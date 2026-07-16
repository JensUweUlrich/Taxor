
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/range/to.hpp>
#include <ranges>
#include <math.h>
#include <algorithm>
#include <tuple>

#include "taxor_profile_configuration.hpp"
#include "taxor_profile.hpp"
#include "search_results.hpp"
#include "parallel_util.hpp"
#include "skim_classifier.hpp"
#include <taxutil.hpp>
#include <profile_output.hpp>

#include <ankerl/unordered_dense.h>

/*!\file taxor_profile.cpp
 * \brief Implements the `taxor profile` subcommand.
 *
 * Turns a `taxor search` results file into a taxonomic abundance profile
 * (CAMI-format report/sequence-abundance/binning files). Two independent
 * read classification strategies are available, selected via
 * `--classification-method`:
 *
 *  - "em" (default): reference filtering to suppress near-identical strain
 *    redundancy (remove_matches_to_nonunique_refs, remove_low_confidence_references,
 *    filter_ref_associations), followed by iterative expectation-maximization
 *    (expectation_maximization) over per-read posterior probabilities.
 *  - "skim": the SKiM binomial-hypothesis-test classifier (see
 *    skim_classifier.hpp), a single-pass per-read statistical test with no
 *    iterative refinement and no cross-read preprocessing.
 *
 * Both paths converge on the same downstream reporting functions
 * (calculate_higher_rank_abundances, calculate_relative_genomic_abundances,
 * and the taxonomy::write_* functions in profile_output.hpp), so their
 * outputs are directly comparable. See tax_profile() for the orchestration
 * of the whole pipeline.
 */

namespace taxor::profile
{

/*!\brief Registers all `taxor profile` command-line options on `parser` and binds them into `config`.
 * \param parser The seqan3 argument parser for the `profile` subcommand.
 * \param config The configuration struct whose fields are bound to the CLI options.
 */
void set_up_subparser_layout(seqan3::argument_parser & parser, taxor::profile::configuration & config)
{
    parser.info.version = "0.2.0";
    parser.info.author = "Jens-Uwe Ulrich";
    parser.info.email = "jens-uwe.ulrich@hpi.de";
    parser.info.short_description = "Taxonomic profiling of a sample by giving read matching results of Taxor search";

    parser.info.description.emplace_back("Taxonomic profiling of the given read set");

    parser.add_subsection("Main options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.search_file,
                      '\0', "search-file", "taxor search file containing results of read querying against the HIXF index",
                      seqan3::option_spec::required);

    parser.add_option(config.report_file, '\0', "cami-report-file", 
                      "output file reporting genomic abundances in CAMI profiling format. "
                      "This is the relative genome abundance in terms of the genome copy number for the respective TAXID in the overall sample. "
                      "Note that this is not identical to the relative abundance in terms of assigned base pairs.",
                      seqan3::option_spec::required);

    parser.add_option(config.sequence_abundance_file, '\0', "seq-abundance-file", 
                      "output file reporting sequence abundance in CAMI profiling format (including unclassified reads). "
                      "This is the relative sequence abundance in terms of sequenced base pairs for the respective TAXID in the overall sample. "
                      "Note that this is not identical to the genomic abundance in terms of genome copy number for the respective TAXID.");

    parser.add_option(config.binning_file, '\0', "binning-file", 
                      "output file reporting read to taxa assignments in CAMI binning format",
                      seqan3::option_spec::required);

    parser.add_option(config.sample_id, '\0', "sample-id", 
                      "Identifier of the analyzed sample",
                      seqan3::option_spec::required);

    parser.add_option(config.threshold, '\0', "min-abundance",
                      "Minimum abundance to report (default: 0.001)",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<double>(0.0), static_cast<double>(1.0)});

    parser.add_option(config.em_steps,
                      '\0', "em-steps",
                      "The number of steps for the expectation maximization (EM) algorithm (default: 100).",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(1000)});

    parser.add_option(config.threads,
                      '\0', "threads",
                      "The number of threads to use.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(32)});

    parser.add_subsection("SKiM classification options (only used with --classification-method skim):");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.classification_method,
                      '\0', "classification-method",
                      "Read classification method: \"em\" (Taxor's own EM-based profiler, default) or "
                      "\"skim\" (the binomial-hypothesis-test classifier from Schneggenburger & Zola, "
                      "\"SKiM: accurately classifying metagenomic ONT reads in limited memory\"). "
                      "Both can be run on the same --search-file (with different output files) to compare them.",
                      seqan3::option_spec::standard,
                      seqan3::value_list_validator{"em", "skim"});

    parser.add_option(config.skim_kmer_size,
                      '\0', "kmer-size",
                      "k-mer size used when building the index (required for --classification-method skim; "
                      "must match the value used for taxor build/search).",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(64)});

    parser.add_option(config.skim_nfixed,
                      '\0', "skim-nfixed",
                      "SKiM's fixed normalization sample size for the binomial test (default: 100, matching the paper).",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(100000)});

    parser.add_option(config.skim_cutoff_exponent,
                      '\0', "skim-cutoff-exponent",
                      "SKiM's p-value cutoff, as a power of ten: cutoff = 10^-e (default: e=12, matching the paper).",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(300)});


    parser.add_flag(config.output_verbose_statistics,
                    '\0', "output-verbose-statistics",
                    "Enable verbose statistics to be "
                    "printed to std::cout. If the flag --determine-best-tmax is not set, this flag is ignored "
                    "and has no effect.",
                    seqan3::option_spec::hidden);

    parser.add_flag(config.debug,
                    '\0', "debug",
                    "Enables debug output in layout file.",
                    seqan3::option_spec::hidden);
}

/*!\brief Splits `str` on every occurrence of `delimiter` into a vector of substrings.
 * \param str The string to split (not modified).
 * \param delimiter The single character to split on.
 * \return The substrings between delimiters, in order; consecutive delimiters produce empty substrings.
 */
std::vector<std::string> str_split(std::string &str, char delimiter)
{

    std::stringstream str_stream(str);
    std::string segment;
    std::vector<std::string> seglist;

    while(std::getline(str_stream, segment, delimiter))
    {
        seglist.push_back(segment);
    }

    return std::move(seglist);
}

/*!\brief Parses a `taxor search` results TSV file into per-read candidate match lists.
 * \param filepath Path to the tab-separated search-results file (see the "search" section of the README for the
 *        10-column format: QUERY_NAME, ACCESSION, REFERENCE_NAME, TAXID, REF_LEN, QUERY_LEN, QHASH_COUNT,
 *        QHASH_MATCH, TAX_STR, TAX_ID_STR; a row with ACCESSION=="-" marks a read with no match at all).
 * \param taxpath[out] Filled with each accession id's (TAX_ID_STR, TAX_STR) lineage pair, the first time that
 *        accession is encountered.
 * \return A map from read id to the list of its candidate reference matches (a single "-" entry if the read has no
 *         match). A read's line(s) are expected to be contiguous in the file; a "-" line is only kept if the read
 *         has no other match recorded yet.
 */
std::map<std::string, std::vector<taxonomy::Search_Result>> parse_search_results(std::string const filepath,
                                                                    std::map<std::string, std::pair<std::string, std::string>> &taxpath)
{
    //std::vector<std::vector<std::string> > tax_file_lines{};
    //taxonomy::read_tsv(filepath, tax_file_lines);
    uint64_t species_counter = 0;
	std::map<std::string, std::vector<taxonomy::Search_Result>> results{};

    size_t idx = 0;

    std::ifstream ifs(filepath);
    if (ifs.fail()) 
    {
        throw std::runtime_error{"Could not open search results file: " + filepath};
    }
    std::string tmpline;

	while (getline(ifs, tmpline)) {

        std::stringstream ss(tmpline);
        std::vector<std::string> line;
        std::string tmp;
        while (getline(ss, tmp, '\t')) 
        {
            line.push_back(tmp);
        }

        if (idx++ == 0) 
            continue;
        std::string read_id = line[0];
        
        if (line[0].find_first_of(" ") != std::string::npos)
            read_id = line[0].substr(0,line[0].find_first_of(" "));
        taxonomy::Search_Result res{};
        res.read_id = read_id;
        if (line[1].compare("-") == 0)
        {
            res.accession_id = "-";
            res.query_len = std::stoull(line[5]);
        }
        else
        {
            res.accession_id = line[1];
            res.tax_id = line[3];
            res.ref_len = std::stoull(line[4]);
            res.query_len = std::stoull(line[5]);
            res.query_hash_count = std::stoull(line[6]);
            res.query_hash_match = std::stoull(line[7]);

            if (!taxpath.contains(res.accession_id))
            {
                std::pair<std::string, std::string> taxpair = std::make_pair(line[9], line[8]);
                taxpath.insert(std::move(std::make_pair(res.accession_id, std::move(taxpair))));
            }
        }

		if (!results.contains(read_id))
        {
            results.insert(std::move(std::make_pair(read_id, std::vector<taxonomy::Search_Result>{})));
        }

        // don't add null result if read already has a reference assignment
        if (results.at(read_id).size() > 0 && res.accession_id.compare("-") == 0)
        {
            continue;
        }
        results.at(read_id).emplace_back(std::move(res));
	}

    return std::move(results);
}


/*!\brief Collects every reference accession that has at least one uniquely-mapping read.
 *
 * "Uniquely mapping" here means a read whose candidate list has exactly one entry (i.e. it
 * matched only this one reference during search) - used as supporting evidence that a
 * reference is a genuine hit and not just along for the ride on another reference's reads.
 *
 * \param search_results Parsed search results (read only).
 * \param threads Number of worker threads to partition the read scan across.
 * \return The set of accession ids with at least one uniquely-mapping read.
 */
ankerl::unordered_dense::set<std::string> get_refs_with_uniquely_mapping_reads(std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results,
                                                                                uint8_t threads = 1u)
{
    auto entries = taxor::util::to_pointer_vector(search_results);
    std::vector<ankerl::unordered_dense::set<std::string>> partials(std::max<size_t>(1, threads));

    auto worker = [&](size_t thread_index, size_t start, size_t end)
    {
        auto & local = partials[thread_index];
        for (size_t idx = start; idx < end; ++idx)
        {
            auto & pair = *entries[idx];
            if (pair.second.size() == 1)
            {
                if (pair.second.at(0).accession_id.compare("-") != 0)
                {
                    local.insert(pair.second.at(0).accession_id);
                }
            }
        }
    };

    taxor::util::parallel_for_indexed(worker, entries.size(), threads);

    // set union across partials; order-independent, so a plain sequential merge is exact
    ankerl::unordered_dense::set<std::string> ref_unique_mappings{};
    for (auto & local : partials)
        for (auto & acc : local)
            ref_unique_mappings.insert(acc);

    return std::move(ref_unique_mappings);
}


/*!\brief Removes ambiguous read-to-reference assignments where the reference has no unique mapping.
 *
 * For each read with more than one candidate reference: if at least one of its candidates is
 * in `ref_unique_mappings` (i.e. it's independently supported by some other, uniquely-mapping
 * read), every candidate *not* in that set is dropped from this read's list, leaving only the
 * independently-supported candidate(s). Reads with no independently-supported candidate at all
 * are left unchanged (still ambiguous). If this empties a read's candidate list, it is replaced
 * with a single "-" (unclassified) placeholder.
 *
 * Each map entry (one read's candidate list) is mutated independently of every other
 * entry, and only reads (never mutates) the shared ref_unique_mappings set, so this is
 * safe to split across threads with no merge step required.
 *
 * \param search_results Search results to filter in place.
 * \param ref_unique_mappings Accession ids known to have at least one uniquely-mapping read (see
 *        get_refs_with_uniquely_mapping_reads); read-only.
 * \param threads Number of worker threads to partition the read scan across.
*/
void remove_matches_to_nonunique_refs(std::map<std::string, std::vector<taxonomy::Search_Result>>& search_results,
                                      ankerl::unordered_dense::set<std::string>& ref_unique_mappings,
                                      uint8_t threads = 1u)
{
    auto entries = taxor::util::to_pointer_vector(search_results);

    auto worker = [&](size_t /*thread_index*/, size_t start, size_t end)
    {
    std::vector<taxonomy::Search_Result>::iterator search_iterator;
    for (size_t idx = start; idx < end; ++idx)
    {
        auto & pair = *entries[idx];
        if (pair.second.size() > 1)
        {

            // first check whether read maps to at least one reference with a uniquely mapping read
            bool unique = false;
            uint64_t query_len = 0;
            for (search_iterator = pair.second.begin(); search_iterator != pair.second.end(); ++search_iterator)
            {
                query_len = (*search_iterator).query_len;
                if (ref_unique_mappings.contains((*search_iterator).accession_id))
                {
                    unique = true;
                    break;
                }
            }

            // only remove ambiguous read-to-reference assignments if at least one ref with a uniquely mapping read matches this read
            if (unique)
            {
                search_iterator = pair.second.begin();
                while (search_iterator != pair.second.end())
                {
                    query_len = (*search_iterator).query_len;
                    if (!ref_unique_mappings.contains((*search_iterator).accession_id))
                        search_iterator = pair.second.erase(search_iterator);
                    else
                        search_iterator++;
                }
            }

            if (pair.second.size() == 0)
            {
                taxonomy::Search_Result res{pair.first,"-",0,query_len,0,0};
                pair.second.emplace_back(std::move(res));
            }
        }
    }
    };

    taxor::util::parallel_for_indexed(worker, entries.size(), threads);
}


/*!\brief Counts, per reference accession, how many reads map to it uniquely vs. ambiguously.
 * \param search_results Parsed search results (read only).
 * \param threads Number of worker threads to partition the read scan across.
 * \return Map from accession id to (unique_read_count, ambiguous_read_count); used by
 *         remove_low_confidence_references to decide which references have enough independent
 *         support to keep.
 */
std::map<std::string,std::pair<uint64_t,uint64_t>> count_unique_ambiguous_mappings_per_reference(
                                std::map<std::string, std::vector<taxonomy::Search_Result>>& search_results,
                                uint8_t threads = 1u)
{
    // <taxid, <unique,ambiguous>>
    auto entries = taxor::util::to_pointer_vector(search_results);
    std::vector<std::map<std::string,std::pair<uint64_t,uint64_t>>> partials(std::max<size_t>(1, threads));

    auto worker = [&](size_t thread_index, size_t start, size_t end)
    {
        auto & map_counts = partials[thread_index];
        for (size_t idx = start; idx < end; ++idx)
        {
            auto & pair = *entries[idx];
            if (pair.second.size() == 1)
            {
                if (pair.second.at(0).accession_id.compare("-") != 0)
                {
                    if (!map_counts.contains(pair.second.at(0).accession_id))
                    {
                        map_counts.insert(std::move(std::make_pair(pair.second.at(0).accession_id, std::move(std::make_pair(0,0)))));
                    }
                    map_counts.at(pair.second.at(0).accession_id).first += 1;
                }
            }
            else
            {
                for (auto & res : pair.second)
                {
                    if (res.accession_id.compare("-") == 0) continue;
                    if (!map_counts.contains(res.accession_id))
                    {
                        map_counts.insert(std::move(std::make_pair(res.accession_id, std::move(std::make_pair(0,0)))));
                    }
                    map_counts.at(res.accession_id).second += 1;
                }
            }
        }
    };

    taxor::util::parallel_for_indexed(worker, entries.size(), threads);

    // additive reduction: commutative/associative, so merge order doesn't affect the result
    std::map<std::string,std::pair<uint64_t,uint64_t>> map_counts{};
    for (auto & partial : partials)
    {
        for (auto & entry : partial)
        {
            auto & dst = map_counts[entry.first];
            dst.first += entry.second.first;
            dst.second += entry.second.second;
        }
    }

    return std::move(map_counts);
}

/*!\brief Keeps only references with enough independently-supported (unique) reads, then re-filters ambiguous reads accordingly.
 *
 * A reference is "accepted" if it has at least `min_unique_mappings` uniquely-mapping reads AND
 * unique reads make up at least `min_fraction_unique` of all its reads (unique + ambiguous).
 * Ambiguous read candidates pointing at a non-accepted reference are then dropped via
 * remove_matches_to_nonunique_refs, reusing the exact same "keep only accepted candidates" logic
 * as the first filtering round.
 *
 * \param search_results Search results to filter in place.
 * \param map_counts Per-reference unique/ambiguous read counts, from count_unique_ambiguous_mappings_per_reference.
 * \param min_unique_mappings Minimum number of uniquely-mapping reads a reference needs to be accepted.
 * \param min_fraction_unique Minimum fraction of a reference's reads that must be unique (not ambiguous) to be accepted.
 * \param threads Number of worker threads to partition the read scan across.
 */
void remove_low_confidence_references(std::map<std::string, std::vector<taxonomy::Search_Result>>& search_results,
                                      std::map<std::string,std::pair<uint64_t,uint64_t>>& map_counts,
                                      uint8_t min_unique_mappings,
                                      float min_fraction_unique,
                                      uint8_t threads = 1u)
{
    ankerl::unordered_dense::set<std::string> accepted_refs{};
    for (auto & ref : map_counts)
    {
        if (ref.second.first >= min_unique_mappings &&
                static_cast<float>(ref.second.first) / static_cast<float>(ref.second.first + ref.second.second) >= min_fraction_unique)
            accepted_refs.insert(ref.first);
    }
    remove_matches_to_nonunique_refs(search_results, accepted_refs, threads);
}


/*!\brief Detects references that "explain" reads better than a near-identical strain, and merges/reassigns accordingly.
 *
 * Filters out suspicious mappings using an approach similar to the two-stage taxonomy
 * assignment algorithm in MegaPath (Leung et al., 2020). In outline:
 *  1. Build, per reference, its unique/total read counts and (for ambiguous reads) which other
 *     references it co-occurs with and how often (Ref_Map_Info::associated_species).
 *  2. For each pair of co-occurring references, decide (based on relative read support and
 *     >=95% mapping overlap) whether one is "explained by" the other, recording that in
 *     `explained_refs`.
 *  3. Flatten transitive "explained by" chains (bounded-pass, cycle-safe - see the loop's own
 *     comment for why an unbounded fixed-point loop is unsafe here).
 *  4. Reassign/drop ambiguous reads' ("size() > 1") ambiguous candidates according to
 *     `explained_refs`, and drop explained references from the returned taxa map.
 *
 * \param search_results Search results to filter/reassign in place.
 * \param threads Number of worker threads to partition the read scan across (used for the
 *        per-read passes; the reference-level "explains" computation is sequential, since it's
 *        bounded by the number of distinct references, not reads).
 * \return Map from surviving accession id to its reference length, for use as the taxa universe
 *         in downstream abundance calculation (e.g. expectation_maximization).
*/
std::map<std::string, size_t> filter_ref_associations(std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results, uint8_t threads)
{

    // taxid, length
    std::map<std::string, size_t> taxa_lengths{};
    std::map<std::string, Ref_Map_Info> ref_associations{};

    // Build ref_associations/taxa_lengths from search_results. Each read is processed
    // independently, but different reads can update the SAME reference's Ref_Map_Info
    // (and taxa_lengths entry), so this is a real reduction: every thread accumulates
    // into its own partial maps, which are then merged (by addition) sequentially below.
    {
        auto entries = taxor::util::to_pointer_vector(search_results);
        std::vector<std::map<std::string, size_t>> partial_lengths(std::max<size_t>(1, threads));
        std::vector<std::map<std::string, Ref_Map_Info>> partial_assoc(std::max<size_t>(1, threads));

        auto worker = [&](size_t thread_index, size_t start, size_t end)
        {
            auto & taxa_lengths_local = partial_lengths[thread_index];
            auto & ref_associations_local = partial_assoc[thread_index];

            for (size_t idx = start; idx < end; ++idx)
            {
                auto & pair = *entries[idx];
                if (pair.second.size() == 0) continue;
                if (pair.second.size() == 1)
                {
                    // if there is one unique mapping between this read and a reference
                    if (pair.second.at(0).accession_id.compare("-") != 0)
                    {
                        if (!ref_associations_local.contains(pair.second.at(0).accession_id))
                            ref_associations_local.insert(std::move(std::make_pair(pair.second.at(0).accession_id, Ref_Map_Info{})));

                        // increment the number of uniquely mapped reads by 1
                        ref_associations_local.at(pair.second.at(0).accession_id).unique_assign_reads += 1;
                        // increment the number of all mapped reads by 1
                        ref_associations_local.at(pair.second.at(0).accession_id).all_assigned_reads += 1;

                        if (!taxa_lengths_local.contains(pair.second.at(0).accession_id))
                            taxa_lengths_local.insert(std::move(std::make_pair(pair.second.at(0).accession_id, pair.second.at(0).ref_len)));
                    }
                }
                else
                {
                    // collect all references assigned to that read
                    std::vector<std::string> acc_ids{};
                    for (auto & res : pair.second)
                    {
                        if (res.accession_id.compare("-") == 0) continue;
                        if (!ref_associations_local.contains(res.accession_id))
                            ref_associations_local.insert(std::move(std::make_pair(res.accession_id, Ref_Map_Info{})));

                        acc_ids.emplace_back(res.accession_id);
                        // increment the number of all mapped reads by 1
                        ref_associations_local.at(res.accession_id).all_assigned_reads += 1;

                        if (!taxa_lengths_local.contains(res.accession_id))
                            taxa_lengths_local.insert(std::move(std::make_pair(res.accession_id, res.ref_len)));
                    }

                    // iterate over all references assigned with that read
                    // (iterate by reference: acc_ids can be large for reads with many ambiguous
                    // matches, and this loop is already O(k^2) per read, so avoid an extra string
                    // copy of every element on top of that)
                    for (std::string const & acc_id1 : acc_ids)
                    {
                        Ref_Map_Info & info1 = ref_associations_local.at(acc_id1);
                        for (std::string const & acc_id2 : acc_ids)
                        {
                            if (acc_id1.compare(acc_id2) == 0) continue;

                            // increment the number of shared reads between ref1 and ref2
                            // (single lookup instead of contains()+insert()+at() against the same key)
                            ++info1.associated_species[acc_id2];
                        }
                    }

                }
            }
        };

        taxor::util::parallel_for_indexed(worker, entries.size(), threads);

        // Merge: all operations below are addition (commutative/associative), so the
        // final numeric result does not depend on partition boundaries or merge order.
        for (auto & partial : partial_lengths)
        {
            for (auto & entry : partial)
            {
                auto it = taxa_lengths.find(entry.first);
                if (it == taxa_lengths.end())
                    taxa_lengths.insert(entry);
                // taxa_lengths is conceptually "first value wins" (ref_len is a property of
                // the reference genome, not of the read), so a differing length for the same
                // accession id indicates malformed/inconsistent input - surface it loudly
                // instead of silently picking whichever partition happened to run first.
                else if (it->second != entry.second)
                    std::cerr << "[TAXOR PROFILE WARNING] inconsistent reference length for accession "
                              << entry.first << ": " << it->second << " vs " << entry.second << std::endl;
            }
        }
        for (auto & partial : partial_assoc)
        {
            for (auto & entry : partial)
            {
                Ref_Map_Info & dst = ref_associations[entry.first];
                dst.unique_assign_reads += entry.second.unique_assign_reads;
                dst.all_assigned_reads += entry.second.all_assigned_reads;
                for (auto & assoc : entry.second.associated_species)
                    dst.associated_species[assoc.first] += assoc.second;
            }
        }
    }

    // first is explained by second => exchange first with second
    std::map<std::string, std::string> explained_refs{};
    // iterate over all found references
    for (auto & ref : ref_associations)
    {
        // find associated references that explain by ref
        for (auto & assoc_ref : ref.second.associated_species)
        {
            // use the ref with more unique mappings as master
            // only use ref if it has more uniquely mapping reads or overall more mapped reads than the associated species
            if (ref.second.unique_assign_reads > ref_associations.at(assoc_ref.first).unique_assign_reads || ref.second.all_assigned_reads > ref_associations.at(assoc_ref.first).all_assigned_reads)
            {
                // if more than 95% of ref's mapped reads also map to current associated ref
                if (ref.second.all_assigned_reads - assoc_ref.second < static_cast<uint64_t>(0.05 * static_cast<double>(ref.second.all_assigned_reads)))
                {
                    // if the number of uniquely mapped reads of ref is less then 5% the number of uniquely mapped reads of associated ref
                    //if (ref.second.unique_assign_reads < static_cast<uint64_t>(0.05 * static_cast<double>(ref_associations.at(assoc_ref.first).unique_assign_reads)))
                    //{
                        explained_refs.insert(std::move(std::make_pair(ref.first, assoc_ref.first)));
                    //}
                }
            }
            // use associated ref as explained ref if it has more uniquely mapping reads and overall more mapped reads
            else
            {
                 // if more than 95% of ref's mapped reads also map to current associated ref
                if (ref_associations.at(assoc_ref.first).all_assigned_reads - ref_associations.at(assoc_ref.first).associated_species.at(ref.first) < static_cast<uint64_t>(0.05 * static_cast<double>(ref_associations.at(assoc_ref.first).all_assigned_reads)))
                {
                    // if the number of uniquely mapped reads of ref is less then 5% the number of uniquely mapped reads of associated ref
                    //if (ref.second.unique_assign_reads < static_cast<uint64_t>(0.05 * static_cast<double>(ref_associations.at(assoc_ref.first).unique_assign_reads)))
                    //{
                        explained_refs.insert(std::move(std::make_pair(assoc_ref.first, ref.first)));
                }
            }
        }
    }
   
    // Flatten transitive "explained by" chains (X explained by Y, Y explained by Z => X explained by Z).
    // The per-pair heuristic above is not guaranteed to be globally acyclic (e.g. three mutually
    // near-identical references can end up with a 3-way cycle A->B->C->A), and a plain
    // "repeat until no changes" loop never converges for such a cycle: two of the three nodes settle
    // into a stable mutual pair while the third perpetually flips between them, looping forever on a
    // single thread. Bound the number of passes: any acyclic chain fully flattens within
    // explained_refs.size() passes, so exceeding that many passes only happens in the presence of a
    // cycle, at which point we stop rather than spin indefinitely.
    bool found = true;
    size_t const max_passes = explained_refs.size() + 1;
    size_t pass = 0;

    while (found && pass < max_passes)
    {
        found = false;
        for (auto & exp_ref : explained_refs)
        {
            if (explained_refs.contains(exp_ref.second) && exp_ref.first.compare(explained_refs.at(exp_ref.second)) != 0)
            {
                exp_ref.second = explained_refs.at(exp_ref.second);
                found = true;
            }
        }
        ++pass;
    }

    
    // iterate over search results and filter results of ambiguous mappings
    // reassign unique mappings of refs explained by another ref => TODO: this should not be done, keep unique mappings
    // in such cases,
    // Independent per read: each thread only mutates its own read's vector, and
    // explained_refs/taxa_lengths are both fully built and read-only from here on, so
    // this is safe to split across threads with no merge step required.
    {
        auto entries = taxor::util::to_pointer_vector(search_results);

        auto worker = [&](size_t /*thread_index*/, size_t start, size_t end)
        {
            for (size_t idx = start; idx < end; ++idx)
            {
                auto & pair = *entries[idx];
                if (pair.second.size() == 0) continue;
                if (pair.second.size() == 1) continue;
/*        {
            // if unique mapping is explained by another ref
            if (explained_refs.contains(pair.second.at(0).accession_id))
            {
                pair.second.at(0).accession_id = explained_refs.at(pair.second.at(0).accession_id);
                pair.second.at(0).ref_len = taxa_lengths.at(pair.second.at(0).accession_id);
            }
        }
*/
//        else
//        {
                // collect all references assigned to that read
                std::set<std::string> acc_ids{};
                for (auto & res : pair.second)
                    acc_ids.emplace(res.accession_id);

                std::vector<taxonomy::Search_Result>::iterator it = pair.second.begin();
                // iterate over search results
                while (it != pair.second.end())
                {
                    // if accession id is explained by another reference
                    if (explained_refs.contains((*it).accession_id))
                    {
                        // if there is another mapping of this read to the reference that explains current accession id
                        // remove this match
                        if (acc_ids.contains(explained_refs.at((*it).accession_id)))
                        {
                            it = pair.second.erase(it);
                            continue;
                        }
                        // reassign otherwise
                        else
                        {
                            (*it).accession_id = explained_refs.at((*it).accession_id);
                            (*it).ref_len = taxa_lengths.at((*it).accession_id);
                        }
                    }
                    it++;
                }
//        }
            }
        };

        taxor::util::parallel_for_indexed(worker, entries.size(), threads);
    }

    std::map<std::string, size_t>::iterator it = taxa_lengths.begin();
    while (it != taxa_lengths.end())
    {
        if (explained_refs.contains((*it).first))
        {
            it = taxa_lengths.erase(it);
            continue;
        }
        it++;
    }
    
    return std::move(taxa_lengths);
}

/*!\brief Initializes each taxon's EM prior to a uniform log-probability (log(1/N) for N taxa).
 * \param taxa The universe of candidate taxa (accession id -> reference length); only the keys are used.
 * \return Map from accession id to its initial log-prior, the EM loop's starting point before any iteration.
 */
std::map<std::string, double> initialize_prior_log_probabilities(std::map<std::string, size_t>& taxa)
{
    std::map<std::string, double> priors {};
    for (auto & taxon : taxa)
    {
        priors.insert(std::move(std::make_pair(taxon.first, log(1.0 / static_cast<double>(taxa.size())))));
    }
    return std::move(priors);
}

/*!\brief Computes each read's log-likelihood of originating from each of its candidate references.
 *
 * For a read with more than one candidate, the likelihood of each candidate is its
 * match-count ratio (query_hash_match/query_hash_count) normalized against the sum of all
 * candidates' ratios, in log space. A read with exactly one (non-"-") candidate gets a
 * likelihood of 0 (i.e. probability 1, log(1)=0), since there's nothing to normalize against.
 * "-" (no-match) candidates never contribute a likelihood entry.
 *
 * \param search_results Parsed search results (read only).
 * \param threads Number of worker threads to partition the read scan across.
 * \return Map from read id to a map of {candidate accession id -> log-likelihood}, consumed by
 *         expectation_maximization's per-read posterior computation.
 */
std::map<std::string, std::map<std::string, double>> calculate_log_likelihoods(std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results,
                                                                                uint8_t threads = 1u)
{
    // Every entry is written under its own read id (pair.first), which is unique per
    // std::map entry by definition, so keys are provably disjoint across partitions -
    // the merge below can never see a real conflict, no reduction logic needed.
    auto entries = taxor::util::to_pointer_vector(search_results);
    std::vector<std::map<std::string, std::map<std::string, double>>> partials(std::max<size_t>(1, threads));

    auto worker = [&](size_t thread_index, size_t start, size_t end)
    {
        auto & likelihoods = partials[thread_index];
        for (size_t idx = start; idx < end; ++idx)
        {
            auto & pair = *entries[idx];
            std::map<std::string, double> read_ref_liklihoods{};
            if (pair.second.size() == 0) continue;
            if (pair.second.size() > 1)
            {
                // calculate sum of matchcount ratios
                double sum_ratio{0.0};
                for (auto & res : pair.second)
                {
                    if (res.accession_id.compare("-") == 0) continue;
                    sum_ratio += static_cast<double>(res.query_hash_match) / static_cast<double>(res.query_hash_count);
                }

                // calculate log likelihoods of the single matches
                for (auto & res : pair.second)
                {
                    if (res.accession_id.compare("-") == 0) continue;
                    double likeL = (log(static_cast<double>(res.query_hash_match)) - log(static_cast<double>(res.query_hash_count))) - log(sum_ratio);
                    read_ref_liklihoods.insert(std::move(std::make_pair(res.accession_id, likeL)));
                }
            }
            else
            {
                // if there is one unique mapping between this read and a reference
                if (pair.second.at(0).accession_id.compare("-") != 0)
                {
                    read_ref_liklihoods.insert(std::move(std::make_pair(pair.second.at(0).accession_id, 0.0)));
                }
            }

            likelihoods.insert(std::move(std::make_pair(pair.first, read_ref_liklihoods)));
        }
    };

    taxor::util::parallel_for_indexed(worker, entries.size(), threads);

    std::map<std::string, std::map<std::string, double>> likelihoods{};
    for (auto & partial : partials)
        likelihoods.merge(partial);

    return std::move(likelihoods);
}

/*!\brief Recomputes each taxon's sequence (nucleotide) abundance from a fixed read classification.
 *
 * For each read currently assigned to a taxon in `profile_results`, credits that taxon with the
 * read's query length; a taxon's log-prior becomes log(its total credited bases / all reads'
 * total bases). This only reads `taxa`'s keys (to pre-seed the accumulator) and `profile_results`
 * - it does not use the *incoming* `log_priors` values in its computation, so it is equally
 * correct whether called once (a one-shot "abundance from a fixed classification", as
 * tax_profile's SKiM path does) or repeatedly inside an EM loop (as expectation_maximization
 * does, where `profile_results` changes between calls).
 *
 * \param log_priors[out] Overwritten (in log space) with each taxon's updated abundance.
 * \param taxa The universe of candidate taxa (accession id -> reference length); only the keys are used.
 * \param profile_results Current read -> best-match assignment (see expectation_maximization / skim::classify_reads).
 * \param threads Number of worker threads to partition the read scan across.
 * \return The (log-space) unclassified abundance, i.e. log(total bases of unclassified reads / all reads' total bases).
 */
double update_log_prior_probabilities(std::map<std::string, double> &log_priors,
                                    std::map<std::string, size_t> & taxa,
                                    std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results,
                                    uint8_t threads = 1u)
{
    std::map<std::string, uint64_t> ref_nts{};
    for (auto & t : taxa)
        ref_nts.insert(std::move(std::make_pair(t.first, 0)));

    // for each taxon sum all lengths of matching reads
    // Different reads can update the SAME taxon's entry in ref_nts, so this is a real
    // reduction: each thread accumulates into its own pre-seeded partial ref_nts map
    // (mirroring the pre-seed-from-taxa invariant .at() below relies on) plus local
    // scalar counters, merged by integer += after joining (exact, no FP concern).
    size_t all_nts{0};
    size_t unclassified_nts{0};
    //std::cout << "Sum nts of matching reads per taxon..." << std::endl << std::flush;
    {
        auto entries = taxor::util::to_pointer_vector(profile_results);
        std::vector<std::map<std::string, uint64_t>> partial_ref_nts(std::max<size_t>(1, threads), ref_nts);
        std::vector<size_t> partial_all_nts(std::max<size_t>(1, threads), 0);
        std::vector<size_t> partial_unclassified_nts(std::max<size_t>(1, threads), 0);

        auto worker = [&](size_t thread_index, size_t start, size_t end)
        {
            auto & ref_nts_local = partial_ref_nts[thread_index];
            for (size_t idx = start; idx < end; ++idx)
            {
                auto & read = *entries[idx];
                if (read.second.size() == 0) continue;
                partial_all_nts[thread_index] += read.second.at(0).query_len;
                if (read.second.at(0).accession_id.compare("-") == 0)
                {
                    partial_unclassified_nts[thread_index] += read.second.at(0).query_len;
                    continue;
                }
                for (auto &ref : read.second)
                {
                    ref_nts_local.at(ref.accession_id) += ref.query_len;
                }
            }
        };

        taxor::util::parallel_for_indexed(worker, entries.size(), threads);

        for (size_t i = 0; i < partial_ref_nts.size(); ++i)
        {
            all_nts += partial_all_nts[i];
            unclassified_nts += partial_unclassified_nts[i];
            for (auto & t : ref_nts)
                t.second += partial_ref_nts[i].at(t.first);
        }
    }
    std::cout << "done" << std::endl;
    // calculate average depth of coverage for each taxon
    // sum of all matched read lengths divided by length of taxon reference sequence
    /*double sum_avg_cov{0.0};
    for (auto & t : ref_nts)
    {
        log_priors.at(t.first) = static_cast<double>(t.second) / static_cast<double>(taxa.at(t.first));
        sum_avg_cov += log_priors.at(t.first);
    }
    */

    // calculate relative genomic abundance for each taxon
    // divide average coverage of each taxon by the sum of average coverage of all taxa
    std::cout << "final update..." << std::endl << std::flush;
    for (auto &t : log_priors)
    {
        // only for genome relative abundance
        //log_priors.at(t.first) = log(t.second + 0.000000000001) - log(sum_avg_cov);

        // 2nd version uses nucleotide abundance
        log_priors.at(t.first) = log(static_cast<double>(ref_nts.at(t.first)) + 0.000000000001) - log(static_cast<double>(all_nts));
    }
    std::cout << "done" << std::endl << std::flush;
    return log(static_cast<double>(unclassified_nts) + 0.000000000001) - log(static_cast<double>(all_nts));

}

// GTDB-style taxon names are prefixed with a 3-character rank marker (e.g. "p__", "s__");
// strip it if present, otherwise return the name unchanged (handles empty/short/malformed entries)
static std::string strip_rank_prefix(std::string const & taxname)
{
    return taxname.size() > 3 ? taxname.substr(3) : std::string{};
}

/*!\brief Rolls up per-taxon abundances to every higher rank (genus, family, ..., superkingdom).
 *
 * For each taxon with nonzero abundance, walks its GTDB-style `;`-separated lineage string
 * (from `taxpath`, keyed by accession id) and, for every rank along that lineage, accumulates
 * the taxon's abundance into that rank's Profile_Output entry (creating it on first encounter).
 * A taxon's rank is read from its own lineage entry's single-character prefix (s/g/f/o/c/p/k),
 * not from its position in the lineage string, so this is robust to a lineage missing a level.
 *
 * \param species_abundances Map from accession id (or "unclassified") to its relative abundance;
 *        entries with abundance 0 are skipped, "unclassified" is passed through as its own entry.
 * \param taxpath Map from accession id to its (TAX_ID_STR, TAX_STR) lineage pair (from parse_search_results).
 * \return Map from taxon id (at every rank) to its rolled-up Profile_Output (rank, lineage strings, percentage).
 */
std::map<std::string, taxonomy::Profile_Output> calculate_higher_rank_abundances(std::map<std::string, double> &species_abundances,
                                                            std::map<std::string, std::pair<std::string, std::string>> &taxpath)
{
    std::map<std::string, taxonomy::Profile_Output> rank_profiles{};
    for (auto &sp : species_abundances)
    {
        //std::cout << sp.first << "\t" << sp.second << std::endl;
        if (sp.second == 0) continue;
        //std::cout << std::string(sp.first) << "\t" << sp.second << std::endl;
        if (sp.first.compare("unclassified") == 0)
        {
            taxonomy::Profile_Output profile{};
            profile.taxid = sp.first;
            profile.percentage = sp.second;
            rank_profiles.insert(std::move(std::pair(sp.first, std::move(profile))));
            continue;
        }
        //std::cout << sp.first << "\t" << sp.second << std::endl;
        std::vector<std::string> taxid_path = str_split(taxpath.at(sp.first).first,';');
        std::vector<std::string> taxname_path = str_split(taxpath.at(sp.first).second,';');

        // taxid_path and taxname_path are expected to be parallel (same length), but the
        // underlying taxonomy input is user-provided and may be malformed/incomplete;
        // guard against index-out-of-range access into taxname_path
        auto tname = [&](size_t i) -> std::string
        {
            return i < taxname_path.size() ? taxname_path[i] : std::string{};
        };

        for (size_t index = 0; index < taxid_path.size();++index)
        {

            if (taxid_path[index].size() < 1) continue;
            if (!rank_profiles.contains(taxid_path[index]))
            {
                taxonomy::Profile_Output profile{};
                profile.taxid = taxid_path[index];
                profile.taxid_string = taxid_path[0];
                profile.taxname_string = strip_rank_prefix(tname(0));
                for (size_t index2 = 1; index2 <= index; ++index2)
                {
                    profile.taxid_string += "|";
                    profile.taxid_string += taxid_path[index2];
                    profile.taxname_string += "|";
                    profile.taxname_string += strip_rank_prefix(tname(index2));
                }

                profile.percentage = 0.0;

                std::string const rank_prefix = tname(index).substr(0, 1);
                if (rank_prefix.compare("s") == 0)
				    profile.rank = "species";
			    else if(rank_prefix.compare("g") == 0)
					profile.rank = "genus";
				else if(rank_prefix.compare("f") == 0)
					profile.rank = "family";
				else if(rank_prefix.compare("o") == 0)
					profile.rank = "order";
				else if(rank_prefix.compare("c") == 0)
					profile.rank = "class";
				else if(rank_prefix.compare("p") == 0)
					profile.rank = "phylum";
				else if(rank_prefix.compare("k") == 0)
					profile.rank = "superkingdom";


                rank_profiles.insert(std::move(std::pair(taxid_path[index], std::move(profile))));
            }

            rank_profiles.at(taxid_path[index]).percentage += species_abundances.at(sp.first);
        }
    }

    return std::move(rank_profiles);
    
}

/*!\brief Taxor's EM-based read classifier: iteratively refines per-taxon abundances and per-read best matches.
 *
 * Each iteration: recomputes per-read log-likelihoods (calculate_log_likelihoods), then for
 * every read combines each candidate's likelihood with its taxon's current log-prior to get a
 * posterior; the candidate(s) with the highest posterior become that read's current best match
 * (profile_results), and the single candidate with the *lowest* posterior is pruned from the
 * read's candidate list (so, over iterations, each read's list shrinks towards its best match).
 * Priors are then updated from the new best-match assignment (update_log_prior_probabilities).
 * Stops after `iterations` steps, or earlier once the total log-likelihood's improvement drops
 * below a fixed convergence threshold.
 *
 * \param iterations Maximum number of EM iterations (see --em-steps).
 * \param taxa The universe of candidate taxa (accession id -> reference length), from filter_ref_associations.
 * \param search_results Parsed (and pre-filtered) search results; each read's candidate list is
 *        progressively pruned in place as the EM loop converges.
 * \param profile_results[out] Final read -> best-match assignment (one entry per read, single
 *        "-" placeholder for reads that never scored against any candidate).
 * \param threads Number of worker threads to partition the read scan across.
 * \return Map from accession id (plus "unclassified") to its final relative sequence abundance.
 */
std::map<std::string, double> expectation_maximization(size_t iterations,
                              std::map<std::string, size_t> & taxa,
                              std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results,
                              std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results,
                              uint8_t threads = 1u)
{
    std::cout << "started" << std::endl;
    std::cout << "Initialize prior probabilities ..." << std::flush;
    std::map<std::string, double> log_priors = initialize_prior_log_probabilities(taxa);
    std::cout << "done" << std::endl;
    double cond_log_likelihood = -__DBL_MAX__;
    size_t iter_step = 0;
    double unclassified_abundance{0.0};

    // Reused across iterations: search_results doesn't change size/structure between EM
    // steps (only individual reads' candidate vectors shrink via erase(worst_match)), so
    // the same pointer partitioning stays valid for the whole loop.
    auto entries = taxor::util::to_pointer_vector(search_results);

    while(iter_step < iterations)
    {
        std::cout << "Starting EM iteration " << iter_step << std::endl << std::flush;
        std::cout << "Calculate Log Likelihoods ..." << std::flush;
        std::map<std::string, std::map<std::string, double>> log_likelihoods = calculate_log_likelihoods(search_results, threads);
        std::cout << "done" << std::endl;
        profile_results.clear();

        // Per-read best/worst-match selection is independent (each thread only mutates
        // its own read's vector). profile_results insertion is disjoint-key (by read id),
        // so it's merged cheaply below. new_cond_log_likelihood is a real reduction: each
        // thread sums into its own slot, added together after joining. Note this changes
        // floating-point summation order vs. the sequential run (~1e-12 relative noise,
        // expected and harmless) - but the per-read best_match/worst_match decisions
        // below only ever depend on log_likelihoods/log_priors, never on this running
        // total, so classification output (--binning-file) is unaffected by the reorder;
        // only the coarse EM convergence check below sees any noise, and its threshold
        // is far too coarse to be affected in practice.
        std::vector<std::map<std::string, std::vector<taxonomy::Search_Result>>> partial_profile_results(std::max<size_t>(1, threads));
        std::vector<double> partial_new_cond_log_likelihood(std::max<size_t>(1, threads), 0.0);

        auto worker = [&](size_t thread_index, size_t start, size_t end)
        {
            auto & profile_results_local = partial_profile_results[thread_index];
            double & new_cond_log_likelihood_local = partial_new_cond_log_likelihood[thread_index];

            for (size_t idx = start; idx < end; ++idx)
            {
                auto & read = *entries[idx];
                if (read.second.size() == 0) continue;
                double max_post = -__DBL_MAX__;
                double min_post = __DBL_MAX__;
                std::vector<taxonomy::Search_Result> best_match{};
                // default-initialize to a valid (dereferenceable) iterator rather than a singular
                // one: if none of the candidates below ever get scored (e.g. all of them were
                // filtered out upstream and are missing from log_likelihoods/log_priors), the
                // erase(worst_match) call further down must still operate on a valid iterator
                std::vector<taxonomy::Search_Result>::iterator worst_match = read.second.begin();
                std::vector<taxonomy::Search_Result>::iterator it = read.second.begin();
                // iterate over search results
                while (it != read.second.end())
                {

                    if ((*it).accession_id.compare("-") == 0)
                    {
                        if (read.second.size() == 1)
                        {
                            best_match.emplace_back((*it));
                            break;
                        }
                        else
                        {
                            // "-" is a no-match placeholder, not a real candidate taxon: mark it
                            // as the (default) worst match and move on without falling through to
                            // the likelihood lookup below, which would otherwise dereference `it`
                            // after it was just advanced (UB / segfault if "-" was the last entry)
                            worst_match = it;
                            ++it;
                            continue;
                        }
                    }

                    double post = 0.0;
                    if (log_likelihoods.contains(read.first) &&
                        log_likelihoods.at(read.first).contains((*it).accession_id) &&
                        log_priors.contains((*it).accession_id))
                    {
                        post = log_likelihoods.at(read.first).at((*it).accession_id) + log_priors.at((*it).accession_id);
                    }
                    else
                    {
                        it++;
                        continue;
                    }

                    new_cond_log_likelihood_local += post;
                    // update best match based on best posterior probability
                    if (post >= max_post)
                    {
                        if (post > max_post)
                        {
                            max_post = post;
                            best_match.clear();
                        }
                        best_match.emplace_back((*it));
                    }
                    // update worst match based on worst posterior probability
                    if (post < min_post)
                    {
                        worst_match = it;
                    }
                    it++;

                }
                profile_results_local.insert(std::move(std::make_pair(read.first, std::move(best_match))));
                // remove reference match with worst posterior probability in each iteration until convergence
                if (read.second.size() > 1)
                    read.second.erase(worst_match);
            }
        };

        taxor::util::parallel_for_indexed(worker, entries.size(), threads);

        double new_cond_log_likelihood = 0.0;
        for (size_t i = 0; i < partial_profile_results.size(); ++i)
        {
            profile_results.merge(partial_profile_results[i]);
            new_cond_log_likelihood += partial_new_cond_log_likelihood[i];
        }

        // update referencs abundances (priors)
        std::cout << "Update prior probabilities ..." << std::flush;
        unclassified_abundance = update_log_prior_probabilities(log_priors, taxa, profile_results, threads);
         std::cout << "done" << std::endl << std::flush;
        double diff = new_cond_log_likelihood - cond_log_likelihood;
        if (diff < abs(log(0.0001)))
            break;

        cond_log_likelihood = new_cond_log_likelihood;
        iter_step++;
    }
    std::cout << "Number of EM steps needed: " << iter_step << std::endl << std::flush;

    log_priors.insert(std::move(std::make_pair("unclassified", unclassified_abundance)));
    for (auto & t : log_priors)
    {
        log_priors.at(t.first) = exp(t.second);
    }

    return std::move(log_priors);
}

/*!\brief Computes each taxon's relative *genomic* abundance (genome copy number), as opposed to sequence abundance.
 *
 * Unlike update_log_prior_probabilities's per-base sequence abundance, this divides each
 * taxon's total credited bases by its own reference length first (an average depth-of-coverage
 * proxy), then normalizes those coverage values against each other - so a short genome
 * covered by few reads can rank comparably to a long genome covered by many, reflecting genome
 * *copies* in the sample rather than raw sequenced bases. Unclassified reads are excluded
 * entirely (there is no "unclassified" entry in the output, unlike the sequence-abundance path).
 *
 * \param log_priors[out] Cleared and refilled with each taxon's relative genomic abundance (linear space).
 * \param taxa The universe of candidate taxa (accession id -> reference length).
 * \param profile_results Final read -> best-match assignment (see expectation_maximization / skim::classify_reads).
 * \param threads Number of worker threads to partition the read scan across.
 */
void calculate_relative_genomic_abundances(std::map<std::string, double> &log_priors,
                                           std::map<std::string, size_t> & taxa,
                                           std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results,
                                           uint8_t threads = 1u)
{
    log_priors.clear();
    std::map<std::string, uint64_t> ref_nts{};
    for (auto & t : taxa)
    {
        ref_nts.insert(std::move(std::make_pair(t.first, 0)));
        log_priors.insert(std::move(std::make_pair(t.first, 0.0)));
    }

    // for each taxon sum all lengths of matching reads
    // Different reads can update the SAME taxon's entry in ref_nts, so this is a real
    // reduction: each thread accumulates into its own pre-seeded partial ref_nts map,
    // merged by integer += after joining (exact, no FP concern) - same pattern as
    // update_log_prior_probabilities.
    {
        auto entries = taxor::util::to_pointer_vector(profile_results);
        std::vector<std::map<std::string, uint64_t>> partial_ref_nts(std::max<size_t>(1, threads), ref_nts);

        auto worker = [&](size_t thread_index, size_t start, size_t end)
        {
            auto & ref_nts_local = partial_ref_nts[thread_index];
            for (size_t idx = start; idx < end; ++idx)
            {
                auto & read = *entries[idx];
                if (read.second.size() == 0) continue;
                if (read.second.at(0).accession_id.compare("-") == 0)
                    continue;
                for (auto &ref : read.second)
                {
                    if (ref_nts_local.contains(ref.accession_id))
                        ref_nts_local.at(ref.accession_id) += ref.query_len;
                }
            }
        };

        taxor::util::parallel_for_indexed(worker, entries.size(), threads);

        for (auto & partial : partial_ref_nts)
            for (auto & t : ref_nts)
                t.second += partial.at(t.first);
    }


    // calculate average depth of coverage for each taxon
    // sum of all matched read lengths divided by length of taxon reference sequence
    double sum_avg_cov{0.0};
    for (auto & t : ref_nts)
    {
        if (log_priors.contains(t.first) && taxa.contains(t.first))
        {
            log_priors.at(t.first) = static_cast<double>(t.second) / static_cast<double>(taxa.at(t.first));
            sum_avg_cov += log_priors.at(t.first);
        }
    }
    

    // calculate relative genomic abundance for each taxon
    // divide average coverage of each taxon by the sum of average coverage of all taxa
    for (auto &t : log_priors)
    {
        // only for genome relative abundance
        log_priors.at(t.first) = log(t.second + 0.000000000001) - log(sum_avg_cov);
    }

    for (auto & t : log_priors)
    {
        log_priors.at(t.first) = exp(t.second);
    }

}

/*!\brief Orchestrates the full `taxor profile` pipeline: parse -> classify -> report.
 *
 * Parses the search-results file once, then branches on `config.classification_method`:
 *  - "skim": skim::classify_reads directly on the parsed results, followed by a single
 *    (non-iterative) call to update_log_prior_probabilities to get sequence abundances - see
 *    skim_classifier.hpp for why one call is mathematically sufficient here.
 *  - "em" (default): the two rounds of reference filtering, filter_ref_associations, then
 *    expectation_maximization.
 *
 * Both branches converge on the same downstream calls (calculate_higher_rank_abundances,
 * calculate_relative_genomic_abundances, and the taxonomy::write_* output functions), so their
 * results are directly comparable.
 *
 * \param config Parsed CLI configuration (input/output file paths and all algorithm parameters).
 */
void tax_profile(taxor::profile::configuration& config)
{
    // <taxid, <taxid_string, taxname_string>>
    std::cout << "Parsing search results..." << std::flush;
    std::map<std::string, std::pair<std::string, std::string>> taxpath{};
    std::map<std::string, std::vector<taxonomy::Search_Result>> search_results = parse_search_results(config.search_file, taxpath);
    std::cout << "done" << std::endl;

    std::map<std::string, size_t> found_taxa{};
    std::map<std::string, std::vector<taxonomy::Search_Result>> profile_results{};
    std::map<std::string, double> tax_abundances{};

    if (config.classification_method == "skim")
    {
        // SKiM's classifier is a direct per-read statistical test with no
        // iterative refinement, so it intentionally skips Taxor's EM-specific
        // preprocessing (reference filtering, filter_ref_associations) and
        // operates straight on the parsed search results - keeping it
        // faithful to the published method is what makes a later
        // side-by-side comparison against the EM path meaningful.
        std::cout << "Classifying reads (SKiM)..." << std::flush;

        skim::parameters skim_params{};
        skim_params.kmer_size = config.skim_kmer_size;
        skim_params.n_fixed = config.skim_nfixed;
        skim_params.cutoff = pow(10.0, -static_cast<double>(config.skim_cutoff_exponent));

        std::tie(profile_results, found_taxa) = skim::classify_reads(search_results, skim_params, config.threads);
        std::cout << "done" << std::endl;

        std::cout << "Calculate sequence abundances..." << std::flush;
        // Mirrors what expectation_maximization does after its own loop: one
        // (non-iterative) call is a mathematically correct "sequence
        // abundance from a fixed classification" computation, since
        // update_log_prior_probabilities doesn't use its incoming log_priors
        // values in the computation, only taxa's keys and profile_results.
        std::map<std::string, double> log_priors = initialize_prior_log_probabilities(found_taxa);
        double const unclassified_abundance = update_log_prior_probabilities(log_priors, found_taxa, profile_results, config.threads);
        log_priors.insert(std::move(std::make_pair("unclassified", unclassified_abundance)));
        for (auto & t : log_priors)
            log_priors.at(t.first) = exp(t.second);
        tax_abundances = std::move(log_priors);
        std::cout << "done" << std::endl;
    }
    else
    {
        std::cout << "Remove matches to nonunique references..." << std::flush;

        // 1st round of reference filtering
        ankerl::unordered_dense::set<std::string> ref_unique_mappings = get_refs_with_uniquely_mapping_reads(search_results, config.threads);
        remove_matches_to_nonunique_refs(search_results, ref_unique_mappings, config.threads);
        ref_unique_mappings.clear();

        std::cout << "done" << std::endl;
        std::cout << "Remove low confidence references..." << std::flush;

        // 2nd round of reference filtering
        std::map<std::string,std::pair<uint64_t,uint64_t>> map_counts = count_unique_ambiguous_mappings_per_reference(search_results, config.threads);
        // at least 3 uniquely mapped reads & at least 10% of all mappings unique
        // TODO: may use different parameters for Illumina data
        remove_low_confidence_references(search_results, map_counts, 3, 0.01, config.threads);
        map_counts.clear();

        std::cout << "done" << std::endl;
        std::cout << "Filter associated references..." << std::flush;

        found_taxa = filter_ref_associations(search_results, config.threads);

        std::cout << "done" << std::endl;
        std::cout << "Start EM algorithm..." << std::flush;

        // returns nucleotide abundances
        tax_abundances = expectation_maximization(config.em_steps, found_taxa, search_results, profile_results, config.threads);

        std::cout << "done" << std::endl;
    }

    std::cout << "Calculate higher rank sequence abundances.." << std::flush;
    
    std::map<std::string, taxonomy::Profile_Output> rank_profiles = calculate_higher_rank_abundances(tax_abundances,taxpath);
    std::cout << "done" << std::endl;
    std::cout << "Write sequence abundances..." << std::flush;

    taxonomy::write_sequence_abundance_file(config.sequence_abundance_file, rank_profiles, config.sample_id, config.threshold);

    std::cout << "done" << std::endl;
    std::cout << "Calculate genomic abundances ..." << std::flush;

    calculate_relative_genomic_abundances(tax_abundances, found_taxa, profile_results, config.threads);

    
    std::cout << "done" << std::endl;
    std::cout << "Write remaining output files ..." << std::flush;


    rank_profiles.clear();
    rank_profiles = calculate_higher_rank_abundances(tax_abundances,taxpath);
    
    taxonomy::write_biobox_profiling_file(config.report_file, rank_profiles, config.sample_id, config.threshold);
    taxonomy::write_biobox_binning_file(config.binning_file, profile_results, config.sample_id);
    std::cout << "done" << std::endl;
}

/*!\brief Entry point for the `taxor profile` subcommand: parses CLI arguments and runs the pipeline.
 * \param parser The seqan3 argument parser for the `profile` subcommand (options registered by set_up_subparser_layout).
 * \return 0 on success; -1 if argument parsing/validation fails or tax_profile() throws.
 */
int execute(seqan3::argument_parser & parser)
{
    taxor::profile::configuration config;
    //chopper::layout::data_store data;

    set_up_subparser_layout(parser, config);

    try
    {
        parser.parse();

        if (config.classification_method == "skim" && config.skim_kmer_size == 0)
            throw seqan3::argument_parser_error{"--kmer-size is required (and must match the value used for "
                                                 "taxor build/search) when --classification-method skim is selected."};

    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[TAXOR PROFILE ERROR] " << ext.what() << '\n';
        return -1;
    }

    try
    {
        tax_profile(config);
    }
    catch (std::exception const & e)
    {
        std::cerr << "[TAXOR PROFILE ERROR] " << e.what() << '\n';
        return -1;
    }

    return 0;
}

}