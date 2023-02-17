



#include "taxor_profile_configuration.hpp"
#include "taxor_profile.hpp"
#include "search_results.hpp"
#include <taxutil.hpp>

#include <ankerl/unordered_dense.h>

namespace taxor::profile
{

void set_up_subparser_layout(seqan3::argument_parser & parser, taxor::profile::configuration & config)
{
    parser.info.version = "0.1.0";
    parser.info.author = "Jens-Uwe Ulrich";
    parser.info.email = "jens-uwe.ulrich@hpi.de";
    parser.info.short_description = "Taxonomic profiling of a sample by giving read matching results of Taxor search";

    parser.info.description.emplace_back("Taxonomic profiling of the given read set");

    parser.add_subsection("Main options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.search_file,
                      '\0', "search-file", "taxor search file containing results of read querying against the HIXF index",
 /*                     "Provide the prefix you used for the output prefix in the chopper count --output-prefix option. "
                      "If you have different means of estimating the k-mer counts of your input data, make sure that a "
                      "file [INPUT-PREFIX].count exists. It needs to be tab-separated and consist of two columns: "
                      "\"[filepath] [tab] [weight/count]\".",*/
                      seqan3::option_spec::required);

    parser.add_option(config.report_file, '\0', "output-file", "A file name for the resulting output.");

    parser.add_option(config.threads,
                      '\0', "threads",
                      "The number of threads to use.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(32)});


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


std::map<std::string, std::vector<taxonomy::Search_Result>> parse_search_results(std::string const filepath)
{
    std::vector<std::vector<std::string> > tax_file_lines{};
    taxonomy::read_tsv(filepath, tax_file_lines);
    uint64_t species_counter = 0;
	std::map<std::string, std::vector<taxonomy::Search_Result>> results{};

    size_t idx = 0;
	for (std::vector<std::string> line : tax_file_lines)
	{
        if (idx++ == 0) 
            continue;
        std::string read_id = line[0];
        
        if (line[0].find_first_of(" ") != std::string::npos)
            read_id = line[0].substr(0,line[0].find_first_of(" "));
        taxonomy::Search_Result res{};
        res.read_id = read_id;
        if (line[1].compare("-") == 0)
        {
            res.taxid = "-";
        }
        else
        {
            res.taxid = line[3];
            res.ref_len = std::stoull(line[4]);
            res.query_len = std::stoull(line[5]);
            res.query_hash_count = std::stoull(line[6]);
            res.query_hash_match = std::stoull(line[7]);
        }

		if (!results.contains(read_id))
        {
            results.insert(std::move(std::make_pair(read_id, std::vector<taxonomy::Search_Result>{})));
        }
        results.at(read_id).emplace_back(std::move(res));
	}

    return std::move(results);
}


ankerl::unordered_dense::set<std::string> get_refs_with_uniquely_mapping_reads(std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results)
{
    ankerl::unordered_dense::set<std::string> ref_unique_mappings{};
    for (auto & pair : search_results)
    {
        if (pair.second.size() == 1)
        {
            if (pair.second.at(0).taxid.compare("-") != 0)
            {
                ref_unique_mappings.insert(pair.second.at(0).taxid);
            }
        }
    }
    return std::move(ref_unique_mappings);
}

void remove_matches_to_nonunique_refs(std::map<std::string, std::vector<taxonomy::Search_Result>>& search_results,
                                      ankerl::unordered_dense::set<std::string>& ref_unique_mappings)
{
    std::vector<taxonomy::Search_Result>::iterator search_iterator;
    for (auto & pair : search_results)
    {
        if (pair.second.size() > 1)
        {   
            size_t count_before = pair.second.size();
            search_iterator = pair.second.begin();
            while (search_iterator != pair.second.end())
            {
                if (!ref_unique_mappings.contains((*search_iterator).taxid))
                    search_iterator = pair.second.erase(search_iterator);
                else
                    search_iterator++;
            }
            if (count_before != pair.second.size())
                std::cout << pair.first << "\t" << count_before << "\t" << pair.second.size() << std::endl;
        }
    }
}

/** 
 * Filtering out suspicious mappings using an approach similar to
 * the two-stage taxonomy assignment algorithm in MegaPath (Leung et al., 2020)
*/
std::map<std::string, size_t> filter_ref_associations(std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results)
{
    std::map<std::string, Ref_Map_Info> ref_associations{};
    // taxid, length
    std::map<std::string, size_t> taxa_lengths{};
    for (auto & pair : search_results)
    {
        if (pair.second.size() == 0) continue;
        if (pair.second.size() == 1)
        {
            // if there is one unique mapping between this read and a reference
            if (pair.second.at(0).taxid.compare("-") != 0)
            {
                if (!ref_associations.contains(pair.second.at(0).taxid))
                    ref_associations.insert(std::move(std::make_pair(pair.second.at(0).taxid, Ref_Map_Info{})));

                // increment the number of uniquely mapped reads by 1
                ref_associations.at(pair.second.at(0).taxid).unique_assign_reads += 1;
                // increment the number of all mapped reads by 1
                ref_associations.at(pair.second.at(0).taxid).all_assigned_reads += 1;

                if (!taxa_lengths.contains(pair.second.at(0).taxid))
                    taxa_lengths.insert(std::move(std::make_pair(pair.second.at(0).taxid, pair.second.at(0).ref_len)));
            }
        }
        else
        {
            // collect all references assigned to that read
            std::vector<std::string> taxids{};
            for (auto & res : pair.second)
            {
                if (!ref_associations.contains(res.taxid))
                    ref_associations.insert(std::move(std::make_pair(res.taxid, Ref_Map_Info{})));
                
                taxids.emplace_back(res.taxid);
                // increment the number of all mapped reads by 1
                ref_associations.at(res.taxid).all_assigned_reads += 1;

                if (!taxa_lengths.contains(res.taxid))
                    taxa_lengths.insert(std::move(std::make_pair(res.taxid, res.ref_len)));
            }

            // iterate over all references assigned with that read
            for (std::string taxid1 : taxids)
            {
                for (std::string taxid2 : taxids)
                {
                    if (taxid1.compare(taxid2) == 0) continue;

                    if (!ref_associations.at(taxid1).associated_species.contains(taxid2))
                        ref_associations.at(taxid1).associated_species.insert(std::move(std::make_pair(taxid2, 0)));

                    // increment the number of shared reads between ref1 and ref2
                    ref_associations.at(taxid1).associated_species.at(taxid2) += 1;
                }
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
            // if more than 95% of ref's mapped reads also map to current associated ref
            if (ref.second.all_assigned_reads - assoc_ref.second < static_cast<uint64_t>(0.05 * static_cast<double>(ref.second.all_assigned_reads)))
            {
                // if the number of uniquely mapped reads of ref is less then 5% the number of uniquely mapped reads of associated ref
                if (ref.second.unique_assign_reads < static_cast<uint64_t>(0.05 * static_cast<double>(ref_associations.at(assoc_ref.first).unique_assign_reads)))
                {
                    explained_refs.insert(std::move(std::make_pair(ref.first, assoc_ref.first)));
                }
            }
        }
    }

    // iterate over search results and filter results of ambiguous mappings
    // reassign unique mappings of refs explained by another ref
    for (auto & pair : search_results)
    {
        if (pair.second.size() == 0) continue;
        if (pair.second.size() == 1)
        {
            // if unique mapping is explained by another ref
            if (explained_refs.contains(pair.second.at(0).taxid))
            {
                pair.second.at(0).taxid = explained_refs.at(pair.second.at(0).taxid);
                pair.second.at(0).ref_len = taxa_lengths.at(pair.second.at(0).taxid);
            }
        }
        else
        {
            // collect all references assigned to that read
            std::set<std::string> taxids{};
            for (auto & res : pair.second)
                taxids.emplace(res.taxid);   
            
            std::vector<taxonomy::Search_Result>::iterator it = pair.second.begin();
            // iterate over search results
            while (it != pair.second.end())
            {
                // if taxid is explained by another reference
                if (explained_refs.contains((*it).taxid))
                {
                    // if there is another mapping of this read to the reference that explains current taxid
                    // remove this match
                    if (taxids.contains(explained_refs.at((*it).taxid)))
                    {
                        it = pair.second.erase(it);
                        continue;
                    }
                    // reassign otherwise
                    else
                    {
                        (*it).taxid = explained_refs.at(pair.second.at(0).taxid);
                        (*it).ref_len = taxa_lengths.at(pair.second.at(0).taxid);
                    }
                }
                it++;
            }


        }
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
    
    return taxa_lengths;
}

std::map<std::string, double> initialize_prior_probabilities(std::map<std::string, size_t>& taxa)
{
    std::map<std::string, double> priors {};
    for (auto & taxon : taxa)
    {
        priors.insert(std::move(std::make_pair(taxon.first, 1.0 / static_cast<double>(taxa.size()))));
        std::cout << priors.at(taxon.first) << "\t" << taxa.size() << std::endl;
    }
    return std::move(priors);
}

void expectation_maximization(size_t iterations, std::map<std::string, size_t> & taxa)
{
    std::map<std::string, double> priors = initialize_prior_probabilities(taxa);
}

void tax_profile(taxor::profile::configuration& config)
{
    
    std::map<std::string, std::vector<taxonomy::Search_Result>> search_results = parse_search_results(config.search_file);
    ankerl::unordered_dense::set<std::string> ref_unique_mappings = get_refs_with_uniquely_mapping_reads(search_results);

    std::cout << ref_unique_mappings.size() << std::endl;
    //remove_matches_to_nonunique_refs(search_results, ref_unique_mappings);
    std::map<std::string, size_t> found_taxa = filter_ref_associations(search_results);
    expectation_maximization(1000, found_taxa);
}

int execute(seqan3::argument_parser & parser)
{
    taxor::profile::configuration config;
    //chopper::layout::data_store data;

    set_up_subparser_layout(parser, config);

    try
    {
        parser.parse();

        // TODO: sanity check of parameters

    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[TAXOR PROFILE ERROR] " << ext.what() << '\n';
        return -1;
    }

    tax_profile(config);

    return 0;
}

}