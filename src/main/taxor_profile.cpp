
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/range/to.hpp>
#include <ranges>
#include <math.h>
#include <algorithm> 

#include "taxor_profile_configuration.hpp"
#include "taxor_profile.hpp"
#include "search_results.hpp"
#include <taxutil.hpp>
#include <profile_output.hpp>

#include <ankerl/unordered_dense.h>

namespace taxor::profile
{

void set_up_subparser_layout(seqan3::argument_parser & parser, taxor::profile::configuration & config)
{
    parser.info.version = "0.1.3";
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


ankerl::unordered_dense::set<std::string> get_refs_with_uniquely_mapping_reads(std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results)
{
    ankerl::unordered_dense::set<std::string> ref_unique_mappings{};
    for (auto & pair : search_results)
    {
        if (pair.second.size() == 1)
        {
            if (pair.second.at(0).accession_id.compare("-") != 0)
            {
                ref_unique_mappings.insert(pair.second.at(0).accession_id);
            }
        }
    }
    return std::move(ref_unique_mappings);
}


/**
 *  Remove all ambiguous read-to-reference assignments where the reference has no unique mapping
*/
void remove_matches_to_nonunique_refs(std::map<std::string, std::vector<taxonomy::Search_Result>>& search_results,
                                      ankerl::unordered_dense::set<std::string>& ref_unique_mappings)
{
    std::vector<taxonomy::Search_Result>::iterator search_iterator;
    for (auto & pair : search_results)
    {
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
}


std::map<std::string,std::pair<uint64_t,uint64_t>> count_unique_ambiguous_mappings_per_reference(
                                std::map<std::string, std::vector<taxonomy::Search_Result>>& search_results)
{
    // <taxid, <unique,ambiguous>>
    std::map<std::string,std::pair<uint64_t,uint64_t>> map_counts{};
    for (auto & pair : search_results)
    {
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
                if (!map_counts.contains(res.accession_id))
                {
                    map_counts.insert(std::move(std::make_pair(res.accession_id, std::move(std::make_pair(0,0)))));
                }
                map_counts.at(res.accession_id).second += 1;
            }
        }
    }

    return std::move(map_counts);
}

void remove_low_confidence_references(std::map<std::string, std::vector<taxonomy::Search_Result>>& search_results,
                                      std::map<std::string,std::pair<uint64_t,uint64_t>>& map_counts,
                                      uint8_t min_unique_mappings,
                                      float min_fraction_unique)
{
    ankerl::unordered_dense::set<std::string> accepted_refs{};
    for (auto & ref : map_counts)
    {   
        if (ref.second.first >= min_unique_mappings &&
                static_cast<float>(ref.second.first) / static_cast<float>(ref.second.first + ref.second.second) >= min_fraction_unique)
            accepted_refs.insert(ref.first);
    }
    remove_matches_to_nonunique_refs(search_results, accepted_refs);
}


/** 
 * Filtering out suspicious mappings using an approach similar to
 * the two-stage taxonomy assignment algorithm in MegaPath (Leung et al., 2020)
*/
std::map<std::string, size_t> filter_ref_associations(std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results, uint8_t threads)
{
    
    // taxid, length
    std::map<std::string, size_t> taxa_lengths{};
    std::map<std::string, Ref_Map_Info> ref_associations{};

    for (auto & pair : search_results)
    {
        if (pair.second.size() == 0) continue;
        if (pair.second.size() == 1)
        {
            // if there is one unique mapping between this read and a reference
            if (pair.second.at(0).accession_id.compare("-") != 0)
            {
                if (!ref_associations.contains(pair.second.at(0).accession_id))
                    ref_associations.insert(std::move(std::make_pair(pair.second.at(0).accession_id, Ref_Map_Info{})));

                // increment the number of uniquely mapped reads by 1
                ref_associations.at(pair.second.at(0).accession_id).unique_assign_reads += 1;
                // increment the number of all mapped reads by 1
                ref_associations.at(pair.second.at(0).accession_id).all_assigned_reads += 1;

                if (!taxa_lengths.contains(pair.second.at(0).accession_id))
                    taxa_lengths.insert(std::move(std::make_pair(pair.second.at(0).accession_id, pair.second.at(0).ref_len)));
            }
        }
        else
        {
            // collect all references assigned to that read
            std::vector<std::string> acc_ids{};
            for (auto & res : pair.second)
            {
                if (!ref_associations.contains(res.accession_id))
                    ref_associations.insert(std::move(std::make_pair(res.accession_id, Ref_Map_Info{})));
                
                acc_ids.emplace_back(res.accession_id);
                // increment the number of all mapped reads by 1
                ref_associations.at(res.accession_id).all_assigned_reads += 1;

                if (!taxa_lengths.contains(res.accession_id))
                    taxa_lengths.insert(std::move(std::make_pair(res.accession_id, res.ref_len)));
            }

            // iterate over all references assigned with that read
            for (std::string acc_id1 : acc_ids)
            {
                for (std::string acc_id2 : acc_ids)
                {
                    if (acc_id1.compare(acc_id2) == 0) continue;

                    if (!ref_associations.at(acc_id1).associated_species.contains(acc_id2))
                        ref_associations.at(acc_id1).associated_species.insert(std::move(std::make_pair(acc_id2, 0)));

                    // increment the number of shared reads between ref1 and ref2
                    ref_associations.at(acc_id1).associated_species.at(acc_id2) += 1;
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
   
    bool found = true;
    
    while (found)
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

    }

    
    // iterate over search results and filter results of ambiguous mappings
    // reassign unique mappings of refs explained by another ref => TODO: this should not be done, keep unique mappings
    // in such cases, 
    for (auto & pair : search_results)
    {
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

std::map<std::string, double> initialize_prior_log_probabilities(std::map<std::string, size_t>& taxa)
{
    std::map<std::string, double> priors {};
    for (auto & taxon : taxa)
    {
        priors.insert(std::move(std::make_pair(taxon.first, log(1.0 / static_cast<double>(taxa.size())))));
    }
    return std::move(priors);
}

std::map<std::string, std::map<std::string, double>> calculate_log_likelihoods(std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results)
{
    std::map<std::string, std::map<std::string, double>> likelihoods{};
    for (auto & pair : search_results)
    {
        std::map<std::string, double> read_ref_liklihoods{};
        if (pair.second.size() == 0) continue;
        if (pair.second.size() > 1)
        {
            // calculate sum of matchcount ratios
            double sum_ratio{0.0};
            for (auto & res : pair.second)
            {
                sum_ratio += static_cast<double>(res.query_hash_match) / static_cast<double>(res.query_hash_count);
            }

            // calculate log likelihoods of the single matches
            for (auto & res : pair.second)
            {
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

    return std::move(likelihoods);
}

double update_log_prior_probabilities(std::map<std::string, double> &log_priors,
                                    std::map<std::string, size_t> & taxa,
                                    std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results)
{
    std::map<std::string, uint64_t> ref_nts{};
    for (auto & t : taxa)
        ref_nts.insert(std::move(std::make_pair(t.first, 0)));

    // for each taxon sum all lengths of matching reads
    size_t all_nts{0};
    size_t unclassified_nts{0};
    std::cout << "Sum nts of matching reads per taxon..." << std::endl << std::flush;
    for (auto &read : profile_results)
    {
        if (read.second.size() == 0) continue;
        all_nts += read.second.at(0).query_len;
        if (read.second.at(0).accession_id.compare("-") == 0)
        {
            unclassified_nts += read.second.at(0).query_len;
            continue;
        }
        for (auto &ref : read.second)
        {
            ref_nts.at(ref.accession_id) += ref.query_len;
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
        for (size_t index = 0; index < taxid_path.size();++index)
        {

            if (!rank_profiles.contains(taxid_path[index]))
            {
                taxonomy::Profile_Output profile{};
                profile.taxid = taxid_path[index];
                profile.taxid_string = taxid_path[0];
                profile.taxname_string = taxname_path[0].substr(3);
                for (size_t index2 = 1; index2 <= index; ++index2)
                {
                    profile.taxid_string += "|";
                    profile.taxid_string += taxid_path[index2];
                    profile.taxname_string += "|";
                    profile.taxname_string += taxname_path[index2].substr(3);
                }

                profile.percentage = 0.0;

                if (taxname_path[index].substr(0, 1).compare("s") == 0)
				    profile.rank = "species";
			    else if(taxname_path[index].substr(0, 1).compare("g") == 0)
					profile.rank = "genus";
				else if(taxname_path[index].substr(0, 1).compare("f") == 0)
					profile.rank = "family";
				else if(taxname_path[index].substr(0, 1).compare("o") == 0)
					profile.rank = "order";
				else if(taxname_path[index].substr(0, 1).compare("c") == 0)
					profile.rank = "class";
				else if(taxname_path[index].substr(0, 1).compare("p") == 0)
					profile.rank = "phylum";
				else if(taxname_path[index].substr(0, 1).compare("k") == 0)
					profile.rank = "superkingdom";


                rank_profiles.insert(std::move(std::pair(taxid_path[index], std::move(profile))));
            }

            rank_profiles.at(taxid_path[index]).percentage += species_abundances.at(sp.first);
        }
    }

    return std::move(rank_profiles);
    
}

std::map<std::string, double> expectation_maximization(size_t iterations, 
                              std::map<std::string, size_t> & taxa,
                              std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results,
                              std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results)
{
    std::cout << "started" << std::endl;
    std::cout << "Initialize prior probabilities ..." << std::flush;
    std::map<std::string, double> log_priors = initialize_prior_log_probabilities(taxa);
    std::cout << "done" << std::endl;
    double cond_log_likelihood = -__DBL_MAX__;
    size_t iter_step = 0;
    double unclassified_abundance{0.0};
    while(iter_step < iterations)
    {
        std::cout << "Starting EM iteration " << iter_step << std::endl << std::flush;
        std::cout << "Calculate Log Likelihoods ..." << std::flush;
        std::map<std::string, std::map<std::string, double>> log_likelihoods = calculate_log_likelihoods(search_results);
        std::cout << "done" << std::endl;
        double new_cond_log_likelihood = 0;
        profile_results.clear();
        for (auto & read : search_results)
        {
            if (read.second.size() == 0) continue;
            double max_post = -__DBL_MAX__;
            double min_post = __DBL_MAX__;
            std::vector<taxonomy::Search_Result> best_match{};
            std::vector<taxonomy::Search_Result>::iterator worst_match{};
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
                        worst_match = it;
                        it++;
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

                new_cond_log_likelihood += post;
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
            profile_results.insert(std::move(std::make_pair(read.first, std::move(best_match))));
            // remove reference match with worst posterior probability in each iteration until convergence
            if (read.second.size() > 1)
                read.second.erase(worst_match);
        }
        // update referencs abundances (priors)
        std::cout << "Update prior probabilities ..." << std::flush;
        unclassified_abundance = update_log_prior_probabilities(log_priors, taxa, profile_results);
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

void calculate_relative_genomic_abundances(std::map<std::string, double> &log_priors,
                                           std::map<std::string, size_t> & taxa,
                                           std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results)
{
    log_priors.clear();
    std::map<std::string, uint64_t> ref_nts{};
    for (auto & t : taxa)
    {
        ref_nts.insert(std::move(std::make_pair(t.first, 0)));
        log_priors.insert(std::move(std::make_pair(t.first, 0.0)));
    }

    // for each taxon sum all lengths of matching reads
    for (auto &read : profile_results)
    {
        if (read.second.size() == 0) continue;
        if (read.second.at(0).accession_id.compare("-") == 0)
            continue;
        for (auto &ref : read.second)
        {
            if (ref_nts.contains(ref.accession_id))
                ref_nts.at(ref.accession_id) += ref.query_len;
        }
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

void tax_profile(taxor::profile::configuration& config)
{
    // <taxid, <taxid_string, taxname_string>>
    std::cout << "Parsing search results..." << std::flush;
    std::map<std::string, std::pair<std::string, std::string>> taxpath{};
    std::map<std::string, std::vector<taxonomy::Search_Result>> search_results = parse_search_results(config.search_file, taxpath);

    std::cout << "done" << std::endl;
    std::cout << "Remove matches to nonunique references..." << std::flush;

    // 1st round of reference filtering
    ankerl::unordered_dense::set<std::string> ref_unique_mappings = get_refs_with_uniquely_mapping_reads(search_results);
    remove_matches_to_nonunique_refs(search_results, ref_unique_mappings);
    ref_unique_mappings.clear();

    std::cout << "done" << std::endl;
    std::cout << "Remove low confidence references..." << std::flush;

    // 2nd round of reference filtering
    std::map<std::string,std::pair<uint64_t,uint64_t>> map_counts = count_unique_ambiguous_mappings_per_reference(search_results);
    // at least 3 uniquely mapped reads & at least 10% of all mappings unique
    // TODO: may use different parameters for Illumina data
    remove_low_confidence_references(search_results, map_counts, 3, 0.01);
    map_counts.clear();
   
    std::cout << "done" << std::endl;
    std::cout << "Filter associated references..." << std::flush;

    std::map<std::string, size_t> found_taxa = filter_ref_associations(search_results, config.threads);
   
    std::cout << "done" << std::endl;
    std::cout << "Start EM algorithm..." << std::flush;

    std::map<std::string, std::vector<taxonomy::Search_Result>> profile_results{};
    // returns nucleotide abundances
    std::map<std::string, double> tax_abundances = expectation_maximization(10, found_taxa, search_results, profile_results);

    std::cout << "done" << std::endl;
    std::cout << "Calculate higher rank sequence abundances.." << std::flush;
    
    std::map<std::string, taxonomy::Profile_Output> rank_profiles = calculate_higher_rank_abundances(tax_abundances,taxpath);
    std::cout << "done" << std::endl;
    std::cout << "Write sequence abundances..." << std::flush;

    taxonomy::write_sequence_abundance_file(config.sequence_abundance_file, rank_profiles, config.sample_id, config.threshold);

    std::cout << "done" << std::endl;
    std::cout << "Calculate genomic abundances ..." << std::flush;

    calculate_relative_genomic_abundances(tax_abundances, found_taxa, profile_results);

    
    std::cout << "done" << std::endl;
    std::cout << "Write remaining output files ..." << std::flush;


    rank_profiles.clear();
    rank_profiles = calculate_higher_rank_abundances(tax_abundances,taxpath);
    
    taxonomy::write_biobox_profiling_file(config.report_file, rank_profiles, config.sample_id, config.threshold);
    taxonomy::write_biobox_binning_file(config.binning_file, profile_results, config.sample_id);
    std::cout << "done" << std::endl;
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