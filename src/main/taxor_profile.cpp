

#include <math.h> 

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


std::map<std::string, std::vector<taxonomy::Search_Result>> parse_search_results(std::string const filepath,
                                                                    std::map<std::string, std::pair<std::string, std::string>> &taxpath)
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

            if (!taxpath.contains(res.taxid))
            {
                std::pair<std::string, std::string> taxpair = std::make_pair(line[9], line[8]);
                taxpath.insert(std::move(std::make_pair(res.taxid, std::move(taxpair))));
            }
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
                read_ref_liklihoods.insert(std::move(std::make_pair(res.taxid, likeL)));
            }
        }
        else
        {
            // if there is one unique mapping between this read and a reference
            if (pair.second.at(0).taxid.compare("-") != 0)
            {
                read_ref_liklihoods.insert(std::move(std::make_pair(pair.second.at(0).taxid, 0.0)));
            }
        }

        likelihoods.insert(std::move(std::make_pair(pair.first, read_ref_liklihoods)));
    }

    return std::move(likelihoods);
}

void update_log_prior_probabilities(std::map<std::string, double> &log_priors,
                                    std::map<std::string, size_t> & taxa,
                                    std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results)
{
    std::map<std::string, uint64_t> ref_nts{};
    for (auto & t : taxa)
        ref_nts.insert(std::move(std::make_pair(t.first, 0)));

    // for each taxon sum all lengths of matching reads
    for (auto &read : profile_results)
    {
        if (read.second.at(0).taxid.compare("-") == 0)
            continue;
        for (auto &ref : read.second)
        {
            ref_nts.at(ref.taxid) += ref.query_len;
        }
    }
    // calculate average depth of coverage for each taxon
    // sum of all matched read lengths divided by length of taxon reference sequence
    double sum_avg_cov{0.0};
    for (auto & t : ref_nts)
    {
        log_priors.at(t.first) = static_cast<double>(t.second) / static_cast<double>(taxa.at(t.first));
        sum_avg_cov += log_priors.at(t.first);
    }

    // calculate relative genomic abundance for each taxon
    // divide average coverage of each taxon by the sum of average coverage of all taxa
    for (auto &t : log_priors)
    {
        log_priors.at(t.first) = log(t.second + 0.000000000001) - log(sum_avg_cov);
    }

}

void calculate_higher_rank_abundances(std::map<std::string, size_t> & taxa,
                                      std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results,
                                      std::map<std::string, std::pair<std::string, std::string>> &taxpath)
{
    
}

std::map<std::string, double> expectation_maximization(size_t iterations, 
                              std::map<std::string, size_t> & taxa,
                              std::map<std::string, std::vector<taxonomy::Search_Result>> &search_results,
                              std::map<std::string, std::vector<taxonomy::Search_Result>> &profile_results)
{
    std::map<std::string, double> log_priors = initialize_prior_log_probabilities(taxa);
    std::map<std::string, std::map<std::string, double>> log_likelihoods = calculate_log_likelihoods(search_results);
    double cond_log_likelihood = -__DBL_MAX__;
    size_t iter_step = 0;
    while(iter_step < iterations)
    {
        double new_cond_log_likelihood = 0;
        profile_results.clear();
        for (auto & read : search_results)
        {
            double max_post = -__DBL_MAX__;
            std::vector<taxonomy::Search_Result> best_match{};
            for (auto &res : read.second)
            {   

                if (read.second.at(0).taxid.compare("-") == 0)
                {
                    best_match.emplace_back(res);
                    break;
                }
               
                double post = log_likelihoods.at(read.first).at(res.taxid) + log_priors.at(res.taxid);

                new_cond_log_likelihood += post;
                if (post >= max_post)
                {
                    if (post > max_post)
                    {
                        max_post = post;
                        best_match.clear();
                    }
                    best_match.emplace_back(res);
                }

            }
            profile_results.insert(std::move(std::make_pair(read.first, std::move(best_match))));
        }
        
        // update referencs abundances (priors)
        update_log_prior_probabilities(log_priors, taxa, profile_results);
        double diff = new_cond_log_likelihood - cond_log_likelihood;
        if (diff < abs(log(0.0001)))
            break;

        cond_log_likelihood = new_cond_log_likelihood;
        iter_step++;
    }
    std::cout << "Number of EM steps needed: " << iter_step << std::endl;

    for (auto & t : log_priors)
    {
        log_priors.at(t.first) = exp(t.second);
    }

    return std::move(log_priors);
}

void tax_profile(taxor::profile::configuration& config)
{
    // <taxid, <taxid_string, taxname_string>>
    std::map<std::string, std::pair<std::string, std::string>> taxpath{};
    std::map<std::string, std::vector<taxonomy::Search_Result>> search_results = parse_search_results(config.search_file, taxpath);
    ankerl::unordered_dense::set<std::string> ref_unique_mappings = get_refs_with_uniquely_mapping_reads(search_results);

    //std::cout << ref_unique_mappings.size() << std::endl;
    //remove_matches_to_nonunique_refs(search_results, ref_unique_mappings);
    std::map<std::string, size_t> found_taxa = filter_ref_associations(search_results);
    std::map<std::string, std::vector<taxonomy::Search_Result>> profile_results{};
    std::map<std::string, double> tax_abundances = expectation_maximization(1000, found_taxa, search_results, profile_results);
    
    for (auto & t: tax_abundances)
    {
        if (t.second >= 0.000001)
        std::cout << t.first << "\t" << t.second << std::endl;
    }

    size_t unclassified = 0;
    for (auto & read : profile_results)
    {
        if (read.second.at(0).taxid.compare("-") == 0)
        {
            unclassified++;
            std::cout << read.second.at(0).query_len << std::endl;
        }
    }
    std::cout << "Unclassified : " << unclassified << " / " << profile_results.size() << std::endl;
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