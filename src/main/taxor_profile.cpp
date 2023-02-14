



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

void tax_profile(taxor::profile::configuration& config)
{
    std::map<std::string, std::vector<taxonomy::Search_Result>> search_results = parse_search_results(config.search_file);
    ankerl::unordered_dense::set<std::string> ref_unique_mappings = get_refs_with_uniquely_mapping_reads(search_results);

    std::cout << ref_unique_mappings.size() << std::endl;
    remove_matches_to_nonunique_refs(search_results, ref_unique_mappings);
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