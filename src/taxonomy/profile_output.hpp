#pragma once

#include <string>

/*!\file profile_output.hpp
 * \brief Output-writing functions for `taxor profile`: turns per-taxon
 *        abundance data and per-read search results into CAMI-format
 *        (biobox) profiling and binning report files.
 */

namespace taxor::taxonomy
{

/*!\brief Aggregated abundance information for one taxon at one taxonomic rank,
 *        as accumulated by `taxor profile` from search results.
 */
struct Profile_Output
{
    //!\brief Taxonomic rank this entry belongs to (e.g. "species", "genus").
    std::string rank;
    //!\brief Taxid of this taxon.
    std::string taxid;
    //!\brief Full lineage of taxids from root to this taxon, semicolon-separated.
    std::string taxid_string;
    //!\brief Full lineage of taxon names from root to this taxon, semicolon-separated.
    std::string taxname_string;
    //!\brief Relative abundance of this taxon, as a fraction in [0, 1].
    double percentage;

};

/*!\brief Format a floating point value as a string with a fixed number of
 *        significant digits.
 * \param f The value to format.
 * \param digits Number of significant digits to keep (passed to
 *               std::ostream::precision()).
 * \return The formatted value as a string.
 */
std::string format(float f, int digits) {
    std::ostringstream ss;
    ss.precision(digits);
    ss << f;
    return ss.str();
}

/*!\brief Write a CAMI-format (biobox) taxonomic profiling report file.
 *
 * Emits one row per taxon in `tax_rank_abundances` whose relative abundance is
 * at or above `threshold`, grouped by rank in the fixed CAMI rank order
 * (superkingdom, phylum, class, order, family, genus, species). Any taxon
 * whose recorded rank is not one of these seven is silently omitted.
 *
 * \param output_file Path of the file to write; overwritten if it already exists.
 * \param tax_rank_abundances Per-taxon abundance entries, keyed by taxid (or
 *                            equivalent identifier used by the caller).
 * \param sample_id Sample identifier written to the report header.
 * \param threshold Minimum relative abundance (fraction in [0, 1]) a taxon must
 *                  have to be included in the output.
 */
void write_biobox_profiling_file(std::string &output_file,
                                 std::map<std::string, Profile_Output> &tax_rank_abundances,
                                 std::string &sample_id,
                                 double &threshold)
{
    std::vector<std::string> tax_ranks{"superkingdom","phylum", "class", "order","family", "genus", "species"};

    std::ofstream fout{output_file};
    fout << "@SampleID:" << sample_id << "\n";
    fout << "@Version:0.10.0\n";
    fout << "@Ranks:superkingdom|phylum|class|order|family|genus|species\n";
    //fout << "@TaxonomyID:ncbi_taxonomy\n";
    fout << "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n";

    for (std::string tr : tax_ranks)
    {
        for (auto & rank : tax_rank_abundances)
        {
            if (rank.second.rank.compare(tr) == 0 && rank.second.percentage >= threshold)
                fout << rank.second.taxid << "\t" << rank.second.rank << "\t" << rank.second.taxid_string << "\t"
                        << rank.second.taxname_string << "\t" << format(rank.second.percentage*100,6) << "\n";
        }
    }
    fout.close();
}

/*!\brief Write a CAMI-format (biobox) sequence abundance report file.
 *
 * Behaves like write_biobox_profiling_file(), with one addition: if
 * `tax_rank_abundances` contains an entry keyed "unclassified", an extra row
 * for it is written first (rank "no rank", taxpath/taxpathsn "-"), regardless
 * of `threshold`.
 *
 * \param output_file Path of the file to write; overwritten if it already exists.
 * \param tax_rank_abundances Per-taxon abundance entries, keyed by taxid (or
 *                            "unclassified" for the unclassified fraction).
 * \param sample_id Sample identifier written to the report header.
 * \param threshold Minimum relative abundance (fraction in [0, 1]) a classified
 *                  taxon must have to be included in the output.
 */
void write_sequence_abundance_file(std::string &output_file,
                                   std::map<std::string, Profile_Output> &tax_rank_abundances,
                                   std::string &sample_id,
                                   double &threshold)
{
    std::vector<std::string> tax_ranks{"superkingdom","phylum", "class", "order","family", "genus", "species"};

    std::ofstream fout{output_file};
    fout << "@SampleID:" << sample_id << "\n";
    fout << "@Version:0.10.0\n";
    fout << "@Ranks:superkingdom|phylum|class|order|family|genus|species\n";
    //fout << "@TaxonomyID:ncbi_taxonomy\n";
    fout << "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n";
    if (tax_rank_abundances.contains("unclassified"))
        fout << "unclassified\tno rank\t-\t-\t" << format(tax_rank_abundances.at("unclassified").percentage*100,6) << "\n";
    for (std::string tr : tax_ranks)
    {
        for (auto & rank : tax_rank_abundances)
        {
            //std::cout << tr << "\t" << rank.second.taxid << std::endl;
            if (rank.second.rank.compare(tr) == 0 && rank.second.percentage >= threshold)
                fout << rank.second.taxid << "\t" << rank.second.rank << "\t" << rank.second.taxid_string << "\t"
                        << rank.second.taxname_string << "\t" << format(rank.second.percentage*100,6) << "\n";
        }
    }
    fout.close();
}

/*!\brief Write a CAMI-format (biobox) binning report file, mapping each read to
 *        the taxid it was classified as.
 *
 * \param output_file Path of the file to write; overwritten if it already exists.
 * \param binning_results Per-read candidate match results, keyed by read id. A
 *                         read is written as unclassified ("-") if its result
 *                         vector is empty OR if its first entry is the "-"
 *                         accession_id sentinel (see Search_Result); both forms
 *                         occur depending on which caller populated the map, so
 *                         both are handled explicitly. Otherwise the tax_id of
 *                         the first (best) candidate is written.
 * \param sample_id Sample identifier written to the report header.
 */
void write_biobox_binning_file(std::string &output_file,
                               std::map<std::string, std::vector<Search_Result>> &binning_results,
                               std::string &sample_id)
{
    std::ofstream fout{output_file};
    fout << "@SampleID:" << sample_id << "\n";
    fout << "@Version:0.10.0\n";
    fout << "@@SEQUENCEID\tTAXID\n";

    for (auto & read : binning_results)
    {
        // Unclassified reads may be represented either as an empty vector, or
        // as a single placeholder entry with accession_id "-" (whose tax_id is
        // never populated by the search-result parser) - print "-" for both,
        // rather than falling through to an empty tax_id in the latter case.
        if (read.second.size() == 0 || read.second.at(0).accession_id.compare("-") == 0)
            fout << read.first << "\t-\n";
        else
            fout << read.first << "\t" << read.second.at(0).tax_id << "\n";
    }


    fout.close();
}

}