#pragma once

#include <string>

namespace taxor::taxonomy
{

struct Profile_Output
{
    std::string rank;
    std::string taxid;
    std::string taxid_string;
    std::string taxname_string;
    double percentage;

};

std::string format(float f, int digits) {
    std::ostringstream ss;
    ss.precision(digits);
    ss << f;
    return ss.str();
}

void write_biobox_profiling_file(std::string &output_file, 
                                 std::map<std::string, Profile_Output> &tax_rank_abundances, 
                                 std::string &sample_id)
                                 //std::string &taxonomy_id)
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
            if (rank.second.rank.compare(tr) == 0)
                fout << rank.second.taxid << "\t" << rank.second.rank << "\t" << rank.second.taxid_string << "\t"
                        << rank.second.taxname_string << "\t" << format(rank.second.percentage*100,6) << "\n";
        }
    }
    fout.close();
}

void write_sequence_abundance_file(std::string &output_file, 
                                   std::map<std::string, Profile_Output> &tax_rank_abundances, 
                                   std::string &sample_id)
                                   //std::string &taxonomy_id)
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
            if (rank.second.rank.compare(tr) == 0)
                fout << rank.second.taxid << "\t" << rank.second.rank << "\t" << rank.second.taxid_string << "\t"
                        << rank.second.taxname_string << "\t" << format(rank.second.percentage*100,6) << "\n";
        }
    }
    fout.close();
}

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
        if (read.second.size() == 0)
            fout << read.first << "\t-\n";
        else
            fout << read.first << "\t" << read.second.at(0).taxid << "\n";
    }


    fout.close();
}

}