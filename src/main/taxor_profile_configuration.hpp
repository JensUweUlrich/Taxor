#pragma once

#include <string>

/*!\file taxor_profile_configuration.hpp
 * \brief The `taxor profile` CLI configuration struct; see set_up_subparser_layout() in
 *        taxor_profile.cpp for how each field is bound to a command-line option.
 */

namespace taxor::profile
{

//!\brief All configurable parameters for a `taxor profile` run.
struct configuration
{
    //!\brief Path to the `taxor search` results TSV file to profile (required).
    std::string search_file{};
    //!\brief Output path for the CAMI-format read-to-taxon binning file (required).
    std::string binning_file{};
    //!\brief Output path for the CAMI-format genomic-abundance report file (required).
    std::string report_file{};
    //!\brief Output path for the CAMI-format sequence-abundance file (includes an "unclassified" row).
    std::string sequence_abundance_file{};
    //!\brief Sample identifier written into the header of every output file.
    std::string sample_id{};
    //!\brief Minimum relative abundance for a taxon to be included in the report/sequence-abundance files.
    double threshold{0.001};
    //double error_rate{0.04};
    //!\brief Number of worker threads used to parallelize the read-level passes of both classification methods.
    uint8_t threads{1u};
    bool output_verbose_statistics{false};
    bool debug{false};
    //!\brief Maximum number of EM iterations (only used by --classification-method em).
    uint16_t em_steps{100};

    //!\brief Read classification method: "em" (default, Taxor's own EM profiler) or
    //!       "skim" (Schneggenburger & Zola's binomial-hypothesis-test classifier).
    std::string classification_method{"em"};
    //!\brief k-mer size used at index build time; required (and validated) only when
    //!       classification_method=="skim", since SKiM's model needs it to derive the
    //!       k-mer/syncmer universe size. 0 is a sentinel meaning "not provided".
    uint8_t skim_kmer_size{0u};
    //!\brief SKiM's fixed normalization sample size for the binomial test (see skim_classifier.hpp).
    uint16_t skim_nfixed{100u};
    //!\brief SKiM's p-value cutoff exponent e, i.e. a read is classified only if its best
    //!       candidate's p-value is below 10^-e.
    uint16_t skim_cutoff_exponent{12u};
};

}
