#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "search_results.hpp"

/*!\file skim_classifier.hpp
 * \brief The SKiM read classifier: a per-read binomial-hypothesis-test alternative to Taxor's
 *        own EM profiler, selectable via `taxor profile --classification-method skim`.
 *
 * Implements the classification algorithm from Schneggenburger & Zola, "SKiM: accurately
 * classifying metagenomic ONT reads in limited memory"
 * (https://pmc.ncbi.nlm.nih.gov/articles/PMC12502918/). See README.md's "SKiM classification
 * method" section for user-facing documentation and the approximation this implementation makes
 * to compute p_i from Taxor's available data.
 */

namespace taxor::profile::skim
{

//!\brief Tunable parameters for classify_reads(), bound to CLI options in taxor_profile.cpp.
struct parameters
{
    //!\brief k-mer size used at index build time (0 is an invalid sentinel; validated by the caller).
    uint8_t kmer_size{0};
    //!\brief Fixed normalization sample size for the binomial test (the paper's default is 100).
    uint16_t n_fixed{100};
    //!\brief p-value cutoff below which a read is classified to its best candidate (the paper's default is 1e-12).
    double cutoff{1e-12};
};

/*!\brief Upper-tail probability P(Bin(n, p) >= x).
 *
 * Computed via direct log-space summation of the binomial PMF (see skim_classifier.cpp for
 * numerical notes on why this, rather than the regularized incomplete beta function, is used).
 * Exposed standalone (rather than only as an internal detail of classify_reads) for testing
 * against reference values, e.g. scipy.stats.binom.sf(x-1, n, p).
 *
 * \param x Number of successes to test for (inclusive lower bound of the tail).
 * \param n Number of binomial trials.
 * \param p Per-trial success probability, in [0, 1].
 * \return P(X >= x) for X ~ Binomial(n, p); 1.0 if x==0, 0.0 if x>n.
 */
double binomial_survival(uint32_t x, uint32_t n, double p);

/*!\brief Classifies every read in `search_results` using the SKiM binomial-hypothesis-test method.
 *
 * For each read, and each of its (non-"-") candidate references: computes p_i (the reference's
 * modeled probability of a random k-mer match, from its reference length and `params.kmer_size`
 * - see skim_classifier.cpp), normalizes the read's match count against `params.n_fixed`, and
 * looks up P(Bin(n_fixed, p_i) >= normalized_match_count). The candidate with the lowest p-value
 * wins (ties broken by lowest accession id) if that p-value clears `params.cutoff`; otherwise the
 * read is unclassified. Unlike Taxor's EM path, this is a direct, single-pass per-read
 * statistical test: no iterative refinement, and no preprocessing shared across reads (it
 * operates directly on parse_search_results's output, not on EM's filtered/merged references).
 *
 * \param search_results Parsed search results. Not mutated (unlike the EM path, SKiM never
 *        prunes candidates from a read's list); taken by non-const reference only so its entries
 *        can be partitioned by pointer for parallel dispatch (see parallel_util.hpp).
 * \param params Classification parameters (k-mer size, n_fixed, cutoff).
 * \param threads Number of worker threads to partition the read scan across.
 * \return A (profile_results, taxa) pair, in exactly the shapes Taxor's existing abundance/report
 *         functions expect:
 *          - profile_results: read id -> a single winning Search_Result if classified (ties
 *            broken by lowest accession id), or a single "-" placeholder entry (matching
 *            parse_search_results's own convention for "no match") if unclassified. Every read
 *            present in `search_results` is represented, since downstream abundance calculations
 *            sum denominators over all reads, including unclassified ones.
 *          - taxa: accession_id -> ref_len, for every accession that won at least one read. Built
 *            as a complete pass before returning, since update_log_prior_probabilities does an
 *            unguarded lookup into a map built this way.
 */
std::pair<std::map<std::string, std::vector<taxonomy::Search_Result>>, std::map<std::string, size_t>>
classify_reads(std::map<std::string, std::vector<taxonomy::Search_Result>> & search_results,
               parameters const & params,
               uint8_t threads);

} // namespace taxor::profile::skim
