#include <algorithm>
#include <cmath>

#include "parallel_util.hpp"
#include "skim_classifier.hpp"

namespace taxor::profile::skim
{

namespace
{

// Builds a table of P(Bin(n_fixed, p) >= x) for x in [0, n_fixed], via a single
// downward recurrence: P(X>=n_fixed) = PMF(n_fixed), then
// P(X>=x) = P(X>=x+1) + PMF(x) for x = n_fixed-1 .. 0. Every summand is
// non-negative, so there is no cancellation to worry about, and PMF is
// computed in log-space (lgamma-based log binomial coefficient) so this stays
// accurate even for the very small p / extreme-tail regime this is used in
// (p_i is often ~1e-6 to 1e-11, and cutoffs go down to 1e-12 or smaller).
std::vector<double> build_survival_table(double p, uint16_t n_fixed)
{
    std::vector<double> table(static_cast<size_t>(n_fixed) + 1, 0.0);

    if (p <= 0.0)
    {
        table[0] = 1.0; // P(X>=0) is always 1; P(X>=x>0) is 0 when p==0
        return table;
    }
    if (p >= 1.0)
    {
        std::fill(table.begin(), table.end(), 1.0);
        return table;
    }

    double const log_p = std::log(p);
    double const log_1mp = std::log(1.0 - p);
    double const log_n_fact = std::lgamma(static_cast<double>(n_fixed) + 1.0);

    double running_sum = 0.0;
    for (int x = static_cast<int>(n_fixed); x >= 0; --x)
    {
        double const log_choose = log_n_fact - std::lgamma(static_cast<double>(x) + 1.0)
                                              - std::lgamma(static_cast<double>(n_fixed - x) + 1.0);
        double const log_pmf = log_choose + static_cast<double>(x) * log_p
                                           + static_cast<double>(n_fixed - x) * log_1mp;
        running_sum += std::exp(log_pmf);
        table[static_cast<size_t>(x)] = running_sum;
    }

    return table;
}

} // namespace

double binomial_survival(uint32_t x, uint32_t n, double p)
{
    if (x == 0)
        return 1.0;
    if (x > n)
        return 0.0;

    std::vector<double> const table = build_survival_table(p, static_cast<uint16_t>(n));
    return table[x];
}

std::pair<std::map<std::string, std::vector<taxonomy::Search_Result>>, std::map<std::string, size_t>>
classify_reads(std::map<std::string, std::vector<taxonomy::Search_Result>> & search_results,
               parameters const & params,
               uint8_t threads)
{
    // Universe size: canonical k-mer space for the given kmer_size. The
    // syncmer-selection density and any FracMinHash-style scaling factor are
    // common to both this and every reference's actually-indexed k-mer count,
    // so they cancel in the ratio p_i = |K_i|/|U| and don't need to be
    // modeled separately here. This approximation assumes roughly uniform
    // syncmer/k-mer density and non-repetitive references; repeat-rich
    // genomes get a systematically conservative (not liberal) bias as a
    // result (their true |K_i| is smaller than ref_len implies, so p_i is
    // overestimated, making classification to them slightly more
    // conservative, not less).
    double const universe_size = std::pow(4.0, static_cast<double>(params.kmer_size)) / 2.0;

    // Precompute a p_i and survival table per distinct reference, sequentially,
    // before any parallel work starts, so the per-read classification loop
    // below is a read-only lookup into an already-built cache - safe to run
    // concurrently with no locking.
    std::map<std::string, std::vector<double>> survival_tables{};
    for (auto const & pair : search_results)
    {
        for (auto const & res : pair.second)
        {
            if (res.accession_id.compare("-") == 0)
                continue;
            if (survival_tables.contains(res.accession_id))
                continue;

            double p_i = static_cast<double>(res.ref_len) / universe_size;
            p_i = std::clamp(p_i, 0.0, 1.0);
            survival_tables.emplace(res.accession_id, build_survival_table(p_i, params.n_fixed));
        }
    }

    auto entries = taxor::util::to_pointer_vector(search_results);
    std::vector<std::map<std::string, std::vector<taxonomy::Search_Result>>> partial_profiles(std::max<size_t>(1, threads));

    auto worker = [&](size_t thread_index, size_t start, size_t end)
    {
        auto & profile_local = partial_profiles[thread_index];
        for (size_t idx = start; idx < end; ++idx)
        {
            auto & pair = *entries[idx];
            taxonomy::Search_Result const * best = nullptr;
            double best_pvalue = 1.0;

            for (auto const & res : pair.second)
            {
                if (res.accession_id.compare("-") == 0)
                    continue;
                if (res.query_hash_count == 0)
                    continue;

                auto const & table = survival_tables.at(res.accession_id);

                uint64_t x_bar;
                if (res.query_hash_count > static_cast<uint64_t>(params.n_fixed))
                    x_bar = (res.query_hash_match * static_cast<uint64_t>(params.n_fixed)) / res.query_hash_count;
                else
                    x_bar = res.query_hash_match;
                x_bar = std::min<uint64_t>(x_bar, params.n_fixed);

                double const pvalue = table[x_bar];

                // lowest p-value wins; ties broken deterministically by lowest
                // accession id (the paper doesn't specify a tie-break rule)
                if (best == nullptr || pvalue < best_pvalue ||
                    (pvalue == best_pvalue && res.accession_id < best->accession_id))
                {
                    best = &res;
                    best_pvalue = pvalue;
                }
            }

            if (best != nullptr && best_pvalue < params.cutoff)
            {
                profile_local.emplace(pair.first, std::vector<taxonomy::Search_Result>{*best});
            }
            else
            {
                // Represent "unclassified" the same way parse_search_results/EM
                // do (a single accession_id=="-" placeholder carrying the
                // read's query_len), not an empty vector: the shared
                // update_log_prior_probabilities/calculate_relative_genomic_abundances
                // functions skip size()==0 entries entirely, which would
                // silently drop unclassified reads' query_len from the
                // abundance denominator. write_biobox_binning_file handles
                // this "-" placeholder form correctly too.
                uint64_t const query_len = pair.second.empty() ? 0 : pair.second.front().query_len;
                taxonomy::Search_Result unclassified{};
                unclassified.read_id = pair.first;
                unclassified.accession_id = "-";
                unclassified.query_len = query_len;
                profile_local.emplace(pair.first, std::vector<taxonomy::Search_Result>{std::move(unclassified)});
            }
        }
    };

    taxor::util::parallel_for_indexed(worker, entries.size(), threads);

    std::map<std::string, std::vector<taxonomy::Search_Result>> profile_results{};
    for (auto & partial : partial_profiles)
        profile_results.merge(partial);

    std::map<std::string, size_t> taxa{};
    for (auto & pair : profile_results)
        for (auto & res : pair.second)
            if (res.accession_id.compare("-") != 0)
                taxa.emplace(res.accession_id, res.ref_len);

    return {std::move(profile_results), std::move(taxa)};
}

} // namespace taxor::profile::skim
