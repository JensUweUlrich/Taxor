
#pragma once

#include <algorithm>
#include <cmath>

#include "threshold_parameters.hpp"
#include "kmer_model.hpp"
#include "fracminhash_model.hpp"
#include "syncmer_model.hpp"

/*!\file threshold.hpp
 * \brief Declares hixf::threshold::threshold, the dispatcher that picks which statistical model
 *        (percentage, k-mer-count, syncmer-match-ratio or FracMinHash) is used to turn a query's number of
 *        matching hashes into a pass/fail decision when searching the Hierarchical Interleaved XOR Filter.
 */
namespace hixf::threshold
{

/*!\brief Selects and evaluates the appropriate statistical threshold model for a search run.
 *
 * Constructed once per search run from a threshold_parameters instance, which chooses one of four
 * threshold_kinds to use for the whole run:
 *  - \c percentage: a fixed fraction of the query's hash count is required to match (user-supplied).
 *  - \c syncmer_model: uses the empirical syncmer-match-ratio lookup table (syncmer_model.hpp).
 *  - \c kmer_model: uses the k-mer mutation-count confidence-interval model (kmer_model.hpp); selected when
 *    every window yields exactly one k-mer (no minimiser subsampling) and FracMinHash is not requested.
 *  - \c fracminhash: uses the FracMinHash containment-index confidence-interval model
 *    (fracminhash_model.hpp); the fallback for windowed/minimised k-mer sets.
 *
 * get() is then called once per query/reference-bin pair (typically from multiple search worker threads) to
 * compute the minimum number of matching hashes required to call the query a hit.
 */
class threshold
{
public:
    threshold() = default;
    threshold(threshold const &) = default;
    threshold & operator=(threshold const &) = default;
    threshold(threshold &&) = default;
    threshold & operator=(threshold &&) = default;
    ~threshold() = default;

    /*!\brief Constructs a threshold evaluator, selecting the statistical model to use for this search run.
     * \param arguments Threshold configuration for this search run (kmer size, error rate, window size, and
     *                   which model to prefer). See threshold_parameters and the class-level documentation
     *                   above for the selection order.
     */
    threshold(threshold_parameters & arguments)
    {
        kmer_size = arguments.kmer_size;
        error_rate = arguments.seq_error_rate;
        size_t kmers_per_window = arguments.window_size - kmer_size + 1;
        //std::cout << arguments.percentage << std::endl << std::flush;
        if (arguments.percentage > 0.0 && arguments.percentage <= 1.0)
        {
            threshold_kind = threshold_kinds::percentage;
            threshold_percentage = arguments.percentage;
            std::cout << "use percentage-model\t" << arguments.percentage << std::endl << std::flush;
        }
        else if (arguments.use_syncmer)
        {
            threshold_kind = threshold_kinds::syncmer_model;
            std::cout << "use syncmer model" << std::endl << std::flush;
        }
        else if (kmers_per_window == 1 && !arguments.fracminhash)
        {
            std::cout << "use kmer-model" << std::endl << std::flush;
            threshold_kind = threshold_kinds::kmer_model;
        }
        else
        {
            threshold_kind = threshold_kinds::fracminhash;
            std::cout << "use frac minhash" << std::endl << std::flush;
        }
    }

    /*!\brief Computes the minimum number of matching hashes required to call a query a hit, using whichever
     *        statistical model was selected in the constructor.
     * \param minimiser_count Total number of query hashes (k-mers/minimisers/syncmers) considered.
     * \param scaling_factor  FracMinHash scaling factor; only used by the \c fracminhash model.
     * \return The match-count threshold: a query/reference-bin pair needs at least this many matching hashes
     *         to be classified as a hit.
     *
     * A small constant false-positive correction (`minimiser_count * 0.0039`) is subtracted for the
     * \c kmer_model and \c fracminhash cases to compensate for the expected XOR-filter false-positive rate.
     * The result is always clamped to be non-negative.
     */
    size_t get(size_t minimiser_count, double scaling_factor) noexcept
    {
        if (minimiser_count == 0)
            return 0;

        double const fp_correction = minimiser_count * 0.0039;
        //std::cout << error_rate << std::endl << std::flush;
        switch (threshold_kind)
        {
            case threshold_kinds::syncmer_model:
            {
                double syncmer_match_ratio = hixf::threshold::get_min_syncmer_match_ratio(kmer_size, error_rate);
                return static_cast<size_t>(minimiser_count * syncmer_match_ratio);
            }
            case threshold_kinds::kmer_model:
            {
                hixf::threshold::TInterval ci = hixf::threshold::calculate_nmut_kmer_CI(error_rate, (size_t) kmer_size, minimiser_count, 0.95);
                double const result = static_cast<double>(minimiser_count) - static_cast<double>(ci.second) - fp_correction;
                return static_cast<size_t>(std::max(0.0, result));
            }
            case threshold_kinds::fracminhash:
            {
                std::pair<double, double> cont_dist_ci = hixf::threshold::calculate_containment_index_CI(error_rate,
                                                                                                         kmer_size,
                                                                                                         minimiser_count,
                                                                                                         scaling_factor,
                                                                                                         0.95);
                double const clow = std::isfinite(cont_dist_ci.first) ? std::clamp(cont_dist_ci.first, 0.0, 1.0) : 0.0;
                double const result = clow * minimiser_count - fp_correction;
                return static_cast<size_t>(std::max(0.0, result));
            }
            default:
            {
                return static_cast<size_t>(minimiser_count * threshold_percentage);
            }
        }
    }

private:
    
    // TODO: add specific syncmer model based on empirical data
    //       minimum number of OCS found for each error rate and each kmer
    //       0.90 <= r <=0.99 and 16 <= k <= 32
    enum class threshold_kinds
    {
        fracminhash,
        percentage,
        kmer_model,
        syncmer_model,
    };

    threshold_kinds threshold_kind{threshold_kinds::percentage};
    uint8_t kmer_size;
    double threshold_percentage{};
    double error_rate{};

    

    
};

} 
