
#pragma once

/*!\file gaussian_inverse.hpp
 * \brief Declares an approximation of the inverse standard normal cumulative distribution function (quantile
 *        function), used to turn a confidence level into a z-score for the statistical threshold models.
 */
namespace hixf::threshold
{

    /*!\brief Approximates the value of the inverse normal cumulative distribution function (probit function).
     * \param p Probability, must be strictly between 0 and 1.
     * \return The z score z such that Phi(z) = p, where Phi is the standard normal CDF.
     *
     * \throws std::invalid_argument if p is not in the open interval (0, 1).
     */
    double NormalCDFInverse(double p);

}