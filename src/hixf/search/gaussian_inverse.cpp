
#include "gaussian_inverse.hpp"
#include <sstream>
#include <cmath>

/*!\file gaussian_inverse.cpp
 * \brief Implements a rational-function approximation of the inverse standard normal CDF (Abramowitz &
 *        Stegun formula 26.2.23), used to convert confidence levels into z-scores for the threshold models.
 */
namespace hixf::threshold
{

    /*!\brief Abramowitz-Stegun rational approximation used as a building block of the inverse normal CDF.
     *
     * Implements Abramowitz and Stegun formula 26.2.23, a well-known rational-function approximation with
     * an absolute error less than 4.5e-4. It is not the inverse normal CDF itself; NormalCDFInverse() calls
     * it with a suitably transformed argument depending on which tail of the distribution \c p falls in.
     * \param t Transformed input, sqrt(-2 * log(x)) for the appropriate tail probability x.
     * \return Approximate value of the standard normal quantile function for that tail.
     */
    double RationalApproximation(double t)
    {
        // Abramowitz and Stegun formula 26.2.23.
        // The absolute value of the error should be less than 4.5 e-4.
        double c[] = { 2.515517, 0.802853, 0.010328 };
        double d[] = { 1.432788, 0.189269, 0.001308 };

        return t - ((c[2] * t + c[1]) * t + c[0]) /
            (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
    }

    /*!\brief Approximates the value of the inverse normal cumulative distribution function (probit function).
     * \param p Probability, must be strictly between 0 and 1.
     * \return The z score z such that Phi(z) = p, where Phi is the standard normal CDF.
     * \throws std::invalid_argument if p is not in the open interval (0, 1).
     */
    double NormalCDFInverse(double p)
    {

        if (p <= 0.0 || p >= 1.0)
        {
            std::stringstream os;
            os << "Invalid input argument (" << p
                << "); must be larger than 0 but less than 1.";
            throw std::invalid_argument(os.str());
        }

        // See article above for explanation of this section.
        if (p < 0.5)
        {
            // F^-1(p) = - G^-1(p)
            return -RationalApproximation(sqrt(-2.0 * log(p)));
        }
        else
        {
            // F^-1(p) = G^-1(1-p)
            //std::cout << "RationalApproximation(sqrt(-2.0 * log(1.0 - p))) is: " << RationalApproximation(sqrt(-2.0 * log(1.0 - p))) << std::endl;
            return RationalApproximation(sqrt(-2.0 * log(1.0 - p)));

        }
    }

}