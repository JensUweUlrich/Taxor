
#pragma once

namespace hixf::threshold
{

    /**
    *   approximates the value of the inverse normal cumulative distribution function
    *   @p      : probability, has to be between 0 and 1
    *   @return : z score
    */
    double NormalCDFInverse(double p);

}