
#pragma once

#include <cstdint>

namespace hixf
{

/*!\file strong_types.hpp
 * \brief Provides small strong (tagged) types used to disambiguate otherwise plain-integer arguments.
 */

//!\brief Strong type for passing the window size.
struct window
{
    uint32_t v{}; //!< The wrapped window size value.
};

}
