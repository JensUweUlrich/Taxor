#pragma once

#include <seqan3/argument_parser/all.hpp>

/*!\file taxor_search.hpp
 * \brief Entry point declaration for the `taxor search` subcommand.
 */

namespace taxor::search
{

    /*!\brief Entry point for the `taxor search` subcommand.
     * \param parser The seqan3::argument_parser instance used to parse the subcommand's options.
     * \return 0 on success, -1 if argument parsing/sanity checks fail or the search itself throws.
     */
    int execute(seqan3::argument_parser & parser);

}