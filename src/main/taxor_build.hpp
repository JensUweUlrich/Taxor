#pragma once

#include <seqan3/argument_parser/all.hpp>

/*!\file taxor_build.hpp
 * \brief Public interface of the `taxor build` subcommand, declaring the
 * entry point invoked by main.cpp; the implementation lives in taxor_build.cpp.
 */

namespace taxor::build
{

    /*!\brief Entry point for the `taxor build` subcommand: parses the build-specific
     * CLI options into a taxor::build::configuration, then drives taxonomy parsing,
     * HIXF layout computation, and index construction.
     * \param parser The `build` sub-parser, already selected by main.cpp's top-level parser.
     * \return 0 on success, -1 on a CLI parsing/validation error, or the propagated
     * error from taxonomy parsing, layout creation, or index building.
     */
    int execute(seqan3::argument_parser & parser);

}