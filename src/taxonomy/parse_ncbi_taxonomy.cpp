
#include <iostream>

#include "parse_ncbi_taxonomy.hpp"
#include "taxutil.hpp"

/*!\file parse_ncbi_taxonomy.cpp
 * \brief Implements parsing of an NCBI/RefSeq-style taxonomy TSV file into
 *        taxor::taxonomy::Species records for `taxor build`.
 */

namespace taxor::taxonomy
{
    /*!\brief Read a tab-separated NCBI/RefSeq taxonomy file and build one Species per row.
     * \param filepath Path to the tab-separated taxonomy input file.
     * \return A vector of Species built from the input rows; rows without a usable
     *         reference file path are skipped (a warning is printed to stderr for each).
     *
     * Each line is expected to hold at least an accession id (column 0) and taxid
     * (column 1); the reference file/ftp path (column 2) is reduced to its final
     * path component and stored as Species::file_stem, and organism name / GTDB-style
     * lineage name and taxid strings (columns 3-5) are copied verbatim if present.
     */
    std::vector<Species> parse_refseq_taxonomy_file(std::string const filepath)
    {
        std::vector<std::vector<std::string> > tax_file_lines{};
	    read_tsv(filepath, tax_file_lines);
	    uint64_t species_counter = 0;
		std::vector<Species> org_list{};
		uint16_t count{0};
		for (std::vector<std::string> line : tax_file_lines)
		{
			
			Species sp{};
			sp.accession_id = line.at(0);
			sp.taxid = line.at(1);

			sp.organism_name = "";
			sp.taxnames_string = "";
			sp.taxid_string = "";

			if (line.size() > 3)
				sp.organism_name = line[3];
			if (line.size() > 4)
				sp.taxnames_string = line[4];
			if (line.size() > 5)
				sp.taxid_string = line[5];

			// Strip trailing path separators first: NCBI's assembly_summary.txt
			// ftp_path (and similar directory-style paths) ends in '/', e.g.
			// ".../GCF_000006685.1_ASM668v1/" - without stripping it,
			// find_last_of() finds that trailing separator and substr(found+1)
			// past it returns an empty string, even though the field clearly
			// does contain a usable path component.
			std::string path_field = line.at(2);
			while (!path_field.empty() && (path_field.back() == '/' || path_field.back() == '\\'))
				path_field.pop_back();

			std::size_t found = path_field.find_last_of("/\\");
			if (found != std::string::npos)
				sp.file_stem = path_field.substr(found+1);
			else
				sp.file_stem = path_field; // no separator at all: the whole field is the stem
			if (sp.file_stem.compare("") == 0 || sp.file_stem.compare(" ") == 0)
			{
				// No usable file path for this accession (e.g. a withdrawn/suppressed
				// NCBI assembly with ftp_path "na" in assembly_summary.txt) - there is
				// no genome to index for it anyway, so skip it instead of aborting the
				// whole build over one bad row.
				std::cerr << "[TAXOR BUILD WARNING] Skipping " << sp.accession_id
						  << ": no file path found in taxonomy input (column 3 has no '/' or '\\').\n";
				continue;
			}

			org_list.emplace_back(std::move(sp));
	 	}
		return std::move(org_list);
    }
    
} // namespace taxor::taxonomy


