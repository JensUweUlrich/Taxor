
#include "parse_ncbi_taxonomy.hpp"
#include "taxutil.hpp"

namespace taxor::taxonomy
{
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
			sp.accession_id = line[0];
			sp.taxid = line[1];

			sp.organism_name = "";
			sp.taxnames_string = "";
			sp.taxid_string = "";

			if (line.size() > 3)
				sp.organism_name = line[3];
			if (line.size() > 4)
				sp.taxnames_string = line[4];
			if (line.size() > 5)
				sp.taxid_string = line[5];
			
			std::size_t found = line[2].find_last_of("/\\");
			if (found != std::string::npos)
				sp.file_stem = line[2].substr(found+1);
			if (sp.file_stem.compare("") == 0 || sp.file_stem.compare(" ") == 0)
				throw std::runtime_error{"No file name found for" + sp.accession_id + " !!!"};
			
			org_list.emplace_back(std::move(sp));
	 	}
		return std::move(org_list);
    }
    
} // namespace taxor::taxonomy


