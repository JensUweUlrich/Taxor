
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

		for (std::vector<std::string> line : tax_file_lines)
		{
			Species sp{};
			sp.accession_id = line[0];
			sp.taxid = line[1];
			sp.organism_name = line[3];
			
			std::string tmp;
			std::stringstream taxonomy(line[4]);
			while (getline(taxonomy, tmp, '\t'))
			{
				if (tmp.substr(0, 1).compare("s") == 0)
					sp.species = tmp.substr(3);
				else if(tmp.substr(0, 1).compare("g") == 0)
					sp.genus = tmp.substr(3);
				else if(tmp.substr(0, 1).compare("f") == 0)
					sp.family = tmp.substr(3);
				else if(tmp.substr(0, 1).compare("o") == 0)
					sp.order = tmp.substr(3);
				else if(tmp.substr(0, 1).compare("c") == 0)
					sp.class_tax = tmp.substr(3);
				else if(tmp.substr(0, 1).compare("p") == 0)
					sp.phylum = tmp.substr(3);
				else if(tmp.substr(0, 1).compare("k") == 0)
					sp.kingdom = tmp.substr(3);
			}
			std::size_t found = line[2].find_last_of("/\\");
			if (found != std::string::npos)
				sp.file_stem = line[2].substr(found+1);
			
			org_list.emplace_back(std::move(sp));
	 	}
		
		return std::move(org_list);
    }
    
} // namespace taxor::taxonomy


