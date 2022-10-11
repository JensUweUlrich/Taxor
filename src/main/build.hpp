

#ifndef build_hpp
#define build_hpp


#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>
#include "multi_interleaved_xor_filter.hpp"

#include "Semaphore.hpp"
#include "SafeMap.hpp"
#include "Species.hpp"
#include <syncmer.hpp>



struct Seqs {
	
	std::string       seqid;
	seqan3::dna5_vector seq;
	
};

typedef std::map<std::string, std::vector<Species>> TGenMap;

// has to be x % 100 = 0 
uint64_t max_bins_per_filter = 4300;
uint64_t breakpoint = 22000;

inline std::string get_seqid( std::string header )
{
    return header.substr( 0, header.find( ' ' ) );
}

/**
   * split underlying sequence based on stretches of N's
   * @seq      : reference sequence as a string
   * @seqlen   : length of the sequence
   * @return   : vector of subsequences, that results from splitting the original sequence
   */
 std::vector<std::string> cutOutNNNs(std::string& seq, uint64_t seqlen)
 {
	 std::vector<std::string> splittedStrings;
	 size_t start = 0;
	 size_t end = 0;

	 while ((start = seq.find_first_not_of("N", end)) != std::string::npos)
	 {
		 end = seq.find("N", start);
		 if (end > seqlen)
		 {
			 std::string s = seq.substr(start, seqlen - start - 1);
			 splittedStrings.push_back(s);
			 break;
		 }
		 std::string s = seq.substr(start, end - start);

		 splittedStrings.push_back(s);
	 }
	 return splittedStrings;
 }


 /**
	   read reference sequences from files and store them in reference sequence queue
	   @queue_refs:    thread safe queue that stores reference sequences to put in the ibf
	   @config:        settings used to build the IBF (reference file path, number of references, kmer size, fragment length)
	   @stats:         statistics (like runtime for specific tasks) for building the ibf
	   @throws:        FileParserException
   */

 void parse_ref_seqs(std::vector< Seqs >& queue_refs, const std::filesystem::path& reference_file)
 {
	 //seqan::SeqFileIn seqFileIn;

	 seqan3::sequence_file_input fin{reference_file};


	 for ( auto& [seq, header, qual] : fin)
	 {
		 const std::string seqid = get_seqid( header );

		 int counter = 1;
		 auto charstring = seq | seqan3::views::to_char;
		 std::string s = std::string(charstring.begin(), charstring.end());

	 	 // remove all Ns from the sequence
		 std::stringstream buf;
		 for (std::string subseq : cutOutNNNs(s, s.size()))
		 {
			 buf << subseq;
		 }
		 s =  buf.str();
		 auto r = s | seqan3::views::char_to<seqan3::dna5>;
		 seqan3::dna5_vector newseq(r.begin(), r.end());

		 queue_refs.emplace_back(Seqs{ seqid, newseq });
					
	 }
			 
 }

 std::vector<std::pair<uint64_t, uint64_t>> compute_bin_number_per_species_batch(const std::vector<Species>& species, uint64_t start_species, uint64_t end_species,
 											   const std::string& gtdb_root, int kmer_size, int sync_size, int t_syncmer)
 {
	uint64_t bin_nr = 0;
	std::vector<std::pair<uint64_t, uint64_t>> sp_index_bins;
	std::set<uint64_t> syncmers{};
	for (uint64_t index = start_species; index <species.size() && index < end_species; ++index)
	{
//			std::cout << sp.name << std::endl;
		auto inpath = std::filesystem::path(gtdb_root) / "gtdb_genomes_reps_r202";
		if (species[index].accession_id.substr(0, 3).compare("GCA") == 0)
		{
			inpath /= "GCA";
		}
		else
		{
			inpath /= "GCF";
		}
		std::string infilename = species[index].accession_id + "_genomic.fna.gz";
		inpath = inpath / species[index].accession_id.substr(4, 3) / species[index].accession_id.substr(7, 3) / species[index].accession_id.substr(10, 3) / infilename;
		std::vector<Seqs> ref_seqs{};
		parse_ref_seqs(ref_seqs, inpath);
					
		
		for (const auto & ref_s : ref_seqs)
		{
			std::vector<uint64_t> strobe_hashes = hashing::seq_to_syncmers(kmer_size, ref_s.seq, sync_size, t_syncmer);
			//randstrobes.insert(randstrobes.end(), strobe_hashes.begin(), strobe_hashes.end());
			//randstrobes.insert(randstrobes.end(),std::make_move_iterator(strobe_hashes.begin()), std::make_move_iterator(strobe_hashes.end()));
			syncmers.insert(std::make_move_iterator(strobe_hashes.begin()), std::make_move_iterator(strobe_hashes.end()));
		}
		
		
		sp_index_bins.push_back(std::move(std::make_pair(index, (syncmers.size() / 100000) + 1)));
		syncmers.clear();
	}
	return sp_index_bins;

 }

std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> compute_bin_numbers(std::vector<Species>& species, const std::string& gtdb_root, 
										  int kmer_size, int sync_size, int t_syncmer)
{
	std::vector<std::future<std::vector<std::pair<uint64_t, uint64_t>>>> future_store{};
	Semaphore max_threads(6);

	uint16_t species_batches = (species.size() / 100) + 1;

	for (uint16_t batch = 0; batch < species_batches; ++batch)
	{
//		std::vector<Species> sp_batch = std::vector(species.begin() + (batch*100), species.begin() + ((batch+1)*100));

		future_store.push_back(std::async(std::launch::async,
               [](const std::vector<Species>& species_vec, Semaphore& maxJobs, uint64_t start_species, uint64_t end_species,
 				const std::string& gtdb_root, int kmer_size, int sync_size, int t_syncmer)
				{
                	 std::scoped_lock w(maxJobs);
                 	return compute_bin_number_per_species_batch(std::ref(species_vec), start_species, end_species,
				  		std::ref(gtdb_root), kmer_size, sync_size, t_syncmer);
               	}, std::ref(species), std::ref(max_threads), batch*100, ((batch+1)*100),
				  std::ref(gtdb_root), kmer_size, sync_size, t_syncmer)
        );
		
	}

	std::vector<std::pair<uint64_t, uint64_t>> finished_batches{};
	while (finished_batches.size() < future_store.size())
	{
		for (auto &future : future_store) 
		{
			std::vector<std::pair<uint64_t, uint64_t>> tmp = future.get();
            finished_batches.insert(finished_batches.end(), tmp.begin(), tmp.end());
        }
	}

	std::cout << "All threads finished" << std::endl;

	uint64_t bin_nr = 0;
	
	uint64_t cum_bin_nr = 0;
	uint16_t f_index = 0;
	uint64_t first_species = 0;
	uint64_t species_counter = 0;
	// <bins in filter, start_species, end_species>
	std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> filter_bins;
	for(auto& result : finished_batches) 
	{
		
		if ((bin_nr + result.second) > max_bins_per_filter)
		{
			filter_bins.push_back(std::make_tuple(bin_nr, first_species, species_counter));
			std::cout << f_index++ << "\t" << bin_nr << "\t" << cum_bin_nr << std::endl;
			species[result.first].filter_index = f_index;
			species[result.first].first_bin = 0;
			bin_nr = result.second;
			species[result.first].last_bin = bin_nr - 1;
			first_species = species_counter + 1;
		}
		else
		{
			species[result.first].filter_index = f_index;
			species[result.first].first_bin = bin_nr;
			bin_nr+= result.second;
			species[result.first].last_bin = bin_nr - 1;
		}
		cum_bin_nr += result.second;
		species_counter++;
	}
	
	if (bin_nr > 0)
	{
		filter_bins.push_back(std::make_tuple(bin_nr,first_species, species_counter-1));
		std::cout << f_index++ << "\t" << bin_nr << "\t" << cum_bin_nr << std::endl;
	}
	double mbytes = (double) (cum_bin_nr * 100000 * 1.23 * 8) / (double) 8388608;
	std::cout << "Expected size of multi-IXF: " << mbytes << " MBytes" << std::endl;
	
	return filter_bins;
}


std::pair<uint16_t, bool> add_species_to_filter(std::vector<Species>& species_vec, const uint64_t start_species_index,
 						   const uint64_t end_species_index, const std::string& gtdb_root, int kmer_size,
						   int sync_size, int t_syncmer, std::mutex& filter_mutex, 
						   seqan3::interleaved_xor_filter<>& ixf,
						   uint16_t filter_index, bool* terminate)
{
	bool success = true;
	uint64_t bin = 0;
	filter_mutex.lock();
	std::cout << filter_index << "\t" << start_species_index << "\t" << end_species_index << std::endl;
	filter_mutex.unlock();
	std::set<uint64_t> syncmers{};
	for (uint64_t spec = start_species_index; spec < species_vec.size() && spec < end_species_index; ++spec)
	{
		if (filter_index == 0)
			std::cout << filter_index << "\t" << spec << "\t" << species_vec[spec].name <<std::endl;
		auto inpath = std::filesystem::path(gtdb_root) / "gtdb_genomes_reps_r202";
		if (species_vec[spec].accession_id.substr(0, 3).compare("GCA") == 0)
		{
			inpath /= "GCA";
		}
		else
		{
			inpath /= "GCF";
		}
		std::string infilename = species_vec[spec].accession_id + "_genomic.fna.gz";
		inpath = inpath / species_vec[spec].accession_id.substr(4, 3) / species_vec[spec].accession_id.substr(7, 3) / 
				species_vec[spec].accession_id.substr(10, 3) / infilename;

		std::vector<Seqs> ref_seqs{};
		parse_ref_seqs(ref_seqs, inpath);

		// hier würfeln wir nur die hashes wild durcheinander durch das einfügen in ein ordered set
		// anschließend teilen wir auch noch die Menge auf
		// wenn wir aufteilen müssen wir di order beibehalten	
		//std::vector<uint64_t> randstrobes{};
		for (const auto & ref_s : ref_seqs)
		{
			std::vector<uint64_t> strobe_hashes = hashing::seq_to_syncmers(kmer_size, ref_s.seq, sync_size, t_syncmer);
			syncmers.insert(std::make_move_iterator(strobe_hashes.begin()), std::make_move_iterator(strobe_hashes.end()));
		}

		std::vector<uint64_t> sync_vec(syncmers.begin(), syncmers.end());
		
		if (spec == 306)
			std::cout << syncmers.size() << std::endl;
		uint64_t first_species_bin = bin;
		for (int j = 0; j < (syncmers.size() / 100000) + 1 ; j++)
		{

			std::vector<uint64_t>::iterator first = sync_vec.begin() + (j*100000);
			std::vector<uint64_t>::iterator last = sync_vec.begin() + ((j+1)*100000);
			if (j == syncmers.size() / 100000)
				last = sync_vec.end();
			std::vector<uint64_t> subvec = std::vector<uint64_t>(first, last);
			
			success = ixf.add_bin_elements(bin++, subvec);
			if (filter_index == 0)
			{
				filter_mutex.lock();
				std::cout << bin << "/" << ixf.bin_count() << "\t" << (syncmers.size() / 100000) + 1 << "\t" << success << std::endl;
				filter_mutex.unlock();
			}
//			std::cout << species.name << "\t" << syncmers.size() << "\t" << bin-1 << "/" << j << "\t" << subvec.size() << "\t" << success << std::endl;
//			filter_mutex.unlock();
			if (!success)
				break;
		}
		syncmers.clear();
		sync_vec.clear();
		
		if (!success)
			break;

		filter_mutex.lock();
		//bin_meta.emplace_back(std::move(std::make_tuple(species.name, filter_index, first_species_bin, bin-1)));
		species_vec[spec].filter_index = filter_index;
		species_vec[spec].first_bin = first_species_bin;
		species_vec[spec].last_bin = bin-1;
		filter_mutex.unlock();

	}
	std::cout << filter_index << "\tBins added to IXF" << std::endl;
	return std::make_pair(filter_index, success);
}

multi_interleaved_xor_filter build_ixf_index(std::vector<Species>& species, const std::string& filter_file, const std::filesystem::path& bin_meta_file,
                     const std::string& gtdb_root, int kmer_size, int sync_size, int t_syncmer)
{	
	//<number of bins, first_species_index, last_species_index>
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> filter_bins = compute_bin_numbers(std::ref(species), gtdb_root,
																			kmer_size, sync_size, t_syncmer);
	std::cout << "bins computed" << std::endl;
	multi_interleaved_xor_filter mixf{};
//	return mixf;
	uint16_t f_index = 0;
	for (auto & fb : filter_bins)
	{
		mixf.add_filter(f_index, std::get<0>(fb), 100000);
		++f_index;
	}
	
	
	Semaphore max_threads(6);
	std::mutex filter_mutex{};
	bool* terminate = new bool(false);
	
	uint64_t species_counter = 0;
	std::vector<bool> finished_filter(filter_bins.size());
	std::fill(finished_filter.begin(), finished_filter.end(), false);
	// repeat until all IXFs have been successfully created
	
	std::queue<std::future<std::pair<uint16_t, bool>>> future_store{};
	// create one thread per filter to add the maximum number of species to the filter
	for (f_index = 0; f_index < filter_bins.size(); ++f_index)
	{
		// don't re-build filter if it already has been created successfully in former iteration round
		future_store.push(std::async(std::launch::async,
            	[](std::vector<Species>& species_vec, Semaphore& maxJobs, const uint64_t start_species_index,
 					   const uint64_t end_species_index, const std::string& gtdb_root, int kmer_size, int sync_size, int t_syncmer,
						   std::mutex& filter_mutex, seqan3::interleaved_xor_filter<>& ixf,
						   uint16_t filter_index, bool* terminate){
                 			std::scoped_lock w(maxJobs);
                 			return add_species_to_filter(std::ref(species_vec), start_species_index, end_species_index, gtdb_root, kmer_size, sync_size, t_syncmer,
						   					filter_mutex, ixf, filter_index, terminate);
               		}, std::ref(species), std::ref(max_threads), std::get<1>(filter_bins[f_index]), 
					   std::get<2>(filter_bins[f_index]), std::ref(gtdb_root), kmer_size, sync_size, t_syncmer, std::ref(filter_mutex), 
					   std::ref(mixf.get(f_index)), f_index, terminate)
        );
	}
		
	uint64_t filter_counter = 0;
	bool success = true;
	int reset_seed = 0;
	// check which IXFs have been created successfully and which not
	while(!future_store.empty())
	{
		auto & future = future_store.front();
		std::pair<uint16_t, bool> result = future.get();
		finished_filter[result.first] = result.second;
		if (result.second)
		{
			std::cout << ++filter_counter << std::endl;
		}
		else
		{
			mixf.get(result.first).clear();
			mixf.get(result.first).set_seed();
			++reset_seed;
			f_index = result.first;
			future_store.push(std::async(std::launch::async,
            	[](std::vector<Species>& species_vec, Semaphore& maxJobs, const uint64_t start_species_index,
 					   const uint64_t end_species_index, const std::string& gtdb_root, int kmer_size, int sync_size, int t_syncmer,
						   std::mutex& filter_mutex, seqan3::interleaved_xor_filter<>& ixf,
						   uint16_t filter_index, bool* terminate){
                 			std::scoped_lock w(maxJobs);
                 			return add_species_to_filter(std::ref(species_vec), start_species_index, end_species_index, gtdb_root, kmer_size, sync_size, t_syncmer,
						   					filter_mutex, ixf, filter_index, terminate);
               		}, std::ref(species), std::ref(max_threads), std::get<1>(filter_bins[f_index]), 
					   std::get<2>(filter_bins[f_index]), std::ref(gtdb_root), kmer_size, sync_size, t_syncmer, std::ref(filter_mutex), 
					   std::ref(mixf.get(f_index)), f_index, terminate)
        	);
				
		}
		future_store.pop();
	}
	
	std::cout << "Reset seed : " << reset_seed << std::endl;	
	mixf.species_vector = std::move(species);	

    std::ofstream os( filter_file, std::ios::binary );
    cereal::BinaryOutputArchive archive( os );
    archive( mixf ); 

	return std::move(mixf);
}

multi_interleaved_xor_filter load_multi_ixf_index(const std::string& filter_file)
{
	multi_interleaved_xor_filter mixf{};
	std::ifstream              is( filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );
    archive( mixf );
	return std::move(mixf);
}

#endif