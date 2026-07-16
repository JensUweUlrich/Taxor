
/*!\file query.hpp
 * \brief Draft/prototype code for a per-read query classification routine (query_read()).
 *
 * The entire remainder of this file, including the `#pragma once` guard and the query_read()
 * function it contains, is wrapped in a single C-style block comment (`/* ... *\/`) and is
 * therefore inert: none of it is compiled or reachable.
 *
 * Usage check performed while documenting this file: grepping the whole src/ tree for both
 * "query.hpp" (as an #include target) and "query_read" (the only symbol this file would define)
 * found no matches outside this file itself. No other file includes query.hpp, and no call sites
 * for query_read exist elsewhere in the codebase as of this writing. The commented-out body
 * sketches a per-read classification loop (hashing a read into strobemer/syncmer hashes, bulk
 * counting them against a hierarchical filter's user bins, and classifying a read based on the
 * best-matching bin's hit ratio) that appears to be an early prototype superseded by the
 * classification logic now implemented in taxor_search.cpp / skim_classifier.hpp.
 */

/*#pragma once

void query_read()
{
    read_hash_vector rhv = seq_to_randstrobes_read(kmer_size, w_min, w_max, seq, sync_size, t_syncmer, q_p, max_dist);
	uint64_t strobe_nr = rhv[0].first.size();
//		TInterval ci = calculateCI(0.1, kmer_size, strobe_nr, 0.95);
//		uint64_t threshold = strobe_nr - ci.second;
		std::vector<uint64_t> max_found(mixf.species_vector.size());
		std::fill(max_found.begin(), max_found.end(), 0);
		for (const auto & p : rhv)
		{
			std::vector<TIXFAgent::counting_vector> count_vectors = mixf.bulk_count(p.first);

			std::vector<uint64_t> result(mixf.species_vector.size());
			for (uint64_t i = 0; i < result.size(); ++i)
			{
				result[i] = 0;
				if (mixf.species_vector[i].filter_index == UINT16_MAX)
					continue;
				for (uint64_t spec_bin = mixf.species_vector[i].first_bin; spec_bin <= mixf.species_vector[i].last_bin; ++spec_bin)
					result[i] += count_vectors[mixf.species_vector[i].filter_index][spec_bin];
			}

			for (int i = 0; i < max_found.size(); ++i)
			{
				if (result[i] > max_found[i])
					max_found[i] = result[i];
			}
		}
		double max_ratio = 0.0;
		uint64_t max_index = UINT64_MAX;
		std::vector<std::pair<uint64_t, double>> potential_indexes{};
		for (int i = 0; i < max_found.size(); ++i)
		{
			double ratio = (double) max_found[i] / (double) strobe_nr;
			sensitivity[i] += ratio;

			// Why does 0.05 work so well?
			if ( ratio > 0.1)
			{
				if (ratio > max_ratio)
					max_ratio = ratio;

				potential_indexes.emplace_back(std::make_pair(i, ratio));
			}
		}

		int spec_count = 0;
		for (const auto & p : potential_indexes)
		{
			if (p.second >= max_ratio)
			{
				classified[p.first] += 1;
				if (++spec_count > 1)
					std::cout << header << std::endl;
			}
		}
}
*/