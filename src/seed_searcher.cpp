/*
	seed_searcher.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "seed_searcher.hpp"
#include "ungapped_extender.hpp"

#include <algorithm>
#include <iostream>
#include <utility>

#include <omp.h>

using namespace std;

void SeedSearcher::Load(std::ifstream &ifs, uint32_t starts_size_check){
	uint32_t num;

	ifs.read((char*)&num, sizeof(uint32_t));
	database_positions.resize(num);
	front_indexes.resize(num);
	back_indexes.resize(num);
	ifs.read((char*)&database_positions[0], sizeof(uint32_t)*num);
	ifs.read((char*)&front_indexes[0], sizeof(uint16_t)*num);
	ifs.read((char*)&back_indexes[0], sizeof(uint16_t)*num);

	ifs.read((char*)&num, sizeof(uint32_t));
	starts.resize(num);
	ifs.read((char*)&starts[0], sizeof(uint32_t)*num);

	if(starts.size() != starts_size_check+1){
		cerr << "[ERROR] database file is bloken" << endl;
		exit(1);
	}
}

void SeedSearcher::Write(std::ofstream &ofs){
	//database_positions
	for(const auto& seed : seeds){
		for(const auto& e : seed){
			database_positions.emplace_back(e.database_position);
		}
	}

	uint32_t num = database_positions.size();
	ofs.write((char*)&num, sizeof(uint32_t));
	ofs.write((char*)&database_positions[0], sizeof(uint32_t)*num);
	database_positions.clear();
	database_positions.shrink_to_fit();

	//front_indexes
	for(const auto& seed : seeds){
		for(const auto& e : seed){
			front_indexes.emplace_back(e.front_index);
		}
	}

	ofs.write((char*)&front_indexes[0], sizeof(uint16_t)*num);
	front_indexes.clear();
	front_indexes.shrink_to_fit();

	//back_indexes
	for(const auto& seed : seeds){
		for(const auto& e : seed){
			back_indexes.emplace_back(e.back_index);
		}
	}

	ofs.write((char*)&back_indexes[0], sizeof(uint16_t)*num);
	back_indexes.clear();
	back_indexes.shrink_to_fit();

	uint32_t start = 0;
	for(const auto& seed : seeds){
		starts.emplace_back(start);
		start += seed.size();
	}
	starts.emplace_back(start);

	num = starts.size();
	ofs.write((char*)&num, sizeof(uint32_t));
	ofs.write((char*)&starts[0], sizeof(uint32_t)*num);
}

void SeedSearcher::BuildParallel(const ReducedAlphabetCoder &seed1_coder, const ReducedAlphabetCoder &seed2_coder,
		const uint32_t seq_length, const int seed1_length, const int seed2_length, const AlphabetCoder::Code *seq,
		const AlphabetCoder::Code *seq_seed1, const AlphabetCoder::Code *seq_seed2){

#if defined(DEBUG) && defined(_OPENMP)
	double s, e;
	s = omp_get_wtime();
	cout << endl;
#endif

	uint32_t max_seed1_index = seed1_coder.GetMaxSeedIndex();
	seeds = vector<vector<Seed> >(max_seed1_index,
					vector<Seed>(0));

	uint32_t number_letters = score_matrix.GetNumberLetters();
	const int *sm = score_matrix.GetMatrix();

	//各threadが管理するseed1 indexの範囲サイズ
	uint32_t thread_handle_index = (max_seed1_index/threads) + 1;

	vector<vector<vector<pair<uint32_t,Seed> > > > threads_seeds =
			vector<vector<vector<pair<uint32_t,Seed> > > >(threads,
			vector<vector<pair<uint32_t, Seed> > >(threads,
			vector<pair<uint32_t, Seed> >(0)));

	#pragma omp parallel num_threads(threads)
	{
		uint32_t start = 0;
		int thread_id = 0;

#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#endif

		int seed_score = 0;
		int pass = 0;

		#pragma omp for //pass計算のため、各threadsの初期位置のみ保存
		for(uint32_t i=seq_length-seed1_length-seed2_length-1; i>(uint32_t)seed2_length; --i){
			if(start == 0) start = i;
		}

		for(uint32_t i=start+seed1_length+seed2_length-1; i>=start-seed2_length+1; --i){
			if(seq_seed1[i] == delimiter){
				pass = seed1_length + 2*seed2_length;
			}
			if(pass){
				--pass;
			}
		}

#if defined(DEBUG) && defined(_OPENMP)
		#pragma omp critical
		{
			cout << "thread " << thread_id << ", start " << start << ", pass " << pass << endl;
		}
#endif

		#pragma omp for
		for(uint32_t i=seq_length-seed1_length-seed2_length-1; i>(uint32_t)seed2_length; --i){
			if(seq_seed1[i-seed2_length] == delimiter){
				pass = seed1_length + 2*seed2_length;
			}

			if(pass == 0){
				seed_score = 0;
				for(int i_seed=0; i_seed<seed1_length; ++i_seed){
					seed_score += sm[seq[i+i_seed] * number_letters + seq[i+i_seed]];
				}
				if(seed_score < seed1_length*4 + seed1_length/2){ //TODO only in BLOSUM62
					continue;
				}

				uint32_t seed_index = seed1_coder.GetSeedIndex(&seq_seed1[i]);
				uint16_t front_index = seed2_coder.GetSeedIndex16(&seq_seed2[i-seed2_length]);
				uint16_t back_index = seed2_coder.GetSeedIndex16(&seq_seed2[i+seed1_length]);

				Seed s = {i, front_index, back_index};

				threads_seeds[thread_id][seed_index/thread_handle_index].emplace_back(make_pair(seed_index, s));
			}else{
				--pass;
				continue;
			}
		}

#if defined(DEBUG) && defined(_OPENMP)
		if(thread_id == 0){
			e = omp_get_wtime();
			cout << endl << "construct threads_seeds: " << e - s << " sec." << endl;
			s = omp_get_wtime();
		}
#endif

		for(int i=0; i<threads; ++i){
			for(const auto& p:threads_seeds[i][thread_id]){
				seeds[p.first].emplace_back(p.second);
			}
		}

#if defined(DEBUG) && defined(_OPENMP)
		if(thread_id == 0){
			e = omp_get_wtime();
			cout << "construct seeds: " << e - s << " sec." << endl;
			s = omp_get_wtime();
		}
#endif

		#pragma omp barrier
		#pragma omp for schedule(dynamic, 10)
		for(uint32_t i=0; i<seeds.size(); ++i){
			std::sort(seeds[i].begin(), seeds[i].end());
		}

#if defined(DEBUG) && defined(_OPENMP)
		if(thread_id == 0){
			e = omp_get_wtime();
			cout << "sort seeds: " << e - s << " sec." << endl << endl;
		}
#endif

	} //end omp parallel
}

void SeedSearcher::Build(const ReducedAlphabetCoder &seed1_coder, const ReducedAlphabetCoder &seed2_coder,
		const uint32_t seq_length, const int seed1_length, const int seed2_length, const AlphabetCoder::Code *seq,
		const AlphabetCoder::Code *seq_seed1, const AlphabetCoder::Code *seq_seed2){

	uint32_t max_seed1_index = seed1_coder.GetMaxSeedIndex();
	seeds = vector<vector<Seed> >(max_seed1_index,
					vector<Seed>(0));

	uint32_t number_letters = score_matrix.GetNumberLetters();
	const int *sm = score_matrix.GetMatrix();
	int seed_score=0;

	int pass = 0;

	for(uint32_t i=seq_length-2; i>=seq_length-seed1_length-seed2_length-seed2_length; --i){
		if(seq_seed1[i] == delimiter){
			pass = seed1_length + 2*seed2_length;
		}
		if(pass){
			--pass;
		}
	}

	for(uint32_t i=seq_length-seed1_length-seed2_length-1; i>(uint32_t)seed2_length; --i){
		if(seq_seed1[i-seed2_length] == delimiter){
			pass = seed1_length + 2*seed2_length;
		}

		if(pass == 0){
			
			seed_score = 0;
			for(int i_seed=0; i_seed<seed1_length; ++i_seed){
				seed_score += sm[seq[i+i_seed] * number_letters + seq[i+i_seed]];
			}
			if(seed_score < seed1_length*4 + seed1_length/2){ //TODO only in BLOSUM62
				continue;
			}

			uint32_t seed_index = seed1_coder.GetSeedIndex(&seq_seed1[i]);
			uint16_t front_index = static_cast<uint16_t>(seed2_coder.GetSeedIndex(&seq_seed2[i-seed2_length]));
			uint16_t back_index = static_cast<uint16_t>(seed2_coder.GetSeedIndex(&seq_seed2[i+seed1_length]));

			Seed s = {i, front_index, back_index};
			seeds[seed_index].emplace_back(s);
		}else{
			--pass;
			continue;
		}
	}

	for(auto& e : seeds){
		std::sort(e.begin(), e.end());
	}
}

vector<SeedSearcher::Hit> SeedSearcher::Search(const std::vector<uint32_t>& indexes_,
		const std::vector<uint16_t>& front_indexes_, const std::vector<uint16_t>& back_indexes_,
		const AlphabetCoder::Code *query, const uint32_t query_position, const AlphabetCoder::Code *database,
		const int ungapped_extension_cutoff, const int gapped_extension_trigger) const{
	UngappedExtender ungapped_extender(score_matrix, delimiter);
	vector<SeedSearcher::Hit> hits;

	for(const uint32_t& index : indexes_){
		for(const uint16_t& front_index : front_indexes_){

			auto itr_s = std::lower_bound(front_indexes.begin()+starts[index],
					front_indexes.begin()+starts[index+1], front_index);

			if(itr_s == front_indexes.begin()+starts[index+1] || *itr_s != front_index){
				continue;
			}

			auto itr_e = std::upper_bound(front_indexes.begin()+starts[index],
					front_indexes.begin()+starts[index+1], front_index);

			size_t search_now = itr_s - front_indexes.begin();
			size_t search_end = itr_e - front_indexes.begin();

			auto back_itr = back_indexes_.begin();

			for(; search_now < search_end; ++search_now){
				if(back_indexes[search_now] < *back_itr){
					
				}else if(back_indexes[search_now] == *back_itr){
					int score = ungapped_extender.UngappedExtend(&query[0], query_position, &database[0],
							database_positions[search_now], ungapped_extension_cutoff, gapped_extension_trigger);
					if(score > gapped_extension_trigger){
						SeedSearcher::Hit hit = {query_position, database_positions[search_now]};
						hits.emplace_back(hit);
					}
				}else{
					++back_itr;
					--search_now;
					if(back_itr == back_indexes_.end()){
						break;
					}
				}
			}
		}
	}
	return hits;
}

void SeedSearcher::Clear(){
	seeds.clear();
	seeds.shrink_to_fit();

	front_indexes.clear();
	back_indexes.clear();
	database_positions.clear();

	starts.clear();
	starts.shrink_to_fit();
}
