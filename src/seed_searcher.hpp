/*
	seed_searcher.hpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#ifndef SEED_SEARCHER_HPP
#define SEED_SEARCHER_HPP

#include "alphabet_coder.hpp"
#include "reduced_alphabet_coder.hpp"
#include "score_matrix.hpp"

#include <fstream>
#include <vector>

struct HitBase{
	uint32_t query_position;
	uint32_t database_position;

	bool operator<(const HitBase& another) const{
		return database_position < another.database_position;
	}
};

struct SeedBase{
	uint32_t database_position;
	uint16_t front_index;
	uint16_t back_index;

	bool operator<(const SeedBase& another) const{
		if(front_index == another.front_index){
			return back_index < another.back_index;
		}
		return front_index < another.front_index;
	}
};

class SeedSearcher{
public:

	typedef struct SeedBase Seed;
	typedef struct HitBase Hit;

	SeedSearcher(const ScoreMatrix &score_matrix_, const AlphabetCoder::Code delimiter_, const int threads_):
		score_matrix(score_matrix_),delimiter(delimiter_),threads(threads_){};
	SeedSearcher(){};
	~SeedSearcher(){};

	void Load(std::ifstream &ifs, const uint32_t starts_size_check);
	void Write(std::ofstream &ofs);
	void BuildParallel(const ReducedAlphabetCoder &seed1_coder, const ReducedAlphabetCoder &seed2_coder,
			const uint32_t seq_length, const int seed1_length, const int seed2_length, const AlphabetCoder::Code *seq,
			const AlphabetCoder::Code *seq_seed1, const AlphabetCoder::Code *seq_seed2);
	void Build(const ReducedAlphabetCoder &seed1_coder, const ReducedAlphabetCoder &seed2_coder,
			const uint32_t seq_length, const int seed1_length, const int seed2_length, const AlphabetCoder::Code *seq,
			const AlphabetCoder::Code *seq_seed1, const AlphabetCoder::Code *seq_seed2);

	//queryのfrontとbackを持ってきて、ヒットするseedのindexの組を返す
	std::vector<SeedSearcher::Hit> Search(const std::vector<uint32_t>& indexes, const std::vector<uint16_t>& front_indexes,
			const std::vector<uint16_t>& back_indexes, const AlphabetCoder::Code *query, const uint32_t query_position,
			const AlphabetCoder::Code *database, const int ungapped_extension_cutoff, const int gapped_extension_trriger) const;

	void Clear();

private:
	ScoreMatrix score_matrix;
	AlphabetCoder::Code delimiter;
	int threads;

	std::vector<std::vector<Seed> > seeds; //only used in db
	std::vector<uint16_t> front_indexes;
	std::vector<uint16_t> back_indexes;
	std::vector<uint32_t> database_positions;
	std::vector<uint32_t> starts;
};

#endif
