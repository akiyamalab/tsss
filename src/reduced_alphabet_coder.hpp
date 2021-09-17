/*
	reduced_alphabet_coder.hpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#ifndef REDUCED_ALPHABET_CODER_HPP
#define REDUCED_ALPHABET_CODER_HPP

#include "alphabet_coder.hpp"

#include <string>
#include <vector>

class ReducedAlphabetCoder{
public:
	ReducedAlphabetCoder(){};
	~ReducedAlphabetCoder(){};
	ReducedAlphabetCoder(const std::string &similarity_sets, int seed_length_, int seed_hamming_);

	std::vector<uint32_t> GetNeighborIndexes(const AlphabetCoder::Code *seed) const ;
	std::vector<uint16_t> GetNeighborIndexes16(const AlphabetCoder::Code *seed) const ;
	void Encode(const AlphabetCoder::Code *before, const unsigned int length, AlphabetCoder::Code *after) const;
	uint32_t GetSeedIndex(const AlphabetCoder::Code *seed) const;
	uint16_t GetSeedIndex16(const AlphabetCoder::Code *seed) const;
	inline uint32_t GetMaxSeedIndex() const{
		return max_seed_index;
	};
	void Profile() const;

private:
	int seed_length;
	int seed_hamming;
	int num_reduced_amino;
	uint32_t max_seed_index;
	std::vector<AlphabetCoder::Code> code_map;
	std::vector<uint32_t> indexer;

	std::vector<std::vector<uint8_t> > neighbor_aminos;
	std::vector<std::vector<uint8_t> > change_positions;
};

#endif
