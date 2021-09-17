/*
	reduced_alphabet_coder.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "reduced_alphabet_coder.hpp"
#include "protein_type.hpp"

#include <iostream>
#include <algorithm>

using namespace std;

namespace{
	void recursive_comb(vector<uint8_t> &position, vector<vector<uint8_t> > &positions, int s, int rest) {
		if (rest == 0) {
			positions.emplace_back(position);
		} else {
			if (s < 0) return;
			recursive_comb(position, positions, s - 1, rest);
			position[rest - 1] = s;
			recursive_comb(position, positions, s - 1, rest - 1);
		}
	}
};

ReducedAlphabetCoder::ReducedAlphabetCoder(const string &similarity_sets, const int seed_length_, int seed_hamming_)
	:seed_length(seed_length_),seed_hamming(seed_hamming_),num_reduced_amino(-1){
	if(seed_length <= seed_hamming){
		cerr << "[ERROR] seed_hamming " << seed_hamming << " must be under seed_length " << seed_length << endl;
		exit(1);
	}

	ProteinType type;
	AlphabetCoder coder(type);

	AlphabetCoder::Code delimiter = coder.GetMaxCode()+1;
	code_map.resize((int)delimiter+1, delimiter); //defなら0-25,初期値25(delimiter)
	indexer.resize(seed_length);

	char rep = '\0';
	for(const auto c:similarity_sets){
		if(c == ' '){
			rep = '\0';
		}else if(rep == '\0'){
			rep = c;
			++num_reduced_amino;
			code_map[coder.Encode(c)] = num_reduced_amino;
		}else{
			code_map[coder.Encode(c)] = num_reduced_amino;
		}
	}
	++num_reduced_amino;

	uint32_t num = 1;
	for(int i=0; i<seed_length; ++i){
		indexer[i] = num;
		num *= num_reduced_amino;
	}
	max_seed_index = num;

	if(seed_hamming > 1){

		vector<uint8_t> neighbor_amino(seed_hamming, 0);

		//repeated permutation
		bool flag = true;
		while(flag){
			neighbor_aminos.emplace_back(neighbor_amino);
			int index = seed_hamming -1;

			while(flag){
				if(++neighbor_amino[index] < num_reduced_amino) break;
				neighbor_amino[index] = 0;
				if(--index < 0) flag = false;
			}
		}

		vector<uint8_t> change_position(seed_hamming, 0);
		recursive_comb(change_position, change_positions,  seed_length-1, seed_hamming);
	}
}

vector<uint32_t> ReducedAlphabetCoder::GetNeighborIndexes(const AlphabetCoder::Code *seed) const{

	vector<uint32_t> seed_indexes;

	if(seed_hamming == 0){
		seed_indexes.emplace_back(ReducedAlphabetCoder::GetSeedIndex(&seed[0]));
		return seed_indexes;
	}

	vector<AlphabetCoder::Code> seed_(seed_length);

	for(int i=0; i<seed_length; ++i){
		seed_[i] = seed[i];
	}

	if(seed_hamming == 1){
		seed_indexes.emplace_back(ReducedAlphabetCoder::GetSeedIndex(&seed[0]));
		for(int i=0; i<seed_length; ++i){
			for(int j=0; j<num_reduced_amino-1; ++j){
				++seed_[i];
				if(seed_[i]==num_reduced_amino){
					seed_[i]=0;
				}
				seed_indexes.emplace_back(ReducedAlphabetCoder::GetSeedIndex(&seed_[0]));
			}
			++seed_[i];
			if(seed_[i]==num_reduced_amino){
				seed_[i] = 0;
			}
		}
		return seed_indexes;
	}

	//hamming >= 2
	for(const auto& pos:change_positions){
		vector<AlphabetCoder::Code> temp_seed(seed_);
		for(const auto& amino:neighbor_aminos){

			for(uint8_t i=0; i<seed_hamming; ++i){
				temp_seed[pos[i]] = amino[i];
			}
			seed_indexes.emplace_back(ReducedAlphabetCoder::GetSeedIndex(&temp_seed[0]));

		}
	}

	std::sort(seed_indexes.begin(), seed_indexes.end());
	seed_indexes.erase(std::unique(seed_indexes.begin(), seed_indexes.end()), seed_indexes.end());
	return seed_indexes;
}

vector<uint16_t> ReducedAlphabetCoder::GetNeighborIndexes16(const AlphabetCoder::Code *seed) const{

	vector<uint16_t> seed_indexes;

	if(seed_hamming == 0){
		seed_indexes.emplace_back(ReducedAlphabetCoder::GetSeedIndex16(&seed[0]));
		return seed_indexes;
	}

	vector<AlphabetCoder::Code> seed_(seed_length);

	for(int i=0; i<seed_length; ++i){
		seed_[i] = seed[i];
	}

	if(seed_hamming == 1){
		seed_indexes.emplace_back(ReducedAlphabetCoder::GetSeedIndex16(&seed[0]));
		for(int i=0; i<seed_length; ++i){
			for(int j=0; j<num_reduced_amino-1; ++j){
				++seed_[i];
				if(seed_[i]==num_reduced_amino){
					seed_[i]=0;
				}
				seed_indexes.emplace_back(ReducedAlphabetCoder::GetSeedIndex16(&seed_[0]));
			}
			++seed_[i];
			if(seed_[i]==num_reduced_amino){
				seed_[i] = 0;
			}
		}
		return seed_indexes;
	}

	//hamming >= 2
	for(const auto& pos:change_positions){
		vector<AlphabetCoder::Code> temp_seed(seed_);
		for(const auto& amino:neighbor_aminos){

			for(uint8_t i=0; i<seed_hamming; ++i){
				temp_seed[pos[i]] = amino[i];
			}
			seed_indexes.emplace_back(ReducedAlphabetCoder::GetSeedIndex16(&temp_seed[0]));

		}
	}

	std::sort(seed_indexes.begin(), seed_indexes.end());
	seed_indexes.erase(std::unique(seed_indexes.begin(), seed_indexes.end()), seed_indexes.end());
	return seed_indexes;
}

void ReducedAlphabetCoder::Encode(const AlphabetCoder::Code *before, const unsigned int length,
		AlphabetCoder::Code *after) const{
	for(uint32_t i=0; i<length; ++i){
		after[i] = code_map[before[i]];
	}
}

uint32_t ReducedAlphabetCoder::GetSeedIndex(const AlphabetCoder::Code *seed) const{
	uint32_t index=0;
	for(int i=0; i<seed_length; ++i){
		index += seed[i]*indexer[i];
	}
	return index;
}

uint16_t ReducedAlphabetCoder::GetSeedIndex16(const AlphabetCoder::Code *seed) const{
	uint16_t index=0;
	for(int i=0; i<seed_length; ++i){
		index += seed[i]*indexer[i];
	}
	return index;
}

void ReducedAlphabetCoder::Profile() const{
	std::cout << "seed length: " << seed_length << ", reduced amino: " << num_reduced_amino;
}
