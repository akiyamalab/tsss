/*
	database_chunk.hpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#ifndef DATABASE_CHUNK_HPP
#define DATABASE_CHUNK_HPP

#include "sequence.hpp"
#include "alphabet_coder.hpp"
#include "reduced_alphabet_coder.hpp"
#include "fasta_sequence_reader.hpp"
#include "seed_searcher.hpp"
#include "score_matrix.hpp"

#include <string>
#include <vector>
#include <fstream>

class DatabaseChunk{
public:

	typedef struct {
		ScoreMatrix score_matrix;
		std::string database_path;
		std::string output_path;
		uint32_t chunk_size;
		int seed1_length;
		int seed1_amino;
		int seed1_hamming;
		int seed2_length;
		int seed2_amino;
		int seed2_hamming;
		int threads;
		AlphabetCoder::Code delimiter;
	}Parameters;

	DatabaseChunk(){};
	~DatabaseChunk(){};
	DatabaseChunk(std::ifstream &in, const Parameters &parameters_);
	DatabaseChunk(const Parameters &parameters_, uint64_t &db_total_length_, uint64_t &db_total_sequences_);

	bool ReadChunk();
	void EncodeChunk();
	void BuildChunk();
	void SaveInfo() const;
	void SaveChunk();
	bool LoadChunk();
	void ClearChunk();
	inline void ResetChunk(){
		chunk_id=0;
	};

	std::vector<SeedSearcher::Hit> SeedSearch(const AlphabetCoder::Code *query, const uint32_t query_length,
			const int ungapped_cutoff, const int gapped_trigger) const;

	uint32_t GetSubjectId(const uint32_t position) const;

	inline AlphabetCoder::Code *GetConcatenatedSeqRef(){
		return &concatenated_sequence[0];
	}

	inline uint32_t GetConcatenatedSeqLength() const{
		return concatenated_seq_length;
	}

	inline int GetChunkId() const{
		return chunk_id;
	}

	inline int GetNumDbChunk() const {
		return num_db_chunk;
	}

	inline std::string GetSubjectName(const uint32_t id) const{
		return sequence_names[id];
	}

	inline uint32_t GetOffset(const uint32_t id) const{
		return offsets[id];
	}

	inline int GetSeedLength() const{
		return parameters.seed1_length + 2 * parameters.seed2_length;
	}

private:
	Parameters parameters;

	FastaSequenceReader reader;
	SeedSearcher searcher;
	AlphabetCoder coder;
	ReducedAlphabetCoder seed1_coder;
	ReducedAlphabetCoder seed2_coder;

	std::vector<Sequence> sequences;
	std::vector<std::string> sequence_names;
	std::vector<AlphabetCoder::Code> concatenated_sequence;
	std::vector<AlphabetCoder::Code> concatenated_seq_seed1;
	std::vector<AlphabetCoder::Code> concatenated_seq_seed2;
	std::vector<uint32_t> offsets;
	uint32_t total_length;
	int concatenated_seq_length;
	int chunk_id;
	int num_db_chunk;
	uint64_t db_total_length;
	uint64_t db_total_sequences;
};

#endif
