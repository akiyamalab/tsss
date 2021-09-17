/*
	aligner.hpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#ifndef ALIGNER_HPP
#define ALIGNER_HPP

#include "query.hpp"
#include "queries.hpp"
#include "score_matrix.hpp"
#include "alphabet_coder.hpp"
#include "alphabet_type.hpp"
#include "edit_blocks.hpp"
#include "database_chunk.hpp"

#include <string>
#include <tr1/memory>

class Aligner{
public:

	enum Outfmt {
		pairwise,
		tabular,
	};

	struct Parameters{
		Outfmt outfmt;
		ScoreMatrix score_matrix;
		AlphabetCoder coder;
		std::string query_path;
		std::string database_path;
		std::string output_path;
		uint32_t chunk_size;
		int num_alignments;
		int threads;
		int seed1_hamming;
		int seed2_hamming;
		AlphabetCoder::Code delimiter;
		std::tr1::shared_ptr<AlphabetType> db_type;
		std::tr1::shared_ptr<AlphabetType> query_type;
		bool filter;
		double evalue;
		int gap_open;
		int gap_extension;
		float normalized_presearched_ungapped_extension_cutoff;
		float normalized_presearched_gapped_extension_trigger;
		float normalized_presearched_gapped_extension_cutoff;
		float normalized_result_gapped_extension_cutoff;
		uint64_t db_total_length;
		uint64_t db_total_sequences;
	}parameters;

	void Align();

	typedef struct {
		uint32_t query_position;
		uint32_t database_position;
	} Coordinate;

	typedef struct {
		uint32_t chunk_id;
		int score;
		Coordinate hit;
		Coordinate start;
		Coordinate end;
	} PresearchedResult;

	typedef struct {
		std::string subject_name;
		int score;
		uint32_t alignment_length;
		uint32_t mismatches;
		uint32_t positives;
		uint32_t gap_openings;
		uint32_t gap_extensions;
		Coordinate hit;
		Coordinate start;
		Coordinate end;
		std::string pairwise_result;
	} Result;

	struct PresearchedResultGreaterScore{
		bool operator()(const PresearchedResult &r1, const PresearchedResult &r2){
			return AlignmentComp(r1.score, r1.start, r1.end, r2.score, r2.start, r2.end) > 0;
		}
	};

	struct ResultGreaterScore{
		bool operator()(const Result &r1, const Result &r2){
			return AlignmentComp(r1.score, r1.start, r1.end, r2.score, r2.start, r2.end) > 0;
		}
	};

private:
	static int AlignmentComp(const int a_score, const Coordinate &a_start, const Coordinate &a_end,
			const int b_score, const Coordinate &b_start, const Coordinate &b_end);
	void BuildQueriesParameters(Queries::Parameters &queries_parameters, const int seed_length);
	void Presearch(Queries &queries, DatabaseChunk& chunk, std::vector<std::vector<PresearchedResult> > &result_list);
	void BuildResults(Queries &queries, DatabaseChunk& chunk, std::vector<std::vector<PresearchedResult> > &presearch_result_list,
			std::vector<std::vector<Result> > &result_list);
	void BuildResult(Query &query, DatabaseChunk &chunk, int score, Coordinate &hit, Coordinate &start,
			Coordinate &end, std::vector<EditBlocks::EditOpType> &edits, Result &result);
	void AddResults(std::vector<PresearchedResult> &added_results, std::vector<PresearchedResult> &results);
	void WriteOutput(std::ofstream &ofs, Queries &queries, std::vector<std::vector<Result> > &results);
};

#endif
