/*
	aligner.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "aligner.hpp"
#include "statistics.hpp"
#include "seed_searcher.hpp"
#include "chain_filter.hpp"
#include "gapped_extender.hpp"
#include "utils.hpp"
#include "protein_type.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>

#include <omp.h>

using namespace std;

namespace {
	bool SameSubject(const AlphabetCoder::Code *seq, const uint32_t before, const uint32_t now,
			const AlphabetCoder::Code delim){
		for(uint32_t i=before; i<now; ++i){
			if(seq[i] == delim){
				return false;
			}
		}
		return true;
	};
};

void Aligner::Align(){
	ofstream ofs;
	utils::FileOpen(ofs, parameters.output_path, ios_base::out);

	DatabaseChunk::Parameters parameters_;
	parameters_.score_matrix  = parameters.score_matrix;
	parameters_.delimiter     = parameters.delimiter;
	parameters_.database_path = parameters.database_path;
	parameters_.seed1_hamming = parameters.seed1_hamming;
	parameters_.seed2_hamming = parameters.seed2_hamming;
	DatabaseChunk chunk(parameters_, parameters.db_total_length, parameters.db_total_sequences);

	Queries::Parameters queries_parameters;
	ifstream queries_is;
	utils::FileOpen(queries_is, parameters.query_path, ios_base::in);
	BuildQueriesParameters(queries_parameters, chunk.GetSeedLength());

	double start,end;

	for(Queries queries(queries_is, queries_parameters); queries.GetNumberSequences() != 0; queries.Next()){
		chunk.ResetChunk();

		cout << "number queries is " << queries.GetNumberSequences() << endl << endl;
		vector<vector<PresearchedResult> > presearch_results_list(queries.GetNumberSequences());
		vector<vector<Result> > results_list(queries.GetNumberSequences());

		cout << "starts presearch" << endl;
		start = omp_get_wtime();
		Presearch(queries, chunk, presearch_results_list);
		end = omp_get_wtime();
		cout << "finish presearch     [time] " << end - start << "sec." << endl << endl;

		chunk.ResetChunk();

		cout << "starts build results" << endl;
		start = omp_get_wtime();
		BuildResults(queries, chunk, presearch_results_list, results_list);
		end = omp_get_wtime();
		cout << "finish build results [time] " << end - start << "sec." << endl << endl;

		cout << "starts write results" << endl;
		start = omp_get_wtime();
		WriteOutput(ofs, queries, results_list);
		end = omp_get_wtime();
		cout << "finish write results [time] " << end - start << "sec." << endl << endl;
	}
	queries_is.close();
	ofs.close();
}

void Aligner::Presearch(Queries &queries, DatabaseChunk& chunk, vector<vector<PresearchedResult> > &result_list){

	Statistics statistics(*(parameters.db_type));
	Statistics::KarlinParameters gapped_karlin_parameters;
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,&gapped_karlin_parameters);
	int gapped_extension_cutoff = Statistics::Normalized2Nominal(
			parameters.normalized_presearched_gapped_extension_cutoff,gapped_karlin_parameters);

	vector<vector<SeedSearcher::Hit> > hits_list(queries.GetNumberSequences());

	std::vector<int> ungapped_extension_cutoffs(queries.GetNumberSequences());
	std::vector<int> gapped_extension_triggers(queries.GetNumberSequences());

	for (uint32_t query_id = 0; query_id < queries.GetNumberSequences(); ++query_id) {
		Query *query = queries.GetQuery(query_id);
		ungapped_extension_cutoffs[query_id] =Statistics::NormalizedCutoff2NominalCutoff(
																					parameters.normalized_presearched_ungapped_extension_cutoff,
																					query->GetUngappedKarlinParameters());
		gapped_extension_triggers[query_id] = Statistics::Normalized2Nominal(
																					parameters.normalized_presearched_gapped_extension_trigger,
																					query->GetUngappedKarlinParameters());
	}

#ifdef DEBUG
	uint32_t sum_seed_search_hit[parameters.threads];
	uint32_t sum_chain_filter_hit[parameters.threads];
	uint32_t sum_presearched_result[parameters.threads];

	double seed_search_time[parameters.threads];
	double chain_filter_time[parameters.threads];
	double gapped_extension_time[parameters.threads];
	double sort_time[parameters.threads];

	for(int i=0; i<parameters.threads; ++i){
		sum_seed_search_hit[i] = 0;
		sum_chain_filter_hit[i] = 0;
		sum_presearched_result[i] = 0;
		seed_search_time[i] = 0.;
		chain_filter_time[i] = 0.;
		gapped_extension_time[i] = 0.;
		sort_time[i] = 0.;
	}
#endif

#ifdef _OPENMP
	omp_set_num_threads(parameters.threads);
#endif

	for(; chunk.LoadChunk(); chunk.ClearChunk()){
		int chunk_id = chunk.GetChunkId();

		cout << "load chunk " << chunk_id << endl;

		const AlphabetCoder::Code *concatenated_seq = chunk.GetConcatenatedSeqRef();
		const uint32_t concatenated_seq_length = chunk.GetConcatenatedSeqLength();

		#pragma omp parallel
		{
			ChainFilter chain_filter(parameters.delimiter, parameters.score_matrix);
			GappedExtender gapped_extender;
			vector<PresearchedResult> temp_results;

#ifdef DEBUG
			double start, end;
			int thread_id = omp_get_thread_num();
#endif

			#pragma omp for schedule(dynamic, 10)
			for(uint32_t query_id = 0; query_id < queries.GetNumberSequences(); ++query_id) {
				Query *query = queries.GetQuery(query_id);

				temp_results.clear();
				PresearchedResult temp_result;
				temp_result.chunk_id = chunk_id;

#ifdef DEBUG
				start = omp_get_wtime();
#endif
				vector<SeedSearcher::Hit> hits = chunk.SeedSearch(query->GetSequence(), query->GetSequenceLength(),
						ungapped_extension_cutoffs[query_id], gapped_extension_triggers[query_id]);

#ifdef DEBUG
				end = omp_get_wtime();
				seed_search_time[thread_id] += end - start;
				start = omp_get_wtime();

				sum_seed_search_hit[thread_id] += hits.size();
#endif

				chain_filter.Filter(query->GetSequence(), query->GetSequenceLength(), ungapped_extension_cutoffs[query_id],
						concatenated_seq, concatenated_seq_length, hits);

#ifdef DEBUG
				end = omp_get_wtime();
				chain_filter_time[thread_id] += end - start;
				start = omp_get_wtime();
				sum_chain_filter_hit[thread_id] += hits.size();
#endif

				for(size_t hit_id=0; hit_id<hits.size(); ++hit_id){
					int score = 0;
					int query_position;
					int database_position;
					temp_result.score = 0;
					temp_result.hit.query_position = hits[hit_id].query_position;
					temp_result.hit.database_position = hits[hit_id].database_position;

					gapped_extender.ExtendOneSide(query->GetSequence()+temp_result.hit.query_position - 1,
							temp_result.hit.query_position - 1, &concatenated_seq[temp_result.hit.database_position - 1],
							parameters.delimiter, true, parameters.score_matrix, parameters.gap_open, parameters.gap_extension,
							gapped_extension_cutoff, &score, &query_position, &database_position, NULL);

					temp_result.score += score;
					temp_result.start.query_position = temp_result.hit.query_position -1 +query_position;
					temp_result.start.database_position = temp_result.hit.database_position -1 +database_position;

					gapped_extender.ExtendOneSide(query->GetSequence()+temp_result.hit.query_position,
							query->GetSequenceLength() - temp_result.hit.query_position,
							&concatenated_seq[temp_result.hit.database_position],
							parameters.delimiter, false, parameters.score_matrix, parameters.gap_open, parameters.gap_extension,
							gapped_extension_cutoff, &score, &query_position, &database_position, NULL);

					temp_result.score += score;
					temp_result.end.query_position += query_position + 1;
					temp_result.end.database_position += database_position + 1;

					if(hit_id > 0 && SameSubject(concatenated_seq, temp_results[temp_results.size()-1].hit.database_position,
								temp_result.hit.database_position, parameters.delimiter)){
						if(temp_results[temp_results.size()-1].score < temp_result.score){
							temp_results.pop_back();
							temp_results.push_back(temp_result);
						}
					}else{
						temp_results.push_back(temp_result);
					}
				}

#ifdef DEBUG
				end = omp_get_wtime();
				gapped_extension_time[thread_id] += end - start;
				start = omp_get_wtime();
				sum_presearched_result[thread_id] += temp_results.size();
#endif

				//上位n件をソート
				if((int)temp_results.size() > parameters.num_alignments){
					partial_sort(temp_results.begin(), temp_results.begin()+parameters.num_alignments, temp_results.end(),
							PresearchedResultGreaterScore());
					temp_results.erase(temp_results.begin() + parameters.num_alignments, temp_results.end());
				}else{
					sort(temp_results.begin(), temp_results.end(), PresearchedResultGreaterScore());
				}

				AddResults(temp_results, result_list[query_id]);

#ifdef DEBUG
				end = omp_get_wtime();
				sort_time[thread_id] += end - start;
#endif

			}
		} // end omp parallel
	}

#ifdef DEBUG
	for(int i=1; i<parameters.threads; ++i){
		sum_seed_search_hit[0] += sum_seed_search_hit[i];
		sum_chain_filter_hit[0] += sum_chain_filter_hit[i];
		sum_presearched_result[0] += sum_presearched_result[i];

		seed_search_time[0] += seed_search_time[i];
		chain_filter_time[0] += chain_filter_time[i];
		gapped_extension_time[0] += gapped_extension_time[i];
		sort_time[0] += sort_time[i];
	}

	cout << "average of seed search hits per query: " << sum_seed_search_hit[0] / queries.GetNumberSequences() << endl;
	cout << "average of chain filter hits per query: " << sum_chain_filter_hit[0] / queries.GetNumberSequences() << endl;
	cout << "average of presearch results per query: " << sum_presearched_result[0] / queries.GetNumberSequences() << endl;

	cout << "total seed search time: " << seed_search_time[0] << " .sec" << endl;
	cout << "total chain filter time: " << chain_filter_time[0] << " .sec" << endl;
	cout << "total gapped extension time: " << gapped_extension_time[0] << " .sec" << endl;
	cout << "total sort time: " << sort_time[0] << " .sec" << endl;
#endif

}

void Aligner::BuildResults(Queries &queries, DatabaseChunk& chunk, vector<vector<PresearchedResult> > &presearch_result_list,
		vector<vector<Result> > &result_list){
	Statistics statistics(*(parameters.db_type));
	Statistics::KarlinParameters gapped_karlin_parameters;
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension, &gapped_karlin_parameters);
	const int gapped_extension_cutoff = Statistics::Normalized2Nominal(
			parameters.normalized_result_gapped_extension_cutoff, gapped_karlin_parameters);

	vector<vector<pair<uint32_t, uint32_t> > > result_ids_list(chunk.GetNumDbChunk());

	if(presearch_result_list.size() < queries.GetNumberSequences()){
		presearch_result_list.resize(queries.GetNumberSequences());
	}
	for(uint32_t query_id = 0; query_id < queries.GetNumberSequences(); ++query_id){
		result_list[query_id].resize(presearch_result_list[query_id].size());
		for(uint32_t result_id = 0; result_id < presearch_result_list[query_id].size(); ++result_id){
			result_ids_list[presearch_result_list[query_id][result_id].chunk_id].push_back(make_pair(query_id, result_id));
		}
	}

#ifdef _OPENMP
	omp_set_num_threads(parameters.threads);
#endif

	for(; chunk.LoadChunk(); chunk.ClearChunk()){
		const AlphabetCoder::Code *concatenated_seq = chunk.GetConcatenatedSeqRef();
		int chunk_id = chunk.GetChunkId();
		vector<pair<uint32_t, uint32_t> > &result_ids = result_ids_list[chunk_id];
		const size_t result_ids_size = result_ids.size();

		//#pragma omp parallel
		//{
		EditBlocks edit_blocks;
		EditBlocks temp_edit_blocks;
		GappedExtender gapped_extender;

		//#pragma omp for schedule(dynamic, 10)
		for(size_t i=0; i<result_ids_size; ++i){
			edit_blocks.Clear();
			uint32_t query_id = result_ids[i].first;
			uint32_t result_id = result_ids[i].second;
			Query *query = queries.GetQuery(query_id);

			int sum_score = 0;
			Coordinate hit;
			Coordinate start;
			Coordinate end;

			PresearchedResult *presearched_result_ptr = &presearch_result_list[query_id][result_id];
			hit.query_position = presearched_result_ptr->hit.query_position;
			hit.database_position = presearched_result_ptr->hit.database_position;
			int score;
			int query_position;
			int database_position;

			gapped_extender.ExtendOneSide(query->GetSequence() + hit.query_position -1, hit.query_position -1,
					concatenated_seq + hit.database_position -1, parameters.delimiter, true,
					parameters.score_matrix, parameters.gap_open, parameters.gap_extension,
					gapped_extension_cutoff, &score, &query_position, &database_position,
					&temp_edit_blocks);

			sum_score += score;
			start.query_position = hit.query_position -1 +query_position;
			start.database_position = hit.database_position -1 +database_position;
			temp_edit_blocks.Reverse();
			edit_blocks.Add(temp_edit_blocks);

			gapped_extender.ExtendOneSide(query->GetSequence() + hit.query_position,
					query->GetSequenceLength() - hit.query_position, concatenated_seq + hit.database_position,
					parameters.delimiter, false, parameters.score_matrix, parameters.gap_open,
					parameters.gap_extension, gapped_extension_cutoff, &score, &query_position,
					&database_position, &temp_edit_blocks);

			sum_score += score;
			end.query_position = hit.query_position + query_position;
			end.database_position = hit.database_position + database_position;
			edit_blocks.Add(temp_edit_blocks);
			vector<EditBlocks::EditOpType> edits = edit_blocks.ToVector();

			BuildResult(*query, chunk, sum_score, hit, start, end, edits, result_list[query_id][result_id]);
		}
		//} // end omp parallel
	}

	for(uint32_t query_id=0; query_id < queries.GetNumberSequences(); ++query_id){
		sort(result_list[query_id].begin(), result_list[query_id].end(), ResultGreaterScore());
	}
}

void Aligner::BuildQueriesParameters(Queries::Parameters &queries_parameters, const int seed_length) {
	Statistics statistics(*parameters.db_type);

	if(parameters.filter){
		queries_parameters.filter = true;
	}else{
		queries_parameters.filter = false;
	}

	queries_parameters.file_sequence_type_ptr = parameters.query_type;
	queries_parameters.aligning_sequence_type_ptr = parameters.db_type;
	statistics.CalculateUngappedIdealKarlinParameters(parameters.score_matrix,
			&queries_parameters.ungapped_karlin_parameters);
	queries_parameters.chunk_size = parameters.chunk_size;
	queries_parameters.score_matrix = parameters.score_matrix;
	queries_parameters.sequence_delimiter = parameters.delimiter;

	if(typeid(*(parameters.query_type)) == typeid(ProteinType)){
		queries_parameters.lower_limit_of_sequence_length = seed_length;
	}else{
		queries_parameters.lower_limit_of_sequence_length = seed_length * 3;
	}
}

void Aligner::AddResults(vector<PresearchedResult> &added_results, vector<PresearchedResult> &results){
	vector<PresearchedResult> new_results(0);

	vector<PresearchedResult>::iterator added_itr = added_results.begin();
	vector<PresearchedResult>::iterator added_itr_e = added_results.end();
	vector<PresearchedResult>::iterator itr = results.begin();
	vector<PresearchedResult>::iterator itr_e = results.end();

	while((int)new_results.size() < parameters.num_alignments && (added_itr != added_itr_e || itr != itr_e)){
		if(added_itr == added_itr_e){
			new_results.push_back(*itr);
			++itr;
		}else if(itr == itr_e || added_itr->score > itr->score){
			new_results.push_back(*added_itr);
			++added_itr;
		}else{
			new_results.push_back(*itr);
			++itr;
		}
	}
	results.clear();
	results.insert(results.begin(), new_results.begin(), new_results.end());
}

void Aligner::WriteOutput(ofstream &ofs, Queries &queries, std::vector<std::vector<Result> > &results) {
	Statistics::KarlinParameters gapped_karlin_parameters;
	double alpha;
	double beta;

	Statistics statistics(*(parameters.db_type));
	statistics.CalculateGappedKarlinParameters(parameters.score_matrix,
			parameters.gap_open, parameters.gap_extension,&gapped_karlin_parameters);
	statistics.CalculateAlphaBeta(parameters.score_matrix, parameters.gap_open,
			parameters.gap_extension, &alpha, &beta);

	if(parameters.outfmt == Outfmt::tabular){

		for (uint32_t i = 0; i < queries.GetNumberSequences(); ++i) {
			Query *query = queries.GetQuery(i);
			string query_name = query->GetName();

			uint64_t search_space = statistics.CalculateSearchSpace(
					query->GetRealSequenceLength(), parameters.db_total_length,
					parameters.db_total_sequences, gapped_karlin_parameters, alpha, beta);

			for (vector<Result>::iterator it = results[i].begin(); it != results[i].end(); ++it) {
				double evalue = Statistics::Nominal2EValue(it->score, search_space, gapped_karlin_parameters);
				if(evalue > parameters.evalue){
					break;
				}

				stringstream ss_i, ss_b, ss_e;
				ss_i << std::fixed << std::setprecision(1)
						 << (1. - (static_cast<float>(it->mismatches + it->gap_openings + it->gap_extensions)
														/ it->alignment_length)) * 100;
				ss_b << std::fixed << std::setprecision(1)
						 << Statistics::Nominal2Normalized(it->score, gapped_karlin_parameters);

				ss_e << std::scientific << std::setprecision(1)
						 << evalue;

				ofs << query_name << "\t";
				ofs << it->subject_name << "\t";
				ofs << ss_i.str() << "\t";
				ofs << it->alignment_length << "\t";
				ofs << it->mismatches << "\t";
				ofs << it->gap_openings << "\t";
				ofs << it->start.query_position << "\t";
				ofs << it->end.query_position << "\t";
				ofs << it->start.database_position << "\t";
				ofs << it->end.database_position << "\t";
				ofs << ss_e.str() << "\t";
				ofs	<< ss_b.str();
				//os << "\t(" << it->hit.query_position << ", "
				//		<< it->hit.database_position << ")";
				ofs << '\n';
			}
		}
	}else if(parameters.outfmt == Outfmt::pairwise){

		for (uint32_t i = 0; i < queries.GetNumberSequences(); ++i) {

#ifdef DEBUG
			ofs << string(60, '#') << '\n';
#endif

			Query *query = queries.GetQuery(i);
			string query_name = query->GetName();
			uint32_t query_length = query->GetSequenceLength();
			ofs << "Query: " << query_name << '\n';
			ofs << "Length: " << query_length << '\n' << '\n';

#ifdef DEBUG
			AlphabetCoder::Code *query_sequence = query->GetSequence();
			string query_data(query_length, ' ');
			parameters.coder.Decode(&query_sequence[0], query_length, &query_data[0]);
			vector<string> data = utils::split(query_data, '#');

			for(const auto& d : data){
				if(d.size() != 0) ofs << d << '\n';
			}
			ofs << '\n';
#endif

			uint64_t search_space = statistics.CalculateSearchSpace(
					query_length, parameters.db_total_length,
					parameters.db_total_sequences, gapped_karlin_parameters, alpha, beta);

			for (vector<Result>::iterator it = results[i].begin(); it != results[i].end(); ++it) {
				double evalue = Statistics::Nominal2EValue(it->score, search_space, gapped_karlin_parameters);
				if(evalue > parameters.evalue){
					break;
				}

				stringstream ss_i, ss_b, ss_e, ss_p, ss_g;
				ss_i << std::fixed << std::setprecision(1)
						 << (1. - (static_cast<float>(it->mismatches + it->gap_openings + it->gap_extensions)
														/ it->alignment_length)) * 100;
				ss_b << std::fixed << std::setprecision(1)
						 << Statistics::Nominal2Normalized(it->score, gapped_karlin_parameters);

				ss_e << std::scientific << std::setprecision(1)
						 << evalue;

				ss_p << std::fixed << std::setprecision(1)
						 << (static_cast<float>(it->positives) / it->alignment_length) * 100;

				ss_g << std::fixed << std::setprecision(1)
						 << (static_cast<float>(it->gap_openings + it->gap_extensions) / it->alignment_length) * 100;

				ofs << "Sbjct: " << it->subject_name << '\n' << '\n';
				ofs	<< "Score = " << ss_b.str() << " bits (" << it->score << "), Expect = " << ss_e.str() << '\n';
				ofs << "Identities = " << it->alignment_length - it->mismatches - it->gap_openings - it->gap_extensions
						<< "/" << it->alignment_length << " (" << ss_i.str() << "%), ";
				ofs << "Positives = " << it->positives << "/" << it->alignment_length
						<< " (" << ss_p.str() << "%), ";
				ofs << "Gaps = " << it->gap_openings + it->gap_extensions << "/" << it->alignment_length
						<< " (" << ss_g.str() << "%)" << '\n' << '\n';
				ofs << it->pairwise_result << '\n';
			}
		}
	}
}

int Aligner::AlignmentComp(const int a_score, const Coordinate &a_start, const Coordinate &a_end,
		const int b_score, const Coordinate &b_start, const Coordinate &b_end) {
	if (a_score > b_score) {
		return 1;
	} else if (a_score < b_score) {
		return -1;
	}

	uint32_t a_q_length = a_end.query_position - a_start.query_position;
	uint32_t b_q_length = b_end.query_position - b_start.query_position;

	if (a_q_length > b_q_length) {
		return 1;
	} else if (a_q_length < b_q_length) {
		return -1;
	}

	uint32_t a_s_length = a_end.database_position - a_start.database_position;
	uint32_t b_s_length = b_end.database_position - b_start.database_position;

	if (a_s_length > b_s_length) {
		return 1;
	} else if (a_s_length < b_s_length) {
		return -1;
	}

	return 0;
}

void Aligner::BuildResult(Query &query, DatabaseChunk &chunk, int score, Coordinate &hit, Coordinate &start,
		Coordinate &end, vector<EditBlocks::EditOpType> &edits, Result &result){
	AlphabetCoder::Code *query_sequence = query.GetSequence();
	AlphabetCoder::Code *database_sequence = chunk.GetConcatenatedSeqRef();
	uint32_t query_position = start.query_position;
	uint32_t database_position = start.database_position;
	uint32_t m = 0;
	uint32_t p = 0;
	uint32_t o = 0;
	uint32_t e = 0;
	EditBlocks::EditOpType prev_op = EditBlocks::kSubstitution;

	uint32_t subject_id = chunk.GetSubjectId(hit.database_position);
	result.subject_name = chunk.GetSubjectName(subject_id);

	if(parameters.outfmt == Outfmt::pairwise){
		uint32_t number_letters = parameters.score_matrix.GetNumberLetters();
		const int *sm = parameters.score_matrix.GetMatrix();

		uint32_t query_length = end.query_position - query_position + 1;
		uint32_t database_length = end.database_position - database_position + 1;
		std::string query_data(query_length, ' ');
		std::string database_data(database_length, ' ');
		parameters.coder.Decode(&query_sequence[query_position], query_length, &query_data[0]);
		parameters.coder.Decode(&database_sequence[database_position], database_length, &database_data[0]);

		stringstream ss_query, ss_match, ss_database;
		ss_query << std::left << setw(7) << "Query" << setw(4) << query.GetRealStart(query_position);
		ss_match << setw(11) << "";
		ss_database << std::left << setw(7) << "Sbjct" << setw(4) << database_position - (chunk.GetOffset(subject_id) - 1);
		uint32_t q=0, d=0;

		for (vector<EditBlocks::EditOpType>::iterator it = edits.begin(); it != edits.end(); ++it) {
			switch (*it) {
				case EditBlocks::kSubstitution:
					ss_query << query_data[q];
					ss_database << database_data[d];
					if (query_data[q] != database_data[d]) {
						++m;
						if (sm[query_data[q] * number_letters + database_data[d]] > 0){
							++p;
							ss_match << "+";
						}else{
							ss_match << " ";
						}
					}else{
						++p;
						ss_match << query_data[q];
					}
					++q;
					++d;
					break;
				case EditBlocks::kGapInSeq0:
					ss_query << "-";
					ss_match << " ";
					ss_database << database_data[d];
					if (prev_op != EditBlocks::kGapInSeq0) {
						++o;
					} else {
						++e;
					}
					++d;
					break;
				case EditBlocks::kGapInSeq1:
					ss_query << query_data[q];
					ss_match << " ";
					ss_database << "-";
					if (prev_op != EditBlocks::kGapInSeq1) {
						++o;
					} else {
						++e;
					}
					++q;
					break;
				default:
					abort();
					break;
			}
			prev_op = *it;
			if(ss_query.str().size() > 11+50){
				ss_query << "  " << query.GetRealEnd(query_position + q - 1);
				ss_database << "  " << database_position + d - chunk.GetOffset(subject_id);
				result.pairwise_result += ss_query.str() + "\n" + ss_match.str() + "\n" + ss_database.str() + "\n\n";
				ss_query.str("");
				ss_match.str("");
				ss_database.str("");
				ss_query.clear(stringstream::goodbit);
				ss_match.clear(stringstream::goodbit);
				ss_database.clear(stringstream::goodbit);
				ss_query << std::left << setw(7) << "Query" << setw(4) << query.GetRealStart(query_position + q);
				ss_match << setw(11) << "";
				ss_database << std::left << setw(7) << "Sbjct"
						<< setw(4) << database_position + d - (chunk.GetOffset(subject_id) - 1);
			}
		}
		ss_query << "  " << query.GetRealEnd(query_position + q - 1);
		ss_database << "  " << database_position + d - chunk.GetOffset(subject_id);
		result.pairwise_result += ss_query.str() + "\n" + ss_match.str() + "\n" + ss_database.str() + "\n";

	}else{
		for (vector<EditBlocks::EditOpType>::iterator it = edits.begin(); it != edits.end(); ++it) {
			switch (*it) {
				case EditBlocks::kSubstitution:
					if (query_sequence[query_position] != database_sequence[database_position]) {
						++m;
					}
					++query_position;
					++database_position;
					break;
				case EditBlocks::kGapInSeq0:
					if (prev_op != EditBlocks::kGapInSeq0) {
						++o;
					} else {
						++e;
					}
					++database_position;
					break;
				case EditBlocks::kGapInSeq1:
					if (prev_op != EditBlocks::kGapInSeq1) {
						++o;
					} else {
						++e;
					}
					++query_position;
					break;
				default:
					abort();
					break;
			}
			prev_op = *it;
		}
	}

	result.score = score;
	result.alignment_length = edits.size();
	result.mismatches = m;
	result.positives = p;
	result.gap_openings = o;
	result.gap_extensions = e;
	result.hit = hit;
	result.start.query_position = query.GetRealStart(start.query_position);
	result.end.query_position = query.GetRealEnd(end.query_position);
	result.start.database_position = start.database_position
			- (chunk.GetOffset(subject_id) - 1);
	result.end.database_position = end.database_position
			- (chunk.GetOffset(subject_id) - 1);
}
