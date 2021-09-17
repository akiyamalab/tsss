/*
	database_chunk.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "aligner.hpp"
#include "database_chunk.hpp"
#include "protein_type.hpp"
#include "alphabet_coder.hpp"
#include "reduced_amino.hpp"
#include "utils.hpp"

#include <iostream>
#include <algorithm>

using namespace std;

DatabaseChunk::DatabaseChunk(std::ifstream &in, const Parameters &parameters_)
	:parameters(parameters_),reader(in),total_length(0),
	chunk_id(0),db_total_length(0),db_total_sequences(0)
{
	ProteinType type;
	coder = AlphabetCoder(type);
	parameters.delimiter = coder.GetMaxCode() + 1;
	searcher = SeedSearcher(parameters.score_matrix, parameters.delimiter, parameters.threads);

	ReducedAmino reduced_amino;

	seed1_coder = ReducedAlphabetCoder(reduced_amino.GetSimilaritySets(parameters.seed1_amino), parameters.seed1_length, 0);
	seed2_coder = ReducedAlphabetCoder(reduced_amino.GetSimilaritySets(parameters.seed2_amino), parameters.seed2_length, 0);

	cout << "[seed 1] ";
	seed1_coder.Profile();
	cout << endl << "[seed 2] ";
	seed2_coder.Profile();
	cout << endl << endl;
}

DatabaseChunk::DatabaseChunk(const Parameters &parameters_, uint64_t &db_total_length_, uint64_t &db_total_sequences_)
	:parameters(parameters_),searcher(parameters_.score_matrix, parameters_.delimiter, 0)
{
	string current_db_path = parameters.database_path;
	ifstream ifs;
	utils::FileOpen(ifs, current_db_path, ios::binary);
	ifs.read((char*)&parameters.seed1_length, sizeof(int));
	ifs.read((char*)&parameters.seed2_length, sizeof(int));
	ifs.read((char*)&parameters.seed1_amino, sizeof(int));
	ifs.read((char*)&parameters.seed2_amino, sizeof(int));
	ifs.read((char*)&num_db_chunk, sizeof(int));
	ifs.read((char*)&db_total_length_, sizeof(uint64_t));
	ifs.read((char*)&db_total_sequences_, sizeof(uint64_t));
	ifs.close();

	ReducedAmino reduced_amino;

	seed1_coder = ReducedAlphabetCoder(reduced_amino.GetSimilaritySets(parameters.seed1_amino),
			parameters.seed1_length, parameters.seed1_hamming);
	seed2_coder = ReducedAlphabetCoder(reduced_amino.GetSimilaritySets(parameters.seed2_amino),
			parameters.seed2_length, parameters.seed2_hamming);

	cout << "[seed 1] ";
	seed1_coder.Profile();
	cout << ", hamming distance: " << parameters.seed1_hamming << endl;
	cout << "[seed 2] ";
	seed2_coder.Profile();
	cout << ", hamming distance: " << parameters.seed2_hamming << endl << endl;
}

bool DatabaseChunk::ReadChunk(){
	string name,sequence;	
	concatenated_seq_length = 1;

	while(true){
		bool check = reader.Read(name, sequence);

		if(check){
			sequences.emplace_back(Sequence(name, sequence));
			total_length += sequence.length();
			concatenated_seq_length += sequence.length() + 1;

			if(total_length >= parameters.chunk_size){
				cout << "read chunk " << chunk_id << endl;
				return true;
			}
		}else{
			if(sequences.empty()){
				return false;
			}else{
				cout << "read chunk " << chunk_id << endl;
				return true;
			}
		}
	}
	return false;
}

void DatabaseChunk::EncodeChunk(){
	concatenated_sequence.resize(concatenated_seq_length);
	concatenated_seq_seed1.resize(concatenated_seq_length);
	concatenated_seq_seed2.resize(concatenated_seq_length);

	uint32_t offset = 1;
	offsets.resize(sequences.size()+1);
	concatenated_sequence[0] = parameters.delimiter;

	for(unsigned int i=0; i<sequences.size(); ++i){
		offsets[i] = offset;
		sequence_names.emplace_back(sequences[i].GetName());
		string sequence = sequences[i].GetSequenceData();
		coder.Encode(&sequence[0], sequence.length(), &concatenated_sequence[offset]);
		offset += sequence.length();
		concatenated_sequence[offset] = parameters.delimiter;
		offset++;
	}
	offsets[sequences.size()] = offset;

	db_total_sequences += sequences.size();
	db_total_length += total_length;

	//concatenated_seq_seed1, seed2 をencode
	seed1_coder.Encode(&concatenated_sequence[0], concatenated_seq_length, &concatenated_seq_seed1[0]);
	seed2_coder.Encode(&concatenated_sequence[0], concatenated_seq_length, &concatenated_seq_seed2[0]);

	sequences.clear();
	sequences.shrink_to_fit();
}

void DatabaseChunk::BuildChunk(){

	if(parameters.threads == 1){
		searcher.Build(seed1_coder, seed2_coder, concatenated_seq_length, parameters.seed1_length,
				parameters.seed2_length, &concatenated_sequence[0], &concatenated_seq_seed1[0], &concatenated_seq_seed2[0]);
	}else{
		searcher.BuildParallel(seed1_coder, seed2_coder, concatenated_seq_length, parameters.seed1_length,
				parameters.seed2_length, &concatenated_sequence[0], &concatenated_seq_seed1[0], &concatenated_seq_seed2[0]);
	}

	concatenated_seq_seed1.clear();
	concatenated_seq_seed2.clear();
	concatenated_seq_seed1.shrink_to_fit();
	concatenated_seq_seed2.shrink_to_fit();
}

void DatabaseChunk::SaveInfo() const{
	string current_output_path = parameters.output_path;
	ofstream ofs;
	utils::FileOpen(ofs, current_output_path, ios::binary);
	ofs.write((char*)&parameters.seed1_length, sizeof(int));
	ofs.write((char*)&parameters.seed2_length, sizeof(int));
	ofs.write((char*)&parameters.seed1_amino, sizeof(int));
	ofs.write((char*)&parameters.seed2_amino, sizeof(int));
	ofs.write((char*)&chunk_id, sizeof(int));
	ofs.write((char*)&db_total_length, sizeof(uint64_t));
	ofs.write((char*)&db_total_sequences, sizeof(uint64_t));
	ofs.close();
}

void DatabaseChunk::SaveChunk(){

	string current_output_path = parameters.output_path + "." + to_string(chunk_id);
	ofstream ofs;
	utils::FileOpen(ofs, current_output_path, ios::binary);

	//seq names
	int num = sequence_names.size();
	ofs.write((char*)&num, sizeof(int));
	for(const auto& name : sequence_names){
		utils::WriteString(ofs, name);
	}
	sequence_names.clear();
	sequence_names.shrink_to_fit();

	//seq data
	ofs.write((char*)&concatenated_seq_length, sizeof(int));
	ofs.write((char*)&concatenated_sequence[0], sizeof(AlphabetCoder::Code)*concatenated_seq_length);
	concatenated_sequence.clear();
	concatenated_sequence.shrink_to_fit();

	//offsets
	num = offsets.size();
	ofs.write((char*)&num, sizeof(int));
	ofs.write((char*)&offsets[0], sizeof(uint32_t)*num);
	offsets.clear();
	offsets.shrink_to_fit();

	//searcher
	searcher.Write(ofs);

	ofs.close();
}

bool DatabaseChunk::LoadChunk(){
	if(chunk_id >= num_db_chunk) return false;

	string current_db_path = parameters.database_path + "." + to_string(chunk_id);
	ifstream ifs;
	utils::FileOpen(ifs, current_db_path, ios::binary);

	//seq names
	int num;
	ifs.read((char*)&num, sizeof(int));
	sequence_names.resize(num);
	for(int i=0; i<num; ++i){
		sequence_names[i] = utils::LoadString(ifs);
	}

	//seq data
	ifs.read((char*)&concatenated_seq_length, sizeof(int));
	concatenated_sequence.resize(concatenated_seq_length);
	ifs.read((char*)&concatenated_sequence[0], sizeof(AlphabetCoder::Code)*concatenated_seq_length);

	//offsets
	ifs.read((char*)&num, sizeof(int));
	offsets.resize(num);
	ifs.read((char*)&offsets[0], sizeof(uint32_t)*num);

	//searcher
	uint32_t starts_size_check = seed1_coder.GetMaxSeedIndex();
	searcher.Load(ifs, starts_size_check);

	ifs.close();
	return true;
}

vector<SeedSearcher::Hit> DatabaseChunk::SeedSearch(const AlphabetCoder::Code *query,
		const uint32_t query_length, const int ungapped_cutoff, const int gapped_trigger) const{
	vector<AlphabetCoder::Code> encoded_query_seed1(query_length);
	vector<AlphabetCoder::Code> encoded_query_seed2(query_length);

	seed1_coder.Encode(&query[0], query_length, &encoded_query_seed1[0]);
	seed2_coder.Encode(&query[0], query_length, &encoded_query_seed2[0]);

	int pass = 0;
	for(int i=1; i<(2*parameters.seed2_length+parameters.seed1_length); ++i){
		if(encoded_query_seed1[i] == parameters.delimiter){
			pass = parameters.seed2_length + parameters.seed2_length + parameters.seed1_length;
		}

		if(pass){
			--pass;
		}
	}

	vector<SeedSearcher::Hit> hits;
	for(uint32_t i=parameters.seed2_length+1; i<query_length-parameters.seed2_length-parameters.seed1_length; ++i){

		if(encoded_query_seed1[i+parameters.seed1_length+parameters.seed2_length-1] == parameters.delimiter){
			pass = parameters.seed2_length + parameters.seed2_length + parameters.seed1_length;
		}
		if(pass){
			--pass;
			continue; //seed+narrow内にdelimiterがある場合
		}

		vector<uint32_t> indexes = seed1_coder.GetNeighborIndexes(&encoded_query_seed1[i]);
		vector<uint16_t> front_indexes = seed2_coder.GetNeighborIndexes16(&encoded_query_seed2[i-parameters.seed2_length]);
		vector<uint16_t> back_indexes = seed2_coder.GetNeighborIndexes16(&encoded_query_seed2[i+parameters.seed1_length]);

		std::sort(back_indexes.begin(), back_indexes.end());

		vector<SeedSearcher::Hit> new_hits = searcher.Search(indexes, front_indexes, back_indexes,
				 &query[0], i, &concatenated_sequence[0], ungapped_cutoff, gapped_trigger);
		std::copy(new_hits.begin(), new_hits.end(), std::back_inserter(hits));
	}

	std::sort(hits.begin(), hits.end());

	return hits;
}

void DatabaseChunk::ClearChunk(){
	searcher.Clear();
	sequences.clear();
	sequence_names.clear();
	concatenated_sequence.clear();
	concatenated_seq_seed1.clear();
	concatenated_seq_seed2.clear();
	offsets.clear();
	total_length = 0;
	concatenated_seq_length = 0;
	chunk_id++;
}

uint32_t DatabaseChunk::GetSubjectId(const uint32_t position) const{
	return distance(offsets.begin(), lower_bound(offsets.begin(), offsets.end(), position)) -1;
}
