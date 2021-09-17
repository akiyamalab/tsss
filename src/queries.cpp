/*
 * query.cpp
 *
 *  Created on: 2009/10/09
 *      Author: shu
 */

#include "fasta_sequence_reader.hpp"
#include "sequence.hpp"
#include "alphabet_coder.hpp"
#include "dna_sequence.hpp"
#include "dna_type.hpp"
#include "sequence_no_filter.hpp"
#include "sequence_codon_filter.hpp"
#include "protein_type.hpp"
#include "statistics.hpp"
#include "translator.hpp"
#include "translated_dna_query.hpp"
#include "protein_query.hpp"
#include "queries.hpp"
#include "logger.hpp"

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <numeric>
#include <typeinfo>
#include <string>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <functional>

using namespace std;

Queries::Queries(istream &in, const Parameters &parameters) :
		parameters_(parameters), max_query_sequence_length_(0), reader_(in), next_query_(
				QueryPtr()) {
	if (parameters_.filter) {
		if (typeid(*(parameters_.aligning_sequence_type_ptr))
				== typeid(ProteinType)) {
			sequence_filter_ = std::tr1::shared_ptr<SequenceFilterInterface>(
					new SequenceCodonFilter());
		} else {
			Logger *logger = Logger::GetInstance();
			logger->ErrorLog("can't use a filter for dna sequences");
			exit(1);
		}
	} else {
		sequence_filter_ = std::tr1::shared_ptr<SequenceFilterInterface>(
				new SequenceNoFilter());
	}
	Next();
}

void Queries::Next() {
  max_query_sequence_length_ = 0;
	SetQueries();
}

bool Queries::SetQueries() {
	queries_.clear();
	max_query_sequence_length_ = 0;
	number_of_skipped_sequences_ = 0;
	unsigned int chunk_size = 0;
	string name;
	string sequence;
	bool reader_ret;
	if (!next_query_) {
		while(1){
			reader_ret = reader_.Read(name, sequence);
			if (reader_ret) {
				Sequence s(name, sequence);
				if(sequence.size() >= parameters_.lower_limit_of_sequence_length){
					next_query_ = BuildQuery(s);
					break;
				}
				else{
					number_of_skipped_sequences_++;
				}
			}
			else{
				break;
			}
		}
	}

	while (1) {
		if (!next_query_) {
			break;
		}
		max_query_sequence_length_ = max(max_query_sequence_length_,
				next_query_->GetSequenceLength());
		unsigned int new_chunk_size = chunk_size
				+ next_query_->GetSequenceLength();
		if (new_chunk_size >= parameters_.chunk_size) {
			break;
		}
		queries_.push_back(next_query_);
		chunk_size = new_chunk_size;
		while(1){
			reader_ret = reader_.Read(name, sequence);
			if (reader_ret) {
				Sequence s(name, sequence);
				if(sequence.size() >= parameters_.lower_limit_of_sequence_length){
					next_query_ = BuildQuery(s);
					break;
				}
				else{
					number_of_skipped_sequences_++;
				}
			} 
			else {
				break;
			}
		}
		if(!reader_ret){
			next_query_ = QueryPtr();
			break;
		}
	}

	if(number_of_skipped_sequences_ > 0){
		cout << "number of queries skipped due to sequence length: " << number_of_skipped_sequences_ << endl;
	}

	return chunk_size;
}

Queries::QueryPtr Queries::BuildQuery(Sequence &sequence) {
	if (typeid(*(parameters_.file_sequence_type_ptr)) == typeid(DnaType)) {
		if (typeid(*(parameters_.aligning_sequence_type_ptr))
				!= typeid(ProteinType)) {
			Logger *logger = Logger::GetInstance();
			logger->ErrorLog("can't use aligner for query dna");
			exit(1);
		} else {
			return Queries::QueryPtr(
					new TranslatedDnaQuery(sequence,
							parameters_.sequence_delimiter, translator_,
							sequence_filter_.get(), parameters_.score_matrix,
							parameters_.ungapped_karlin_parameters));
		}
	} else {
		return Queries::QueryPtr(
				new ProteinQuery(sequence, parameters_.sequence_delimiter,
						sequence_filter_.get(), parameters_.score_matrix,
						parameters_.ungapped_karlin_parameters));
	}
	return QueryPtr();
}

