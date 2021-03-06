/*
 * queries.h
 *
 *  Created on: 2009/10/09
 *      Author: shu
 */

#ifndef QUERIES_HPP
#define QUERIES_HPP

#include "sequence.hpp"
#include "fasta_sequence_reader.hpp"
#include "sequence_filter_interface.hpp"
#include "query.hpp"
#include "statistics.hpp"
#include "translator.hpp"

#include <vector>
#include <string>
#include <stdint.h>
#include <tr1/memory>

class Queries {
public:
  typedef struct {
	  bool filter;
    std::tr1::shared_ptr<AlphabetType> file_sequence_type_ptr;
    std::tr1::shared_ptr<AlphabetType> aligning_sequence_type_ptr;
    Statistics::KarlinParameters ungapped_karlin_parameters;
    unsigned int chunk_size;
    ScoreMatrix score_matrix;
    AlphabetCoder::Code sequence_delimiter;
    uint32_t lower_limit_of_sequence_length;
  } Parameters;

  Queries(std::istream &in, const Parameters &parameters);

  virtual ~Queries()
  {}

  void Next();

  Query *GetQuery(uint32_t id) {
    return queries_[id].get();
  }
  uint32_t GetNumberSequences() {
    return queries_.size();
  }

  uint32_t GetMaxQuerySequenceLength() {
    return max_query_sequence_length_;
  }

  uint32_t GetNumberSkippedSequences() {
    return number_of_skipped_sequences_;
  }

private:
  typedef std::tr1::shared_ptr<Query> QueryPtr;
  bool SetQueries();
  QueryPtr BuildQuery(Sequence &sequence);

  Parameters parameters_;
  std::vector<QueryPtr> queries_;
  uint32_t max_query_sequence_length_;
  uint32_t number_of_skipped_sequences_;
  Translator translator_;
  std::tr1::shared_ptr<SequenceFilterInterface> sequence_filter_;
  FastaSequenceReader reader_;
  QueryPtr next_query_;
};

#endif /* QUERIES_HPP */
