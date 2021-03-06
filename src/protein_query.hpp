/*
 * protein_query.h
 *
 *  Created on: 2011/04/13
 *      Author: shu
 */

#ifndef NORMAL_QUERY_HPP
#define NORMAL_QUERY_HPP

#include "query.hpp"
#include "alphabet_coder.hpp"
#include "statistics.hpp"
#include "sequence_filter_interface.hpp"

#include <vector>

using namespace std;

class ProteinQuery: public Query {
public:
	static uint32_t GetSequenceLength(uint32_t length) {
		if (length == 0) {
			return 0;
		} else {
			return length + 2;
		}
	}
	ProteinQuery(const Sequence &sequence,
			AlphabetCoder::Code sequence_delimiter,
			SequenceFilterInterface *filter, ScoreMatrix score_matrix,
			Statistics::KarlinParameters &ungapped_ideal_karlin_parameters);

	virtual uint32_t GetRealSequenceLength() {
		return sequence_.size() - 2;
	}

	Statistics::KarlinParameters GetUngappedKarlinParameters() {
		return ungapped_karlin_parameters_;
	}

private:
	void SetName(const Sequence &sequence);
	void SetEncodedSequence(const Sequence &sequence,
			const AlphabetCoder &coder,  SequenceFilterInterface *filter);
	void SetStaticalParameters(ScoreMatrix &score_matrix,
			Statistics::KarlinParameters &ungapped_ideal_karlin_parameters);
};

#endif /* NORMAL_QUERY_HPP */
