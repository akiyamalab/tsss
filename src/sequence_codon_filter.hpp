/*
 * sequence_no_filter.h
 *
 *  Created on: 2020/3/11
 *      Author: Motohiro Akikawa
 */

#ifndef SEQUENCE_CODON_FILTER_HPP
#define SEQUENCE_CODON_FILTER_HPP

#include "sequence_filter_interface.hpp"

class SequenceCodonFilter : public SequenceFilterInterface{
public:
	SequenceCodonFilter() {
	}

	~SequenceCodonFilter() {

	}
	int Filter(const std::string &seq, std::string *masked_seq);
};

#endif /* SEQUENCE_Codon_FILTER_H_ */
