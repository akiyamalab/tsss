/*
 * protein_sequence.h
 *
 *  Created on: 2011/04/11
 *      Author: shu
 */

#ifndef PROTEIN_SEQUENCE_HPP
#define PROTEIN_SEQUENCE_HPP

#include "sequence.hpp"
#include <string>

using namespace std;

class ProteinSequence : public Sequence{
public:
  ProteinSequence(const std::string &name, const std::string &sequence_data)
  : Sequence(name, sequence_data)
  {}

  virtual ~ProteinSequence()
  {}
};

#endif /* PROTEIN_SEQUENCE_HPP */
