/*
 * translator.h
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#ifndef TRANSLATOR_HPP
#define TRANSLATOR_HPP

#include <string>
#include <vector>
#include <map>
#include <stdint.h>

class DnaSequence;
class ProteinSequence;


class Translator {
public:
  Translator();
  void Translate(const DnaSequence &dna, std::vector<ProteinSequence> &proteins);

private:
  std::vector<std::string> start_codons_;
  std::map<std::string, std::string>codon_table_;
};

#endif /* TRANSLATOR_HPP */
