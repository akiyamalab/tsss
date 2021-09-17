/*
 * score_matrix_reader.h
 *
 *  Created on: 2010/06/22
 *      Author: shu
 */

#ifndef SCORE_MATRIX_READER_HPP
#define SCORE_MATRIX_READER_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include "alphabet_coder.hpp"
#include "alphabet_type.hpp"


class ScoreMatrix;

class ScoreMatrixReader {
public:
  void Read(std::istream &in, AlphabetType &type, std::vector<int> &matrix, unsigned int &number_letters);

private:
  std::list<std::string> Split(std::string str, std::string delim);
  int GetNumberLetters(const AlphabetCoder &coder);
  int GetMatrixSize(const AlphabetCoder &coder);
  void Parse(std::istream &in, const AlphabetCoder &coder, int *matrix);
};

#endif /* SCORE_MATRIX_READER_HPP */
