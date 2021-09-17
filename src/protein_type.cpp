/*
 * protein.cpp
 *
 *  Created on: 2010/09/14
 *      Author: shu
 */

#include "protein_type.hpp"

using namespace std;

const string ProteinType::kRegularLetters = "ARNDCQEGHILKMFPSTWYV";
const string ProteinType::kAmbiguousLetters = "BZ*";
const char ProteinType::kUnknownLetter = 'X';
