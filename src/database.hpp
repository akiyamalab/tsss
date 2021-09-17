/*
	database.hpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#ifndef DATABASE_HPP
#define DATABASE_HPP

#include "score_matrix.hpp"

#include <string>

class Database{
public:
	void BuildDatabase();

	struct Parameters{
		ScoreMatrix score_matrix;
		std::string database_path;
		std::string output_path;
		uint32_t chunk_size;
		int threads;
		int seed1_length;
		int seed1_amino;
		int seed2_length;
		int seed2_amino;
	}parameters;

private:
};

#endif
