/*
	database_build.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "database_build.hpp"
#include "score_matrix_reader.hpp"
#include "protein_type.hpp"

#include <iostream>
#include <string>
#include <climits>
#include <cmath>

#include <omp.h>

using namespace std;

void DatabaseBuild::Run(po::variables_map& vm){
	Database database;
	SetParameters(vm, database);
	database.BuildDatabase();

	cout << "finished" << endl << endl;
}

void DatabaseBuild::SetParameters(po::variables_map& vm, Database &database) const{

	const string default_matrix_name = "BLOSUM62";
	const string default_matrix ="   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

	ScoreMatrixReader score_matrix_reader;
	vector<int> matrix;
	unsigned int number_letters;
	ProteinType type;

	istringstream default_score_matrix_is(default_matrix);
	score_matrix_reader.Read(default_score_matrix_is, type, matrix, number_letters);
	database.parameters.score_matrix = ScoreMatrix(default_matrix_name, &matrix[0], number_letters);

	database.parameters.database_path = vm["database"].as<string>();
	database.parameters.output_path = vm["output"].as<string>();
	database.parameters.chunk_size = vm["chunk_size"].as<uint32_t>();


	int threads = vm["threads"].as<int>();
	if(threads > 0){

#ifdef _OPENMP
		int max_threads = omp_get_max_threads();

		if(threads < max_threads){
			database.parameters.threads = threads;
		}else{
			database.parameters.threads = max_threads;
		}
#else
		database.parameters.threads = 1;
#endif

	}else{
		cerr << "[ERROR] please set number of threads > 0" << endl;
		exit(1);
	}

	cout << "[database path] " << database.parameters.database_path << endl;
	cout << "[output path]   " << database.parameters.output_path << endl << endl;
	cout << "[CPU threads]   " << database.parameters.threads << endl << endl;

	int s1_l = vm["seed1_length"].as<int>();
	int s1_a = vm["seed1_amino"].as<int>();
	int s2_l = vm["seed2_length"].as<int>();
	int s2_a = vm["seed2_amino"].as<int>();

	if(pow(s1_a, s1_l) > UINT_MAX){
		cerr << "[ERROR] please set (seed1_amino ^ seed1_length) <= UINT_MAX \n"
				 << "current seed1_length is " << s1_l << ", seed1_amino is " << s1_a << endl;
		exit(1);
	}

	if(pow(s2_a, s2_l) > USHRT_MAX){
		cerr << "[ERROR] please set (seed2_amino ^ seed2_length) <= USHRT_MAX \n"
				 << "current seed2_length is " << s2_l << ", seed2_amino is " << s2_a << endl;
		exit(1);
	}

	database.parameters.seed1_length = s1_l;
	database.parameters.seed1_amino = s1_a;
	database.parameters.seed2_length = s2_l;
	database.parameters.seed2_amino = s2_a;
}
