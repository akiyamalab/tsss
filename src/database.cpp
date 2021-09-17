/*
	database.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "database.hpp"
#include "database_chunk.hpp"

#include <iomanip>
#include <iostream>
#include <fstream>

#include <omp.h>

using namespace std;

void Database::BuildDatabase(){
	ifstream in(parameters.database_path);

	DatabaseChunk::Parameters parameters_;
	parameters_.score_matrix  = parameters.score_matrix;
	parameters_.database_path = parameters.database_path;
	parameters_.output_path   = parameters.output_path;
	parameters_.chunk_size    = parameters.chunk_size;
	parameters_.seed1_length  = parameters.seed1_length;
	parameters_.seed2_length  = parameters.seed2_length;
	parameters_.seed1_amino   = parameters.seed1_amino;
	parameters_.seed2_amino   = parameters.seed2_amino;
	parameters_.threads       = parameters.threads;

	double omp_s, omp_e;
	omp_s = omp_get_wtime();

	cout << std::fixed << std::setprecision(2);

	DatabaseChunk chunk(in, parameters_);
	for(; chunk.ReadChunk(); ){

		omp_e = omp_get_wtime();
		cout << "read   database chunk\t\t" << "[time] " << omp_e - omp_s << " sec." << endl;
		omp_s = omp_get_wtime();

		chunk.EncodeChunk();

		omp_e = omp_get_wtime();
		cout << "encode database chunk\t\t" << "[time] " << omp_e - omp_s << " sec." << endl;
		omp_s = omp_get_wtime();

		chunk.BuildChunk();

		omp_e = omp_get_wtime();
		cout << "build  database chunk\t\t" << "[time] " << omp_e - omp_s << " sec." << endl;
		omp_s = omp_get_wtime();

		chunk.SaveChunk();

		omp_e = omp_get_wtime();
		cout << "write  database chunk\t\t" << "[time] " << omp_e - omp_s << " sec." << endl;
		omp_s = omp_get_wtime();

		chunk.ClearChunk();

		omp_e = omp_get_wtime();
		cout << "clear  database chunk\t\t" << "[time] " << omp_e - omp_s << " sec.\n" << endl;
		omp_s = omp_get_wtime();

	}

	chunk.SaveInfo();
}
