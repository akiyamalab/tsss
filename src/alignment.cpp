/*
	alignment.hpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "alignment.hpp"
#include "aligner.hpp"
#include "alphabet_coder.hpp"
#include "protein_type.hpp"
#include "dna_type.hpp"
#include "score_matrix_reader.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <tr1/memory>

#include <omp.h>

using namespace std;

void Alignment::Run(po::variables_map& vm){
	Aligner aligner;
	SetParameters(vm, aligner);
	aligner.Align();

	cout << "finished" << endl << endl;
}

void Alignment::SetParameters(po::variables_map& vm, Aligner &aligner) const{
	const string default_matrix_name = "BLOSUM62";
	const string default_matrix ="   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

	ScoreMatrixReader score_matrix_reader;
	vector<int> matrix;
	unsigned int number_letters;
	ProteinType type;

	istringstream default_score_matrix_is(default_matrix);
	score_matrix_reader.Read(default_score_matrix_is, type, matrix, number_letters);
	aligner.parameters.score_matrix = ScoreMatrix(default_matrix_name, &matrix[0], number_letters);

	aligner.parameters.query_path = vm["query"].as<string>();
	aligner.parameters.database_path = vm["database"].as<string>();
	aligner.parameters.output_path = vm["output"].as<string>();
	aligner.parameters.chunk_size = vm["chunk_size"].as<uint32_t>();
	aligner.parameters.num_alignments = vm["alignments"].as<int>();
	aligner.parameters.seed1_hamming = vm["seed1_hamming"].as<int>();
	aligner.parameters.seed2_hamming = vm["seed2_hamming"].as<int>();

	int threads = vm["threads"].as<int>();
	if(threads > 0){

#ifdef _OPENMP
		int max_threads = omp_get_max_threads();

		if(threads < max_threads){
			aligner.parameters.threads = threads;
		}else{
			aligner.parameters.threads = max_threads;
		}
#else
		aligner.parameters.threads = 1;
#endif

	}else{
		cerr << "[ERROR] please set number of threads > 0" << endl;
		exit(1);
	}

	Aligner::Outfmt fmt = static_cast<Aligner::Outfmt>(vm["outfmt"].as<int>());
	if(fmt == Aligner::Outfmt::pairwise || fmt == Aligner::Outfmt::tabular){
		aligner.parameters.outfmt = fmt;
	}else{
		cerr << "[ERROR] undefined output format: " << fmt << endl;
		exit(1);
	}

	cout << "[query path]    " << aligner.parameters.query_path << endl;
	cout << "[output path]   " << aligner.parameters.output_path << endl;
	cout << "[database path] " << aligner.parameters.database_path << endl << endl;
	cout << "[CPU threads]   " << aligner.parameters.threads << endl << endl;

	AlphabetCoder coder(type);
	aligner.parameters.coder = coder;
	aligner.parameters.delimiter = coder.GetMaxCode()+1;

	aligner.parameters.db_type = std::tr1::shared_ptr<AlphabetType>(new ProteinType());

	char query_type = vm["query_type"].as<char>();
	if(query_type == 'p'){
		aligner.parameters.query_type = std::tr1::shared_ptr<AlphabetType>(new ProteinType());
	}else if(query_type == 'd'){
		aligner.parameters.query_type = std::tr1::shared_ptr<AlphabetType>(new DnaType());
	}else{
		cerr << "[ERROR] undefined sequence type: " << query_type << endl;
		exit(1);
	}

	aligner.parameters.filter = vm["query_filter"].as<bool>();
	aligner.parameters.evalue = vm["evalue"].as<double>();

	aligner.parameters.gap_open = 11;
	aligner.parameters.gap_extension = 1;
	aligner.parameters.normalized_presearched_ungapped_extension_cutoff = 7.0;
	aligner.parameters.normalized_presearched_gapped_extension_trigger = 22.0;
	aligner.parameters.normalized_presearched_gapped_extension_cutoff = 15.0;
	aligner.parameters.normalized_result_gapped_extension_cutoff = 25.0;
}
