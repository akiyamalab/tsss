/*
	main.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <utility>

#include "database_build.hpp"
#include "alignment.hpp"
#include "utils.hpp"

#include <omp.h>

namespace po = boost::program_options;

namespace {
	po::options_description DefineDbOption(){
		po::options_description opt;
		po::options_description req_opt("(Required)");
		po::options_description opt_opt("(Optional)");

		req_opt.add_options()
			("database,i", po::value<std::string>(), "database file")
			("output,o", po::value<std::string>(), "output file");

		opt_opt.add_options()
			("chunk_size,c", po::value<uint32_t>()->default_value(1 << 30), "chunk size (bytes) [1073741824 (=1GB)]")
			("threads,t", po::value<int>()->default_value(1), "num of threads [1]")
			("seed1_length,l", po::value<int>()->default_value(3), "length of first seed [3]")
			("seed1_amino,a", po::value<int>()->default_value(14), "type of amino in first seed [14]")
			("seed2_length,L", po::value<int>()->default_value(5), "length of second seed [5]")
			("seed2_amino,A", po::value<int>()->default_value(8), "type of amino in second seed [8]");

		opt.add(req_opt).add(opt_opt);
		return opt;
	}

	po::options_description DefineAlnOption(){
		po::options_description opt;
		po::options_description req_opt("(Required)");
		po::options_description opt_opt("(Optional)");

		req_opt.add_options()
			("query,i", po::value<std::string>(), "query file")
			("database,d", po::value<std::string>(), "database file")
			("output,o", po::value<std::string>(), "output file");

		opt_opt.add_options()
			("query_type,q", po::value<char>()->default_value('p'), "query sequence type, p(protein) or d(dna) [p]")
			("query_filter,f", po::value<bool>()->default_value(1), "filter query sequence, 1(enable) or 0(disable) [1]")
      ("alignments,n", po::value<int>()->default_value(10), "number of outputs for each query [10]")
      ("threads,t", po::value<int>()->default_value(1), "number of threads [1]]")
			("seed1_hamming,h", po::value<int>()->default_value(0), "allowed hamming distance of first seed in seed search [0]")
			("seed2_hamming,H", po::value<int>()->default_value(1), "allowed hamming distance of second seed in seed search [1]")
			("outfmt,F", po::value<int>()->default_value(1), "output format, 0(pairwise) or 1(tabular) [1]")
			("evalue,e", po::value<double>()->default_value(10), "evalue threshold [10]")
			("chunk_size,c", po::value<uint32_t>()->default_value(1<<27), "query chunk size [128MB]");      

		opt.add(req_opt).add(opt_opt);
		return opt;
	}

	po::variables_map CommandParse(int argc, char *argv[], po::options_description opt){
		po::variables_map vm;
		try{
			po::store(po::parse_command_line(argc, argv, opt), vm);
		}catch(const po::error_with_option_name& e){
			std::cerr << e.what() << std::endl;
			exit(1);
		}
		po::notify(vm);
		return vm;
	}
}

int main(int argc, char *argv[]){

	double start, end;
	start = omp_get_wtime();

	po::options_description db_opt = DefineDbOption();
	po::options_description aln_opt = DefineAlnOption();

	if (argc > 1 && std::strcmp(argv[1], "db") == 0){
		po::variables_map vm = CommandParse(argc-1, argv+1, db_opt);

		if (!vm.count("output") || !vm.count("database")){
			std::cout << "usage: tsss db [<options>]\n"
								<< db_opt << std::endl;
			exit(0);
		}
		DatabaseBuild db;
		db.Run(vm);

	}else if (argc > 1 && std::strcmp(argv[1], "aln") == 0){
		po::variables_map vm = CommandParse(argc-1, argv+1, aln_opt);

		if (!vm.count("query") || !vm.count("database") || !vm.count("output")){
			std::cout << "usage: tsss aln [<options>]\n"
								<< aln_opt << std::endl;
			exit(0);
		}
		Alignment align;
		align.Run(vm);

	}else{
		std::cout << "********** TSSS **********\n"
							<< "There are two TSSS sub-commands:\n\n"
							<< "db    convert a Protein fasta file to a TSSS format database file\n"
							<< "usage: tsss db [options]\n"
							<< db_opt << "\n"
							<< "aln   search homologues of queries from database\n"
							<< "usage: tsss aln [options]\n"
							<< aln_opt << std::endl;
		exit(0);
	}

	end = omp_get_wtime();
	std::cout << std::fixed << std::setprecision(2) <<
			"real time: " << end - start << " sec." << std::endl;
	std::cout << "peak memory usage: " << utils::Get_max_memory_consumption() << " GB" << std::endl;

	return 0;
}
