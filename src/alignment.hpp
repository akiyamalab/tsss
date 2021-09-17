/*
	alignment.hpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include "aligner.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

class Alignment{
public:
	void Run(po::variables_map& vm);

private:
	void SetParameters(po::variables_map& vm, Aligner &aligner) const;
};

#endif
