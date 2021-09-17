/*
	database_build.hpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#ifndef DATABASE_BUILD_HPP
#define DATABASE_BUILD_HPP

#include "database.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

class DatabaseBuild{
public:
	void Run(po::variables_map& vm) ;

private:
	void SetParameters(po::variables_map& vm, Database &database) const;
};

#endif
