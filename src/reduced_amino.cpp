/*
	reduced_amino.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "reduced_amino.hpp"

#include <iostream>

using namespace std;

string ReducedAmino::GetSimilaritySets(const int num) const{
	string sets;

	switch(num){
		case 22:
			sets = "A S T R K Q E Z N D B H G P C L I V M F Y W";
			break;
		case 21:
			sets = "A S T R K Q EZ N D B H G P C L I V M F Y W";
			break;
		case 20:
			sets = "A S T R K Q EZ N DB H G P C L I V M F Y W";
			break;
		case 19:
			sets = "A S T R K Q EZ N DB H G P C L IV M F Y W";
			break;
		case 18:
			sets = "A S T R K QEZ N DB H G P C L IV M F Y W";
			break;
		case 17:
			sets = "A S T R K QEZ N DB H G P C LIV M F Y W";
			break;
		case 16:
			sets = "A S T RK QEZ N DB H G P C LIV M F Y W";
			break;
		case 15:
			sets = "A S T RK QEZ N DB H G P C LIV M FY W";
			break;
		case 14:
			sets = "A S T RK QEZ NDB H G P C LIV M FY W";
			break;
		case 13:
			sets = "A S T RK QEZ NDB H G P C LIVM FY W";
			break;
		case 12:
			sets = "A S T RKQEZ NDB H G P C LIVM FY W";
			break;
		case 11:
			sets = "A ST RKQEZ NDB H G P C LIVM FY W";
			break;
		case 10:
			sets = "A ST RKQEZ NDB H G P C LIVM FYW";
			break;
		case 9:
			sets = "AST RKQEZ NDB H G P C LIVM FYW";
			break;
		case 8:
			sets = "AST RKQEZNDB H G P C LIVM FYW";
			break;
		case 7:
			sets = "AST RKQEZNDBH G P C LIVM FYW";
			break;
		case 6:
			sets = "AST RKQEZNDBH G P CLIVM FYW";
			break;
		case 5:
			sets = "AST RKQEZNDBHG P CLIVM FYW";
			break;
		case 4:
			sets = "ASTRKQEZNDBHG P CLIVM FYW";
			break;
		case 3:
			sets = "ASTRKQEZNDBHGP CLIVM FYW";
			break;
		case 2:
			sets = "ASTRKQEZNDBHGP CLIVMFYW";
			break;
		default:
			std::cerr << "[ERROR] reduced amino type should be from 2 to 22" << endl;
			exit(1);
	}

	return sets;
}
