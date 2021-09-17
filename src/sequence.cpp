/*
 * sequence.cpp
 *
 *  Created on: 2009/06/19
 *      Author: shu
 */

#include "sequence.hpp"

using namespace std;

Sequence::Sequence(const string &name, const string &sequence_data)
	:sequence_data_(sequence_data){
	auto pos = name.find(" ");
	if(pos != string::npos){
		name_ = name.substr(0, pos);
	}else{
		name_ = name;
	}
}
