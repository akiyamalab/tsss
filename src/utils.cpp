/*
	utils.cpp
	Author: kazuki takabatake
	Updated: 2019/12/14
*/

#include "utils.hpp"

#include <algorithm>
#include <sstream>
#include <sys/unistd.h>

using namespace std;

namespace utils{
  double Get_max_memory_consumption(){
    char command[128];
    char field[51];
    FILE *p_file;

    sprintf(command, "grep VmHWM /proc/%d/status | cut -c 7- | tr -d ' '", getpid());
    p_file = popen(command, "r");
    fgets(field, 50, p_file);
    pclose(p_file);
    return std::atof(field)/(1024*1024);
  }

	std::vector<std::string> split(const std::string &str, const char delim){
		std::vector<std::string> elems;
		std::string item;
		for(char ch: str){
			if(ch == delim){
				elems.push_back(item);
				item.clear();
			}else{
				item += ch;
			}
		}

		elems.push_back(item);
		return elems;
	}

	void WriteString(ofstream &ofs, const string &s){
		int size = s.size()+1;
		const char *s_char = s.c_str();

		ofs.write((char*)&size, sizeof(int));
		ofs.write((char*)s_char, sizeof(char)*size);
	}

	string LoadString(ifstream &ifs){
		int size;
		ifs.read((char*)&size, sizeof(int));
		char s_char[size];

		ifs.read((char*)&s_char, sizeof(char)*size);
		return string(s_char);
	}
}
