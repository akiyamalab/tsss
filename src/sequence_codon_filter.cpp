/*
 * sequence_codon_filter.cpp
 *
 *  Created on: 2020/3/11
 *      Author: Motohiro Akikawa
 */

#include "sequence_codon_filter.hpp"
#include <iostream>
#include <stdio.h>
int SequenceCodonFilter::Filter(const std::string &seq, std::string *masked_seq) {
  unsigned int start=0;
  unsigned int end=0;
  unsigned int i = 0;
    *masked_seq=seq; 

  while(start<=seq.length()){
    i = start;
    while(seq[i]!='*' && i<=seq.length()){
       i++;

      }
    if(i-start<20){
      end=i;
      for(unsigned int j = start; j<end;j++){
	(*masked_seq)[j]='X';
      }
    }

    start=i+1;
  }
  //for debug
  //std::cout << seq << std::endl;
  //std::cout << *masked_seq << std::endl;
  //
  return 0;
}
