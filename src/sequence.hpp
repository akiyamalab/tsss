/*
 * sequence.h
 *
 *  Created on: 2009/06/19
 *      Author: shu
 */

#ifndef SEQUENCE_HPP
#define SEQUENCE_HPP

#include <string>

class Sequence {
public:
  Sequence(const std::string &name, const std::string &sequence_data);

  virtual ~Sequence()
  {}

  std::string GetName() const {
    return name_;
  }

  std::string GetSequenceData() const {
    return sequence_data_;
  }




private:
  std::string name_;
  std::string sequence_data_;

};

#endif /* SEQUENCE_HPP */
