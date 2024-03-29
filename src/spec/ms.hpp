//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef PROT_SPEC_MS_HPP_
#define PROT_SPEC_MS_HPP_

#include <sstream>

#include "spec/ms_header.hpp"

namespace prot {

template <class T>
class Ms {
 public:
  Ms() {};

  Ms(MsHeaderPtr header_ptr) {header_ptr_ = header_ptr;}

  Ms(MsHeaderPtr header_ptr, const std::vector<T> &peak_ptr_list) {
    header_ptr_ = header_ptr;
    peak_ptr_list_ = peak_ptr_list;
  }

  std::string toString() {
    std::string header_str = header_ptr_->toString();
    std::stringstream tmp;
    for (size_t i = 0; i < peak_ptr_list_.size(); i++) {
      tmp << i << " " << peak_ptr_list_[i]->getPosition() 
          << " " << peak_ptr_list_[i]->getIntensity() << "\n";
    }
    return header_str + tmp.str();
  }

  MsHeaderPtr getMsHeaderPtr() {return header_ptr_;}

  void setHeaderPtr(MsHeaderPtr header_ptr) {header_ptr_ = header_ptr;}

  size_t size() {return peak_ptr_list_.size();}

  T getPeakPtr(int i) {return peak_ptr_list_[i];}
  
  std::vector<T> getPeakPtrVec() {return peak_ptr_list_;}

 private:
  MsHeaderPtr header_ptr_;
  std::vector<T> peak_ptr_list_;
};

}
#endif
