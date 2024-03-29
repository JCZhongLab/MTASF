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



#ifndef PROT_ZERO_PTM_FILTER_MNG_HPP_
#define PROT_ZERO_PTM_FILTER_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class ZeroPtmFilterMng {
 public:
  ZeroPtmFilterMng(PrsmParaPtr prsm_para_ptr,
                   int thread_num,
                   const std::string & output_file_ext):
      prsm_para_ptr_(prsm_para_ptr),
      thread_num_(thread_num),
      output_file_ext_(output_file_ext) {}

  PrsmParaPtr prsm_para_ptr_;

  /** parameters for fast filteration */
  int max_proteoform_mass_ = 100000;

  //Candidate protein number for each spectrum
  unsigned int comp_num_ = 5;
  unsigned int pref_suff_num_ = 5;
  unsigned int inte_num_ = 10;
  int filter_scale_ = 100;

  int thread_num_ = 1;

  std::string output_file_ext_;
};

typedef std::shared_ptr<ZeroPtmFilterMng> ZeroPtmFilterMngPtr;

} /* namespace tools */

#endif /* ZERO_PTM_FILTER_MNG_HPP_ */
