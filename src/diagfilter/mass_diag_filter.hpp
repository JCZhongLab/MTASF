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

#ifndef PROT_MASS_DIAG_FILTER_H_
#define PROT_MASS_DIAG_FILTER_H_

#include <vector>

#include "base/proteoform.hpp"
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/mass_match.hpp"
#include "diagfilter/diag_filter_mng.hpp"

namespace prot {

class MassDiagFilter {
 public:
  MassDiagFilter(const ProteoformPtrVec &proteo_ptrs, DiagFilterMngPtr mng_ptr);

  SimplePrsmPtrVec getBestMatch(const PrmMsPtrVec &ms_ptr_vec);

 private:
  DiagFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  MassMatchPtr index_ptr_;

  SimplePrsmPtrVec compute(const PrmMsPtrVec &ms_ptr_vec);
};

typedef std::shared_ptr<MassDiagFilter> MassDiagFilterPtr;
} /* namespace prot */

#endif /* PROT_DIAG_FILTER_H_ */
