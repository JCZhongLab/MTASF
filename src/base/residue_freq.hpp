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


#ifndef PROT_BASE_RESIDUE_FREQ_HPP_
#define PROT_BASE_RESIDUE_FREQ_HPP_

#include <vector>

#include "base/residue.hpp"

namespace prot {

class ResidueFreq: public Residue {
 public:
  ResidueFreq(AminoAcidPtr acid_ptr, PtmPtr ptm_ptr, double freq):
      Residue(acid_ptr, ptm_ptr),
      freq_(freq) {}

  double getFreq() {return freq_;}

 private:
  double freq_;
};

typedef std::shared_ptr<ResidueFreq> ResFreqPtr;
typedef std::vector<ResFreqPtr> ResFreqPtrVec;

}  // namespace prot
#endif