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


#ifndef PROT_ONE_PTM_SEARCH_HPP_
#define PROT_ONE_PTM_SEARCH_HPP_

#include "spec/spectrum_set.hpp"
#include "prsm/prsm.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"

namespace prot {

class OnePtmSearchProcessor {
 public:
  explicit OnePtmSearchProcessor(PtmSearchMngPtr mng_ptr): mng_ptr_(mng_ptr) {}
  void process();

 private:
  PrsmPtrVec onePtmSearchOneSpec(SpectrumSetPtr spec_set_ptr,
                                 const SimplePrsmPtrVec &simple_prsm_ptr_vec,
                                 FastaIndexReaderPtr reader_ptr,
                                 PtmSearchMngPtr mng_ptr,
                                 AlignTypePtr type_ptr);

  PtmSearchMngPtr mng_ptr_;
};

typedef std::shared_ptr<OnePtmSearchProcessor> OnePtmSearchProcessorPtr;

}  // namespace prot

#endif
