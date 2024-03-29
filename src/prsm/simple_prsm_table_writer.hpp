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

#ifndef PROT_PRSM_SIMPLE_TABLE_WRITER_HPP_
#define PROT_PRSM_SIMPLE_TABLE_WRITER_HPP_

#include <string>

#include "prsm/prsm_para.hpp"
#include "prsm/simple_prsm.hpp"

namespace prot {

class SimplePrsmTableWriter {
 public:
  SimplePrsmTableWriter(PrsmParaPtr prsm_para_ptr,
                        const std::string &input_file_ext,
                        const std::string &output_file_ext):
    prsm_para_ptr_(prsm_para_ptr),
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext) {}

  void write();

 private:
  PrsmParaPtr prsm_para_ptr_;
  std::string input_file_ext_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<SimplePrsmTableWriter> SimplePrsmTableWriterPtr;

}  // namespace prot

#endif /* SIMPLE_TABLE_WRITER_HPP_ */
