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


#ifndef PROT_BASE_MOD_BASE_HPP_
#define PROT_BASE_MOD_BASE_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/mod.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class ModBase {
 public:
  static void initBase(const std::string &file_name);

  static const ModPtrVec& getBaseModPtrVec() {return mod_ptr_vec_;}

  static ModPtr getNoneModPtr() {return none_mod_ptr_;}

  static ModPtr getC57ModPtr() {return c57_mod_ptr_;}

  static ModPtr getC58ModPtr() {return c58_mod_ptr_;}

  static ModPtr getBaseModPtr(ModPtr mod_ptr);

  static ModPtr getBaseModPtr(ResiduePtr ori_residue, ResiduePtr mod_residue);

  static bool isNoneModPtr(ModPtr mod_ptr) {return mod_ptr == none_mod_ptr_;}

  static ModPtr getModPtrFromXml(xercesc::DOMElement * element);

 private:
  static ModPtrVec mod_ptr_vec_;
  static ModPtr none_mod_ptr_;
  //C57
  static ModPtr c57_mod_ptr_;
  //C58
  static ModPtr c58_mod_ptr_;
};

}

#endif

