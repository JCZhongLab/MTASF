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


#ifndef PROT_BASE_MOD_HPP_
#define PROT_BASE_MOD_HPP_

#include <string>
#include <vector>

#include "base/residue.hpp"
#include "base/residue_base.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Mod;
typedef std::shared_ptr<Mod> ModPtr;

class Mod {
 public:
  Mod(ResiduePtr ori_residue_ptr, ResiduePtr mod_residue_ptr):
    ori_residue_ptr_(ori_residue_ptr),
    mod_residue_ptr_(mod_residue_ptr) {}

  explicit Mod(xercesc::DOMElement* element);

  ResiduePtr getOriResiduePtr() { return ori_residue_ptr_;}

  ResiduePtr getModResiduePtr() { return mod_residue_ptr_;}

  bool isSame(ModPtr mod_ptr) {
    return ori_residue_ptr_ == mod_ptr->getOriResiduePtr()
        && mod_residue_ptr_ == mod_ptr->getModResiduePtr();
  }

  double getShift() {
    return mod_residue_ptr_->getMass() - ori_residue_ptr_->getMass();
  }

  void appendToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "mod";}

 private:
  ResiduePtr ori_residue_ptr_;
  ResiduePtr mod_residue_ptr_;
};

typedef std::vector<ModPtr> ModPtrVec;

}  // namespace prot
#endif
