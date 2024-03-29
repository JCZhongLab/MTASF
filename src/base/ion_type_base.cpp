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

#include <string>

#include "base/ion_type_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

IonTypePtrVec IonTypeBase::ion_type_ptr_vec_;

IonTypePtr IonTypeBase::ion_type_ptr_B_;

IonTypePtr IonTypeBase::ion_type_ptr_PREC_;

void IonTypeBase::initBase(const std::string &file_name) {
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = IonType::getXmlElementName();
    int ion_type_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < ion_type_num; i++) {
      xercesc::DOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      IonTypePtr ion_type_ptr = std::make_shared<IonType>(element);
      ion_type_ptr_vec_.push_back(ion_type_ptr);
      if (ion_type_ptr->getName() == getName_B()) {
        ion_type_ptr_B_ = ion_type_ptr;
      }
      if (ion_type_ptr->getName() == getName_PREC()) {
        ion_type_ptr_PREC_ = ion_type_ptr;
      }
    }
  }
}

IonTypePtr IonTypeBase::getIonTypePtrByName(const std::string &name) {
  for (size_t i = 0; i < ion_type_ptr_vec_.size(); i++) {
    std::string n = ion_type_ptr_vec_[i]->getName();
    if (n == name) {
      return ion_type_ptr_vec_[i];
    }
  }
  return IonTypePtr(nullptr);
}

}  // namespace prot
