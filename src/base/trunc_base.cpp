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

#include "base/logger.hpp"
#include "base/trunc_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

TruncPtrVec TruncBase::trunc_ptr_vec_;

void TruncBase::initBase(const std::string &file_name) {
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Trunc::getXmlElementName();
    int trunc_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < trunc_num; i++) {
      xercesc::DOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      TruncPtr trunc_ptr = std::make_shared<Trunc>(element);
      trunc_ptr_vec_.push_back(trunc_ptr);
    }
  }
}

TruncPtr TruncBase::getTruncPtrByName(const std::string &name) {
  for (size_t i = 0; i < trunc_ptr_vec_.size(); i++) {
    std::string n = trunc_ptr_vec_[i]->getName();
    if (n == name) {
      return trunc_ptr_vec_[i];
    }
  }
  return TruncPtr(nullptr);
}

TruncPtr TruncBase::getTruncPtrFromXml(xercesc::DOMElement * element) {
  std::string name = Trunc::getNameFromXml(element);
  TruncPtr trunc_ptr = getTruncPtrByName(name);
  return trunc_ptr;
}

}  // namespace prot

