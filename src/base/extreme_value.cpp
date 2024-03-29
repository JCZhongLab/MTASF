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


#include <cmath>
#include <limits>
#include <string>
#include <algorithm>

#include "base/base_data.hpp"
#include "base/extreme_value.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

void ExtremeValue::setOneProtProb(double one_prot_prob) {
  one_prot_prob_ = one_prot_prob;
  init();
}

ExtremeValue::ExtremeValue(xercesc::DOMElement* element) {
  one_prot_prob_ = xml_dom_util::getScientificChildValue(element, "one_protein_probability", 0);
  test_num_ = xml_dom_util::getScientificChildValue(element, "test_number", 0);
  adjust_factor_ = xml_dom_util::getDoubleChildValue(element, "adjust_factor", 0);
  init();
}

void ExtremeValue::init() {
  e_value_ = one_prot_prob_ * test_num_ * adjust_factor_;
  if (one_prot_prob_ >= 1 || test_num_ == base_data::getMaxDouble()) {
    p_value_  = 1.0;
  } else {
    double n = std::max(test_num_ * adjust_factor_, 1.0);
    // approximation of 1 - (1- one_prot_prob)^n
    p_value_ =  n * one_prot_prob_
        - (n * (n - 1)) / 2 * one_prot_prob_ * one_prot_prob_
        + (n * (n - 1) * (n - 2)) / 6 * std::pow(one_prot_prob_, 3);
    if (p_value_ > 1.0) {
      p_value_ = 1.0;
    }
  }
}

void ExtremeValue::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = string_util::convertToScientificStr(one_prot_prob_, 4);
  xml_doc->addElement(element, "one_protein_probability", str.c_str());
  str = string_util::convertToScientificStr(test_num_, 4);
  xml_doc->addElement(element, "test_number", str.c_str());
  str = string_util::convertToString(adjust_factor_);
  xml_doc->addElement(element, "adjust_factor", str.c_str());
  str = string_util::convertToScientificStr(p_value_, 4);
  xml_doc->addElement(element, "p_value", str.c_str());
  str = string_util::convertToScientificStr(e_value_, 4);
  xml_doc->addElement(element, "e_value", str.c_str());
  parent->appendChild(element);
}

ExtremeValuePtr ExtremeValue::getMaxEvaluePtr() {
  ExtremeValuePtr evalue_ptr
      = std::make_shared<ExtremeValue>(1.0, base_data::getMaxDouble(), 1.0);
  return evalue_ptr;
}

}  // namespace prot
