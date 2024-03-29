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


#ifndef PROT_BASE_EXTREME_VALUE_HPP_
#define PROT_BASE_EXTREME_VALUE_HPP_

#include <memory>
#include <vector>
#include <string>

#include "base/xml_dom_document.hpp"

namespace prot {

class ExtremeValue;
typedef std::shared_ptr<ExtremeValue> ExtremeValuePtr;

class ExtremeValue {
 public:
  ExtremeValue(double one_prot_prob, double test_num, double adjust_factor):
      one_prot_prob_(one_prot_prob),
      test_num_(test_num),
      adjust_factor_(adjust_factor) {
        init();
      }

  explicit ExtremeValue(xercesc::DOMElement* element);

  double getPValue() {return p_value_;}

  double getEValue() {return e_value_;}

  double getOneProtProb() {return one_prot_prob_;}

  double getTestNum() {return test_num_;}

  double getAdjustFactor() { return adjust_factor_;}

  void setOneProtProb(double one_prot_prob);

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "extreme_value";}

  static ExtremeValuePtr getMaxEvaluePtr();

 private:
  // one_prot_prob is the probability that the spectrum and a randem problem
  // have a protein-spectrum-match with a score no less than the threshold
  double one_prot_prob_;

  double test_num_;

  double adjust_factor_;

  double p_value_;

  double e_value_;

  void init();
};

typedef std::vector<ExtremeValuePtr> ExtremeValuePtrVec;

}  // namespace prot

#endif
