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
#include <vector>
#include <string>

#include "boost/algorithm/string.hpp"

#include "base/activation_base.hpp"
#include "base/string_util.hpp"
#include "spec/msalign_reader.hpp"

namespace prot {

std::vector<std::string> MsAlignReader::readOneSpectrum() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    line = string_util::trim(line);
    if (line == "BEGIN IONS") {
      line_list.push_back(line);
    } 
	else if (line == "END IONS") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    } 
	else if (line == "" || line[0] == '#') {
      continue;
    } 
	else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

void MsAlignReader::readNext() {
  deconv_ms_ptr_ = nullptr;
  spectrum_str_vec_ = readOneSpectrum();
  /*
  从光谱中读取了一段字符串BEGIN IONS开始，END IONS结束，每一行为数组的一项。
  */
  if (spectrum_str_vec_.size() == 0) {
    input_.close();
    return;
  }
  int id = -1;
  int prec_id = 0;
  std::string scans;
  double retention_time = -1;
  std::string activation;
  std::string title;
  int level = 2;
  int ms_one_id = -1;
  int ms_one_scan = -1;
  double prec_mass = -1;
  int prec_charge = -1;
  double prec_inte = -1;
  int feature_id = -1;
  double feature_inte = -1;
  std::vector<std::string> strs;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0, 1);
    if (letter >= "A" && letter <= "Z") {
      strs = string_util::split(spectrum_str_vec_[i], '=');
      if (strs[0] == "ID") {
        id = std::stoi(strs[1]);
      }

      if (strs[0] == "PRECURSOR_ID") {  //无该字段， 默认为0
        prec_id = std::stoi(strs[1]);
      } else if (strs[0] == "SCANS") {
        scans = strs[1];
      } else if (strs[0] == "RETENTION_TIME") {
        retention_time = std::stod(strs[1]);
      } else if (strs[0] == "ACTIVATION") {
        activation = strs[1];
      } else if (strs[0] == "TITLE") { //无该字段， 默认为空字符串
        title = strs[1];
      } else if (strs[0] == "LEVEL") { //无该字段， 默认为2
        level = std::stoi(strs[1]);
      } else if (strs[0] == "MS_ONE_ID") {
        ms_one_id = std::stod(strs[1]);
      } else if (strs[0] == "MS_ONE_SCAN") {
        ms_one_scan = std::stod(strs[1]);
      } else if (strs[0] == "PRECURSOR_MASS") {
        prec_mass = std::stod(strs[1]);
      } else if (strs[0] == "PRECURSOR_CHARGE") {
        prec_charge = std::stoi(strs[1]);
      } else if (strs[0] == "PRECURSOR_INTENSITY") {
        prec_inte = std::stod(strs[1]);
      } else if (strs[0] == "FEATURE_ID") { //无该字段， 默认为-1
        feature_id = std::stoi(strs[1]);
      } else if (strs[0] == "FEATURE_INTENSITY") { //无该字段， 默认为-1
        feature_inte = std::stod(strs[1]);
      }
    }
  }
  if (id < 0 || prec_charge < 0 || prec_mass < 0) {
    LOG_WARN("Input file format error: sp id " << id << " prec_chrg "
             << prec_charge << " prec mass " << prec_mass);
  }

  MsHeaderPtr header_ptr = std::make_shared<MsHeader>();
  header_ptr->setFileName(file_name_);
  header_ptr->setId(id);       //  光谱中的ID字段
  header_ptr->setPrecId(prec_id); // 无， 0

  if (scans != "") {
    header_ptr->setScans(scans);  //  光谱中的SCANS字段
  } else {
    header_ptr->setScans("");
  }
  header_ptr->setRetentionTime(retention_time); //  光谱中的RETENTION_TIME字段 
  // LOG_DEBUG("retention time " << retention_time);

  if (title != "") {
    std::stringstream ss;
    ss << "sp_" << id;
    header_ptr->setTitle(ss.str());
  } else {
    header_ptr->setTitle(title);  // 无 ， 空字符串
  }

  if (activation_ptr_ != nullptr) {
    header_ptr->setActivationPtr(activation_ptr_);
  } else if (activation != "") { //  光谱中的ACTIVATION字段 
    ActivationPtr activation_ptr = ActivationBase::getActivationPtrByName(activation); 
	// 通过avtivation.xml构建的数据结构数组，通过name查找，取出对应的activationPtr

    header_ptr->setActivationPtr(activation_ptr);  
  }
  header_ptr->setMsLevel(level); // 无 ， 2

  header_ptr->setMsOneId(ms_one_id); //  光谱中的MS_ONE_ID字段

  header_ptr->setMsOneScan(ms_one_scan); //  光谱中的MS_ONE_SCAN字段

  header_ptr->setPrecMonoMz(prec_mass /prec_charge + mass_constant::getProtonMass());
  //  PRECURSOR_MASS / PRECURSOR_CHARGE + 1.007276

  header_ptr->setPrecCharge(prec_charge); //  光谱中的PRECURSOR_CHARGE字段

  header_ptr->setPrecInte(prec_inte); //  光谱中的PRECURSOR_INTENSITY字段

  header_ptr->setFeatureId(feature_id); // 无 ， -1

  header_ptr->setFeatureInte(feature_inte); // 无 ， -1

  std::vector<DeconvPeakPtr> peak_ptr_list;
  int idx = 0;
  for (size_t i = 1; i < spectrum_str_vec_.size() - 1; i++) {
    std::string letter = spectrum_str_vec_[i].substr(0, 1);
    if (letter >= "0" && letter <= "9") {
      boost::split(strs, spectrum_str_vec_[i], boost::is_any_of("\t "));
      double mass = std::stod(strs[0]);
      double inte = std::stod(strs[1]);
      int charge = std::stoi(strs[2]);
      DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(idx, mass, inte, charge);
      peak_ptr_list.push_back(peak_ptr);
      idx++;
      if (static_cast<int>(peak_ptr_list.size()) >= this->peak_num_limit_) break;
    }
  }

  deconv_ms_ptr_ = std::make_shared<Ms<DeconvPeakPtr> >(header_ptr, peak_ptr_list);

  current_++;
}

DeconvMsPtr MsAlignReader::getNextMs() {
  readNext();
  while (deconv_ms_ptr_ != nullptr
         && skip_list_.find(deconv_ms_ptr_->getMsHeaderPtr()->getScansString()) != skip_list_.end()) {
    readNext();
  }
  return deconv_ms_ptr_;
}

std::vector<SpectrumSetPtr> MsAlignReader::getNextSpectrumSet(SpParaPtr sp_para_ptr) {
  std::vector<SpectrumSetPtr> spec_set_vec;
  DeconvMsPtrVec deconv_ms_ptr_vec;
  for (int i = 0; i < group_spec_num_; i++) {
    readNext();
    while (deconv_ms_ptr_ != nullptr
           && skip_list_.find(deconv_ms_ptr_->getMsHeaderPtr()->getScansString()) != skip_list_.end()) {
      readNext();
    }
    if (deconv_ms_ptr_ == nullptr) {
      spec_set_vec.push_back(nullptr);
      return spec_set_vec;
    }
    deconv_ms_ptr_vec.push_back(deconv_ms_ptr_);
  }
  double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  // prec_mono_mass = PRECURSOR_MASS / PRECURSOR_CHARGE + 1.007276

  // LOG_DEBUG("prec mass " << prec_mono_mass);
  int count = 1;
  for (int i = 1; i < group_spec_num_; i++) // group_spec_num_ 为1 不进入循环
  {
    double new_mass = deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getPrecMonoMass();
    if (std::abs(prec_mono_mass - new_mass) < 0.5) {
      prec_mono_mass = (prec_mono_mass * count + new_mass)/ (count+1);
      count++;
    }
  }
  // LOG_DEBUG("prec mass result " << prec_mono_mass);
  std::vector<double> prec_errors;
  prec_errors.push_back(0);
  for (int i = 1; i <= sp_para_ptr->prec_error_; i++) {
    prec_errors.push_back(- i * mass_constant::getIsotopeMass());
    prec_errors.push_back(i * mass_constant::getIsotopeMass());
  }
  //prec_errors (0, -1.00235, 1.00235)
  for (size_t i = 0; i< prec_errors.size(); i++) {
    spec_set_vec.push_back(std::make_shared<SpectrumSet>(deconv_ms_ptr_vec,
                                                         sp_para_ptr,
                                                         prec_mono_mass + prec_errors[i]));
  }
  /*
  spec_set_vec (
            SpectrumSet(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass + 0),
			SpectrumSet(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass - 1.00235),
			SpectrumSet(deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass + 1.00235)
  )
  */
  return spec_set_vec;
}

void MsAlignReader::close() {
  input_.close();
}

}  // namespace prot
