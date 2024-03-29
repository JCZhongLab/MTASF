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


#ifndef PROT_FEATURE_DECONV_DATA_HPP_
#define PROT_FEATURE_DECONV_DATA_HPP_

#include <memory>
#include <vector>

#include "spec/peak.hpp"

namespace prot {

class DeconvData {
 public:
  DeconvData(PeakPtrVec &peak_list, double max_mass, int max_charge, double win_size); 

	int getBgnPeak(int i) {return win_bgn_peaks_[i];}

	int getEndPeak(int i) {return win_end_peaks_[i];}

	int getIntervalPeakNum(int i) {return win_end_peaks_[i] - win_bgn_peaks_[i] + 1;}

	double getMaxMass() {return max_mass_;}

	int getMaxCharge() {return max_charge_;}

	PeakPtrVec getPeakList() {return peak_list_;}

	int getWinId(int i) {return win_ids_[i];}

	int getWinNum() {return win_num_;}

	void setMaxCharge(int max_chrg) {max_charge_ = max_chrg;}

	void setMaxMass(double mass) {max_mass_ = mass;}

 private:
	PeakPtrVec peak_list_;
  double max_mass_;
	int max_charge_;

	// the number of windows 
  int win_num_;
	// the length of each window 
  double win_size_;

	// the window id for each peak 
  std::vector<int> win_ids_;
	// the first peak of each window 
  std::vector<int> win_bgn_peaks_;
	// the last peak of each window 
  std::vector<int> win_end_peaks_;

  void initWinId();

  void initWinBgnEnd();
};

typedef std::shared_ptr<DeconvData> DeconvDataPtr;

}

#endif
