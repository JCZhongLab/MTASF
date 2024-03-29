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
#include "base/amino_acid_base.hpp"
#include "base/ptm_base.hpp"
#include "base/fasta_seq.hpp"
#include "base/residue_base.hpp"
#include "base/residue_util.hpp"

namespace prot {

namespace residue_util {

ResiduePtrVec convertStrToResiduePtrVec(const std::string &seq) {
  ResiduePtrVec residue_ptr_vec;
  std::string seq2 = FastaSeq::rmChar(seq);
  for (size_t i = 0; i < seq2.length(); i++) {
    AminoAcidPtr acid_ptr = AminoAcidBase::getAminoAcidPtrByOneLetter(seq2.substr(i, 1));
    PtmPtr ptm_ptr = PtmBase::getEmptyPtmPtr();
    ResiduePtr residue_ptr = ResidueBase::getBaseResiduePtr(acid_ptr, ptm_ptr);
    residue_ptr_vec.push_back(residue_ptr);
  }
  return residue_ptr_vec;
}

void applyFixedMod(ResiduePtrVec &residue_ptrs, const ModPtrVec &fix_mod_ptr_vec) {
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    for (size_t j = 0; j < fix_mod_ptr_vec.size(); j++) {
      if (residue_ptrs[i] == fix_mod_ptr_vec[j]->getOriResiduePtr()) {
        residue_ptrs[i] = fix_mod_ptr_vec[j]->getModResiduePtr();
        break;
      }
    }
  }
}

ResiduePtrVec convertStrToResiduePtrVec(const std::string & seq,
                                        const ModPtrVec &fix_mod_ptr_vec) {
  ResiduePtrVec residue_ptrs = convertStrToResiduePtrVec(seq);
  applyFixedMod(residue_ptrs, fix_mod_ptr_vec);
  return residue_ptrs;
}

ResiduePtrVec convertStrToResiduePtrVec(const StringPairVec &string_pair_vec) {
  ResiduePtrVec residue_ptr_vec;
  for (size_t i = 0; i < string_pair_vec.size(); i++) {
    std::string acid_str = string_pair_vec[i].first;
    AminoAcidPtr acid_ptr = AminoAcidBase::getAminoAcidPtrByOneLetter(acid_str);
    std::string ptm_str = string_pair_vec[i].second;
    PtmPtr ptm_ptr = PtmBase::getPtmPtrByAbbrName(ptm_str);
    ResiduePtr residue_ptr = ResidueBase::getBaseResiduePtr(acid_ptr, ptm_ptr);
    residue_ptr_vec.push_back(residue_ptr);
  }
  return residue_ptr_vec;
}

ResiduePtrVec convertStrToResiduePtrVec(const StringPairVec &string_pair_vec,
                                                     const ModPtrVec &fix_mod_ptr_vec) {
  ResiduePtrVec residue_ptrs = convertStrToResiduePtrVec(string_pair_vec);
  applyFixedMod(residue_ptrs, fix_mod_ptr_vec);
  return residue_ptrs;
}

int findResidue(const ResiduePtrVec &residue_list, ResiduePtr residue_ptr) {
  for (size_t i = 0; i < residue_list.size(); i++) {
    if (residue_list[i] == residue_ptr) {
      return i;
    }
  }
  return -1;
}

double compResiduePtrVecMass(const ResiduePtrVec &ptr_vec) {
  double mass = 0;
  for (size_t i = 0; i < ptr_vec.size(); i++) {
    mass += ptr_vec[i]->getMass();
  }
  return mass;
}

double compResiduePtrVecMass(const std::string & seq) {
  return compResiduePtrVecMass(convertStrToResiduePtrVec(seq));
}

double compResiduePtrVecMass(const std::string & seq,
                             const ModPtrVec &fix_mod_ptr_vec) {
  return compResiduePtrVecMass(convertStrToResiduePtrVec(seq, fix_mod_ptr_vec));
}

}  // namespace residue_util

}  // namespace prot
