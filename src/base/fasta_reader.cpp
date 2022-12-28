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
#include "base/fasta_reader.hpp"

namespace prot {

FastaReader::FastaReader(const std::string &file_name) {
  input_.open(file_name.c_str(), std::ios::in);
  if (!input_.is_open()) {
    LOG_ERROR("fasta file  " << file_name << " does not exist.");
    throw "fasta file does not exist.";
  }
  std::getline(input_, ori_name_);

  while (ori_name_.length() <= 1) {
    std::getline(input_, ori_name_);
  }
}


/*
 读取一段数据库中的序列，包括name, des, seq
 以及对原始序列进行变化的compAcidPtmPairVec();构建acid_ptm_pair_vec_
  初始的fasta数据中的序列没有[]，，其中数据结构全为
  （字母， No PTM）
  （字母， No PTM）
  （字母， No PTM）
  （字母， No PTM）
  字母变化：
  B - D
  Z - E
  X - A
  J - I
  其他字母不变。
  */
FastaSeqPtr FastaReader::getNextSeq() {
  if (!input_.is_open()) {
    return FastaSeqPtr(nullptr);
  }

  // get the letters of sequence
  std::string ori_seq;
  ori_name_ = string_util::trim(ori_name_);
  std::string prot_name = ori_name_.substr(1, ori_name_.size() - 1);
  std::string line;
  while (std::getline(input_, line)) {
    if (line.length() >= 1 && line.substr(0, 1) == ">") {
      ori_name_ = string_util::trim(line);
      return std::make_shared<FastaSeq>(prot_name, ori_seq);
    }
    line = string_util::trim(line);
    ori_seq = ori_seq + line;
    if (ori_seq.size() >= 1000000) {
      LOG_ERROR("Protein sequences are too long! Incorrect fasta file!");
      throw("fasta file error");
    }
  }
  input_.close();
  return std::make_shared<FastaSeq>(prot_name, ori_seq);
}

void FastaReader::close() {
  input_.close();
}

}  // namespace prot
