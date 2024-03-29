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


#include <map>
#include <string>
#include <vector>
#include <ctime>

#include "base/base_data.hpp"
#include "base/file_util.hpp"
#include "base/string_util.hpp"
#include "base/fasta_util.hpp"

#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_sample_merge.hpp"
#include "prsm/prsm_stat.hpp"

#include "console/topmerge_argument.hpp"
#include "console/topmerge_process.hpp"

namespace prot {

int topMergeProcess(std::map<std::string, std::string> &arguments,
                    std::vector<std::string> &proteo_file_list) {

  Argument::outputArguments(std::cout, arguments);
  std::string resource_dir = arguments["resourceDir"];
  base_data::init(resource_dir);

  std::string ori_db_file_name = arguments["databaseFileName"];
  std::string db_file_name = ori_db_file_name + "_target";
  fasta_util::dbSimplePreprocess(ori_db_file_name, db_file_name);

  double error_tole = std::stod(arguments["errorTolerance"]);

  std::string base_path = file_util::absoluteDir(proteo_file_list[0]);
  std::string output_file_name = base_path + file_util::getFileSeparator() 
      + arguments["mergedOutputFileName"];
  LOG_DEBUG("Output file name " << output_file_name);
  std::string fixed_mod = arguments["fixedMod"];

  std::cout << "Merging files - started." << std::endl;
  PrsmSampleMergePtr merge_ptr 
      = std::make_shared<PrsmSampleMerge>(db_file_name, 
                                          proteo_file_list,
                                          output_file_name,
                                          fixed_mod,
                                          error_tole);
  merge_ptr->process();
  merge_ptr = nullptr;
  std::cout << "Merging files - finished." << std::endl;
  return 0;
}

}
