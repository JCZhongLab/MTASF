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


#ifndef PROT_FEATURE_ENV_RESCORE_HPP_
#define PROT_FEATURE_ENV_RESCORE_HPP_

#include <vector>

#include "feature/match_env.hpp"

namespace prot {
namespace EnvRescore {
void rescore(MatchEnvPtr2D &match_envs, const std::vector<std::vector<double> > para);
}  // namespace EnvRescore
}  // namespace prot
#endif
