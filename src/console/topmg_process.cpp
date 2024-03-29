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
#include <algorithm>
#include <vector>

#include "base/version.hpp"
#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"
#include "base/file_util.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/feature_util.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_cluster.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_str_combine.hpp"

#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"

#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"

#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_processor.hpp"

#include "graph/graph.hpp"
#include "graph/proteo_graph.hpp"
#include "graphalign/graph_align_mng.hpp"
#include "graphalign/graph_align_processor.hpp"
#include "graphalign/graph_post_processor.hpp"

#include "mcmc/mcmc_mng.hpp"
#include "mcmc/mcmc_dpr_processor.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "console/topmg_argument.hpp"
#include "console/topmg_process.hpp"
#include "console/generate_ms_tag.hpp"
#include "console/filter_db.hpp"

#include "spec/msalign_writer.hpp"
#include "base/amino_acid_base.hpp"
#include "base/mass_constant.hpp"

namespace prot {
	using namespace std;
	const std::string cstr_ori_database_backup = "ori_database";
	const std::string cstr_ori_database_sub = "database_sub_I";
	const std::string cstr_ori_spfile_backup = "ori_spfile";

	void copy_file_fun(const std::string& ori_file_name,
		const std::string& st_file_name) {
		std::ifstream ori_db(ori_file_name);
		std::string line;
		std::ofstream standard_db;
		standard_db.open(st_file_name.c_str(), std::ios::out);

		while (std::getline(ori_db, line)) {
			standard_db << line << std::endl;
		}
		ori_db.close();
		standard_db.close();
	}

	void merge_function(std::string sp_file_name, int sp_nums, std::string result_name) {
		std::string sp_file_name_pre = sp_file_name.substr(0, sp_file_name.size() - 7);

		std::string line;
		std::ofstream result_file;
		result_file.open(std::string(sp_file_name_pre + result_name).c_str(), std::ios::out);

		std::string file_name = std::string(sp_file_name_pre + result_name + string_util::convertToString(0));
		std::ifstream ori_db(file_name);
		while (std::getline(ori_db, line)) {
			if (line == "</simple_prsm_list>")
				break;
			result_file << line << std::endl;
		}
		ori_db.close();
		remove(file_name.c_str());


		for (int i = 1; i < sp_nums; i++) {
			int cnt_index = 0;
			file_name = std::string(sp_file_name_pre + result_name + string_util::convertToString(i));
			ori_db.open(file_name, std::ios::in);
			while (std::getline(ori_db, line)) {
				if (cnt_index < 2)
				{
					++cnt_index;
					continue;
				}
				if (line == "</simple_prsm_list>")
					break;
				result_file << line << std::endl;
			}
			ori_db.close();
			remove(file_name.c_str());
		}
		result_file << "</simple_prsm_list>" << std::endl;
		result_file.close();
	}

	void merge_result_data(std::string sp_file_name, int sp_nums, std::string result_name) {
		merge_function(sp_file_name, sp_nums, result_name);
	}

	int TopMG_identify(std::map<std::string, std::string>& arguments) {
		try {
			std::vector<int> tag_nums;
			std::vector<double> mod_num;

			size_t db_size;
			double err;
			int tag_len;
			double trunc_ratio;


			std::cout << "Please enter terminal truncation ratio:" << std::endl;
			std::cin >> trunc_ratio;

			std::cout << "Please enter tag len:" << std::endl;
			std::cin >> tag_len;

			std::cout << "Please enter tag_tolerance:" << std::endl;
			std::cin >> err;

			std::cout << "Please enter filter database size:" << std::endl;
			std::cin >> db_size;

			std::time_t start = time(nullptr);
			char buf[50];
			std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

			arguments["startTime"] = buf;
			Argument::outputArguments(std::cout, arguments);

			std::string resource_dir = arguments["resourceDir"];

			base_data::init(resource_dir);

			LOG_DEBUG("Init base data completed");

			std::string db_file_name = arguments["databaseFileName"];
			std::string sp_file_name = arguments["spectrumFileName"];
			std::string ori_db_file_name = arguments["oriDatabaseFileName"];
			std::string var_mod_file_name = arguments["varModFileName"];
			
			std::string feature_file_name = sp_file_name.substr(0, sp_file_name.length() - 12) + ".feature";

			if (arguments["useFeatureFile"] == "true") {
				if (!boost::filesystem::exists(feature_file_name)) {
					LOG_ERROR("TopFD feature file does not exist!. Please use -x option in command line or select 'Missing MS1 feature file in GUI'.");
					return 1;
				}
			}

			int ptm_num = std::stoi(arguments["ptmNumber"]);
			LOG_DEBUG("num of unknown shfit " << ptm_num);
			int thread_num = std::stoi(arguments["threadNumber"]); 
			int filter_result_num = std::stoi(arguments["filteringResultNumber"]);
			double max_ptm_mass = std::stod(arguments["maxPtmMass"]);

			bool decoy = false;
			if (arguments["searchType"] == "TARGET+DECOY") {
				decoy = true;
			}
			LOG_DEBUG("Decoy " << decoy);
			LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
			int db_block_size = std::stoi(arguments["databaseBlockSize"]);

			PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

			int sp_num = msalign_util::geneSpIndex(sp_file_name, prsm_para_ptr->getSpParaPtr());
			copy_file_fun(sp_file_name, cstr_ori_spfile_backup);
			int max_mod_num = std::stoi(arguments["varPtmNumber"]); 
			int gap = std::stoi(arguments["proteoGraphGap"]);

			fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
			copy_file_fun(ori_db_file_name + std::string("_target_decoy"), cstr_ori_database_backup);

			int spectrum_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
			{
				MsAlignReader reader(cstr_ori_spfile_backup, 1, nullptr, prsm_para_ptr->getSpParaPtr()->getSkipList());
				DeconvMsPtr deconv_ms_ptr;
				int cnt = 0;
				std::string mod_file_name;
				std::vector<double> protein_mass;
				std::vector<std::vector<double>*> all_Trunc_mass;
				
				double temp_mass = 0;
				//building proteoform mass set
				std::cout << "constructing proteoform mass set..." << std::endl;
				get_protein_mass(cstr_ori_database_backup, protein_mass, all_Trunc_mass,trunc_ratio);

				std::vector<double> mod_list = get_modlist(var_mod_file_name, max_mod_num);
				std::vector<std::vector<double>> mod_protein = get_mod_protein_mass(protein_mass, mod_list);

				std::cout << "protein num:" << protein_mass.size() << std::endl;
				std::cout << "all_Trunc_mass num:" << all_Trunc_mass.size() << std::endl;
				std::cout << "mod protein num:" << mod_protein.size() << std::endl;

				TagGenerate tag_gen(tag_len, err);
				//std::cout << "DataBase filtering - started." << std::endl;
				while ((deconv_ms_ptr = reader.getNextMs()) != nullptr){
					std::cout << std::flush << "MTASF filtering - processing " << cnt
						<< " of " << spectrum_num << " spectra.\r";
					remove(sp_file_name.c_str());
					std::ofstream* output = new std::ofstream(sp_file_name);
					msalign_writer::write(*output, deconv_ms_ptr);
					output->close();
					//MF+TF
					get_filter_db(deconv_ms_ptr, mod_protein, all_Trunc_mass,tag_gen, sp_file_name, ori_db_file_name + std::string("_target_decoy"), tag_len, db_size, err, max_mod_num);
		
					fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size, false);

					//std::cout << "ASF-One PTM filtering - started." << std::endl;
					OnePtmFilterMngPtr one_ptm_filter_mng_ptr =
						std::make_shared<OnePtmFilterMng>(prsm_para_ptr, "topmg_one_filter", thread_num,
							var_mod_file_name, 1);
					
					//OnePtmFilterMng->output_file_ext_ = "topmg_one_filter"
					one_ptm_filter_mng_ptr->inte_num_ = 4;
					one_ptm_filter_mng_ptr->pref_suff_num_ = 4;
					one_ptm_filter_mng_ptr->comp_num_ = 4;
					OnePtmFilterProcessorPtr one_filter_processor =
						std::make_shared<OnePtmFilterProcessor>(one_ptm_filter_mng_ptr);
					one_filter_processor->process();
					one_filter_processor = nullptr;
					//std::cout << "ASF-One PTM filtering - finished." << std::endl;

					std::vector<std::string> input_exts;
					input_exts.emplace_back("topmg_one_filter_complete");
					input_exts.emplace_back("topmg_one_filter_prefix");
					input_exts.emplace_back("topmg_one_filter_suffix");
					input_exts.emplace_back("topmg_one_filter_internal");

					if (arguments["useAsfDiag"] == "true") {
						std::cout << "ASF-Diagonal PTM filtering - started." << std::endl;
						filter_result_num = 15;
						DiagFilterMngPtr diag_filter_mng_ptr
							= std::make_shared<DiagFilterMng>(prsm_para_ptr, filter_result_num,
								thread_num, "topmg_multi_filter",
								var_mod_file_name, 1);
						DiagFilterProcessorPtr diag_filter_processor
							= std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
						diag_filter_processor->process();
						diag_filter_processor = nullptr;
						std::cout << "ASF-Diagonal filtering - finished." << std::endl;

						input_exts.emplace_back("topmg_multi_filter");
					}

					//std::cout << "Combining filtering results - started." << std::endl;
					SimplePrsmStrCombinePtr asf_filter_combiner
						= std::make_shared<SimplePrsmStrCombine>(sp_file_name,
							input_exts,
							"topmg_graph_filter" + string_util::convertToString(cnt), 20 * input_exts.size());
					asf_filter_combiner->process();
					asf_filter_combiner = nullptr;
					//std::cout << "Combining filtering results - finished." << std::endl;
					cnt++;
				}
				std::cout << "MTASF filtering - processing " << cnt << " of " << spectrum_num << " spectra.\r";
				std::cout << "MTASF filtering - finished " << ".\r";
				reader.close();
				merge_result_data(sp_file_name, sp_num, "topmg_graph_filter");
				copy_file_fun(cstr_ori_spfile_backup, sp_file_name);
				remove(cstr_ori_spfile_backup.c_str());
			}

			////fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
			copy_file_fun(cstr_ori_database_backup, ori_db_file_name + std::string("_target_decoy"));
			fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size, false);
			int var_ptm_in_gap = std::min(std::stoi(arguments["varPtmNumInGap"]), max_mod_num);
			GraphAlignMngPtr ga_mng_ptr
				= std::make_shared<GraphAlignMng>(prsm_para_ptr,
					var_mod_file_name,
					ptm_num, max_mod_num,
					gap, var_ptm_in_gap, max_ptm_mass,
					thread_num, "topmg_graph_filter", "topmg_graph_align");
			std::cout << "Graph alignment - started." << std::endl;
			GraphAlignProcessorPtr ga_processor_ptr = std::make_shared<GraphAlignProcessor>(ga_mng_ptr);
			ga_processor_ptr->process();
			ga_processor_ptr = nullptr;
			std::cout << "Graph alignment - finished." << std::endl;

			std::cout << "Graph alignment post-processing - started." << std::endl;
			GraphPostProcessorPtr ga_post_processor_ptr
				= std::make_shared<GraphPostProcessor>(ga_mng_ptr, "topmg_graph_align", "topmg_graph_post");
			ga_post_processor_ptr->process();
			ga_post_processor_ptr = nullptr;
			std::cout << "Graph alignment post-processing - finished." << std::endl;
			std::cout << "E-value computation using MCMC - started." << std::endl;
			MCMCMngPtr mcmc_mng_ptr
				= std::make_shared<MCMCMng>(prsm_para_ptr, "topmg_graph_post", "topmg_evalue",
					var_mod_file_name, max_mod_num, thread_num);
			DprProcessorPtr processor = std::make_shared<DprProcessor>(mcmc_mng_ptr);
			processor->process();
			processor = nullptr;

			std::cout << "E-value computation using MCMC - finished." << std::endl;
			TopMG_post(arguments);
		}
		catch (const char* e) {
			std::cout << "[Exception]" << std::endl;
			std::cout << e << std::endl;
		}

		return 0;
	}

	int TopMG_post(std::map<std::string, std::string>& arguments) {
		try {

			std::string db_file_name = arguments["databaseFileName"];
			std::string sp_file_name = arguments["spectrumFileName"];
			std::string ori_db_file_name = arguments["oriDatabaseFileName"];
			std::string var_mod_file_name = arguments["varModFileName"];
			std::string resource_dir = arguments["resourceDir"];

			base_data::init(resource_dir);

			LOG_DEBUG("Initialization completed");

			PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
			int n_top = std::stoi(arguments["numOfTopPrsms"]);


			std::cout << "Finding PrSM clusters - started." << std::endl;
			if (arguments["useFeatureFile"] == "true") {
				// TopFD msalign file with feature ID
				double prec_error_tole = 1.2;
				ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
				PrsmFeatureClusterPtr prsm_clusters
					= std::make_shared<PrsmFeatureCluster>(db_file_name,
						sp_file_name,
						"topmg_evalue",
						"topmg_cluster",
						fix_mod_list,
						prec_error_tole,
						prsm_para_ptr);
				prsm_clusters->process();
				prsm_clusters = nullptr;
			}
			else {
				double ppo;
				std::istringstream(arguments["errorTolerance"]) >> ppo;
				ppo = ppo / 1000000.0;  //  15 / 1000000.0
				PrsmClusterPtr prsm_clusters
					= std::make_shared<PrsmCluster>(db_file_name, sp_file_name,
						"topmg_evalue", prsm_para_ptr->getFixModPtrVec(),  // fix_mod_list 为空
						"topmg_cluster", ppo);
				prsm_clusters->process();
				prsm_clusters = nullptr;
			}
			std::cout << "Finding PrSM clusters - finished." << std::endl;
			if (arguments["searchType"] == "TARGET") {
				std::cout << "Top PrSM selecting - started" << std::endl;
				PrsmTopSelectorPtr selector
					= std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "topmg_cluster", "topmg_top", n_top);
				selector->process();
				selector = nullptr;
				std::cout << "Top PrSM selecting - finished." << std::endl;
			}
			else {
				std::cout << "Top PrSM selecting - started " << std::endl;
				PrsmTopSelectorPtr selector
					= std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "topmg_cluster", "topmg_top_pre", n_top);
				selector->process();  //选出一个光谱匹配到所有序列中E-value值最小的一条序列
				selector = nullptr;
				std::cout << "Top PrSM selecting - finished." << std::endl;
				std::cout << "FDR computation - started. " << std::endl;
				PrsmFdrPtr fdr = std::make_shared<PrsmFdr>(db_file_name, sp_file_name, "topmg_top_pre", "topmg_top");
				fdr->process();
				fdr = nullptr;
				std::cout << "FDR computation - finished." << std::endl;
			}
			std::string cutoff_type = arguments["cutoffSpectralType"];
			std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
			double cutoff_value = std::stod(arguments["cutoffSpectralValue"]);
			PrsmCutoffSelectorPtr cutoff_selector
				= std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "topmg_top",
					"topmg_prsm_cutoff", cutoff_type, cutoff_value);
			cutoff_selector->process();
			cutoff_selector = nullptr;
			std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

			std::time_t end = time(nullptr);
			char buf[50];
			std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&end));
			arguments["endTime"] = buf;

			std::string argu_str = Argument::outputCsvArguments(arguments);

			std::cout << "Outputting PrSM table - started." << std::endl;
			PrsmTableWriterPtr table_out
				= std::make_shared<PrsmTableWriter>(prsm_para_ptr, argu_str, "topmg_prsm_cutoff",
					std::string("_topmg_prsm.csv"));
			table_out->write();
			table_out = nullptr;
			std::cout << "Outputting PrSM table - finished." << std::endl;

			std::cout << "Generating PrSM xml files - started." << std::endl;
			XmlGeneratorPtr xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, "topmg_prsm_cutoff", "topmg_prsm_cutoff");

			xml_gene->process();
			xml_gene = nullptr;
			std::cout << "Generating PrSM xml files - finished." << std::endl;

			std::cout << "Converting PrSM xml files to html files - started." << std::endl;
			translate(arguments, "topmg_prsm_cutoff");
			std::cout << "Converting PrSM xml files to html files - finished." << std::endl;

			cutoff_type = (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR" : "EVALUE";
			std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
			std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
			cutoff_selector = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "topmg_top",
				"topmg_form_cutoff", cutoff_type, cutoff_value);
			cutoff_selector->process();
			cutoff_selector = nullptr;
			std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

			std::cout << "Selecting top PrSMs for proteoforms - started." << std::endl;
			PrsmFormFilterPtr form_filter
				= std::make_shared<PrsmFormFilter>(db_file_name, sp_file_name, "topmg_form_cutoff",
					"topmg_form_cutoff_form");
			form_filter->process();
			form_filter = nullptr;
			std::cout << "Selecting top PrSMs for proteoforms - finished." << std::endl;

			std::cout << "Outputting proteoform table - started." << std::endl;
			PrsmTableWriterPtr form_out
				= std::make_shared<PrsmTableWriter>(prsm_para_ptr, argu_str, "topmg_form_cutoff_form",
					"_topmg_proteoform.csv");
			form_out->write();
			form_out = nullptr;
			std::cout << "Outputting proteoform table - finished." << std::endl;
			exit(1);

			std::cout << "Generating proteoform xml files - started." << std::endl;
			xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, "topmg_form_cutoff", "topmg_proteoform_cutoff");
			xml_gene->process();
			xml_gene = nullptr;
			std::cout << "Generating proteoform xml files - finished." << std::endl;

			std::cout << "Converting proteoform xml files to html files - started." << std::endl;
			translate(arguments, "topmg_proteoform_cutoff");
			std::cout << "Converting proteoform xml files to html files - finished." << std::endl;

		}
		catch (const char* e) {
			std::cout << "[Exception]" << std::endl;
			std::cout << e << std::endl;
		}
		return 0;
	}

	int TopMGProcess(std::map<std::string, std::string>& arguments) {
		if (TopMG_identify(arguments) != 0) {
			return 1;
		}
		/*修改*/
		return 0;
		/**/
		//return TopMG_post(arguments);
	}

	int TopMGProgress_multi_file(std::map<std::string, std::string>& arguments,
		const std::vector<std::string>& spec_file_lst) {

		std::string base_path = file_util::absoluteDir(spec_file_lst[0]);
		std::string base_name = base_path + file_util::getFileSeparator()
			+ arguments["combinedOutputName"];

		std::time_t start = time(nullptr);
		char buf[50];
		std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
		std::string combined_start_time = buf;

		std::cout << "TopMG " << prot::version_number << std::endl;

		for (size_t k = 0; k < spec_file_lst.size(); k++) {
			std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
			std::string start_time = buf;
			arguments["startTime"] = start_time;
			arguments["spectrumFileName"] = spec_file_lst[k];
			if (prot::TopMGProcess(arguments) != 0) {
				return 1;
			}
		}

		if (spec_file_lst.size() > 1 && arguments["combinedOutputName"] != "") {
			std::cout << "Merging files - started." << std::endl;
			int N = 1000000;
			// merge msalign files
			prot::msalign_util::mergeMsalignFiles(spec_file_lst, N, base_name + "_ms2.msalign");
			// merge feature files
			std::vector<std::string> feature_file_lst(spec_file_lst.size());
			for (size_t i = 0; i < spec_file_lst.size(); i++) {
				std::string sp_file_name = spec_file_lst[i];
				feature_file_lst[i] = sp_file_name.substr(0, sp_file_name.length() - 12) + ".feature";
			}
			prot::feature_util::mergeFeatureFiles(feature_file_lst, N, base_name + ".feature");
			// merge EVALUE files
			std::vector<std::string> prsm_file_lst(spec_file_lst.size());
			for (size_t i = 0; i < spec_file_lst.size(); i++) {
				prsm_file_lst[i] = prot::file_util::basename(spec_file_lst[i]) + ".topmg_evalue";
			}
			prot::prsm_util::mergePrsmFiles(prsm_file_lst, N, base_name + "_ms2.topmg_evalue");
			std::cout << "Merging files - finished." << std::endl;

			std::string sp_file_name = base_name + "_ms2.msalign";
			arguments["spectrumFileName"] = sp_file_name;
			arguments["startTime"] = combined_start_time;
			prot::TopMG_post(arguments);
		}

		if (arguments["keepTempFiles"] != "true") {
			std::cout << "Deleting temporary files - started." << std::endl;
			std::string ori_db_file_name = arguments["oriDatabaseFileName"];

			for (size_t k = 0; k < spec_file_lst.size(); k++) {
				std::string sp_file_name = spec_file_lst[k];
				prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_topmg_proteoform_cutoff_xml");
				prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_topmg_prsm_cutoff_xml");
				prot::file_util::cleanTopmgDir(ori_db_file_name, sp_file_name);
			}

			if (spec_file_lst.size() > 1 && arguments["combinedOutputName"] != "") {
				std::string sp_file_name = base_name + "_ms2.msalign";
				prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_topmg_proteoform_cutoff_xml");
				prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_topmg_prsm_cutoff_xml");
				prot::file_util::cleanTopmgDir(ori_db_file_name, sp_file_name);
			}
			std::cout << "Deleting temporary files - finished." << std::endl;
		}

		std::cout << "TopMG finished." << std::endl << std::flush;
		return 0;
	}

}  // namespace prot
