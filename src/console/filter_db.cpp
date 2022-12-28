#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <map>

#include <ctime>

#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"
#include "base/version.hpp"
#include "base/file_util.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/feature_util.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_cluster.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/simple_prsm_str_combine.hpp"

#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"
#include "zeroptmsearch/zero_ptm_search_mng.hpp"
#include "zeroptmsearch/zero_ptm_search_processor.hpp"

#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
#include "oneptmsearch/one_ptm_search_processor.hpp"

#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_processor.hpp"

#include "ptmsearch/ptm_search_processor.hpp"

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/evalue_processor.hpp"

#include "local/local_mng.hpp"
#include "local/local_processor.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "console/toppic_argument.hpp"

#include "base/amino_acid_base.hpp"
#include "base/mass_constant.hpp"
#include "spec/msalign_writer.hpp"

#include "console/toppic_process.hpp"
#include "console/filter_db.hpp"
#include <exception>
#include<numeric>
#include<algorithm>
#include<cmath>
namespace prot {
	double get_Ms_seq(const DeconvMsPtr deconv_ms_ptr, std::vector<double>& d_vector)
	{
		MsHeaderPtr header_ptr = deconv_ms_ptr->getMsHeaderPtr();
		std::vector<double> x_vec;
		std::vector<double> m_minus_x_vec;
		std::vector<double> x_minus_water_vec;
		std::vector<double> m_minus_x_minus_water_vec;

		double preMass = header_ptr->getPrecMonoMass();
		return preMass;

	}


	//Combination modification(PTM)
	void Combination(std::vector<double>& a, std::vector<double>& b, std::vector<double>& mod_mass, int l, int m, int M) {

		int N = a.size();
		if (m == 0) {
			//for (auto i : b)
			//	cout << i << ' ';
			double mass = accumulate(b.begin(), b.end(), 0.);
			mass = round(mass*1000.0)/1000.0;

			mod_mass.push_back(mass);
			//cout << endl;
			return;
		}
		for (int i = l; i < N; i++) {
			b[M - m] = a[i];
			Combination(a, b, mod_mass, i + 1, m - 1, M);
		}
	}

	//get combination modification(PTM) mass set
	std::vector<double> get_modlist(std::string var_mod_name, int max_mod) {
		std::vector<double> mod;
		std::vector<double> temp_mod;
		std::vector<double> mod_mass;
		std::vector<double> all_mod;
		std::ifstream var_file(var_mod_name);
		if (!var_file.is_open())
			return mod_mass;
		std::string line;
		while (std::getline(var_file, line)) {
			if (line[0] == '#') continue;
			if (line == "") continue;
			mod.push_back(atof(string_util::split(line, ',')[1].c_str()));
		}
		for (int i = 1; i <= max_mod; i++) {
			//vector<double>().swap(temp_mod);
			temp_mod.resize(i);
			all_mod.insert(all_mod.end(), mod.begin(), mod.end());
			Combination(all_mod, temp_mod, mod_mass, 0, i, i);
		}
		sort(mod_mass.begin(), mod_mass.end());
		mod_mass.erase(unique(mod_mass.begin(), mod_mass.end()), mod_mass.end());
		return mod_mass;
	}


	void get_One_Ms_vec(const DeconvMsPtr deconv_ms_ptr, std::vector<double>& ms_vector) {
		MsHeaderPtr header_ptr = deconv_ms_ptr->getMsHeaderPtr();
		double preMass = header_ptr->getPrecMonoMass();
		double waterMass = prot::mass_constant::getWaterMass();
		for (size_t i = 0; i < deconv_ms_ptr->size(); i++) {
			double x = deconv_ms_ptr->getPeakPtr(i)->getPosition();
			ms_vector.emplace_back(x);
		}
		std::sort(ms_vector.begin(), ms_vector.end(), std::less<double>());
	}

	void get_protein_mass(const std::string cstr_ori_database_backup, std::vector<double>& protein_mass, std::vector<std::vector<double>*>& all_Trunc_mass, double trunc) {
		double temp_mass;
		double N_mass;
		double C_mass;

		std::vector<double> trunc_N;
		std::vector<double> trunc_C;
		int trunc_num;
		FastaReader reader(cstr_ori_database_backup);
		FastaSeqPtr seq_info = reader.getNextSeq();
		double total_mass;
		size_t prt_str_size;
		while (seq_info != nullptr) {
			StringPairVec str_pair_vec = seq_info->getAcidPtmPairVec();
			total_mass = 0;
			prt_str_size = str_pair_vec.size();
			for (size_t i = 0; i < prt_str_size; i++) {
				std::string acid_str = str_pair_vec[i].first;
				total_mass += AminoAcidBase::getAminoAcidPtrByOneLetter(acid_str)->getMonoMass();
			}
			protein_mass.emplace_back(total_mass);

			//terminal truncation mass
			trunc_num = prt_str_size * trunc;

			N_mass = 0;
			C_mass = 0;
			std::vector<double>* trunc_vec = new std::vector<double>();
			trunc_vec->emplace_back(0);
			if (trunc_num >= 1) {
				for (size_t i = 0; i < trunc_num; i++) {
					N_mass += AminoAcidBase::getAminoAcidPtrByOneLetter(str_pair_vec[i].first)->getMonoMass();
					trunc_N.emplace_back(N_mass);
				}


				for (size_t i = 0; i < trunc_num; i++) {
					C_mass += AminoAcidBase::getAminoAcidPtrByOneLetter(str_pair_vec[prt_str_size - i - 1].first)->getMonoMass();
					trunc_C.emplace_back(C_mass);
				}

				//C and N truncation
				for (size_t i = 0; i < trunc_num / 2; i++) {
					for (size_t j = 0; j < trunc_num / 2; j++) {
						trunc_vec->emplace_back(trunc_N[i] + trunc_C[j]);
					}
				}

				trunc_vec->insert(trunc_vec->end(), trunc_N.begin(), trunc_N.end());
				trunc_vec->insert(trunc_vec->end(), trunc_C.begin(), trunc_C.end());
				sort(trunc_vec->begin(), trunc_vec->end());
				trunc_vec->erase(unique(trunc_vec->begin(), trunc_vec->end()), trunc_vec->end());
			}
			all_Trunc_mass.push_back(trunc_vec);
			std::vector<double>().swap(trunc_C);
			std::vector<double>().swap(trunc_N);

			seq_info = reader.getNextSeq();
		}
	}

	std::vector<std::vector<double>> get_mod_protein_mass(const std::vector<double>& protein_mass, std::vector<double>& mod_list) {
		std::vector<std::vector<double>> mod_protein_mass;
		std::vector<double> one_protein_mass;

		for (size_t i = 0; i < protein_mass.size(); i++) {
			one_protein_mass.emplace_back(protein_mass[i]);
			for (size_t k = 0; k < mod_list.size(); k++)
				one_protein_mass.emplace_back(protein_mass[i] + mod_list[k]);
			sort(one_protein_mass.begin(), one_protein_mass.end());
			//one_protein_mass.erase(unique(one_protein_mass.begin(), one_protein_mass.begin()), one_protein_mass.end());
			mod_protein_mass.push_back(one_protein_mass);
			std::vector<double>().swap(one_protein_mass);
		}

		return mod_protein_mass;
	}

	void gen_pro_database(const std::string ori_db_file_name, const std::string sp_file_name, const std::vector<std::pair<double, FastaSeqPtr>>& pro_score_vec,
		size_t db_size) {

		std::ofstream output;
		std::string sp_file_name_pre = sp_file_name.substr(0, sp_file_name.size() - 7);
		output.open(ori_db_file_name.c_str(), std::ios::out);

		for (size_t i = 0; i < db_size; i++) {
			if (i > 0 && pro_score_vec[i].first == 0)
				break;
			output << ">" << pro_score_vec[i].second->getName() << " " << pro_score_vec[i].second->getDesc() << std::endl;
			output << pro_score_vec[i].second->getRawSeq() << std::endl;
		}
		output.close();
	}

	int binarySearch(const std::vector<double>& v, double target) {
		int low = 0, n = v.size(), middle = 0;
		int high = n - 1;
		/*if ((target < v[0] - 1) && (target > v[n - 1] + 1))
			return 0;*/
		while (low < high) {
			middle = (low + high) / 2;
			if (abs(target - v[middle]) <= 1.) {
				return middle;
			}
			else if (target < v[middle] + 1) {
				high = middle;
			}
			else if (target > v[middle] + 1) {
				low = middle + 1;
			}
		}
		return 0;

	}
	//已经构造 好的截断
	double specu_error(const std::vector<double>& mod_protein, double pre_Mass) {
		double best_error = 65535;
		int pos = binarySearch(mod_protein, pre_Mass - 18.0152);
		for (size_t i = pos; i >= 0; i--)
		{
			if (std::abs(mod_protein[i] + 18.0152 - pre_Mass) < 1.) {
				if (best_error > std::abs(mod_protein[i] + 18.0152 - pre_Mass))
					best_error = std::abs(mod_protein[i] + 18.0152 - pre_Mass);
			}
			else
				break;

		}
		for (size_t i = pos; i < mod_protein.size(); i++) {
			if (std::abs(mod_protein[i] + 18.0152 - pre_Mass) < 1.) {
				if (best_error > std::abs(mod_protein[i] + 18.0152 - pre_Mass))
					best_error = std::abs(mod_protein[i] + 18.0152 - pre_Mass);
			}
			else
				break;
		}
		return best_error;
	}

	std::string& replace_all(std::string& src, const std::string& old_value, const std::string& new_value) {
		// 每次重新定位起始位置，防止上轮替换后的字符串形成新的old_value
		for (std::string::size_type pos(0); pos != std::string::npos; pos += new_value.length()) {
			if ((pos = src.find(old_value, pos)) != std::string::npos) {
				src.replace(pos, old_value.length(), new_value);
			}
			else break;
		}
		return src;
	}


	void compute_score(const FastaSeqPtr seq_info, const double pre_Mass, std::vector<std::pair<double, FastaSeqPtr>>& pro_score_vec,
		const std::vector<double>& mod_protein_mass, const std::vector<double>& trunc_mass, const size_t tag_len, const Tag_vec& ms_tag, const double err, int max_mod) {

		double prt_score = 0;
		double mod_score = 0;
		double tag_score = 0;
		StringPairVec str_pair_vec = seq_info->getAcidPtmPairVec();
		std::vector<double> ori_db_seq_vec;
		std::vector<double> mod_trunc_mass;
		int prt_str_size = str_pair_vec.size();
		int mod_size = mod_protein_mass.size();
		int trunc_size = trunc_mass.size();
		double max_prt, min_prt;
		max_prt = mod_protein_mass[mod_size - 1] + 18.0152 + 1;
		min_prt = mod_protein_mass[0] - trunc_mass[trunc_size - 1] + 18.0152 - 1;


		if ((pre_Mass < max_prt) && (pre_Mass > min_prt)) {

			//截断数目,修饰表,根据修饰和截断过滤
			for (size_t i = 0; i < trunc_size; i++)
				for (size_t j = 0; j < mod_size; j++)
					mod_trunc_mass.emplace_back(mod_protein_mass[j] - trunc_mass[i]);
			sort(mod_trunc_mass.begin(), mod_trunc_mass.end());
			double mod_error;
			double score = 0;
			double temp_mass = 0;

			mod_error = specu_error(mod_trunc_mass, pre_Mass);
			if (mod_error < 65535)
			{
				mod_score = 1 - mod_error;
			}
		}

		if (prt_str_size > tag_len) {
			if (!ms_tag.empty()) {
				std::string seq = seq_info->getString(str_pair_vec);
				std::string* re_seq = &seq;
				reverse(re_seq->begin(), re_seq->end());

				for (size_t t = 0; t < ms_tag.size(); t++)
				{
					if (strstr(seq.c_str(), ms_tag[t].first.c_str()) || strstr(re_seq->c_str(), ms_tag[t].first.c_str()))
						tag_score += (1 - ms_tag[t].second);
				}
				tag_score = tag_score * log2(prt_str_size) / prt_str_size;
			}
		}
		prt_score = mod_score + tag_score;
		pro_score_vec.emplace_back(std::make_pair(prt_score, seq_info));
	}

	bool GreaterSort(std::pair<double, FastaSeqPtr> a, std::pair<double, FastaSeqPtr> b) { return (a.first > b.first); }

	void get_filter_db(const DeconvMsPtr deconv_ms_ptr, const std::vector<std::vector<double>>& mod_protein_mass, std::vector<std::vector<double>*>& all_Trunc_mass,
		TagGenerate& tag_gen, const std::string sp_file_name, const std::string ori_db_file_name,
		const int tag_len, const size_t db_size, const double err, int max_mod)
	{
		std::vector<std::pair<double, FastaSeqPtr>> pro_score_vec;
		std::vector<double> d_ms_seq_vector;
		//质谱质量列表
		MsHeaderPtr header_ptr = deconv_ms_ptr->getMsHeaderPtr();
		double pre_Mass = header_ptr->getPrecMonoMass();
		int scan_id = header_ptr->getFirstScanNum();
		int sp_id = header_ptr->getId();

		std::vector<double> ms_vector;
		get_One_Ms_vec(deconv_ms_ptr, ms_vector);
		//生成质谱tag
		Tag_vec ms_tag = tag_gen.Gen_tag(ms_vector);

		//std::cout << "get_Ms_seq - ended. d_ms_tag_vector size : " << ms_tag.size() << std::endl;
		std::string cstr_ori_database_backup = "ori_database";
		FastaReader reader(cstr_ori_database_backup);
		FastaSeqPtr seq_info = reader.getNextSeq();
		int temp = 0;
		//i定位protein_mass
		int i = 0;
		while (seq_info != nullptr) {
			compute_score(seq_info, pre_Mass, pro_score_vec, mod_protein_mass[i], *(all_Trunc_mass[i]), tag_len, ms_tag, err, max_mod);
			seq_info = reader.getNextSeq();
			i++;
		}

		//将所有蛋白质序列的E-score得分进行排序
		std::sort(pro_score_vec.begin(), pro_score_vec.end(), GreaterSort);

		//得到ASF方法的输入蛋白质序列数据库
		gen_pro_database(ori_db_file_name, sp_file_name, pro_score_vec, db_size);
	}
}

	