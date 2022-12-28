#pragma once
#include <string>
#include <map>
#include <vector>
#include "console/toppic_process.hpp"
#include "console/generate_ms_tag.hpp"

#ifndef FILTER_DB_H
#define FILTER_DB_H
namespace prot {

	//get precursor
	double get_Ms_seq(const DeconvMsPtr deconv_ms_ptr, std::vector<double>& d_vector);

	//get the protein theoretical mass
	void get_protein_mass(const std::string cstr_ori_database_backup, std::vector<double>& protein_mass,std::vector<std::vector<double>*>& all_Trunc_mass, double trunc_ratio);

	//get modification mass
	std::vector<double> get_modlist(std::string mod_name, int max_mod);

	//get proteoform mass set
	std::vector<std::vector<double>> get_mod_protein_mass(const std::vector<double>& protein_mass,std::vector<double> &mod_list);

	//get ms_vec
	void get_One_Ms_vec(const DeconvMsPtr deconv_ms_ptr, std::vector<double>& ms_vector);

	//write filter_file
	void gen_pro_database(const std::string ori_db_file_name,  const std::string sp_file_name, const std::vector<std::pair<double, FastaSeqPtr>> &pro_score_vec,
		size_t db_size);
	
	// MF_TF score
	void compute_score(const FastaSeqPtr seq_info, const double pre_Mass,std::vector<std::pair<double, FastaSeqPtr>>& pro_score_vec,
		const std::vector<double>& mod_protein_mass, const std::vector<double>& trunc_mass,
		const size_t tag_len, const Tag_vec& ms_tag,const double err, int max_mod);

	//MF_TF
	void get_filter_db(const DeconvMsPtr deconv_ms_ptr, const std::vector<std::vector<double>>& mod_protein_mass, std::vector<std::vector<double>*>& all_Trunc_mass,
		TagGenerate& tag_gen, const std::string sp_file_name, const std::string ori_db_file_name,
		const int tag_len, const size_t db_size, const double err, int max_mod);
}
#endif