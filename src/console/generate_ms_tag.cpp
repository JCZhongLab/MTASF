#include "generate_ms_tag.hpp"



namespace prot {
    TagGenerate::TagGenerate(int tag_len, double error) {
        this->error = error;
        this->tag_len = tag_len;
        this->Amino2Mass.push_back(std::pair<char, double>('B', 0.0));
        this->Amino2Mass.push_back(std::pair<char, double>('G', 57.02146));
        this->Amino2Mass.push_back(std::pair<char, double>('A', 71.03711));
        this->Amino2Mass.push_back(std::pair<char, double>('S', 87.03203));
        this->Amino2Mass.push_back(std::pair<char, double>('P', 97.05276));
        this->Amino2Mass.push_back(std::pair<char, double>('V', 99.06841));
        this->Amino2Mass.push_back(std::pair<char, double>('T', 101.0477));
        this->Amino2Mass.push_back(std::pair<char, double>('C', 103.0092));
        this->Amino2Mass.push_back(std::pair<char, double>('I', 113.0841));
        this->Amino2Mass.push_back(std::pair<char, double>('L', 113.0841));
        this->Amino2Mass.push_back(std::pair<char, double>('N', 114.0429));
        this->Amino2Mass.push_back(std::pair<char, double>('D', 115.0296));
        this->Amino2Mass.push_back(std::pair<char, double>('Q', 128.0586));
        this->Amino2Mass.push_back(std::pair<char, double>('K', 128.095));
        this->Amino2Mass.push_back(std::pair<char, double>('E', 129.0426));
        this->Amino2Mass.push_back(std::pair<char, double>('M', 131.0405));
        this->Amino2Mass.push_back(std::pair<char, double>('H', 137.0589));
        this->Amino2Mass.push_back(std::pair<char, double>('F', 147.0684));
        this->Amino2Mass.push_back(std::pair<char, double>('R', 156.101));
        this->Amino2Mass.push_back(std::pair<char, double>('Y', 163.0633));
        this->Amino2Mass.push_back(std::pair<char, double>('W', 186.0793));
    }

    bool cmp(std::pair<std::string, double>& p1, std::pair<std::string, double>& p2) {
        return p1.second < p2.second;
    }

    Tag_vec TagGenerate::Gen_tag(std::vector<double>& ms_vec)
    {
        //std::vector<std::vector<Edge_acid *>> edges;
        std::vector<Edge_acid*> edge_list;
        Tag_vec tag_vec;

        if (!this->tag_list.empty())
            std::vector<Tag>().swap(this->tag_list);

        if (!this->edges.empty())
            std::vector<std::vector<Edge_acid*>>().swap(this->edges);

        if (!this->edges_topeak.empty())
            std::vector<std::vector<Edge_acid*>>().swap(this->edges_topeak);

        this->edges.resize(ms_vec.size());
        this->edges_topeak.resize(ms_vec.size());
        double mass_dif;
        if(ms_vec.size() < 2)
            return tag_vec;
        for (int i = 0; i < ms_vec.size() - 1; i++) {
            for (int j = i + 1; j < ms_vec.size(); j++)
            {
                mass_dif = ms_vec[j] - ms_vec[i];
                if (mass_dif > MAX_EDGE_MASS)
                    break;
                //Amino2Mass的最大值
                if (mass_dif < (this->Amino2Mass.end() - 1)->second + 1) {
                    edge_list = this->mass_difference_to_amino_acid(mass_dif, ms_vec[j], i, j);
                    if (!edge_list.empty())
                        for (std::vector<Edge_acid*>::iterator e = edge_list.begin(); e != edge_list.end(); e++) {
                            this->edges[i].push_back(*e);
                            this->edges_topeak[j].push_back(*e);
                        }
                }
            }
        }

        for (size_t i = 0; i < ms_vec.size(); i++) {
            Tag tag(this->tag_len);
            this->dfs_search(i, tag);
        }
        //去重
        this->remove_same_tag(tag_vec);

        return tag_vec;
    }

    //remove duplicate tag
    void TagGenerate::remove_same_tag(Tag_vec& tag_vec) {

        bool flag_isotope;
        Tag_vec tag_vec_1;
        int tag_num = 0;
        std::map<std::string, bool> temp_map;
        for (std::vector<Tag>::iterator t = this->tag_list.begin(); t != this->tag_list.end();)
        {
            flag_isotope = false;
            std::string temp_s = "";
            for (int i = 0; i < t->sequence_letter.size(); i++)
            {
                if (t->sequence_letter[i] == 'B') {
                    flag_isotope = true;
                    t = this->tag_list.erase(t);
                    break;
                }
            }
            if (!flag_isotope) {
                tag_vec_1.push_back(std::make_pair(t->sequence_letter, t->tag_error));
                t++;
            }
        }
        if (tag_vec_1.empty()) return;

        sort(tag_vec_1.begin(), tag_vec_1.end(), cmp);
        for (size_t i = 0; i < tag_vec_1.size(); i++) 
            if (temp_map.count(tag_vec_1[i].first) == 0) {
                temp_map[tag_vec_1[i].first] = true;
                tag_vec.push_back(tag_vec_1[i]);
                tag_num++;
            }            
    }

    std::vector<Edge_acid*> TagGenerate::mass_difference_to_amino_acid(double mass_dif, double this_peak, int from_peak, int to_peak) {
        std::vector<Edge_acid*> edge_list;
        double a_diff;
        double temp_e= this->error;
        char temp_c;
        for (std::vector<std::pair<char, double>>::iterator acid = this->Amino2Mass.begin(); acid != this->Amino2Mass.end(); acid++) {
            a_diff = abs(mass_dif - acid->second);

            //if(a_diff/ acid->second * 1000000<=  this->error)
            if (a_diff < this->error)
            {
                if (a_diff < temp_e) {
                    temp_c = acid->first;
                    temp_e = a_diff;
                }    
                //break;
            }
            //if (mass_dif+1 < acid->second) break;
        }
        if (temp_e < this->error) {
            Edge_acid* e = new(Edge_acid);
            e->amino_acid = temp_c;
            e->error = temp_e;
            e->from_peak = from_peak;
            e->to_peak = to_peak;
            edge_list.push_back(e);
            if (temp_c == 'I') {
                Edge_acid* e = new(Edge_acid);
                e->amino_acid = 'L';
                e->error = temp_e;
                e->from_peak = from_peak;
                e->to_peak = to_peak;
                edge_list.push_back(e);
            }
        }
        return edge_list;
    }



    void TagGenerate::dfs_search(int current_pos, Tag tag) {
        if (tag.length == this->tag_len) {
            tag.tag_error /= this->tag_len;
            this->tag_list.push_back(tag);
            return;
        }

        for (std::vector<Edge_acid*>::iterator e = this->edges[current_pos].begin(); e != this->edges[current_pos].end(); e++) {
            tag.sequence_letter[tag.length] = (*e)->amino_acid;
            tag.length += 1;
            tag.tag_error += (*e)->error;
            this->dfs_search((*e)->to_peak, tag);
            tag.length -= 1;
            tag.tag_error -= (*e)->error;
        }
    }
}