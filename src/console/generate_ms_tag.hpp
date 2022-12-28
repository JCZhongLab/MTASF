
#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<vector>
#include<map>
#include<algorithm>
#define MAX_EDGE_MASS 200

#ifndef GENERATE_MS_TAG_H
#define GENERATE_MS_TAG_H
namespace prot {
    class Tag {
    public:
        Tag(int len = 5) {
            for (int i = 0; i < len; i++)
                this->sequence_letter.push_back('A');//Initial
            this->length = 0;
            this->tag_error = 0;
        }
        std::string sequence_letter;
        int length;
        double tag_error;
    };


    typedef struct Edge_ACID
    {
        char amino_acid = 'A';
        double error = 0.0;
        int from_peak = 0;
        int to_peak = 0;
    }Edge_acid;

    typedef std::vector < std::pair<std::string, double> > Tag_vec;
    typedef std::pair< Tag_vec, Tag_vec> Pair_tag_comb;


    class TagGenerate {

    public:
        TagGenerate(int tag_len, double error);
        Tag_vec Gen_tag(std::vector<double>& ms_vec);
        std::vector<Edge_acid*> mass_difference_to_amino_acid(double mass_dif, double this_peak, int from_peak, int to_peak);
        void dfs_search(int current_pos, Tag tag);
        void remove_same_tag(Tag_vec& tag_vec);
        double error;
        int tag_len;
        std::vector<Tag> tag_list;
        std::vector<std::pair<char, double>> Amino2Mass;
        std::vector<std::vector<Edge_acid*>> edges;
        std::vector<std::vector<Edge_acid*>> edges_topeak;
    };
}
#endif