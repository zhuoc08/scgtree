#pragma once


#include<vector>
#include<map>
#include<unordered_map>
#include"Matrix.h"



class TreeNode
{
    public:
        std::vector<int> borders;
        std::vector<int> children;
        bool isleaf;
        std::vector<int> leafnodes;
        int father;
        int in_father_index;
        int level;	
        int vertexs_num;  

        std::vector<int> border_in_father, border_in_union_border;
        std::unordered_map<int, int> vertex_in_leafnodes;
        std::unordered_map<int, int> border_in_borders_map, union_border_in_union_borders_map;

        std::vector<int> union_borders; 
        std::vector<int> mind; 

        Matrix dist;
        std::unordered_map<int, int> short_cuts;
        std::vector<int> nonleafinvlist;
        std::vector<int> leafinvlist;
        std::vector<std::pair<int, int> > short_cuts_dis;
};

