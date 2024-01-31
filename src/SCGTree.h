#pragma once


#include<vector>
#include<queue>

#include"Graph.h"
#include"TreeNode.h"
#include"AuxiliaryData.h"
#include"Matrix.h"


class SCGTree
{
    public:

        Graph graph;

        int fanout, tau;
        std::vector<TreeNode> tree;
        int treesize=0, treedepth=0;

		std::vector<bool> isborder_set;
		std::vector<std::vector<int> > treepath_set;

        std::vector<int> border_max_level, border_max_level_index;

        int flevel, nflevel;
        std::vector<std::pair<int, int> > short_cut_pairs_in_level;
        std::vector<Matrix> short_cut_matrix_in_level;

        SCGTree(Graph g);


        void build(int fanout, int tau, int lamda);

        void build_tree_skeleton(int fanout, int tau);

        void calculate_distance_matrix_up_and_down();

        void calculate_shortcut_matrix_with_given_level(int level);

        void calculate_shortcut_matrix_with_given_level(int level, 	std::vector< std::vector<int> > &all_minds, std::vector<std::pair<int, int> > &all_pairs);

        int cal_shortcut_value(int tnode1, int tnode2);

        void cal_temp_dist(int tnode, std::map<int, std::vector< std::vector<int> > > &temp_map);

        void calculate_shortcut_matrix_with_lamda(int lamda);

        void cal_short_cuts_dis();

        void calculate_auxiliary_data();

        int spsp_query(int s, int d);

        void cal_candidate_nodes(std::vector<int> &objects);

        inline void expand_search_area_knn(int locid, int &expandable, int &tn, int &tmin, std::priority_queue<QueryStatus, std::vector<QueryStatus >, query_status_competor > &pq, std::unordered_map<int, std::vector<int> > &itm, std::vector<bool> &shortcuts_visited, std::vector<int> &shortcuts_real_distance);

        std::vector<ResultSet> knn_query( int locid, int K);

        void expand_search_area_range(int locid, int range, int &expandable, int &tn, int &tmin, std::queue<QueryStatus> &q, std::unordered_map<int, std::vector<int> > &itm, std::vector<bool> &shortcuts_visited, std::vector<int> &shortcuts_real_distance, std::vector<ResultSet> rstset);

        std::vector<ResultSet> range_query( int locid, int range);

        void save(const char *scgtree_file, bool compress, bool dump_path);

        void load(const char *scgtree_file, bool compress, bool dump_path);

        bool edge_weight_update(int s, int d, int weight, bool show_info, bool with_shortcut);

        bool update_shortcut(std::set<int> &updated_shortcuts);

        int find_lca(int s, int d);

        int find_lca_level(int s, int d);

        SpspResult spsp_path_query(int s, int d, bool with_shortcut);

        std::vector<int> tree_path_recovery(int s, int d, SpspResult sr);

        std::vector<int> path_recovery(int v1, int v2, int d, std::vector<int> &path);

        PathRecoveryResult find_border(int v1, int v2,  int d, int tn);

};