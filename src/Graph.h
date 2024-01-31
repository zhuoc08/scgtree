#pragma once

#include<vector>
#include<set>
#include<unordered_map>

#include<metis.h>

#include"Vertex.h"
#include"Matrix.h"



class Graph
{
    public:
        int type=0; // graph type, 0-undirected , 1-directed

        int nov=0, noe=0; // number of vertex/edge
        std::vector<Vertex> vertices, rvertices; // vertices:out-edges, rvertices:in-edges

    public:
        Graph();

        Graph(int n);

        bool load(const char* coordinate_file, const char* edge_file);

        void load_coordinate(const char* coordinate_file);

        void load_edge(const char* edge_file);

        bool edge_exist(int s, int d);

        void make_rvertices();

        std::vector<int> dijkstra( int s, std::vector<int> &cands);

        int dijkstra_visited_edges( int s, std::vector<int> &cands);

        std::vector<int> dijkstra_with_path(int s, int d, std::vector<int> &path);

        int bi_dijkstra( int s, int d);

        Matrix floyd();

        double rad(double d);

        double heuristic_function(int s, int d);

        int a_star_algorithm(int s, int d);

        bool edge_weight_update(int s, int t, int weight);


    public:
        bool adjweight_set_to_all_one = true;

        // use for metis
        // idx_t = int64_t / real_t = double
        idx_t nvtxs; // |vertices|
        idx_t ncon; // number of weight per vertex
        idx_t* xadj; // array of adjacency of indices
        idx_t* adjncy; // array of adjacency nodes
        idx_t* vwgt; // array of weight of nodes
        idx_t* adjwgt; // array of weight of edges in adjncy
        idx_t nparts; // number of parts to partition
        idx_t objval; // edge cut for partitioning solution
        idx_t* part; // array of partition vector
        idx_t options[METIS_NOPTIONS]; // option array

        void metis_options_setting();

        std::unordered_map<int,int> graph_partition(int fanout, std::set<int> &nset );



};
