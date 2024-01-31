#include<iostream>
#include<sys/time.h>
#include<fstream>
#include<sstream>
#include<math.h>
#include<algorithm>

#include"SCGTree.h"
#include"Graph.h"




void build_and_save_test(std::string dist_file, std::string co_file, std::string scg_file, int fanout, int tau, int lamda) {

    struct timeval tv;
    long long ts, te;
    long long t1, t2, total_time;

    Graph graph;
    graph.load(co_file.c_str(), dist_file.c_str());

    std::cout << "load graph done!" << std::endl;

    SCGTree tree(graph);


    gettimeofday( &tv, NULL );
    t1 = tv.tv_sec;

	tree.build(fanout, tau, lamda); 

    gettimeofday( &tv, NULL );
    t2 = tv.tv_sec;
    total_time = t2 - t1;
    std::cout << "build scgtree time: " << total_time << " s." <<std::endl;


	tree.save(scg_file.c_str(), false, false);


    // spsp query test
    // int s, d, result;
    // while(std::cin >> s >> d) {
	// 	gettimeofday( &tv, NULL );
	// 	t1 = tv.tv_sec * 1000000 + tv.tv_usec ;

    //     result = tree.spsp_query(s, d);

	// 	gettimeofday( &tv, NULL );
	// 	t2 = tv.tv_sec * 1000000 + tv.tv_usec ;
	// 	total_time = t2 - t1;
    //     std::cout << "scgtree result: " << result << ", time: " << total_time << " us." <<std::endl;
    // }



}



void spsp_query_test(std::string dist_file, std::string co_file, std::string scg_file, int query_time){

    Graph graph;
    graph.load(co_file.c_str(), dist_file.c_str());
    // std::cout << "load graph done!" << std::endl;

    SCGTree tree(graph);
    tree.load(scg_file.c_str(), false, false);
    // std::cout << "load tree done!" << std::endl;

	int s, d;
    std::vector<int> svertexs, dvertexs;

    while(svertexs.size()<query_time){
        s = rand()%tree.graph.nov;
        d = rand()%tree.graph.nov;

        svertexs.push_back(s);
        dvertexs.push_back(d);
    }

    struct timeval tv;
    long long ts, te, t1, t2, total_time;

    gettimeofday( &tv, NULL );
    t1 = tv.tv_sec * 1000000 + tv.tv_usec ;

    for(int i=0; i<query_time; i++){
        s = svertexs[i];
        d = dvertexs[i];
        tree.spsp_query(s, d);
    }

    gettimeofday( &tv, NULL );
    t2 = tv.tv_sec * 1000000 + tv.tv_usec ;
    total_time = t2 - t1;

    std::cout << "avg time of " << query_time << " times spsp query: " << ((long double)total_time)/query_time << "(us)" << std::endl;

}



void knn_query_test(std::string dist_file, std::string co_file, std::string scg_file, int query_time){


    Graph graph;
    graph.load(co_file.c_str(), dist_file.c_str());
    // std::cout << "load graph done!" << std::endl;

    SCGTree tree(graph);
    tree.load(scg_file.c_str(), false, false);
    // std::cout << "load tree done!" << std::endl;

    int k=10, vq, c, C=0.001*graph.nov;
	std::vector<int> vqs, objects;

    // v_q
    vqs.clear();
    for(int i=0; i<query_time; i++) {
        vq = rand()%graph.nov;
        vqs.push_back(vq);
    }

    // candidate vertices
    while(objects.size() != C){
        c = rand()%graph.nov;
        if(find(objects.begin(), objects.end(), c) == objects.end()){
            objects.push_back(c);
        }
    }
    tree.cal_candidate_nodes(objects); 


    struct timeval tv;
    long long ts, te, t1, t2, total_time;

    gettimeofday( &tv, NULL );
    t1 = tv.tv_sec * 1000000 + tv.tv_usec ;

    for(int i=0; i<query_time; i++) {
        tree.knn_query(vqs[i], k);
    }

    gettimeofday( &tv, NULL );
    t2 = tv.tv_sec * 1000000 + tv.tv_usec ;
    total_time = t2 - t1;

    std::cout << "avg time of " << query_time << " times knn query: " << ((long double)total_time)/query_time << "(us)" << std::endl;

}


void range_query_test(std::string dist_file, std::string co_file, std::string scg_file, int query_time){


    Graph graph;
    graph.load(co_file.c_str(), dist_file.c_str());
    // std::cout << "load graph done!" << std::endl;

    SCGTree tree(graph);
    tree.load(scg_file.c_str(), false, false);
    // std::cout << "load tree done!" << std::endl;

    int d=5000, vq, c, C=0.001*graph.nov;
	std::vector<int> vqs, objects;

    // v_q
    vqs.clear();
    for(int i=0; i<query_time; i++) {
        vq = rand()%graph.nov;
        vqs.push_back(vq);
    }

    // candidate vertices
    while(objects.size() != C){
        c = rand()%graph.nov;
        if(find(objects.begin(), objects.end(), c) == objects.end()){
            objects.push_back(c);
        }
    }
    tree.cal_candidate_nodes(objects); 

    struct timeval tv;
    long long ts, te, t1, t2, total_time;

    gettimeofday( &tv, NULL );
    t1 = tv.tv_sec * 1000000 + tv.tv_usec ;

    for(int i=0; i<query_time; i++) {
        tree.range_query(vqs[i], d);
    }

    gettimeofday( &tv, NULL );
    t2 = tv.tv_sec * 1000000 + tv.tv_usec ;
    total_time = t2 - t1;

    std::cout << "avg time of " << query_time << " times range query: " << ((long double)total_time)/query_time << "(us)" << std::endl;

}


int main(){

	std::string dist_file="BAY.d", co_file="BAY.co", scg_file="BAY.scg";

    int fanout = 4;
    int tau = 128;
    int lamda = 20; // ratio of total size of distance matirces
    int query_time = 100000;

    build_and_save_test(dist_file, co_file, scg_file, fanout, tau, lamda);

    spsp_query_test(dist_file, co_file, scg_file, query_time);

    knn_query_test(dist_file, co_file, scg_file, query_time);

    range_query_test(dist_file, co_file, scg_file, query_time);

    return 0;
}