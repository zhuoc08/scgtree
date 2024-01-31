#pragma once


#include<set>
#include<unordered_map>
#include<utility>




class SpspResult{
	public:
		int dis;
		std::vector<int> path;
		std::vector<int> path_dis;
};


// knn_query result
class ResultSet{
	public:
		int id;
		int dis;
};



// init search node
class QueryStatus{
	public:
		int id;
		bool isvertex;
		int dis;
};



// init status struct
class Status{
	public:
		int tnid; // tree node id
		std::set<int> nset; // node set
};




class PathRecoveryResult{
	public:
		int status; // 0:unfind, 1:find, 2:v1/v2 is not in union border list
		int v1, v2;
		int b;
		int sp_dis1, sp_dis2;
};





bool cmp(std::pair<long long, int> a, std::pair<long long, int> b);



class query_status_competor{
	public:
	bool operator()( const QueryStatus& l, const QueryStatus& r ){
		return l.dis > r.dis;
	}
};

