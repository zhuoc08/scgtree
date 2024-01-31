#include<iostream>
#include<metis.h>
#include<unordered_map>
#include<queue>
#include<set>
#include<map>
#include<algorithm>
#include<limits.h>
#include<math.h>


#include"Vertex.h"
#include"Graph.h"
#include"MinHeap.h"
#include"Matrix.h"


Graph::Graph(){
    metis_options_setting();

}


Graph::Graph(int n){
    metis_options_setting();
	vertices.resize(n);
	rvertices.resize(n);
}


bool Graph::load(const char* coordinate_file, const char* edge_file){


	FILE *fin;
	fin = fopen(coordinate_file, "r");
	if(!fin) {
		std::cout << "coordinate file not exist!" << std::endl;
		return false;
	}

	char tag;
	int vertexNum, num;
	
	double xpos, ypos, x=0, y=0;

	// read first line
	fscanf(fin, "%d", &vertexNum );

	nov = vertexNum;
	vertices.resize(nov, Vertex());
	rvertices.resize(nov, Vertex());

	while( fscanf(fin, "%d %lf %lf", &num, &xpos, &ypos ) != EOF ){
		vertices[num-1].x = xpos;
		vertices[num-1].y = ypos;

		rvertices[num-1].x = xpos;
		rvertices[num-1].y = ypos;
	}


	fclose(fin);


	// load edge
	fin = fopen(edge_file, "r");
	if(!fin) {
		std::cout << "edge file not exist!" << std::endl;
		return false;
	}

	// init nodes
	int edgeNum;

	// read first line	
	fscanf(fin, "%d %d", &vertexNum, &edgeNum );

	noe = edgeNum;

	int eid;
	int snid, enid;
	double weight;
	int iweight;

	while( fscanf(fin,"%d %d %lf", &snid, &enid, &weight ) == 3 ){
		iweight = (int) weight;
		vertices[snid-1].adjnodes.push_back( enid-1 );
		vertices[snid-1].adjweight.push_back(iweight );

		rvertices[enid-1].adjnodes.push_back( snid-1 );
		rvertices[enid-1].adjweight.push_back(iweight );

	}

	fclose(fin);


	//sort adjnodes
	std::unordered_map<int, int> edge_weight;
	for(int i=0; i<vertexNum; i++){
		edge_weight.clear();
		for(int j=0; j<vertices[i].adjnodes.size(); j++){
			edge_weight[vertices[i].adjnodes[j]] = vertices[i].adjweight[j];

		}
		std::sort(vertices[i].adjnodes.begin(), vertices[i].adjnodes.end());

		vertices[i].adjweight.clear();
		for(int j=0; j<vertices[i].adjnodes.size(); j++){
			vertices[i].adjweight.push_back(edge_weight[vertices[i].adjnodes[j]]);
		}

	}



	return true;

}


void Graph::load_coordinate(const char* coordinate_file){



	if(noe==0){
		std::cout << "edges should be loaded first!\n" << std::endl;
		return;
	}


	FILE *fin;
	fin = fopen(coordinate_file, "r");
	if(!fin) {
		std::cout << "coordinate file not exist!" << std::endl;
		return;
	}

	std::cout << "LOADING NODE COORDINATE... ";

	char tag;
	int vertexNum, num;
	
	double xpos, ypos, x=0, y=0;

	// read first line
	fscanf(fin, "%d", &vertexNum );

	if(nov == 0) {
		nov = vertexNum;
		for(int i=0; i<vertexNum; i++){
			Vertex vertex = { x, y };
			vertices.push_back(vertex);
		}

	} else if(vertexNum != nov) {
		std::cout << "node file and edge file are not matched!" << std::endl;
		return;
	}



	while( fscanf(fin, "%d %lf %lf", &num, &xpos, &ypos ) != EOF ){
		vertices[num-1].x = xpos;
		vertices[num-1].y = ypos;
	}


	fclose(fin);

	std::cout << " done!" << std::endl;


}


void Graph::load_edge(const char* edge_file){

	FILE *fin;


	fin = fopen(edge_file, "r");
	if(!fin) {
		std::cout << "edge file not exist!" << std::endl;
		return;
	}

	// init nodes
	int vertexNum, edgeNum;
	double x=0, y=0;

	// read first line	
	fscanf(fin, "%d %d", &vertexNum, &edgeNum );
	if(nov == 0) {

		nov = vertexNum;
		for(int i=0; i<vertexNum; i++){
			Vertex vertex = { x, y };
			vertices.push_back(vertex);
		}

	} else if(vertexNum != nov) {
		std::cout << "node file and edge file are not matched!" << std::endl;
		return;
	}


	if(noe != 0) {
		std::cout << "edge are already loaded!" << std::endl;
		return;
	} else {
		noe = edgeNum;
	}


	std::cout << "NODE_COUNT=" << nov << std::endl;
	std::cout << "EDGE_COUNT=" << noe << std::endl;




	// load edge
	std::cout << "LOADING EDGE... ";
	int eid;
	int snid, enid;
	double weight;
	int iweight;

	while( fscanf(fin,"%d %d %lf", &snid, &enid, &weight ) == 3 ){
		iweight = (int) weight;
		vertices[snid-1].adjnodes.push_back( enid-1 );
		vertices[snid-1].adjweight.push_back(iweight );
	}

	fclose(fin);

	std::cout << " done!" << std::endl;


}


// check if edge exist
bool Graph::edge_exist(int s, int d){
	return std::binary_search(vertices[s].adjnodes.begin(), vertices[s].adjnodes.end(), d);
}


// make rvertices
void Graph::make_rvertices(){
	nov = vertices.size();
	noe = 0;
	int adjv, adjw;

	for(int i=0; i<rvertices.size(); i++) {
		rvertices[i].adjnodes.clear();
		rvertices[i].adjweight.clear();
	}

	for(int i=0; i<vertices.size(); i++){
		for(int j=0; j<vertices[i].adjnodes.size(); j++){
			adjv = vertices[i].adjnodes[j];
			adjw = vertices[i].adjweight[j];
			rvertices[adjv].adjnodes.push_back(i);
			rvertices[adjv].adjweight.push_back(adjw);
			noe ++;
		}
	}
}


// dijkstra algorithm, s = source node, cands = candidate nodes
std::vector<int> Graph::dijkstra( int s, std::vector<int> &cands){
	// init
	std::set<int> todo;
	todo.clear();
	todo.insert(cands.begin(), cands.end());

	std::unordered_map<int,int> result;
	result.clear();
	std::set<int> visited;
	visited.clear();

	MinHeap Q;

	Q.insert(std::make_pair(s, 0));
	int adjnode, weight;
	while(!todo.empty() && !Q.empty()) {
		std::pair<int, int> top = Q.extract_top();

		// put min to result, add to visited
		result[top.first] = top.second;
		visited.insert( top.first );

		// unvisited vertex
		if ( todo.find( top.first ) != todo.end() ){
			todo.erase( top.first );
		}

		// expand
		for ( int i = 0; i < vertices[top.first].adjnodes.size(); i++ ){
			adjnode = vertices[top.first].adjnodes[i];
			if ( visited.find( adjnode ) != visited.end() ){
				continue;
			}

			weight = vertices[top.first].adjweight[i];

			int pos = Q.find(adjnode);
			if ( pos == -1 ){
				Q.insert(std::make_pair(adjnode, top.second+weight));
			} else{
				if ( top.second+weight < Q.get(pos).second ){
					Q.decrease(pos, std::make_pair( adjnode, top.second+weight));
				}
			}
		}
	}

	std::vector<int> output;
	for ( int i = 0; i < cands.size(); i++ ){
		output.push_back( result[cands[i]] );
	}

	return output;
}


int Graph::dijkstra_visited_edges( int s, std::vector<int> &cands){

	int visited_edges = 0;

	// init
	std::set<int> todo;
	todo.clear();
	todo.insert(cands.begin(), cands.end());

	std::unordered_map<int,int> result;
	result.clear();
	std::set<int> visited;
	visited.clear();

	MinHeap Q;

	Q.insert(std::make_pair(s, 0));
	int adjnode, weight;
	while(!todo.empty() && !Q.empty()) {
		std::pair<int, int> top = Q.extract_top();
		
		// put min to result, add to visited
		result[top.first] = top.second;
		visited.insert( top.first );

		if ( todo.find( top.first ) != todo.end() ){
			todo.erase( top.first );
		}

		// expand
		for ( int i = 0; i < vertices[top.first].adjnodes.size(); i++ ){
			adjnode = vertices[top.first].adjnodes[i];
			if ( visited.find( adjnode ) != visited.end() ){
				continue;
			}

			visited_edges += 1;

			weight = vertices[top.first].adjweight[i];

			int pos = Q.find(adjnode);

			if ( pos == -1 ){
				Q.insert(std::make_pair(adjnode, top.second+weight));
			} else{
				if ( top.second+weight < Q.get(pos).second ){
					Q.decrease(pos, std::make_pair( adjnode, top.second+weight));
				}
			}
		}
	}

	return visited_edges;
}


// obtain shortest path from s to d
std::vector<int> Graph::dijkstra_with_path(int s, int d, std::vector<int> &path){
	// init
	if(s==d) {
		return std::vector<int>{s, d};
	}

	std::unordered_map<int,int> result, pred;
	result.clear();
	std::set<int> visited;
	visited.clear();

	MinHeap Q;
	pred[s] = -1;
	Q.insert(std::make_pair(s, 0));
	int adjnode, weight;
	while(!Q.empty()) {
		std::pair<int, int> top = Q.extract_top();

		if(top.first==d){
			break;
		}
		
		// put min to result, add to visited
		result[top.first] = top.second;
		visited.insert( top.first );


		// expand
		for ( int i = 0; i < vertices[top.first].adjnodes.size(); i++ ){
			adjnode = vertices[top.first].adjnodes[i];
			if ( visited.find( adjnode ) != visited.end() ){
				continue;
			}

			weight = vertices[top.first].adjweight[i];
			
			int pos = Q.find(adjnode);

			if ( pos == -1 ){
				Q.insert(std::make_pair(adjnode, top.second+weight));
				pred[adjnode] = top.first;
			} else{
				if ( top.second+weight < Q.get(pos).second ){
					Q.decrease(pos, std::make_pair( adjnode, top.second+weight));
					pred[adjnode] = top.first;
				}
			}
		}
	}


	std::vector<int> path_temp;
	path_temp.push_back(d);

	int v_pred = pred[d];

	while(v_pred != -1){
		path_temp.push_back(v_pred);
		v_pred = pred[v_pred];
	}

	std::reverse(path_temp.begin(), path_temp.end());
	path.insert(path.end(), path_temp.begin(), path_temp.end());

	return path;
}


// bi_dijkstra
int Graph::bi_dijkstra( int s, int d){

	if(s==d) {
		return 0;
	}

	std::unordered_map<int,int> fresult, bresult;
	std::set<int> fvisited, bvisited;

	int miu = INT_MAX;

	MinHeap qf, qb;

	qf.insert(std::make_pair(s, 0));
	qb.insert(std::make_pair(d, 0));

	int adjnode, weight, relax_weight;
	while(!qf.empty() && !qb.empty()) {

		std::pair<int, int> ftop = qf.extract_top();
		std::pair<int, int> btop = qb.extract_top();

		fresult[ftop.first] = ftop.second;
		bresult[btop.first] = btop.second;

		fvisited.insert( ftop.first );
		bvisited.insert( btop.first );

		// expand and relax in qf
		for ( int i = 0; i < vertices[ftop.first].adjnodes.size(); i++ ){
			adjnode = vertices[ftop.first].adjnodes[i];
			if ( fvisited.find( adjnode ) != fvisited.end() ){
				continue;
			}

			weight = vertices[ftop.first].adjweight[i];
			
			int pos = qf.find(adjnode);
			if ( pos == -1 ){
				qf.insert(std::make_pair(adjnode, ftop.second+weight));
				relax_weight = ftop.second+weight;
			} else{
				if ( ftop.second+weight < qf.get(pos).second ){
					qf.decrease(pos, std::make_pair( adjnode, ftop.second+weight));

					relax_weight = ftop.second+weight;
				} else {
					relax_weight = qf.get(pos).second;
				}
			}

			if(bvisited.find(adjnode) != bvisited.end() && miu > relax_weight + bresult[adjnode]) {
				miu = relax_weight + bresult[adjnode];
			}
		}

		// expand and relax in qb
		for ( int i = 0; i < rvertices[btop.first].adjnodes.size(); i++ ){
			adjnode = rvertices[btop.first].adjnodes[i];
			if ( bvisited.find( adjnode ) != bvisited.end() ){
				continue;
			}

			weight = rvertices[btop.first].adjweight[i];

			int pos = qb.find(adjnode);

			if ( pos == -1 ){
				qb.insert(std::make_pair(adjnode, btop.second+weight));
				relax_weight = btop.second+weight;
			} else{
				if ( btop.second+weight < qb.get(pos).second ){
					qb.decrease(pos, std::make_pair( adjnode, btop.second+weight));
					relax_weight = btop.second+weight;
				} else {
					relax_weight = qb.get(pos).second;
				}
			}

			if(fvisited.find(adjnode) != fvisited.end() && miu > relax_weight + fresult[adjnode]) {
				miu = relax_weight + fresult[adjnode];
			}
		}

		if(qf.top().second+qb.top().second > miu) {
			break;
		}
	}

	return miu;

}


// floyd algorithm
Matrix Graph::floyd(){

	// initial matrix
	Matrix dist;
	dist.init(nov, nov, std::vector<int>(nov*nov, INT_MAX));
	for(int i=0; i<nov; i++) {
		dist.a[i][i] = 0;
	}

	for(int i=0; i<nov; i++) {
		for(int j=0; j<vertices[i].adjnodes.size(); j++) {
			dist.a[i][vertices[i].adjnodes[j]] = vertices[i].adjweight[j];
		}
	}

	// do floyd
	for(int k=0; k<nov; k++) {
		for(int i=0; i<nov; i++) {
			for(int j=0; j<nov; j++) {
				if(dist.a[i][k]==INT_MAX || dist.a[k][j]==INT_MAX) {
					continue;
				}
				if(dist.a[i][j] > dist.a[i][k] + dist.a[k][j]) {
					dist.a[i][j] = dist.a[i][k] + dist.a[k][j];
				}
			}
		}
	}

	return dist;


}


// cal rad for A* algorithm
double Graph::rad(double d) {
	double pi = 3.1415926535897;
	return d * pi / 180.0;
}


// heuristic function for A* algorithm
double Graph::heuristic_function(int s, int d) {
	double x1, y1, x2, y2;
	double dis;
	double earth_radius = 6378.137;

	x1 = rad(vertices[s].x);
	y1 = rad(vertices[s].y);
	x2 = rad(vertices[d].x);
	y2 = rad(vertices[d].y);

	dis = 2*asin(sqrt(pow( sin((x1-x2)/2), 2) + cos(x1)*cos(x2)*pow( sin(abs(y1-y2)/2), 2))) * earth_radius * 1000;

	return dis;
}


// A* algorithm
int Graph::a_star_algorithm(int s, int d){

	int result;
	int x, y;
	double tentative_g_score;
	std::pair<double, int > p;
	std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int> >, std::greater<std::pair<double, int> > > openset;
	std::set<int> closedset;
	std::map<int, int> came_from;
	std::map<int, int> g_score;

	came_from[s] = -1;
	g_score[s] = 0;


	openset.push(std::make_pair(0, s));

	while(!openset.empty()) {
		p = openset.top();
		openset.pop();	

		x = p.second;
		closedset.insert(x);  

		if(x==d) {
			return g_score[x];
		}

		for(int i=0; i<vertices[x].adjnodes.size(); i++) { 
			y = vertices[x].adjnodes[i];
			if(closedset.find(y) != closedset.end()) { 
				continue;
			}

			tentative_g_score =  g_score[x] + vertices[x].adjweight[i]; 
			
			if(g_score.find(y)==g_score.end() || g_score[y]>tentative_g_score) { 
				g_score[y] = tentative_g_score;
				openset.push(std::make_pair(tentative_g_score + heuristic_function(y, d), y));
				came_from[y] = x;
			}
		}
	}


}


bool Graph::edge_weight_update(int s, int t, int weight){

	bool updated = true;

	// if undirected graph
	if(type == 1) {
		if(std::find(vertices[s].adjnodes.begin(), vertices[s].adjnodes.end(), t) != vertices[s].adjnodes.end()
		&& std::find(vertices[t].adjnodes.begin(), vertices[t].adjnodes.end(), s) != vertices[t].adjnodes.end()) {
			int index;
			index = std::find(vertices[s].adjnodes.begin(), vertices[s].adjnodes.end(), t) - vertices[s].adjnodes.begin();
			vertices[s].adjweight[index] = weight;

			index = std::find(vertices[t].adjnodes.begin(), vertices[t].adjnodes.end(), s) - vertices[t].adjnodes.begin();
			vertices[t].adjweight[index] = weight;

		} else {
			updated = false;
			std::cout << "edge is not exist!" << std::endl;
		}

	// if undirected graph
	} else {
		if(std::find(vertices[s].adjnodes.begin(), vertices[s].adjnodes.end(), t) != vertices[s].adjnodes.end()) {
			int index = std::find(vertices[s].adjnodes.begin(), vertices[s].adjnodes.end(), t) - vertices[s].adjnodes.begin();
			vertices[s].adjweight[index] = weight;
		} else {
			updated = false;
			std::cout << "edge is not exist!" << std::endl;
		}

	}

	return updated;

}


// METIS setting options
void Graph::metis_options_setting(){
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // _RB
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // _VOL
	options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // _RM
	options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM; // _GROW _EDGE _NODE
	options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM; // _GREEDY _SEP2SIDED _SEP1SIDED
	// options[METIS_OPTION_NCUTS] = 1;
	// options[METIS_OPTION_NITER] = 10;
	/* balance factor, used to be 500 */
	options[METIS_OPTION_UFACTOR] = 500;
	// options[METIS_OPTION_MINCONN];
	options[METIS_OPTION_CONTIG] = 1;
	// options[METIS_OPTION_SEED];
	options[METIS_OPTION_NUMBERING] = 0;
	// options[METIS_OPTION_DBGLVL] = 0;
}



// graph partition
// input: nset = a set of node id
// output: <node, node belong to partition id>
std::unordered_map<int,int> Graph::graph_partition(int fanout, std::set<int> &nset){
	std::unordered_map<int,int> result;

	xadj = new idx_t[nset.size() + 1];
	adjncy = new idx_t[noe];
	adjwgt = new idx_t[noe];

	// transform data to metis
	// nvtxs, ncon
	nvtxs = nset.size();
	ncon = 1;

	int xadj_pos = 1;
	int xadj_accum = 0;
	int adjncy_pos = 0;

	// xadj, adjncy, adjwgt
	std::unordered_map<int,int> nodemap; // mapping: real node id -> new node id begin with zero
	nodemap.clear();

	xadj[0] = 0;
	int i = 0;
	for ( std::set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++ ){
		// init node map
		nodemap[*it] = i;

		int nid = *it;
		int f = vertices[nid].adjnodes.size();
		for ( int j = 0; j < f; j++ ){
			int enid = vertices[nid].adjnodes[j];
			if ( nset.find( enid ) != nset.end() ){
				xadj_accum ++;

				adjncy[adjncy_pos] = enid;
				adjwgt[adjncy_pos] = vertices[nid].adjweight[j];
				adjncy_pos ++;
			}
		}
		xadj[xadj_pos++] = xadj_accum;
	}

	// adjust nodes number started by 0
	for ( int i = 0; i < adjncy_pos; i++ ){
		adjncy[i] = nodemap[adjncy[i]];
	}

	// adjwgt -> 1
	if (adjweight_set_to_all_one){
		for ( int i = 0; i < adjncy_pos; i++ ){
			adjwgt[i] = 1;
		}
	}

	// nparts
	nparts = fanout;
	part = new idx_t[nset.size()];

	// partition, result -> part
	// k way partition
	METIS_PartGraphKway(
        &nvtxs,
        &ncon,
        xadj,
        adjncy,
        NULL,
        NULL,
        adjwgt,
        &nparts,
        NULL,
        NULL,
        options,
        &objval,
        part
    );	

	// push to result, key is real node id, value is part
	result.clear();
	i = 0;
	for ( std::set<int>::iterator it = nset.begin(); it != nset.end(); it++, i++ ){
		result[*it] = part[i];
	}

	// finalize
    delete xadj;
    delete adjncy;
    delete adjwgt;
    delete part;

	return result;

}