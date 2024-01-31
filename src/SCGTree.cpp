#include<iostream>
#include<vector>
#include<stack>
#include<unordered_map>
#include<set>
#include<algorithm>
#include<limits.h>
#include<cassert>
#include<sys/time.h>

#include"SCGTree.h"
#include"Graph.h"
#include"TreeNode.h"
#include"AuxiliaryData.h"

#define NDEBUG


SCGTree::SCGTree(Graph g): graph(g) {
	isborder_set.resize(graph.nov, false);
	treepath_set.resize(graph.nov, std::vector<int>());
}
 

void SCGTree::build(int fanout, int tau, int lamda){

	this->fanout = fanout;
	this->tau = tau;

	// std::cout << "build_tree_skeleton...";
	build_tree_skeleton(fanout, tau); 
	// std::cout << "done!" << std::endl;

	// std::cout << "calculate_auxiliary_data...";
	calculate_auxiliary_data();
	// std::cout << "done!" << std::endl;

	// std::cout << "calculate_distance_matrix_up_and_down...";
	calculate_distance_matrix_up_and_down();
	// std::cout << "done!" << std::endl;

	int matrix_size = 0;
    for(int i=0; i<tree.size(); i++){
        matrix_size += tree[i].dist.r * tree[i].dist.c;
    }
	int lamda_size = lamda * matrix_size / 100;

	// std::cout << "calculate_shortcut_matrix_with_lamda...";
	calculate_shortcut_matrix_with_lamda(lamda_size);
	// std::cout << "done!" << std::endl;

	// std::cout << "cal_short_cuts_dis...";
	cal_short_cuts_dis();
	// std::cout << "done!" << std::endl;

}


// graph partition and build tree skeleton
void SCGTree::build_tree_skeleton(int fanout, int tau){
    // init root
	TreeNode root;
	root.isleaf = false;
	root.father = -1;
	tree.push_back(root);

	// init stack
	std::stack<Status> buildstack;
	Status rootstatus;
	rootstatus.tnid = 0;
	rootstatus.nset.clear();
	for ( int i = 0; i < graph.nov; i++ ){
		rootstatus.nset.insert(i);
	}
	buildstack.push( rootstatus );

	// start to build
	std::unordered_map<int,int> presult;
	std::set<int> childset[fanout];

	while( buildstack.size() > 0 ){
		Status current = buildstack.top();
		buildstack.pop();

		// update treepath
		for ( std::set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
			treepath_set[*it].push_back( current.tnid );
		}

		// check limitation
		if ( current.nset.size() <= tau ){
			// build leaf node
			tree[current.tnid].isleaf = true;
			tree[current.tnid].leafnodes.clear();
			for ( std::set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
				tree[current.tnid].leafnodes.push_back( *it );
			}
			continue;
		}

		// partition
		presult = graph.graph_partition(fanout, current.nset );

		for ( int i = 0; i < fanout; i++ ){
			childset[i].clear();
		}
		int slot;
		for ( std::set<int>::iterator it = current.nset.begin(); it != current.nset.end(); it++ ){
			slot = presult[*it]; 
			childset[slot].insert(*it);
		}

		// create children nodes
		int childpos;
		for ( int i = 0; i < fanout; i++ ){
			TreeNode tnode;
			tnode.isleaf = false;
			tnode.father = current.tnid;
			
            tree.push_back(tnode);
			childpos = tree.size() - 1;

			tree[current.tnid].children.push_back( childpos );

			// identify border nodes
			tree[childpos].borders.clear();
			for ( std::set<int>::iterator it = childset[i].begin(); it != childset[i].end(); it++ ){

				bool isborder = false;
				for ( int j = 0; j < graph.vertices[*it].adjnodes.size(); j++ ){
					if ( childset[i].find( graph.vertices[*it].adjnodes[j] ) == childset[i].end() ){
						isborder = true;
						break;
					}
				}
				if ( isborder ){
					tree[childpos].borders.push_back(*it);
					isborder_set[*it] = true;
				}
			}

			Status ongoingstatus;
			ongoingstatus.tnid = childpos;
			ongoingstatus.nset = childset[i];
			buildstack.push(ongoingstatus);

		}
	}
}


// calculate the distance matrix, up_and_down method
void SCGTree::calculate_distance_matrix_up_and_down(){
	// initialize matrices
	for(int i=0; i<treesize; i++) {
		int vs;
		if(tree[i].isleaf) {
			vs = tree[i].leafnodes.size();
		} else {
			vs = tree[i].union_borders.size();
		}

		tree[i].dist.init(vs);
	
	}

	// add edges
	for(int i=0; i<graph.nov; i++) {
		int s, d, w, lca, lca_level, spos, dpos, tn;
		s = i;
		for(int j=0; j<graph.vertices[i].adjnodes.size(); j++) {
			d = graph.vertices[i].adjnodes[j];
			w = graph.vertices[i].adjweight[j];
			lca = find_lca(s, d);
			lca_level = find_lca_level(s, d);

			if(tree[lca].isleaf) {
				spos = tree[lca].vertex_in_leafnodes[s];
				dpos = tree[lca].vertex_in_leafnodes[d];
			} else {
				tn = treepath_set[s][lca_level+1];
				spos = tree[tn].border_in_father[tree[tn].border_in_borders_map[s]];
				
				tn = treepath_set[d][lca_level+1];
				dpos = tree[tn].border_in_father[tree[tn].border_in_borders_map[d]];

				assert(spos<tree[lca].dist.r && dpos<tree[lca].dist.r);
			}

			tree[lca].dist.a[spos][dpos] = w; 
			tree[lca].dist.midv[spos][dpos] = -1; 
		}
	}

	// up calculate
	for(int i=treesize-1; i>=0; i--){
		int tn, fn;
		tn = i;
		tree[tn].dist.do_floyd();

		if(tn!=0){
			// push up
			fn = tree[tn].father;
			int b1, b2, xpos, xpos1, ypos, ypos1;
			for(int j=0; j<tree[tn].borders.size(); j++) {
				b1 = tree[tn].borders[j];
				xpos = tree[tn].border_in_father[tree[tn].border_in_borders_map[b1]];// b1
				
				for(int k=0; k<tree[tn].borders.size(); k++) {
					b2 = tree[tn].borders[k];
					ypos = tree[tn].border_in_father[tree[tn].border_in_borders_map[b2]];// b2

					if(tree[tn].isleaf) { 
						xpos1 = tree[tn].vertex_in_leafnodes[b1]; // b1
						ypos1 = tree[tn].vertex_in_leafnodes[b2];    // b2

					} else {
						xpos1 = tree[tn].border_in_union_border[j]; // b1
						ypos1 = tree[tn].border_in_union_border[k]; // b2
					}

					if(tree[fn].dist.a[xpos][ypos] > tree[tn].dist.a[xpos1][ypos1]) { // push up
						tree[fn].dist.a[xpos][ypos] = tree[tn].dist.a[xpos1][ypos1]; 
						tree[fn].dist.midv[xpos][ypos] = -2;
					}	


				}
			}
		}
	}


	// down calculate
	for(int i=0; i<treesize; i++) {
		int tn=i, cn;
		if(!tree[i].isleaf) {
			for(int j=0; j<tree[i].children.size(); j++) {
				cn = tree[i].children[j];

				bool dirty = false;
				int b1, b2, xpos, xpos1, ypos, ypos1;
				for(int i=0; i<tree[cn].borders.size(); i++) {
					b1 = tree[cn].borders[i];
					xpos = tree[cn].border_in_father[tree[cn].border_in_borders_map[b1]]; // b1
					for(int j=0; j<tree[cn].borders.size(); j++) {
						b2 = tree[cn].borders[j];
						ypos = tree[cn].border_in_father[tree[cn].border_in_borders_map[b2]]; // b2

						if(tree[cn].isleaf) {
							xpos1 = tree[cn].vertex_in_leafnodes[b1]; // b1
							ypos1 = tree[cn].vertex_in_leafnodes[b2]; // b2
						} else {
							xpos1 = tree[cn].border_in_union_border[i]; // b1
							ypos1 = tree[cn].border_in_union_border[j]; // b2
						}

						if(tree[cn].dist.a[xpos1][ypos1] != tree[tn].dist.a[xpos][ypos]) {
							tree[cn].dist.a[xpos1][ypos1] = tree[tn].dist.a[xpos][ypos];
							tree[cn].dist.midv[xpos1][ypos1] = -3;
							dirty = true;
						}
					}
				}

				if(dirty) {
					// do floyd
					tree[cn].dist.do_floyd();
				}
			}
		}
	}

	// 重构leaf node matrix
	for(int i=0; i<treesize; i++) {
		if(tree[i].isleaf) {
			std::vector<int> minds, mid;
			for(int j=0; j<tree[i].borders.size(); j++) {
				int b = tree[i].borders[j];
				int bindex = tree[i].vertex_in_leafnodes[b];
				for(int k=0; k<tree[i].leafnodes.size(); k++) {
					minds.push_back(tree[i].dist.a[bindex][k]);
				}
			}

			for(int j=0; j<tree[i].leafnodes.size(); j++){
				for(int k=0; k<tree[i].leafnodes.size(); k++){
					mid.push_back(tree[i].dist.midv[j][k]);
				}
			}
			tree[i].dist.init_with_midv(tree[i].borders.size(), tree[i].leafnodes.size(), minds, mid);
		}
	}
}

 
// calculate auxiliary data
void SCGTree::calculate_auxiliary_data(){
	// tree nodes size
    treesize = tree.size();

	// union_borders
	std::vector<int> cands;
	std::set<int> nset;

	for(int tn=0; tn<treesize; tn++) {

		cands.clear();
		if ( tree[tn].isleaf ){
			cands = tree[tn].leafnodes;
			tree[tn].union_borders = tree[tn].borders;
		}
		else{
			nset.clear();
			for ( int k = 0; k < tree[tn].children.size(); k++ ){
				int cid = tree[tn].children[k];
				nset.insert( tree[cid].borders.begin(), tree[cid].borders.end() );
			}

			cands.clear();
			for ( std::set<int>::iterator it = nset.begin(); it != nset.end(); it ++ ){
				cands.push_back( *it );
			}
			tree[tn].union_borders = cands;
		}
	}

	std::vector<int> borders;
	bool is_border;
	border_max_level.clear();
	border_max_level_index.clear();
	for(int i=0; i< graph.nov; i++) {
		is_border = false;
		for(int j=1; j<treepath_set[i].size(); j++) {
			borders = tree[treepath_set[i][j]].borders;
			for(int k=0; k<borders.size(); k++) {
				if(borders[k] == i) {
					border_max_level.push_back(j); 
					border_max_level_index.push_back(k); 
					
					is_border = true;
					isborder_set[i] = is_border;
					break;
				}
			}
			if(is_border) {
				break;
			}
		}
		if(!is_border) {
			border_max_level.push_back(treepath_set[i].size());
			border_max_level_index.push_back(-1);
			isborder_set[i] = is_border;
		}
	}

	for(int i=1; i<treesize; i++) {
		std::vector<int> borders = tree[i].borders;
		for(int b : borders) {
			std::vector<int> father_union_borders = tree[tree[i].father].union_borders;
			for(int j=0; j<father_union_borders.size(); j++) {
				if(b==father_union_borders[j]) {
					tree[i].border_in_father.push_back(j);
				}
			}
		}
	}
		
	for(int i=1; i<treesize; i++) {
		std::vector<int> borders = tree[i].borders;
		std::vector<int> union_borders = tree[i].union_borders;
		for(int b : borders) {
			for(int j=0; j<union_borders.size(); j++) {
				if(b==union_borders[j]) {
					tree[i].border_in_union_border.push_back(j);
				}
			}
		}
	}

	for(int i=1; i<treesize; i++) {
		if(tree[i].isleaf) {
			std::vector<int> leafnodes = tree[i].leafnodes;
			for(int v : leafnodes) {
				for(int j=0; j<leafnodes.size(); j++) {
					if(v==leafnodes[j]) {
						tree[i].vertex_in_leafnodes[v] = j;
					}
				}
			}
		}
	}


	tree[0].level=0;
	for(int i=0; i<treesize; i++) {
		for(int j=0; j<tree[i].children.size(); j++) {
			tree[tree[i].children[j]].level = tree[i].level+1;
		}
	}

    for(int i=treesize-1; i>=0; i--) {
        if(tree[i].isleaf) {
            tree[i].vertexs_num = tree[i].leafnodes.size();
        } else {
            int vertexs_num=0;
            for(int j=0; j<tree[i].children.size(); j++) {
                vertexs_num += tree[tree[i].children[j]].vertexs_num;
            }
            tree[i].vertexs_num = vertexs_num;
        }
    }


	for(int i=1; i<treesize; i++) {
		std::vector<int> borders = tree[i].borders;
		for(int j=0; j<borders.size(); j++) {
			tree[i].border_in_borders_map[borders[j]] = j;
		}
	}


	tree[0].in_father_index = -1;
	for(int i=0; i<treesize; i++) {
		for(int j=0; j<tree[i].children.size(); j++) {
			int child = tree[i].children[j];
			tree[child].in_father_index = j;
		}
	}


	for(int i=0; i<treesize;i++) {
		if(tree[i].level >treedepth) {
			treedepth = tree[i].level;
		}
	}



}


void SCGTree::cal_short_cuts_dis(){

	int x, y;
	std::vector<std::vector<std::pair<int, int> > > short_cuts_dis(tree.size());
	for(int i=0; i<short_cut_pairs_in_level.size(); i++) {
		x = short_cut_pairs_in_level[i].first;
		y = short_cut_pairs_in_level[i].second;

		int min = short_cut_matrix_in_level[i].vmin;

		short_cuts_dis[x].push_back(std::make_pair(min, y));
		short_cuts_dis[y].push_back(std::make_pair(min, x));
	}

	bool first = true;
	for(int i=0; i<short_cuts_dis.size(); i++) {
		if(short_cuts_dis[i].empty()) continue;
		sort(short_cuts_dis[i].begin(), short_cuts_dis[i].end());

		for(int j=0; j<short_cuts_dis[i].size(); j++) {
			tree[i].short_cuts_dis.push_back(short_cuts_dis[i][j]);
		}
	}
}


// calculate all shortcuts at given level
void SCGTree::calculate_shortcut_matrix_with_given_level(int level){
	
	flevel=level;
	nflevel=level+1;
    std::vector<std::vector<int> > level_tnode;
	int total_border_num=0;

	level_tnode.push_back(std::vector<int>(1, 0));
	for(int l=0; l<level; l++) {
		std::vector<int> nds;
		for(int i=0; i<level_tnode[l].size(); i++) {
			for(int child : tree[level_tnode[l][i]].children) {
				nds.push_back(child);
			}
		}
		level_tnode.push_back(nds);
	}

	std::unordered_map<int, int> all_borders;
	int index = 0;
	for(int i=0; i<level_tnode[level].size(); i++) {
		int tn = level_tnode[level][i];
		for(int j=0; j<tree[tn].borders.size(); j++) {
			all_borders[tree[tn].borders[j]] = index;
			index ++;
		}
	}


	std::vector<std::vector<int> > dist(index, std::vector<int>(index, INT_MAX));
	int tn, x, y, ix, iy;
	for(int l=0; l<level_tnode.size()-1; l++) {
		for(int n=0; n<level_tnode[l].size(); n++) {
			tn = level_tnode[l][n];
			for(int u=0; u<tree[tn].union_borders.size(); u++) {
				for(int k=u+1; k<tree[tn].union_borders.size(); k++ ) {
					x = tree[tn].union_borders[u];
					y = tree[tn].union_borders[k];
					ix = all_borders[x];
					iy = all_borders[y];

					dist[ix][iy] = tree[tn].dist.a[u][k];
					dist[iy][ix] = tree[tn].dist.a[u][k];
				}
				
			}
		}
	}

	std::cout << "all borders size: " << all_borders.size() << std::endl;

	for(int i=0; i<index; i++) {
		dist[i][i] = 0;
	}

	for(int k=0; k<index; k++) {
		for(int i=0; i<index; i++) {
			for(int j=0; j<index; j++) {
				if(dist[i][k]==INT_MAX || dist[k][j]==INT_MAX) {
					continue;
				}
				if(dist[i][j] > dist[i][k] + dist[k][j]) {
					dist[i][j] = dist[i][k] + dist[k][j];
				}
			}
		}
	}

	short_cut_matrix_in_level.clear();
	short_cut_pairs_in_level.clear();
	for(int i=0; i<tree.size(); i++) {
		tree[i].short_cuts.clear();
	}

	std::vector<int> flevel = level_tnode[level];

	for(int i=0; i<flevel.size(); i++) {
        for(int j=i+1; j<flevel.size(); j++) {
            if(tree[flevel[i]].father == tree[flevel[j]].father) {
                continue;
            }
            std::pair<int, int> p = std::make_pair(flevel[i], flevel[j]);
            short_cut_pairs_in_level.push_back(p);
        }
    }
    short_cut_matrix_in_level.resize(short_cut_pairs_in_level.size());

	std::pair<int, int> p;
	std::vector<int> minds;
	int b1, b2, ib1, ib2, total_short_cut_size=0;
	for(int i=0; i<short_cut_pairs_in_level.size(); i++) {
		minds.clear();
		p = short_cut_pairs_in_level[i];

        tree[p.first].short_cuts[p.second] = i+1;
        tree[p.second].short_cuts[p.first] = -(i+1);

		for(int j=0; j<tree[p.first].borders.size(); j++) {
			b1 = tree[p.first].borders[j];
			assert(all_borders.find(b1) != all_borders.end());
			ib1 = all_borders[b1];
			for(int k=0; k<tree[p.second].borders.size(); k++) {
				b2 = tree[p.second].borders[k];
				assert(all_borders.find(b2) != all_borders.end());
				ib2 = all_borders[b2];
				minds.push_back(dist[ib1][ib2]);
			}
		}

		short_cut_matrix_in_level[i].init(tree[p.first].borders.size(), tree[p.second].borders.size(), minds);
		total_short_cut_size += minds.size();
	}

	printf("full filled shortcuts in level %d done!, total shortcuts size: %d \n", level, total_short_cut_size);

	std::cout << "calculate shortcut distance." << std::endl;

	std::vector<std::vector<std::pair<int, int> > > short_cuts_dis(tree.size());
	for(int i=0; i<short_cut_pairs_in_level.size(); i++) {
		x = short_cut_pairs_in_level[i].first;
		y = short_cut_pairs_in_level[i].second;

		int min = INT_MAX;
		for(int j=0; j<tree[x].borders.size(); j++) {
			for(int k=0; k<tree[y].borders.size(); k++) {
				if(min > short_cut_matrix_in_level[i].a[j][k]) {
					min = short_cut_matrix_in_level[i].a[j][k];
				}
			}
		}

		short_cuts_dis[x].push_back(std::make_pair(min, y));
		short_cuts_dis[y].push_back(std::make_pair(min, x));
	}

	bool first = true;
	for(int i=0; i<short_cuts_dis.size(); i++) {
		if(short_cuts_dis[i].empty()) continue;
		sort(short_cuts_dis[i].begin(), short_cuts_dis[i].end());

		for(int j=0; j<short_cuts_dis[i].size(); j++) {
			tree[i].short_cuts_dis.push_back(short_cuts_dis[i][j]);
		}
	}

    return;

}


// calculate all shortcuts at given level, save results
void SCGTree::calculate_shortcut_matrix_with_given_level(int level, std::vector< std::vector<int> > &all_minds, std::vector<std::pair<int, int> > &all_pairs) {

    std::vector<std::vector<int> > level_tnode;
	int total_border_num=0;

	level_tnode.push_back(std::vector<int>(1, 0));
	for(int l=0; l<level; l++) {
		std::vector<int> nds;
		for(int i=0; i<level_tnode[l].size(); i++) {
			for(int child : tree[level_tnode[l][i]].children) {
				nds.push_back(child);
			}
		}
		level_tnode.push_back(nds);

	}

	std::unordered_map<int, int> all_borders;
	int index = 0;
	for(int i=0; i<level_tnode[level].size(); i++) {
		int tn = level_tnode[level][i];
		for(int j=0; j<tree[tn].borders.size(); j++) {
			all_borders[tree[tn].borders[j]] = index;
			index ++;
		}
	}

	std::vector<std::vector<int> > dist(index, std::vector<int>(index, INT_MAX));

	for(int l=0; l<level_tnode.size()-1; l++) {
		for(int n=0; n<level_tnode[l].size(); n++) {
			int tn = level_tnode[l][n];
			for(int u=0; u<tree[tn].union_borders.size(); u++) {
				for(int k=u+1; k<tree[tn].union_borders.size(); k++ ) {
					int x = tree[tn].union_borders[u];
					int y = tree[tn].union_borders[k];
					int ix = all_borders[x];
					int iy = all_borders[y];

					dist[ix][iy] = tree[tn].dist.a[u][k];
					dist[iy][ix] = tree[tn].dist.a[u][k];
				}
				
			}
		}
	}

	for(int i=0; i<index; i++) {
		dist[i][i] = 0;
	}

	for(int k=0; k<index; k++) {
		for(int i=0; i<index; i++) {
			for(int j=0; j<index; j++) {
				if(dist[i][k]==INT_MAX || dist[k][j]==INT_MAX) {
					continue;
				}
				if(dist[i][j] > dist[i][k] + dist[k][j]) {
					dist[i][j] = dist[i][k] + dist[k][j];
				}
			}
		}
	}

	short_cut_pairs_in_level.clear();
	std::vector<int> flevel = level_tnode[level];

	for(int i=0; i<flevel.size(); i++) {
        for(int j=i+1; j<flevel.size(); j++) {
            if(tree[flevel[i]].father == tree[flevel[j]].father) {
                continue;
            }
            std::pair<int, int> p = std::make_pair(flevel[i], flevel[j]);
            short_cut_pairs_in_level.push_back(p);
        }
    }


	std::pair<int, int> p;
	std::vector<int> minds;
	int b1, b2, ib1, ib2;
	for(int i=0; i<short_cut_pairs_in_level.size(); i++) {
		minds.clear();
		p = short_cut_pairs_in_level[i];

		for(int j=0; j<tree[p.first].borders.size(); j++) {
			b1 = tree[p.first].borders[j];
			assert(all_borders.find(b1) != all_borders.end());
			ib1 = all_borders[b1];
			for(int k=0; k<tree[p.second].borders.size(); k++) {
				
				b2 = tree[p.second].borders[k];
				assert(all_borders.find(b2) != all_borders.end());
				ib2 = all_borders[b2];
				minds.push_back(dist[ib1][ib2]);
			}
		}

		assert(tree[p.first].borders.size()*tree[p.second].borders.size() == minds.size());
		all_minds.push_back(minds);
		all_pairs.push_back(std::make_pair(p.first, p.second));
	}

    return;
}


int SCGTree::cal_shortcut_value(int tnode1, int tnode2) {

    int level1, level2, bs1, bs2;
    int value1=0, value2=0, tot_value=0;

    if(tree[tnode1].level < tree[tnode2].level) std::swap(tnode1, tnode2);
    level1 = tree[tnode1].level;
    level2 = tree[tnode2].level;


    while(level1 < level2) {
        bs1 = tree[tnode1].borders.size();
        tnode1 = tree[tnode1].father;
        bs2 = tree[tnode1].borders.size();
        value1 += bs1 * bs2;
    }

    while(tree[tnode1].father != tree[tnode2].father) {
        bs1 = tree[tnode1].borders.size();
        tnode1 = tree[tnode1].father;
        bs2 = tree[tnode1].borders.size();
        value1 += bs1 * bs2;

        bs1 = tree[tnode2].borders.size();
        tnode2 = tree[tnode2].father;
        bs2 = tree[tnode2].borders.size();
        value2 += bs1 * bs2;

    }
    bs1 = tree[tnode1].borders.size();
    bs2 = tree[tnode2].borders.size();

    tot_value = value1 + value2 + bs1 * bs2;

    return tot_value;
}


void SCGTree::cal_temp_dist(int tnode, std::map<int, std::vector< std::vector<int> > > &temp_map) {
	std::vector< std::vector<int> > temp_vec;
	std::vector<int> minds, minds1;

	assert(tree[tnode].borders.size()>0);
	int ix, iy, fnode, ffnode;

	// 对tnode的每个border结点
	for(int j=0; j<tree[tnode].borders.size(); j++) {
		temp_vec.clear();
		minds.clear();
		assert(j< tree[tnode].border_in_father.size());

		ix = tree[tnode].border_in_father[j];
		fnode = tree[tnode].father;

		// tnode 的border 到父结点所有border的最短距离
		for(int k=0; k<tree[fnode].borders.size(); k++) {
			iy = tree[fnode].border_in_union_border[k];
			assert(ix<tree[fnode].dist.r && iy<tree[fnode].dist.c);
			minds.push_back(tree[fnode].dist.a[ix][iy]);
		}
		temp_vec.push_back(minds);


		while(tree[fnode].father != 0) {
			int ffnode = tree[fnode].father;
			minds1.clear();
			minds1.resize(tree[ffnode].borders.size(), INT_MAX);
			for(int k=0; k<minds1.size(); k++) {
				ix = tree[ffnode].border_in_union_border[k];
				for(int l=0; l<minds.size(); l++) {
					iy = tree[fnode].border_in_father[l];
					int dis = tree[ffnode].dist.a[ix][iy];
					if(minds1[k] > minds[l] + dis) {
						minds1[k] = minds[l] + dis;
					}
				}
			}
			assert(minds1.size()>0);
			temp_vec.push_back(minds1);
			minds = minds1;
			fnode = ffnode;
		}
		assert(temp_vec.size()>0);
		temp_map[tree[tnode].borders[j]] = temp_vec;
		assert(temp_map.size()>0);
	}


}


// calculate all shortcuts with space constraint
void SCGTree::calculate_shortcut_matrix_with_lamda(int lamda){

	for(int i=0; i<tree.size(); i++) {
		tree[i].short_cuts.clear();
	}
	short_cut_matrix_in_level.clear();
	short_cut_pairs_in_level.clear();

    std::vector<int> level_tnode, level_size_arr;
	int total_short_cut_size=0;

	level_size_arr.resize(treedepth, -1);
	int level_size=0;
	int slevel, flag=0;

	for(int l=2; l<treedepth; l++) {
		level_size = 0;

		level_tnode.clear();
		short_cut_pairs_in_level.clear();

		for(int i=0; i<tree.size(); i++) {
			if(tree[i].level==l ) {
				level_tnode.push_back(i);
			}
		}

		for(int i=0; i<level_tnode.size(); i++) {
			for(int j=i+1; j<level_tnode.size(); j++) {
				if(tree[level_tnode[i]].father == tree[level_tnode[j]].father) {
					continue;
				}
				level_size += tree[level_tnode[i]].borders.size() * tree[level_tnode[j]].borders.size();
			}
		}

		level_size_arr[l] = level_size;

		if(l==2) {
			if(lamda < level_size) {
				slevel = l;
				flag = 1;
				break;
			}
		} else {
			if(lamda < level_size) {
				slevel = l-1;
				break;
			}
		}
	}

	std::vector< std::vector<int> > all_minds; 
	std::vector<std::pair<int, int> > all_pairs; 
	if(flag==1) {
		// printf("non-full filled level: %d, \n", slevel);
		flevel = slevel-1;
		nflevel = slevel;

		int remain_lamda = lamda;

		// printf("full filled level done! \nfull filled level size: %d, non-full filled level size: %d\n", level_size_arr[slevel], remain_lamda);

		level_tnode.clear();
		short_cut_pairs_in_level.clear();
		std::vector<std::pair<int, int> > shortcut_value_in_level;

		for(int i=0; i<tree.size(); i++) {
			if(tree[i].level==slevel+1 ) {
				level_tnode.push_back(i);
			}
		}

		int k=0;
		for(int i=0; i<level_tnode.size(); i++) {
			for(int j=i+1; j<level_tnode.size(); j++) {
				if(tree[level_tnode[i]].father == tree[level_tnode[j]].father) {
					continue;
				}
				std::pair<int, int> p = std::make_pair(level_tnode[i], level_tnode[j]);
				short_cut_pairs_in_level.push_back(p);

				int vertex_set1=0, vertex_set2=0;
				vertex_set1 = tree[level_tnode[i]].vertexs_num;
				vertex_set2 = tree[level_tnode[j]].vertexs_num;

				int short_cut_size;
				short_cut_size = tree[level_tnode[i]].borders.size() * tree[level_tnode[j]].borders.size();
				int value = cal_shortcut_value(level_tnode[i], level_tnode[j]);

				int unit_value = (float(vertex_set1) /graph.nov) * (float(vertex_set2)/graph.nov) * value / short_cut_size;
				shortcut_value_in_level.push_back(std::make_pair(unit_value, k++));
			}
		}

		std::sort(shortcut_value_in_level.begin(), shortcut_value_in_level.end(), cmp);

		std::vector<std::pair<int, int> > candidate_shortcut_pairs_in_level;

		assert(shortcut_value_in_level.size() == short_cut_pairs_in_level.size());
		for(int i=0; i<shortcut_value_in_level.size(); i++) {
			std::pair<int, int> p1 = shortcut_value_in_level[i];
			std::pair<int, int> p = short_cut_pairs_in_level[p1.second];
			if(total_short_cut_size + tree[p.first].borders.size()*tree[p.second].borders.size() > remain_lamda) {
				break;
			}
			total_short_cut_size += tree[p.first].borders.size()*tree[p.second].borders.size();
			candidate_shortcut_pairs_in_level.push_back(p);
		}

		std::map<int, std::vector< std::vector<int> > > temp_map;
		std::set<int> temp_saved;
		int snode, dnode;
		std::vector<int> minds, sminds, dminds; 

		for(int i=0; i<candidate_shortcut_pairs_in_level.size(); i++) {
			minds.clear();
			std::pair<int, int> p = candidate_shortcut_pairs_in_level[i];
			snode = p.first;
			dnode = p.second;

			if(temp_saved.find(snode) == temp_saved.end()) {
				cal_temp_dist(snode, temp_map);
				temp_saved.insert(snode);

			}

			if(temp_saved.find(dnode) == temp_saved.end()) {
				cal_temp_dist(dnode, temp_map);
				temp_saved.insert(dnode);

			}

			int flevel=-1, ix, iy;

			while(tree[snode].father != tree[dnode].father) {
				snode = tree[snode].father;
				dnode = tree[dnode].father;
				flevel += 1;
			}
			int lca = tree[snode].father;

			for(int j=0; j<tree[p.first].borders.size(); j++) {
				int s = tree[p.first].borders[j];
				assert(temp_map.find(s) != temp_map.end() && flevel<temp_map[s].size());

				sminds = temp_map[s][flevel];
				for(int k=0; k<tree[p.second].borders.size(); k++) {
				
					int d = tree[p.second].borders[k];
					assert(temp_map.find(d) != temp_map.end() && flevel<temp_map[d].size());
					dminds = temp_map[d][flevel];

					int res = INT_MAX;
					for(int m=0; m<sminds.size(); m++) {
						ix = tree[snode].border_in_father[m];
						for(int n=0; n<dminds.size(); n++) {
							iy = tree[dnode].border_in_father[n];
							int dis = tree[lca].dist.a[ix][iy];
							if (res > sminds[m] + dis + dminds[n]) {
								res = sminds[m] + dis + dminds[n];
							}
						}
					}
					minds.push_back(res);
				}
			}

			all_minds.push_back(minds);
			all_pairs.push_back(std::make_pair(p.first, p.second) );
		}



	} else {

		// printf("full filled level: %d, non-full filld level: %d \n", slevel, slevel+1);
		flevel = slevel;
		nflevel = slevel+1;
		calculate_shortcut_matrix_with_given_level(slevel, all_minds, all_pairs);

		int remain_lamda = lamda - level_size_arr[slevel];

		// printf("full filled level done! \nfull filled level size: %d, non-full filled level size: %d\n", level_size_arr[slevel], remain_lamda);

		level_tnode.clear();
		short_cut_pairs_in_level.clear();
		std::vector<std::pair<int, int> > shortcut_value_in_level;

		for(int i=0; i<tree.size(); i++) {
			if(tree[i].level==slevel+1 ) {
				level_tnode.push_back(i);
			}
		}

		int k=0;
		for(int i=0; i<level_tnode.size(); i++) {
			for(int j=i+1; j<level_tnode.size(); j++) {
				if(tree[level_tnode[i]].father == tree[level_tnode[j]].father) {
					continue;
				}
				std::pair<int, int> p = std::make_pair(level_tnode[i], level_tnode[j]);
				short_cut_pairs_in_level.push_back(p);

				int vertex_set1=0, vertex_set2=0;
				vertex_set1 = tree[level_tnode[i]].vertexs_num;
				vertex_set2 = tree[level_tnode[j]].vertexs_num;

				int short_cut_size;
				short_cut_size = tree[level_tnode[i]].borders.size() * tree[level_tnode[j]].borders.size();
				int value = cal_shortcut_value(level_tnode[i], level_tnode[j]);

				int unit_value = (float(vertex_set1) /graph.nov) * (float(vertex_set2)/graph.nov) * value / short_cut_size;
				shortcut_value_in_level.push_back(std::make_pair(unit_value, k++));

			}
		}

		std::sort(shortcut_value_in_level.begin(), shortcut_value_in_level.end(), cmp);

		std::vector<std::pair<int, int> > candidate_shortcut_pairs_in_level;

		assert(shortcut_value_in_level.size() == short_cut_pairs_in_level.size());
		for(int i=0; i<shortcut_value_in_level.size(); i++) {
			std::pair<int, int> p1 = shortcut_value_in_level[i];
			std::pair<int, int> p = short_cut_pairs_in_level[p1.second];
			if(total_short_cut_size + tree[p.first].borders.size()*tree[p.second].borders.size() > remain_lamda) {
				break;
			}
			total_short_cut_size += tree[p.first].borders.size()*tree[p.second].borders.size();
			candidate_shortcut_pairs_in_level.push_back(p);
		}

		std::map<int, std::vector< std::vector<int> > > temp_map;
		std::set<int> temp_saved;
		int snode, dnode;
		std::vector<int> minds, sminds, dminds; 

		for(int i=0; i<candidate_shortcut_pairs_in_level.size(); i++) {
			minds.clear();
			std::pair<int, int> p = candidate_shortcut_pairs_in_level[i];
			snode = p.first;
			dnode = p.second;

			if(temp_saved.find(snode) == temp_saved.end()) {
				cal_temp_dist(snode, temp_map);
				temp_saved.insert(snode);
			}

			if(temp_saved.find(dnode) == temp_saved.end()) {
				cal_temp_dist(dnode, temp_map);
				temp_saved.insert(dnode);
			}

			int flevel=-1, ix, iy;

			while(tree[snode].father != tree[dnode].father) {
				snode = tree[snode].father;
				dnode = tree[dnode].father;
				flevel += 1;
			}
			int lca = tree[snode].father;

			for(int j=0; j<tree[p.first].borders.size(); j++) {
				int s = tree[p.first].borders[j];
				assert(temp_map.find(s) != temp_map.end() && flevel<temp_map[s].size());

				sminds = temp_map[s][flevel];
				for(int k=0; k<tree[p.second].borders.size(); k++) {
				
					int d = tree[p.second].borders[k];
					assert(temp_map.find(d) != temp_map.end() && flevel<temp_map[d].size());
					dminds = temp_map[d][flevel];

					int res = INT_MAX;
					for(int m=0; m<sminds.size(); m++) {
						ix = tree[snode].border_in_father[m];
						for(int n=0; n<dminds.size(); n++) {
							iy = tree[dnode].border_in_father[n];
							int dis = tree[lca].dist.a[ix][iy];
							if (res > sminds[m] + dis + dminds[n]) {
								res = sminds[m] + dis + dminds[n];
							}
						}
					}
					minds.push_back(res);
				}
			}
			all_minds.push_back(minds);
			all_pairs.push_back(std::make_pair(p.first, p.second) );
		}
	}

	assert(all_minds.size()==all_pairs.size());
	short_cut_matrix_in_level.resize(all_minds.size());
	short_cut_pairs_in_level.clear();

	for(int i=0; i<all_minds.size(); i++) {
		std::pair<int, int> p = all_pairs[i];
		tree[p.first].short_cuts[p.second] = i+1;
		tree[p.second].short_cuts[p.first] = -(i+1);

		assert(tree[p.first].borders.size() * tree[p.second].borders.size() == all_minds[i].size());
		short_cut_matrix_in_level[i].init(tree[p.first].borders.size(), tree[p.second].borders.size(), all_minds[i]);
		short_cut_pairs_in_level.push_back(p);
	}

	return;

}


// spsp query
int SCGTree::spsp_query(int s, int d){

	if(s==d) {
		return 0;
	}

	std::vector<int> spath = treepath_set[s];
	std::vector<int> dpath = treepath_set[d];
    bool use_short_cut = false;

    if(spath.back()==dpath.back()) {
		return graph.a_star_algorithm( s, d);
    }  else {

		int same_level=-1;
		if(spath.size()>dpath.size()){
			swap(spath, dpath);
		}
        // find lca
		for(int i=0; i<spath.size(); i++) {
			if(spath[i]==dpath[i]){
				same_level=i;
			}
		}
		
		if(tree[spath[flevel]].short_cuts.find(dpath[flevel]) != tree[spath[flevel]].short_cuts.end()) {
				same_level = flevel-1;
				use_short_cut = true;
		}

		if(nflevel != -1 && tree[spath[nflevel]].short_cuts.find(dpath[nflevel]) != tree[spath[nflevel]].short_cuts.end()) {
				same_level = nflevel-1;
				use_short_cut = true;
		}

		spath = treepath_set[s];
		dpath = treepath_set[d];

        std::vector<int> s_search_node_path; // s -> same level tree node path
		std::vector<int> d_search_node_path; // d -> same level tree node path

		// search_node_path
		s_search_node_path.push_back(-1);
        for(int l=border_max_level[s]-1; l>same_level; l--) {
			s_search_node_path.push_back(spath[l]);
        }

		d_search_node_path.push_back(-1);
        for(int l=border_max_level[d]-1; l>same_level; l--) {
			d_search_node_path.push_back(dpath[l]);
        }

		std::vector<int> s_dist, d_dist, dist;
		int x_border_size, y_border_size, ix, iy, weight, index;
		bool s_path_only_one_level, d_path_only_one_level;

		// s_dp_result
		s_dist.clear();
		s_dist.push_back(0);

		if(s_search_node_path.size()>1) {
			s_dist.clear();
			s_dist.resize(tree[s_search_node_path[1]].borders.size());
			if(tree[s_search_node_path[1]].isleaf) {
				index = tree[s_search_node_path[1]].vertex_in_leafnodes[s];
				for(int i=0; i<tree[s_search_node_path[1]].borders.size(); i++){
					s_dist[i] = tree[s_search_node_path[1]].dist.a[i][index];
				}
			} else {
				ix = tree[spath[border_max_level[s]]].border_in_father[border_max_level_index[s]];
				for(int i=0; i<tree[s_search_node_path[1]].borders.size(); i++){
					iy = tree[s_search_node_path[1]].border_in_union_border[i];
					s_dist[i] = tree[s_search_node_path[1]].dist.a[ix][iy];
				}
			}
		}

		if(s_search_node_path.size()>2) {
			for(int i=2; i<s_search_node_path.size(); i++) {
				x_border_size = tree[s_search_node_path[i-1]].borders.size();
				y_border_size = tree[s_search_node_path[i]].borders.size();
				dist.clear();
				dist.resize(y_border_size, INT_MAX);
				
				for(int j=0; j<y_border_size; j++){
					for(int k=0; k<x_border_size; k++) {
						ix = tree[s_search_node_path[i-1]].border_in_father[k];
						iy = tree[s_search_node_path[i]].border_in_union_border[j];
						weight = tree[s_search_node_path[i]].dist.a[ix][iy];
						if(dist[j] > s_dist[k] + weight){
							dist[j] = s_dist[k] + weight;
						}
					}
				}
				s_dist = dist;
			}
		}

		// d_dp_result
		d_dist.clear();
		d_dist.push_back(0);

		if(d_search_node_path.size()>1) {
			d_dist.clear();
			d_dist.resize(tree[d_search_node_path[1]].borders.size());
			if(tree[d_search_node_path[1]].isleaf) {
				index = tree[d_search_node_path[1]].vertex_in_leafnodes[d];
				for(int i=0; i<tree[d_search_node_path[1]].borders.size(); i++){
					d_dist[i] = tree[d_search_node_path[1]].dist.a[i][index];
				}
			} else { 
				ix = tree[dpath[border_max_level[d]]].border_in_father[border_max_level_index[d]];
				for(int i=0; i<tree[d_search_node_path[1]].borders.size(); i++){
					iy = tree[d_search_node_path[1]].border_in_union_border[i];
					d_dist[i] = tree[d_search_node_path[1]].dist.a[ix][iy];
				}
			}
		}
		
		if(d_search_node_path.size()>2) {
			for(int i=2; i<d_search_node_path.size(); i++) {
				x_border_size = tree[d_search_node_path[i-1]].borders.size();
				y_border_size = tree[d_search_node_path[i]].borders.size();
				dist.clear();
				dist.resize(y_border_size, INT_MAX);
				for(int j=0; j<y_border_size; j++){
					for(int k=0; k<x_border_size; k++) {
						ix = tree[d_search_node_path[i-1]].border_in_father[k];
						iy = tree[d_search_node_path[i]].border_in_union_border[j];
						weight = tree[d_search_node_path[i]].dist.a[ix][iy];
						if(dist[j] > d_dist[k] + weight){
							dist[j] = d_dist[k] + weight;
						}
					}
				}
				d_dist = dist;
			}
		}

		int min_dis = INT_MAX;
		x_border_size = s_dist.size();
		y_border_size = d_dist.size();

		s_path_only_one_level = s_search_node_path.back()== -1 ? true : false;
		d_path_only_one_level = d_search_node_path.back()== -1 ? true : false;

        if(use_short_cut) {
            int matrix_index = tree[spath[same_level+1]].short_cuts[dpath[same_level+1]];
            if(!s_path_only_one_level && !d_path_only_one_level) {
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        ix = i;
                        iy = j;
                        if(matrix_index>0) {
                            weight = short_cut_matrix_in_level[matrix_index-1].a[ix][iy];
                        } else {
                            weight = short_cut_matrix_in_level[-matrix_index-1].a[iy][ix];
                        }
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
                        }
                    }
                }
            } else if (s_path_only_one_level && !d_path_only_one_level) {
				ix = tree[spath[same_level+1]].border_in_borders_map[s];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        iy = j;
                        if(matrix_index>0) {
                            weight = short_cut_matrix_in_level[matrix_index-1].a[ix][iy];
                        } else {
                            weight = short_cut_matrix_in_level[-matrix_index-1].a[iy][ix];
                        }
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
                        }
                    }
                }
            } else if (!s_path_only_one_level && d_path_only_one_level) {
				iy = tree[dpath[same_level+1]].border_in_borders_map[d];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        ix = i;
                        if(matrix_index>0) {
                            weight = short_cut_matrix_in_level[matrix_index-1].a[ix][iy];
                        } else {
                            weight = short_cut_matrix_in_level[-matrix_index-1].a[iy][ix];
                        }
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
                        }
                    }
                }
            } else if (s_path_only_one_level && d_path_only_one_level) {
				ix = tree[spath[same_level+1]].border_in_borders_map[s];
				iy = tree[dpath[same_level+1]].border_in_borders_map[d];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        if(matrix_index>0) {
                            weight = short_cut_matrix_in_level[matrix_index-1].a[ix][iy];
                        } else {
                            weight = short_cut_matrix_in_level[-matrix_index-1].a[iy][ix];
                        }
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
                        }
                    }
                }
            }
        } else {
            if(!s_path_only_one_level && !d_path_only_one_level) {
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        ix = tree[s_search_node_path.back()].border_in_father[i];
                        iy = tree[d_search_node_path.back()].border_in_father[j];
                        weight = tree[dpath[same_level]].dist.a[ix][iy];
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
                        }
                    }
                }
            } else if (s_path_only_one_level && !d_path_only_one_level) {
				ix = tree[spath[same_level+1]].border_in_father[tree[spath[same_level+1]].border_in_borders_map[s]];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        iy = tree[d_search_node_path.back()].border_in_father[j];
                        weight = tree[dpath[same_level]].dist.a[ix][iy];
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
                        }
                    }
                }
            } else if (!s_path_only_one_level && d_path_only_one_level) {
				iy = tree[dpath[same_level+1]].border_in_father[tree[dpath[same_level+1]].border_in_borders_map[d]];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        ix = tree[s_search_node_path.back()].border_in_father[i];
                        weight = tree[dpath[same_level]].dist.a[ix][iy];
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
                        }
                    }
                }
            } else if (s_path_only_one_level && d_path_only_one_level) {
				ix = tree[spath[same_level+1]].border_in_father[tree[spath[same_level+1]].border_in_borders_map[s]];
				iy = tree[dpath[same_level+1]].border_in_father[tree[dpath[same_level+1]].border_in_borders_map[d]];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        weight = tree[dpath[same_level]].dist.a[ix][iy];
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
                        }
                    }
                }
            }
        }
		return min_dis;
    }

}


// identify candidate nodes
void SCGTree::cal_candidate_nodes(std::vector<int> &objects){
	for ( int i = 0; i < tree.size(); i++ ){
		tree[i].nonleafinvlist.clear();
		tree[i].leafinvlist.clear();
	}

	for ( int i = 0; i < objects.size(); i++ ){
		int oid = objects[i];
		assert(oid>=0 && oid<graph.nov);
		int max_level, current, pos, child;
		if(isborder_set[oid]){
			max_level = border_max_level[oid];
			if(max_level<flevel) {
				current = treepath_set[oid][flevel];
				pos = tree[current].border_in_borders_map[oid];
			} else {
				current = treepath_set[oid][max_level];
				pos = border_max_level_index[oid];
			}

			tree[current].leafinvlist.push_back(pos);
			while( current != -1 ){
				child = current;
				current = tree[current].father;
				if ( current == -1 ) break;
				if ( std::find(tree[current].nonleafinvlist.begin(), tree[current].nonleafinvlist.end(), child ) == tree[current].nonleafinvlist.end() ){
				tree[current].nonleafinvlist.push_back(child);
				}
			}

		} else {
			current = treepath_set[oid].back();

			pos = tree[current].vertex_in_leafnodes[oid];
			tree[current].leafinvlist.push_back(pos);		

			while( current != -1 ){
				child = current;
				current = tree[current].father;
				if ( current == -1 ) break;
				if ( std::find(tree[current].nonleafinvlist.begin(), tree[current].nonleafinvlist.end(), child ) == tree[current].nonleafinvlist.end() ){
					tree[current].nonleafinvlist.push_back(child);
				}
			}
		}
	}
}


// expansion for knn query
void SCGTree::expand_search_area_knn(int locid, int &expandable, int &tn, int &tmin, std::priority_queue<QueryStatus, 
					std::vector<QueryStatus >, query_status_competor > &pq, 
					std::unordered_map<int, std::vector<int> > &itm, std::vector<bool> &shortcuts_visited, 
					std::vector<int> &shortcuts_real_distance){

	int father, cid, bid, posa, posb, min, dis, allmin;
	father = tree[tn].father;
	int stn = tn;

	if(tree[father].level < flevel-1) { // full filled level
		tn = treepath_set[locid][flevel];
		int start_sc=-1, sc_tn, sc_index, x_y_dis;
		for(int i=0; i<shortcuts_visited.size(); i++) {
			if(!shortcuts_visited[i]) {
				sc_tn = tree[tn].short_cuts_dis[i].second;
				if(tree[sc_tn].nonleafinvlist.empty() && tree[sc_tn].leafinvlist.empty()){
					shortcuts_visited[i] = true;
					continue;
				}
				sc_index = tree[tn].short_cuts[sc_tn];

				if(start_sc==-1) {
					start_sc = i;
					if(shortcuts_real_distance[i]==-1) {
						allmin = INT_MAX;
						for(int i=0; i<tree[sc_tn].borders.size(); i++) {
							min = INT_MAX;
							for(int j=0; j<tree[tn].borders.size(); j++) {
								if( sc_index>0 ){
										x_y_dis = short_cut_matrix_in_level[sc_index-1].a[j][i];
									} else {
										x_y_dis = short_cut_matrix_in_level[-sc_index-1].a[i][j];
									}
								dis = x_y_dis + itm[tn][j];
								if(min > dis) {
									min = dis;
								}
							}
							itm[sc_tn].push_back(min);
							if(allmin > min) {
								allmin = min;
							}
						}
						shortcuts_real_distance[i] = allmin;

					} else {
						allmin = shortcuts_real_distance[i];
					}

					tmin = allmin;
					QueryStatus status = { sc_tn, false, allmin };
					pq.push(status);
					shortcuts_visited[i] = true;
					if(i==shortcuts_visited.size()-1) {
						expandable=-1;
						tmin = INT_MAX;
					}

					for(int i=0; i<tree[sc_tn].leafinvlist.size(); i++){
						posa = tree[sc_tn].leafinvlist[i];
						bid = tree[sc_tn].borders[posa];
						allmin = itm[sc_tn][posa];
						
						QueryStatus status = { bid, true, allmin };
						pq.push(status);
					}
					
				} else {
					if(tree[tn].short_cuts_dis[i].first <= tmin) {
						if(shortcuts_real_distance[i]==-1) {
							itm[sc_tn].clear();
							allmin = INT_MAX;
							for(int i=0; i<tree[sc_tn].borders.size(); i++) {
								min = INT_MAX;
								for(int j=0; j<tree[tn].borders.size(); j++) {
									if( sc_index>0 ){
											x_y_dis = short_cut_matrix_in_level[sc_index-1].a[j][i];
										} else {
											x_y_dis = short_cut_matrix_in_level[-sc_index-1].a[i][j];
										}
									dis = x_y_dis + itm[tn][j];
									if(min > dis) {
										min = dis;
									}
								}
								itm[sc_tn].push_back(min);
								if(allmin > min) {
									allmin = min;
								}
							}
							shortcuts_real_distance[i] = allmin;
						} else {
							allmin = shortcuts_real_distance[i];
						}

						if(allmin <= tmin) {
							QueryStatus status = { sc_tn, false, allmin };
							pq.push(status);
							shortcuts_visited[i] = true;
							if(i==shortcuts_visited.size()-1) {
								expandable=-1;
								tmin = INT_MAX;
							}

							for(int i=0; i<tree[sc_tn].leafinvlist.size(); i++){
								posa = tree[sc_tn].leafinvlist[i];
								bid = tree[sc_tn].borders[posa];
								allmin = itm[sc_tn][posa];
								QueryStatus status = { bid, true, allmin };
								pq.push(status);
							}
						}
					} else {
						break;
					}
				}	
			}
		}

		if(start_sc==-1) {
			expandable=-1;
			tmin = INT_MAX;
		}

		tn = stn;

	} else {
		itm[father].clear();
		tmin = INT_MAX;
		for(int i=0; i<tree[father].borders.size(); i++) {
			min = INT_MAX;
			posa = tree[father].border_in_union_border[i];
			for(int j=0; j<itm[tn].size(); j++) {
				posb = tree[tn].border_in_father[j];
				dis = tree[father].dist.a[posa][posb] + itm[tn][j];
				if(min > dis) {
					min = dis;
				}
			}
			itm[father].push_back(min);
			if(tmin > min) {
				tmin = min;
			}
		}
		for(int i=0; i<tree[father].nonleafinvlist.size(); i++) {
			cid = tree[father].nonleafinvlist[i];
			if( cid == tn) continue;
			itm[cid].clear();
			allmin = INT_MAX;
			for ( int j = 0; j<tree[cid].borders.size(); j++ ){
				min = INT_MAX;
				posa = tree[cid].border_in_father[j];
				for(int k=0; k<tree[tn].borders.size(); k++) {
					posb = tree[tn].border_in_father[k];
					dis = tree[father].dist.a[posa][posb] + itm[tn][k];
					if(min > dis ) {
						min = dis;
					}
				}	
				itm[cid].push_back(min);
				if(allmin > min) {
					allmin = min;
				}
			}
			QueryStatus status = { cid, false, allmin };
			pq.push(status);
		}

		for(int i=0; i<tree[father].leafinvlist.size(); i++) {
			posa = tree[father].leafinvlist[i];
			bid = tree[father].borders[posa];
			allmin = itm[father][posa];
			
			QueryStatus status = { bid, true, allmin };
			pq.push(status);
		}
		tn = father;
	}

}


// kNN query
std::vector<ResultSet> SCGTree::knn_query( int locid, int K){
	std::priority_queue<QueryStatus, std::vector<QueryStatus >, query_status_competor > pq;
	std::vector<ResultSet> rstset;
	rstset.clear();

	std::unordered_map<int, std::vector<int> > itm;
	
	itm.clear();
	int tn, tmin, cid, bid, posa, posb, min, dis, allmin;

	tn = treepath_set[locid].back();
	itm[tn].clear();

	tmin = INT_MAX;
	posa = tree[tn].vertex_in_leafnodes[locid];
	for ( int j = 0; j < tree[tn].borders.size(); j++ ){
		dis = tree[tn].dist.a[j][posa];
		itm[tn].push_back(dis);
		if(tmin > dis) {
			tmin = dis;
		}
	}

	std::vector<int> cands, result;
	cands.clear();
	for ( int i = 0; i < tree[tn].leafinvlist.size(); i++ ){
		cands.push_back( tree[tn].leafnodes[tree[tn].leafinvlist[i]] );
	}
	result = graph.dijkstra( locid, cands );

	for ( int i = 0; i < cands.size(); i++ ){
		QueryStatus status = { cands[i], true, result[i] };
		pq.push(status);
	}

	std::vector<bool> shortcuts_visited(tree[treepath_set[locid][flevel]].short_cuts.size(), false);
	std::vector<int> shortcuts_real_distance(tree[treepath_set[locid][flevel]].short_cuts.size(), -1);
	int expandable=0;
	// do search
	while(rstset.size()<K && (!pq.empty() || expandable!=-1)) {
		while(pq.empty()) {
			expand_search_area_knn(locid, expandable, tn, tmin, pq, itm, shortcuts_visited, shortcuts_real_distance);
		}
		QueryStatus top = pq.top();
		pq.pop();
		if(top.dis>tmin) {
			expand_search_area_knn(locid, expandable, tn, tmin, pq, itm, shortcuts_visited, shortcuts_real_distance);
			pq.push(top);
		} else if(top.isvertex) {
			ResultSet rs = { top.id, top.dis };
			rstset.push_back(rs);
		} else if(!top.isvertex) {
			if(tree[top.id].isleaf) {
				for(int i=0; i<tree[top.id].leafinvlist.size(); i++){
					min = INT_MAX;
					posa = tree[top.id].leafinvlist[i];
					for(int j=0; j<tree[top.id].borders.size(); j++) {
						posb = j;
						dis = tree[top.id].dist.a[posb][posa] + itm[top.id][j];
						if(min > dis) {
							min = dis;
						}
					}

					QueryStatus status = { tree[top.id].leafnodes[tree[top.id].leafinvlist[i]], true, min };
					pq.push(status);
				}
			} else {
				for(int i=0; i<tree[top.id].nonleafinvlist.size(); i++){
					cid = tree[top.id].nonleafinvlist[i];
					allmin = INT_MAX;
					for(int j=0; j<tree[cid].borders.size(); j++) {
						posa = tree[cid].border_in_father[j];
						min = INT_MAX;
						for(int k=0; k<tree[top.id].borders.size(); k++) {
							posb = tree[top.id].border_in_union_border[k];
							dis = tree[top.id].dist.a[posa][posb] + itm[top.id][k];
							if(min > dis) {
								min = dis;
							}
						}
						itm[cid].push_back(min);
						if(allmin > min) {
							allmin = min;
						}
					}
					QueryStatus status = { cid, false, allmin };
					pq.push(status);
				}
				for(int i=0; i<tree[top.id].leafinvlist.size(); i++){					
					posa = tree[top.id].leafinvlist[i];
					bid = tree[top.id].borders[posa];
					allmin = itm[top.id][posa];
					QueryStatus status = { bid, true, allmin };
					pq.push(status);
				}
			}
		}
	}

	return rstset;

}


// expansion for range query
void SCGTree::expand_search_area_range(int locid, int range, int &expandable, int &tn, int &tmin, std::queue<QueryStatus> &q, 
								std::unordered_map<int, std::vector<int> > &itm, std::vector<bool> &shortcuts_visited, 
								std::vector<int> &shortcuts_real_distance, std::vector<ResultSet> rstset){


	int father, cid, bid, posa, posb, min, dis, allmin;
	father = tree[tn].father;
	int stn = tn;

	if(tree[father].level < flevel-1) { // full filled level
		tn = treepath_set[locid][flevel];
		int start_sc=-1, sc_tn, sc_index, x_y_dis;
		for(int i=0; i<shortcuts_visited.size(); i++) {
			if(!shortcuts_visited[i]) {
				sc_tn = tree[tn].short_cuts_dis[i].second;
				if(tree[sc_tn].nonleafinvlist.empty() && tree[sc_tn].leafinvlist.empty()){
					shortcuts_visited[i] = true;
					continue;
				}
				sc_index = tree[tn].short_cuts[sc_tn];

				if(start_sc==-1) {
					start_sc = i;
					if(shortcuts_real_distance[i] == -1) {
						itm[sc_tn].clear();
						allmin = INT_MAX;
						for(int i=0; i<tree[sc_tn].borders.size(); i++) {
							min = INT_MAX;
							for(int j=0; j<tree[tn].borders.size(); j++) {
								if( sc_index>0 ){
										x_y_dis = short_cut_matrix_in_level[sc_index-1].a[j][i];
									} else {
										x_y_dis = short_cut_matrix_in_level[-sc_index-1].a[i][j];
									}
								dis = x_y_dis + itm[tn][j];
								if(min > dis) {
									min = dis;
								}
							}
							itm[sc_tn].push_back(min);
							if(allmin > min) {
								allmin = min;
							}
						}
						shortcuts_real_distance[i] = allmin;
					} else {
						allmin = shortcuts_real_distance[i];
					}
					tmin = allmin;

					shortcuts_visited[i] = true;
					if(i==shortcuts_visited.size()-1) {
						expandable=-1;
					}

					if(allmin>range) {
						start_sc = -1;
						continue;
					}
					QueryStatus status = { sc_tn, true, allmin };
					q.push(status);

					for(int i=0; i<tree[sc_tn].leafinvlist.size(); i++){
						posa = tree[sc_tn].leafinvlist[i];
						bid = tree[sc_tn].borders[posa];
						min = itm[sc_tn][posa];
						if(min<=range){
							ResultSet rs = { bid, min };
							rstset.push_back(rs);
						}
					}
					
				} else {
					if(tree[tn].short_cuts_dis[i].first <= tmin) {
						if(shortcuts_real_distance[i] == -1) {
							itm[sc_tn].clear();
							allmin = INT_MAX;
							for(int i=0; i<tree[sc_tn].borders.size(); i++) {
								min = INT_MAX;
								for(int j=0; j<tree[tn].borders.size(); j++) {
									if( sc_index>0 ){
											x_y_dis = short_cut_matrix_in_level[sc_index-1].a[j][i];
										} else {
											x_y_dis = short_cut_matrix_in_level[-sc_index-1].a[i][j];
										}
									dis = x_y_dis + itm[tn][j];
									if(min > dis) {
										min = dis;
									}
								}
								itm[sc_tn].push_back(min);
								if(allmin > min) {
									allmin = min;
								}
							}
							shortcuts_real_distance[i] = allmin;

						} else {
							allmin = shortcuts_real_distance[i];
						}

						shortcuts_visited[i] = true;

						if(i==shortcuts_visited.size()-1) {
							expandable=-1;
						}

						if(allmin <= range) {
							QueryStatus status = { sc_tn, true, allmin };
							q.push(status);
							for(int i=0; i<tree[sc_tn].leafinvlist.size(); i++){
								posa = tree[sc_tn].leafinvlist[i];
								bid = tree[sc_tn].borders[posa];
								min = itm[sc_tn][posa];
								if(min<=range){
									ResultSet rs = { bid, min };
									rstset.push_back(rs);
								}
							}
						}
					} else {
						break;
					}
				}	
			}
		}
		if(start_sc==-1) {
			expandable=-1;
		}
		tn = stn;

	} else {
		itm[father].clear();
		tmin = INT_MAX;
		for(int i=0; i<tree[father].borders.size(); i++) {
			min = INT_MAX;
			posa = tree[father].border_in_union_border[i];
			for(int j=0; j<itm[tn].size(); j++) {
				posb = tree[tn].border_in_father[j];
				dis = tree[father].dist.a[posa][posb] + itm[tn][j];
				if(min > dis) {
					min = dis;
				}
			}
			itm[father].push_back(min);
			if(tmin > min) {
				tmin = min;
			}
		}

		for(int i=0; i<tree[father].nonleafinvlist.size(); i++) {
			cid = tree[father].nonleafinvlist[i];
			if( cid == tn) continue;
			itm[cid].clear();
			allmin = INT_MAX;
			for ( int j = 0; j<tree[cid].borders.size(); j++ ){
				min = INT_MAX;
				posa = tree[cid].border_in_father[j];
				for(int k=0; k<tree[tn].borders.size(); k++) {
					posb = tree[tn].border_in_father[k];
					dis = tree[father].dist.a[posa][posb] + itm[tn][k];
					if(min > dis ) {
						min = dis;
					}
				}	
				itm[cid].push_back(min);
				if(allmin > min) {
					allmin = min;
				}
			}
			if(allmin <= range){
				QueryStatus status = { cid, true, allmin };
				q.push(status);
			}
		}

		for(int i=0; i<tree[father].leafinvlist.size(); i++) {
			posa = tree[father].leafinvlist[i];
			bid = tree[father].borders[posa];
			min = itm[father][posa];

			if(min<=range){
				ResultSet rs = { bid, min };
				rstset.push_back(rs);
			}
		}
		tn = father;
	}
}


// range query
std::vector<ResultSet> SCGTree::range_query( int locid, int range){

	std::queue<QueryStatus> q;
	std::vector<ResultSet> rstset;
	rstset.clear();
	std::unordered_map<int, std::vector<int> > itm;
	
	itm.clear();
	int tn, tmin, cid, bid, posa, posb, min, dis, allmin;

	tn = treepath_set[locid].back();
	itm[tn].clear();

	tmin = INT_MAX;
	posa = tree[tn].vertex_in_leafnodes[locid];
	for ( int j = 0; j < tree[tn].borders.size(); j++ ){
		dis = tree[tn].dist.a[j][posa];

		itm[tn].push_back(dis);
		if(tmin > dis) {
			tmin = dis;
		}
	}

	std::vector<int> cands, result;

	cands.clear();
	for ( int i = 0; i < tree[tn].leafinvlist.size(); i++ ){
		cands.push_back( tree[tn].leafnodes[tree[tn].leafinvlist[i]] );
	}
	result = graph.dijkstra( locid, cands );
	for ( int i = 0; i < cands.size(); i++ ){
		if(result[i]<=range){
			ResultSet rs = { cands[i], result[i] };
			rstset.push_back(rs);
		}
	}

	std::vector<bool> shortcuts_visited(tree[treepath_set[locid][flevel]].short_cuts.size(), false);
	std::vector<int> shortcuts_real_distance(tree[treepath_set[locid][flevel]].short_cuts.size(), -1);
	int expandable=0;

	while(!q.empty() || (expandable!=-1 && tmin<=range)) {
		if(q.empty()){
			while(q.empty() && (expandable!=-1 && tmin<=range)) {
				expand_search_area_range(locid, range, expandable, tn, tmin, q, itm, 
										shortcuts_visited, shortcuts_real_distance, rstset);
			}
			if(q.empty()) break;
		}
		assert(!q.empty());

		QueryStatus top = q.front();
		q.pop();

		if(tree[top.id].isleaf) {
			for(int i=0; i<tree[top.id].leafinvlist.size(); i++){
				min = INT_MAX;
				posa = tree[top.id].leafinvlist[i];
				for(int j=0; j<tree[top.id].borders.size(); j++) {
					posb = j;
					dis = tree[top.id].dist.a[posb][posa] + itm[top.id][j];
					if(min > dis) {
						min = dis;
					}
				}
				if(min<=range){
					ResultSet rs = { tree[top.id].leafnodes[tree[top.id].leafinvlist[i]], min };
					rstset.push_back(rs);
				}
			}
		} else {
			for(int i=0; i<tree[top.id].nonleafinvlist.size(); i++){
				cid = tree[top.id].nonleafinvlist[i];
				allmin = INT_MAX;
				for(int j=0; j<tree[cid].borders.size(); j++) {
					posa = tree[cid].border_in_father[j];
					min = INT_MAX;
					for(int k=0; k<tree[top.id].borders.size(); k++) {
						posb = tree[top.id].border_in_union_border[k];
						dis = tree[top.id].dist.a[posa][posb] + itm[top.id][k];
						if(min > dis) {
							min = dis;
						}
					}
					itm[cid].push_back(min);
					if(allmin > min) {
						allmin = min;
					}
				}

				if(allmin<=range){
					QueryStatus status = { cid, true, allmin };
					q.push(status);
				}
			}
			for(int i=0; i<tree[top.id].leafinvlist.size(); i++){					
				posa = tree[top.id].leafinvlist[i];
				bid = tree[top.id].borders[posa];
				min = itm[top.id][posa];
				if(min<=range){
					ResultSet rs = { bid, min };
					rstset.push_back(rs);
				}
			}
		}

	}

	return rstset;
}


// dump tree
void SCGTree::save(const char *scgtree_file, bool compress=false, bool dump_path=false){

	FILE *fout = fopen( scgtree_file, "wb" );
	int *buf = new int[ graph.nov ];

	int tree_size = treesize;
	fwrite( &tree_size, sizeof(int), 1, fout );
	int nodes_size = graph.nov;
	fwrite( &nodes_size, sizeof(int), 1, fout );

	for ( int i = 0; i < treesize; i++ ){
		// borders
		int count_borders = tree[i].borders.size();
		fwrite( &count_borders, sizeof(int), 1, fout );
		std::copy( tree[i].borders.begin(), tree[i].borders.end(), buf );
		fwrite( buf, sizeof(int), count_borders, fout );
		// children
		int count_children = tree[i].children.size();
		fwrite( &count_children, sizeof(int), 1, fout );
		std::copy( tree[i].children.begin(), tree[i].children.end(), buf );
		fwrite( buf, sizeof(int), count_children, fout );
		// isleaf
		fwrite( &tree[i].isleaf, sizeof(bool), 1, fout );
		// leafnodes
		int count_leafnodes = tree[i].leafnodes.size();
		fwrite( &count_leafnodes, sizeof(int), 1, fout );
		std::copy( tree[i].leafnodes.begin(), tree[i].leafnodes.end(), buf );
		fwrite( buf, sizeof(int), count_leafnodes, fout );
		// father
		fwrite( &tree[i].father, sizeof(int), 1, fout );
	}

	// path
	if(dump_path) {
		for ( int i = 0; i < graph.nov; i++ ){
			int count = treepath_set[i].size();
			fwrite( &count, sizeof(int), 1, fout );
			std::copy( treepath_set[i].begin(), treepath_set[i].end(), buf );
			fwrite( buf, sizeof(int), count, fout );
		}
	}



	int count;

	for ( int i = 0; i < tree.size(); i++ ){
		// union borders
		count = tree[i].union_borders.size();
		fwrite( &count, sizeof(int), 1, fout );
		buf = new int[count];
		std::copy( tree[i].union_borders.begin(), tree[i].union_borders.end(), buf );
		fwrite( buf, sizeof(int), count, fout );
		delete[] buf;

		if(compress) {
			if(tree[i].isleaf) {
				// mind
				count = tree[i].mind.size();
				fwrite( &count, sizeof(int), 1, fout );
				buf = new int[count];
				std::copy( tree[i].mind.begin(), tree[i].mind.end(), buf );
				fwrite( buf, sizeof(int), count, fout );
				delete[] buf;
			} else {
				// mind
				std::vector<int> compress_minds;
				for(int j=0; j<count; j++) {
					for(int k=0; k<count; k++) {
						if(j<k) {
							compress_minds.push_back(tree[i].mind[j*count+k]);
						}
					}
				}
				count = compress_minds.size();

				fwrite( &count, sizeof(int), 1, fout );
				buf = new int[count];
				copy( compress_minds.begin(), compress_minds.end(), buf );
				fwrite( buf, sizeof(int), count, fout );
				delete[] buf;
			}
		} else {
			count = tree[i].dist.r * tree[i].dist.c;
			fwrite( &count, sizeof(int), 1, fout );
			for(int j=0; j<tree[i].dist.r; j++) {
				fwrite( tree[i].dist.a[j], sizeof(int), tree[i].dist.c, fout );
			}
		}
	}

	int short_cut_size = short_cut_pairs_in_level.size();
	fwrite( &short_cut_size, sizeof(int), 1, fout );

	for ( int i = 0; i < short_cut_pairs_in_level.size(); i++ ){

		int x = short_cut_pairs_in_level[i].first;
		int y = short_cut_pairs_in_level[i].second;
		
		fwrite( &x, sizeof(int), 1, fout );
		fwrite( &y, sizeof(int), 1, fout );

		for(int j=0; j<tree[x].borders.size(); j++) {
			fwrite( short_cut_matrix_in_level[i].a[j], sizeof(int), tree[y].borders.size(), fout );
		}
	}
	fwrite( &flevel, sizeof(int), 1, fout );
	fwrite( &nflevel, sizeof(int), 1, fout );

	fclose(fout);

}


// load tree
void SCGTree::load(const char *scgtree_file, bool compress=false, bool dump_path=false){

	FILE *fin = fopen( scgtree_file, "rb" );
	int *buf = new int[ graph.nov ];
	int count_borders, count_children, count_leafnodes;
	bool isleaf;
	int father;

	tree.clear();

	int tree_size;
	int nodes_size;
	fread( &tree_size, sizeof(int), 1, fin );
	fread( &nodes_size, sizeof(int), 1, fin );

	tree.resize(tree_size);

	for(int i=0; i<tree_size; i++) {
		fread( &count_borders, sizeof(int), 1, fin );
		tree[i].borders.clear();
		fread( buf, sizeof(int), count_borders, fin );
		for ( int j = 0; j < count_borders; j++ ){
			tree[i].borders.push_back(buf[j]);
		}
		// children
		fread( &count_children, sizeof(int), 1, fin );
		fread( buf, sizeof(int), count_children, fin );

		tree[i].children.clear();
		for ( int j = 0; j < count_children; j++ ){
			tree[i].children.push_back(buf[j]);
		}
		// isleaf
		fread( &isleaf, sizeof(bool), 1, fin );
		tree[i].isleaf = isleaf;
		// leafnodes
		fread( &count_leafnodes, sizeof(int), 1, fin );
		fread( buf, sizeof(int), count_leafnodes, fin );
		tree[i].leafnodes.clear();
		for ( int j = 0; j < count_leafnodes; j++ ){
			tree[i].leafnodes.push_back(buf[j]);
		}
		// father
		fread( &father, sizeof(int), 1, fin );
		tree[i].father = father;
	}

	int count;
	if(dump_path) {
		for(int pos=0; pos<nodes_size; pos++) {

			fread( &count, sizeof(int), 1, fin );
			fread( buf, sizeof(int), count, fin );

			treepath_set[pos].clear();
			for ( int i = 0; i < count; i++ ){
				treepath_set[pos].push_back( buf[i] );
			}
		}
	} else {
		std::vector<int> leaf;
		for(int i=0; i<tree.size(); i++) {
			if(tree[i].isleaf) {
				leaf.push_back(i);
			}
		}

		int v, tnode;
		for(int i=0; i<leaf.size(); i++) {
			tnode = leaf[i];
			while(tnode != -1) {
				for(int j=0; j<tree[leaf[i]].leafnodes.size(); j++) {
					v = tree[leaf[i]].leafnodes[j];
					treepath_set[v].push_back(tnode);
				}
				tnode = tree[tnode].father;
			}

			for(int j=0; j<tree[leaf[i]].leafnodes.size(); j++) {
				v = tree[leaf[i]].leafnodes[j];
				std::reverse(treepath_set[v].begin(), treepath_set[v].end());
			}
		}
	}

	for(int pos=0; pos<tree_size; pos++) {
		fread( &count, sizeof(int), 1, fin );

		// union borders
		buf = new int[count];
		fread( buf, sizeof(int), count, fin );
		tree[pos].union_borders.clear();
		for ( int i = 0; i < count; i++ ){
			tree[pos].union_borders.push_back(buf[i]);
		}
		delete[] buf;

		if(compress) {
			if(tree[pos].isleaf) {
				// mind
				fread( &count, sizeof(int), 1, fin );
				buf = new int[count];
				fread( buf, sizeof(int), count, fin );
				tree[pos].mind.clear();
				for ( int i = 0; i < count; i++ ){
					tree[pos].mind.push_back(buf[i]);
				}

				delete[] buf;
			} else {

				// mind
				fread( &count, sizeof(int), 1, fin );
				buf = new int[count];
				fread( buf, sizeof(int), count, fin );
				tree[pos].mind.clear();

				int ubs = tree[pos].union_borders.size();
				std::vector<std::vector<int> > arr(ubs, std::vector<int>(ubs, 0));


				int index = 0;
				for(int i=0; i<ubs; i++) {
					for(int j=0; j<ubs; j++) { 
						if(i<j) {
							arr[i][j] = buf[index];
							arr[j][i] = buf[index];
							index ++;
						}
					}
				}

				for(int i=0; i<ubs; i++) {
					for(int j=0; j<ubs; j++) { 
						tree[pos].mind.push_back(arr[i][j]);
					}
				}

				delete[] buf;
			}
		} else {
			// mind
			fread( &count, sizeof(int), 1, fin );
			buf = new int[count];
			fread( buf, sizeof(int), count, fin );
			tree[pos].mind.clear();
			for ( int i = 0; i < count; i++ ){
				tree[pos].mind.push_back(buf[i]);
			}
			delete[] buf;
		}
	}

	int short_cut_size, x, y;

	fread( &short_cut_size, sizeof(int), 1, fin );
	short_cut_matrix_in_level.resize(short_cut_size);

	for(int i=0; i<short_cut_size; i++  ){
		fread( &x, sizeof(int), 1, fin );
		fread( &y, sizeof(int), 1, fin );
		short_cut_pairs_in_level.push_back(std::make_pair(x, y));
		// borders
		std::vector<int> minds;
		
		for(int j=0; j<tree[x].borders.size(); j++) {
			buf = new int[ tree[y].borders.size() ];
			fread( buf, sizeof(int), tree[y].borders.size(), fin );
			for(int k=0;k<tree[y].borders.size(); k++) {
				minds.push_back(buf[k]);
			}
			delete [] buf;
		}
		short_cut_matrix_in_level[i].init(tree[x].borders.size(), tree[y].borders.size(), minds);

		tree[x].short_cuts[y] = i+1;
        tree[y].short_cuts[x] = -(i+1);
	}

	fread( &flevel, sizeof(int), 1, fin );
	fread( &nflevel, sizeof(int), 1, fin );

	fclose(fin);

	// Distance Matrix
	for(int i=0; i<tree.size(); i++) {
		if(tree[i].isleaf) {
			tree[i].dist.init( tree[i].borders.size(),  tree[i].leafnodes.size(), tree[i].mind);
		} else {
			tree[i].dist.init(tree[i].union_borders.size(), tree[i].union_borders.size(), tree[i].mind);
		}
	}
	calculate_auxiliary_data();
	cal_short_cuts_dis();
}


// edge update
bool SCGTree::edge_weight_update(int s, int d, int weight, bool show_info, bool with_shortcut){

    struct timeval tv;
    long long t1, t2, update_time, sc_update_time;

    gettimeofday( &tv, NULL );
    t1 = tv.tv_sec * 1000000 + tv.tv_usec ;

	std::vector<int> ulist, dlist;
	std::vector<bool> updated_tn = std::vector<bool>(treesize, false);
	std::vector<bool> add_to_uq=updated_tn, add_to_dq=updated_tn;

	std::set<int> updated_shortcuts; // updated shortcuts

	std::priority_queue<int, std::vector<int >, std::less<int> > uq; // push up
	std::priority_queue<int, std::vector<int >, std::greater<int> > dq; // push down

	// update graph and re-initialize matrices
	bool graph_updated = graph.edge_weight_update(s, d, weight);

	if(!graph_updated) {
		return false;
	}

	// re-initialize matrices
	int lca_level = find_lca_level(s, d);
	int lca = treepath_set[s][lca_level];

	if(!tree[lca].isleaf){

		int stn, dtn, xpos, ypos, ubs, cn;
		ubs = tree[lca].union_borders.size();

		for(int j=0; j<ubs; j++) {
			for(int k=0; k<ubs; k++) {
				if(j==k) tree[lca].dist.a[j][k] = 0;
				else tree[lca].dist.a[j][k] = INT_MAX;
			}
		}

		for(int i=0; i<tree[lca].children.size(); i++) {
			cn = tree[lca].children[i];
			int b1, b2, xpos, xpos1, ypos, ypos1;
			for(int i=0; i<tree[cn].borders.size(); i++) {
				b1 = tree[cn].borders[i];
				xpos = tree[cn].border_in_father[tree[cn].border_in_borders_map[b1]];
				for(int j=0; j<tree[cn].borders.size(); j++) {
					b2 = tree[cn].borders[j];
					ypos = tree[cn].border_in_father[tree[cn].border_in_borders_map[b2]];
					if(tree[cn].isleaf) {
						xpos1 = tree[cn].border_in_union_border[i];
						ypos1 = tree[cn].vertex_in_leafnodes[b2];
						tree[lca].dist.a[xpos][ypos] = tree[cn].dist.a[xpos1][ypos1];
					} else {
						xpos1 = tree[cn].border_in_union_border[i];
						ypos1 = tree[cn].border_in_union_border[j];
						tree[lca].dist.a[xpos][ypos] = tree[cn].dist.a[xpos1][ypos1];
					}
				}
			}
		}

		stn = treepath_set[s][lca_level+1];
		dtn = treepath_set[d][lca_level+1];

		xpos = tree[stn].border_in_father[tree[stn].border_in_borders_map[s]];
		ypos = tree[dtn].border_in_father[tree[dtn].border_in_borders_map[d]];

		if(graph.type==0) {
			tree[lca].dist.a[xpos][ypos] = weight;
			tree[lca].dist.a[ypos][xpos] = weight;
		} else {
			tree[lca].dist.a[xpos][ypos] = weight;
		}
	}

	if(!add_to_uq[lca]){
		uq.push(lca);
	}

	int tn, fn, cn, cnc;
	// push up
	while(!uq.empty()) {
		tn = uq.top();
		uq.pop();

		ulist.push_back(tn);
		updated_tn[tn] = true;

		if(tree[tn].level<nflevel) {
			for(int i=0; i<tree[tn].children.size(); i++) {
				cn = tree[tn].children[i];
				for(auto it : tree[cn].short_cuts) {
					if(it.second>0) {
						updated_shortcuts.insert(it.second-1);
					} else {
						updated_shortcuts.insert(-it.second-1);
					}
				}
			}
		}


		if(tree[tn].isleaf) {
			std::vector<int> minds, results;
			for(int i=0; i<tree[tn].borders.size(); i++) {
				results = graph.dijkstra(tree[tn].borders[i], tree[tn].leafnodes);
				minds.insert(minds.end(), results.begin(), results.end());
			}
			assert(tree[tn].dist.r==tree[tn].borders.size() && tree[tn].dist.c==tree[tn].leafnodes.size() );
			tree[tn].dist.init(tree[tn].borders.size(), tree[tn].leafnodes.size(), minds);
		} else {
			tree[tn].dist.do_floyd();
		}

		fn = tree[tn].father;
		if(fn != -1){
			bool updated = false;
			int b1, b2, xpos, xpos1, ypos, ypos1;
			for(int i=0; i<tree[tn].borders.size(); i++) {
				b1 = tree[tn].borders[i];
				xpos = tree[tn].border_in_father[tree[tn].border_in_borders_map[b1]]; // b1
				for(int j=0; j<tree[tn].borders.size(); j++) {
					b2 = tree[tn].borders[j];
					ypos = tree[tn].border_in_father[tree[tn].border_in_borders_map[b2]]; // b2
					if(tree[tn].isleaf) {
						xpos1 = tree[tn].border_in_union_border[i]; // b1
						ypos1 = tree[tn].vertex_in_leafnodes[b2]; // b2
						if(tree[fn].dist.a[xpos][ypos] != tree[tn].dist.a[xpos1][ypos1]) {
							tree[fn].dist.a[xpos][ypos] = tree[tn].dist.a[xpos1][ypos1];
							updated = true;
						}
					} else {
						xpos1 = tree[tn].border_in_union_border[i]; // b1
						ypos1 = tree[tn].border_in_union_border[j]; // b2
						if(tree[fn].dist.a[xpos][ypos] != tree[tn].dist.a[xpos1][ypos1]) {
							tree[fn].dist.a[xpos][ypos] = tree[tn].dist.a[xpos1][ypos1];
							updated = true;
						}
					}
				}
			}

			if(updated) {
				if(!add_to_uq[tn]){
					uq.push(tn);
				}
				if(!add_to_dq[tn]){
					dq.push(tn);
				}
			}
		} else {
			if(!add_to_dq[tn]){
				dq.push(tn);
			}
		}
	}

	// push down
	while(!dq.empty()) {
		tn = dq.top();
		dq.pop();

		dlist.push_back(tn);

		if(!tree[tn].isleaf) {
			for(int i=0; i<tree[tn].children.size(); i++) {
				cn = tree[tn].children[i];
				bool dirty = false;
				int b1, b2, xpos, xpos1, ypos, ypos1;
				for(int i=0; i<tree[cn].borders.size(); i++) {
					b1 = tree[cn].borders[i];
					xpos = tree[cn].border_in_father[tree[cn].border_in_borders_map[b1]]; // b1
					for(int j=0; j<tree[cn].borders.size(); j++) {
						b2 = tree[cn].borders[j];
						ypos = tree[cn].border_in_father[tree[cn].border_in_borders_map[b2]]; // b2
						if(tree[cn].isleaf) {
							xpos1 = tree[cn].border_in_union_border[i];  // b1
							ypos1 = tree[cn].vertex_in_leafnodes[b2];  // b2
							if(tree[cn].dist.a[xpos1][ypos1] != tree[tn].dist.a[xpos][ypos]) {
								tree[cn].dist.a[xpos1][ypos1] = tree[tn].dist.a[xpos][ypos];
								dirty = true;
							}
						} else {
							xpos1 = tree[cn].border_in_union_border[i];  // b1
							ypos1 = tree[cn].border_in_union_border[j];  // b2
							if(tree[cn].dist.a[xpos1][ypos1] != tree[tn].dist.a[xpos][ypos]) {
								tree[cn].dist.a[xpos1][ypos1] = tree[tn].dist.a[xpos][ypos];
								dirty = true;
							}
						}
					}
				}

				if(dirty) {
					updated_tn[cn] = true;
					if(tree[cn].level<nflevel) {
						for(int i=0; i<tree[cn].children.size(); i++) {
							cnc = tree[cn].children[i];
							for(auto it : tree[cnc].short_cuts) {
								if(it.second>0) {
									updated_shortcuts.insert(it.second-1);
								} else {
									updated_shortcuts.insert(-it.second-1);
								}
							}
						}
					}

					if(tree[cn].isleaf){
						std::vector<int> minds, results;
						for(int i=0; i<tree[cn].borders.size(); i++) {
							results = graph.dijkstra(tree[cn].borders[i], tree[cn].leafnodes);
							minds.insert(minds.end(), results.begin(), results.end());
						}
						assert(tree[cn].dist.r==tree[cn].borders.size() && tree[cn].dist.c==tree[cn].leafnodes.size() );
						tree[cn].dist.init(tree[cn].borders.size(), tree[cn].leafnodes.size(), minds);

					} else {
						tree[cn].dist.do_floyd();
						if(!add_to_dq[cn]){
							dq.push(cn);
						}
					}
				}
			}
		}
	}
	
	gettimeofday( &tv, NULL );
    t2 = tv.tv_sec * 1000000 + tv.tv_usec ;
    update_time = t2 - t1;


	if(!with_shortcut){
		if(show_info) {
			int updated_tn_num = std::count(updated_tn.begin(), updated_tn.end(), true);
			int updated_st_num = updated_shortcuts.size();
			std::cout << "total update tree node ratio: " << ((float) updated_tn_num) / treesize*100 << "%.";
			std::cout << " use time: " << update_time << " us." << std::endl;
		}
		return true;
	}

    gettimeofday( &tv, NULL );
    t1 = tv.tv_sec * 1000000 + tv.tv_usec ;

	update_shortcut(updated_shortcuts);

	gettimeofday( &tv, NULL );
    t2 = tv.tv_sec * 1000000 + tv.tv_usec ;
	sc_update_time = t2 - t1;

	if(show_info) {
		int updated_tn_num = std::count(updated_tn.begin(), updated_tn.end(), true);
		int updated_st_num = updated_shortcuts.size();
		std::cout << "total update tree node ratio: " << ((float) updated_tn_num) / treesize*100 << "%.";
		std::cout << " use time: " << update_time << " us." << std::endl;
		std::cout << "total update shortcut ratio: " << ((float) updated_st_num) / short_cut_matrix_in_level.size()*100 << "%." ;
		std::cout << " use time: " << sc_update_time << " us." << std::endl;
	}
    
	return true;

}


bool SCGTree::update_shortcut(std::set<int> &updated_shortcuts){

	std::map<int, std::vector< std::vector<int> > > temp_map;
	std::set<int> temp_saved;
	int snode, dnode;
	std::vector<int> minds, sminds, dminds; 

	for(std::set<int>::iterator it=updated_shortcuts.begin(); it!=updated_shortcuts.end(); it++) {
		minds.clear();
		std::pair<int, int> p = short_cut_pairs_in_level[*it];
		snode = p.first;
		dnode = p.second;
		if(temp_saved.find(snode) == temp_saved.end()) {
			cal_temp_dist(snode, temp_map);
			temp_saved.insert(snode);
		}

		if(temp_saved.find(dnode) == temp_saved.end()) {
			cal_temp_dist(dnode, temp_map);
			temp_saved.insert(dnode);
		}

		int dlevel=-1, ix, iy;

		// find lca
		while(tree[snode].father != tree[dnode].father) {
			snode = tree[snode].father;
			dnode = tree[dnode].father;
			dlevel += 1;
		}
		int lca = tree[snode].father;

		for(int j=0; j<tree[p.first].borders.size(); j++) {
			int s = tree[p.first].borders[j];
			assert(temp_map.find(s) != temp_map.end() && dlevel<temp_map[s].size());
			sminds = temp_map[s][dlevel];
			for(int k=0; k<tree[p.second].borders.size(); k++) {
				int d = tree[p.second].borders[k];
				assert(temp_map.find(d) != temp_map.end() && dlevel<temp_map[d].size());
				dminds = temp_map[d][dlevel];
				int res = INT_MAX;
				for(int m=0; m<sminds.size(); m++) {
					ix = tree[snode].border_in_father[m];
					for(int n=0; n<dminds.size(); n++) {
						iy = tree[dnode].border_in_father[n];
						int dis = tree[lca].dist.a[ix][iy];
						if (res > sminds[m] + dis + dminds[n]) {
							res = sminds[m] + dis + dminds[n];
						}
					}
				}
				minds.push_back(res);
			}
		}

		short_cut_matrix_in_level[*it].init(tree[p.first].borders.size(), tree[p.second].borders.size(), minds);
	}
	cal_short_cuts_dis();

	return true;

}


// find lca's level of s and d
int SCGTree::find_lca_level(int s, int d){

	int same_level=-1;
	std::vector<int> spath = treepath_set[s];
	std::vector<int> dpath = treepath_set[d];

	if(spath.size()>dpath.size()){
		swap(spath, dpath);
	}

	for(int i=0; i<spath.size(); i++) {
		if(spath[i]==dpath[i]){
			same_level=i;
		}
	}

	return same_level;
}


// find lca's level of s and d
int SCGTree::find_lca(int s, int d){

	int same_level=-1;
	std::vector<int> spath = treepath_set[s];
	std::vector<int> dpath = treepath_set[d];

	if(spath.size()>dpath.size()){
		swap(spath, dpath);
	}

	for(int i=0; i<spath.size(); i++) {
		if(spath[i]==dpath[i]){
			same_level=i;
		}
	}

	return treepath_set[s][same_level];
}


SpspResult SCGTree::spsp_path_query(int s, int d, bool with_shortcut){
	if(s==d) {
		SpspResult sr{0, std::vector<int>{s, d}, std::vector<int>{0}};
		return sr;
	}

	std::vector<int> spath = treepath_set[s];
	std::vector<int> dpath = treepath_set[d];
    bool use_short_cut = false;

    if(spath.back()==dpath.back()) {
        std::vector<int> cands;
        cands.push_back(d);
        std::vector<int> result = graph.dijkstra( s, cands);
		SpspResult sr{result[0], std::vector<int>{s, d}, std::vector<int>{result[0]}};

		return sr;

    }  else {

		int same_level=-1;

		if(spath.size()>dpath.size()){
			swap(spath, dpath);
		}

		for(int i=0; i<spath.size(); i++) {
			if(spath[i]==dpath[i]){
				same_level=i;
			}
		}
		
		if(with_shortcut){

			if(tree[spath[flevel]].short_cuts.find(dpath[flevel]) != tree[spath[flevel]].short_cuts.end()) {
					same_level = flevel-1;
					use_short_cut = true;
			}

			if(nflevel != -1 && tree[spath[nflevel]].short_cuts.find(dpath[nflevel]) != tree[spath[nflevel]].short_cuts.end()) {
					same_level = nflevel-1;
					use_short_cut = true;
			}
		}

		spath = treepath_set[s];
		dpath = treepath_set[d];

        std::vector<int> s_search_node_path; // s -> same level tree node path
		std::vector<int> d_search_node_path; // d -> same level tree node path

		s_search_node_path.push_back(-1);
        for(int l=border_max_level[s]-1; l>same_level; l--) {
			s_search_node_path.push_back(spath[l]);
        }

		d_search_node_path.push_back(-1);
        for(int l=border_max_level[d]-1; l>same_level; l--) {
			d_search_node_path.push_back(dpath[l]);
        }

		std::vector<int> s_dist, d_dist, dist;
		int x_border_size, y_border_size, ix, iy, weight, index;
		bool s_path_only_one_level, d_path_only_one_level;

		//record path
		std::vector<std::vector<int> > s_dp_path, d_dp_path;
		std::vector<int> s_dp_path_temp, d_dp_path_temp;

		// s_dp_result
		s_dist.clear();
		s_dist.push_back(0);

		s_dp_path.push_back(std::vector<int>{s});

		if(s_search_node_path.size()>1) {
			s_dist.clear();
			s_dist.resize(tree[s_search_node_path[1]].borders.size());

			// 添加路径
			s_dp_path_temp.clear();
			s_dp_path_temp.resize(tree[s_search_node_path[1]].borders.size(), 0);
			s_dp_path.push_back(s_dp_path_temp);

			if(tree[s_search_node_path[1]].isleaf) {
				index = tree[s_search_node_path[1]].vertex_in_leafnodes[s];
				for(int i=0; i<tree[s_search_node_path[1]].borders.size(); i++){
					s_dist[i] = tree[s_search_node_path[1]].dist.a[i][index];
				}
			} else {
				ix = tree[spath[border_max_level[s]]].border_in_father[border_max_level_index[s]];
				for(int i=0; i<tree[s_search_node_path[1]].borders.size(); i++){
					iy = tree[s_search_node_path[1]].border_in_union_border[i];
					s_dist[i] = tree[s_search_node_path[1]].dist.a[ix][iy];
				}
			}
		}

		if(s_search_node_path.size()>2) {
			for(int i=2; i<s_search_node_path.size(); i++) {
				x_border_size = tree[s_search_node_path[i-1]].borders.size();
				y_border_size = tree[s_search_node_path[i]].borders.size();

				dist.clear();
				dist.resize(y_border_size, INT_MAX);

				s_dp_path_temp.clear();
				s_dp_path_temp.resize(y_border_size);

				int up;
				for(int j=0; j<y_border_size; j++){
					for(int k=0; k<x_border_size; k++) {
						ix = tree[s_search_node_path[i-1]].border_in_father[k];
						iy = tree[s_search_node_path[i]].border_in_union_border[j];
						weight = tree[s_search_node_path[i]].dist.a[ix][iy];
						if(dist[j] > s_dist[k] + weight){
							dist[j] = s_dist[k] + weight;
							up = k;
						}
					}
					s_dp_path_temp[j] = up;
				}
				s_dist = dist;
				s_dp_path.push_back(s_dp_path_temp);
			}
		}

		d_dist.clear();
		d_dist.push_back(0);
		d_dp_path.push_back(std::vector<int>{d});

		if(d_search_node_path.size()>1) {
			d_dist.clear();
			d_dist.resize(tree[d_search_node_path[1]].borders.size());

			d_dp_path_temp.clear();
			d_dp_path_temp.resize(tree[d_search_node_path[1]].borders.size(), 0);
			d_dp_path.push_back(d_dp_path_temp);



			if(tree[d_search_node_path[1]].isleaf) {
				index = tree[d_search_node_path[1]].vertex_in_leafnodes[d];
				for(int i=0; i<tree[d_search_node_path[1]].borders.size(); i++){
					d_dist[i] = tree[d_search_node_path[1]].dist.a[i][index];
				}
			} else { 
				ix = tree[dpath[border_max_level[d]]].border_in_father[border_max_level_index[d]];
				for(int i=0; i<tree[d_search_node_path[1]].borders.size(); i++){
					iy = tree[d_search_node_path[1]].border_in_union_border[i];
					d_dist[i] = tree[d_search_node_path[1]].dist.a[ix][iy];
				}
			}
		}

		if(d_search_node_path.size()>2) {
			for(int i=2; i<d_search_node_path.size(); i++) {
				
				x_border_size = tree[d_search_node_path[i-1]].borders.size();
				y_border_size = tree[d_search_node_path[i]].borders.size();

				dist.clear();
				dist.resize(y_border_size, INT_MAX);

				d_dp_path_temp.clear();
				d_dp_path_temp.resize(y_border_size);

				int up;
				for(int j=0; j<y_border_size; j++){
					for(int k=0; k<x_border_size; k++) {
						ix = tree[d_search_node_path[i-1]].border_in_father[k];
						iy = tree[d_search_node_path[i]].border_in_union_border[j];

						weight = tree[d_search_node_path[i]].dist.a[ix][iy];


						if(dist[j] > d_dist[k] + weight){
							dist[j] = d_dist[k] + weight;
							up = k;
						}
					}

					// d_dp_path_temp.push_back(up);
					d_dp_path_temp[j] = up;

				}
				d_dist = dist;
				d_dp_path.push_back(d_dp_path_temp);
			}
		}
		int min_dis = INT_MAX;
		x_border_size = s_dist.size();
		y_border_size = d_dist.size();

		int p1, p2, w;
		s_path_only_one_level = s_search_node_path.back()== -1 ? true : false;
		d_path_only_one_level = d_search_node_path.back()== -1 ? true : false;

        if(use_short_cut) {
            int matrix_index = tree[spath[same_level+1]].short_cuts[dpath[same_level+1]];
            if(!s_path_only_one_level && !d_path_only_one_level) {
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        ix = i;
                        iy = j;
                        if(matrix_index>0) {
                            weight = short_cut_matrix_in_level[matrix_index-1].a[ix][iy];
                        } else {
                            weight = short_cut_matrix_in_level[-matrix_index-1].a[iy][ix];
                        }
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
							p1 = i;
							p2 = j;
							w = weight;
                        }
                    }
                }
            } else if (s_path_only_one_level && !d_path_only_one_level) {
				ix = tree[spath[same_level+1]].border_in_borders_map[s];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        iy = j;
                        if(matrix_index>0) {
                            weight = short_cut_matrix_in_level[matrix_index-1].a[ix][iy];
                        } else {
                            weight = short_cut_matrix_in_level[-matrix_index-1].a[iy][ix];
                        }
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
							p1 = i;
							p2 = j;
							w = weight;
                        }
                    }
                }
            } else if (!s_path_only_one_level && d_path_only_one_level) {
				iy = tree[dpath[same_level+1]].border_in_borders_map[d];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        ix = i;
                        if(matrix_index>0) {
                            weight = short_cut_matrix_in_level[matrix_index-1].a[ix][iy];
                        } else {
                            weight = short_cut_matrix_in_level[-matrix_index-1].a[iy][ix];
                        }

                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
							p1 = i;
							p2 = j;
							w = weight;
                        }
                    }
                }
            } else if (s_path_only_one_level && d_path_only_one_level) {
				ix = tree[spath[same_level+1]].border_in_borders_map[s];
				iy = tree[dpath[same_level+1]].border_in_borders_map[d];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        if(matrix_index>0) {
                            weight = short_cut_matrix_in_level[matrix_index-1].a[ix][iy];
                        } else {
                            weight = short_cut_matrix_in_level[-matrix_index-1].a[iy][ix];
                        }
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
							p1 = i;
							p2 = j;
							w = weight;
                        }
                    }
                }
            }
        } else {
            if(!s_path_only_one_level && !d_path_only_one_level) {
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        ix = tree[s_search_node_path.back()].border_in_father[i];
                        iy = tree[d_search_node_path.back()].border_in_father[j];
                        weight = tree[dpath[same_level]].dist.a[ix][iy];
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
							p1 = i;
							p2 = j;
							w = weight;
                        }
                    }
                }
            } else if (s_path_only_one_level && !d_path_only_one_level) {
				ix = tree[spath[same_level+1]].border_in_father[tree[spath[same_level+1]].border_in_borders_map[s]];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        iy = tree[d_search_node_path.back()].border_in_father[j];
                        weight = tree[dpath[same_level]].dist.a[ix][iy];
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
							p1 = i;
							p2 = j;
							w = weight;
                        }
                    }
                }
            } else if (!s_path_only_one_level && d_path_only_one_level) {
				iy = tree[dpath[same_level+1]].border_in_father[tree[dpath[same_level+1]].border_in_borders_map[d]];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        ix = tree[s_search_node_path.back()].border_in_father[i];
                        weight = tree[dpath[same_level]].dist.a[ix][iy];
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
							p1 = i;
							p2 = j;
							w = weight;
                        }
                    }
                }
            } else if (s_path_only_one_level && d_path_only_one_level) {
				ix = tree[spath[same_level+1]].border_in_father[tree[spath[same_level+1]].border_in_borders_map[s]];
				iy = tree[dpath[same_level+1]].border_in_father[tree[dpath[same_level+1]].border_in_borders_map[d]];
                for(int i=0; i<x_border_size; i++) {
                    for(int j=0; j<y_border_size; j++) {
                        weight = tree[dpath[same_level]].dist.a[ix][iy];
                        if(min_dis > s_dist[i] + d_dist[j] + weight) {
                            min_dis = s_dist[i] + d_dist[j] + weight;
							p1 = i;
							p2 = j;
							w = weight;
                        }
                    }
                }
            }
        }

		SpspResult sr;
		sr.dis = min_dis;

		std::vector<int> s_path, s_dis;
		std::vector<int> d_path, d_dis;

		assert(s_search_node_path.size() == s_dp_path.size());

		if(s_search_node_path.size()>1) {
			s_path.push_back(tree[s_search_node_path.back()].borders[p1]);
			for(int i=s_search_node_path.size()-1; i>1; i--) {
				p1 = s_dp_path[i][p1];
				s_path.push_back(tree[s_search_node_path[i-1]].borders[p1]);
			}
		}

		s_path.push_back(s);
		std::reverse(s_path.begin(), s_path.end());

		if(s_search_node_path.size()>1) {
			int bpos, xpos, ypos, tn;
			for(int i=s_search_node_path.size()-1; i>0; i--) {
				if(i>=2){
					bpos = tree[s_search_node_path[i]].border_in_borders_map[s_path[i]];
					xpos = tree[s_search_node_path[i]].border_in_union_border[bpos];
					bpos = tree[s_search_node_path[i-1]].border_in_borders_map[s_path[i-1]];
					ypos = tree[s_search_node_path[i-1]].border_in_father[bpos];
					s_dis.push_back(tree[s_search_node_path[i]].dist.a[ypos][xpos]);
				} else if(i==1) {
					tn = s_search_node_path[i];
					if(tree[tn].isleaf){
						int xpos, ypos;
						xpos = tree[tn].border_in_borders_map[s_path[i]];
						ypos = tree[tn].vertex_in_leafnodes[s];
						s_dis.push_back(tree[tn].dist.a[xpos][ypos]);
					} else {
						int bpos, xpos, ypos;
						bpos = tree[tn].border_in_borders_map[s_path[i]];
						xpos = tree[tn].border_in_union_border[bpos];
						ypos = tree[spath[border_max_level[s]]].border_in_father[border_max_level_index[s]];
						s_dis.push_back(tree[tn].dist.a[ypos][xpos]);
					}
				}
			}
		}

		std::reverse(s_dis.begin(), s_dis.end());

		assert(d_search_node_path.size() == d_dp_path.size());

		if(d_search_node_path.size()>1) {
			d_path.push_back(tree[d_search_node_path.back()].borders[p2]);
			for(int i=d_search_node_path.size()-1; i>1; i--) {
				p2 = d_dp_path[i][p2];
				d_path.push_back(tree[d_search_node_path[i-1]].borders[p2]);
			}
		}
		d_path.push_back(d);
		std::reverse(d_path.begin(), d_path.end());

		if(d_search_node_path.size()>1) {
			int bpos, xpos, ypos, tn;
			for(int i=d_search_node_path.size()-1; i>0; i--) {
				if(i>=2){
					bpos = tree[d_search_node_path[i]].border_in_borders_map[d_path[i]];
					xpos = tree[d_search_node_path[i]].border_in_union_border[bpos];
					bpos = tree[d_search_node_path[i-1]].border_in_borders_map[d_path[i-1]];
					ypos = tree[d_search_node_path[i-1]].border_in_father[bpos];
					d_dis.push_back(tree[d_search_node_path[i]].dist.a[ypos][xpos]);
				} else if(i==1){
					tn = d_search_node_path[i];
					if(tree[tn].isleaf){
						int xpos, ypos;
						xpos = tree[tn].border_in_borders_map[d_path[i]];
						ypos = tree[tn].vertex_in_leafnodes[d];
						d_dis.push_back(tree[tn].dist.a[xpos][ypos]);
					} else {
						int bpos, xpos, ypos;
						bpos = tree[tn].border_in_borders_map[d_path[i]];
						xpos = tree[tn].border_in_union_border[bpos];
						ypos = tree[dpath[border_max_level[d]]].border_in_father[border_max_level_index[d]];
						d_dis.push_back(tree[tn].dist.a[ypos][xpos]);
					}
				}
			}
		}

		std::reverse(d_path.begin(), d_path.end());

		sr.path.insert(sr.path.begin(), s_path.begin(), s_path.end());
		sr.path.insert(sr.path.end(), d_path.begin(), d_path.end());

		sr.path_dis.insert(sr.path_dis.begin(), s_dis.begin(), s_dis.end());
		sr.path_dis.insert(sr.path_dis.end(), w);
		sr.path_dis.insert(sr.path_dis.end(), d_dis.begin(), d_dis.end());

		return sr;
    }

}


// tree path recovery
std::vector<int> SCGTree::tree_path_recovery(int s, int d, SpspResult sr){
	std::vector<int> path;
	for(int i=0; i<sr.path.size()-1; i++) {
		path_recovery(sr.path[i], sr.path[i+1], sr.path_dis[i], path);
	}

	return path;
}


// path recovery
std::vector<int> SCGTree::path_recovery(int v1, int v2, int sp_dis, std::vector<int> &path){
	if(v1==v2){
		path.push_back(v1);
		return path;
	}

	if(graph.edge_exist(v1, v2)){
		path.push_back(v1);
		path.push_back(v2);
		return path;
	}

	int ln1 = treepath_set[v1].back();
	int ln2 = treepath_set[v2].back();

	if(ln1==ln2){
		int lca = find_lca(v1, v2);
		PathRecoveryResult prr = find_border(v1, v2, sp_dis, lca);
		while(prr.status==0){
			lca = tree[lca].father;
			if(lca==-1){
				break;
			}
			prr = find_border(v1, v2, sp_dis, lca);
		}

		if(prr.status==1){
			path_recovery(v1, prr.b, prr.sp_dis1, path);
			path_recovery(prr.b, v2, prr.sp_dis2, path);
			return path;

		} else if (prr.status==0) {
			graph.dijkstra_with_path(v1, v2, path);
			return path;
		}
	} else {
		int lca = find_lca(v1, v2);
		PathRecoveryResult prr = find_border(v1, v2, sp_dis, lca);
		while(prr.status==0){
			lca = tree[lca].father;
			if(lca==-1){
				break;
			}
			prr = find_border(v1, v2, sp_dis, lca);
		}
 
		if(prr.status==1){
			path_recovery(v1, prr.b, prr.sp_dis1, path);
			path_recovery(prr.b, v2, prr.sp_dis2, path);
			return path;
		} else {
			path.push_back(v1);
			path.push_back(v2);

			SpspResult sub_result;
			std::vector<int> sub_path;

			sub_result = spsp_path_query(v1, v2, false);
			sub_path = tree_path_recovery(v1, v2, sub_result);
			path.insert(path.end(), sub_path.begin(), sub_path.end());

			return path;
		}
	}
}


// find border in tn, satisfy sp(v1,v2) = sp(v1,b) + sp(b,v2)
PathRecoveryResult SCGTree::find_border(int v1, int v2, int sp_dis, int tn){

	int vi1, vi2;

	std::vector<int>::iterator it = std::find(tree[tn].union_borders.begin(), tree[tn].union_borders.end(), v1);

	if(it!=tree[tn].union_borders.end() && *it==v1){
		vi1 = it - tree[tn].union_borders.begin();
	} else {
		PathRecoveryResult prr;
		prr.status = 0;
		return prr;
	}

	it = std::find(tree[tn].union_borders.begin(), tree[tn].union_borders.end(), v2);
	if(it!=tree[tn].union_borders.end() && *it==v2){
		vi2 = it - tree[tn].union_borders.begin();
	} else {
		PathRecoveryResult prr;
		prr.status = 0;
		return prr;
	}

	int dis, b;
	for(int i=0; i<tree[tn].union_borders.size(); i++){
		b = tree[tn].union_borders[i];
		if(b!=v1 && b!=v2){
			if(sp_dis == tree[tn].dist.a[vi1][i]+tree[tn].dist.a[i][vi2]){
				PathRecoveryResult prr = {1, v1, v2, b, tree[tn].dist.a[vi1][i], tree[tn].dist.a[i][vi2]};
				return prr;
			}
		}
	}

	PathRecoveryResult prr;
	prr.status = 0;
	return prr;
}
