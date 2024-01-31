#pragma once

#include<vector>
#include<unordered_map>


class MinHeap 
{
	public:
		MinHeap();

		int n; 
		std::vector<std::pair<int,int> > a;
		std::unordered_map<int, int> vertinheap;

		std::pair<int, int> get(int i);

		int find(int id);

		bool empty();


		std::pair<int, int> top();

		std::pair<int, int> extract_top();

		void fix_down(int i);

		void insert(std::pair<int,int> key);

		void decrease(int i, std::pair<int, int> key);

		void init();

		void show();

};
