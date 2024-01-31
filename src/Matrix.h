#pragma once

#include<vector>


class Matrix
{
	public:

		int r;
		int c;
		int mn;
		int **a;
		int **midv; 
		
		int vmin;

		Matrix();
        
		~Matrix();

		void init(int row, int col, std::vector<int> mind);

		void init_with_midv(int row, int col, std::vector<int> mind, std::vector<int> mid);

		void init(int n);

		void clear();

		Matrix(const Matrix &om);

		Matrix& operator=(const Matrix &om);

		void update_value(int x, int y, int value);

		void do_floyd();

		int min_value();

};

