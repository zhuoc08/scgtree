#pragma once

#include<vector>


class Vertex
{
	public:
		double x,y; 
		std::vector<int> adjnodes; 
		std::vector<int> adjweight;

		Vertex();
		Vertex(double x, double y);

};


