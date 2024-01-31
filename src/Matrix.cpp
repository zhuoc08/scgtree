#include<iostream>
#include<cassert>
#include<vector>
#include<set>
#include<unordered_map>
#include<limits.h>
#include <malloc.h>

#include"Matrix.h"
#include"MinHeap.h"


Matrix::Matrix():r(1),c(1),mn(1){

    a = new int*[r];
    for(int i=0;i<r;i++)a[i]=new int[c];

    midv = new int*[mn];
    for(int i=0;i<mn;i++)midv[i]=new int[mn];

}


Matrix::~Matrix(){
    clear();
}


void Matrix::clear()
{
    if(a && r>0){
        for(int i=0;i<r;i++)delete [] a[i];
        delete [] a;
    }

    if(midv && mn>0){
        for(int i=0;i<mn;i++)delete [] midv[i];
        delete [] midv;
    }

    return;
}


void Matrix::init(int row, int col, std::vector<int> mind)
{

    if(row*col == mind.size()) {
        if(a){
            for(int i=0;i<r;i++)delete [] a[i];
            delete [] a;
        }

        r = row;
        c = col;

        a = new int*[r];
        for(int i=0;i<r;i++)a[i]=new int[c];

        mn = 1;
        midv = new int*[mn];
        for(int i=0;i<mn;i++)midv[i]=new int[mn];

        int min=INT_MAX;
        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                a[i][j]=mind[i*c + j];
                if(a[i][j]!=0 && a[i][j]<min) {
                    min = a[i][j];
                }
            }
        }
        vmin = min;
    } else {
        std::cout << "matrix size is not equal!" << std::endl;
    }

    return;

}


void Matrix::init_with_midv(int row, int col, std::vector<int> mind, std::vector<int> mid)
{

    assert(row*col == mind.size());
    if(row*col == mind.size()) {

        if(a){
            for(int i=0;i<r;i++)delete [] a[i];
            delete [] a;
        }
        if(midv){
            for(int i=0;i<mn;i++)delete [] midv[i];
            delete [] midv;
        }

        r = row;
        c = col;

        a = new int*[r];
        for(int i=0;i<r;i++)a[i]=new int[c];

        int min=INT_MAX;
        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                a[i][j]=mind[i*c + j];
                if(a[i][j]!=0 && a[i][j]<min) {
                    min = a[i][j];
                }
            }
        }

        mn = col;
        midv = new int*[mn];
        for(int i=0;i<mn;i++)midv[i]=new int[mn];
        for(int i=0;i<mn;i++){
            for(int j=0;j<mn;j++){
                midv[i][j] = mid[i*mn + j];
            }
        }

        vmin = min;
    } else {
        std::cout << "matrix size is not equal!" << std::endl;
    }

    return;

}


void Matrix::init(int n){

    if(a){
        for(int i=0;i<r;i++)delete [] a[i];
        delete [] a;
    }
    if(midv){
        for(int i=0;i<mn;i++)delete [] midv[i];
        delete [] midv;
    }

    r = n;
    c = n;
    a = new int*[r];
    for(int i=0;i<r;i++)a[i]=new int[c];


    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            if(i==j){
                a[i][j]=0;
            } else {
                a[i][j] = INT_MAX;
            }
        }
    }

    mn = n;
    midv = new int*[mn];
    for(int i=0;i<mn;i++)midv[i]=new int[mn];

    for(int i=0; i<mn; i++) {
        for(int j=0; j<mn; j++) {
            if(i==j){
                midv[i][j]= -1;
            } else {
                midv[i][j] = -4;
            }
        }
    }
    return;
}


Matrix::Matrix(const Matrix &om){
    r = om.r;
    c = om.c;
    a = new int*[r];
    for(int i=0;i<r;i++)a[i]=new int[c];

    for(int i=0; i<r; i++) {
        for(int j=0; j<c; j++) {
            a[i][j] = om.a[i][j];
        }
    }

    mn = om.c;
    midv = new int*[mn];
    for(int i=0;i<mn;i++)midv[i]=new int[mn];

    for(int i=0; i<mn; i++) {
        for(int j=0; j<mn; j++) {
            midv[i][j] = om.midv[i][j];
        }
    }

}


Matrix& Matrix::operator=(const Matrix &om){
    r = om.r;
    c = om.c;
    a = new int*[r];
    for(int i=0;i<r;i++)a[i]=new int[c];

    for(int i=0; i<r; i++) {
        for(int j=0; j<c; j++) {
            a[i][j] = om.a[i][j];
        }
    }

    mn = om.c;
    midv = new int*[mn];
    for(int i=0;i<mn;i++)midv[i]=new int[mn];

    for(int i=0; i<mn; i++) {
        for(int j=0; j<mn; j++) {
            midv[i][j] = om.midv[i][j];
        }
    }

    return *this;
}


void Matrix::update_value(int x, int y, int value){
    if(r<c && y<c){
        a[x][y] = value;
    } else {
        std::cout << "position error!" << std::endl;
    }
}


void Matrix::do_floyd(){

    if(r!=c) {
        std::cout << "matrix row and column not equal!" << std::endl;
        return;
    }
	for(int k=0; k<r; k++) {
		for(int i=0; i<r; i++) {
			for(int j=0; j<r; j++) {
				if(a[i][k]==INT_MAX || a[k][j]==INT_MAX) {
					continue;
				}
				if(a[i][j] > a[i][k] + a[k][j]) {
					a[i][j] = a[i][k] + a[k][j]; 
                    midv[i][j] = k; 
				}
			}
		}
	}
}


int Matrix::min_value(){
    return vmin;
}
