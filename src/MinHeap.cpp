#include<iostream>
#include<cassert>
#include<vector>
#include<unordered_map>
#include<limits.h>
#include<math.h>


#include"MinHeap.h"



MinHeap::MinHeap(){init();}


std::pair<int, int> MinHeap::get(int i) {
    return a[i];
}

int MinHeap::find(int id) {
    if(vertinheap.find(id) != vertinheap.end()) {
        return vertinheap[id];
    } else {
        return -1;
    }
}


bool MinHeap::empty() {
    if(n==0) {
        return true;
    } else {
        return false;
    }
}


std::pair<int, int> MinHeap::top(){
    return a[1];
}


std::pair<int, int> MinHeap::extract_top() {
    assert(n>0);
    std::pair<int, int> min = a[1];
    a[1] = a[n];
    a.pop_back();
    vertinheap.erase(min.first);
    vertinheap[a[1].first] = 1;
    n --;
    fix_down(1);
    return min;
}


void MinHeap::fix_down(int i) {
    int l = i*2;
    int r = i*2 + 1;
    int next = i;
    if(l<=n && a[l].second < a[i].second){
        next = l;
    }
    if(r<=n && a[r].second < a[next].second){
        next = r;
    }
    if(next != i) {
        swap(a[i], a[next]);
        vertinheap[a[i].first] = i;
        vertinheap[a[next].first] = next;
        fix_down(next);
    }
}


void MinHeap::insert(std::pair<int,int> key) {
    n ++;
    a.push_back(std::make_pair(key.first, INT_MAX));
    vertinheap[key.first] = n;
    decrease(n, key);
}


void MinHeap::decrease(int i, std::pair<int, int> key) {
    assert(key.first == a[i].first);
    assert(key.second < a[i].second);
    a[i] = key;

    while(i>1 && (a[i].second < a[floor(i/2)].second)) {
        swap(a[i], a[floor(i/2)]);
        vertinheap[a[i].first] = i;
        vertinheap[a[floor(i/2)].first] = floor(i/2);
        i = floor(i/2);
    }

}


void MinHeap::init()
{
    a.clear();
    a.push_back(std::make_pair(0, 0));
    n = 0;
    vertinheap.clear();
}


void MinHeap::show() {
    std::cout << "length:" << n << std::endl;
    for(int i=1; i<n+1; i++) {
        std::cout << a[i].first << ": " << a[i].second << ", " ;
    }
    std::cout << std::endl;
}

