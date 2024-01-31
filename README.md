# SCG-tree

Shortcut Enhanced Graph Hierarchy Tree for Efficient Spatial Queries on Massive Road Networks

###  Datasets

- [DIMACS Datasets](https://www.diag.uniroma1.it//challenge9/download.shtml)
- An example of inputs
  - BAY.co (vertices) : first line contains the total number of vertices, following lines contain three elements (vertex id, Latitude, Longitude)
  - BAY.d (edges) : first line contains the total number of vertices and edges, following lines contain three elements (start vertex id, end vertex id, distance)

### Graph Partition

- Graph partition process is performed by [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/download), please install Metis first.

### Quick Test

- You can quickly construct and save BAY.scg by using BAY.co and BAY.d, and perform SPSP, kNN, and range queries on BAY.scg.

make test

./test

