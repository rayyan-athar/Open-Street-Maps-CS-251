// graph.h 
// Rayyan Athar
//
// Basic graph class using adjacency list representation.  
//
// Adam T Koehler, PhD
// University of Illinois Chicago
// CS 251, Fall 2023
//
// Project Original Variartion By:
// Joe Hummel, PhD
// University of Illinois at Chicago

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <algorithm>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
  // Adjacency list representation of the graph. With vertices and edges.
  vector<VertexT> Vertices;
  map<VertexT, map<VertexT, WeightT> > AdjacencyList;
  map<VertexT, set<VertexT> > NeighboursList;
  int numEdges;
 
 public:

  // Constructor: Constructs an empty graph 
  graph() {
    numEdges = 0;
  }

  // Clear: Deletes all vertices and edges
  void clear() {
    Vertices.clear();
    AdjacencyList.clear();
  }

  // NumVertices: Returns the # of vertices currently in the graph.
  int NumVertices() const {
    return static_cast<int>(this->Vertices.size());
  }

  // NumEdges: Returns the # of edges currently in the graph.
  int NumEdges() const {
    return numEdges;
  }

  // addVertex: Adds the vertex v to the graph returns true.  If the  vertex already exists in the graph, then false is returned.
  bool addVertex(VertexT v) {
    // Check to see if the vertex already exists in the graph and return false if it does
    if(AdjacencyList.count(v) == 1) {
      return false;
    }

    // Add the vertex to the graph and return true, and the vertex to the set of vertices
    set<VertexT> S;
    NeighboursList.emplace(v, S);
    Vertices.push_back(v);
    AdjacencyList[v];
    return true;
  }

  // addEdge: Adds the edge (from, to, weight) to the graph, and returns true.  If the edge already exists in the graph, then false is returned.
  bool addEdge(VertexT from, VertexT to, WeightT weight) {

    // Check to see if the from or to vertex already exist in the graph and return false if it does not exist
    if(AdjacencyList.count(from) == 0 || AdjacencyList.count(to) == 0) {
      return false;
    }

    // Add the edge to the graph and return true
    if(AdjacencyList[from].count(to) == 0) { 
      AdjacencyList.at(from).emplace(to, weight);
      NeighboursList.at(from).emplace(to);
      numEdges++;
    }
    else {
      AdjacencyList.at(from).at(to) = weight;
    }
    return true;
  }

  // getWeight: Returns the weight associated with a given edge.  If the edge exists, the weight is returned via the reference parameter and true is returned.  If the edge does not exist, the weight parameter is unchanged and false is returned.
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {

    // Check to see if the from or to vertex already exist in the graph and return false if it does not exist
    if(AdjacencyList.count(from) == 0 || AdjacencyList.count(to) == 0) {
      return false;
    }

    // If the edge exists, update the weight and return true
    if(AdjacencyList.at(from).count(to) == 1) { 
      weight = AdjacencyList.at(from).at(to);
      return true;
    }

    return false;    
  }

  // neighbors: Returns a set containing the neighbors of v, i.e. all vertices that can be reached from v along one edge. Since a set is returned, the neighbors are returned in sorted order; use foreach to iterate through the set.
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> S;

    if(AdjacencyList.count(v) == 0) {
      return S;
    }

    return NeighboursList.at(v);
  }

  // getVertices: Returns a vector containing all the vertices currently in the graph.
  vector<VertexT> getVertices() const {
    return this->Vertices;  
  }

  // dump: Dumps the internal state of the graph for debugging purposes.
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    for (int i = 0; i < this->NumVertices(); ++i) {
      output << " " << i << ". " << this->Vertices[i] << endl;
    }
    
    output << endl;
    output << "**Edges:" << endl;
    for (auto &x : AdjacencyList) {
        output << " row " << x.first << ": ";
        for (auto &y : Vertices) {
          if(x.second.count(y) == 1) {
            output << "(T," << x.second.at(y) << ") ";
          } else {
            output << "F ";
          }
        }
        output << endl;
    }

  }
};
