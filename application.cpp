// application.cpp 
// Rayyan Athar
//
//
// Adam T Koehler, PhD
// University of Illinois Chicago
// CS 251, Fall 2023
//
// Project Original Variartion By:
// Joe Hummel, PhD
// University of Illinois at Chicago
//
// 
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include <stack>
#include <limits>
#include "tinyxml2.h"
#include "dist.h"
#include "graph.h"
#include "osm.h"


using namespace std;
using namespace tinyxml2;

// simply copy and paste this line as a global variable
const double INF = numeric_limits<double>::max();

//
// Implement your standard application here
//
 
class prioritize  
{
public:
  bool operator()(const pair<long long,double>& p1, const pair<long long,double>& p2) const
  {
    return p1.second > p2.second; 
  }
};

// Function checks if a vertex is visited or not and returns true if it is visited or false otherwise
bool isVisited(vector<long long> visited, long long currentV) {
  int size = visited.size();
	for (int i = 0; i < size; i++) {
  	if (visited[i] == currentV)
      return true;
  }
  return false;
}

// Function checks if a string of a building exists in the vector of buildings
BuildingInfo searchBuilding(string query, vector<BuildingInfo>& Buildings) {
  // An empty building object is created, and two for loops are created to loop through the building vectors returning the building the matches the name or abbreviation, otherwise if not found the empty building is returned
  BuildingInfo result;
  result.Fullname = "";
  result.Abbrev = "";

  for(auto& building : Buildings) {
    if(building.Abbrev == query) {
      return building;
    }
  }
  for(auto& building : Buildings) { 
    if (building.Fullname.find(query) != string::npos) {
      return building;
    }
  }
  return result;
}

// Function returns the coordinates of the nearest node to the specified building 
Coordinates nearestNode(vector<FootwayInfo> Footways, map<long long, Coordinates> Nodes, BuildingInfo building) {
  Coordinates nearest;
  double min = INF;

  // Two nested for loops are created and the distance between a node of a footway and the building is calculated, and if it is less that the minimum distance, the minimum distance and the coordinates are updated accordingly
  for(auto& footway : Footways) {
    for(int i = 0; i < footway.Nodes.size(); i++) {
      double dist = distBetween2Points(Nodes.at(footway.Nodes.at(i)).Lat, Nodes.at(footway.Nodes.at(i)).Lon, building.Coords.Lat, building.Coords.Lon);
      if(dist < min) {
        min = dist;
        nearest.Lat = Nodes.at(footway.Nodes.at(i)).Lat;
        nearest.Lon = Nodes.at(footway.Nodes.at(i)).Lon;
        nearest.ID = Nodes.at(footway.Nodes.at(i)).ID;
      }
    }    
  }
  return nearest;
}

// Function performs dijkstra's algorithm on a given point and calculates the shortest distance along with all the predecessors
void dijkstraAlgorithm(graph<long long, double> &G, long long start, map<long long, double> &distances, map<long long, long long> &pred) {
  priority_queue<pair<long long, double>, vector<pair<long long, double>>, prioritize> unvisitedQueue;
  vector<long long> visited;

  // For loop sets all distances to infinity and predecessors to 0, pushes all vertices to the priority queue with infinite distances
  for (auto &currV : G.getVertices()) {
    distances[currV] = INF;
    pred[currV] = 0;
    unvisitedQueue.push(make_pair(currV, INF));
  }

  distances[start] = 0;
  unvisitedQueue.push(make_pair(start, 0));

  // Main loop: Continue until all vertices are visited or the priority queue is empty
  while (!unvisitedQueue.empty()) {
    // Extract the vertex with the smallest tentative distance from the priority queue
    long long currentV = unvisitedQueue.top().first;
    unvisitedQueue.pop();

    // If the distance to the current vertex is still infinity, break the loop (no more reachable vertices)
    if (distances[currentV] == INF) {
      break;
    } 
    // If the current vertex is already visited, skip to the next iteration
    else if (isVisited(visited, currentV)) {
      continue;
    } 
    else {
      visited.push_back(currentV);
    }

    // For loop loops through the neighbors of the current vertex
    for (auto &adjV : G.neighbors(currentV)) {
      // Get the weight of the edge between the current vertex and its neighbor
      double edgeWeight = 0.0;
      G.getWeight(currentV, adjV, edgeWeight);

      // Calculate the alternative distance to the neighbor through the current vertex
      double alternativePathDistance = distances[currentV] + edgeWeight;

      // If the alternative distance is shorter than the known distance, update the distance and predecessor
      if (alternativePathDistance < distances[adjV]) {
        distances[adjV] = alternativePathDistance;
        pred[adjV] = currentV;
        // Push the neighbor to the priority queue with the updated distance
        unvisitedQueue.push(make_pair(adjV, alternativePathDistance));
      }
    }
  }
}


// Function gets the path from the start to end of the predecessor and outputs the path
void getPath(map<long long, long long> &pred, long long start, long long end) {
  long long currV = pred[end];
  stack<long long> path;
  path.push(end);

  while(currV != 0){
    path.push(currV);
    currV = pred[currV];
  }

  cout << path.top();
  path.pop();
  while(!path.empty()) {
    cout << " -> " << path.top();
    path.pop();
  }

  cout << endl << endl;
}

void application(
    map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings, graph<long long, double>& G) {
  string person1Building, person2Building;

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);

    // Searches Buildings 1 and 2 and if not found continues to ask input and search again
    BuildingInfo building1 = searchBuilding(person1Building, Buildings);
    BuildingInfo building2 = searchBuilding(person2Building, Buildings);

    if(building1.Fullname == ""){
      cout << "Person 1's building not found" << endl;
      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }
    
    if(building2.Fullname == ""){
      cout << "Person 2's building not found" << endl;
      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }

    cout << "\nPerson 1's point:\n " << building1.Fullname << "\n (" << building1.Coords.Lat << ", " << building1.Coords.Lon << ")" << endl;
    cout << "Person 2's point:\n " << building2.Fullname << "\n ("  << building2.Coords.Lat << ", " << building2.Coords.Lon << ")" << endl;

    // Calculates the midpoint between both building1 and building2
    BuildingInfo buildingCenter;
    Coordinates midpoint = centerBetween2Points(building1.Coords.Lat, building1.Coords.Lon, building2.Coords.Lat, building2.Coords.Lon);
    
    // A for loop loops through all the buildings and caluclates the distance between them and the midpoint, and if the ditance is less than the minimum, the minimum and the center building are updated accordingly
    double min = INF;
    for(auto& building: Buildings) {
      double dist = distBetween2Points(building.Coords.Lat, building.Coords.Lon, midpoint.Lat, midpoint.Lon);
      if(dist < min) {
        min = dist;
        buildingCenter = building;
      }
    }

    cout << "Destination Building:\n " << buildingCenter.Fullname << "\n ("  << buildingCenter.Coords.Lat << ", " << buildingCenter.Coords.Lon << ")" << endl;
    cout << endl;

    // The nearest nodes to building1, building2 and the center building are calculated
    Coordinates nearestNodeBuilding1 = nearestNode(Footways, Nodes, building1);
    Coordinates nearestNodeBuilding2 = nearestNode(Footways, Nodes, building2);
    Coordinates nearestNodeBuildingCenter = nearestNode(Footways, Nodes, buildingCenter);


    cout << "Nearest P1 node:\n " << nearestNodeBuilding1.ID << "\n (";
    cout << nearestNodeBuilding1.Lat << ", " << nearestNodeBuilding1.Lon << ")" << endl;

    cout << "Nearest P2 node:\n " << nearestNodeBuilding2.ID << "\n (";
    cout << nearestNodeBuilding2.Lat << ", " << nearestNodeBuilding2.Lon << ")" << endl;

    cout << "Nearest destination node:\n " << nearestNodeBuildingCenter.ID << "\n (";
    cout << nearestNodeBuildingCenter.Lat << ", " << nearestNodeBuildingCenter.Lon << ")" << endl;


    map<long long, double> distances1;
    map<long long, long long> pred1;

    map<long long, double> distances2;
    map<long long, long long> pred2;

    long long start1 = nearestNodeBuilding1.ID;
    long long start2 = nearestNodeBuilding2.ID;
    long long start3 = nearestNodeBuildingCenter.ID;


    // Dijkstra's algorithm is called on the first and second nearest nodes to building1 and building2 respectively, to calculate the shortest path between them and the center
    dijkstraAlgorithm(G, start1, distances1, pred1);
    dijkstraAlgorithm(G, start2, distances2, pred2);


    // If Dijkstra's algorithm failed to calculate a shortest path between either the first two starting points and the center than the user inputs both buildings and resets the whole process
    if(distances1[start3] == INF) {
      cout << "Sorry, destination unreachable." << endl;
      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }
    else if (distances2[start3] == INF) {
      cout << "Sorry, destination unreachable." << endl;
      cout << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }
    // Otherwise the shortest distance and path between the first two starting points and the center are displayed
    else {
      cout << "Person 1's distance to dest: " << distances1[start3];
      cout << " miles" << endl;
      cout << "Path: ";
      getPath(pred1, start1, start3);

      cout << "Person 2's distance to dest: " << distances2[start3];
      cout << " miles" << endl;
      cout << "Path: ";
      getPath(pred2, start2, start3);
    }

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }    
}

int main() {
  // graph<long long, double> G;

  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  graph<long long, double> G;

  // A for loop adds the nodes to the graph
  for(auto& node : Nodes) {
    G.addVertex(node.first);
  }
  

  // Two nested loops add the nodes to the graph
  for(auto& footway : Footways) {
    for(int i = 0; i < footway.Nodes.size()-1; i++) {
      double dist = distBetween2Points(Nodes.at(footway.Nodes.at(i)).Lat, Nodes.at(footway.Nodes.at(i)).Lon, Nodes.at(footway.Nodes.at(i+1)).Lat, Nodes.at(footway.Nodes.at(i+1)).Lon);
      G.addEdge(footway.Nodes.at(i), footway.Nodes.at(i+1), dist);
      G.addEdge(footway.Nodes.at(i+1), footway.Nodes.at(i), dist);
    }
  }
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;


  cout << endl;

  // Execute Application
  application(Nodes, Footways, Buildings, G);

  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
