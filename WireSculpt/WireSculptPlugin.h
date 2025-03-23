#pragma once
#include <string>
#include "Vertex.h"
#include "Edge.h"
#include <vector>

class WireSculptPlugin
{
public:
	WireSculptPlugin() {};

	bool ProcessFile(std::string filePath);
	std::vector<Vertex>* GetVerticies();
	void GetExtremePoints(const std::string& filePath);
	
	
	// Path finding procedure:
	// 1. Run FindTspPath on landmark vertices, using regular distance calculation for nearest neighbors
	// 2. Run FindPath for each pair of consecutive vertices in the path
	// Note: We can also use FindPath in place of the tsp distance calculation, but this may increase computation heavily
	std::vector<Vertex*> FindPath(std::vector<Vertex>& verticies, Vertex* start, Vertex* goal, int vertexCount);
	std::vector<Vertex*> FindTspPath(std::vector<Vertex>& landmarks, Vertex* start);


private:
	std::vector<Vertex> verticies;
	std::vector<Edge> edges; 

	std::vector<std::string> SplitString(const std::string& input, char delimiter);
	bool GetFileExtension(const std::string& filePath);

	// Tsp Nearest Neighbor helper functions
	float findNearestNeighbor(int indexOfVertex, std::vector<Vertex>& landmarks);
	int pickUnvisitedCity(std::vector<int> used);
	int findMinTriangularDistanceEdge(int newVertex, std::vector<int> tour, std::vector<Vertex>& landmarks);
};



