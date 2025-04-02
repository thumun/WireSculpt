#pragma once
#include <string>
#include "Vertex.h"
#include "Edge.h"
#include <vector>

#include <Eigen/Core>

struct vec3f {
	float x, y, z;
	vec3f(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
};

class WireSculptPlugin
{
public:
	WireSculptPlugin() {};

	bool ProcessFile(std::string filePath);
	std::vector<Vertex>* GetVerticies();
	std::vector<int> GetExtremePoints(const std::string& filePath);
	
	
	// Path finding procedure:
	// 1. Run FindTspPath on landmark vertices, using regular distance calculation for nearest neighbors
	// 2. Run FindPath for each pair of consecutive vertices in the path
	// Note: We can also use FindPath in place of the tsp distance calculation, but this may increase computation heavily
	std::vector<Vertex*> FindPath(std::vector<Vertex>& verticies, Vertex* start, Vertex* goal, int vertexCount);
	std::vector<int> FindTspPath(std::vector<Vertex*> landmarks, int start);	// outputs naive tsp path
	std::vector<int> TwoOptTspPath(std::vector<Vertex*> landmarks, int start, int maxIters);	// optimize naive tsp path

	// Suggestive contours
	std::vector<std::pair<vec3f, vec3f>> GetContours(float fovChoice, int viewChoice, int contoursChoice, float testSCChoice, const char* filename);

private:
	std::vector<Vertex> verticies;
	std::vector<Edge> edges; 

	std::vector<std::string> SplitString(const std::string& input, char delimiter);
	bool GetFileExtension(const std::string& filePath);

	// Tsp Nearest Neighbor helper functions
	float computeLMDistance(Vertex& lm1, Vertex& lm2);
	float insertionCost(Vertex& lm1, Vertex& lm2, Vertex& lm3);
	int findNearestNeighbor(int indexOfVertex, std::vector<Vertex*> landmarks, std::vector<int> used);
	int pickUnvisitedCity(std::vector<int> used);
	int findMinTriangularDistanceEdge(int newLM, std::vector<int> tour, std::vector<Vertex*> landmarks);
	std::vector<int> swapEdge(std::vector<int> tour, int i, int j);
	
};



