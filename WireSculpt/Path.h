#pragma once
#include "Vertex.h"
#include <vector>

class Path {
public:
	Path();
	~Path();

	// Vector of vertices that make up path
	std::vector<const Vertex*> vertices;	// potentially unneeded

	// Vector of path segments 
	// where each element is a vector of vertices in a single path segment
	std::vector<std::vector<const Vertex*>> pathSegments;

};