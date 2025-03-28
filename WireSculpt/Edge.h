#pragma once
#include "Vertex.h"
#include <set>

class Edge {
public:
	Edge(const Vertex* v1, const Vertex* v2);
	~Edge();

	// Two Vertex endpoints
	std::set<const Vertex*> endpoints;

	// Warped edge length between two vertices
	float warpedLength;
	int id;
	
protected:
	// Original edge length between two vertices
	float length;
};