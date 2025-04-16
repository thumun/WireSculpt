#pragma once
#include "Vertex.h"
#include <set>

class Edge {
public:
	Edge(const Vertex* v1, const Vertex* v2);
	~Edge();

	// Two Vertex endpoints
	std::pair<const Vertex*, const Vertex*> endpoints;

	static int lastId;

	// Warped edge length between two vertices
	float featureLength;	// length after warping by feature attraction
	float warpedLength;		// accumulated warping (applying path repulsion)
	int id;

	float getLength(); 
	
protected:
	// Original edge length between two vertices
	float length;
};