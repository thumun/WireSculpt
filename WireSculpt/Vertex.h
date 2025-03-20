#pragma once
//#include "Edge.h"
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <map>
#include <vector>
#include <queue>
#include <cmath>

class Edge;

struct NonLandmarkData {
	// Feature attraction weight w, and path repulsion weight w*
	float wAttract;
	float wRepel;

	// Geodesic distance from feature lines
	float geoDistance;
};

// Represents the vertices of the input 3D model 
class Vertex {
public:
	Vertex(const MPoint& position, int idNum, bool landmark = false);
	~Vertex();

	MPoint mPosition;
	//MVector mNormal;
	//bool isLandmark;
	std::pair<bool, NonLandmarkData> isLandmark; 
	std::map<Vertex*, Edge*> neighbors;
	int id;

	// A* path traversal
	float f;	// f = g + h
	float g;	// movement cost from starting vertex to given vertex
	float h;	// estimated movement cost from given vertex to target vertex

	// Overload comparison operators for priority queue
	bool operator>(const Vertex& other) const;
	bool operator==(const Vertex& other) const;

	void resetFGH();

	void setLandMarkStatus(bool landMark);
	void setNonLandMarkData(float attract, float repel, float dist); 
	void setNeighbor(Vertex *, Edge *);

	bool getLandMarkStatus(); 


protected:
};