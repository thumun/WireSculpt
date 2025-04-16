#pragma once
//#include "Edge.h"
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <map>
#include <vector>
#include <queue>
#include <cmath>

class Edge;

//struct NonLandmarkData {
//	// Feature attraction weight w, and path repulsion weight w*
//	float wAttract;
//	float wRepel;
//
//	// Geodesic distance from feature lines
//	float geoDistance;
//};

// Represents the vertices of the input 3D model 
class Vertex {
public:
	typedef std::pair<Vertex*, Edge*> Ve;

	Vertex(const MPoint& position, bool landmark = false);
	//Vertex(const Vertex& v); 
	~Vertex();

	MPoint mPosition;
	//MVector mNormal;
	bool isLandmark;
	//std::pair<bool, NonLandmarkData> isLandmark; 
	std::map<Vertex*, Edge*> neighbors;
	int id;

	static int lastId;

	// Feature attraction weight w, and path repulsion weight w*
	float wAttract;
	float wRepel;

	// Geodesic distance from feature lines
	float geoDistance;

	// A* path traversal
	float f;	// f = g + h
	float g;	// movement cost from starting vertex to given vertex
	float h;	// estimated movement cost from given vertex to target vertex

	// Overload comparison operators for priority queue
	bool operator==(const Vertex& other) const;

	void resetFGH();

	void setLandMarkStatus(bool landMark);
	void setNonLandMarkData(float attract, float repel, float dist); 
	void setNeighbor(Vertex *v, Edge *e);

	bool getLandMarkStatus(); 


protected:

};

struct VertexPtrCompare
{
	bool operator()(const Vertex* lhs, const Vertex* rhs) const {
		return lhs->f > rhs->f; // min-heap: lower f is higher priority
	}
};