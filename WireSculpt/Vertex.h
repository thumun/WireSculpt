#pragma once
//#include "Edge.h"
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <map>

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
	Vertex(const MPoint& position, bool isLandmark = false);
	~Vertex();

	MPoint mPosition;
	MVector mNormal;
	//bool isLandmark;
	std::pair<bool, NonLandmarkData> isLandmark; 
	std::map<const Vertex*, const Edge*> neighbors;

	// Graph Traversal
	float f;
	float g;
	float h;

	void setLandMarkStatus(bool landMark);
	void setNonLandMarkData(float attract, float repel, float dist); 
	void setNeighbor(Vertex *, const Edge *);

	bool getLandMarkStatus(); 


protected:
};