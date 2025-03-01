#pragma once
//#include "Edge.h"
#include <maya/MPoint.h>
#include <maya/MVector.h>
#include <map>

class Edge;

// Represents the vertices of the input 3D model 
class Vertex {
public:
	Vertex(const MPoint& position, const MVector& normal, bool isLandmark = false);
	~Vertex();

	MPoint mPosition;
	MVector mNormal;
	bool isLandmark;
	std::map<const Vertex*, const Edge*> neighbors;

protected:
};