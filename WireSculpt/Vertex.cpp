#include "Vertex.h"

Vertex::Vertex(const MPoint& position, const MVector& normal, bool isLandmark) : 
	mPosition(position), mNormal(normal) {
	this->isLandmark = isLandmark;
	//this->neighbors = {};
}

Vertex::~Vertex() {}