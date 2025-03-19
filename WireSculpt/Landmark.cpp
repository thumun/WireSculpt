#include "Landmark.h"

Landmark::Landmark(const MPoint& position, const MVector& normal, bool isLandmark, 
	const Vertex* specialEndpt) : Vertex(position, isLandmark) {
	this->specialEndpt = specialEndpt;
}

Landmark::~Landmark() {}