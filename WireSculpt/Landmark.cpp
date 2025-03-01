#include "Landmark.h"

Landmark::Landmark(const MPoint& position, const MVector& normal, bool isLandmark, 
	const Vertex* specialEndpt) : Vertex(position, normal, isLandmark) {
	this->specialEndpt = specialEndpt;
}

Landmark::~Landmark() {}