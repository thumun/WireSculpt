#include "NonLandmark.h"

NonLandmark::NonLandmark(const MPoint& position, const MVector& normal, bool isLandmark, 
	float wAttract, float wRepel, float geoDistance) : Vertex(position, isLandmark) {
	this->wAttract = wAttract;
	this->wRepel = wRepel;
	this->geoDistance = geoDistance;
}

NonLandmark::~NonLandmark() {}