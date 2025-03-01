#pragma once
#include "Vertex.h"

// Represents the non-landmark sites (not selected as extreme points from Au et al. 2011)
class NonLandmark : public Vertex {
public:
	NonLandmark(const MPoint& position, const MVector& normal, bool isLandmark = false, float wAttract = 0, float wRepel = 0, float geoDistance = 0);
	~NonLandmark();

	// Feature attraction weight w, and path repulsion weight w*
	float wAttract;
	float wRepel;

	// Geodesic distance from feature lines
	float geoDistance;
};