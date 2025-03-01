#pragma once
#include "Vertex.h"

// Represents the landmark sites (extreme points from Au et al. 2011)
class Landmark : public Vertex {
public:
	Landmark(const MPoint& position, const MVector& normal, bool isLandmark = true, const Vertex* specialEndpt = nullptr);
	~Landmark();

	// Reference to the vertex endpoint of the special edge it represents
	const Vertex* specialEndpt;
};