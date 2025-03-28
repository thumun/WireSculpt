#include "Edge.h"

Edge::Edge(const Vertex* v1, const Vertex* v2) {
	this->endpoints = { v1, v2 };
	this->length = (v2->mPosition - v1->mPosition).length();
	this->warpedLength = length; // temporary
}

Edge::~Edge() {}