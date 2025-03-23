#include "Edge.h"

Edge::Edge(int idVal, const Vertex* v1, const Vertex* v2) {
	this->id = idVal;
	this->endpoints = { v1, v2 };
	this->length = (v2->mPosition - v1->mPosition).length();
	this->warpedLength = 1;
}

Edge::~Edge() {}