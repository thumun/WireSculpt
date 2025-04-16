#include "Edge.h"

int Edge::lastId = 0;

Edge::Edge(const Vertex* v1, const Vertex* v2): id(lastId) {
	this->endpoints = std::make_pair(v1, v2);
	this->length = (v2->mPosition - v1->mPosition).length();
	this->featureLength = length; // temporary
	this->warpedLength = length; // temporary
	lastId++;
}

Edge::~Edge() {}

float Edge::getLength() {
	return this->length;
}