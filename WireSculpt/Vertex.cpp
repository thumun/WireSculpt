#include "Vertex.h"

Vertex::Vertex(const MPoint& position, bool isLandmark) : 
	mPosition(position) {
	this->isLandmark.first = isLandmark;
	//this->neighbors = {};
}

Vertex::~Vertex() {}

void Vertex::setNeighbor(Vertex* v, const Edge* e) {
	this->neighbors.insert({ v, e });
	v->neighbors.insert({this, e});
}