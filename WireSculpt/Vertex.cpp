#include "Vertex.h"
#include "Edge.h"
#include <algorithm>
#include <iostream>
#include <maya/MGlobal.h>

Vertex::Vertex(const MPoint& position, int idNum, bool landmark) : 
	mPosition(position), id(idNum) {
	this->isLandmark.first = landmark;
	resetFGH();
	this->neighbors = {};
}

Vertex::~Vertex() {}

void Vertex::setNeighbor(Vertex* v, Edge* e) {
	this->neighbors.insert(std::make_pair( v, e ));
}

void Vertex::resetFGH() {
	f = 999999;
	g = 999999;
	h = 999999;
}

