#include "Vertex.h"
#include "Edge.h"
#include <algorithm>
#include <iostream>
#include <maya/MGlobal.h>

int Vertex::lastId = 0;

Vertex::Vertex(const MPoint& position, bool landmark) : 
	mPosition(position), id(lastId) {

	//this->isLandmark.first = landmark;
	this->isLandmark = landmark;
	this->wAttract = 0;
	this->wRepel = 0;
	resetFGH();
	this->neighbors = {};
	lastId++;
}

//Vertex::Vertex(const Vertex& v) : mPosition(v.mPosition), 
//	id(lastId++), isLandmark(v.isLandmark), neighbors(v.neighbors) {
//
//	resetFGH();
//}

Vertex::~Vertex() {}

void Vertex::setNeighbor(Vertex* v, Edge* e) {
	this->neighbors.insert(std::make_pair( v, e ));
}

bool Vertex::operator==(const Vertex& other) const
{
	return (other.id == this->id);
}

void Vertex::resetFGH() {
	f = 999999;
	g = 999999;
	h = 999999;
}

