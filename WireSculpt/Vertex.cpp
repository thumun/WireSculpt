#include "Vertex.h"
#include "Edge.h"
#include <algorithm>
#include <iostream>

Vertex::Vertex(const MPoint& position, int idNum, bool landmark) : 
	mPosition(position), id(idNum) {
	this->isLandmark.first = landmark;
	resetFGH();
	this->neighbors = {};
}

Vertex::~Vertex() {}

void Vertex::setNeighbor(Vertex* v, Edge* e) {
	this->neighbors.insert({ v, e });
	v->neighbors.insert({this, e});
}

// A* path traversal
bool Vertex::operator>(const Vertex& other) const
{
	return f > other.f;
}

//bool Vertex::operator==(const Vertex& other) const
//{
//	// ignoring normals for now
//	return mPosition[0] == other.mPosition[0] 
//		&& mPosition[1] == other.mPosition[1]
//		&& mPosition[2] == other.mPosition[2];
//}

void Vertex::resetFGH() {
	f = 999999;
	g = 999999;
	h = 999999;
}

