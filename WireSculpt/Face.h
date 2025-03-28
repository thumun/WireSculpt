#pragma once
#include "Vertex.h"

class Vertex; 

class Face
{
public: 
	Face(Vertex * v1, Vertex * v2, Vertex * v3, std::vector<float> norm);
	~Face() {};

	// make private later
	Vertex * v1;
	Vertex * v2;
	Vertex * v3;
	std::vector<float> normal; 

private: 
	

};

