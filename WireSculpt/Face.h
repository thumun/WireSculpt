#pragma once
#include "Vertex.h"

class Vertex; 

class Face
{
public: 
	Face(Vertex * v1, Vertex * v2, Vertex * v3, std::vector<float> norm);
	~Face() {};

	static int lastId;

	// make private later
	/*Vertex * v1;
	Vertex * v2;
	Vertex * v3;*/
	std::vector<Vertex*> verticies; 
	std::vector<float> normal; 
	int id;


private: 
	

};

