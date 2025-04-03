#include "Face.h"

#include <vector>
#include "Vertex.h"

int Face::lastId = 0;

Face::Face(Vertex * v1, Vertex * v2, Vertex * v3, std::vector<float> norm): normal(norm), id(lastId) {

	verticies.push_back(v1);
	verticies.push_back(v2);
	verticies.push_back(v3);

	lastId++; 
}