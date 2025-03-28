#include "Face.h"

#include <vector>
#include "Vertex.h"

Face::Face(Vertex * v1, Vertex * v2, Vertex * v3, std::vector<float> norm): 
	v1(v1), v2(v2), v3(v3), normal(norm) {}