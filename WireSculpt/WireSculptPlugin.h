#pragma once
#include <string>
#include "Vertex.h"
#include "Edge.h"
#include <vector>

#include <Eigen/Core>

class WireSculptPlugin
{
public:
	WireSculptPlugin() {};

	bool ProcessFile(std::string filePath);
	std::vector<Vertex>* GetVerticies();
	void GetExtremePoints(const std::string& filePath);

private:
	std::vector<Vertex> verticies;
	std::vector<Edge> edges; 

	std::vector<std::string> SplitString(const std::string& input, char delimiter);
	bool GetFileExtension(const std::string& filePath);
};



