#pragma once
#include <string>
#include "Vertex.h"
#include <vector>

class WireSculptPlugin
{
public:
	WireSculptPlugin() {};

	bool ProcessFile(std::string filePath);
	std::vector<Vertex>* GetVerticies();

private:
	std::vector<Vertex> verticies;

	std::vector<std::string> SplitString(const std::string& input, char delimiter);
	bool GetFileExtension(const std::string& filePath);
};



