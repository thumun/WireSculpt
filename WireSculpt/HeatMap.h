#pragma once
#include <vector>
#include <string>

// Implement Heat Map method
class HeatMap {
public:
	HeatMap() {};

	bool OpenObjFile(std::string filePath);


private:
	std::vector<std::vector<float>> vertexCoords;
	std::vector<std::vector<int>> faceData;
};