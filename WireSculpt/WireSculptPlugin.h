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

	void compute_all_features(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXi& E, Eigen::MatrixXd& N,
		std::vector<std::vector<int>>& VF,
		std::vector<std::vector<int>>& VFi, Eigen::MatrixXi& IF, Eigen::MatrixXi& OV,
		Eigen::MatrixXd& FC,
		Eigen::MatrixXd& FN, Eigen::MatrixXd& DA, Eigen::MatrixXd& D, Eigen::MatrixXd& L,
		Eigen::MatrixXd& G, Eigen::MatrixXd& dblA);
	void compute_laplacian(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXi& E, Eigen::MatrixXd& G,
		Eigen::MatrixXd& N, Eigen::MatrixXd& L, Eigen::MatrixXd& vertex_is_concave, double beta,
		double eps, double sigma, int clip_bound, int lap_weighting, double filter1_thresh);
};



