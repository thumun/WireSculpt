#pragma once

#include <vector>
#include <Eigen/Core>

void compute_all_features(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXi& E, Eigen::MatrixXd& N,
	std::vector<std::vector<int>>& VF,
	std::vector<std::vector<int>>& VFi, Eigen::MatrixXi& IF, Eigen::MatrixXi& OV,
	Eigen::MatrixXd& FC,
	Eigen::MatrixXd& FN, Eigen::MatrixXd& DA, Eigen::MatrixXd& D, Eigen::MatrixXd& L,
	Eigen::MatrixXd& G, Eigen::MatrixXd& dblA);
void compute_laplacian(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXi& E, Eigen::MatrixXd& G,
	Eigen::MatrixXd& N, Eigen::MatrixXd& L, Eigen::MatrixXd& vertex_is_concave, double beta,
	double eps, double sigma, int clip_bound, int lap_weighting, double filter1_thresh);
std::vector<int> get_extreme_points(Eigen::MatrixXi& F, Eigen::MatrixXd& V, Eigen::MatrixXd& L, int index_given,
	Eigen::MatrixXi& E, double proximity, double maxVal);
