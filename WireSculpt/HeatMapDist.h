#pragma once

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>

#include "WireSculptPlugin.h"

#include <unordered_map>
#include <unordered_set>


class HeatMapDist
{
public: 
	HeatMapDist(WireSculptPlugin ws);
    ~HeatMapDist() {};

    void heatDiffusion(int sInput);
    void computePhi(int sInput, WireSculptPlugin& ws);
    void colorScheme(WireSculptPlugin ws, char c);

    std::vector<double> colors;
    double MaxColor;
    int s;
    double t;
    Eigen::MatrixXd A;
    Eigen::MatrixXd Lc;
    Eigen::SparseMatrix<double> M;
    Eigen::FullPivLU<Eigen::MatrixXd> lu;
    Eigen::VectorXd L;
    Eigen::VectorXd phi;
    Eigen::FullPivLU<Eigen::MatrixXd> lulc;

private:
    void computeA(WireSculptPlugin& ws);
    void computeLc(WireSculptPlugin& ws);
    void computeM(WireSculptPlugin& ws);

    void computeFaceArea(std::vector<Face>* faces, std::unordered_map<Face*, double>* mp);
    double computeVertexArea(Vertex* vert, std::unordered_map<Face*, double>* mp, std::vector<Face>* faces);

    std::unordered_set<Vertex*> getNeighbor(const Vertex& vertex);
    double computeCotan(const Vertex * v1, const Vertex * v2, std::vector<Face>* faces); 

    Eigen::VectorXd computeB(int s, WireSculptPlugin& ws);
    double computeTime(WireSculptPlugin& ws);

    Eigen::Vector3d gradientFace(const Face& f, const std::unordered_map<Vertex, double> vu);
    double computeDeltaXu(Vertex* u, std::unordered_map<Face, Eigen::Vector3d> fv, WireSculptPlugin& ws);
    double computeDeltaXuFace(Vertex* curr, Vertex* v1, Vertex* v2);

    double computeAngle(Vertex * u, Vertex * v);
};

