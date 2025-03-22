#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <maya/MPoint.h>
#include <maya/MVector.h>

#include "WireSculptPlugin.h"
#include "ExtremePoints.h"

#include <igl/read_triangle_mesh.h>
#include <Eigen/Core>
#include <igl/is_vertex_manifold.h>


using namespace std;

// checks to make sure file is obj 
bool WireSculptPlugin::GetFileExtension(const string& filePath) {
    size_t dotPosition = filePath.find_last_of(".");
    return (filePath.substr(dotPosition + 1) == "obj");
}

// split string helper 
vector<string> WireSculptPlugin::SplitString(const string& input, char delimiter) {
    vector<string> tokens;
    size_t start = 0;
    size_t end = input.find(delimiter);

    while (end != string::npos) {
        tokens.push_back(input.substr(start, end - start));
        start = end + 1;
        end = input.find(delimiter, start);
    }

    tokens.push_back(input.substr(start));

    return tokens;
}

// processes the obj file 
// returns false: issue with file (contents or type) 
// returns true: success 
bool WireSculptPlugin::ProcessFile(std::string filePath) {

    // incorrect obj type 
    if (!GetFileExtension(filePath)) {
        return false; 
    }

    // failure to open 
    ifstream fin;
    fin.open(filePath);
    if (!fin.is_open()) {
        return false;
    }

    string line;

    // temp storage 
    vector<MVector> pos; 
    vector<MVector> vertexNormals;

    while (getline(fin, line)) {

        vector<string> fields = SplitString(line, ' ');

        if (line == "-1") {
            break;
        }
        // adding vert positions to temp vector
        else if (fields[0] == "v") {
            pos.push_back(MVector(stof(fields[1]), stof(fields[2]), stof(fields[3])));
            verticies.push_back(MPoint(pos[pos.size()-1]));
        }

        // calculating normals by going through faces 
        else if (fields[0] == "f") {
            // temp storage to calc normals w/ vert pos
            vector<int> faceVerts; 

            // going through verticies that make up a face 
            for (int i = 1; i < fields.size(); i++) {
                vector<string> triplets = SplitString(fields[i], '/');

                faceVerts.push_back(stoi(triplets[0]));
            }

            // draw edges with verts to get cross prod
            for (int i = 0; i < faceVerts.size(); i++) {
                
                int indx = faceVerts[i] - 1;

                MVector current = pos[indx];
                MVector prev = pos[(indx - 1)%pos.size()];
                MVector next = pos[(indx + 1)%pos.size()];

                MVector crossProd = (current - prev) ^ (next-current);
                crossProd.normalize();

                // for each calc can make Vertex obj and add to list
                //verticies.push_back(Vertex(MPoint(pos[0]), crossProd));
            }

            // setting up neighbors 
            for (int i = 0; i < faceVerts.size(); i++) {
                // check if edge already exists before creating 

                auto vertOne = faceVerts[i] - 1; 
                auto vertTwo = faceVerts[(i + 1)%faceVerts.size()] - 1;

                bool edgeFound = false; 
                int edgeIndx = -1; 

                // want to do while but how to loop at same time
                if (edges.size() > 0) {
                    for (int j = 0; j < edges.size(); j++) {

                        if ((edges[j].endpoints.find(&verticies[vertOne]) != edges[j].endpoints.end()) &&
                            (edges[j].endpoints.find(&verticies[vertTwo]) != edges[j].endpoints.end())) {
                            edgeFound = true; 
                            edgeIndx = j; 
                            break;
                        }
                    }
                }

                if (!edgeFound) {
                    edges.push_back(Edge(&verticies[vertOne], &verticies[vertTwo]));
                    edgeIndx = edges.size() - 1;
                }

                verticies[vertOne].setNeighbor(&verticies[vertTwo], &edges[edgeIndx]);
            }
        }
    }

    fin.close();
    return true;
}

void WireSculptPlugin::GetExtremePoints(const std::string& filePath) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Eigen::MatrixXi E;    // Edges
    Eigen::MatrixXd N;    // Normals
    std::vector<std::vector<int>> VF;
    std::vector<std::vector<int>> VFi;
    Eigen::MatrixXi IF;   // Incident Faces
    Eigen::MatrixXi OV;   // Opposite Vertices
    Eigen::MatrixXd FC;   //
    Eigen::MatrixXd FN;   // Face Normals
    Eigen::MatrixXd DA;   // Dihedral Angle
    Eigen::MatrixXd D;    // Distance
    Eigen::MatrixXd L;    // Laplacian
    Eigen::MatrixXd G;    // Gaussian Curvature
    Eigen::MatrixXd dblA; // Area

    double beta = 0.1; 
    double eps = 0.000001;
    double sigma = 0.001;
    int lap_weighting = 0;
    Eigen::MatrixXd vertex_is_concave;
    int index = 0;
    std::vector<int> extreme_points;

    // no idea what this value is 
    int clip_bound = 0;

    igl::read_triangle_mesh(filePath, V, F);

    Eigen::MatrixXi B;
    igl::is_vertex_manifold(F, B);
    // link this to plug in 
    if (B.minCoeff() == 0) {
        std::cerr << ">> The loaded mesh is not manifold.\n";
    }

    // debugged - this seems to be working 
    compute_all_features(V, F, E, N, VF, VFi, IF, OV, FC, FN, DA, D, L, G, dblA);

    compute_laplacian(V, F, E, G, N, L, vertex_is_concave, beta, eps, sigma, clip_bound, lap_weighting, 0.4);

    extreme_points = get_extreme_points(F, V, L, index, E);
}

std::vector<Vertex>* WireSculptPlugin::GetVerticies() {
    return &(this->verticies);
}


#if EXEDEBUG
int main() {
    WireSculptPlugin ws = WireSculptPlugin();
    //ws.ProcessFile("D:/CGGT/AdvTopics/WireSculpt/testobj/cow.obj");
    ws.GetExtremePoints("D:/CGGT/AdvTopics/WireSculpt/testobj/cube.obj");
}
#endif // EXEDEBUG