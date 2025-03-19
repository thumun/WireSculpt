#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <maya/MPoint.h>
#include <maya/MVector.h>

#include "WireSculptPlugin.h"

#include <igl/read_triangle_mesh.h>
#include <Eigen/Core>
#include <igl/is_vertex_manifold.h>

// compute laplacian files
#include <igl/per_vertex_normals.h>
#include <igl/edges.h>
#include <igl/gaussian_curvature.h>
#include <igl/adjacency_list.h>


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

// everything below is code from au et al 

// compute all features code 
void compute_incident_faces(Eigen::MatrixXi& E, Eigen::MatrixXi& F, Eigen::MatrixXi& IF) {
    std::map<std::set<int>, std::vector<int>> edge_to_faces_map;

    int not_assigned = F.rows();
    IF = Eigen::MatrixXi::Ones(E.rows(), 2);
    for (int i = 0; i < F.rows(); i++) {
        /// Register this face in all three edges if not already done
        edge_to_faces_map[{F(i, 0), F(i, 1)}].push_back(i);
        edge_to_faces_map[{F(i, 1), F(i, 2)}].push_back(i);
        edge_to_faces_map[{F(i, 2), F(i, 0)}].push_back(i);
    }

    for (int i = 0; i < E.rows(); i++) {
        IF(i, 0) = edge_to_faces_map[{E(i, 0), E(i, 1)}][0];
        IF(i, 1) = edge_to_faces_map[{E(i, 0), E(i, 1)}][1];
    }
}

void compute_opposite_vertices(Eigen::MatrixXi& E, Eigen::MatrixXi& F, Eigen::MatrixXi& OV) {
    std::map<std::set<int>, std::vector<int>> edge_to_opposite_vertices_map;

    int not_assigned = F.rows();
    OV = Eigen::MatrixXi::Ones(E.rows(), 2);
    for (int i = 0; i < F.rows(); i++) {
        /// Register this face in all three edges if not already done
        edge_to_opposite_vertices_map[{F(i, 0), F(i, 1)}].push_back(F(i, 2));
        edge_to_opposite_vertices_map[{F(i, 1), F(i, 2)}].push_back(F(i, 0));
        edge_to_opposite_vertices_map[{F(i, 2), F(i, 0)}].push_back(F(i, 1));
    }

    for (int i = 0; i < E.rows(); i++) {
        OV(i, 0) = edge_to_opposite_vertices_map[{E(i, 0), E(i, 1)}][0];
        OV(i, 1) = edge_to_opposite_vertices_map[{E(i, 0), E(i, 1)}][1];
    }
}

void compute_dihedral_angle(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXi& E, Eigen::MatrixXi& IF,
    Eigen::MatrixXi& OV, Eigen::MatrixXd& FN, Eigen::MatrixXd& DA) {
    DA.resize(E.rows(), 1);
    std::map<std::set<int>, int> vertices_to_face_index_map;
    for (int i = 0; i < F.rows(); i++) {
        vertices_to_face_index_map[{F(i, 0), F(i, 1), F(i, 2)}] = i;
    }

    for (int i = 0; i < E.rows(); i++) {
        Eigen::Vector3d p0 = V.row(E(i, 0));
        Eigen::Vector3d p1 = V.row(E(i, 1));
        Eigen::Vector3d p2 = V.row(OV(i, 0));
        Eigen::Vector3d p3 = V.row(OV(i, 1));
        int f1 = vertices_to_face_index_map[{E(i, 0), E(i, 1), OV(i, 0)}];
        int f2 = vertices_to_face_index_map[{E(i, 0), E(i, 1), OV(i, 1)}];
        Eigen::Vector3d n1 = FN.row(f1);
        Eigen::Vector3d n2 = FN.row(f2);
        double normal_angle = atan2((n1.cross(n2)).norm(), n1.dot(n2));

        if ((p3 - p2).dot(n1 - n2) < 0)
            normal_angle = -normal_angle;

        DA(i) = normal_angle;
    }
}

void compute_distance(Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXd& D) {
    D.resize(E.rows(), 1);
    for (int i = 0; i < E.rows(); i++) {
        D(i) = sqrt((V(E(i, 0), 0) - V(E(i, 1), 0)) * (V(E(i, 0), 0) - V(E(i, 1), 0)) +
            (V(E(i, 0), 1) - V(E(i, 1), 1)) * (V(E(i, 0), 1) - V(E(i, 1), 1)) +
            (V(E(i, 0), 2) - V(E(i, 1), 2)) * (V(E(i, 0), 2) - V(E(i, 1), 2)));
    }
}

void WireSculptPlugin::compute_all_features(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXi& E, Eigen::MatrixXd& N,
    std::vector<std::vector<int>>& VF,
    std::vector<std::vector<int>>& VFi, Eigen::MatrixXi& IF, Eigen::MatrixXi& OV,
    Eigen::MatrixXd& FC,
    Eigen::MatrixXd& FN, Eigen::MatrixXd& DA, Eigen::MatrixXd& D, Eigen::MatrixXd& L,
    Eigen::MatrixXd& G, Eigen::MatrixXd& dblA) {

    igl::edges(F, E);

    igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_DEFAULT, N);

    igl::doublearea(V, F, dblA);

    Eigen::MatrixXd PD1, PD2, PV1, PV2;

    compute_incident_faces(E, F, IF);
    compute_opposite_vertices(E, F, OV);

    igl::per_face_normals(V, F, FN);

    compute_dihedral_angle(V, F, E, IF, OV, FN, DA);
    compute_distance(V, E, D);

    igl::gaussian_curvature(V, F, G);

    igl::per_face_normals(V, F, FN);

    igl::per_vertex_normals(V, F, FN, N);
}

// compute laplacian

bool is_concave(Eigen::MatrixXd& V, Eigen::MatrixXd& N, std::vector<std::vector<int>>& A, int i) {
    bool concave = false;
    for (int j = 0; j < A[i].size(); j++) {
        if ((((V.row(i) - V.row(A[i][j])) / ((V.row(i) - V.row(A[i][j])).norm())).dot(N.row(A[i][j]) - N.row(i))) >
            0.01) {
            concave = true;
        }
    }
    return concave;
}

double getAngle3D(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const bool in_degree) {
    // Compute the actual angle
    double rad = v1.normalized().dot(v2.normalized());
    if (rad < -1.0)
        rad = -1.0;
    else if (rad > 1.0)
        rad = 1.0;
    return (in_degree ? acos(rad) * 180.0 / M_PI : acos(rad));
}

bool connIsConvex(const Eigen::Vector3d& source_centroid, const Eigen::Vector3d& target_centroid,
    const Eigen::Vector3d& source_normal, const Eigen::Vector3d& target_normal, double& normal_angle) {
    bool is_convex = true;
    bool is_smooth = true;

    normal_angle = getAngle3D(source_normal, target_normal, true);
    //  Geometric comparisons
    Eigen::Vector3d vec_t_to_s, vec_s_to_t;

    vec_t_to_s = source_centroid - target_centroid;
    vec_s_to_t = -vec_t_to_s;

    Eigen::Vector3d ncross;
    ncross = source_normal.cross(target_normal);

    // vec_t_to_s is the reference direction for angle measurements
    // Convexity Criterion: Check if connection of patches is convex. If this is the case the two supervoxels should be merged.
    if ((getAngle3D(vec_t_to_s, source_normal, true) - getAngle3D(vec_t_to_s, target_normal, true)) <= 0) {
        normal_angle = -normal_angle;
        is_convex &= true;  // connection convex
    }
    else {
        is_convex &= (normal_angle < 0.0);  // concave connections will be accepted  if difference of normals is small
    }
    return (is_convex && is_smooth);
}

std::set<int> two_ring_neighbors(int id, std::vector<std::vector<int>>& A) {
    std::set<int> neighbors;
    for (int j = 0; j < A[id].size(); j++) {
        int neighbor = A[id][j];
        for (int k = 0; k < A[neighbor].size(); k++) {
            neighbors.insert(A[neighbor][k]);
        }
    }
    return neighbors;
}

bool is_connection_concave(Eigen::MatrixXd& V, Eigen::MatrixXd& N, int i, int j) {
    bool concave = false;
    if ((((V.row(i) - V.row(j)) / ((V.row(i) - V.row(j)).norm())).dot(N.row(j) - N.row(i))) > 0.01) {
        concave = true;
    }
    return concave;
}

bool filter_1b(Eigen::MatrixXd& V, Eigen::MatrixXd& N, int i, std::vector<std::vector<int>>& A, double ratio) {
    std::set<int> neighbors = two_ring_neighbors(i, A);
    neighbors.erase(i);
    double concave_neighbors = 0.0;
    for (int vertex_id : neighbors) {
        if (is_connection_concave(V, N, vertex_id, i)) {
            concave_neighbors += 1.0;
        }
    }
    return (concave_neighbors / neighbors.size() > ratio);
}

void covariance_matrix(Eigen::MatrixXd& mat, Eigen::MatrixXd& cov) {
    Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
    cov = Eigen::MatrixXd::Zero(mat.cols(), mat.cols());
    cov = (centered.adjoint() * centered) / double(mat.rows() - 1);
}

void eigenvalues(Eigen::MatrixXd& mat, Eigen::MatrixXd& ev) {
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(mat.cols(), mat.cols());
    covariance_matrix(mat, cov);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(cov);
    ev = eigensolver.eigenvalues();
}

bool filter_2(Eigen::MatrixXd& V, int i, std::vector<std::vector<int>>& A) {
    //Eigen::MatrixXd one_ring = Eigen::MatrixXd::Zero(A[i].size(), 3);
    Eigen::MatrixXd one_ring = Eigen::MatrixXd::Zero(A[i].size() + 1, 3);
    for (int l = 0; l < A[i].size(); l++) {
        int j = A[i][l];
        one_ring.row(l) = V.row(j);
    }
    one_ring.row(A[i].size()) = V.row(i);
    Eigen::MatrixXd ev;
    eigenvalues(one_ring, ev);
    return (ev(0, 0) / ev.sum() > 0.001);
}

void WireSculptPlugin::compute_laplacian(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXi& E, Eigen::MatrixXd& G,
    Eigen::MatrixXd& N, Eigen::MatrixXd& L, Eigen::MatrixXd& vertex_is_concave, double beta,
    double eps, double sigma, int clip_bound, int lap_weighting, double filter1_thresh) {

    Eigen::MatrixXd E2 = Eigen::MatrixXd::Zero(E.rows(), 1);
    for (int i = 0; i < E.rows(); i++) {
        E2(i, 0) = (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
    }
    double mean_dist = E2.mean();
    Eigen::MatrixXd::Index maxRow, maxCol;
    E2.maxCoeff(&maxRow, &maxCol);
    double max_dist = E2(maxRow, 0);

    std::vector<std::vector<int>> VF;
    std::vector<std::vector<int>> VI;
    igl::vertex_triangle_adjacency(F.maxCoeff() + 1, F, VF, VI);
    Eigen::MatrixXd area;
    igl::doublearea(V, F, area);

    std::vector<double> connections;
    for (int i = 0; i < E.rows(); i++) {
        connections.push_back((abs(G(E(i, 0))) + abs(G(E(i, 1), 0))));
    }

    double max_val = *std::max_element(connections.begin(), connections.end());
    G = G.array() / max_val * 1.0;

    std::vector<std::vector<int>> A;
    igl::adjacency_list(F, A);
    Eigen::MatrixXd V_is_concave = Eigen::MatrixXd::Zero(V.rows(), 1);

    for (int i = 0; i < V.rows(); i++) {
        if (is_concave(V, N, A, i) && filter_1b(V, N, i, A, 0.3) && filter_2(V, i, A)) {
            V_is_concave(i, 0) = 1.0;
        }
    }

    vertex_is_concave = V_is_concave;
    Eigen::MatrixXd w = Eigen::MatrixXd::Zero(V.rows(), V.rows());

    for (int i = 0; i < V.rows(); i++) {
        double sum_w_ik = 0.0;
        for (int l = 0; l < A[i].size(); l++) {
            int j = A[i][l];
            double normal_angle;
            connIsConvex(V.row(i), V.row(j), N.row(i), N.row(j), normal_angle);
            ///Convex
            if (V_is_concave(i, 0) == 1.0 || V_is_concave(j, 0) == 1.0) {
                w(i, j) = (1 - ((V.row(i) - V.row(j)).norm() / max_dist)) * (1 - (abs(G(i, 0)) + abs(G(j, 0)))) * beta +
                    eps;
            }
            else {
                w(i, j) = (1 - ((V.row(i) - V.row(j)).norm() / max_dist)) * (1 - (abs(G(i, 0)) + abs(G(j, 0)))) + eps;
            }
            sum_w_ik += w(i, j);
        }


        for (int l = 0; l < A[i].size(); l++) {
            int j = A[i][l];
            w(i, j) = w(i, j) / sum_w_ik;
        }
    }

    for (int i = 0; i < w.rows(); i++) {
        w(i, i) = -1.0;
    }

    L = w;
}

// extreme pts calc 

/*
std::vector<int> get_extreme_points(Eigen::MatrixXi& F, Eigen::MatrixXd& V, Eigen::MatrixXd& L, int index_given,
    Eigen::MatrixXi& E) {

    std::vector<int> extreme_points;
    std::set<int> extreme_point_set;
    std::ranlux48 gen;
    std::uniform_int_distribution<int> uniform_0_255(0, V.rows() - 1);
    int prev_previous_num_points = 0;
    int previous_num_points = 0;
    int current_num_points = 0;
    int max_iters = 2;
    int iters = 0;

    std::vector<std::vector<int>> A;
    igl::adjacency_list(F, A);

    std::vector<T> basic_tripletList;
    basic_tripletList.reserve(E.rows() * 2 + V.rows() + 20);
    for (int i = 0; i < L.rows(); i++) {
        for (int j = 0; j < L.cols(); j++) {
            if (L(i, j) != 0.0)
                basic_tripletList.push_back(T(i, j, L(i, j)));
        }
    }

    //std::cout << ">> Finding extreme points...\n";
    do {

        iters++;
        int index = uniform_0_255(gen);
        std::cout << ">> Starting with random point " << index << "\n";

        int q = get_starting_point_fast(F, V, L, index, A, E, basic_tripletList);

        index = q;

        std::cout << ">> Selecting new point " << index << "\n";

        Eigen::SparseMatrix<double> LC_sparse(V.rows() + A[index].size() + 1, V.rows());

        std::vector<T> tripletList = basic_tripletList;

        for (int i = 0; i < A[index].size(); i++) {
            tripletList.push_back(T(V.rows() + i, A[index][i], 1.0));
        }
        tripletList.push_back(T(V.rows() + A[index].size(), index, 1.0));
        LC_sparse.setFromTriplets(tripletList.begin(), tripletList.end());

        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(A[index].size() + 1, V.rows());
        for (int i = 0; i < A[index].size(); i++) {
            C(i, A[index][i]) = 1.0;
        }
        C(A[index].size(), index) = 1.0;
        Eigen::VectorXd B = Eigen::VectorXd::Ones(A[index].size() + 1);
        B(A[index].size(), 0) = 0.0;
        Eigen::VectorXd nulls = Eigen::VectorXd::Zero(V.rows());

        Eigen::VectorXd Mb(L.rows() + C.rows());
        Mb << nulls,
            B;

        Eigen::SPQR<Eigen::SparseMatrix<double>> solver;
        solver.compute(LC_sparse);

        Eigen::MatrixXd x = solver.solve(Mb);

        Eigen::MatrixXd::Index maxRow, maxCol;
        x.maxCoeff(&maxRow, &maxCol);

        Eigen::MatrixXd::Index minRow, minCol;
        x.minCoeff(&minRow, &minCol);

        //std::cout << ">> Field varies from " << x(minRow,0) << " to " << x(maxRow,0) << std::endl;

        std::vector<int> candidate_extreme_points;

        std::vector<std::pair<double, int>> extreme_point_queue;

        for (int i = 0; i < V.rows(); i++) {
            double max_val = std::numeric_limits<double>::lowest();
            double min_val = std::numeric_limits<double>::max();
            for (int adj_vert : A[i]) {
                max_val = std::max(x(adj_vert, 0), max_val);
                min_val = std::min(x(adj_vert, 0), min_val);
            }
            if (x(i) >= max_val || x(i) <= min_val) {
                //std::cout << ">> Found extreme value: " << i << "\n";
                if (extreme_point_set.size() < 2) {
                    candidate_extreme_points.push_back(i);
                }
                else {
                    double max_dist = max_geodesic_dist(extreme_point_set, V);
                    if (!in_proximity_to(V.row(i), extreme_point_set, V, max_dist * 0.15)) {
                        double neighbors = 0.0;
                        for (int j = 0; j < A[i].size(); j++) {
                            neighbors = neighbors + x(A[i][j]);
                        }
                        neighbors = abs(x(i) - (neighbors / A[i].size()));
                        //std::cout << "Neighbor_values: " << neighbors << "\n";
                        extreme_point_queue.push_back(std::pair<double, int>(neighbors, i));
                    }
                }
            }
        }

        if (extreme_point_queue.size() > 0 && extreme_point_set.size() >= 2) {
            std::sort(extreme_point_queue.begin(), extreme_point_queue.end(),
                boost::bind(&std::pair<double, int>::first, _1) <
                boost::bind(&std::pair<double, int>::first, _2));

            for (int i = extreme_point_queue.size() - 1; i >= 0; i--) {
                double max_dist = max_geodesic_dist(extreme_point_set, V);
                if (!in_proximity_to(V.row(extreme_point_queue[i].second), extreme_point_set, V, max_dist * 0.15)) {
                    extreme_point_set.insert(extreme_point_queue[i].second);
                    //std::cout << ">> Added extreme point " << extreme_point_queue[i].second << std::endl;
                }
                else {
                    //std::cout << ">> Rejected extreme point " << extreme_point_queue[i].second << std::endl;
                }

            }
        }

        if (candidate_extreme_points.size() > 0) {
            std::vector<double> pair_dist;
            std::vector<int> candidate1;
            std::vector<int> candidate2;
            int selected_1 = 0;
            int selected_2 = 0;
            /// Add pair, that is the farthest apart first
            for (int i = 0; i < candidate_extreme_points.size(); i++) {
                for (int j = 0; j < i; j++) {
                    pair_dist.push_back(
                        (V.row(candidate_extreme_points[i]) - V.row(candidate_extreme_points[j])).norm());
                    candidate1.push_back(candidate_extreme_points[i]);
                    candidate2.push_back(candidate_extreme_points[j]);
                }
            }
            double max_dist = 0.0;
            for (int i = 0; i < pair_dist.size(); i++) {
                if (max_dist < pair_dist[i]) {
                    max_dist = pair_dist[i];
                    selected_1 = candidate1[i];
                    selected_2 = candidate2[i];
                }
            }
            extreme_point_set.insert(selected_1);
            extreme_point_set.insert(selected_2);
        }

    } while (iters < max_iters);

    extreme_points.clear();
    for (auto val : extreme_point_set) {
        //std::cout << ">> Final Extreme Point: " << val << std::endl;
        extreme_points.push_back(val);
    }
    return extreme_points;
}
*/

#if EXEDEBUG
int main() {
    WireSculptPlugin ws = WireSculptPlugin();
    //ws.ProcessFile("D:/CGGT/AdvTopics/WireSculpt/testobj/cow.obj");
    ws.GetExtremePoints("D:/CGGT/AdvTopics/WireSculpt/testobj/cube.obj");
}
#endif // EXEDEBUG