#include "HeatMapDist.h"

#include <unordered_map>
#include <unordered_set>
#include <maya/MGlobal.h>

HeatMapDist::HeatMapDist(WireSculptPlugin ws) {
    colors.resize(ws.faces.size());
    t = 100 * computeTime(ws);
    computeA(ws);
    computeLc(ws);
    computeM(ws);

    lulc.compute(Lc);
    lu.compute(M);
}

double HeatMapDist::computeTime(WireSculptPlugin& ws) {
    double sum = 0.0; 

    for (auto e : ws.edges) {
        sum += e.getLength();
    }

    return (sum/ws.edges.size()) * (sum / ws.edges.size());
}

void HeatMapDist::computeA(WireSculptPlugin& ws) {
    std::vector<Vertex>& vertices = ws.verticies;
    Eigen::SparseMatrix<double> ACalc(vertices.size(), vertices.size()); // Create M |V| x |V| matrix
    std::unordered_map<Face*, double> mp;
    computeFaceArea(&ws.faces, &mp);

    for (int i = 0; i < vertices.size(); ++i) {
        ACalc.insert(i, i) = computeVertexArea(&vertices[i], &mp, &ws.faces);
    }

    ACalc.makeCompressed();

    this->A = ACalc;
}

void HeatMapDist::computeFaceArea(std::vector<Face> * faces, std::unordered_map<Face*, double>* mp) {
    // Implement your vertex area calculation here.
    // Calculate the sum of one-third of the areas of incident faces.

    for (auto f : *faces) {

        Eigen::Vector3d p0 = {f.verticies[0]->mPosition.x, f.verticies[0]->mPosition.y, f.verticies[0]->mPosition.z};
        Eigen::Vector3d p1 = {f.verticies[1]->mPosition.x, f.verticies[1]->mPosition.y, f.verticies[1]->mPosition.z};
        Eigen::Vector3d p2 = { f.verticies[2]->mPosition.x, f.verticies[2]->mPosition.y, f.verticies[2]->mPosition.z};

        Eigen::Vector3d edge1 = p1 - p0;
        Eigen::Vector3d edge2 = p2 - p0;

        double area = 0.5 * edge1.cross(edge2).norm();

        /*double area = 0.5 * (f.verticies[0]->mPosition.x * (f.verticies[1]->mPosition.y - f.verticies[2]->mPosition.y) +
                            f.verticies[1]->mPosition.x * (f.verticies[2]->mPosition.y - f.verticies[0]->mPosition.y) +
                            f.verticies[2]->mPosition.x * (f.verticies[0]->mPosition.y - f.verticies[1]->mPosition.y));*/

        mp->insert({&f, area});
    }
}

double HeatMapDist::computeVertexArea(Vertex* vert, std::unordered_map<Face*, double> * mp, std::vector<Face> * faces) {
    // face area divided by 3!
    // find face with this vertex 

    double area = 0.0;
    for (auto f : *faces) {
        if (&f.verticies[0] == &vert || &f.verticies[1] == &vert || &f.verticies[2] == &vert) {
            area += (mp->at(&f) / 3.0f);
        }
    }

    return area;
}

void HeatMapDist::computeLc(WireSculptPlugin& ws) {
    const std::vector<Vertex>& vertices = ws.verticies;
    Lc.resize(vertices.size(), vertices.size());
    Eigen::SparseMatrix<double> LcCalc(vertices.size(), vertices.size()); // Create M |V|x|V| matrix

    for (size_t i = 0; i < vertices.size(); ++i) {
        std::unordered_set<Vertex*> neighbors = getNeighbor(vertices[i]); // Get the neighbor of the vertex vi
        double sum = 0.0;

        for (Vertex* v : neighbors) {
            double cotan = computeCotan(&vertices[i], v, &ws.faces); // Compute cotan formula
            LcCalc.coeffRef(i, v->id) = -1.0 * cotan;
            sum += cotan;
        }

        LcCalc.insert(i, i) = sum;
    }

    LcCalc.makeCompressed();
    this->Lc = Eigen::MatrixXd(LcCalc);
}

std::unordered_set<Vertex*> HeatMapDist::getNeighbor(const Vertex& vertex) {
    std::unordered_set<Vertex*> neighbors;

    for (auto n: vertex.neighbors) {
        neighbors.insert(n.first);
    }

    return neighbors;
}

double HeatMapDist::computeCotan(const Vertex * v1, const Vertex * v2, std::vector<Face>* faces) {
    Eigen::Vector3d p_i = {v1->mPosition.x, v1->mPosition.y, v1->mPosition.z};
    Eigen::Vector3d p_j = { v2->mPosition.x, v2->mPosition.y, v2->mPosition.z };
    
    // Find the two adjacent faces (if they exist)
    Vertex* v_alpha = nullptr; // Third vertex in the first adjacent face
    Vertex* v_beta = nullptr;  // Third vertex in the second adjacent face
    
    for (const Face& f : *faces) {
        bool has_vi = false, has_vj = false;
        Vertex* third_vertex = nullptr;
    
        for (auto v : f.verticies) {
            if (v->id == v1->id) {
                has_vi = true;
            }
            else if (v->id == v2->id) {
                has_vj = true;
            }
            else {
                third_vertex = v;
            }
        }
    
        if (has_vi && has_vj) {
            if (!v_alpha) {
                v_alpha = third_vertex;
            }
            else if (!v_beta) {
                v_beta = third_vertex;
            }
        }
    }
    
    double cot_alpha = 0.0, cot_beta = 0.0;
    
    // Compute cotangent for the first adjacent face (if it exists)
    if (v_alpha) {
        Eigen::Vector3d p_alpha = { v_alpha->mPosition.x, v_alpha->mPosition.y, v_alpha->mPosition.z };
        Eigen::Vector3d e1 = p_i - p_alpha;
        Eigen::Vector3d e2 = p_j - p_alpha;
        double dot = e1.dot(e2);
        double cross_norm = e1.cross(e2).norm();
        if (cross_norm > 1e-8) { // Avoid division by zero
            cot_alpha = dot / cross_norm;
        }
    }
    
    // Compute cotangent for the second adjacent face (if it exists)
    if (v_beta) {
        Eigen::Vector3d p_beta = { v_beta->mPosition.x, v_beta->mPosition.y, v_beta->mPosition.z };
        Eigen::Vector3d e1 = p_i - p_beta;
        Eigen::Vector3d e2 = p_j - p_beta;
        double dot = e1.dot(e2);
        double cross_norm = e1.cross(e2).norm();
        if (cross_norm > 1e-8) { // Avoid division by zero
            cot_beta = dot / cross_norm;
        }
    }
    
    // Return the average cotangent weight
    if (!v_beta) {
        return cot_alpha;
    } 
    else {
        return 0.5 * (cot_alpha + cot_beta);
    }

}

void HeatMapDist::computeM(WireSculptPlugin& ws) {
    Eigen::SparseMatrix<double> MCalc(ws.verticies.size(), ws.verticies.size()); // Create M |V| x |V| matrix

    for (int i = 0; i < A.rows(); ++i) {

        std::unordered_set<Vertex*> s = getNeighbor(ws.verticies[i]);
        double tmp = A(i, i) - t * Lc.coeff(i, i); // Use coeff() for sparse matrix access

        MCalc.insert(i, i) = tmp;

        //int indx = 0; 

        for (Vertex* v : s) {
            tmp = -t * Lc.coeff(i, v->id); // Use coeff() for sparse matrix access
            if (std::abs(tmp) != 0.0) {
                MCalc.insert(i, v->id) = tmp;
            }
            //indx++; 
        }
    }

    MCalc.makeCompressed();
    this->M = MCalc;
}

void HeatMapDist::heatDiffusion(int sInput) {
    s = sInput;
    Eigen::VectorXd K = Eigen::VectorXd::Zero(A.rows());

    K(s) = 1.0f; 

    this->L = this->lu.solve(K);
}

void HeatMapDist::heatDiffusionFromPath(const std::vector<Vertex*>& sourceVertexIndices) {
    Eigen::VectorXd K = Eigen::VectorXd::Zero(A.rows());

    double heatPerVertex = 1.0 / sourceVertexIndices.size();

    for (auto v : sourceVertexIndices) {
        K(v->id) = heatPerVertex;
    }

    this->L = this->lu.solve(K);
}

void HeatMapDist::computePhi(int sInput, WireSculptPlugin& ws) {
    s = sInput;
    Eigen::VectorXd b = computeB(s, ws);

    this->phi = this->lulc.solve(b);
    
    double mini = std::numeric_limits<double>::infinity();
    
    for (int i = 0; i < phi.size(); ++i) {
        mini = std::min(phi(i), mini);
    }

    if (mini < 0) {
        for (int i = 0; i < phi.size(); ++i) {
            phi(i) = phi(i) + std::abs(mini);
        }
    }
    
}

Eigen::VectorXd HeatMapDist::computeB(int s, WireSculptPlugin& ws) {
    std::unordered_map<int, double> vu;
    std::vector<Vertex> vertices = ws.verticies;

    heatDiffusion(s); // compute the heat diffusion

    // store in a map the heat associated with each vector
    for (size_t i = 0; i < vertices.size(); ++i) {
        vu[vertices[i].id] = L(i);
    }

    // compute the gradient vector for each face and store it in a map.
    std::vector<Face> faces = ws.faces;
    std::unordered_map<int, Eigen::Vector3d> fv;

    for (size_t i = 0; i < faces.size(); ++i) {
        fv[faces[i].id] = gradientFace(faces[i], vu);
    }

    // compute delta X
    Eigen::VectorXd X(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        X(i) = -1.0 * computeDeltaXu(&vertices[i], fv, ws);
    }
    return X;
}

Eigen::Vector3d HeatMapDist::gradientFace(const Face& f, const std::unordered_map<int, double> vu) {
    Eigen::Vector3d p0 = {f.verticies[0]->mPosition.x, f.verticies[0]->mPosition.y, f.verticies[0]->mPosition.z};
    Eigen::Vector3d p1 = {f.verticies[1]->mPosition.x, f.verticies[1]->mPosition.y, f.verticies[1]->mPosition.z};
    Eigen::Vector3d p2 = { f.verticies[2]->mPosition.x, f.verticies[2]->mPosition.y, f.verticies[2]->mPosition.z};

    auto a = p1 - p0;
    auto b = p2 - p0;

    auto N = a.cross(b);

    double Af = 0.5 * N.norm();

    if (Af < 1e-8) {
        return Eigen::Vector3d::Zero();
    }

    double u0 = vu.at(f.verticies[0]->id);
    double u1 = vu.at(f.verticies[1]->id);
    double u2 = vu.at(f.verticies[2]->id);

    Eigen::Vector3d gradient =
        (u0 * (p2 - p1) +
            u1 * (p0 - p2) +
            u2 * (p1 - p0)) / (2.0 * Af);

    return gradient;
}

double HeatMapDist::computeDeltaXu(Vertex* u, std::unordered_map<int, Eigen::Vector3d> fv, WireSculptPlugin& ws) {
    // get each f on v & divide by num faces 
    
    double delta = -1.0f; 

    bool v1 = false;  
    bool v2 = false; 
    bool v3 = false; 

    for (auto f : ws.faces) {
        if (f.verticies[0]->id == u->id) {
            v1 = true; 
            delta += computeDeltaXuFace(u, f.verticies[1], f.verticies[2], fv[f.id]);
        }
        else if (f.verticies[1]->id == u->id) {
            v2 = true;
            delta += computeDeltaXuFace(u, f.verticies[0], f.verticies[2], fv[f.id]);
        }
        else if (f.verticies[2]->id == u->id) {
            v3 = true;
            delta += computeDeltaXuFace(u, f.verticies[0], f.verticies[1], fv[f.id]);
        }
    }

    return 0.5 * delta;
}

double HeatMapDist::computeDeltaXuFace(Vertex* curr, Vertex* v1, Vertex* v2, Eigen::Vector3d Xj) {
    /*double theta2 = computeAngle(v2, v1);
    double theta1 = computeAngle(curr, v2);*/
    double theta1 = computeAngle(curr, v2);
    double theta2 = computeAngle(curr, v1);

    Eigen::Vector3d E2 = { v2->mPosition.x - curr->mPosition.x, v2->mPosition.y - curr->mPosition.y, v2->mPosition.z - curr->mPosition.z };
    Eigen::Vector3d E1 = { v1->mPosition.x - curr->mPosition.x, v1->mPosition.y - curr->mPosition.y, v1->mPosition.z - curr->mPosition.z };

    double inner1 = E1.dot(Xj);
    double inner2 = E2.dot(Xj);

    double cot1 = 1 / std::tan(theta1);
    double cot2 = 1 / std::tan(theta2);

    return (cot1 * inner1 + cot2 * inner2);
    
}

double HeatMapDist::computeAngle(Vertex * v1, Vertex * v2) {

    Eigen::Vector3d u = { v1->mPosition[0], v1->mPosition[1], v1->mPosition[2] };
    Eigen::Vector3d v = { v2->mPosition[0], v2->mPosition[1], v2->mPosition[2] };

    Eigen::Vector3d e1 = u - v;

    // For typical usage, you'll want to pass a second adjacent vertex
    // to compute the angle between two edges
    Eigen::Vector3d e2 = v - u; 

    double dot = e1.dot(e2);
    double norm_e1 = e1.norm();
    double norm_e2 = e2.norm();

    // Clamp to avoid numerical issues with acos()
    double cos = dot / (norm_e1 * norm_e2);
    //cos = std::clamp(cos, -1.0, 1.0);
    if (cos > 1.0) {
        cos = 1.0; 
    }
    else if (cos < -1.0) {
        cos = -1.0;
    }

    /*double innerProduct = u.dot(v); 

    double lu = std::sqrt(u.dot(u));
    double lv = std::sqrt(v.dot(v));

    double cos = innerProduct / (lu * lv);*/

    return std::acos(cos);
}

std::unordered_map<Vertex*, float> HeatMapDist::colorScheme(WireSculptPlugin& ws, char c) {
    Eigen::VectorXd cs;

    if (c == 'h') {
        cs = L;
    }
    else if (c == 'd') {
        cs = phi;
    }

    std::unordered_map<Vertex*, float> lv; 
    std::vector<Vertex>* vertices = ws.GetVerticies();
    
    for (int i = 0; i < vertices->size(); i++) {
        Vertex* v = &((*vertices)[i]);
        lv.insert({ v, cs(i) });
    }

    return lv;
}