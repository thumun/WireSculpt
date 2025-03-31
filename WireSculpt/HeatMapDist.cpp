#include "HeatMapDist.h"

#include <unordered_map>
#include <unordered_set>

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
    double area = 0.0;
    // Implement your vertex area calculation here.
    // Calculate the sum of one-third of the areas of incident faces.

    for (auto f : *faces) {
        double area = 0.5 * (f.v1->mPosition.x * (f.v2->mPosition.y - f.v3->mPosition.y) +
                            f.v2->mPosition.x * (f.v3->mPosition.y - f.v1->mPosition.y) +
                            f.v3->mPosition.x * (f.v1->mPosition.y - f.v2->mPosition.y));

        mp->insert({&f, area});
    }
}

double HeatMapDist::computeVertexArea(Vertex* vert, std::unordered_map<Face*, double> * mp, std::vector<Face> * faces) {
    // face area divided by 3!
    // find face with this vertex 

    double area = 0.0;
    for (auto f : *faces) {
        if (&f.v1 == &vert || &f.v2 == &vert || &f.v3 == &vert) {
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
    Eigen::Vector3d a, b;
    Vertex* v3 = nullptr;
    bool found = false;

    for (auto f : *faces) {

        if ((*f.v1 == *v1 && *f.v2 == *v2) || (*f.v1 == *v2 && *f.v2 == *v1) 
            || (*f.v2 == *v1 && *f.v3 == *v2) || (*f.v3 == *v1 && *f.v1 == *v2) 
            || (*f.v3 == *v1 && *f.v1 == *v2) || (*f.v1 == *v1 && *f.v3 == *v2)) {

            v3 = (f.v1 != v1 && f.v1 != v2) ? f.v1 :
                (f.v2 != v1 && f.v2 != v2) ? f.v2 : f.v3;

            found = true;
            break;
        }

        //if ((*f.v1 == *v1) && (*f.v2 == *v2)) {
        //    
        //    //auto v3 = f.v3; 
        //    
        //    Eigen::Vector3d p1 = Eigen::Vector3d(v1->mPosition.x, v1->mPosition.y, v1->mPosition.z);
        //    Eigen::Vector3d p2 = Eigen::Vector3d(v2->mPosition.x, v2->mPosition.y, v2->mPosition.z);
        //    Eigen::Vector3d p3 = Eigen::Vector3d(v3->mPosition.x, v3->mPosition.y, v3->mPosition.z);

        //    a = p1 - p3;
        //    b = p2 - p3;

        //    found = true;
        //    break;
        //}
    }

    // error val
    if (!found) {
        return 0.0;
    }
    else {

        // cos = a . b / magnitude 
        // cot = cos / sqrt(1 - cos ^ 2)

        /*a.normalize();
        b.normalize();

        auto cos = a.dot(b);
        auto cot = cos / sqrt(1 - cos * cos);

        return cot; */

        Eigen::Vector3d p1 = Eigen::Vector3d(v1->mPosition.x, v1->mPosition.y, v1->mPosition.z);
        Eigen::Vector3d p2 = Eigen::Vector3d(v2->mPosition.x, v2->mPosition.y, v2->mPosition.z);
        Eigen::Vector3d p3 = Eigen::Vector3d(v3->mPosition.x, v3->mPosition.y, v3->mPosition.z);

        a = p1 - p3;
        b = p2 - p3;

        double dotProduct = a.dot(b);
        double crossProductMagnitude = a.cross(b).norm();

        if (crossProductMagnitude < 1e-8) { // Handle degenerate triangles
            return 0.0;
        }

        return dotProduct / crossProductMagnitude;
    }
}

void HeatMapDist::computeM(WireSculptPlugin& ws) {
    for (int i = 0; i < A.rows(); ++i) {
        std::unordered_set<Vertex*> s = getNeighbor(ws.verticies[i]);
        double tmp = A(i, i) - t * Lc.coeff(i, i); // Use coeff() for sparse matrix access

        this->M.insert(i, i) = tmp;

        //int indx = 0; 

        for (Vertex* v : s) {
            tmp = -t * Lc.coeff(i, v->id); // Use coeff() for sparse matrix access
            if (std::abs(tmp) != 0.0) {
                this->M.insert(i, v->id) = tmp;
            }
            //indx++; 
        }
    }

    this->M.makeCompressed();
    //this->M = M; 
}

void HeatMapDist::heatDiffusion(int sInput) {
    s = sInput;
    Eigen::VectorXd K = Eigen::VectorXd::Zero(A.rows());

    K(s) = 1.0f; 

    this->L = this->lu.solve(K);
}

void HeatMapDist::computePhi(int sInput) {
    s = sInput;
    Eigen::VectorXd b = computeB(s);

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

Eigen::VectorXd HeatMapDist::computeB(int s) {
    Eigen::VectorXd b(A.rows());
    for (int i = 0; i < b.rows(); i++) {
        b(i) = 0.0;
    }
    return b;
}


void HeatMapDist::colorScheme(WireSculptPlugin ws, char c) {
    Eigen::VectorXd cs;

    if (c == 'h') {
        cs = L;
    }
    else if (c == 'd') {
        cs = phi;
    }
    else {
        return;
    }

    std::unordered_map<Vertex*, float> lv; 
    std::vector<Vertex> vertices = ws.verticies;
    
    for (int i = 0; i < vertices.size(); i++) {
        lv.insert({&vertices[i], cs(i)});
    }
}