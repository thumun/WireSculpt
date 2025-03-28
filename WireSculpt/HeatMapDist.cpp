#include "HeatMapDist.h"

#include <unordered_map>
#include <unordered_set>

HeatMapDist::HeatMapDist(WireSculptPlugin ws) {
    colors.resize(ws.faces.size());
    computeA(ws);
    computeLc(ws);
    computeM(ws);

    lulc.compute(Lc);
    //lu.compute(M);
}

void HeatMapDist::computeA(WireSculptPlugin& ws) {
    std::vector<Vertex>& vertices = ws.verticies;
    Eigen::SparseMatrix<double> M(vertices.size(), vertices.size()); // Create M |V| x |V| matrix
    std::unordered_map<Face*, double> mp;
    computeFaceArea(&ws.faces, &mp);

    for (int i = 0; i < vertices.size(); ++i) {
        M.insert(i, i) = computeVertexArea(&vertices[i], &mp, &ws.faces);
    }

    M.makeCompressed();

    this->A = M;
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
    //Eigen::SparseMatrix<double> M(vertices.size(), vertices.size()); // Create M |V|x|V| matrix

    for (size_t i = 0; i < vertices.size(); ++i) {
        std::unordered_set<Vertex*> s = getNeighbor(vertices[i]); // Get the neighbor of the vertex vi
        double sum = 0.0;

        int indx = 0; 

        for (Vertex* v : s) {
            double cotan = computeCotan(vertices[i], *v, &ws.faces); // Compute cotan formula
            Lc.insert(i, indx) = -1.0 * cotan;
            sum += cotan;
            indx++; 
        }

        Lc.insert(i, i) = sum;
    }

    Lc.makeCompressed();

    //this->Lc = M; 
}

std::unordered_set<Vertex*> HeatMapDist::getNeighbor(const Vertex& vertex) {
    std::unordered_set<Vertex*> neighbors;

    for (auto n: vertex.neighbors) {
        neighbors.insert(n.first);
    }

    return neighbors;
}

double HeatMapDist::computeCotan(const Vertex& v1, const Vertex& v2, std::vector<Face>* faces) {
    Eigen::Vector3d a, b;
    bool found = false;

    for (auto f : *faces) {

        if (f.v1 == &v1 && f.v2 == &v2) {
            
            auto v3 = f.v3; 
            
            Eigen::Vector3d p1 = Eigen::Vector3d(v1.mPosition.x, v1.mPosition.y, v1.mPosition.z);
            Eigen::Vector3d p2 = Eigen::Vector3d(v2.mPosition.x, v2.mPosition.y, v2.mPosition.z);
            Eigen::Vector3d p3 = Eigen::Vector3d(v3->mPosition.x, v3->mPosition.y, v3->mPosition.z);

            a = p1 - p3;
            b = p2 - p3;

            found = true;
            break;
                
        }
    }

    if (!found) {
        return 0.0; //or handle the error in an appropriate way.
    }

    double dotProduct = a.dot(b);
    double crossProductMagnitude = a.cross(b).norm();

    // Handle the case where the cross product magnitude is zero (degenerate triangle).
    if (crossProductMagnitude < 1e-8) {
        return 0.0; // Or handle it in a way appropriate for your application.
    }

    return dotProduct / crossProductMagnitude;    
}

void HeatMapDist::computeM(WireSculptPlugin& ws) {
    for (int i = 0; i < A.rows(); ++i) {
        std::unordered_set<Vertex*> s = getNeighbor(ws.verticies[i]);
        double tmp = A(i, i) - t * Lc.coeff(i, i); // Use coeff() for sparse matrix access

        M.insert(i, i) = tmp;

        int indx = 0; 

        for (Vertex* v : s) {
            tmp = -t * Lc.coeff(i, indx); // Use coeff() for sparse matrix access
            if (std::abs(tmp) != 0.0) {
                M.insert(i, indx) = tmp;
            }
            indx++; 
        }
    }

    M.makeCompressed();
    this->M = M; 
}

void HeatMapDist::heatDiffusion(int sInput) {
    s = sInput;
    Eigen::VectorXd K = Eigen::VectorXd::Zero(A.rows());

    K(s) = 1.0f; 

    //L = lu.solve(K);
}

void HeatMapDist::computePhi(int sInput) {
    s = sInput;
    Eigen::VectorXd b = computeB(s);

    //phi = lulc.solve(b);
    
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