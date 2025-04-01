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
        double area = 0.5 * (f.verticies[0]->mPosition.x * (f.verticies[1]->mPosition.y - f.verticies[2]->mPosition.y) +
                            f.verticies[1]->mPosition.x * (f.verticies[2]->mPosition.y - f.verticies[0]->mPosition.y) +
                            f.verticies[2]->mPosition.x * (f.verticies[0]->mPosition.y - f.verticies[1]->mPosition.y));

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
    Eigen::Vector3d a, b;
    Vertex* v3 = nullptr;
    bool found = false;

    //std::cout << "v1 input:" << v1->id;
    //std::cout << "\tv2 input:" << v2->id << std::endl;

    for (auto f : *faces) {

        //std::cout << "\nface verts:" << f.verticies[0]->id << ", " << f.verticies[1]->id << ", " << f.verticies[2]->id << std::endl; 

        bool foundV1 = false; 
        bool foundV2 = false; 

        for (auto v : f.verticies) {
            if (*v == *v1) {
                foundV1 = true; 
            }
            else if (*v == *v2) {
                foundV2 = true; 
            }
            else {
                v3 = v; 
            }
        }

        if (foundV1 && foundV2) {
            found = true; 
            break;
        }

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
    /*Eigen::VectorXd b(A.rows());
    for (int i = 0; i < b.rows(); i++) {
        b(i) = 0.0;
    }
    return b;*/

    std::unordered_map<Vertex, double> vu;
    std::vector<Vertex> vertices = ws.verticies;

    heatDiffusion(s); // compute the heat diffusion

    // store in a map the heat associated with each vector
    for (size_t i = 0; i < vertices.size(); ++i) {
        vu[vertices[i]] = L(i);
    }

    // compute the gradient vector for each face and store it in a map.
    std::vector<Face> faces = ws.faces;
    std::unordered_map<Face, Eigen::Vector3d> fv;

    for (size_t i = 0; i < faces.size(); ++i) {
        Eigen::Vector3d x = gradientFace(faces[i], vu);
        // normalize x
        auto sum = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        x[0] /= sum;
        x[1] /= sum;
        x[2] /= sum;
        //x = Utilities::normalize(x);
        fv[faces[i]] = x;
    }

    // compute delta X
    Eigen::VectorXd X(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        X(i) = -1.0 * computeDeltaXu(vertices[i], fv);
    }
    return X;
}

Eigen::Vector3d HeatMapDist::gradientFace(const Face& f, const std::unordered_map<Vertex, double> vu) {
    double Af = 0.5 * (f.verticies[0]->mPosition.x * (f.verticies[1]->mPosition.y - f.verticies[2]->mPosition.y) +
        f.verticies[1]->mPosition.x * (f.verticies[2]->mPosition.y - f.verticies[0]->mPosition.y) +
        f.verticies[2]->mPosition.x * (f.verticies[0]->mPosition.y - f.verticies[1]->mPosition.y));

    auto v1 = f.verticies[0];
    auto v2 = f.verticies[1];
    auto v3 = f.verticies[2];

    Eigen::Vector3d p1 = Eigen::Vector3d(v1->mPosition.x, v1->mPosition.y, v1->mPosition.z);
    Eigen::Vector3d p2 = Eigen::Vector3d(v2->mPosition.x, v2->mPosition.y, v2->mPosition.z);
    Eigen::Vector3d p3 = Eigen::Vector3d(v3->mPosition.x, v3->mPosition.y, v3->mPosition.z);

    auto a = p1 - p3;
    auto b = p2 - p3;

    double N = a.cross(b);

    double u = vu.at(*v2);
    Eigen::Vector3d c1 = N.crossProduct(v1);
    c1 = c1 * u;

    u = vu.at(*v3);
    Eigen::Vector3d c2 = N.crossProduct(v2);
    c2 = c2 * u;

    u = vu.at(*v1);
    Eigen::Vector3d c3 = N.crossProduct(v3);
    c3 = c3 * u;

    double delta_x = c1.x + c2.x + c3.x;
    double delta_y = c1.y + c2.y + c3.y;
    double delta_z = c1.z + c2.z + c3.z;

    double area = 1.0 / (2.0 * Af);
    delta_x = delta_x * area;
    delta_y = delta_y * area;
    delta_z = delta_z * area;

    return Eigen::Vector3d(delta_x, delta_y, delta_z);
}

double computeDeltaXu(Vertex* u, std::unordered_map<Face, Eigen::Vector3d> fv, WireSculptPlugin& ws) {
    // get each f on v & divide by num faces 
    
    double delta = -1.0f; 

    bool v1 = false;  
    bool v2 = false; 
    bool v3 = false; 

    for (auto f : ws.faces) {
        if (*f.verticies[0] == *u) {
            v1 = true; 
            delta += computeDeltaXuFace(u, f.verticies[1], f.verticies[2]);
        }
        else if (*f.verticies[1] == *u) {
            v2 = true;
            delta += computeDeltaXuFace(u, f.verticies[0], f.verticies[2]);
        }
        else if (*f.verticies[2] == *u) {
            v3 = true;
            delta += computeDeltaXuFace(u, f.verticies[0], f.verticies[1]);
        }
    }

    //double delta = computeDeltaXuFace(h, fv.at(h.getFace()));
    /*Halfedge<Point_3> current = h.next.opposite;
    while (current != h) {
        delta += computeDeltaXuFace(current, fv.at(current.getFace()));
        current = current.next.opposite;
    }*/

    return 0.5 * delta;
}

// what is this
double computeDeltaXuFace(Vertex * curr, Vertex * v1, Vertex * v2) {
    
    auto ahh = curr->mPosition - curr->mPosition;
    auto argh = v1->mPosition - curr->mPosition;

    Eigen::Vector3d c = {ahh[0], ahh[1], ahh[2]};
    Eigen::Vector3d d = {argh[0], argh[1], argh[2]};
    double dotProduct = c.dot(d);
    double magnitudeC = c.norm();
    double magnitudeD = d.norm();

    double cosTheta = dotProduct / (magnitudeC * magnitudeD);

    cosTheta = std::max(std::min(cosTheta, 1.0), -1.0);
    double theta = acos(cosTheta);

    double theta2 = Utilities.computeAngele(v3, v2);
    v3 = Utilities.getVector(e3.opposite);
    double theta1 = Utilities.computeAngele(v1, v3);

    double inner1 = E1.innerProduct(Xj).doubleValue();
    double inner2 = E2.innerProduct(Xj).doubleValue();

    double cot1 = 1 / std::tan(theta1);
    double cot2 = 1 / std::tan(theta2);

    return (cot1 * inner1 + cot2 * inner2);
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