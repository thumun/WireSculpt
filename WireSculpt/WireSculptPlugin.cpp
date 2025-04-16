#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <maya/MPoint.h>
#include <maya/MVector.h>

#include "WireSculptPlugin.h"
#include "ExtremePoints.h"
#include "Contours.h"

#include <igl/read_triangle_mesh.h>
#include <Eigen/Core>
#include <igl/is_vertex_manifold.h>
#include <maya/MGlobal.h>

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

void WireSculptPlugin::resetIDs() {
    Vertex::lastId = 0; 
    Edge::lastId = 0; 
    Face::lastId = 0; 
}

WireSculptPlugin::WireSculptPlugin()
{
    geodesicData = std::make_unique<igl::HeatGeodesicsData<double>>();
}

// processes the obj file 
// returns false: issue with file (contents or type) 
// returns true: success 
bool WireSculptPlugin::ProcessFile(std::string filePath) {

    resetIDs();

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
            verticies.push_back(Vertex(MPoint(pos[pos.size() - 1]), false));
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

            //if (faceVerts.size() < 3) return;

            // Compute face normal
            MVector v0 = pos[faceVerts[0]-1];
            MVector v1 = pos[faceVerts[1]-1];
            MVector v2 = pos[faceVerts[2]-1];
            MVector crossProd = (v1 - v0) ^ (v2 - v0);
            crossProd.normalize();

            std::vector<float> norm;
            norm.push_back(crossProd.x); 
            norm.push_back(crossProd.y);
            norm.push_back(crossProd.z);

            // Triangulate the face (assumes convex polygons)
            for (size_t i = 1; i + 1 < faceVerts.size(); i++) {
                faces.push_back(Face(&verticies[faceVerts[0]-1],
                    &verticies[faceVerts[i]-1],
                    &verticies[faceVerts[i + 1]-1],
                    norm));
            }

            edges.reserve(verticies.size() * (verticies.size() - 1) * 0.5f);

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

std::vector<int> WireSculptPlugin::GetExtremePoints(const std::string& filePath) {
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

    igl::heat_geodesics_precompute(V, F, *geodesicData.get());

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

    return extreme_points;
}

// Optimized Tour
std::vector<int> WireSculptPlugin::TwoOptTspPath(std::vector<Vertex*> landmarks, int start, int maxIters) {
    std::vector<int> tour = FindTspPath(landmarks, start);

    int iters = 0;
    int numLMs = tour.size();

    while (iters < maxIters) {
        for (int i = 1; i < numLMs; i++) {
            for (int j = i + 1; j < numLMs - 1; j++) {
                Vertex& lm1 = *landmarks[i - 1];
                Vertex& lm2 = *landmarks[i];
                Vertex& lm3 = *landmarks[j];
                Vertex& lm4 = *landmarks[j + 1];
                if (computeLMDistance(lm1, lm3) + computeLMDistance(lm2, lm4)
                    < computeLMDistance(lm1, lm2) + computeLMDistance(lm3, lm4)) {
                    tour = swapEdge(tour, i, j);
                }
            }
        }
        iters++;
    }
    return tour;
}

std::vector<int> WireSculptPlugin::swapEdge(std::vector<int> tour, int i, int j) {
    std::vector<int> newTour;

    for (int p = 0; p < i; p++) {
        newTour.push_back(tour[p]);
    }
    
    int q = 0;
    for (int p = i; p < j + 1; p++) {
        newTour.push_back(tour[j - q]);
        q++;
    }

    for (int p = j + 1; p < tour.size(); p++) {
        newTour.push_back(tour[p]);
    }

    return newTour;
}


// Nearest Insertion Tour
std::vector<int> WireSculptPlugin::FindTspPath(std::vector<Vertex*> landmarks, int start) {
    int numLandmarks = landmarks.size();
    std::vector<int> tour;
    std::vector<int> used(numLandmarks);

    tour.push_back(start);
    used[start] = 1;

    int newLM;
    int indexNewLM;
    int numVisitedLM = 1;

    while (numVisitedLM != numLandmarks) {
        if (numVisitedLM == 1) {
            newLM = findNearestNeighbor(start, landmarks, used);
            tour.push_back(newLM);
        }
        else {
            // Selection
            newLM = pickUnvisitedCity(used);

            // Insertion
            indexNewLM = findMinTriangularDistanceEdge(newLM, tour, landmarks);
            if (indexNewLM > tour.size()) {
                MGlobal::displayInfo("index is too big");
            }
            else {
                tour.insert(tour.begin() + indexNewLM, newLM);

            }
        }
        used[newLM] = 1;
        numVisitedLM++;
    }

    return tour;
}

// replace this with FindPath if doesn't perform well
int WireSculptPlugin::findNearestNeighbor(int indexOfVertex, std::vector<Vertex*> landmarks, std::vector<int> used) {
    int bestNeighborIndex = -1;
    for (unsigned int i = 0; i < landmarks.size(); i++) {
        if (used[i] == 0) {
            if (bestNeighborIndex == -1) {
                bestNeighborIndex = i;
            }
            else if (computeLMDistance(*landmarks[indexOfVertex], *landmarks[i]) < 
                computeLMDistance(*landmarks[indexOfVertex], *landmarks[bestNeighborIndex])) {
                bestNeighborIndex = i;
            }
        }
    }
    
    return bestNeighborIndex;
}

int WireSculptPlugin::pickUnvisitedCity(std::vector<int> used) {
    int numLMs = used.size();
    int index = rand() % numLMs;
    while (used[index] == 1) {
        index = rand() % numLMs;
    }
    return index;
}

// return the index in which newLM will be inserted in tour 
// that will minimize the insertion cost for the tour
int WireSculptPlugin::findMinTriangularDistanceEdge(int newLM, std::vector<int> tour, std::vector<Vertex*> landmarks) {
    int indexToInsert = 0;
    float minCost = 99999;
    float currCost;

    for (int i = 0; i < tour.size() - 1; i++) {
        currCost = insertionCost(*landmarks[tour[i]], *landmarks[tour[i + 1]], *landmarks[newLM]);
        if (currCost < minCost) {
            minCost = currCost;
            indexToInsert = i + 1;
        }
    }
    return indexToInsert;
}
float WireSculptPlugin::insertionCost(Vertex& lm1, Vertex& lm2, Vertex& lm3) {
    return computeLMDistance(lm1, lm3) + computeLMDistance(lm3, lm2) - computeLMDistance(lm1, lm2);
}

float WireSculptPlugin::computeLMDistance(Vertex& lm1, Vertex& lm2) {
    MPoint start = lm1.mPosition;
    MPoint end = lm2.mPosition;
    return start.distanceTo(end);
}


// A* Path Finding Algorithm
std::vector<Vertex*> WireSculptPlugin::FindPath(std::vector<Vertex>& verticies, Vertex* start, Vertex* goal, int vertexCount)
{
    // Initialize the open and closed lists
    std::priority_queue<Vertex*, std::vector<Vertex*>, VertexPtrCompare> openList;
    std::vector<bool> inOpenList(vertexCount, false);

    // Keep track of parents
    std::vector<int> parentList(vertexCount, -1);

    // Reset all f,g,h values of all verticies
    for (auto& v : verticies) {
        v.resetFGH();
    }

    start->f = 0;
    start->g = 0;
    start->h = 0;
    
    // Start vertex
    openList.push(start);
    inOpenList[start->id] = true;

    // Main loop
    while (!openList.empty()) {

        // Get vertex with lowest f value from open list
        Vertex* current = openList.top();
        openList.pop();
        inOpenList[current->id] = false;

        // Check if current vertex is goal
        if (current->id == goal->id) {

            // Reconstruct the path
            std::vector<Vertex*> path;
            while (!(current->id == start->id))
            {
                path.push_back(current);

                if (parentList[current->id] >= 0) {
                    current = &verticies[parentList[current->id]];	// assuming that verticies are stored in order of id

                }
                else {
                    return path;
                }
            }

            path.push_back(start);
            reverse(path.begin(), path.end());
            return path;
            
        }

        // Explore neighbors
        for (auto& neighbor : current->neighbors) {
            Vertex* nVert = neighbor.first;
            Edge* nEdge = neighbor.second;

            float newG = current->g + nEdge->warpedLength;
            if (newG < nVert->g) {
                float newH = 0; //(goal->mPosition - nVert->mPosition).length();	// for now - the distance from n to goal
                float newF = newG + newH;
                nVert->g = newG;
                nVert->h = newH;
                nVert->f = newF;
                parentList[nVert->id] = current->id;
                        
                if (!inOpenList[nVert->id]) {
                    openList.push(nVert);
                    inOpenList[nVert->id] = true;
                }
            }

        }
    }
    return std::vector<Vertex*>();
}

std::unordered_map<Vertex*, float> WireSculptPlugin::GetHeatMapDistance(WireSculptPlugin& ws) {

    Eigen::VectorXi gamma (1); // init w/ size 1 
    gamma << 0; // source is vert 0 

    //Eigen::VectorXi gamma = { 0 };

    Eigen::VectorXd D;

    igl::heat_geodesics_solve(*geodesicData.get(), gamma, D);

    std::unordered_map<Vertex*, float> outHeatInfo; 

    // converting D to be in our prev format
    for (int i = 0; i < D.size(); i++) {
        outHeatInfo.insert({&verticies[i], D(i)});
    }

    /*HeatMapDist dist = HeatMapDist(ws);
    dist.heatDiffusion(0);
    dist.computePhi(0, ws);
    return dist.colorScheme(ws, 'd');*/

    return outHeatInfo; 
}

std::unordered_map<Vertex*, float> WireSculptPlugin::GetHeatMapDistance(WireSculptPlugin& ws, std::vector<int>* segments) {
    Eigen::VectorXi gamma(segments->size()); // init w/ size 1 
    //gamma << 0; // source is vert 0 

    for (int i = 0; i < segments->size(); ++i) {
        gamma[i] = (*segments)[i];
    }

    Eigen::VectorXd D;

    igl::heat_geodesics_solve(*geodesicData.get(), gamma, D);

    std::unordered_map<Vertex*, float> outHeatInfo;

    // converting D to be in our prev format
    for (int i = 0; i < D.size(); i++) {
        outHeatInfo.insert({ &verticies[i], D(i) });
    }

    return outHeatInfo;
}

std::vector<std::pair<vec3f, vec3f>> WireSculptPlugin::GetContours(float fovChoice, int viewChoice, int contoursChoice, float testSCChoice, const char* filename) {
    Contours contour(fovChoice, viewChoice, contoursChoice, testSCChoice, filename);

    std::vector<std::pair<vec3f, vec3f>> featureSegments;

    if (contour.featureLines.size() == 0) {
        MGlobal::displayInfo("No feature lines");
        return std::vector<std::pair<vec3f, vec3f>>();
    }

    for (int index = 0; index < contour.featureLines.size(); index++) {
        std::vector<std::vector<float>> featurePoints = contour.featureLines[index];
        for (int i = 0; i < featurePoints.size() - 1; i += 2) {
            std::vector<float> p1 = featurePoints[i];
            std::vector<float> p2 = featurePoints[i + 1];

            vec3f start(p1[0], p1[1], p1[2]);
            vec3f end(p2[0], p2[1], p2[2]);

            featureSegments.push_back({ start, end });
        }
    }
    return featureSegments;
}

int WireSculptPlugin::findClosestVertex(float x, float y, float z) {
    Eigen::Vector3d vertPos = { x, y, z};
    int storeIndex = -1;
    double dist = std::numeric_limits<double>::max();

    // find based on pos
    int index = 0;
    for (auto& v : this->verticies) {
        Eigen::Vector3d comparePos = { v.mPosition.x, v.mPosition.y, v.mPosition.z };
        auto distCompare = (comparePos - vertPos).norm();

        if (distCompare < dist) {
            storeIndex = index;
            dist = distCompare;
        }
        index += 1;
    }
    return storeIndex;
}
std::vector<int> WireSculptPlugin::processSegments(std::vector<std::pair<vec3f, vec3f>>* segments) {
    std::vector<int> vertIndices;
    std::vector<bool> inIndexList(this->verticies.size(), false);

    // Find index of closest vertex for start and end point
    for (int i = 0; i < segments->size(); i++) {
        vec3f start = (*segments)[i].first;   
        int startIndex = findClosestVertex(start.x, start.y, start.z);
        if (!inIndexList[startIndex]) {
            inIndexList[startIndex] = true;
            vertIndices.push_back(startIndex);
        }

        vec3f end = (*segments)[i].second;  
        int endIndex = findClosestVertex(end.x, end.y, end.z);
        if (!inIndexList[endIndex]) {
            inIndexList[endIndex] = true;
            vertIndices.push_back(endIndex);
        }
    }

    return vertIndices;
}

std::vector<Vertex>* WireSculptPlugin::GetVerticies() {
    return &(this->verticies);
}

std::vector<Edge>* WireSculptPlugin::GetEdges() {
    return &(this->edges);
}

#if EXEDEBUG
int main() {
    WireSculptPlugin ws = WireSculptPlugin();
    ws.ProcessFile("D:/CGGT/AdvTopics/WireSculpt/testobj/subSphere.obj");
    //ws.GetExtremePoints("D:/CGGT/AdvTopics/WireSculpt/testobj/cow.obj");
    //ws.ProcessFile("C:/Users/53cla/Documents/Penn/CIS_6600/Authoring Tool/WireSculpt/Test objs/suzanne.obj");
    //std::vector<Vertex>* verticies = ws.GetVerticies();
    //Vertex* source = &verticies[2];
    //Vertex* goal = &verticies[5];   // arbitrary
    //std::vector<Vertex*> path = ws.FindPath((*verticies), source, goal, (*verticies).size());

    ws.GetExtremePoints("D:/CGGT/AdvTopics/WireSculpt/testobj/subSphere.obj");

    auto featSeg = ws.GetContours(0.7f, 0, 1, 0.1f, "D:/CGGT/AdvTopics/WireSculpt/testobj/subSphere.obj");
    auto test = ws.processSegments(&featSeg);

    //ws.GetHeatMapDistance(ws, &test);
    ws.GetHeatMapDistance(ws);

}
#endif // EXEDEBUG