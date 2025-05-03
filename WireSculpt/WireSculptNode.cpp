#include "WireSculptNode.h"
#include "cylinder.h"
#include "sphere.h"
#include "WireSculptPlugin.h"
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNumericData.h>
#include <maya/MFnData.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnMesh.h>
#include <maya/MGlobal.h>
#include <random>


#include <maya/MFnLambertShader.h>
#include <maya/MFnSet.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MItMeshVertex.h>

MTypeId WireSculptNode::id(0x0007F015);
MObject WireSculptNode::inMeshFile;
MObject WireSculptNode::aAttract;
MObject WireSculptNode::bAttract;
MObject WireSculptNode::aRepel;
MObject WireSculptNode::bRepel;
//MObject WireSculptNode::lambda;
//MObject WireSculptNode::K;
//MObject WireSculptNode::M;
MObject WireSculptNode::isAbstract;
MObject WireSculptNode::thickness;
MObject WireSculptNode::outGeom;
MObject WireSculptNode::fov;
MObject WireSculptNode::view;
MObject WireSculptNode::contour;
MObject WireSculptNode::testSC;
MObject  WireSculptNode::proximityThresh;
MObject  WireSculptNode::filterThresh;
MObject  WireSculptNode::maxValThresh;

WireSculptNode::WireSculptNode() : MPxNode()
{
}

WireSculptNode::~WireSculptNode()
{
}

MStatus WireSculptNode::initialize()
{
    MFnNumericAttribute numAttr;
    MFnTypedAttribute typedAttr;
    MFnUnitAttribute unitAttr;

    MStatus returnStatus;

    /** Creating + Adding Attributes */

    // Mesh
    WireSculptNode::inMeshFile = typedAttr.create("inputMeshFile", "im", MFnNumericData::kString, MObject::kNullObj, &returnStatus);
    typedAttr.setUsedAsFilename(true);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating mesh file attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::inMeshFile);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding mesh file attribute\n");
        return returnStatus;
    }
    
    // Range & Steepness Attraction Parameters
    WireSculptNode::aAttract = numAttr.create("rangeAttraction", "ra", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating range attract attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::aAttract);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding range attract attribute\n");
        return returnStatus;
    }

    WireSculptNode::bAttract = numAttr.create("steepAttraction", "sa", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating steepness attract attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::bAttract);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding steepness attract attribute\n");
        return returnStatus;
    }

    // Range & Steepness Repulsion Parameters
    WireSculptNode::aRepel = numAttr.create("rangeRepulsion", "rr", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating range repulsion attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::aRepel);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding range repulsion attribute\n");
        return returnStatus;
    }

    WireSculptNode::bRepel = numAttr.create("steepRepulsion", "sr", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating steepness repulsion attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::bRepel);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding steepness repulsion attribute\n");
        return returnStatus;
    }

    // thickness of the wire path
    WireSculptNode::thickness = numAttr.create("wireThickness", "wt", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating thickness attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::thickness);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding thickness attribute\n");
        return returnStatus;
    }

    /* Feature Lines parameters */

    // Fov
    WireSculptNode::fov = numAttr.create("fieldOfView", "fov", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating fov attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::fov);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding fov attribute\n");
        return returnStatus;
    }

    // View
    WireSculptNode::view = numAttr.create("viewChoice", "view", MFnNumericData::kInt, 0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating view attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::view);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding view attribute\n");
        return returnStatus;
    }

    // Contour Choice
    WireSculptNode::contour = numAttr.create("contourChoice", "contour", MFnNumericData::kInt, 0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating contour attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::contour);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding contour attribute\n");
        return returnStatus;
    }

    // Test SC Value
    WireSculptNode::testSC = numAttr.create("testSC", "tsc", MFnNumericData::kDouble, 0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating testSC attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::testSC);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding testSC attribute\n");
        return returnStatus;
    }

    // Extreme Points Params
    WireSculptNode::isAbstract = numAttr.create("isAbstract", "isa", MFnNumericData::kBoolean, true, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating isAbstract attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::isAbstract);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding isAbstract attribute\n");
        return returnStatus;
    }
    WireSculptNode::proximityThresh = numAttr.create("proximityThresh", "pt", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating extreme proximity threshold attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::proximityThresh);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding extreme proximity threshold attribute\n");
        return returnStatus;
    }

    WireSculptNode::filterThresh = numAttr.create("filterThresh", "ft", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating extreme filter threshold attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::filterThresh);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding extreme filter threshold attribute\n");
        return returnStatus;
    }

    WireSculptNode::maxValThresh = numAttr.create("maxValThresh", "mt", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating extreme max threshold attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::maxValThresh);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding extreme max threshold attribute\n");
        return returnStatus;
    }

    // Output Geometry
    WireSculptNode::outGeom = typedAttr.create("geometry", "geo", MFnData::kMesh, MObject::kNullObj, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating geometry attribute\n");
        return returnStatus;
    }
    typedAttr.setWritable(false);
    typedAttr.setStorable(false);
    returnStatus = addAttribute(WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding geometry attribute\n");
        return returnStatus;
    }


    /** Adding Attribute Affects */
    returnStatus = attributeAffects(WireSculptNode::inMeshFile,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::aAttract,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::bAttract,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::aRepel,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::bRepel,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::thickness,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::fov,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::view,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::contour,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::testSC,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::isAbstract,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::proximityThresh,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::filterThresh,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::maxValThresh,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    return MS::kSuccess;
}

// Visualizing Mesh with Spheres
void WireSculptNode::createWireframeMesh(const double& radius, std::vector<Vertex>& verticies,
    MColorArray* colors, MColor color) {
    for (auto vertex : verticies) {

        MPoint start = vertex.mPosition;

        MPointArray currPoints;
        MIntArray currFaceCounts;
        MIntArray currFaceConnects;

        SphereMesh sphere(start, radius);
        sphere.getMesh(currPoints, currFaceCounts, currFaceConnects);
        sphere.appendToMesh(points, faceCounts, faceConnects);
        int numVerticesThisSphere = currPoints.length();

        for (unsigned int i = 0; i < numVerticesThisSphere; ++i) {
            (*colors).append(color);
        }
    }
}

// Visualizing Contours Lines
void WireSculptNode::createContoursMesh(const double& radius, std::vector<std::pair<vec3f, vec3f>> featureSegments,
    MColorArray* colors, MColor color) {
    for (int index = 0; index < featureSegments.size(); index++) {
        std::pair<vec3f, vec3f> line = featureSegments[index];
        MPoint start(line.first.x, line.first.y, line.first.z);
        MPoint end(line.second.x, line.second.y, line.second.z);

        MPointArray currPoints;
        MIntArray currFaceCounts;
        MIntArray currFaceConnects;

        CylinderMesh cylinder(start, end, radius);
        cylinder.getMesh(currPoints, currFaceCounts, currFaceConnects);
        cylinder.appendToMesh(points, faceCounts, faceConnects);
        int numVerticesThisSphere = currPoints.length();

        for (unsigned int i = 0; i < numVerticesThisSphere; ++i) {
            (*colors).append(color);
        }
    }
}

// Visualizing Mesh Contours Vertices
void WireSculptNode::createFeatureVertsMesh(const double& radius, std::vector<Vertex>& verticies, 
    std::vector<int> featureVertices, MColorArray* colors, MColor color) {
    for (int index = 0; index < featureVertices.size(); index++) {
        Vertex* featureVert = &verticies[featureVertices[index]];
        MPoint start = (*featureVert).mPosition;

        MPointArray currPoints;
        MIntArray currFaceCounts;
        MIntArray currFaceConnects;

        SphereMesh sphere(start, radius * 1.1);
        sphere.getMesh(currPoints, currFaceCounts, currFaceConnects);
        sphere.appendToMesh(points, faceCounts, faceConnects);
        int numVerticesThisSphere = currPoints.length();

        for (unsigned int i = 0; i < numVerticesThisSphere; ++i) {
            (*colors).append(color);
        }
    }
}

// Visualize edge weights after applying feature attraction mesh
void WireSculptNode::createEdgeWeightsMesh(const double& radius, std::vector<Edge>& edges, MColorArray* colors) {
    float maxDist = -1;
    float shortDist = 99999;
    for (auto e : edges) {
        if (e.featureLength > maxDist) {
            maxDist = e.featureLength;
        }
        if (e.featureLength < shortDist) {
            shortDist = e.featureLength;
        }
    }
    for (auto e : edges) {
        float length = e.featureLength;
        const Vertex* v1 = e.endpoints.first;
        const Vertex* v2 = e.endpoints.second;
        MPointArray currPoints;
        MIntArray currFaceCounts;
        MIntArray currFaceConnects;

        MPoint start(v1->mPosition);
        MPoint end(v2->mPosition);

        CylinderMesh cylinder(start, end, radius * 0.8);
        cylinder.getMesh(currPoints, currFaceCounts, currFaceConnects);
        cylinder.appendToMesh(points, faceCounts, faceConnects);

        int numVerticesThisSphere = currPoints.length();
        float r;
        float b;
        for (unsigned int c = 0; c < numVerticesThisSphere; ++c) {
            r = ((float)length - shortDist) / (maxDist - shortDist);
            b = 1.0 - r;
            MColor color(r, 0, b, 0.6f);
            (*colors).append(color);
        }
    }
}

//void WireSculptNode::mapToColors(std::unordered_map<Vertex*, float> colorScheme) {
//    float maxDist = -1;
//    for (auto vc : colorScheme) {
//        if (vc.second > maxDist) {
//            maxDist = vc.second;
//        }
//    }
//    for (auto vc : colorScheme) {
//        float distance = colorScheme[(vc.first)] / maxDist;
//        float r = std::min(std::max(distance, 0.0f), 1.0f);
//        vc.second = 1.0 - r;
//
//        if (r < 0) {
//            MGlobal::displayInfo("Invalid r value");
//        }
//    }
//}

//void WireSculptNode::remapFeatureLengths(std::vector<Edge>& edges, float gamma = 0.5f) {
//    float maxDist = -1;
//    float shortDist = 99999;
//    for (auto e : edges) {
//        if (e.featureLength > maxDist) {
//            maxDist = e.featureLength;
//        }
//        if (e.featureLength < shortDist) {
//            shortDist = e.featureLength;
//        }
//    }
//    for (auto e : edges) {
//        float normalized = (e.featureLength - shortDist) / (maxDist - shortDist);  // maps to [0, 1]
//        float curved = pow(normalized, gamma);
//        float temp = e.featureLength;
//        e.featureLength = curved;
//        
//        MGlobal::displayInfo("Remap: Curved F edge value: " + MString() + curved + ": Normalized F: " + normalized + "; Original: " + temp);
//    }
//}

// Visualizing Heatmap Mesh
void WireSculptNode::createHeatMapMesh(const double& radius, std::unordered_map<Vertex*, float> colorScheme, MColorArray* colors) {
    float r = 1.0;
    float b = 0.0;
    float g = 0.0;
    float maxDist = -1;
    for (auto vc : colorScheme) {
        if (vc.second > maxDist) {
            maxDist = vc.second;
        }
    }
    for (auto vc : colorScheme) {
        MPoint start = (vc.first)->mPosition;

        MPointArray currPoints;
        MIntArray currFaceCounts;
        MIntArray currFaceConnects;

        SphereMesh sphere(start, radius);
        sphere.getMesh(currPoints, currFaceConnects, currFaceConnects);
        int numVerticesThisSphere = currPoints.length();

        for (unsigned int i = 0; i < numVerticesThisSphere; ++i) {
            float distance = colorScheme[(vc.first)];
            r = distance / maxDist; //(maxDist > 1e-6f) ? (distance / maxDist) : 0.0f;
            r = std::min(std::max(r, 0.0f), 1.0f);
            b = 1.0 - r;
           /* r = 1.0 * distance / maxDist;
            b = 1.0 - 1.0 * distance / maxDist;*/

            MColor color(r, g, b);
            (*colors).append(color);
            //MGlobal::displayInfo("Color: " + MString() + "; r: " + r + "; g: " + MString() + g + "; b: " + MString() + b);
        }

        sphere.appendToMesh(points, faceCounts, faceConnects);
    }

}

MObject WireSculptNode::createMesh(const double& radius, const double& aAttract, const double& bAttract, 
    const double& aRepel, const double& bRepel, const double& fovVal, const int& viewChoice,
    const int& contourChoice, const double& testSCVal, const double& proximity, const double& filter,
    const double& maxVal, const bool& isAbstract,
    WireSculptPlugin& ws, const std::string& filePath, 
    std::vector<Vertex>& verticies, std::vector<Edge>& edges, MObject& outData, MStatus& status) {
    
    /* Color variables intitialization */
    MColorArray colorsContours;
    MColorArray colorsHeatMap;

    MColor gray(0.5, 0.5, 0.5);
    MColor red(1.0, 0, 0);

    // Making sphere wireframe
    // createWireframeMesh(radius, verticies, &colorsContours, gray);

    /* Step 1 - Extract Extreme Points */

    /*MGlobal::displayInfo("Extreme points parameters: ");
    MGlobal::displayInfo("Proximity: " + MString() + proximity);
    MGlobal::displayInfo("Filter: " + MString() + filter);
    MGlobal::displayInfo("Max Val: " + MString() + maxVal);*/

    float proxNum = 0.15;
    float filterNum = 0.05;
    float maxNum = 1.0;
    if (isAbstract) {
        proxNum = proximity;
        filterNum = filter;
        maxNum = maxVal; 
    }
    std::vector<int> extremePoints = ws.GetExtremePoints(filePath, proxNum, filterNum, maxNum);
    MGlobal::displayInfo("Extreme points: ");
    for (int i : extremePoints) {
        MGlobal::displayInfo("point: " + MString() + i);

    }
    // Set landmark vertices to be vertices of indexes extremePts chose:
    std::vector<Vertex*> landmarks;
    std::set<Vertex*> landmarksSet;
    int colorIndex = 0;
    for (int index : extremePoints) {
        if (index >= verticies.size()) {
            MGlobal::displayInfo("Error: index too large");
        }
        else {
            //landmarks.push_back(&verticies[index]);
            landmarksSet.insert(&verticies[index]);
            MGlobal::displayInfo("landmark pushed back");

            // Draw each Landmark Vertex
            MPointArray currPoints;
            MIntArray currFaceCounts;
            MIntArray currFaceConnects;

            SphereMesh sphere(verticies[index].mPosition, radius * 2);
            sphere.getMesh(currPoints, currFaceCounts, currFaceConnects);
            sphere.appendToMesh(points, faceCounts, faceConnects);
            int numVerticesThisSphere = currPoints.length();

            for (unsigned int i = 0; i < numVerticesThisSphere; ++i) {
               /* float r = ((float) colorIndex) / (extremePoints.size() - 1.0);
                MColor color(r, 0.0, 1.0 - r);*/
                colorsHeatMap.append(gray);
            }
            colorIndex += 1;
        }
    }

    /* Step 2 - Extract Feature Lines */
    // Draw Contours
    std::vector<std::pair<vec3f, vec3f>> featureSegments = ws.GetContours(fovVal, viewChoice, contourChoice, testSCVal, filePath.c_str());
    // createContoursMesh(radius, featureSegments, &colorsContours, gray);

    std::vector<int> featureVertices = ws.processSegments(&featureSegments);
    MColor red_alpha(1, 0, 0, 0.7);
    //createFeatureVertsMesh(radius, verticies, featureVertices, &colorsHeatMap, red_alpha);

    MGlobal::displayInfo("Finished: set up contours");

    // Add feature verts into landmarks
    float fVThreshold = 0.2;
    //int seed = 42;  // Seed
    //std::mt19937 gen(seed); // Mersenne Twister engine
    std::uniform_real_distribution<> dis(0.0, 1.0); // Range [0, 1)
    for (int index : featureVertices) {
        std::mt19937 gen(index); // Mersenne Twister engine
        double rand = dis(gen); // Generate random number
        //std::mt19937 gen(index); // Seed with index
        //uint32_t randInt = gen(); // Get deterministic integer
        //double rand = static_cast<double>(randInt) / gen.max();
        if (rand <= fVThreshold) {
            landmarksSet.insert(&verticies[index]);
        }
    }
    
    for (Vertex* v : landmarksSet) {
        landmarks.push_back(v);
    }


    /* Step 3 - Build Heat Map from Feature Lines Vertices */
    std::unordered_map<Vertex*, float> contourHeatMap = ws.GetHeatMapDistance(ws, &featureVertices);
    //createHeatMapMesh(radius, contourHeatMap, &colorsHeatMap);

    /* Step 4 - Compute Feature Attraction Weights */
    float lBar = 0; 
    for (auto e : edges) {      // find average edge length
        lBar += e.getLength();
    }
    lBar /= edges.size();

    for (int i = 0; i < verticies.size(); i++) {    // compute feature attraction weight for each vertex
        Vertex* vert = &verticies[i];
        float distance = contourHeatMap[vert] < 0 ? 0 : contourHeatMap[vert];   // clamping - getting rid of negative values
        vert->wAttract = (aAttract / (1.0 + std::exp(-bAttract * distance / lBar))) + (1 - aAttract);
    }

    for (auto& e: edges) {
        const Vertex* vi = e.endpoints.first;
        const Vertex* vj = e.endpoints.second;
        e.featureLength = 0.5 * (vi->wAttract + vj->wAttract) * e.getLength();
        e.warpedLength = e.featureLength;
    }

    // Begin visualize - colors of contours heat map on vertices and warped feature weights on edges
    std::unordered_map<Vertex*, float> featureWtsHeatMap;
    for (int i = 0; i < verticies.size(); i++) {
        Vertex* vert = &verticies[i];
        featureWtsHeatMap[vert] = vert->wAttract;
    }
    //createHeatMapMesh(radius, featureWtsHeatMap, &colorsHeatMap);
    //createEdgeWeightsMesh(radius, edges, &colorsHeatMap);
    // End visualize

    /* Step 5 - Run TSP Optimized Nearest Neighbors on landmark vertices */ 
    std::vector<int> tour = ws.TwoOptTspPath(landmarks, 0, 20);     // max 20 iterations
    std::vector<int> wirePath;                                      // accumulated path
    
    for (int t = 0; t < tour.size(); t++) {
        int index1;
        int index2;

        // if last element, cycle back to first node
        if (t == tour.size() - 1) {
            index1 = tour.size() - 1;
            index2 = 0;
        }
        else {
            index1 = t;
            index2 = t + 1;
        }

        Vertex* source = landmarks[tour[index1]];
        Vertex* goal = landmarks[tour[index2]];

        /* Step 5a - Adjust edge weights via path repulsion */ 
        // Call heatmap on existing path
        if (wirePath.size() > 0) {
            std::unordered_map<Vertex*, float> pathHeatMap = ws.GetHeatMapDistance(ws, &wirePath);
            //createHeatMapMesh(radius, pathHeatMap, &colorsHeatMap);
            //mapToColors(pathHeatMap);

            // Update warpedLengths
            for (int i = 0; i < verticies.size(); i++) {    // compute feature attraction weight for each vertex
                Vertex* vert = &verticies[i];
                float distance = pathHeatMap[vert] < 0 ? 0 : pathHeatMap[vert];
                vert->wRepel = (aRepel / (1.0 + std::exp(-bRepel * distance / lBar))) + (1.0 - aRepel);   // issue ? - doubles in float math!!
            }
            for (auto& e : edges) {
                const Vertex* vi = e.endpoints.first;
                const Vertex* vj = e.endpoints.second;
                e.warpedLength = (vi->wAttract + vj->wAttract) * e.getLength() / (vi->wRepel + vj->wRepel);   // IS THIS FEATURE LENGTH OR ORIGINAL LENGTH??
            }
        }
        
        // Run A* between each of the vertices
        std::vector<int> path = ws.FindPath(verticies, source, goal, verticies.size());
        wirePath.insert(wirePath.end(), path.begin(), path.end());  // concatenate current path to accumulated wirePath

        if (path.size() == 0) {
            MGlobal::displayInfo("No path found");
        }
        else {
            // Draw path
            for (int i = 0; i < path.size() - 1; i++) {
                Vertex* v1 = &verticies[path[i]];
                Vertex* v2 = &verticies[path[i + 1]];
                MPoint start = v1->mPosition;
                MPoint end = v2->mPosition;

                MPointArray currPoints;
                MIntArray currFaceCounts;
                MIntArray currFaceConnects;

                CylinderMesh cylinder(start, end, radius * 0.5);
                cylinder.getMesh(currPoints, currFaceCounts, currFaceConnects);
                cylinder.appendToMesh(points, faceCounts, faceConnects);

                int numVerticesThisSphere = currPoints.length();
                for (unsigned int c = 0; c < numVerticesThisSphere; ++c) {
                    float r = ((float) i) / (path.size() - 2);
                    float b = 1.0 - r;
                    //MColor color(r, 0, b); //(0, g, b);
                    colorsHeatMap.append(gray);
                }
            }
        }
    }
    
    MFnMesh meshFn;
    MObject meshObject = meshFn.create(points.length(), faceCounts.length(), points, faceCounts, faceConnects, outData, &status);

    /* Color Data */
    // Get all vertex indices
    MItMeshVertex itVertex(meshObject, &status);

    MIntArray vertexIndices;
    while (!itVertex.isDone()) {
        vertexIndices.append(itVertex.index());
        itVertex.next();
    }

    // Create or use the default color set
    MString colorSetName = "colorSet1";
    meshFn.createColorSetWithName(colorSetName, nullptr, nullptr, &status);

    // Assign colors to vertices
    status = meshFn.setVertexColors(colorsHeatMap, vertexIndices);
    //status = meshFn.setVertexColors(colorsContours, vertexIndices);

    return meshObject;
}

MStatus WireSculptNode::compute(const MPlug& plug, MDataBlock& data) {
    MStatus returnStatus;

    // Check if plug is geometry
    // TODO: test how to make things change based on the input that got changed
    if (plug == outGeom) {

        /* Clear points and face data */
        points.clear();
        faceCounts.clear();
        faceConnects.clear();

        /* Get inputs */
        // Input Mesh
        MDataHandle fileData = data.inputValue(inMeshFile, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting grammar data handle\n");
            return returnStatus;
        }
        MString meshFilePath = fileData.asString();
        std::string meshFilePathStr = meshFilePath.asChar();
        MString objFilePath = MString() + meshFilePathStr.c_str();
        MGlobal::displayInfo("File path: " + objFilePath);

        // Range Attract
        MDataHandle aAttractData = data.inputValue(aAttract, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting range attract data handle\n");
            return returnStatus;
        }
        double aAttractVal = aAttractData.asDouble();

        // Steepness Attract
        MDataHandle bAttractData = data.inputValue(bAttract, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting steep attract data handle\n");
            return returnStatus;
        }
        double bAttractVal = bAttractData.asDouble();

        // Range Repel
        MDataHandle aRepelData = data.inputValue(aRepel, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting range repel data handle\n");
            return returnStatus;
        }
        double aRepelVal = aRepelData.asDouble();

        // Steepness Repel
        MDataHandle bRepelData = data.inputValue(bRepel, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting steep repel data handle\n");
            return returnStatus;
        }
        double bRepelVal = bRepelData.asDouble();

        // Thickness
        MDataHandle thicknessData = data.inputValue(thickness, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting thickness data handle\n");
            return returnStatus;
        }
        double thicknessVal = thicknessData.asDouble() * 0.1f;

        // FOV
        MDataHandle fovData = data.inputValue(fov, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting fov data handle\n");
            return returnStatus;
        }
        double fovVal = fovData.asDouble();

        // View Choice
        MDataHandle viewData = data.inputValue(view, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting view data handle\n");
            return returnStatus;
        }
        double viewVal = viewData.asInt();

        // Contour Choice
        MDataHandle contourData = data.inputValue(contour, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting contour data handle\n");
            return returnStatus;
        }
        double contourVal = contourData.asInt();

        // Contour Choice
        MDataHandle testSCData = data.inputValue(testSC, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting testSC data handle\n");
            return returnStatus;
        }
        double testSCVal = testSCData.asDouble() * 0.1f;

        MDataHandle isAbstractData = data.inputValue(isAbstract, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting proximityThresh data handle\n");
            return returnStatus;
        }
        bool isAbstractVal = isAbstractData.asBool();

        MDataHandle proximityThreshData = data.inputValue(proximityThresh, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting proximityThresh data handle\n");
            return returnStatus;
        }
        double proximityThreshVal = proximityThreshData.asDouble() * 0.1f;

        MDataHandle filterThreshData = data.inputValue(filterThresh, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting filterThresh data handle\n");
            return returnStatus;
        }
        double filterThreshVal = filterThreshData.asDouble() * 0.1f;

        MDataHandle maxValThreshData = data.inputValue(maxValThresh, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting maxValThresh data handle\n");
            return returnStatus;
        }
        double maxValThreshVal = maxValThreshData.asDouble() * 0.1f;

        /* Process file */
        WireSculptPlugin ws = WireSculptPlugin();
        bool returnVal = ws.ProcessFile(meshFilePathStr);
        MString errorMsg = "\n";
        if (returnVal) {
            errorMsg += "file processed successfully";
            MGlobal::displayInfo(errorMsg);
            MGlobal::displayInfo(MString() + ws.GetVerticies()->size());
        }
        else {
            errorMsg += "issue with file format or contents";
            MGlobal::displayInfo(errorMsg);
        }

        /* Get output object - geometry */
        MDataHandle outGeometry = data.outputValue(outGeom, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting geometry data handle\n");
            return returnStatus;
        }
        MFnMeshData dataCreator;
        MObject newOutputData = dataCreator.create(&returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR creating output geometry data\n");
            return returnStatus;
        }

        // Create new geometry
        createMesh(thicknessVal, aAttractVal, bAttractVal, aRepelVal, bRepelVal, fovVal, viewVal, contourVal, testSCVal, 
            proximityThreshVal, filterThreshVal, maxValThreshVal, isAbstractVal,
            ws, meshFilePathStr, *(ws.GetVerticies()), *(ws.GetEdges()), newOutputData, returnStatus);

        if (!returnStatus) {
            returnStatus.perror("ERROR creating new mesh\n");
            return returnStatus;
        }
        outGeometry.set(newOutputData);
        data.setClean(plug);
    }
    else {
        return MS::kUnknownParameter;
    }
    return MS::kSuccess;

}