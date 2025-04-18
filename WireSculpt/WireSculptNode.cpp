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
MObject WireSculptNode::lambda;
MObject WireSculptNode::K;
MObject WireSculptNode::M;
MObject WireSculptNode::thickness;
MObject WireSculptNode::outGeom;
MObject WireSculptNode::fov;
MObject WireSculptNode::view;
MObject WireSculptNode::contour;
MObject WireSculptNode::testSC;

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

    // Lambda - error
    WireSculptNode::lambda = numAttr.create("lambdaError", "le", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating lambda attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::lambda);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding lambda error attribute\n");
        return returnStatus;
    }

    // K - total paths to expand in iteration
    WireSculptNode::K = numAttr.create("kValue", "kV", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating K Value attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::K);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding K Value attribute\n");
        return returnStatus;
    }

    // M - total number of paths to keep in an iteration
    WireSculptNode::M = numAttr.create("mValue", "mV", MFnNumericData::kDouble, 0.0, &returnStatus);
    if (!returnStatus) {
        returnStatus.perror("ERROR creating M Value attribute\n");
        return returnStatus;
    }
    returnStatus = addAttribute(WireSculptNode::M);
    if (!returnStatus) {
        returnStatus.perror("ERROR adding M Value attribute\n");
        return returnStatus;
    }

    // M - total number of paths to keep in an iteration
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

    returnStatus = attributeAffects(WireSculptNode::lambda,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::K,
        WireSculptNode::outGeom);
    if (!returnStatus) {
        returnStatus.perror("ERROR in attributeAffects\n");
        return returnStatus;
    }

    returnStatus = attributeAffects(WireSculptNode::M,
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

        SphereMesh sphere(start, radius);
        sphere.getMesh(currPoints, currFaceCounts, currFaceConnects);
        sphere.appendToMesh(points, faceCounts, faceConnects);
        int numVerticesThisSphere = currPoints.length();

        for (unsigned int i = 0; i < numVerticesThisSphere; ++i) {
            (*colors).append(color);
        }
    }
}
void WireSculptNode::mapToColors(std::unordered_map<Vertex*, float> colorScheme) {
    float maxDist = -1;
    for (auto vc : colorScheme) {
        if (vc.second > maxDist) {
            maxDist = vc.second;
        }
    }
    for (auto vc : colorScheme) {
        float distance = colorScheme[(vc.first)];
        float r = std::min(std::max(distance, 0.0f), 1.0f);
        vc.second = r;
    }
}

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
            r = distance; //(maxDist > 1e-6f) ? (distance / maxDist) : 0.0f;
            r = std::min(std::max(r, 0.0f), 1.0f);
            b = 1.0 - r;
           /* r = 1.0 * distance / maxDist;
            b = 1.0 - 1.0 * distance / maxDist;*/

            MColor color(r, g, b);
            (*colors).append(color);
        }

        sphere.appendToMesh(points, faceCounts, faceConnects);
    }

}

MObject WireSculptNode::createMesh(const double& radius, const double& aAttract, const double& bAttract, 
    const double& aRepel, const double& bRepel, const double& fovVal, const int& viewChoice,
    const int& contourChoice, const double& testSCVal, 
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
    std::vector<int> extremePoints = ws.GetExtremePoints(filePath);
    MGlobal::displayInfo("Extreme points: ");
    for (int i : extremePoints) {
        MGlobal::displayInfo("point: " + MString() + i);

    }
    // Set landmark vertices to be vertices of indexes extremePts chose:
    std::vector<Vertex*> landmarks;

    for (int index : extremePoints) {
        if (index >= verticies.size()) {
            MGlobal::displayInfo("Error: index too large");
        }
        else {
            landmarks.push_back(&verticies[index]);
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
                colorsHeatMap.append(gray);
            }

        }
    }

    /* Step 2 - Extract Feature Lines */
    // Draw Contours
    std::vector<std::pair<vec3f, vec3f>> featureSegments = ws.GetContours(fovVal, viewChoice, contourChoice, testSCVal, filePath.c_str());
    // createContoursMesh(radius, featureSegments, &colorsContours, gray);

    std::vector<int> featureVertices = ws.processSegments(&featureSegments);
    // createFeatureVertsMesh(radius, verticies, featureVertices, &colorsContours, red);

    MGlobal::displayInfo("Finished: set up contours");


    /* Step 3 - Build Heat Map from Feature Lines Vertices */
    //std::unordered_map<Vertex*, float> colorScheme = ws.GetHeatMapDistance(ws);
    std::unordered_map<Vertex*, float> contourHeatMap = ws.GetHeatMapDistance(ws, &featureVertices);
    createHeatMapMesh(radius, contourHeatMap, &colorsHeatMap);
    mapToColors(contourHeatMap);

    /* Step 4 - Compute Feature Attraction Weights */
    MGlobal::displayInfo("edge size: " + MString() + edges.size());

    float lBar = 0; 
    for (auto e : edges) {      // find average edge length
        lBar += e.getLength();
    }
    lBar /= edges.size();

    for (int i = 0; i < verticies.size(); i++) {    // compute feature attraction weight for each vertex
        Vertex* vert = &verticies[i];
        float distance = contourHeatMap[vert];
        vert->wAttract = (aAttract / (1.0 + std::exp(-bAttract * distance / lBar))) + (1 - aAttract);   // issue ? - doubles in float math!!
        MGlobal::displayInfo("Wattract: " + MString() + vert->wAttract);
    }

    for (auto e: edges) {
        const Vertex* vi = e.endpoints.first;
        const Vertex* vj = e.endpoints.second;
        e.featureLength = 0.5 * (vi->wAttract + vj->wAttract) * e.getLength();   // does this change the correct edges? (pointer issues?)
        e.warpedLength = e.featureLength;   // does this change the correct edges? (pointer issues?)
    }

    // print warped edge lengths for each vertex to check
    /*for (auto v : verticies) {
        for (auto ve : v.neighbors) {
            Edge* e = ve.second;
            MGlobal::displayInfo("Warped edge feature lengths: " + MString() + e->featureLength);
        }
    }*/

    /* Step 5 - Run TSP Optimized Nearest Neighbors on landmark vertices */ 
    std::vector<int> tour = ws.TwoOptTspPath(landmarks, 0, 20); // max 20 iterations
    std::vector<int> wirePath;
    
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
        }
        std::unordered_map<Vertex*, float> pathHeatMap = ws.GetHeatMapDistance(ws, &featureVertices);
        //createHeatMapMesh(radius, pathHeatMap, &colorsHeatMap);
        mapToColors(pathHeatMap);

        // Update warpedLengths
        for (int i = 0; i < verticies.size(); i++) {    // compute feature attraction weight for each vertex
            Vertex* vert = &verticies[i];
            float distance = pathHeatMap[vert];
            vert->wRepel = (aRepel / (1.0 + std::exp(-bRepel * distance / lBar))) + (1 - aRepel);   // issue ? - doubles in float math!!
        }
        for (auto e : edges) {
            const Vertex* vi = e.endpoints.first;
            const Vertex* vj = e.endpoints.second;
            e.warpedLength = (vi->wAttract + vj->wAttract) * e.featureLength / (vi->wRepel + vj->wRepel);   // IS THIS FEATURE LENGTH OR ORIGINAL LENGTH??
        }

        // Run A* between each of the vertices
        std::vector<int> path = ws.FindPath(verticies, source, goal, verticies.size());
        wirePath.insert(wirePath.end(), path.begin(), path.end());  // concatenate current path to accumulated wirePath

        if (path.size() == 0) {
            MGlobal::displayInfo("No path found");
        }
        else {
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
                for (unsigned int i = 0; i < numVerticesThisSphere; ++i) {
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

        // Lambda
        MDataHandle lambdaData = data.inputValue(lambda, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting lambda data handle\n");
            return returnStatus;
        }
        double lambdaVal = lambdaData.asDouble();

        // K
        MDataHandle kData = data.inputValue(K, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting k data handle\n");
            return returnStatus;
        }
        double kVal = kData.asDouble();

        // M
        MDataHandle mData = data.inputValue(M, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting m data handle\n");
            return returnStatus;
        }
        double mVal = mData.asDouble();

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