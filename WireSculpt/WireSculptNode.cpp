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

    return MS::kSuccess;
}

MObject WireSculptNode::createMesh(const double& radius, WireSculptPlugin& ws, std::vector<Vertex>& verticies, MObject& outData, MStatus& status) {

    // Making sphere wireframe
    for (auto vertex : verticies) {
      
        MPoint start = vertex.mPosition;

        MPointArray currPoints;
        MIntArray currFaceCounts;
        MIntArray currFaceConnects;

        SphereMesh sphere(start, radius);
        sphere.getMesh(currPoints, currFaceConnects, currFaceConnects);
        sphere.appendToMesh(points, faceCounts, faceConnects);
    }

    // Running A* or TSP and drawing cylinders to map out path
    int testTSP = 1;
    if (testTSP) {
        /* TSP */
        // Choose five random vertices
        
        std::vector<Vertex*> landmarks;
        MPointArray currPoints;
        MIntArray currFaceCounts;
        MIntArray currFaceConnects;

        std::vector<int> chosen(verticies.size());
        int numLMs = chosen.size();
        int index = rand() % numLMs;

        for (int i = 0; i < 5; i++) {
            while (chosen[index] == 1) {
                index = rand() % numLMs;
            }
            chosen[index] = 1;
            landmarks.push_back(&verticies[index]);
                
            // Draw each Landmark Vertex
            SphereMesh sphere(verticies[index].mPosition, radius * 2);
            sphere.getMesh(currPoints, currFaceConnects, currFaceConnects);
            sphere.appendToMesh(points, faceCounts, faceConnects);

        }

        // Run TSP Nearest Neighbors on five vertices
        std::vector<int> tour = ws.FindTspPath(landmarks, 0);

        // Make cylinder meshes to connect vertices found
        for (int t = 0; t < tour.size() - 1; t++) {
            MPoint start = landmarks[tour[t]]->mPosition;
            MPoint end = landmarks[tour[t + 1]]->mPosition;

            CylinderMesh cylinder(start, end, radius * 0.7);
            cylinder.getMesh(currPoints, currFaceConnects, currFaceConnects);
            cylinder.appendToMesh(points, faceCounts, faceConnects);
        }
        MPoint start = landmarks[tour[0]]->mPosition;
        MPoint end = landmarks[tour[tour.size() - 1]]->mPosition;

        CylinderMesh cylinder(start, end, radius * 0.7);
        cylinder.getMesh(currPoints, currFaceConnects, currFaceConnects);
        cylinder.appendToMesh(points, faceCounts, faceConnects);
    }
    else {
        /* A* Path */
        // Choose random vertices to be source and goal
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0, verticies.size() - 1);

        Vertex* source = &verticies[dist6(rng)];
        Vertex* goal = &verticies[dist6(rng)];   // arbitrary

        // Draw the source and goal
        MPointArray currPoints;
        MIntArray currFaceCounts;
        MIntArray currFaceConnects;

        SphereMesh sphere(source->mPosition, radius * 1.7);
        sphere.getMesh(currPoints, currFaceConnects, currFaceConnects);
        sphere.appendToMesh(points, faceCounts, faceConnects);

        SphereMesh sphere2(goal->mPosition, radius * 3);
        sphere2.getMesh(currPoints, currFaceConnects, currFaceConnects);
        sphere2.appendToMesh(points, faceCounts, faceConnects);

        // Run A*
        std::vector<Vertex*> path = ws.FindPath(verticies, source, goal, verticies.size());
        if (path.size() == 0) {
            MGlobal::displayInfo("No path found");
        }
        else {
            for (int i = 0; i < path.size() - 1; i++) {
                MPoint start = path[i]->mPosition;
                MPoint end = path[i + 1]->mPosition;

                MPointArray currPoints;
                MIntArray currFaceCounts;
                MIntArray currFaceConnects;

                CylinderMesh cylinder(start, end, radius * 0.5);
                cylinder.getMesh(currPoints, currFaceConnects, currFaceConnects);
                cylinder.appendToMesh(points, faceCounts, faceConnects);
            }
        }
    }

    // Display wireframe:
    /*for (auto v : verticies) {
        MPoint start = v.mPosition;
        for (auto n : v.neighbors) {
            MPoint end = n.first->mPosition;

            MPointArray currPoints;
            MIntArray currFaceCounts;
            MIntArray currFaceConnects;

            CylinderMesh cylinder(start, end, radius * 0.5);
            cylinder.getMesh(currPoints, currFaceConnects, currFaceConnects);
            cylinder.appendToMesh(points, faceCounts, faceConnects);
        }
    }*/
    
    
    MFnMesh meshFS;
    MObject meshObject = meshFS.create(points.length(), faceCounts.length(), points, faceCounts, faceConnects, outData, &status);

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
        double thicknessVal = thicknessData.asDouble();

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
        createMesh(thicknessVal, ws, *(ws.GetVerticies()), newOutputData, returnStatus);

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