#include "WireSculptNode.h"
#include "cylinder.h"
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNumericData.h>
#include <maya/MFnData.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnMesh.h>
#include <maya/MGlobal.h>


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

MObject WireSculptNode::createMesh(const double& radius, MObject& outData, MStatus& status) {
    MPoint start({0, 0, 0 });
    MPoint end({ 1, 1, 1 });
    MPointArray currPoints;
    MIntArray currFaceCounts;
    MIntArray currFaceConnects;
    CylinderMesh cylinder(start, end, 0.12);
    cylinder.getMesh(currPoints, currFaceConnects, currFaceConnects);
    cylinder.appendToMesh(points, faceCounts, faceConnects);
    
    MFnMesh meshFS;
    return meshFS.create(points.length(), faceCounts.length(), points, faceCounts, faceConnects, outData, &status);
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
        MDataHandle grammarData = data.inputValue(inMeshFile, &returnStatus);
        if (!returnStatus) {
            returnStatus.perror("ERROR getting grammar data handle\n");
            return returnStatus;
        }
        MString meshFilePath = grammarData.asString();
        std::string meshFilePathStr = meshFilePath.asChar();
        MString grammar = MString() + meshFilePathStr.c_str();

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

        // TODO - idk why i didnt add the other parameters... i must have forgot??

        // Need to create new geometry
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

        createMesh(aAttractVal, newOutputData, returnStatus);

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