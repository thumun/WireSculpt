
#ifndef CreateWireSculptNode_H_
#define CreateWireSculptNode_H_
#include "Vertex.h"
#include "WireSculptPlugin.h"
#include <maya/MPxNode.h>
#include <maya/MPointArray.h>
#include <vector>

class WireSculptNode : public MPxNode
{
public:
    WireSculptNode();
    virtual ~WireSculptNode();
    static void* creator() { return new WireSculptNode; }
    static MStatus initialize();

    MStatus compute(const MPlug& plug, MDataBlock& data) override;

    static MTypeId  id;
    static MObject  inMeshFile;
    static MObject  aAttract;
    static MObject  bAttract;
    static MObject  aRepel;
    static MObject  bRepel;
    static MObject  lambda;
    static MObject  K;
    static MObject  M;
    static MObject  thickness;
    static MObject  outGeom;
    static MObject  fov;
    static MObject  view;
    static MObject  contour;
    static MObject  testSC;

protected:
    MObject createMesh(const double& radius, const double& fovVal, const int& viewChoice, 
        const int& contourChoice, const double& testSCVal, 
        WireSculptPlugin& ws, const std::string& filePath, std::vector<Vertex>& verticies, 
        MObject& outData, MStatus& status);

    // CreateMesh helpers for visualization
    void createWireframeMesh(const double& radius, std::vector<Vertex>& verticies,
        MColorArray* colors, MColor color);
    void createContoursMesh(const double& radius, std::vector<std::pair<vec3f, vec3f>> featureSegments,
        MColorArray* colors, MColor color);
    void createFeatureVertsMesh(const double& radius, std::vector<Vertex>& verticies, 
        std::vector<int> featureVertices, MColorArray* colors, MColor color);
    void createHeatMapMesh(const double& radius, std::unordered_map<Vertex*, float> colorScheme, MColorArray* colors);
    
    MPointArray points;
    MIntArray faceCounts;
    MIntArray faceConnects;

};

#endif