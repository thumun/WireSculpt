
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
    static MObject  thickness;
    static MObject  outGeom;
    static MObject  fov;
    static MObject  view;
    static MObject  contour;
    static MObject  testSC;
    static MObject  isAbstract;
    static MObject  proximityThresh;
    static MObject  filterThresh;
    static MObject  maxValThresh;

protected:
    MObject createMesh(const double& radius, const double& aAttract, const double& bAttract, 
        const double& aRepel, const double& bRepel, const double& fovVal, const int& viewChoice,
        const int& contourChoice, const double& testSCVal, const double& proximity, const double& filter, 
        const double& maxVal, const bool& isAbstract,
        WireSculptPlugin& ws, const std::string& filePath, std::vector<Vertex>& verticies, 
        std::vector<Edge>& edges, MObject& outData, MStatus& status);

    // map to colors - red channel
    /*void mapToColors(std::unordered_map<Vertex*, float> colorScheme);
    void remapFeatureLengths(std::vector<Edge>& edges, float gamma);*/

    /* CreateMesh helpers for visualization */
    void createWireframeMesh(const double& radius, std::vector<Vertex>& verticies,
        MColorArray* colors, MColor color);
    void createContoursMesh(const double& radius, std::vector<std::pair<vec3f, vec3f>> featureSegments,
        MColorArray* colors, MColor color);

    // visualize corresponding vertices belonging to contours
    void createFeatureVertsMesh(const double& radius, std::vector<Vertex>& verticies, 
        std::vector<int> featureVertices, MColorArray* colors, MColor color);
    
    // visualize heat map on vertices with given colorScheme
    void createHeatMapMesh(const double& radius, std::unordered_map<Vertex*, float> colorScheme, MColorArray* colors);
    
    // visualize feature length warped edge weights as colors -- normalizes distances and maps from range 0-1 to fit rgb
    void createEdgeWeightsMesh(const double& radius, std::vector<Edge>& edges,
        MColorArray* colors);


    MPointArray points;
    MIntArray faceCounts;
    MIntArray faceConnects;

};

#endif