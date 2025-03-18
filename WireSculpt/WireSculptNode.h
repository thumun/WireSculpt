
#ifndef CreateWireSculptNode_H_
#define CreateWireSculptNode_H_

#include <maya/MPxNode.h>
#include <maya/MPointArray.h>

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

protected:
    MObject createMesh(const double& radius, MObject& outData, MStatus& status);
    
    MPointArray points;
    MIntArray faceCounts;
    MIntArray faceConnects;

};

#endif