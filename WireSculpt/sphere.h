#ifndef sphere_H_
#define sphere_H_

#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MVector.h>
#include <maya/MVectorArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>


class SphereMesh 
{
public:
    SphereMesh(const MPoint& pos, double r = 0.25);
    ~SphereMesh();

    void getMesh(
        MPointArray& points,
        MIntArray& faceCounts,
        MIntArray& faceConnects);

    void appendToMesh(
        MPointArray& points,
        MIntArray& faceCounts,
        MIntArray& faceConnects);

protected:
    void transform(MPointArray& points, MVectorArray& normals);
    MPoint mPos;
    double r;

    // Creates a unit sphere at (0,0,0) with radius r
    static void initSphereMesh(double r);
    static MPointArray gPoints;
    static MVectorArray gNormals;
    static MIntArray gFaceCounts;
    static MIntArray gFaceConnects;
};
#endif