#include "sphere.h"
#include <maya/MMatrix.h>
#include <math.h>

MPointArray SphereMesh::gPoints;
MVectorArray SphereMesh::gNormals;
MIntArray SphereMesh::gFaceCounts;
MIntArray SphereMesh::gFaceConnects;

SphereMesh::SphereMesh(
    const MPoint& pos, double _r) :
    mPos(pos), r(_r)
{
    initSphereMesh(r);
}

SphereMesh::~SphereMesh() {}

void SphereMesh::transform(MPointArray& points, MVectorArray& normals)
{
    for (int i = 0; i < gPoints.length(); i++)
    {
        MPoint p = gPoints[i] + mPos;
        points.append(p);

        MVector n = gNormals[i];
        normals.append(n);
    }
}

void SphereMesh::appendToMesh(
    MPointArray& points,
    MIntArray& faceCounts,
    MIntArray& faceConnects)
{
    MPointArray cpoints;
    MVectorArray cnormals;
    transform(cpoints, cnormals);

    int startIndex = points.length(); // offset for indexes
    for (int i = 0; i < cpoints.length(); i++)
    {
        points.append(cpoints[i]);
    }
    for (int i = 0; i < gFaceCounts.length(); i++)
    {
        faceCounts.append(gFaceCounts[i]);
    }

    for (int i = 0; i < gFaceConnects.length(); i++)
    {
        faceConnects.append(gFaceConnects[i] + startIndex);
    }
}

void SphereMesh::getMesh(
    MPointArray& points,
    MIntArray& faceCounts,
    MIntArray& faceConnects)
{
    MVectorArray cnormals;
    transform(points, cnormals);
    faceCounts = gFaceCounts;
    faceConnects = gFaceConnects;
}

void SphereMesh::initSphereMesh(double r)
{
    int slices = 10;
    int stacks = 10;

    // Add points and normals
    gPoints.clear();
    gNormals.clear();
    gFaceCounts.clear();
    gFaceConnects.clear();

    // Loop through stacks
    for (int i = 0; i <= stacks; i++) {
        float V = (float)i / (float)stacks;
        float phi = V * M_PI;

        // Loop through slices
        for (int j = 0; j <= slices; j++) {
            float U = (float)j / (float)slices;
            float theta = U * (M_PI * 2);

            // use spherical coordinates to calculate the positions.
            float x = r * cos(theta) * sin(phi);
            float y = r * cos(phi);
            float z = r * sin(theta) * sin(phi);

            gPoints.append(MPoint(x, y, z));
            gNormals.append(MVector(x, y, z));  // normals even needed? need to normalize?
        }
    }

    // Set indices for sphere
    for (int i = 0; i < slices * stacks + slices; i++) {
        gFaceCounts.append(3);
        gFaceConnects.append(i);
        gFaceConnects.append(i + slices + 1);
        gFaceConnects.append(i + slices);

        gFaceCounts.append(3);
        gFaceConnects.append(i + slices + 1);
        gFaceConnects.append(i);
        gFaceConnects.append(i + 1);
    }
}
