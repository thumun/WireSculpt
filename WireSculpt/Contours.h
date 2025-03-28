#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
//#include "GLCamera.h"
#include "timestamp.h"
//#include "GL/glui.h"
#ifndef DARWIN
//#include <GL/glext.h>
#endif
#include <algorithm>

using namespace trimesh;
using namespace std;

class Contours {
public:
	Contours(float fovVal);

	// Variables
	TriMesh* themesh;
	xform xf;
	float fov = 0.7f;
	char* xffilename; // Filename where we look for "home" position
	point viewpos;    // Current view position
	int draw_c = 1, draw_sc = 1;
	int test_sc = 1;
	float sug_thresh = 0.01;

	// Other miscellaneous variables
	float feature_size;	// Used to make thresholds dimensionless
	vec currcolor;


	void compute_perview(vector<float>& ndotv, vector<float>& kr,
		vector<float>& sctest_num, vector<float>& sctest_den,
		vector<float>& shtest_num, vector<float>& q1,
		vector<vec2>& t1, vector<float>& Dt1q1,
		bool extra_sin2theta = false);

	vec gradkr(int i);

	float find_zero_hermite(int v0, int v1, float val0, float val1,
		const vec& grad0, const vec& grad1);
	void draw_face_isoline2(int v0, int v1, int v2,
		const vector<float>& val,
		const vector<float>& test_num,
		const vector<float>& test_den,
		bool do_hermite, bool do_test, float fade);
	void draw_face_isoline(int v0, int v1, int v2,
		const vector<float>& val,
		const vector<float>& test_num,
		const vector<float>& test_den,
		const vector<float>& ndotv,
		bool do_bfcull, bool do_hermite,
		bool do_test, float fade);

	void draw_isolines(const vector<float>& val,
		const vector<float>& test_num,
		const vector<float>& test_den,
		const vector<float>& ndotv,
		bool do_bfcull, bool do_hermite,
		bool do_test, float fade);
	void draw_mesh();
	void compute_feature_size();
	int main(int argc, char* argv[]);
	void resetview();
	void redraw();
};