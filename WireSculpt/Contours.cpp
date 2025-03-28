/*
Authors:
  Szymon Rusinkiewicz, Princeton University
  Doug DeCarlo, Rutgers University

With contributions by:
  Xiaofeng Mi, Rutgers University
  Tilke Judd, MIT

rtsc.cc
Real-time suggestive contours - these days, it also draws many other lines.
*/

//#include <stdio.h>
//#include <stdlib.h>
//#include "TriMesh.h"
//#include "TriMesh_algo.h"
//#include "XForm.h"
////#include "GLCamera.h"
//#include "timestamp.h"
////#include "GL/glui.h"
//#ifndef DARWIN
////#include <GL/glext.h>
//#endif
//#include <algorithm>

#include "Contours.h"
#include <maya/MGlobal.h>
using namespace trimesh;
using namespace std;


// Set to false for hardware that has problems with display lists
const bool use_dlists = true;
// Set to false for hardware that has problems with supplying 3D texture coords
const bool use_3dtexc = false;
/*

// Globals: mesh...
TriMesh* themesh;

// Two cameras: the primary one, and an alternate one to fix the lines
// and see them from a different direction
//GLCamera camera;
xform xf;
float fov = 0.7f;
char* xffilename; // Filename where we look for "home" position
point viewpos;    // Current view position

// Toggles for drawing various lines
int draw_c = 1, draw_sc = 1;

// Toggles for tests we perform
int test_sc = 1;
float sug_thresh = 0.01;

// Toggles for style
//int draw_faded = 0;

//int lighting_style = 0;
////GLUI_Rotation* lightdir = NULL;
//float lightdir_matrix[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
//int light_wrt_camera = true;

// Other miscellaneous variables
float feature_size;	// Used to make thresholds dimensionless
vec currcolor;		// Current line color
*/
//void Contours::setUpContours(const char* filename) {
//	themesh = TriMesh::read(filename);
//	if (!themesh)
//		MGlobal::displayInfo("Contours: Error reading file");
//
//	xffilename = new char[strlen(filename) + 4];
//	strcpy(xffilename, filename);
//	char* dot = strrchr(xffilename, '.');
//	if (!dot)
//		dot = strrchr(xffilename, 0);
//	strcpy(dot, ".xf");
//
//	themesh->need_tstrips();
//	themesh->need_bsphere();
//	themesh->need_normals();
//	themesh->need_curvatures();
//	themesh->need_dcurv();
//	compute_feature_size();
//
//	redraw();
//
//	resetview();
//}

Contours::Contours(float fovVal) {
	this->fov = fovVal;
}

// Compute per-vertex n dot l, n dot v, radial curvature, and
// derivative of curvature for the current view
void Contours::compute_perview(vector<float>& ndotv, vector<float>& kr,
	vector<float>& sctest_num, vector<float>& sctest_den,
	vector<float>& shtest_num, vector<float>& q1,
	vector<vec2>& t1, vector<float>& Dt1q1,
	bool extra_sin2theta)
{

	int nv = themesh->vertices.size();

	float scthresh = sug_thresh / sqr(feature_size);
	bool need_DwKr = draw_sc;

	ndotv.resize(nv);
	kr.resize(nv);

	if (need_DwKr) {
		sctest_num.resize(nv);
		sctest_den.resize(nv);
	}

	// Compute quantities at each vertex
#pragma omp parallel for
	for (int i = 0; i < nv; i++) {
		// Compute n DOT v
		vec viewdir = viewpos - themesh->vertices[i];
		float rlv = 1.0f / len(viewdir);
		viewdir *= rlv;
		ndotv[i] = viewdir DOT themesh->normals[i];

		float u = viewdir DOT themesh->pdir1[i], u2 = u * u;
		float v = viewdir DOT themesh->pdir2[i], v2 = v * v;

		// Note:  this is actually Kr * sin^2 theta
		kr[i] = themesh->curv1[i] * u2 + themesh->curv2[i] * v2;

		if (!need_DwKr)
			continue;

		// Use DwKr * sin(theta) / cos(theta) for cutoff test
		sctest_num[i] = u2 * (u * themesh->dcurv[i][0] +
			3.0f * v * themesh->dcurv[i][1]) +
			v2 * (3.0f * u * themesh->dcurv[i][2] +
				v * themesh->dcurv[i][3]);
		float csc2theta = 1.0f / (u2 + v2);
		sctest_num[i] *= csc2theta;
		float tr = (themesh->curv2[i] - themesh->curv1[i]) *
			u * v * csc2theta;
		sctest_num[i] -= 2.0f * ndotv[i] * sqr(tr);
		if (extra_sin2theta)
			sctest_num[i] *= u2 + v2;

		sctest_den[i] = ndotv[i];

		sctest_num[i] -= scthresh * sctest_den[i];
	}
}


// Compute gradient of (kr * sin^2 theta) at vertex i
vec Contours::gradkr(int i)
{
	vec viewdir = viewpos - themesh->vertices[i];
	float rlen_viewdir = 1.0f / len(viewdir);
	viewdir *= rlen_viewdir;

	float ndotv = viewdir DOT themesh->normals[i];
	float sintheta = sqrt(1.0f - sqr(ndotv));
	float csctheta = 1.0f / sintheta;
	float u = (viewdir DOT themesh->pdir1[i]) * csctheta;
	float v = (viewdir DOT themesh->pdir2[i]) * csctheta;
	float kr = themesh->curv1[i] * u * u + themesh->curv2[i] * v * v;
	float tr = u * v * (themesh->curv2[i] - themesh->curv1[i]);
	float kt = themesh->curv1[i] * (1.0f - u * u) +
		themesh->curv2[i] * (1.0f - v * v);
	vec w = u * themesh->pdir1[i] + v * themesh->pdir2[i];
	vec wperp = u * themesh->pdir2[i] - v * themesh->pdir1[i];
	const Vec<4>& C = themesh->dcurv[i];

	vec g = themesh->pdir1[i] * (u * u * C[0] + 2.0f * u * v * C[1] + v * v * C[2]) +
		themesh->pdir2[i] * (u * u * C[1] + 2.0f * u * v * C[2] + v * v * C[3]) -
		2.0f * csctheta * tr * (rlen_viewdir * wperp +
			ndotv * (tr * w + kt * wperp));
	g *= (1.0f - sqr(ndotv));
	g -= 2.0f * kr * sintheta * ndotv * (kr * w + tr * wperp);
	return g;
}


// Find a zero crossing between val0 and val1 by linear interpolation
// Returns 0 if zero crossing is at val0, 1 if at val1, etc.
static inline float find_zero_linear(float val0, float val1)
{
	return val0 / (val0 - val1);
}


// Find a zero crossing using Hermite interpolation
float Contours::find_zero_hermite(int v0, int v1, float val0, float val1,
	const vec& grad0, const vec& grad1)
{
	if (unlikely(val0 == val1))
		return 0.5f;

	// Find derivatives along edge (of interpolation parameter in [0,1]
	// which means that e01 doesn't get normalized)
	vec e01 = themesh->vertices[v1] - themesh->vertices[v0];
	float d0 = e01 DOT grad0, d1 = e01 DOT grad1;

	// This next line would reduce val to linear interpolation
	//d0 = d1 = (val1 - val0);

	// Use hermite interpolation:
	//   val(s) = h1(s)*val0 + h2(s)*val1 + h3(s)*d0 + h4(s)*d1
	// where
	//  h1(s) = 2*s^3 - 3*s^2 + 1
	//  h2(s) = 3*s^2 - 2*s^3
	//  h3(s) = s^3 - 2*s^2 + s
	//  h4(s) = s^3 - s^2
	//
	//  val(s)  = [2(val0-val1) +d0+d1]*s^3 +
	//            [3(val1-val0)-2d0-d1]*s^2 + d0*s + val0
	// where
	//
	//  val(0) = val0; val(1) = val1; val'(0) = d0; val'(1) = d1
	//

	// Coeffs of cubic a*s^3 + b*s^2 + c*s + d
	float a = 2 * (val0 - val1) + d0 + d1;
	float b = 3 * (val1 - val0) - 2 * d0 - d1;
	float c = d0, d = val0;

	// -- Find a root by bisection
	// (as Newton can wander out of desired interval)

	// Start with entire [0,1] interval
	float sl = 0.0f, sr = 1.0f, valsl = val0, valsr = val1;

	// Check if we're in a (somewhat uncommon) 3-root situation, and pick
	// the middle root if it happens (given we aren't drawing curvy lines,
	// seems the best approach..)
	//
	// Find extrema of derivative (a -> 3a; b -> 2b, c -> c),
	// and check if they're both in [0,1] and have different signs
	float disc = 4 * b - 12 * a * c;
	if (disc > 0 && a != 0) {
		disc = sqrt(disc);
		float r1 = (-2 * b + disc) / (6 * a);
		float r2 = (-2 * b - disc) / (6 * a);
		if (r1 >= 0 && r1 <= 1 && r2 >= 0 && r2 <= 1) {
			float vr1 = (((a * r1 + b) * r1 + c) * r1) + d;
			float vr2 = (((a * r2 + b) * r2 + c) * r2) + d;
			// When extrema have different signs inside an
			// interval with endpoints with different signs,
			// the middle root is in between the two extrema
			if ((vr1 < 0.0f && vr2 >= 0.0f) ||
				(vr1 > 0.0f && vr2 <= 0.0f)) {
				// 3 roots
				if (r1 < r2) {
					sl = r1;
					valsl = vr1;
					sr = r2;
					valsr = vr2;
				}
				else {
					sl = r2;
					valsl = vr2;
					sr = r1;
					valsr = vr1;
				}
			}
		}
	}

	// Bisection method (constant number of interations)
	for (int iter = 0; iter < 10; iter++) {
		float sbi = (sl + sr) / 2.0f;
		float valsbi = (((a * sbi + b) * sbi) + c) * sbi + d;

		// Keep the half which has different signs
		if ((valsl < 0.0f && valsbi >= 0.0f) ||
			(valsl > 0.0f && valsbi <= 0.0f)) {
			sr = sbi;
			valsr = valsbi;
		}
		else {
			sl = sbi;
			valsl = valsbi;
		}
	}

	return 0.5f * (sl + sr);
}


// Draw part of a zero-crossing curve on one triangle face, but only if
// "test_num/test_den" is positive.  v0,v1,v2 are the indices of the 3
// vertices, "val" are the values of the scalar field whose zero
// crossings we are finding, and "test_*" are the values we are testing
// to make sure they are positive.  This function assumes that val0 has
// opposite sign from val1 and val2 - the following function is the
// general one that figures out which one actually has the different sign.
void Contours::draw_face_isoline2(int v0, int v1, int v2,
	const vector<float>& val,
	const vector<float>& test_num,
	const vector<float>& test_den,
	bool do_hermite, bool do_test, float fade)
{
	// How far along each edge?
	float w10 = do_hermite ?
		find_zero_hermite(v0, v1, val[v0], val[v1],
			gradkr(v0), gradkr(v1)) :
		find_zero_linear(val[v0], val[v1]);
	float w01 = 1.0f - w10;
	float w20 = do_hermite ?
		find_zero_hermite(v0, v2, val[v0], val[v2],
			gradkr(v0), gradkr(v2)) :
		find_zero_linear(val[v0], val[v2]);
	float w02 = 1.0f - w20;

	// Points along edges
	point p1 = w01 * themesh->vertices[v0] + w10 * themesh->vertices[v1];
	point p2 = w02 * themesh->vertices[v0] + w20 * themesh->vertices[v2];

	float test_num1 = 1.0f, test_num2 = 1.0f;
	float test_den1 = 1.0f, test_den2 = 1.0f;
	float z1 = 0.0f, z2 = 0.0f;
	bool valid1 = true;
	if (do_test) {
		// Interpolate to find value of test at p1, p2
		test_num1 = w01 * test_num[v0] + w10 * test_num[v1];
		test_num2 = w02 * test_num[v0] + w20 * test_num[v2];
		if (!test_den.empty()) {
			test_den1 = w01 * test_den[v0] + w10 * test_den[v1];
			test_den2 = w02 * test_den[v0] + w20 * test_den[v2];
		}
		// First point is valid iff num1/den1 is positive,
		// i.e. the num and den have the same sign
		valid1 = ((test_num1 >= 0.0f) == (test_den1 >= 0.0f));
		// There are two possible zero crossings of the test,
		// corresponding to zeros of the num and den
		if ((test_num1 >= 0.0f) != (test_num2 >= 0.0f))
			z1 = test_num1 / (test_num1 - test_num2);
		if ((test_den1 >= 0.0f) != (test_den2 >= 0.0f))
			z2 = test_den1 / (test_den1 - test_den2);
		// Sort and order the zero crossings
		if (z1 == 0.0f)
			z1 = z2, z2 = 0.0f;
		else if (z2 < z1)
			swap(z1, z2);
	}

	// If the beginning of the segment was not valid, and
	// no zero crossings, then whole segment invalid
	if (!valid1 && !z1 && !z2)
		return;

	// Draw the valid piece(s)
	int npts = 0;
	if (valid1) {
		/*glColor4f(currcolor[0], currcolor[1], currcolor[2],
			test_num1 / (test_den1 * fade + test_num1));
		glVertex3fv(p1);*/
		npts++;
	}
	if (z1) {
		float num = (1.0f - z1) * test_num1 + z1 * test_num2;
		float den = (1.0f - z1) * test_den1 + z1 * test_den2;
		/*glColor4f(currcolor[0], currcolor[1], currcolor[2],
			num / (den * fade + num));
		glVertex3fv((1.0f - z1) * p1 + z1 * p2);*/
		npts++;
	}
	if (z2) {
		float num = (1.0f - z2) * test_num1 + z2 * test_num2;
		float den = (1.0f - z2) * test_den1 + z2 * test_den2;
		/*glColor4f(currcolor[0], currcolor[1], currcolor[2],
			num / (den * fade + num));
		glVertex3fv((1.0f - z2) * p1 + z2 * p2);*/
		npts++;
	}
	if (npts != 2) {
		/*glColor4f(currcolor[0], currcolor[1], currcolor[2],
			test_num2 / (test_den2 * fade + test_num2));
		glVertex3fv(p2);*/
	}

}


// See above.  This is the driver function that figures out which of
// v0, v1, v2 has a different sign from the others.
void Contours::draw_face_isoline(int v0, int v1, int v2,
	const vector<float>& val,
	const vector<float>& test_num,
	const vector<float>& test_den,
	const vector<float>& ndotv,
	bool do_bfcull, bool do_hermite,
	bool do_test, float fade)
{
	// Backface culling
	if (likely(do_bfcull && ndotv[v0] <= 0.0f &&
		ndotv[v1] <= 0.0f && ndotv[v2] <= 0.0f))
		return;

	// Quick reject if derivs are negative
	if (do_test) {
		if (test_den.empty()) {
			if (test_num[v0] <= 0.0f &&
				test_num[v1] <= 0.0f &&
				test_num[v2] <= 0.0f)
				return;
		}
		else {
			if (test_num[v0] <= 0.0f && test_den[v0] >= 0.0f &&
				test_num[v1] <= 0.0f && test_den[v1] >= 0.0f &&
				test_num[v2] <= 0.0f && test_den[v2] >= 0.0f)
				return;
			if (test_num[v0] >= 0.0f && test_den[v0] <= 0.0f &&
				test_num[v1] >= 0.0f && test_den[v1] <= 0.0f &&
				test_num[v2] >= 0.0f && test_den[v2] <= 0.0f)
				return;
		}
	}

	// Figure out which val has different sign, and draw
	if ((val[v0] < 0.0f && val[v1] >= 0.0f && val[v2] >= 0.0f) ||
		(val[v0] > 0.0f && val[v1] <= 0.0f && val[v2] <= 0.0f))
		draw_face_isoline2(v0, v1, v2,
			val, test_num, test_den,
			do_hermite, do_test, fade);
	else if ((val[v1] < 0.0f && val[v2] >= 0.0f && val[v0] >= 0.0f) ||
		(val[v1] > 0.0f && val[v2] <= 0.0f && val[v0] <= 0.0f))
		draw_face_isoline2(v1, v2, v0,
			val, test_num, test_den,
			do_hermite, do_test, fade);
	else if ((val[v2] < 0.0f && val[v0] >= 0.0f && val[v1] >= 0.0f) ||
		(val[v2] > 0.0f && val[v0] <= 0.0f && val[v1] <= 0.0f))
		draw_face_isoline2(v2, v0, v1,
			val, test_num, test_den,
			do_hermite, do_test, fade);
}


// Takes a scalar field and renders the zero crossings, but only where
// test_num/test_den is greater than 0.
void Contours::draw_isolines(const vector<float>& val,
	const vector<float>& test_num,
	const vector<float>& test_den,
	const vector<float>& ndotv,
	bool do_bfcull, bool do_hermite,
	bool do_test, float fade)
{
	const int* t = &themesh->tstrips[0];
	const int* stripend = t;
	const int* end = t + themesh->tstrips.size();

	// Walk through triangle strips
	while (1) {
		if (unlikely(t >= stripend)) {
			if (unlikely(t >= end))
				return;
			// New strip: each strip is stored as
			// length followed by indices
			stripend = t + 1 + *t;
			// Skip over length plus first two indices of
			// first face
			t += 3;
		}
		// Draw a line if, among the values in this triangle,
		// at least one is positive and one is negative
		const float& v0 = val[*t], & v1 = val[*(t - 1)], & v2 = val[*(t - 2)];
		if (unlikely((v0 > 0.0f || v1 > 0.0f || v2 > 0.0f) &&
			(v0 < 0.0f || v1 < 0.0f || v2 < 0.0f)))
			draw_face_isoline(*(t - 2), *(t - 1), *t,
				val, test_num, test_den, ndotv,
				do_bfcull, do_hermite, do_test, fade);
		t++;
	}
}

// Draw the mesh, possibly including a bunch of lines
void Contours::draw_mesh()
{
	// These are static so the memory isn't reallocated on every frame
	static vector<float> ndotv, kr;
	static vector<float> sctest_num, sctest_den, shtest_num;
	static vector<float> q1, Dt1q1;
	static vector<vec2> t1;
	compute_perview(ndotv, kr, sctest_num, sctest_den, shtest_num,
		q1, t1, Dt1q1);
	int nv = themesh->vertices.size();

	// Enable antialiased lines
	/*glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);*/

	// Kr = 0 loops
	if (draw_sc && !test_sc) {
		currcolor = vec(0.6, 0.6, 0.6);
		/*glLineWidth(1.5);
		glBegin(GL_LINES);*/
		draw_isolines(kr, sctest_num, sctest_den, ndotv,
			true, false, false, 0.0f);
		/*glEnd();*/
		currcolor = vec(0.0, 0.0, 0.0);
	}

	// Suggestive contours and contours
	if (draw_sc) {
		float fade = 0.0f;
		/*glLineWidth(2.5);
		glBegin(GL_LINES);*/
		draw_isolines(kr, sctest_num, sctest_den, ndotv,
			true, false, true, fade);
		/*glEnd();*/
	}
	if (draw_c) {
		/*glLineWidth(2.5);
		glBegin(GL_LINES);*/
		draw_isolines(ndotv, kr, vector<float>(), ndotv,
			false, false, true, 0.0f);
		/*glEnd();*/
	}

	/*glDisable(GL_LINE_SMOOTH);
	glDisable(GL_POINT_SMOOTH);
	glDisable(GL_BLEND);
	glDepthMask(GL_TRUE);*/
}


// Signal a redraw
void need_redraw()
{
	/*glutPostRedisplay();*/
}


// Clear the screen and reset OpenGL modes to something sane
void cls()
{
	/*glDisable(GL_DITHER);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	glClearColor(1, 1, 1, 0);
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);*/
}


// Set up viewport and scissoring for the subwindow, and optionally draw
// a box around it (actually, just clears a rectangle one pixel bigger
// to black).  Assumes current viewport is set up for the whole window.
void set_subwindow_viewport(bool draw_box = false)
{
	/*GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	GLint x = V[0], y = V[1], w = V[2], h = V[3];*/
	/*int boxsize = min(w, h) / 3;

	x += w - boxsize * 11 / 10;
	y += h - boxsize * 11 / 10;
	w = h = boxsize;*/

	if (draw_box) {
		/*glViewport(x - 1, y - 1, w + 2, h + 2);
		glScissor(x - 1, y - 1, w + 2, h + 2);
		glClearColor(0, 0, 0, 0);
		glEnable(GL_SCISSOR_TEST);
		glClear(GL_COLOR_BUFFER_BIT);
		glScissor(x, y, w, h);*/
	}

	/*glViewport(x, y, w, h);*/
}


// Draw the scene
void Contours::redraw()
{
	timestamp t = now();
	viewpos = inv(xf) * point(0, 0, 0);
	/*GLUI_Master.auto_set_viewport();*/

	//camera.setupGL(xf * themesh->bsphere.center, themesh->bsphere.r);

	cls();

	// Transform and draw
	/*glPushMatrix();
	glMultMatrixd((double*)xf);*/
	draw_mesh();
	/*glPopMatrix();*/

	/*glDisable(GL_SCISSOR_TEST);
	glutSwapBuffers();*/
	printf("\rElapsed time: %.2f msec.", 1000.0f * (now() - t));
	fflush(stdout);

	// See if we need to autospin the camera(s)
	/*if (camera.autospin(xf))
		need_redraw();*/
}


// Set the view to look at the middle of the mesh, from reasonably far away
void Contours::resetview()
{
	//camera.stopspin();

	if (!xf.read(xffilename))
		xf = xform::trans(0, 0, -3.5f / fov * themesh->bsphere.r) *
		xform::trans(-themesh->bsphere.center);

	// Reset light position too
	//lightdir->reset();
}


// Compute a "feature size" for the mesh: computed as 1% of
// the reciprocal of the 10-th percentile curvature
void Contours::compute_feature_size()
{
	int nv = themesh->curv1.size();
	int nsamp = min(nv, 500);

	vector<float> samples;
	samples.reserve(nsamp * 2);

	for (int i = 0; i < nsamp; i++) {
		// Quick 'n dirty portable random number generator
		static unsigned randq = 0;
		randq = unsigned(1664525) * randq + unsigned(1013904223);

		int ind = randq % nv;
		samples.push_back(fabs(themesh->curv1[ind]));
		samples.push_back(fabs(themesh->curv2[ind]));
	}

	const float frac = 0.1f;
	const float mult = 0.01f;
	themesh->need_bsphere();
	float max_feature_size = 0.05f * themesh->bsphere.r;

	int which = int(frac * samples.size());
	nth_element(samples.begin(), samples.begin() + which, samples.end());

	feature_size = min(mult / samples[which], max_feature_size);
}


// Handle mouse button and motion events
static unsigned buttonstate = 0;
static const unsigned ctrl_pressed = 1 << 30;

void mousemotionfunc(int x, int y)
{
	// Ctrl+mouse = relight
	if (buttonstate & ctrl_pressed) {
		/*GLUI_Master.auto_set_viewport();
		GLint V[4];
		glGetIntegerv(GL_VIEWPORT, V);*/
		//y = V[1] + V[3] - 1 - y; // Adjust for top-left vs. bottom-left
		//float xx = 2.0f * float(x - V[0]) / float(V[2]) - 1.0f;
		//float yy = 2.0f * float(y - V[1]) / float(V[3]) - 1.0f;
		//float theta = M_PI * min(sqrt(xx * xx + yy * yy), 1.0f);
		//float phi = atan2(yy, xx);
		//XForm<float> lightxf = lightxf.rot(phi, 0, 0, 1) *
		//	lightxf.rot(theta, 0, 1, 0);
		/*lightdir->set_float_array_val((float*)lightxf);*/
		need_redraw();
		return;
	}

	/*static const Mouse::button physical_to_logical_map[] = {
		Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
		Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
	};
	Mouse::button b = Mouse::NONE;
	if (buttonstate & (1 << 3))
		b = Mouse::WHEELUP;
	else if (buttonstate & (1 << 4))
		b = Mouse::WHEELDOWN;
	else
		b = physical_to_logical_map[buttonstate & 7];

	camera.mouse(x, y, b,
		xf * themesh->bsphere.center,
		themesh->bsphere.r, xf);*/

	need_redraw();
	/*GLUI_Master.sync_live_all();*/
}

void mousebuttonfunc(int button, int state, int x, int y)
{
	/*if (state == GLUT_DOWN)
		buttonstate |= (1 << button);
	else
		buttonstate &= ~(1 << button);

	if (glutGetModifiers() & GLUT_ACTIVE_CTRL)
		buttonstate |= ctrl_pressed;
	else
		buttonstate &= ~ctrl_pressed;

	mousemotionfunc(x, y);*/
}


#define Ctrl (1-'a')

// Reshape the window.  We clear the window here to possibly avoid some
// weird problems.  Yuck.
void reshape(int x, int y)
{
	/*GLUI_Master.auto_set_viewport();
	cls();
	glutSwapBuffers();
	need_redraw();*/
}


//void usage(const char* myname)
//{
//	fprintf(stderr, "Usage: %s [-options] infile\n", myname);
//	exit(1);
//}


int Contours::main(int argc, char* argv[])
{

	int wwid = 820, wht = 700;
	/*for (int j = 1; j < argc; j++) {
		if (argv[j][0] == '+') {
			sscanf(argv[j] + 1, "%d,%d,%f,%f", &wwid, &wht,
				&sug_thresh, &ph_thresh);
		}
	}*/

	/*glutInitWindowSize(wwid, wht);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInit(&argc, argv);*/

	//if (argc < 2)
		//usage(argv[0]);

	// Skip over any parameter beginning with a '-' or '+'
	int i = 1;
	while (i < argc - 1 && argv[i][0] == '-' && argv[i][0] == '+') {
		i++;
		if (!strcmp(argv[i - 1], "--"))
			break;
	}
	const char* filename = argv[i];

	themesh = TriMesh::read(filename);
	//if (!themesh)
		//usage(argv[0]);

	xffilename = new char[strlen(filename) + 4];
	strcpy(xffilename, filename);
	char* dot = strrchr(xffilename, '.');
	if (!dot)
		dot = strrchr(xffilename, 0);
	strcpy(dot, ".xf");

	themesh->need_tstrips();
	themesh->need_bsphere();
	themesh->need_normals();
	themesh->need_curvatures();
	themesh->need_dcurv();
	compute_feature_size();
	//currsmooth = 0.5f * themesh->feature_size();

	char windowname[255];
	sprintf(windowname, "RTSC - %s", filename);
	//int main_win = glutCreateWindow(windowname);

	//glutDisplayFunc(redraw);
	//GLUI_Master.set_glutMouseFunc(mousebuttonfunc);
	//glutMotionFunc(mousemotionfunc);
	//GLUI_Master.set_glutReshapeFunc(reshape);

	//GLUI* glui = GLUI_Master.create_glui_subwindow(main_win, GLUI_SUBWINDOW_BOTTOM);
	//glui->set_main_gfx_window(main_win);
	//GLUI_Rollout* g = glui->add_rollout("Options", false);
	//glui->add_statictext_to_panel(g, "Lines:");
	//glui->add_checkbox_to_panel(g, "Occluding contours", &draw_c);
	//glui->add_checkbox_to_panel(g, "Suggestive contours", &draw_sc);

	//glui->add_column_to_panel(g, false);
	//glui->add_statictext_to_panel(g, "Line tests:");
	//glui->add_checkbox_to_panel(g, "Trim SC", &test_sc);
	//glui->add_slider_to_panel(g, "SC thresh", GLUI_SLIDER_FLOAT,
	//	0.0, 0.1, &sug_thresh);

	//glui->add_statictext_to_panel(g, " ");
	//glui->add_statictext_to_panel(g, "Mesh style:");
	//GLUI_RadioGroup* r = glui->add_radiogroup_to_panel(g, 0 /*&color_style*/);
	//glui->add_radiobutton_to_group(r, "White");
	//if (!themesh->colors.empty())
	//	glui->add_radiobutton_to_group(r, "Mesh colors");

	//glui->add_column_to_panel(g, false);
	//glui->add_statictext_to_panel(g, "Lighting:");
	//r = glui->add_radiogroup_to_panel(g, &lighting_style);
	//glui->add_radiobutton_to_group(r, "None");
	//lightdir = glui->add_rotation_to_panel(g, "Direction",
	//	(float*)&lightdir_matrix);
	//glui->add_checkbox_to_panel(g, "On camera", &light_wrt_camera);
	//lightdir->reset();

	//glui->add_column_to_panel(g, false);
	//glui->add_button_to_panel(g, "Exit", 0, exit);

	//resetview();

	//glutMainLoop();
	return 0;

}
