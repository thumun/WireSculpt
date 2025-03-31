#include "Contours.h"
#include <maya/MGlobal.h>

using namespace trimesh;
using namespace std;


// Set to false for hardware that has problems with display lists
const bool use_dlists = true;
// Set to false for hardware that has problems with supplying 3D texture coords
const bool use_3dtexc = false;


Contours::Contours(float fovVal, const char* filename) {

	this->fov = fovVal;
	this->feature_size = 0;
	this->featurePoints = {};

	this->themesh = TriMesh::read(filename);
	if (!this->themesh)
		MGlobal::displayInfo("Contours: Error reading file");
	

	this->xffilename = new char[strlen(filename) + 4];
	strcpy(this->xffilename, filename);
	char* dot = strrchr(this->xffilename, '.');
	if (!dot)
		dot = strrchr(this->xffilename, 0);
	strcpy(dot, ".xf");
	
}

// Camera helper functions
MPoint Contours::getCameraPosition()
{
	M3dView view = M3dView::active3dView();
	MDagPath cameraDagPath;
	view.getCamera(cameraDagPath);
	MFnCamera cameraFn(cameraDagPath);

	return cameraFn.eyePoint(MSpace::kWorld);
}

MVector Contours::getCameraDirection()
{
	M3dView view = M3dView::active3dView();
	MDagPath cameraDagPath;
	view.getCamera(cameraDagPath);
	MFnCamera cameraFn(cameraDagPath);

	return cameraFn.viewDirection(MSpace::kWorld);
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
		point fp = p1;
		//MGlobal::displayInfo("p1 is: " + MString() + fp.x + MString() + fp.y + MString() + fp.z);
		featurePoints.push_back({ static_cast<float>(fp.x),
								  static_cast<float>(fp.y),
								  static_cast<float>(fp.z) });
		//featurePoints.push_back(std::vector<float>({ fp.x, fp.y, fp.z }));
		npts++;
	}
	if (z1) {
		float num = (1.0f - z1) * test_num1 + z1 * test_num2;
		float den = (1.0f - z1) * test_den1 + z1 * test_den2;
		/*glColor4f(currcolor[0], currcolor[1], currcolor[2],
			num / (den * fade + num));
		glVertex3fv((1.0f - z1) * p1 + z1 * p2);*/
		point fp = (1.0f - z1) * p1 + z1 * p2;
		//MGlobal::displayInfo("p1 2 is: " + MString() + fp.x + MString() + fp.y + MString() + fp.z);
		featurePoints.push_back({ static_cast<float>(fp.x),
								  static_cast<float>(fp.y),
								  static_cast<float>(fp.z) });
		//featurePoints.push_back(std::vector<float>({ fp.x, fp.y, fp.z }));
		npts++;
	}
	if (z2) {
		float num = (1.0f - z2) * test_num1 + z2 * test_num2;
		float den = (1.0f - z2) * test_den1 + z2 * test_den2;
		/*glColor4f(currcolor[0], currcolor[1], currcolor[2],
			num / (den * fade + num));
		glVertex3fv((1.0f - z2) * p1 + z2 * p2);*/
		point fp = (1.0f - z2) * p1 + z2 * p2;
		//MGlobal::displayInfo("p1 3 is: " + MString() + fp.x + MString() + fp.y + MString() + fp.z);
		featurePoints.push_back({ static_cast<float>(fp.x),
								  static_cast<float>(fp.y),
								  static_cast<float>(fp.z) });
		//featurePoints.push_back(std::vector<float>({ fp.x, fp.y, fp.z }));
		npts++;
	}
	if (npts != 2) {
		/*glColor4f(currcolor[0], currcolor[1], currcolor[2],
			test_num2 / (test_den2 * fade + test_num2));
		glVertex3fv(p2);*/
		point fp = p2;
		//MGlobal::displayInfo("p1 4  is: " + MString() + fp.x + MString() + fp.y + MString() + fp.z);
		featurePoints.push_back({ static_cast<float>(fp.x),
								  static_cast<float>(fp.y),
								  static_cast<float>(fp.z) });
		//featurePoints.push_back(std::vector<float>({ fp.x, fp.y, fp.z }));
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

	//MGlobal::displayInfo("Contours drawface isoline entered");

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

	//MGlobal::displayInfo("Contours drawface isoline exited");

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
}


// Draw the scene
void Contours::redraw()
{
	timestamp t = now();
	viewpos = inv(xf) * point(0, 0, 0);

	//camera.setupGL(xf * themesh->bsphere.center, themesh->bsphere.r);


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
	if (!xf.read(xffilename))
		xf = xform::trans(0, 0, -3.5f / fov * themesh->bsphere.r) *
		xform::trans(-themesh->bsphere.center);
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


