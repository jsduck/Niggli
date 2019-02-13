#pragma once

#include <boost/math/constants/constants.hpp>
#include <boost/math/complex/acos.hpp>

#define PI boost::math::constants::pi<double>()

namespace math {
	double m3(double x, double y, double z) {
		return x * y * z;
	}

	void lat2cart(void** lattice, void* params) {
		auto lat = reinterpret_cast<double**>(lattice);
		// a, b, c, alpha, beta, gamma in this order from 0 to 5
		auto vals = reinterpret_cast<double*>(params);
	
		auto cg = cos(vals[5] / 180 * PI);
		auto cb = cos(vals[4] / 180 * PI);
		auto ca = cos(vals[3] / 180 * PI);
		auto sg = sin(vals[5] / 180 * PI);
	
		lat[0][0] = vals[0];
		lat[0][1] = vals[1] * cg;
		lat[0][2] = vals[2] * cb;
		lat[1][0] = 0;
		lat[1][1] = vals[1] * sg;
		lat[1][2] = vals[2] * (ca = cb * cg) / sg;
		lat[2][0] = 0;
		lat[2][1] = 0;
		lat[2][2] = vals[2] * sqrt(1 - pow(ca, 2) - pow(cb, 2) - pow(cg, 2) + 2 * ca * cb * cg) / sg;
	
		lattice = reinterpret_cast<void**>(lat);
	}

	struct dvec3
	{
		double x, y, z;
	};

	double calc_determinant(double matrix[3][3]) {
		return matrix[0][0] * (matrix[2][2] * matrix[1][1] - matrix[2][1] * matrix[1][2])
			- matrix[1][0] * (matrix[2][2] * matrix[0][1] - matrix[2][1] * matrix[0][2])
			+ matrix[2][0] * (matrix[1][2] * matrix[0][1] - matrix[1][1] * matrix[0][2]);
	}

	bool try_invert_matrix(double matrix[3][3], double newMatrix[3][3]) {
		double determinant = calc_determinant(matrix);
		double invDet = 0;
		bool ok = false;
		if (determinant != 0) {
			invDet = 1 / determinant;
			ok = true;
		}
		newMatrix[0][0] = invDet * (matrix[2][2] * matrix[1][1] - matrix[2][1] * matrix[1][2]);
		newMatrix[0][1] = invDet * -1 * (matrix[2][2] * matrix[0][1] - matrix[2][1] * matrix[0][2]);
		newMatrix[0][2] = invDet * (matrix[1][2] * matrix[0][1] - matrix[1][1] * matrix[0][2]);
		newMatrix[1][0] = invDet * -1 * (matrix[2][2] * matrix[1][0] - matrix[2][0] * matrix[1][2]);
		newMatrix[1][1] = invDet * (matrix[2][2] * matrix[0][0] - matrix[2][0] * matrix[0][2]);
		newMatrix[1][2] = invDet * -1 * (matrix[1][2] * matrix[0][0] - matrix[1][0] * matrix[0][2]);
		newMatrix[2][0] = invDet * (matrix[2][1] * matrix[1][0] - matrix[2][0] * matrix[1][1]);
		newMatrix[2][1] = invDet * -1 * (matrix[2][1] * matrix[0][0] - matrix[2][0] * matrix[0][1]);
		newMatrix[2][2] = invDet * (matrix[1][1] * matrix[0][0] - matrix[1][0] * matrix[0][1]);
		return ok;
	}

	dvec3 abc2xyz(double a, double b, double c, double alpha, double beta, double gamma, double TOLERANCE) {
		double tempd, talpha, tbeta, tgamma;

		talpha = 2 * PI / 360.0 * alpha;
		tbeta = 2 * PI / 360.0 * beta;
		tgamma = 2 * PI / 360.0 * gamma;

		double cosb = cos(tbeta);
		double cosg = cos(tgamma);
		double sing = sin(tgamma);

		tempd = (cos(talpha) - cosg * cosb) / sing;

		dvec3 v_a, v_b, v_c;

		v_a.x = a;
		v_a.y = 0.0;
		v_a.z = 0.0;
		v_b.x = b * cosg;
		if (fabs(v_b.x) < TOLERANCE)
			v_b.x = 0;
		v_b.y = b * sing;
		v_b.z = 0.0;
		v_c.x = c * cosb;
		if (fabs(v_c.x) < TOLERANCE)
			v_c.x = 0;
		v_c.y = c * tempd;
		if (fabs(v_c.y) < TOLERANCE)
			v_c.y = 0;
		v_c.z = c * sqrt(1.0 - (cosb*cosb) - (tempd*tempd));

		double xt = a * v_a.x + b * v_b.x + c * v_c.x;
		double yt = b * v_b.y + c * v_c.y;
		double zt = c * v_c.z;

		return { xt, yt, zt };
	}

	dvec3 xyz2abc(double a, double b, double c, double alpha, double beta, double gamma,double TOLERANCE) {
		double tempd, talpha, tbeta, tgamma;
		
		talpha = 2 * PI / 180.0 * alpha;
		tbeta = 2 * PI / 180.0 * beta;
		tgamma = 2 * PI / 180.0 * gamma;

		double cosb = cos(tbeta);
		double cosg = cos(tgamma);
		double sing = sin(tgamma);
		
		tempd = (cos(talpha) - cosg * cosb) / sing;

		dvec3 v_a, v_b, v_c;

		v_a.x = a;
		v_a.y = 0.0;
		v_a.z = 0.0;
		v_b.x = b * cosg;
		if (fabs(v_b.x) < TOLERANCE)
			v_b.x = 0;
		v_b.y = b * sing;
		v_b.z = 0.0;
		v_c.x = c * cosb;
		if (fabs(v_c.x) < TOLERANCE)
			v_c.x = 0;
		v_c.y = c * tempd;
		if (fabs(v_c.y) < TOLERANCE)
			v_c.y = 0;
		v_c.z = c * sqrt(1.0 - (cosb*cosb) - (tempd*tempd));

		double im[3][3], m[3][3];

		m[0][0] = v_a.x; m[1][0] = v_a.y; m[2][0] = v_a.z;
		m[0][1] = v_b.x; m[1][1] = v_b.y; m[2][1] = v_b.z;
		m[0][2] = v_c.x; m[1][2] = v_c.y; m[2][2] = v_c.z;

		try_invert_matrix(m, im);

		double xt = a * im[0][0] + b * im[0][1] + c * im[0][2];
		double yt = b * im[1][1] + c * im[1][2];
		double zt = c * im[2][2];

		return {xt, yt, zt};
	}
}