#pragma once

#include "math.h"

#include "unit_cell.h"

#include <boost/units/cmath.hpp>
#include <boost/math/special_functions.hpp>

struct niggli
{
	bool find_point(unit_cell uc) {
		using namespace boost::math;
		using namespace math;

		int br = 0;

		auto a = uc.a; auto b = uc.b; auto c = uc.c;
		auto alpha = uc.alpha; auto beta = uc.beta; auto gamma = uc.gamma;

		auto ca = cos(uc.alpha / 180 * PI);
		auto cb = cos(uc.beta / 180 * PI);
		auto cg = cos(uc.gamma / 180 * PI);

		// STEP 1
		step1:
		if (br > 50) {
			printf("NIGGLI POINT NOT FOUND\n");
			return false;
		}
		br++;
		if (abs(a) > abs(b) || (abs(a) == abs(b) && abs(m3(b, c, ca) > abs(m3(a, c, cb))))) {
			std::swap(a, b);
		}
		// STEP 2
		if (abs(b) > abs(c) || (abs(b) == abs(c) && abs(m3(a, c, cb) > abs(m3(a, b, cg))))) {
			std::swap(b, c);
			goto step1;
		}
		// STEP 3
		step3:
		if (!((m3(a, b, cg) > 0 && m3(a, c, cb) > 0 && m3(b, c, ca) > 0) || 
			(m3(a, b, cg) <= 0 && m3(a, c, cb) <= 0 && m3(b, c, ca) <= 0))) {
			a = -a;
		}
		// STEP 4
		if (!((m3(a, b, cg) > 0 && m3(a, c, cb) > 0 && m3(b, c, ca) > 0) ||
			(m3(a, b, cg) <= 0 && m3(a, c, cb) <= 0 && m3(b, c, ca) <= 0))) {
			b = -b;
			goto step3;
		}
		// STEP 5
		if (2 * abs(m3(b, c, ca)) > pow(b, 2) ||
			(2 * m3(b, c, ca) == pow(b, 2) && 2 * m3(a, c, cb) < m3(a, b, cg)) || 
			(2 * m3(b, c, ca) == -pow(b, 2) && m3(a, b, cg) < 0)) {
			c = c - sign(m3(b, c, ca)) * b;
			goto step1;
		}
		// STEP 6
		if (2 * abs(m3(a, c, cb)) > pow(a, 2) ||
			(2 * m3(a, c, cb) == pow(a, 2) && 2 * m3(b, c, ca) < m3(a, b, cg)) ||
			(2 * m3(a, c, cb) == -pow(a, 2) && m3(a, b, cg) < 0)) {
			c = c - sign(m3(a, c, ca)) * a;
			goto step1;
		}
		// STEP 7
		if (2 * abs(m3(a, b, cg)) > pow(c, 2) ||
			(2 * m3(a, b, cg) == pow(c, 2) && 2 * m3(b, c, ca) < m3(a, c, cb)) ||
			(2 * m3(a, b, cg) == -pow(c, 2) && m3(a, c, cb) < 0)) {
			b = b - sign(m3(a, b, cg)) * a;
			goto step1;
		}
		// STEP 8
		auto K = pow(a, 2) + pow(b, 2) + 2 * m3(a, b, cg) + 2 * m3(a, c, cb) + 2 * m3(b, c, ca);
		// STEP 9
		if (K < 0 || (K == 0 && pow(a, 2) + 2 * m3(a, c, cb) > pow(b, 2) + 2 * m3(b, c, ca))) {
			c = a + b + c;
			goto step1;
		}
		// STEP 10
		auto u = pow(a, 2) / pow(c, 2);
		auto v = pow(b, 2) / pow(c, 2);
		auto x = 2 * b * c * ca / pow(uc.c, 2);
		auto y = 2 * a * c * cb / pow(uc.c, 2);
		auto z = 2 * a * b * cg / pow(uc.c, 2);

		printf("NIGGLI POINT FOUND :\n u: %.4f v: %.4f\n x: %.4f y: %.4f z: %.4f\n", u, v, x, y, z);
		printf(" a: %.4f b: %.4f c: %.4f\n", a, b, c);

		return true;
	}

	bool verify_point(unit_cell uc) {
		using namespace boost::math;

		auto u = pow(uc.a, 2) / pow(uc.c, 2);
		auto v = pow(uc.b, 2) / pow(uc.c, 2);

		auto x = 2 * uc.b * uc.c * cos(uc.alpha / 180 * PI) / pow(uc.c, 2);
		auto y = 2 * uc.a * uc.c * cos(uc.beta / 180 * PI) / pow(uc.c, 2);
		auto z = 2 * uc.a * uc.b * cos(uc.gamma / 180 * PI) / pow(uc.c, 2);

		if (!(u <=v && v <= 1))
			return false;

		if (!((x > 0 && y > 0 && z > 0) || (x <= 0 && y <= 0 && z <= 0)))
			return false;

		if (u == v && !(abs(x) <= abs(y)))
			return false;

		printf("NIGGLI POINT (v):\n u: %.4f v: %.4f\n x: %.4f y: %.4f z: %.4f\n", u, v, x, y, z);

		return true;
	}

	bool reduce(unit_cell uc, double epsilon) {
		using namespace boost::math;

		int br = 0;

		auto a = uc.a; auto b = uc.b; auto c = uc.c;
		auto alpha = uc.alpha; auto beta = uc.beta; auto gamma = uc.gamma;

		//printf("INPUT:\n a: %.4f b: %.4f c: %.4f \n alpha: %.4f beta: %.4f gamma: %.4f\n", a, b, c, alpha, beta, gamma);
		//find_point({ a, b, c, alpha, beta, gamma });

		double A = pow(a, 2);
		double B = pow(b, 2);
		double C = pow(c, 2);

		//boost::math::cos

		double xi = 2 * b * c * cos_pi(alpha / 180.f);
		double eta = 2 * a * c * cos_pi(beta / 180.f);
		double zeta = 2 * a * b * cos_pi(gamma / 180.f);

		// TEST VALUES
		// A = 9; B = 27; C = 4;
		// xi = -5; eta = -4; zeta = -22;

		// Implementation of the "Unified Algorithm for Determining the Reduced (Niggli) Cell" by I. KRIVY and B. GRUBER - 1975
		// STEP 1
		step1:
		if (br > 50) {
			//printf("NIGGLI CELL NOT REDUCED\n");
			return false;
		}
		br++;
		if (A > B + epsilon || (!(abs(A - B) > epsilon) && abs(xi) > abs(eta) + epsilon)) {
			std::swap(A, B);
			std::swap(xi, eta);
			//br = 0;
		}
		// STEP 2
		if (B > C + epsilon || (!(abs(B - C) > epsilon) && abs(eta) > abs(zeta) + epsilon)) {
			std::swap(B, C);
			std::swap(eta, zeta);
			goto step1;
		}
		// STEP 3
		if (xi * eta * zeta > 0) {
			xi = abs(xi);
			eta = abs(eta);
			zeta = abs(zeta);
			//br = 0;
		}
		// STEP 4
		if (xi * eta * zeta <= 0) {
			xi = -abs(xi);
			eta = -abs(eta);
			zeta = -abs(zeta);
			//br = 0;
		}
		// STEP 5
		if (abs(xi) > B + epsilon || (!(abs(B - xi) > epsilon) && 2 * eta < zeta - epsilon) || (!(abs(B + xi) > epsilon) && zeta < -epsilon)) {
			C = B + C - xi * sign(xi);
			eta = eta - zeta * sign(xi);
			xi = xi - 2 * B * sign(xi);
			goto step1;
		}
		// STEP 6
		if (abs(eta) > A + epsilon || (!(abs(A - eta) > epsilon) && 2 * xi < zeta - epsilon) || (!(abs(A + eta) > epsilon) && zeta < -epsilon)) {
			C = A + C - eta * sign(eta);
			xi = xi - zeta * sign(eta);
			eta = eta - 2 * A * sign(eta);
			goto step1;
		}
		// STEP 7
		if (abs(zeta) > A + epsilon || (!(abs(A - zeta) > epsilon) && 2 * xi < eta - epsilon) || (!(abs(A + zeta) > epsilon) && eta < -epsilon)) {
			B = A + B - zeta * sign(zeta);
			xi = xi - eta * sign(zeta);
			zeta = zeta - 2 * A * sign(zeta);
			goto step1;
		}
		// STEP 8
		if (xi + eta + zeta + A + B < -epsilon || (!(abs(xi + eta + zeta + A + B) > epsilon) && 2 * (A + eta) + zeta > epsilon)) {
			C = A + B + C + xi + eta + zeta;
			xi = 2 * B + xi + zeta;
			eta = 2 * A + eta + zeta;
			goto step1;
		}

		//printf("NIGGLI CELL REDUCED\n");
		
		// Divide xi, eta and zeta by 2 to obtain Niggli form
		xi /= 2;
		eta /= 2;
		zeta /= 2;

		//printf("RAW:\n A: %.4f B: %.4f C: %.4f \n xi: %.4f eta: %.4f zeta: %.4f\n", A, B, C, xi, eta, zeta);

		//find_point({ A, B, C, xi, eta, zeta });
		//verify_point({ A, B, C, xi, eta, zeta });

		// Normalise the results to describe the Niggli cell in the usual way

		A = sqrt(A);
		B = sqrt(B);
		C = sqrt(C);

		xi = acos(xi / (B * C)) / PI * 180;
		eta = acos(eta / (A * C)) / PI * 180;
		zeta = acos(zeta / (A * B)) / PI * 180;

		//printf("NORM:\n a': %.4f b': %.4f c': %.4f \n alpha': %.4f beta': %.4f gamma': %.4f\n", A, B, C, xi, eta, zeta);

		//find_point({ A, B, C, xi, eta, zeta });
		//verify_point({ A, B, C, xi, eta, zeta });

		uc.a = A;
		uc.b = B;
		uc.c = C;
		uc.alpha = xi;
		uc.beta = eta;
		uc.gamma = zeta;

		return true;
	}
};

