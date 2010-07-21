#ifndef GAUSSNEWTONSOLVER_H
#define GAUSSNEWTONSOLVER_H

#include <vector>
#include <cmath>
#include <cassert>

#include "math/vector3d.h"
#include "math/point3d.h"
#include "math/matrix3d.h"
#include "math/SparseMatrix.h"

namespace meshtalent {

class DeformationGraph;

class GaussNewtonSolver {
public:
	typedef meshtalent::math::Vector3d<double> V3d;
	typedef meshtalent::math::Point3d<double>  P3d;
	typedef meshtalent::math::Matrix3d<double> M3d;
public:
	GaussNewtonSolver(const DeformationGraph* _dg, double wrot, double wreg, double wcon, double wlen) 
		: dg(_dg), wrotsqrt(sqrt(wrot)), wregsqrt(sqrt(wreg)), wconsqrt(sqrt(wcon)), wlensqrt(sqrt(wlen)) {}
	double F(const std::vector<double>& x) const;
	void f(const std::vector<double>& x, std::vector<double>& result) const;
	/* compute the inverse of the Jacobi Matrix, the size is n * m, n is num of unknown number, m is num of functions.
	   in : x, out : J */
	void JacobiInv(const std::vector<double>& x, SparseMatrix<double, int>& J) const;
	void gradsofF(const std::vector<double>& x, std::vector<double>& result) const;
	void iterateOnce(const std::vector<double>& xk, std::vector<double>& xkp1) const;
	bool isConvergent(const std::vector<double>& xkm1, const std::vector<double>& xk, const std::vector<double>& xkp1) const;
	void solve(const std::vector<double>& x, const double eps) const; // eps is epsilon.
public:
	static const int step = 12;
private:
	// compute f's components.
	double computeCC(const std::vector<double>& x, int i, int lhs, int rhs) const;
	double computeCCm1(const std::vector<double>& x, int i, int lhs) const;
	V3d computeNeighborReg(const std::vector<double>& x, int j, int k) const;
	V3d computeConstraints(const std::vector<double>& x, int i) const; // index is a index of dg->dmesh.handleIDs.
	V3d computeLength(const std::vector<double>& x, int j, int k) const; // add length component for f.

	// compute JacobiInv's components.
	void addCC(const std::vector<double>& x, SparseMatrix<double, int>& J, int nowindex) const;
	void addCCm1(const std::vector<double>& x, SparseMatrix<double, int>& J, int nowindex) const;
	void addNeighborReg(const std::vector<double>& x, SparseMatrix<double, int>& J,
		int nowindexj, int nowindexk, const V3d& gk_gj) const;
	void addConstraints(const std::vector<double>& x, SparseMatrix<double, int>& J, int i) const;
	void addLength(const std::vector<double>& x, SparseMatrix<double, int>& J,
		int nowindexj, int nowindexk) const; // add length component for JacobiInv.

private:
	const DeformationGraph* dg;
	const double wrotsqrt;
	const double wregsqrt;
	const double wconsqrt;
	const double wlensqrt;
};

} // end of namespace meshtalent

#endif // GAUSSNEWTONSOLVER_H
