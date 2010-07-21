#include <cmath>
#include <cassert>
#include <numeric>

#include <cholmod.h>

#include "GaussNewtonSolver.h"
#include "DeformableMesh3d.h"
#include "DeformationGraph.h"
#include "debugtrace.h"

namespace meshtalent {

double GaussNewtonSolver::F(const std::vector<double>& x) const
{
	std::vector<double> result;
	f(x, result);
	return inner_product(result.begin(), result.end(), result.begin(), 0.0);
}

double GaussNewtonSolver::computeCC(const std::vector<double>& x, int i, int lhs, int rhs) const
{
	assert(lhs != rhs && lhs < 3 && rhs < 3);
	const double* xl = &x[i*step+lhs];
	const double* xr = &x[i*step+rhs];
	return xl[0] * xr[0] + xl[3] * xr[3] + xl[6] * xr[6];
}
double GaussNewtonSolver::computeCCm1(const std::vector<double>& x, int i, int lhs) const
{
	assert(lhs < 3);
	const double* xl = &x[i*step+lhs];
	return xl[0] * xl[0] + xl[3] * xl[3] + xl[6] * xl[6] - 1.0;
}

GaussNewtonSolver::V3d GaussNewtonSolver::computeNeighborReg(const std::vector<double>& x, int j, int k) const
{
	const double* xj = &x[j*step];
	const double* xk = &x[k*step];
	const P3d& gk = dg->edges[k].first.g;
	const P3d& gj = dg->edges[j].first.g;
	V3d gk_gj = gk - gj;
	V3d v;
	// compute each components of v separately, just for speed, easy to get coding error.
	v[0] = (xj[0] * gk_gj[0] + xj[1] * gk_gj[1] + xj[2] * gk_gj[2]) + 
		gj[0] + xj[9] - (gk[0] + xk[9]);
	v[1] = (xj[3] * gk_gj[0] + xj[4] * gk_gj[1] + xj[5] * gk_gj[2]) + 
		gj[1] + xj[10] - (gk[1] + xk[10]);
	v[2] = (xj[6] * gk_gj[0] + xj[7] * gk_gj[1] + xj[8] * gk_gj[2]) + 
		gj[2] + xj[11] - (gk[2] + xk[11]);
	return v;
}

GaussNewtonSolver::V3d GaussNewtonSolver::computeConstraints(const std::vector<double>& x, int i) const
{
	typedef DeformableMesh3d::InterMesh InterMesh;
	DeformableMesh3d& dmesh = dg->dmesh;
	InterMesh::VertexHandle index = dmesh.handleIDs[i]; // index is a handle.
	std::vector<DeformableMesh3d::IndexWeightPair>& nearnodes = 
		dmesh.pMesh->property(dmesh.nearnodesArr, index);
	assert(nearnodes.size() == static_cast<size_t>(dg->relatenum));
#ifndef NDEBUG
	// assert nearnodes is sorted by index.
	for (int ii = 0; ii < static_cast<int>(nearnodes.size()-1); ++ii) {
		assert(nearnodes[ii].first < nearnodes[ii+1].first);
	}
#endif
	InterMesh::Point& vi = dmesh.pMesh->point(index);
	P3d vip(0.0, 0.0, 0.0);
	for (int k = 0; k < static_cast<int>(nearnodes.size()); ++k) {
		int j = nearnodes[k].first; // the graph node's index in dg->edges.
		assert(j < static_cast<int>(dg->edges.size()) && j >= 0);
		double w = nearnodes[k].second;
		const double* xj = &x[j*step];
		const P3d& gj = dg->edges[j].first.g;
		V3d tmpV(0.0, 0.0, 0.0);
		tmpV[0] = vi[0] - gj[0];
		tmpV[1] = vi[1] - gj[1];
		tmpV[2] = vi[2] - gj[2];
		vip[0] += (xj[0] * tmpV[0] + xj[1] * tmpV[1] + xj[2] * tmpV[2] + gj[0] + xj[9]) * w;
		vip[1] += (xj[3] * tmpV[0] + xj[4] * tmpV[1] + xj[5] * tmpV[2] + gj[1] + xj[10]) * w;
		vip[2] += (xj[6] * tmpV[0] + xj[7] * tmpV[1] + xj[8] * tmpV[2] + gj[2] + xj[11]) * w;
	}
	// now vip is the vi's new position using x's transformation.
	P3d q = dmesh.specifiedposArr[i];
	return vip-q;
}

GaussNewtonSolver::V3d GaussNewtonSolver::computeLength(const std::vector<double>& x, int j, int k) const
{
	const double* xj = &x[j*step];
	const double* xk = &x[k*step];
	V3d v(xj[9] - xk[9], xj[10] - xk[10], xj[11] - xk[11]);
	return v;
}

void GaussNewtonSolver::f(const std::vector<double>& x, std::vector<double>& result) const
{
	/* assume every 12 elements of x is ordered like this:
	   M[0][0], M[0][1], M[0][2], M[1][0], M[1][1], M[1][2], M[2][0], M[2][1], M[2][2], t[0], t[1], t[2] */
	const int nodesize = static_cast<int>(dg->edges.size());
	assert(static_cast<int>(x.size()) == nodesize * step);
	assert(result.size() == 0);
	result.reserve(x.size() * 2); // need to recompute ! 
	DeformableMesh3d& dmesh = dg->dmesh;
	// compute Erot's components.
	for (int i = 0; i < nodesize; ++i) {
		result.push_back(wrotsqrt * computeCC(x, i, 0, 1));
		result.push_back(wrotsqrt * computeCC(x, i, 0, 2));
		result.push_back(wrotsqrt * computeCC(x, i, 1, 2));
		result.push_back(wrotsqrt * computeCCm1(x, i, 0));
		result.push_back(wrotsqrt * computeCCm1(x, i, 1));
		result.push_back(wrotsqrt * computeCCm1(x, i, 2));
	}
	// compute Ereg's components.
	for (int j = 0; j < nodesize; ++j) {
		const std::vector<int>& edgelist = dg->edges[j].second;
		const int edgelistsize = static_cast<int>(edgelist.size());
		for (int i = 0; i < edgelistsize; ++i) {
			int k = edgelist[i];
			V3d v = computeNeighborReg(x, j, k);
			result.push_back(wregsqrt * v[0]);
			result.push_back(wregsqrt * v[1]);
			result.push_back(wregsqrt * v[2]);
		}
	}
	// compute Econ's components.
	int selectsize = static_cast<int>(dmesh.handleIDs.size());
	assert(static_cast<size_t>(selectsize) == dmesh.specifiedposArr.size());
	
	for (int i = 0; i < selectsize; ++i) {
		V3d vip_q = computeConstraints(x, i);
		result.push_back(wconsqrt * vip_q[0]);
		result.push_back(wconsqrt * vip_q[1]);
		result.push_back(wconsqrt * vip_q[2]);
	}
#ifdef ELEN
	// compute Elen's components.
	for (int j = 0; j < nodesize; ++j) {
		const std::vector<int>& edgelist = dg->edges[j].second;
		const int edgelistsize = static_cast<int>(edgelist.size());
		for (int i = 0; i < edgelistsize; ++i) {
			int k = edgelist[i];
			V3d v = computeLength(x, j, k);
			result.push_back(wlensqrt * v[0]);
			result.push_back(wlensqrt * v[1]);
			result.push_back(wlensqrt * v[2]);
		}
	}
#endif

	// assert num of function is more than num of unknown numbers.
	assert(result.size() > x.size());
}

void GaussNewtonSolver::addCC(const std::vector<double>& x, SparseMatrix<double, int>& J, int nowindex) const
{
	typedef std::pair<double, int> PAIR;
	std::vector<PAIR> acol;
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+1], nowindex + 0));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+0], nowindex + 1));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+4], nowindex + 3));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+3], nowindex + 4));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+7], nowindex + 6));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+6], nowindex + 7));
	J.addACol(acol);
	acol.clear();

	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+2], nowindex + 0));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+0], nowindex + 2));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+5], nowindex + 3));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+3], nowindex + 5));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+8], nowindex + 6));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+6], nowindex + 8));
	J.addACol(acol);
	acol.clear();

	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+2], nowindex + 1));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+1], nowindex + 2));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+5], nowindex + 4));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+4], nowindex + 5));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+8], nowindex + 7));
	acol.push_back(std::make_pair(wrotsqrt * x[nowindex+7], nowindex + 8));
	J.addACol(acol);

}

void GaussNewtonSolver::addCCm1(const std::vector<double>& x, SparseMatrix<double, int>& J, int nowindex) const
{
	typedef std::pair<double, int> PAIR;
	std::vector<PAIR> acol;
	for (int i = 0; i < 3; ++i) {
		acol.clear();
		acol.push_back(std::make_pair(wrotsqrt * 2 * x[nowindex+0+i], nowindex + 0 + i));
		acol.push_back(std::make_pair(wrotsqrt * 2 * x[nowindex+3+i], nowindex + 3 + i));
		acol.push_back(std::make_pair(wrotsqrt * 2 * x[nowindex+6+i], nowindex + 6 + i));
		J.addACol(acol);
	}
}

void GaussNewtonSolver::addNeighborReg(const std::vector<double>& x, SparseMatrix<double, int>& J,
									   int nowindexj, int nowindexk, const V3d& gk_gj) const
{
	typedef std::pair<double, int> PAIR;
	std::vector<PAIR> acol;
	// which should be push_back first depends on j' and k' value which is larger.
	assert(nowindexj != nowindexk);
	if (nowindexj > nowindexk) {
		for (int i = 0; i < 3; ++i) {
			acol.clear();
			acol.push_back(std::make_pair(-1 * wregsqrt, nowindexk + 9 + i));
			acol.push_back(std::make_pair(gk_gj[0] * wregsqrt, nowindexj + 0 + i*3));
			acol.push_back(std::make_pair(gk_gj[1] * wregsqrt, nowindexj + 1 + i*3));
			acol.push_back(std::make_pair(gk_gj[2] * wregsqrt, nowindexj + 2 + i*3));
			acol.push_back(std::make_pair(1 * wregsqrt, nowindexj + 9 + i));
			J.addACol(acol);
		}
	} else {
		for (int i = 0; i < 3; ++i) {
			acol.clear();
			acol.push_back(std::make_pair(gk_gj[0] * wregsqrt, nowindexj + 0 + i*3));
			acol.push_back(std::make_pair(gk_gj[1] * wregsqrt, nowindexj + 1 + i*3));
			acol.push_back(std::make_pair(gk_gj[2] * wregsqrt, nowindexj + 2 + i*3));
			acol.push_back(std::make_pair(1 * wregsqrt, nowindexj + 9 + i));
			acol.push_back(std::make_pair(-1 * wregsqrt, nowindexk + 9 + i));
			J.addACol(acol);
		}
	}
}

void GaussNewtonSolver::addConstraints(const std::vector<double>& x, SparseMatrix<double, int>& J, int i) const
{
	typedef DeformableMesh3d::InterMesh InterMesh;
	DeformableMesh3d& dmesh = dg->dmesh;
	InterMesh::VertexHandle index = dmesh.handleIDs[i]; // index is a handle.
	std::vector<DeformableMesh3d::IndexWeightPair>& nearnodes = 
		dmesh.pMesh->property(dmesh.nearnodesArr, index);
	assert(nearnodes.size() == static_cast<size_t>(dg->relatenum));
#ifndef NDEBUG
	// assert nearnodes is sorted by index.
	for (int i = 0; i < static_cast<int>(nearnodes.size()-1); ++i) {
		assert(nearnodes[i].first < nearnodes[i+1].first);
	}
#endif
	InterMesh::Point& vi = dmesh.pMesh->point(index);
	typedef std::pair<double, int> PAIR;
	std::vector<PAIR> acol;
	for (int ii = 0; ii < 3; ++ii) {
		acol.clear();
		for (int k = 0; k < static_cast<int>(nearnodes.size()); ++k) {
			int j = nearnodes[k].first; // the graph node's index in dg->edges.
			int nowindexj = j * step;
			double w = nearnodes[k].second;
			const P3d& gj = dg->edges[j].first.g;
			V3d vi_gj(0.0, 0.0, 0.0);
			vi_gj[0] = vi[0] - gj[0];
			vi_gj[1] = vi[1] - gj[1];
			vi_gj[2] = vi[2] - gj[2];
			// R.
			acol.push_back(std::make_pair(wconsqrt * w * vi_gj[0], nowindexj + 0 + ii*3));
			acol.push_back(std::make_pair(wconsqrt * w * vi_gj[1], nowindexj + 1 + ii*3));
			acol.push_back(std::make_pair(wconsqrt * w * vi_gj[2], nowindexj + 2 + ii*3));
			// t.
			acol.push_back(std::make_pair(wconsqrt * w, nowindexj + 9 + ii));
		}
		J.addACol(acol);
	}
}

void GaussNewtonSolver::addLength(const std::vector<double>& x, SparseMatrix<double, int>& J, int nowindexj, int nowindexk) const
{
	typedef std::pair<double, int> PAIR;
	std::vector<PAIR> acol;
	assert(nowindexj != nowindexk);
	if (nowindexj > nowindexk) {
		for (int i = 0; i < 3; ++i) {
			acol.clear();
			acol.push_back(std::make_pair(-1.0 * wlensqrt, nowindexk + 9 + i));
			acol.push_back(std::make_pair(1.0 * wlensqrt, nowindexj + 9 + i));
			J.addACol(acol);
		}
	} else {
		for (int i = 0; i < 3; ++i) {
			acol.clear();
			acol.push_back(std::make_pair(1.0 * wlensqrt, nowindexj + 9 + i));
			acol.push_back(std::make_pair(-1.0 * wlensqrt, nowindexk + 9 + i));
			J.addACol(acol);
		}
	}
}

void GaussNewtonSolver::JacobiInv(const std::vector<double>& x, SparseMatrix<double, int>& J) const
{
	const int nodesize = static_cast<int>(dg->edges.size());
	J.clear();
	assert(nodesize * step == static_cast<int>(x.size()));
	
	// part of Erot.
	for (int i = 0; i < nodesize; ++i) {
		int nowindex = i * step; // the ith node. x[nowindex] is R[0][0] of the ith node.
		addCC(x, J, nowindex);
		addCCm1(x, J, nowindex);
	}
	// part of Ereg.
#ifndef NDEBUG
	// assert all edges in one node's edgelist is sorted.
	for (int i = 0; i < nodesize; ++i) {
		const std::vector<int>& edgelist = dg->edges[i].second;
		for (int j = 0; j < static_cast<int>(edgelist.size() - 1); ++j) {
			assert(edgelist[j] < edgelist[j+1]);
		}
	}
#endif
	// compute.
	for (int j = 0; j < nodesize; ++j) {
		const std::vector<int>& edgelist = dg->edges[j].second;
		for (int i = 0; i < static_cast<int>(edgelist.size()); ++i) {
			int k = edgelist[i];
			const P3d& gk = dg->edges[k].first.g;
			const P3d& gj = dg->edges[j].first.g;
			V3d gk_gj = gk - gj;
			int indexnowj = j * step;
			int indexnowk = k * step;
			addNeighborReg(x, J, indexnowj, indexnowk, gk_gj);
		}
	}

	// part of Econ.
	int selectsize = static_cast<int>(dg->dmesh.handleIDs.size());
	assert(static_cast<size_t>(selectsize) == dg->dmesh.specifiedposArr.size());
	for (int i = 0; i < selectsize; ++i) {
		addConstraints(x, J, i);
	}

#ifdef ELEN
	// part of Elen.
	for (int j = 0; j < nodesize; ++j) {
		const std::vector<int>& edgelist = dg->edges[j].second;
		for (int i = 0; i < static_cast<int>(edgelist.size()); ++i) {
			int k = edgelist[i];
			int indexnowj = j * step;
			int indexnowk = k * step;
			addLength(x, J, indexnowj, indexnowk);
		}
	}
#endif
	
	// should assert colptr.size() is equal to result in f()'s size.
}

// unwritten.
//void GaussNewtonSolver::gradsofF(const std::vector<double>& x, std::vector<double>& result) const
//{
//	const int nodesize = static_cast<int>(dg->edges.size());
//	assert(x.size() == nodesize * step);
//	assert(result.size() == 0);
//	result.reserve(x.size());
//	
//}

void GaussNewtonSolver::iterateOnce(const std::vector<double>& xk, std::vector<double>& xkp1) const
{
	const int nodesize = static_cast<int>(dg->edges.size());
	assert(static_cast<int>(xk.size()) == nodesize * step);
	assert(xkp1.size() == 0);

	// compute f(xk).
	std::vector<double> fxk;
	f(xk, fxk);
	const int M = static_cast<int>(fxk.size());
	const int N = static_cast<int>(xk.size());

	// compute JT(xk).
	SparseMatrix<double, int> J;
	JacobiInv(xk, J);

	assert(static_cast<int>(J.colsize()) == M); // the column num should be the same with M.

	// get access with the internal of the JacobiInv.
	std::vector<int>& colptr = J.getColptr();
	std::vector<int>& rowind = J.getRowind();
	std::vector<double>& values = J.getValues();
	
#ifdef DEBUGTRACE
	meshtalent::DebugTrace dtj("./Jacobi.log");
	dtj.Trace("colptr :\n");
	for (int i = 0; i < static_cast<int>(colptr.size()); ++i) {
		dtj.Trace("%d\n", colptr[i]);
	}
	dtj.Trace("rowind :\n");
	for (int i = 0; i < static_cast<int>(rowind.size()); ++i) {
		dtj.Trace("%d\n", rowind[i]);
	}
	dtj.Trace("values :\n");
	for (int i = 0; i < static_cast<int>(values.size()); ++i) {
		dtj.Trace("%lf\n", values[i]);
	}
	dtj.Trace(".......\n");
	dtj.Trace(".......\n");
#endif

	// construct JT and JTJ.
	cholmod_common cm;
	cholmod_start(&cm);
	cholmod_sparse JT; // N * M.
	cholmod_sparse *pJTJ; // N * N.
	JT.nrow = N;
	JT.ncol = M;
	JT.nzmax = values.size();
	JT.p = &colptr[0];
	JT.i = &rowind[0];
	JT.nz = NULL;
	JT.x = &values[0];
	JT.z = NULL;
	JT.stype = 0;
	JT.itype = CHOLMOD_INT;
	JT.xtype = CHOLMOD_REAL;
	JT.dtype = CHOLMOD_DOUBLE;
	JT.sorted = 1;
	JT.packed = 1;
	pJTJ = cholmod_aat(&JT, NULL, 0, 1, &cm);
	pJTJ->stype = 1; // use upper triangle.

	cholmod_dense fxkdense;
	fxkdense.nrow = M;
	fxkdense.ncol = 1;
	fxkdense.nzmax = M;
	fxkdense.d = M;
	fxkdense.x = &fxk[0];
	fxkdense.xtype = CHOLMOD_REAL;
	fxkdense.dtype = CHOLMOD_DOUBLE;

	cholmod_dense *pNegJTfxk;
	pNegJTfxk = cholmod_zeros(N, 1, CHOLMOD_REAL, &cm);
	double coone[2] = { -1, 0 };
	double cozero[2] = { 0, 0 };
	
	/*meshtalent::DebugTrace dt("./cholmod_sdmult.log");
	dt.Trace("JT info : \n");
	dt.Trace("colptr total %d :\n", JT.ncol + 1);
	for (int i = 0; i <= JT.ncol; ++i) {
		dt.Trace("%d : %d\n", i, ((int*)(JT.p))[i]);
	}
	dt.Trace("rowind total %d : \n", JT.nzmax);
	for (int i = 0; i < JT.nzmax; ++i) {
		dt.Trace("%d : %d\n", i, ((int*)(JT.i))[i]);
	}
	dt.Trace("values total %d : \n", JT.nzmax);
	for (int i = 0; i < JT.nzmax; ++i) {
		dt.Trace("%d : %lf\n", i, ((double*)(JT.x))[i]);
	}
	dt.Trace("fxkdense info :\n");
	dt.Trace("x total %d :\n", fxkdense.nzmax);
	for (int i = 0; i < fxkdense.nzmax; ++i) {
		dt.Trace("%d : %lf\n", i, ((double*)(fxkdense.x))[i]);
	}*/
	int r = cholmod_sdmult(&JT, 0, coone, cozero, &fxkdense, pNegJTfxk, &cm);

#ifdef DEBUGTRACE
	// check pNegJTfxk if it is a zero vector.
	meshtalent::DebugTrace dt("./pNegJTfxk.log");
	for (int i = 0; i < N; ++i) {
		dt.Trace("%d : %lf\n", i, ((double*)(pNegJTfxk->x))[i]);
	}
	dt.Trace("r = %d\n", r);
#endif

	//solve equations.
	cholmod_factor *pL;
	cholmod_dense  *pfxkp1dense;
	pL = cholmod_analyze(&JT, &cm);
	cholmod_factorize(&JT, pL, &cm);
	pfxkp1dense = cholmod_solve(CHOLMOD_A, pL, pNegJTfxk, &cm);
	const double* pinner = (const double*)pfxkp1dense->x;
	xkp1.reserve(N);
	for (int i = 0; i < N; ++i) {
		xkp1.push_back(pinner[i] + xk[i]);
	}

	cholmod_free_sparse(&pJTJ, &cm);
	cholmod_free_dense(&pNegJTfxk, &cm);
	cholmod_free_factor(&pL, &cm);
	cholmod_free_dense(&pfxkp1dense, &cm);
	cholmod_finish(&cm);
}

bool GaussNewtonSolver::isConvergent(const std::vector<double>& xkm1, const std::vector<double>& xk, const std::vector<double>& xkp1) const
{
	const int nodesize = static_cast<int>(dg->edges.size());
	const int N = nodesize * step;
	const double eps = 1.0e-6;
	assert(static_cast<int>(xkm1.size()) == N &&
		   	static_cast<int>(xk.size()) == N &&
		   	static_cast<int>(xkp1.size()) == N);
	bool b1, b2, b3;
	// compute |Fk-Fkp1| < eps * (1+Fk).
	std::vector<double> fxkm1, fxk;
	f(xkm1, fxkm1);
	f(xk, fxk);
	double Fxkm1 = inner_product(fxkm1.begin(), fxkm1.end(), fxkm1.begin(), 0.0);
	double Fxk = inner_product(fxk.begin(), fxk.end(), fxk.begin(), 0.0);
	b1 = fabs(Fxk-Fxkm1) < eps * (1 + Fxk);
	// compute ||gradofFxk||endlessnorm < cubic_root(eps) * (1+Fk).
	b2 = true;
	// compute ||deltak||endlessnorm < square_root(eps) * (1+||deltak||endlessnorm).
	std::vector<double> deltak;
	deltak.reserve(N);
	for (int i = 0; i < N; ++i) {
		deltak.push_back(xkp1[i] - xk[i]);
	}
	double endlessnormOfdeltak = *max_element(deltak.begin(), deltak.end());
	b3 = endlessnormOfdeltak < pow(eps, 1.0/3) * (1 + endlessnormOfdeltak);
	// result.
	return b1 && b2 && b3;
}

} // end of namespace meshtalent
