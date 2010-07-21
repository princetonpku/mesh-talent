#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cfloat>

#include "DeformableMesh3d.h"
#include "DeformationGraph.h"
#include "GaussNewtonSolver.h"
#include "debugtrace.h"

namespace meshtalent {

template <typename FUN>
void DeformableMesh3d::specifyposArr(FUN fun)
{
	std::vector<InterMesh::VertexHandle>& svIDs = pMesh->property(selectedVertices);
	sort(svIDs.begin(), svIDs.end());

#ifndef NDEBUG
	// assert handleIDs is sorted.
	for (int i = 0; i < static_cast<int>(handleIDs.size() - 1); ++i) { 
		assert(handleIDs[i] < handleIDs[i+1]);
	}
#endif

	specifiedposArr.clear();
	int handleIDsize = static_cast<int>(handleIDs.size());
	specifiedposArr.reserve(handleIDsize);
	for (int i = 0; i < handleIDsize; ++i) {
		InterMesh::Point& p = pMesh->point(handleIDs[i]); // get one handle.
		P3d pnow(p[0], p[1], p[2]);
		bool b = binary_search(svIDs.begin(), svIDs.end(), handleIDs[i]); // we can also use sequential search to speed up.
		if (b) { // this handle should move
			pnow = fun(pnow);
		}
		specifiedposArr.push_back(pnow);
	}
}

void DeformableMesh3d::updateVerticesUsingGraph()
{
	InterMesh::VertexIter v_it, v_end(pMesh->vertices_end());
	for (v_it = pMesh->vertices_begin(); v_it != v_end; ++v_it) {
		// get point pos.
		InterMesh::Point& vi = pMesh->point(v_it);
		P3d vi_p3d(vi[0], vi[1], vi[2]);
		P3d vihat_p3d(0.0, 0.0, 0.0);
		double wsum = 0.0;
		// get nearnodes.
		std::vector<IndexWeightPair>& nearnodes = pMesh->property(nearnodesArr, v_it);
		int nearnodessize = nearnodes.size();
		assert(nearnodessize == pGraph->relatenum);
		for (int k = 0; k < nearnodessize; ++k) {
			int j = nearnodes[k].first;
			double w = nearnodes[k].second;
			//assert(w < 1.0 && w > 0.0);
			wsum += w;
			DeformationGraph::GraphNode& gnode = pGraph->edges[j].first;
			P3d tmp = gnode.g + gnode.R * (vi_p3d - gnode.g) + gnode.t;
			vihat_p3d[0] += tmp[0] * w;
			vihat_p3d[1] += tmp[1] * w;
			vihat_p3d[2] += tmp[2] * w;
		}
		assert(fabs(1.0 - wsum) < FLT_EPSILON);
		vi = InterMesh::Point(vihat_p3d[0], vihat_p3d[1], vihat_p3d[2]);
	}
}

void DeformableMesh3d::deform()
{
	// try to iterate times from that each node has a init value of identity matrix and zero vector.
	GaussNewtonSolver gns(pGraph, pGraph->wrot, pGraph->wreg, pGraph->wcon, pGraph->wlen);
	std::vector<double> x[2];
	const int step = gns.step;
	const int nodesize = static_cast<int>(pGraph->edges.size());
	const double initarray[gns.step] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
	x[0].reserve(nodesize*step);
	for (int i = 0; i < nodesize; ++i) {
		for (int j = 0; j < step; ++j) {
			x[0].push_back(initarray[j]);
		}
	}
#ifdef DEBUGTRACE
	meshtalent::DebugTrace dt("./iterate.log");
	clock_t beg, end;
	beg = clock();
#endif
	int times = 0;
	gns.iterateOnce(x[0], x[1]);
	double result[2];
	result[0] = gns.F(x[0]);
	result[1] = gns.F(x[1]);
	int aswitch = 0;
#ifdef DEBUGTRACE
	dt.Trace("the %d time: result is %lf\n", times+1, result[aswitch]);
#endif
	while (result[aswitch] > result[1-aswitch]) {
		++times;
		aswitch = 1 - aswitch;
		x[1-aswitch].clear();
		gns.iterateOnce(x[aswitch], x[1-aswitch]);
		result[1-aswitch] = gns.F(x[1-aswitch]);
#ifdef DEBUGTRACE
		dt.Trace("the %d time: result is %lf\n", times+1, result[aswitch]);
#endif
	}
#ifdef DEBUGTRACE
	end = clock();
	dt.Trace("The time of iterating is %lf\n", (double)(end-beg) / CLOCKS_PER_SEC);
	dt.Trace("All iterated %d times\n.", times);
#endif
	// apply the right x to pGraph, compute all new positions of the mesh vertices, then compute new pos of graph nodes.
	pGraph->getFromx(x[aswitch]);
	updateVerticesUsingGraph();
	pGraph->updateNodesPos();
}

void DeformableMesh3d::translate(const V3d& t)
{
	specifyposArr(translateFun(t));
	deform();
}

void DeformableMesh3d::scale(double lamada)
{
	P3d pcenter = centerOfSV();
	specifyposArr(scaleFun(pcenter, lamada));
	deform();
}

void DeformableMesh3d::rotate(const V3d& t, const P3d& p, double beta)
{
	specifyposArr(rotateFun(p, t.normalize(), beta));
	deform();
}

DeformableMesh3d::P3d DeformableMesh3d::centerOfSV() const
{
	InterMesh::Point psum(0.0, 0.0, 0.0);
	InterMesh::VertexIter v_it, v_end(pMesh->vertices_end());
	for (v_it = pMesh->vertices_begin(); v_it != v_end; ++v_it) {
		psum += pMesh->point(v_it);
	}
	psum /= static_cast<double>(pMesh->n_vertices());
	P3d p3dsum(psum[0], psum[1], psum[2]);
	return p3dsum;
}

class CmpDist {
public:
	typedef DeformableMesh3d::P3d P3d;
public:
	CmpDist(const P3d& _vp, 
		const std::vector<P3d>& _nodecordcoll) : vp(_vp), nodecordcoll(_nodecordcoll) {}
	bool operator() (int lhs, int rhs) {
		return (nodecordcoll[lhs] - vp).magsqr() < (nodecordcoll[rhs] - vp).magsqr();
	}
private:
	const P3d vp;
	const std::vector<P3d>& nodecordcoll;
};

class CmpFirst {
public:
	template <typename T1, typename T2>
	bool operator() (const std::pair<T1,T2>& lhs, const std::pair<T1,T2>& rhs) {
		return (lhs.first < rhs.first) || (lhs.first == rhs.first && lhs.second < rhs.second);
	}
};

void DeformableMesh3d::InitDatas(DeformationGraph* _pGraph)
{
	setGraph(_pGraph);

	// copy all the node's coordinate into "nodecordcoll".
	int rnumplusone = pGraph->relatenum + 1;
	int nodenum = pGraph->nodenum;
	assert(rnumplusone < nodenum);

	std::vector<DeformableMesh3d::P3d> nodecordcoll; 
	nodecordcoll.reserve(nodenum);
	assert(nodenum == static_cast<int>(pGraph->edges.size()));
	for (int i = 0; i < nodenum; ++i) {
		nodecordcoll.push_back(pGraph->edges[i].first.g);
	}

	std::vector<int> nodes(nodenum);
	InterMesh::VertexIter v_it, v_end(pMesh->vertices_end());
	int k;
	for (v_it = pMesh->vertices_begin(), k = 0; v_it != v_end; ++v_it, ++k) {
		// init nodes, nodes will be partial sorted.
		for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
			nodes[i] = i;
		}
		InterMesh::Point& vp = pMesh->point(v_it);
		P3d tp(vp[0], vp[1], vp[2]);
		std::vector<int>::iterator sortend = nodes.begin() + rnumplusone;
		partial_sort(nodes.begin(), sortend, nodes.end(), CmpDist(tp, nodecordcoll));
		double maxdist = (nodecordcoll[nodes[rnumplusone-1]] - tp).mag();

		// Construct nearnodesArr and compute weights.
		std::vector<IndexWeightPair>& nearnodesnow = pMesh->property(nearnodesArr, v_it);
		for (int i = 0; i < rnumplusone - 1; ++i) {
			double d = (nodecordcoll[nodes[i]] - tp).mag();
			//assert(d < maxdist && d >= 0.0);
			double w = (1 - d / maxdist) * (1 - d / maxdist);
			nearnodesnow.push_back(std::make_pair(nodes[i], w));
			// put this vertex index in the node.
			pGraph->edges[nodes[i]].first.vertices.push_back(k);
		}
		sort(nearnodesnow.begin(), nearnodesnow.end(), CmpFirst());
		// normalize weight in nearnodesArr[k].
		double dsum = 0.0;
		assert(rnumplusone - 1 == static_cast<int>(nearnodesnow.size()));
		for (int i = 0; i < rnumplusone - 1; ++i) {
			dsum += nearnodesnow[i].second;
		}
		for (int i = 0; i < rnumplusone - 1; ++i) {
			nearnodesnow[i].second /= dsum;
		}
	} // end of for 
}

void DeformableMesh3d::gethandles()
{
	handleIDs.clear();
	std::vector<InterMesh::VertexHandle>& svIDs = pMesh->property(selectedVertices);
	copy(svIDs.begin(), svIDs.end(), back_inserter(handleIDs));
	sort(handleIDs.begin(), handleIDs.end());
}

double DeformableMesh3d::boundingboxVolume() const
{
	double coordmax[3] = { DBL_MIN, DBL_MIN, DBL_MIN };
	double coordmin[3] = { DBL_MAX, DBL_MAX, DBL_MAX };
	InterMesh::VertexIter v_it, v_end(pMesh->vertices_end());
	for (v_it = pMesh->vertices_begin(); v_it != v_end; ++v_it) {
		InterMesh::Point& p = pMesh->point(v_it);
		for (int j = 0; j < 3; ++j) {
			if (coordmax[j] < p[j]) {
				coordmax[j] = p[j];
			}
			if (coordmin[j] > p[j]) {
				coordmin[j] = p[j];
			}
		}
	}
	double result = (coordmax[0]-coordmin[0]) * (coordmax[1]-coordmin[1]) * (coordmax[2]-coordmin[2]);
	assert(result > 0);
	return result;
}

double DeformableMesh3d::surfaceArea() const
{
	double result = 0.0;
	InterMesh::FaceIter f_it, f_end(pMesh->faces_end());
	for (f_it = pMesh->faces_begin(); f_it != f_end; ++f_it) {
		P3d p3d[3];
		//assert((*f_it).is_triangle());
		InterMesh::FaceVertexIter fv_it(pMesh->fv_iter(f_it.handle()));
		int j = 0;
		for (; fv_it; ++fv_it, ++j) {
			InterMesh::Point& p = pMesh->point(fv_it);
			p3d[j] = P3d(p[0], p[1], p[2]);
		}
		V3d v3d = (p3d[1]-p3d[0]) % (p3d[2]-p3d[0]);
		result += v3d.mag() / 2;
	}
	return result;
}

} // end of namespace meshtalent
