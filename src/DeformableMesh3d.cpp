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
void DeformableMesh3d::specifyposArr(const std::vector<int>& selectedHandles, FUN fun)
{
#ifndef NDEBUG
	// assert selectedHandles is sorted.
	for (int i = 0; i < static_cast<int>(selectedHandles.size() - 1); ++i) {
		assert(selectedHandles[i] < selectedHandles[i+1]);
	}
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
		specifiedposArr.push_back(pnow);
	}

	for (int i = 0; i < selectedHandles.size(); ++i) {
		int index = selectedHandles[i];
		specifiedposArr[index] = fun(specifiedposArr[index]);
	}
}

DeformableMesh3d::DeformableMesh3d(InterMesh* _pMesh) : pMesh(_pMesh), pGraph(NULL) 
{
	pMesh->add_property(selectedVertices);
	pMesh->add_property(nearnodesArr);

	// for voronoi.
	pMesh->add_property(voroInfo);
	pMesh->add_property(halfEdgeDis);
	pMesh->add_property(overlapped_voroInfo);
	// compute halfEdgeDis.
	InterMesh::HalfedgeIter he_it(pMesh->halfedges_begin()), he_end(pMesh->halfedges_end());
	for (; he_it != he_end; ++he_it) {
		InterMesh::HalfedgeHandle heh = he_it.handle();
		InterMesh::VertexHandle phfrom = pMesh->from_vertex_handle(heh);
		InterMesh::VertexHandle phto = pMesh->to_vertex_handle(heh);
		InterMesh::Point pf = pMesh->point(phfrom);
		InterMesh::Point pt = pMesh->point(phto);
		P3d p3dfrom(pf[0], pf[1], pf[2]);
		P3d p3dto(pt[0], pt[1], pt[2]);

		pMesh->property(halfEdgeDis, he_it) = (p3dto - p3dfrom).mag();
	}
}

DeformableMesh3d::~DeformableMesh3d() 
{ 
	// for voronoi.
	pMesh->remove_property(overlapped_voroInfo);
	pMesh->remove_property(halfEdgeDis);
	pMesh->remove_property(voroInfo);

	pMesh->remove_property(nearnodesArr); 
	pMesh->remove_property(selectedVertices);
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

void DeformableMesh3d::translate(const std::vector<int>& selectedHandles, const V3d& t)
{
	specifyposArr(selectedHandles, translateFun(t));
	deform();
}

void DeformableMesh3d::scale(const std::vector<int>& selectedHandles, double lamada)
{
	P3d pcenter = centerOfSV();
	specifyposArr(selectedHandles, scaleFun(pcenter, lamada));
	deform();
}

void DeformableMesh3d::rotate(const std::vector<int>& selectedHandles, const V3d& t, const P3d& p, double beta)
{
	specifyposArr(selectedHandles, rotateFun(p, t.normalize(), beta));
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

	
	//genVoronoi();
	genOverlappedVoronoi();

	// init the map from vertexID to nodeID.
	std::map<int, int> nodesmap;
	for (int i = 0; i < pGraph->edges.size(); ++i) {
		nodesmap.insert(std::make_pair(pGraph->edges[i].first.vertexID, i));
	}

	// init nearnodes of every vertex.
	InterMesh::VertexIter v_it, v_end(pMesh->vertices_end());
	for (v_it = pMesh->vertices_begin(); v_it != v_end; ++v_it) {
		const OverlappedVoronoiInfo& OLVINow = pMesh->property(overlapped_voroInfo, v_it.handle());
		assert(OLVINow.indexNow == OverlappedVoronoiInfo::nearNodesNum);

		double maxdist = OLVINow.vis[OLVINow.indexNow-1].dis;
		std::vector<IndexWeightPair>& nearnodesnow = pMesh->property(nearnodesArr, v_it);

		for (int i = 0; i < OverlappedVoronoiInfo::nearNodesNum - 1; ++i) {
			double d = OLVINow.vis[i].dis;
			double w = (1 - d / maxdist) * (1 - d / maxdist);
			int indexOfNode = nodesmap[OLVINow.vis[i].handleRoot.idx()];
			nearnodesnow.push_back(std::make_pair(indexOfNode, w));
			// put this vertex index in the node.
			pGraph->edges[indexOfNode].first.vertices.push_back(v_it.handle().idx());

		}
		sort(nearnodesnow.begin(), nearnodesnow.end(), CmpFirst());
		// normalize weight in nearnodesArr[k].
		double dsum = 0.0;
		assert(OverlappedVoronoiInfo::nearNodesNum - 1 == static_cast<int>(nearnodesnow.size()));
		for (int i = 0; i < OverlappedVoronoiInfo::nearNodesNum - 1; ++i) {
			dsum += nearnodesnow[i].second;
		}
		for (int i = 0; i < OverlappedVoronoiInfo::nearNodesNum - 1; ++i) {
			nearnodesnow[i].second /= dsum;
		}
	}
	
	


	/*
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
	*/
	

}

void DeformableMesh3d::gethandles()
{
	handleIDs.clear();
	std::set<InterMesh::VertexHandle>& svIDs = pMesh->property(selectedVertices);
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

struct ElementInPQ {
	typedef DeformableMesh3d::InterMesh::VertexHandle VertexHandle;
	explicit ElementInPQ(VertexHandle _vh = VertexHandle(), double _d = 10e8) : vh(_vh), d(_d) {}

	VertexHandle vh;
	double d;
};

void DeformableMesh3d::genVoronoi()
{
	// init all vertex's voronoi information.
	InterMesh::VertexIter v_it = pMesh->vertices_begin(), v_end = pMesh->vertices_end();
	for (; v_it != v_end; ++v_it) {
		pMesh->property(voroInfo, v_it.handle()) = VoronoiInfo();
	}
	// init all node's voronoi information.
	int nodeSize = pGraph->edges.size();
	for (int i = 0; i < nodeSize; ++i) {
		InterMesh::VertexHandle vh(pGraph->edges[i].first.vertexID);
		pMesh->property(voroInfo, vh) = VoronoiInfo(0.0, vh);
	}

	// add a property to vertices to tell this vertex's pos in PQ.
	OpenMesh::VPropHandleT<int> posInPQ;
	pMesh->add_property(posInPQ);

	// init the PQ.
	std::vector<ElementInPQ> PQ;
	PQ.reserve(pMesh->n_vertices() + 1);
	PQ.push_back(ElementInPQ());
	for (int i = 0; i < nodeSize; ++i) {
		InterMesh::VertexHandle vh(pGraph->edges[i].first.vertexID);
		pMesh->property(posInPQ, vh) = PQ.size();
		PQ.push_back(ElementInPQ(vh, 0.0));
	}
	for (v_it = pMesh->vertices_begin(); v_it != v_end; ++v_it) {
		if (pMesh->property(voroInfo, v_it.handle()).dis != 0.0) { // it is a node vertex.
			pMesh->property(posInPQ, v_it.handle()) = PQ.size();
			PQ.push_back(ElementInPQ(v_it.handle()));
		}
	}

	// dijkstra.
	while (PQ.size() > 1) {
		// extract min.
		ElementInPQ e = PQ[1];
		PQ[1] = PQ.back();
		pMesh->property(posInPQ, PQ[1].vh) = 1;
		PQ.pop_back();
		int index = 1;
		int indexChild = index * 2;
		while (indexChild < PQ.size()) {
			if (indexChild + 1 < PQ.size() && PQ[indexChild+1].d < PQ[indexChild].d) { indexChild++; }
			if (PQ[index].d <= PQ[indexChild].d) { break; }
			// swap.
			ElementInPQ te = PQ[index];
			PQ[index] = PQ[indexChild];
			PQ[indexChild] = te;
			pMesh->property(posInPQ, PQ[index].vh) = index;
			pMesh->property(posInPQ, PQ[indexChild].vh) = indexChild;

			index = indexChild;
			indexChild = index * 2;
		}

		// relax.
		InterMesh::VertexHandle vh = e.vh;
		const VoronoiInfo& VIParent = pMesh->property(voroInfo, vh);
		InterMesh::VertexOHalfedgeIter voh_it(*pMesh, vh);
		for (; voh_it; ++voh_it) {
			double disOfEdge = pMesh->property(halfEdgeDis, voh_it.handle());
			InterMesh::VertexHandle vhto = pMesh->to_vertex_handle(voh_it.handle());
			VoronoiInfo& VINow = pMesh->property(voroInfo, vhto);
			if (VIParent.dis + disOfEdge < VINow.dis) {
				// should relax here.
				VINow.dis = VIParent.dis + disOfEdge;
				VINow.handleRoot = VIParent.handleRoot;
				int indexInPQ = pMesh->property(posInPQ, vhto);
				PQ[indexInPQ].d = VINow.dis;
				int indexNow = indexInPQ;

				while (indexNow/2 > 0 && PQ[indexNow].d < PQ[indexNow/2].d) {
					// swap.
					ElementInPQ te = PQ[indexNow];
					PQ[indexNow] = PQ[indexNow/2];
					PQ[indexNow/2] = te;
					pMesh->property(posInPQ, PQ[indexNow].vh) = indexNow;
					pMesh->property(posInPQ, PQ[indexNow/2].vh) = indexNow/2;
					indexNow /= 2;
				} // end of while.
			} // end of if.
		} // end of for.
	} // end of while.

	pMesh->remove_property(posInPQ);

	// set color.
	for (v_it = pMesh->vertices_begin(); v_it != v_end; ++v_it) {
		const VoronoiInfo& viNow = pMesh->property(voroInfo, v_it.handle());
		int idx = viNow.handleRoot.idx();
		pMesh->set_color(v_it, InterMesh::Color(idx, idx*idx, 255-idx));
	}
}

int DeformableMesh3d::OverlappedVoronoiInfo::nearNodesNum = 4;

void DeformableMesh3d::genOverlappedVoronoi()
{
	// init all vertex's overlapped voronoi information.
	const int nearNodesNum = pGraph->relatenum + 1;
	OverlappedVoronoiInfo::setNearNodesNum(nearNodesNum);
	InterMesh::VertexIter v_it = pMesh->vertices_begin(), v_end = pMesh->vertices_end();
	for (; v_it != v_end; ++v_it) {
		pMesh->property(overlapped_voroInfo, v_it.handle()) = OverlappedVoronoiInfo();
	}
	// init all node's overlapped voronoi information.
	int nodeSize = pGraph->edges.size();
	for (int i = 0; i < nodeSize; ++i) {
		InterMesh::VertexHandle vh(pGraph->edges[i].first.vertexID);
		OverlappedVoronoiInfo& olviNow = pMesh->property(overlapped_voroInfo, vh);
		olviNow.vis[0] = VoronoiInfo(0.0, vh);
	}

	// add a property to vertices to tell this vertex's pos in PQ.
	OpenMesh::VPropHandleT<int> posInPQ;
	pMesh->add_property(posInPQ);

	// init the PQ, the same with voronoi.
	std::vector<ElementInPQ> PQ;
	PQ.reserve(pMesh->n_vertices() + 1);
	PQ.push_back(ElementInPQ());
	for (int i = 0; i < nodeSize; ++i) {
		InterMesh::VertexHandle vh(pGraph->edges[i].first.vertexID);
		pMesh->property(posInPQ, vh) = PQ.size();
		PQ.push_back(ElementInPQ(vh, 0.0));
	}
	for (v_it = pMesh->vertices_begin(); v_it != v_end; ++v_it) {
		if (pMesh->property(overlapped_voroInfo, v_it.handle()).vis[0].dis != 0.0) { // it is not a node vertex.
			pMesh->property(posInPQ, v_it.handle()) = PQ.size();
			PQ.push_back(ElementInPQ(v_it.handle()));
		}
	}

	// dijkstra.
	while (PQ.size() > 1) {
		// extract min.
		ElementInPQ e = PQ[1];
		PQ[1] = PQ.back();
		pMesh->property(posInPQ, PQ[1].vh) = 1;
		PQ.pop_back();
		int index = 1;
		int indexChild = index * 2;
		while (indexChild < PQ.size()) {
			if (indexChild + 1 < PQ.size() && PQ[indexChild+1].d < PQ[indexChild].d) { indexChild++; }
			if (PQ[index].d <= PQ[indexChild].d) { break; }
			// swap.
			ElementInPQ te = PQ[index];
			PQ[index] = PQ[indexChild];
			PQ[indexChild] = te;
			pMesh->property(posInPQ, PQ[index].vh) = index;
			pMesh->property(posInPQ, PQ[indexChild].vh) = indexChild;

			index = indexChild;
			indexChild = index * 2;
		}

		// relax.
		InterMesh::VertexHandle vh = e.vh;
		OverlappedVoronoiInfo& OLVIParent = pMesh->property(overlapped_voroInfo, vh);
		InterMesh::VertexOHalfedgeIter voh_it(*pMesh, vh);
		for (; voh_it; ++voh_it) {
			double disOfEdge = pMesh->property(halfEdgeDis, voh_it.handle());
			InterMesh::VertexHandle vhto = pMesh->to_vertex_handle(voh_it.handle());

			OverlappedVoronoiInfo& OLVINow = pMesh->property(overlapped_voroInfo, vhto);
			bool keyUpdated = overlappedVoronoiRelax(OLVIParent, OLVINow, disOfEdge);

			if (keyUpdated) {
				int indexInPQ = pMesh->property(posInPQ, vhto);
				PQ[indexInPQ].d = OLVINow.vis[OLVINow.indexNow].dis;
				int indexNow = indexInPQ;

				while (indexNow/2 > 0 && PQ[indexNow].d < PQ[indexNow/2].d) {
					// swap.
					ElementInPQ te = PQ[indexNow];
					PQ[indexNow] = PQ[indexNow/2];
					PQ[indexNow/2] = te;
					pMesh->property(posInPQ, PQ[indexNow].vh) = indexNow;
					pMesh->property(posInPQ, PQ[indexNow/2].vh) = indexNow/2;
					indexNow /= 2;
				} // end of while.
			} // end of if.
		} // end of for.
		
		// if the times of this Element being popped is not enough, push it back to PQ.
		if (++OLVIParent.indexNow < OverlappedVoronoiInfo::nearNodesNum) {
			assert(e.d <= OLVIParent.vis[OLVIParent.indexNow].dis);
			e.d = OLVIParent.vis[OLVIParent.indexNow].dis;
			// push back and sift up.
			PQ.push_back(e);
			int indexNow = PQ.size() - 1;
			while (indexNow/2 > 0 && PQ[indexNow].d < PQ[indexNow/2].d) {
				// swap.
				ElementInPQ te = PQ[indexNow];
				PQ[indexNow] = PQ[indexNow/2];
				PQ[indexNow/2] = te;
				pMesh->property(posInPQ, PQ[indexNow].vh) = indexNow;
				pMesh->property(posInPQ, PQ[indexNow/2].vh) = indexNow/2;
				indexNow /= 2;
			} // end of while.
		} // end of if.
		
	} // end of while.

	pMesh->remove_property(posInPQ);
}
bool DeformableMesh3d::overlappedVoronoiRelax(const OverlappedVoronoiInfo& OLVIParent, OverlappedVoronoiInfo& OLVINow, double disOfEdge)
{
	VoronoiInfo buffer[OverlappedVoronoiInfo::MAXSIZE * 2];
	int indexInBuffer = 0;

	assert(OLVIParent.indexNow < OverlappedVoronoiInfo::nearNodesNum);

	if (OLVINow.indexNow >= OverlappedVoronoiInfo::nearNodesNum) { return false; }

	// put parent's vi in buffer.
	for (int i = OLVIParent.indexNow; i < OverlappedVoronoiInfo::nearNodesNum; ++i) {
		VoronoiInfo vi = OLVIParent.vis[i];
		bool inNowPoped = false;
		for (int j = 0; j < OLVINow.indexNow; ++j) {
			if (vi.handleRoot == OLVINow.vis[j].handleRoot) {
				assert(OLVINow.vis[j].dis <= vi.dis + disOfEdge);
				inNowPoped = true;
				break;
			}
		}
		if (!inNowPoped) { 
			vi.dis += disOfEdge;
			buffer[indexInBuffer++] = vi;
	   	}
	}
	// put now's vi in buffer.
	for (int i = OLVINow.indexNow; i < OverlappedVoronoiInfo::nearNodesNum; ++i) {
		buffer[indexInBuffer++] = OLVINow.vis[i];
	}
	//root duplicate remove and sort.
	std::sort(buffer, buffer + indexInBuffer, VoronoiInfo::DisCmp());
	int indexInOLVINow = OLVINow.indexNow;
	for (int i = 0; i < indexInBuffer; ++i) {
		bool alreadyHasThisRoot = false;
		for (int j = OLVINow.indexNow; j < indexInOLVINow; ++j) {
			if (OLVINow.vis[j].handleRoot == buffer[i].handleRoot) {
				alreadyHasThisRoot = true;
				break;
			}
		}
		if (!alreadyHasThisRoot) {
			OLVINow.vis[indexInOLVINow++] = buffer[i];
		}
		if (indexInOLVINow >= OverlappedVoronoiInfo::nearNodesNum) { // has already filled all vis in OLVINow.
			break;
		}
	} // end of for.
	return true;
}


} // end of namespace meshtalent
