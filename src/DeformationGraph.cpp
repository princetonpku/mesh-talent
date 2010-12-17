#include <list>
#include <vector>
#include <ctime>
#include <cstdlib>
// ------------------------------
#include "DeformationGraph.h"
#include "DeformableMesh3d.h"
#include "debugtrace.h"
// ------------------------------

#ifdef WIN32
static void srand48(unsigned int s) 
{
	srand(s);
}
static long lrand48()
{
	return rand() * (RAND_MAX+1) + rand();
}
static double drand48() 
{
	return lrand48() / (static_cast<double>(RAND_MAX+1) * static_cast<double>(RAND_MAX+1));
}
#endif

namespace meshtalent {

class RandDoubleGenerator {
public:
	static void srand(unsigned int s) {
		::srand48(s);
	}
	double operator() (double low, double high) const {
		return low + (high - low) * ::drand48();
	}
};

class RandUIntGenerator {
public:
	static void srand(unsigned int s) {
		::srand48(s);
	}
	long operator() () const {
		return lrand48();
	}
};

DeformationGraph::DeformationGraph(DeformableMesh3d& _dmesh, int _nodenum, int _relatenum, int _samplescale, double _deletearearate, 
								   double _wrot, double _wreg, double _wcon, double _wlen)
: dmesh(_dmesh), nodenum(_nodenum), relatenum(_relatenum), samplescale(_samplescale), deletearearate(_deletearearate), 
  wrot(_wrot), wreg(_wreg), wcon(_wcon), wlen(_wlen)
{
	BuildGraph();
}

void DeformationGraph::BuildGraph()
{
	GenerateRandomNodes();
	GenerateEdges();
}

typedef std::pair<DeformationGraph::P3d, int> PIPair;
typedef std::pair<double, int> DIPair;

class DistCmp {
public:
	typedef DeformationGraph::PDPair PDPair;
	bool operator() (const PDPair& lhs, const PDPair& rhs) {
		return lhs.second < rhs.second;
	}
};

class WithinCircle {
public:
	typedef DeformationGraph::P3d P3d;
	WithinCircle(const P3d& p, double r) : pcenter(p), rsquare(r*r) {}
	bool operator() (const DeformationGraph::PointWithID& pointwithid) {
		return (pointwithid.p - pcenter).magsqr() < rsquare;
	}
private:
	P3d pcenter;
	double rsquare;
};

void DeformationGraph::GenerateDensePointSet(std::list<PointWithID>& pointcoll, double delta)
{
#ifdef DEBUGTRACE
	meshtalent::DebugTrace dt("./sample.log");
#endif
	assert(pointcoll.size() == 0);

	typedef DeformableMesh3d::InterMesh InterMesh;
	InterMesh* pMesh = dmesh.pMesh;
	int vertexnum = pMesh->n_vertices();
	//int edgenum = pMesh->n_edges();
	//int facetnum = pMesh->n_faces();
	// sample points on vertices.
	InterMesh::VertexIter v_it, v_end(pMesh->vertices_end());
	for (v_it = pMesh->vertices_begin(); v_it != v_end; ++v_it) {
		InterMesh::Point& p = pMesh->point(v_it);
		P3d p3d(p[0], p[1], p[2]);
		pointcoll.push_back(PointWithID(p3d, PointWithID::VERTEX, v_it.handle().idx()));
	}
#ifdef DEBUGTRACE
	dt.Trace("There %d vertex_points sampled.\n", vertexnum);
#endif
	// sample points on edges.
	int epcount = 0;
	InterMesh::EdgeIter e_it, e_end(pMesh->edges_end());
	for (e_it = pMesh->edges_begin(); e_it != e_end; ++e_it) {
		InterMesh::EdgeHandle eh = pMesh->handle(*e_it);
		InterMesh::HalfedgeHandle heh = pMesh->halfedge_handle(eh, 0); // always select the first edge, but if it's a boundary??
		InterMesh::VertexHandle phfrom = pMesh->from_vertex_handle(heh);
		InterMesh::VertexHandle phto = pMesh->to_vertex_handle(heh);
		InterMesh::Point pf = pMesh->point(phfrom);
		InterMesh::Point pt = pMesh->point(phto);
		P3d p3dfrom(pf[0], pf[1], pf[2]);
		P3d p3dto(pt[0], pt[1], pt[2]);
		double dist = (p3dto - p3dfrom).mag();
		//compute delta vector and sample on this edge.
		V3d deltaV = (p3dto - p3dfrom) * delta;
		int times = static_cast<int>(dist / delta) - 1;
		P3d p3dnow = p3dfrom;
		for (int j = 0; j < times; ++j) {
			p3dnow += deltaV;
			pointcoll.push_back(PointWithID(p3dnow, PointWithID::EDGE, e_it.handle().idx()));
		}
		epcount += (times < 0 ? 0 : times);
	}
#ifdef DEBUGTRACE
	dt.Trace("There are %d edge_points sampled.\n", epcount);
#endif
	// sample points on facet.
	int fpcount = 0;
	InterMesh::FaceIter f_it, f_end(pMesh->faces_end());
	for (f_it = pMesh->faces_begin(); f_it != f_end; ++f_it) {
		P3d trivertices[3];
		//assert((*f_it).is_triangle());
		InterMesh::FaceVertexIter fv_it(pMesh->fv_iter(f_it.handle()));
		int j = 0;
		for (; fv_it; ++fv_it, ++j) {
			InterMesh::Point& p = pMesh->point(fv_it);
			trivertices[j] = P3d(p[0], p[1], p[2]);
		}
		// compute num should be sampled.
		V3d v1 = trivertices[1] - trivertices[0];
		V3d v2 = trivertices[2] - trivertices[0];
		//double lenv1 = v1.mag();
		//double lenv2 = v2.mag();
		double area = (v1 % v2).mag() / 2;
		int samplenum = static_cast<int>(area / (delta * delta));
		fpcount += samplenum;
		RandDoubleGenerator drand;
		drand.srand(time(NULL));
		while (samplenum > 0) {
			double lamada1 = drand(0, 1);
			double lamada2 = drand(0, 1);
			if (lamada1 + lamada2 >= 1) { // judge if this point is not in triangle.
				continue;
			}
			P3d tp3d = trivertices[0] + (v1 * lamada1 + v2 * lamada2);
			pointcoll.push_back(PointWithID(tp3d, PointWithID::FACE, f_it.handle().idx()));
			--samplenum;
		} // end of while.
	} // end of for.
#ifdef DEBUGTRACE
	dt.Trace("There are %d facet_points sampled.\n", fpcount);
	dt.Trace("There are %d points sampled totally.\n", pointcoll.size());
#endif
}

void DeformationGraph::moveFromEdgeToVertex(GraphNode* pnode, const PointWithID& pwid)
{
	typedef DeformableMesh3d::InterMesh InterMesh;
	InterMesh* pMesh = dmesh.pMesh;

	InterMesh::EdgeHandle eh(pwid.id);
	InterMesh::HalfedgeHandle heh = pMesh->halfedge_handle(eh, 0); // always select the first edge, but if it's a boundary??
	InterMesh::VertexHandle phfrom = pMesh->from_vertex_handle(heh);
	InterMesh::VertexHandle phto = pMesh->to_vertex_handle(heh);
	InterMesh::Point pf = pMesh->point(phfrom);
	InterMesh::Point pt = pMesh->point(phto);
	P3d p3dfrom(pf[0], pf[1], pf[2]);
	P3d p3dto(pt[0], pt[1], pt[2]);
	double disToFrom = (p3dfrom - pwid.p).mag();
	double disToTo = (p3dto - pwid.p).mag();

	pnode->vertexID = disToFrom < disToTo ? phfrom.idx() : phto.idx();
	pnode->g = disToFrom < disToTo ? p3dfrom : p3dto;
}

void DeformationGraph::moveFromFaceToVertex(GraphNode* pnode, const PointWithID& pwid)
{
	typedef DeformableMesh3d::InterMesh InterMesh;
	InterMesh* pMesh = dmesh.pMesh;

	InterMesh::FaceHandle fh(pwid.id);
	InterMesh::FaceVertexIter fv_it(pMesh->fv_iter(fh));
	double minDis = 10e8;
	for (; fv_it; ++fv_it) {
		InterMesh::Point& p = pMesh->point(fv_it);
		P3d p3d = P3d(p[0], p[1], p[2]);
		double d = (p3d - pwid.p).mag();
		if (d < minDis) {
			pnode->g = p3d;
			pnode->vertexID = fv_it.handle().idx();
			minDis = d;
		}
	}
}

void DeformationGraph::GenerateRandomNodes()
{
	std::list<PointWithID> pointcoll;
	// compute the sampledelta.
	const double samplelamada = 0.62; // this is a magic num, should not be modified.
	double surfacearea = dmesh.surfaceArea();
	int sizeOfVertices = dmesh.pMesh->n_vertices();
	double sampledelta = sqrt(surfacearea / (samplelamada * samplescale * sizeOfVertices));
	// sample points on surface.
	GenerateDensePointSet(pointcoll, sampledelta);
	// compute the deleteradius.
	const double pi = acos(-1.0);
	const double deleteradius = sqrt(surfacearea / pi / deletearearate);
	// generate nodes.
	edges.reserve(nodenum);
#ifdef DEBUGTRACE
	meshtalent::DebugTrace dt("./pointleft.log");
#endif
	RandUIntGenerator bigrand; 
	bigrand.srand(time(NULL));
	for (int i = 0; i < nodenum; ++i) {
		// adjust nodenum automatically.
		if (pointcoll.empty()) {
			const_cast<int&>(nodenum) = i;
			break;
		}
#ifdef DEBUGTRACE
		dt.Trace("%d\n", pointcoll.size());
#endif
		// randomly select a point.
		int ranindex = bigrand() % static_cast<long>(pointcoll.size());
		std::list<PointWithID>::iterator it = pointcoll.begin();
		advance(it, ranindex);
		// add to Graph.
		GraphNode tnode(it->p, i);
		switch (it->from) {
		case PointWithID::VERTEX:
			tnode.vertexID = it->id;
			break;
		case PointWithID::EDGE:
			moveFromEdgeToVertex(&tnode, *it);
			break;
		case PointWithID::FACE:
			moveFromFaceToVertex(&tnode, *it);
			break;
		default:
			assert(false);
			break;
		}
		edges.push_back(make_pair(tnode, std::vector<int>()));
		pointcoll.remove_if(WithinCircle(it->p, deleteradius)); // remove nearby points.
	}

#ifndef NDEBUG
	// assert the distance between each two Graph nodes is larger than radius.
	for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
		for (int j = i + 1; j < static_cast<int>(edges.size()); ++j) {
			//assert((edges[i].first.g - edges[j].first.g).mag() >= deleteradius);
		}
	}
#endif
	
}

// invoked by GenerateEdges,check if two sorted vector share same element.
template <typename Iter>
bool mergeCheck(Iter lhsBeg, Iter lhsEnd, Iter rhsBeg, Iter rhsEnd)
{
	bool result = false;
	while (lhsBeg != lhsEnd && rhsBeg != rhsEnd) {
		if (*lhsBeg == *rhsBeg) {
			result = true;
			break;
		} else if (*lhsBeg < *rhsBeg) {
			++lhsBeg;
		} else {
			++rhsBeg;
		}
	}
	return result;
}

void DeformationGraph::GenerateEdges()
{
	dmesh.InitDatas(this);
	assert(nodenum == static_cast<int>(edges.size()));
#ifndef NDEBUG
	// assert sorted.
	for (int i = 0; i < nodenum; ++i) {
		std::vector<int>& vs = edges[i].first.vertices;
		for (int j = 0; j < static_cast<int>(vs.size()) - 1; ++j) {
			assert(vs[j] <= vs[j+1]);
		}
	}
	// assert there is no edges information now.
	for (int i = 0; i < nodenum; ++i) {
		assert(edges[i].second.size() == 0);
	}
#endif
	// for each two nodes, check if they has impact on the same vertex.
	for (int i = 0; i < nodenum; ++i) {
		for (int j = i + 1; j < nodenum; ++j) {
			std::vector<int>& vs1 = edges[i].first.vertices;
			std::vector<int>& vs2 = edges[j].first.vertices;
			if (mergeCheck(vs1.begin(), vs1.end(), vs2.begin(), vs2.end())) {
				edges[i].second.push_back(j);
				edges[j].second.push_back(i);
			}
		}
	}
	// sort the indices in edgelist for every node.
	for (int i = 0; i < nodenum; ++i) {
		sort(edges[i].second.begin(), edges[i].second.end());
	}
	// watch how many edges each node has.
	/*
	meshtalent::DebugTrace dt("E:\\howmanyedges.log");
		for (int i = 0; i < nodenum; ++i) {
			dt.Trace("Node num%d: %d edges.\n", i, edges[i].second.size());
		}*/
	
}

void DeformationGraph::reportTox(std::vector<double>& x) const
{
	assert(x.size() == 0);
	// we have a size of ((3*3+3)*nodenum) unknown numbers.
	const int dim = 3;
	x.reserve((dim*dim+dim)*edges.size());
	for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
		// report M3d : R
		for (int j = 0; j < dim; ++j) {
			for (int k = 0; k < dim; ++k) {
				x.push_back(edges[i].first.R[j][k]);
			}
		}
		// report V3d : t
		for (int j = 0; j < dim; ++j) {
			x.push_back(edges[i].first.t[j]);
		}
	}
	assert(x.size() == (dim*dim+dim)*edges.size());
}

void DeformationGraph::getFromx(const std::vector<double>& x)
{
	const int dim = 3;
	assert(x.size() == (dim*dim+dim)*edges.size());
	int index = 0;
	for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
		// get M3d : R
		for (int j = 0; j < dim; ++j) {
			for (int k = 0; k < dim; ++k) {
				edges[i].first.R[j][k] = x[index++];
			}
		}
		// get V3d : t
		for (int j = 0; j < dim; ++j) {
			edges[i].first.t[j] = x[index++];
		}
	}
	assert(index == static_cast<int>(x.size()));
}

void DeformationGraph::updateNodesPos()
{
	assert(static_cast<int>(edges.size()) == nodenum);
	for (int i = 0; i < nodenum; ++i) {
		edges[i].first.g += edges[i].first.t;
	}
}

} // end of meshtalent
