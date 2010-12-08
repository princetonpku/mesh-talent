#ifndef DEFORMATIONGRAPH_H
#define DEFORMATIONGRAPH_H

#include <utility>
#include <vector>
#include <list>
#include <algorithm>
#include <ostream>
#include <istream>
// --------------------- math utility
#include "math/vector3d.h"
#include "math/matrix3d.h"
#include "math/point3d.h"
// ---------------------

namespace meshtalent {

class DeformableMesh3d;
class GaussNewtonSolver;

class DeformationGraph {
public:
	friend class DeformableMesh3d;
	friend class GaussNewtonSolver;
public:
	typedef math::Vector3d<double> V3d;
	typedef math::Matrix3d<double> M3d;
	typedef math::Point3d<double> P3d;
	typedef std::pair<P3d, double> PDPair;
public:
	DeformationGraph(DeformableMesh3d& _dmesh, int _nodenum, int _relatenum, int _samplescale, double _deleteradius, 
					 double _wrot, double _wreg, double _wcon, double _wlen);
	// for numeric solve.
	void reportTox(std::vector<double>& x) const;
	void getFromx(const std::vector<double>& x);
	// graph nodes should be updated after updating mesh vertices's pos.
	void updateNodesPos();
public:
	friend std::ostream& operator<< (std::ostream& os, const DeformationGraph& dg) {
		int nodesize = static_cast<int>(dg.edges.size());
		int edgesize = 0;
		for (int i = 0; i < nodesize; ++i) {
			edgesize += static_cast<int>(dg.edges[i].second.size());
		}
		os << nodesize << ' ' << edgesize << '\n';
		for (int i = 0; i < nodesize; ++i) {
			const P3d& g = dg.edges[i].first.g;
			os << g[0] << ' ' << g[1] << ' ' << g[2] << '\n';
		}
		for (int i = 0; i < nodesize; ++i) {
			for (int j = 0; j < static_cast<int>(dg.edges[i].second.size()); ++j) {
				os << i << ' ' << dg.edges[i].second[j] << '\n';
			}
		}
		return os;
	}
public:
	struct PointWithID {
		P3d p;
		enum From { VERTEX, EDGE, FACE };
		enum From from;
		int id;
		PointWithID(const P3d& _p, enum From _from, int _id) : p(_p), from(_from), id(_id) {}
	};
public:
	struct GraphNode {
		GraphNode(const P3d& _g, int _index, int _vertexID = -1) : g(_g), index(_index), vertexID(_vertexID) {}
		P3d g; // the space coordinate of this graph node.
		V3d t; // the translation transformation vector.
		M3d R; // the linear transformation matrix.
		int index; // the index in "edges".
		std::vector<int> vertices; // the coll of vertex whose deformation is correspond to this node.
		int vertexID; // the id of vertex it belongs.
	};
	typedef std::pair<GraphNode, std::vector<int> > Link; // node and adjacent edges
private:
	void BuildGraph();
	// called by GenerateRandomNodes().
	void moveFromEdgeToVertex(GraphNode* pnode, const PointWithID& pwid);
	void moveFromFaceToVertex(GraphNode* pnode, const PointWithID& pwid);
	void GenerateRandomNodes();
	void GenerateEdges();
	void GenerateDensePointSet(std::list<PointWithID>& pointcoll, double delta);// sample points on the surface.

private:
	DeformableMesh3d& dmesh;
	const int nodenum; // the num of nodes of this graph
	const int relatenum; // the num of related nodes of each mesh vertex
	const int samplescale;
	const double deletearearate;
	std::vector<Link> edges;
public:
	std::vector<Link>& getEdges() { return edges; }
	const std::vector<Link>& getEdges() const { return edges; }
private:
	const double wrot;
	const double wreg;
	const double wcon;
	const double wlen;
};

} // end of namespace meshtalent

#endif // DEFORMATIONGRAPH_H
