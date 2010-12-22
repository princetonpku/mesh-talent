#ifndef DEFORMABLEMESH3D_H
#define DEFORMABLEMESH3D_H

#include <vector>
#include <utility>
#include <set>
// ---------------------OpenMesh
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
// ---------------------math utility
#include "math/vector3d.h"
#include "math/matrix3d.h"
#include "math/point3d.h"
// ---------------------

namespace meshtalent {

class DeformationGraph;
class GaussNewtonSolver;

class DeformableMesh3d {
public:
	friend class DeformationGraph;
	friend class GaussNewtonSolver;
public:
	typedef math::Vector3d<double> V3d;
	typedef math::Point3d<double> P3d;
	typedef math::Matrix3d<double> M3d;
	typedef std::pair<int, double> IndexWeightPair; // index of GraphNode, and weight of the impact of deformation.
public:
	struct InterTraits : public OpenMesh::DefaultTraits {
		typedef OpenMesh::Vec3d Point;
		typedef OpenMesh::Vec3d Normal;
	};
	typedef OpenMesh::TriMesh_ArrayKernelT<InterTraits> InterMesh;
public:
	// constructor.
	DeformableMesh3d(InterMesh* _pMesh);
	// destructor.
	~DeformableMesh3d();

	InterMesh* getpMesh() const {
		return pMesh;
	}
	
	// make all the selected vertex translate by t.
	void translate(const std::vector<int>& selectedHandles, const V3d& t);
	// make all the selected vertex scaling by their center.
	void scale(const std::vector<int>& selectedHandles, double lamada);
	// make all the selected vertex rotate round a given vector starting at a given point, the angle is beta.
	void rotate(const std::vector<int>& selectedHandles, const V3d& t, const P3d& p, double beta);

	void gethandles();
	size_t handleSize() const {
		return handleIDs.size();
	}
	// bounding box volume & surface area.
public:
	double boundingboxVolume() const;
	double surfaceArea() const;
public:
	P3d centerOfSV() const; // compute every time.
private:
	// should be called
	void InitDatas(DeformationGraph* _pGraph); // should be invoked before the edges of graph were generated.
	void setGraph(DeformationGraph* _pGraph) {
		pGraph = _pGraph;
	}
private:
	template <typename FUN>
	void specifyposArr(const std::vector<int>& selectedHandles, FUN fun); // fill specifiedposArr with the pos that user specified.
	class translateFun {
	public:
		translateFun(const V3d& _t) : t(_t) {}
		const P3d operator() (const P3d& p) {
			return p + t;
		}
	private:
		const V3d t;
	};
	class scaleFun {
	public:
		scaleFun(const P3d& _pcenter, double _lamada) : pcenter(_pcenter), lamada(_lamada) {}
		const P3d operator() (const P3d& p) {
			return pcenter + (p - pcenter) * lamada;
		}
	private:
		const P3d pcenter;
		const double lamada;
	};
	class rotateFun {
	public:
		rotateFun(const P3d& _pcenter, const V3d& _u, double _beta) : pcenter(_pcenter), u(_u), beta(_beta) {
			assert(u.mag() == 1.0);
			const double pi = acos(-1.0);
			double radianbeta = beta * pi / 180.0;
			double c = cos(radianbeta);
			double s = sin(radianbeta);
			// compute rotate matrix.
			Ru[0][0] = c + (1-c) * u[0] * u[0];
			Ru[0][1] = (1-c) * u[1] * u[0] - s * u[2];
			Ru[0][2] = (1-c) * u[2] * u[0] + s * u[1];
			Ru[1][0] = (1-c) * u[0] * u[1] + s * u[2];
			Ru[1][1] = c + (1-c) * u[1] * u[1];
			Ru[1][2] = (1-c) * u[2] * u[1] - s * u[0];
			Ru[2][0] = (1-c) * u[0] * u[2] - s * u[1];
			Ru[2][1] = (1-c) * u[1] * u[2] + s * u[0];
			Ru[2][2] = c + (1-c) * u[2] * u[2];
#ifndef NDEBUG
			{
				// check if this matrix is a rotate matrix.
				for (int i = 0; i < 3; ++i) {
					assert(Ru[i].mag() == 1.0);
				}
				assert(Ru.determinant() == 1.0);
			}
#endif
		}
		const P3d operator() (const P3d& p) {
			return pcenter + (Ru * (p - pcenter));
		}
	private:
		const P3d pcenter;
		const V3d u;
		double beta;
		M3d Ru;
	};
	void deform(); // invoked by translate, scale, rotate.
	void updateVerticesUsingGraph(); // invoked by deform.
private:
	InterMesh* pMesh;
	DeformationGraph* pGraph;
	// nearnodesArr's every element is correspond to pMesh's Vertex of the same index.
	// std::vector<IndexWeightPair> is a coll of nodes that has impact of a Vertex, always has size of "pGraph->relatenum".
	OpenMesh::VPropHandleT<std::vector<IndexWeightPair> > nearnodesArr;
	// handleIDs's elements tell which vertices are specified as control handles.
	std::vector<InterMesh::VertexHandle> handleIDs;
	/* specifiedposArr's every element is correspond to the vertex with a index in handleIDs. */
	std::vector<P3d> specifiedposArr;
public:
	std::vector<InterMesh::VertexHandle>& getHandleIDs() {
		return handleIDs;
	}
	const std::vector<InterMesh::VertexHandle>& getHandleIDs() const {
		return handleIDs;
	}
private:
	// Every Vertex which is selected, its handle will be stored in selectedVertices.
	OpenMesh::MPropHandleT<std::set<InterMesh::VertexHandle> > selectedVertices;
public:
	std::set<InterMesh::VertexHandle>& getSVSet() {
		return pMesh->property(selectedVertices);
	}
	const std::set<InterMesh::VertexHandle>& getSVSet() const {
		return pMesh->property(selectedVertices);
	}

// for computing voronoi of the mesh.
public:
	struct VoronoiInfo {
		double dis;
		InterMesh::VertexHandle handleRoot;
		VoronoiInfo(double _dis = 10e8, InterMesh::VertexHandle _handleRoot = InterMesh::VertexHandle()) : dis(_dis), handleRoot(_handleRoot) {}
		struct DisCmp {
			bool operator() (const VoronoiInfo& lhs, const VoronoiInfo& rhs) {
				return lhs.dis < rhs.dis;
			}
		};
	};
private:
	OpenMesh::VPropHandleT<VoronoiInfo> voroInfo;
	OpenMesh::HPropHandleT<double> halfEdgeDis;
private:
	void genVoronoi(); // called by InitDatas().

// for computing overlapped voronoi.
public:
	struct OverlappedVoronoiInfo {
	public:
		static const int MAXSIZE = 8;
	public:
		VoronoiInfo vis[MAXSIZE]; // always ordered by it's dis, and handleRoot are unique.
		int indexNow; // have computed out how many nearest nodes.
		OverlappedVoronoiInfo() : indexNow(0) {}
	public:
		static int getNearNodesNum() { return nearNodesNum; }
		static void setNearNodesNum(int k) { nearNodesNum = k; }
	public:
		static int nearNodesNum; // find nearNodesNum nearest nodes for every vertex, no more than MAXSIZE.
	};
private:
	OpenMesh::VPropHandleT<OverlappedVoronoiInfo> overlapped_voroInfo;
private:
	void genOverlappedVoronoi();
	bool overlappedVoronoiRelax(const OverlappedVoronoiInfo& OLVIParent, OverlappedVoronoiInfo& OLVINow, double disOfEdge);
};

} // end of namespace meshtalent

#endif // DEFORMABLEMESH3D_H
