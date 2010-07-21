#include <iostream>
#include <fstream>
#include <cstdlib>
// ---------------------OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
// ---------------------meshtalent
#include "DeformableMesh3d.h"
#include "DeformationGraph.h"

int main(int argc, char** argv)
{
	using namespace meshtalent;
	typedef DeformableMesh3d::InterMesh MyMesh;
	MyMesh mesh;

	if (argc != 3) {
		std::cerr << "Usage: " << argv[0] << " inputfile outputfile\n";
		return EXIT_FAILURE;
	}

	if (!OpenMesh::IO::read_mesh(mesh, argv[1])) {
		std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
		return EXIT_FAILURE;
	}

	DeformableMesh3d* pDM = new DeformableMesh3d(&mesh);
	DeformationGraph* pDG = new DeformationGraph(*pDM, 300, 3, 20, 150, 1, 10, 100, 20);

	typedef DeformableMesh3d::InterMesh::VertexHandle VH;
	std::vector<VH>& svs = pDM->getSVArr();
	svs.push_back(VH(148));
	svs.push_back(VH(1268));
	svs.push_back(VH(2196));
	svs.push_back(VH(2249));
	svs.push_back(VH(3128));
	pDM->gethandles();

	svs.clear();
	svs.push_back(VH(148));
	typedef DeformableMesh3d::V3d V3d;
	V3d v(20, 0, -15);
	pDM->translate(v);

	delete pDM;
	delete pDG;
	if (!OpenMesh::IO::write_mesh(mesh, argv[2])) {
		std::cerr << "Error: cannot write mesh to " << argv[2] << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
