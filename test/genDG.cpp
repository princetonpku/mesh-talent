#include <iostream>
#include <fstream>
#include <cstdlib>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

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

	std::ofstream DGInfo(argv[2]);
	DGInfo << *pDG << std::endl;
	return EXIT_SUCCESS;
}
