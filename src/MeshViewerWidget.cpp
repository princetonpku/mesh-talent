#include <iostream>

#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "MeshViewerWidget.h"

using namespace OpenMesh;

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
 : QGLViewerWidget(parent), pdmesh_(NULL), showlist_start_(0)
{
}

MeshViewerWidget::~MeshViewerWidget()
{
	delete pdmesh_;
}

bool MeshViewerWidget::openMesh(const char* filename)
{
	mesh_.request_face_normals();
	//mesh_.request_face_colors();
	mesh_.request_vertex_normals();
	//mesh_.request_vertex_colors();

	if (OpenMesh::IO::read_mesh(mesh_, filename)) {
		mesh_.update_face_normals();
		mesh_.update_vertex_normals();

		// bounding box
		InterMesh::ConstVertexIter vIt(mesh_.vertices_begin());
		InterMesh::ConstVertexIter vEnd(mesh_.vertices_end());      
		
		typedef InterMesh::Point Point;
		using OpenMesh::Vec3d;
		
		Vec3d bbMin, bbMax;
		
		bbMin = bbMax = OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt));
		
		for (size_t count = 0; vIt != vEnd; ++vIt, ++count) {
			bbMin.minimize(OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt)));
			bbMax.maximize(OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt)));
		}
    
		// set center and radius
		set_scene_pos( (bbMin+bbMax)*0.5, (bbMin-bbMax).norm()*0.5 );

		// create show list
		init_showlist();

#if defined(WIN32)
		updateGL();
#endif

		setWindowTitle(QFileInfo(filename).fileName());

		// loading done
		return true;
	}
	return false;
}

void MeshViewerWidget::init_showlist()
{
	showlist_start_ = glGenLists(N_DRAW_MODES);
	
	// wireframe.
	glNewList(showlist_start_ + WIRE_FRAME, GL_COMPILE);
	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	draw_mesh_wireframe();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEndList();

	// solidflat.
	glNewList(showlist_start_ + SOLID_FLAT, GL_COMPILE);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);
	draw_mesh_solidflat();
	glEndList();

	// solidsmooth.
	glNewList(showlist_start_ + SOLID_SMOOTH, GL_COMPILE);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	draw_mesh_solidsmooth();
	glEndList();

	// pointset.
	glNewList(showlist_start_ + POINT_SET, GL_COMPILE);
	glDisable(GL_LIGHTING);
	draw_mesh_pointset();
	glEndList();
}

void MeshViewerWidget::open_mesh_gui(QString fname)
{
	if ( fname.isEmpty() || !openMesh(fname.toLocal8Bit()) )
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical( NULL, windowTitle(), msg);
	}
}

void MeshViewerWidget::draw_scene(int drawmode)
{
	if (!mesh_.n_vertices()) { return; }
	assert(drawmode < N_DRAW_MODES);
	glCallList(showlist_start_ + drawmode);
}

void MeshViewerWidget::draw_mesh_wireframe() const
{
	InterMesh::ConstFaceIter fIt(mesh_.faces_begin()),
							 fEnd(mesh_.faces_end());
	InterMesh::ConstFaceVertexIter fvIt;

	glBegin(GL_TRIANGLES);
	for (; fIt != fEnd; ++fIt) {
		fvIt = mesh_.cfv_iter(fIt.handle()); 
		glVertex3dv(&mesh_.point(fvIt)[0]);
		++fvIt;
		glVertex3dv(&mesh_.point(fvIt)[0]);
		++fvIt;
		glVertex3dv(&mesh_.point(fvIt)[0]);
	}
	glEnd();
}

void MeshViewerWidget::draw_mesh_solidflat() const
{
	InterMesh::ConstFaceIter fIt(mesh_.faces_begin()),
							 fEnd(mesh_.faces_end());
	InterMesh::ConstFaceVertexIter fvIt;

	glBegin(GL_TRIANGLES);
	for (; fIt != fEnd; ++fIt) {
		glNormal3dv(&mesh_.normal(fIt)[0]);

		fvIt = mesh_.cfv_iter(fIt.handle());
		glVertex3dv(&mesh_.point(fvIt)[0]);
		++fvIt;
		glVertex3dv(&mesh_.point(fvIt)[0]);
		++fvIt;
		glVertex3dv(&mesh_.point(fvIt)[0]);
	}
	glEnd();
}

void MeshViewerWidget::draw_mesh_solidsmooth() const 
{
	InterMesh::ConstFaceIter fIt(mesh_.faces_begin()),
							 fEnd(mesh_.faces_end());
	InterMesh::ConstFaceVertexIter fvIt;

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh_.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh_.vertex_normals());

	glBegin(GL_TRIANGLES);
	for (; fIt != fEnd; ++fIt) {
		fvIt = mesh_.cfv_iter(fIt.handle());
		glArrayElement(fvIt.handle().idx());
		++fvIt;
		glArrayElement(fvIt.handle().idx());
		++fvIt;
		glArrayElement(fvIt.handle().idx());
	}
	glEnd();

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewerWidget::draw_mesh_pointset() const 
{
	InterMesh::ConstFaceIter fIt(mesh_.faces_begin()),
							 fEnd(mesh_.faces_end());
	InterMesh::ConstFaceVertexIter fvIt;

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh_.points());
	glDrawArrays(GL_POINTS, 0, mesh_.n_vertices());
	glDisableClientState(GL_VERTEX_ARRAY);
}
