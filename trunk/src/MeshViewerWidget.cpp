#include <iostream>
#include <algorithm>
#include <cmath>

#include <qapplication.h>
#include <QMouseEvent>
#include <QLineEdit>

#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "MeshViewerWidget.h"
#include "genGraphDialog.h"
#include "DeformationGraph.h"

using namespace OpenMesh;
using namespace Qt;

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
 : QGLViewerWidget(parent), pdmesh_(NULL), pdgraph_(NULL)
{
}

MeshViewerWidget::~MeshViewerWidget()
{
	delete pdmesh_;
}

bool MeshViewerWidget::openMesh(const char* filename)
{
	mesh_.request_face_normals();
	mesh_.request_face_colors();
	mesh_.request_vertex_normals();
	mesh_.request_vertex_colors();

	if (OpenMesh::IO::read_mesh(mesh_, filename)) {
		
		mesh_.update_face_normals();
		mesh_.update_vertex_normals();

		InterMesh::VertexIter vIt = mesh_.vertices_begin();
		InterMesh::VertexIter vEnd = mesh_.vertices_end();
		for (; vIt != vEnd; ++vIt) {
			mesh_.set_color(vIt, InterMesh::Color(127, 127, 0));
		}

		// bounding box
		typedef InterMesh::Point Point;
		using OpenMesh::Vec3d;
		
		Vec3d bbMin, bbMax;
		
		vIt = mesh_.vertices_begin();
		bbMin = bbMax = OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt));
		
		size_t count = 0;
		for (; vIt != vEnd; ++vIt, ++count) {
			bbMin.minimize(OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt)));
			bbMax.maximize(OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt)));
		}
    
		// set center and radius
		set_scene_pos( (bbMin+bbMax)*0.5, (bbMin-bbMax).norm()*0.5 );

#if defined(WIN32)
		updateGL();
#endif

		// set deformable mesh.
		pdmesh_ = new DeformableMesh3d(&mesh_);

		setWindowTitle(QFileInfo(filename).fileName());

		// loading done
		return true;
	}
	return false;
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
	//assert(drawmode < N_DRAW_MODES);
	switch (drawmode) {
	case WIRE_FRAME:
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_mesh_wireframe();
		draw_mesh_pointset();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case SOLID_FLAT:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		draw_mesh_pointset();
		break;
	case SOLID_SMOOTH:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		draw_mesh_solidsmooth();
		draw_mesh_pointset();
		break;
	case POINT_SET:
		glDisable(GL_LIGHTING);
		draw_mesh_pointset();
		break;
	case DRAW_GRAPH:
		glDisable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_graph();
		break;
	default:
		break;
	}
	draw_select_boxes();
}

void MeshViewerWidget::draw_mesh_wireframe() const
{
	InterMesh::ConstFaceIter fIt(mesh_.faces_begin()),
							 fEnd(mesh_.faces_end());
	InterMesh::ConstFaceVertexIter fvIt;

	glLoadName(mesh_.n_vertices());

	glColor3f(1.0, 1.0, 0.0);
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

	glLoadName(mesh_.n_vertices());

	glBegin(GL_TRIANGLES);
	for (int i = 0; fIt != fEnd; ++fIt, ++i) {
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

	glLoadName(mesh_.n_vertices());

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
	InterMesh::ConstVertexIter vIt(mesh_.vertices_begin()),
							vEnd(mesh_.vertices_end());
	for (int i = 0; vIt != vEnd; ++i, ++vIt) {
		glLoadName(i);
		glBegin(GL_POINTS);
		glColor3bv((const GLbyte*)&mesh_.color(vIt)[0]);
		glVertex3dv(&mesh_.point(vIt)[0]);
		glEnd();
	}
}

void MeshViewerWidget::draw_graph() const 
{
	assert(pdgraph_ != NULL);

	const std::vector<DeformationGraph::Link>& edges = pdgraph_->getEdges();
	const int size = edges.size();
	// graph nodes.
	glBegin(GL_POINTS);
	for (int i = 0; i < size; ++i) {
		glVertex3dv(&edges[i].first.g[0]);
	}
	glEnd();
	// graph edges.
	glBegin(GL_LINES);
	for (int i = 0; i < size; ++i) {
		for (size_t j = 0; j < edges[i].second.size(); ++j) {
			glVertex3dv(&edges[i].first.g[0]);
			glVertex3dv(&edges[edges[i].second[j]].first.g[0]);
		}
	}
	glEnd();
}
void MeshViewerWidget::mousePressEvent(QMouseEvent* _event)
{
	if (mesh_.n_vertices() == 0) return;
	// rotate, translate, scale, pick, or deform.
	if (_event->buttons() == LeftButton) {
		int mm = mouse_mode();
		if (mm < N_MOUSE_MODES) { // rotate, translate, scale.
			QGLViewerWidget::mousePressEvent(_event);
		} else {                  // pick, deform.
			switch (mm) {
			case MOUSE_PICK:
				// record the press pos.
				pick_press_point_.rx() = _event->x();
				pick_press_point_.ry() = _event->y();
				break;
			case MOUSE_DEFORM:
				break;
			}
		} // end of else
	} // end of if
}

void MeshViewerWidget::mouseMoveEvent(QMouseEvent* _event)
{
	if (mesh_.n_vertices() == 0) return;
	// rotate, translate, scale, pick, or deform.
	if (_event->buttons() == LeftButton) {
		int mm = mouse_mode();
		if (mm < N_MOUSE_MODES) { // rotate, translate, scale.
			QGLViewerWidget::mouseMoveEvent(_event);
		} else {                  // pick, deform.
			switch (mm) {
			case MOUSE_PICK:
				break;
			case MOUSE_DEFORM:
				break;
			}
		} // end of else
	} // end of if
}

void MeshViewerWidget::mouseReleaseEvent(QMouseEvent* _event)
{
	if (mesh_.n_vertices() == 0) return;
	// rotate, translate, scale, pick, or deform.
	int mm = mouse_mode();
	if (mm < N_MOUSE_MODES) { // rotate, translate, scale.
		QGLViewerWidget::mouseReleaseEvent(_event);
	} else {                  // pick, deform.
		switch (mm) {
		case MOUSE_PICK:
			processPick(_event);
			break;
		case MOUSE_DEFORM:
			break;
		}
	} // end of else
}

void MeshViewerWidget::processHits(GLint hits, GLuint* buffer, bool controled)
{
	int n_vertices = mesh_.n_vertices();
	if (!n_vertices) return;

	if (!controled) {
		pdmesh_->getSVSet().clear();
	}

	GLuint* ptr = buffer;
	std::cout << "hits = " << hits << std::endl;
	for (int i = 0; i < hits; ++i) {
		int names = *ptr;
		assert(names == 1);
		*ptr++;
		std::cout << "   z1 is " << *ptr;
		*ptr++;
		std::cout << "   z2 is " << *ptr;
		*ptr++;
		std::cout << "   index is " << *ptr << std::endl;;
		assert(*ptr <= n_vertices);
		if ((int)(*ptr) < n_vertices) {
			// set vertex color
			//mesh_.set_color(InterMesh::VertexHandle(*ptr), InterMesh::Color(0, 0, 127));
			pdmesh_->getSVSet().insert(InterMesh::VertexHandle(*ptr));
		}
		*ptr++;
	}
}

void MeshViewerWidget::processPick(QMouseEvent* _event)
{
	std::vector<GLuint> selectBuf(mesh_.n_vertices()); // assert enough space.
	GLint hits;
	GLint viewport[4];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glSelectBuffer(mesh_.n_vertices(), &selectBuf[0]);
	glRenderMode(GL_SELECT);
	
	glInitNames();
	glPushName(0);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	if (_event->x() == pick_press_point_.x() && _event->y() == pick_press_point_.y()) { // point select.
		gluPickMatrix((GLdouble)(_event->x()), (GLdouble)(viewport[3]-_event->y()), 5.0, 5.0, viewport);
	} else { // rect select.
		// compute the pick rect.
		GLdouble x1 = (GLdouble)(_event->x());
		GLdouble y1 = (GLdouble)(viewport[3]-(_event->y()));
		GLdouble x2 = (GLdouble)(pick_press_point_.x());
		GLdouble y2 = (GLdouble)(viewport[3]-(pick_press_point_.y()));
		GLdouble detX = fabs(x1-x2);
		GLdouble detY = fabs(y1-y2);
		GLdouble cenX = (x1+x2)/2;
		GLdouble cenY = (y1+y2)/2;
		std::cout << "The rect is " << cenX << " " << cenY << " " << detX << " " << detY << std::endl;
		gluPickMatrix(cenX, cenY, detX, detY, viewport);
	}
	glMultMatrixd(projection_matrix());
	draw_scene(draw_mode());

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glFlush();

	hits = glRenderMode(GL_RENDER);
	processHits(hits, &selectBuf[0], _event->modifiers() == ControlModifier);
	updateGL();
}

void MeshViewerWidget::drawSelectBox(const Vec3d& center, double radius)
{
	double xmin, ymin, zmin, xmax, ymax, zmax;
	xmin = center[0] - radius;
	ymin = center[1] - radius;
	zmin = center[2] - radius;
	xmax = center[0] + radius;
	ymax = center[1] + radius;
	zmax = center[2] + radius;
	glBegin(GL_QUADS);

	glVertex3d(xmin, ymin, zmin);
	glVertex3d(xmin, ymin, zmax);
	glVertex3d(xmax, ymin, zmax);
	glVertex3d(xmax, ymin, zmin);
	
	glVertex3d(xmin, ymax, zmin);
	glVertex3d(xmin, ymax, zmax);
	glVertex3d(xmax, ymax, zmax);
	glVertex3d(xmax, ymax, zmin);

	glVertex3d(xmin, ymin, zmin);
	glVertex3d(xmin, ymax, zmin);
	glVertex3d(xmin, ymax, zmax);
	glVertex3d(xmin, ymin, zmax);

	glVertex3d(xmax, ymin, zmin);
	glVertex3d(xmax, ymax, zmin);
	glVertex3d(xmax, ymax, zmax);
	glVertex3d(xmax, ymin, zmax);

	glVertex3d(xmin, ymin, zmin);
	glVertex3d(xmin, ymax, zmin);
	glVertex3d(xmax, ymax, zmin);
	glVertex3d(xmax, ymin, zmin);

	glVertex3d(xmin, ymin, zmax);
	glVertex3d(xmin, ymax, zmax);
	glVertex3d(xmax, ymax, zmax);
	glVertex3d(xmax, ymin, zmax);

	glEnd();
}

void MeshViewerWidget::draw_select_boxes()
{
	const std::set<InterMesh::VertexHandle>& selectedVertices = 
		pdmesh_->getSVSet();
	glColor3f(1.0, 0.0, 1.0);
	typedef std::set<InterMesh::VertexHandle>::const_iterator CITER;
	for (CITER cit = selectedVertices.begin(); cit != selectedVertices.end(); ++cit) {
		Vec3d center = mesh_.point(*cit);
		double radius = 0.5;
		drawSelectBox(center, radius);
	}
}

void MeshViewerWidget::gen_graph_query()
{
	if (pdgraph_ || !mesh_.n_vertices()) { 
		QMessageBox(QMessageBox::NoIcon, "Alert", "Condition not matched");
		return;
   	}

	GenGraphDialog dialog(this);
	if (dialog.exec()) {
		pdgraph_ = new DeformationGraph(*pdmesh_, dialog.nodeNumLineEdit->text().toDouble(), 
				dialog.relateNumLineEdit->text().toDouble(), 
				dialog.sampleScaleLineEdit->text().toDouble(), 
				dialog.deleteRadiusLineEdit->text().toDouble(), 
				dialog.ErotLineEdit->text().toDouble(), 
				dialog.EregLineEdit->text().toDouble(), 
				dialog.EconLineEdit->text().toDouble(),
				20);
	}
}
