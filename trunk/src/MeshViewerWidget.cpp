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
 : QGLViewerWidget(parent), pdmesh_(NULL), pdgraph_(NULL), box_radius_(0.0)
{
}

MeshViewerWidget::~MeshViewerWidget()
{
	delete pdmesh_;
	delete pdgraph_;
}

void MeshViewerWidget::updateMeshCenter()
{
	typedef InterMesh::Point Point;
	using OpenMesh::Vec3d;
	Vec3d bbMin, bbMax;

	InterMesh::VertexIter vIt = mesh_.vertices_begin();
	InterMesh::VertexIter vEnd = mesh_.vertices_end();
	bbMin = bbMax = OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt));

	size_t count = 0;
	for (; vIt != vEnd; ++vIt, ++count) {
		bbMin.minimize(OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt)));
		bbMax.maximize(OpenMesh::vector_cast<Vec3d>(mesh_.point(vIt)));
	}

	// set center and radius and box's radius.
	set_scene_pos( (bbMin+bbMax)*0.5, (bbMin-bbMax).norm()*0.5 );
	box_radius_ = 0.005 * (bbMin-bbMax).norm();
}

void MeshViewerWidget::updateMeshNormals()
{
	mesh_.update_face_normals();
	mesh_.update_vertex_normals();
}

bool MeshViewerWidget::openMesh(const char* filename)
{
	mesh_.request_face_normals();
	mesh_.request_face_colors();
	mesh_.request_vertex_normals();
	mesh_.request_vertex_colors();

	if (OpenMesh::IO::read_mesh(mesh_, filename)) {

		InterMesh::VertexIter vIt = mesh_.vertices_begin();
		InterMesh::VertexIter vEnd = mesh_.vertices_end();
		for (; vIt != vEnd; ++vIt) {
			mesh_.set_color(vIt, InterMesh::Color(127, 127, 0));
		}

		update_mesh();

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
	if (fname.isEmpty() || !openMesh(fname.toLocal8Bit())) {
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}

void MeshViewerWidget::save_mesh_gui(QString fname)
{
	if (fname.isEmpty() || !saveMesh(fname.toLocal8Bit())) {
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
	
}

bool MeshViewerWidget::saveMesh(const char* filename)
{
	return OpenMesh::IO::write_mesh(mesh_, filename);
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
		draw_select_boxes();
		break;
	case SOLID_FLAT:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		draw_mesh_pointset();
		glDisable(GL_LIGHTING);
		draw_select_boxes();
		break;
	case SOLID_SMOOTH:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		draw_mesh_solidsmooth();
		draw_mesh_pointset();
		glDisable(GL_LIGHTING);
		draw_select_boxes();
		break;
	case POINT_SET:
		glDisable(GL_LIGHTING);
		draw_mesh_pointset();
		draw_select_boxes();
		break;
	case VORONOI_DIAGRAM:
		glDisable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_voronoidiagram();
		break;
	case DRAW_GRAPH:
		glDisable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_graph();
		break;
	default:
		break;
	}
}

void MeshViewerWidget::draw_mesh_wireframe() const
{
	InterMesh::ConstFaceIter fIt(mesh_.faces_begin()),
							 fEnd(mesh_.faces_end());
	InterMesh::ConstFaceVertexIter fvIt;

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
		if (mouse_mode() == MOUSE_PICK) {
			glLoadName(i);
		}
		glBegin(GL_POINTS);
		glColor3bv((const GLbyte*)&mesh_.color(vIt)[0]);
		glNormal3dv(&mesh_.normal(vIt)[0]);
		glVertex3dv(&mesh_.point(vIt)[0]);
		glEnd();
	}
	if (mouse_mode() == MOUSE_PICK) {
		glLoadName(mesh_.n_vertices());
	}
	
}
void MeshViewerWidget::draw_voronoidiagram() const
{
	InterMesh::ConstFaceIter fIt(mesh_.faces_begin()),
							fEnd(mesh_.faces_end());
	for (; fIt != fEnd; ++fIt) {
		InterMesh::Point endpoints[3];
		InterMesh::Point midpoints[3];
		InterMesh::Point center;
		InterMesh::Color colors[3];
		//assert((*f_it).is_triangle());
		InterMesh::ConstFaceVertexIter fv_it(mesh_, fIt.handle());
		int j = 0;
		for (; fv_it; ++fv_it, ++j) {
			endpoints[j] = mesh_.point(fv_it);
			colors[j] = mesh_.color(fv_it);
		}
		for (j = 0; j < 3; ++j) {
			midpoints[j] = (endpoints[(j+1)%3] + endpoints[(j+2)%3]) / 2;
		}
		center = (endpoints[0] + endpoints[1] + endpoints[2]) / 3;

		for (j = 0; j < 3; ++j) {
			glColor3bv((const GLbyte*)&colors[j]);
			glBegin(GL_QUADS);
				glVertex3dv(&endpoints[j][0]);
				glVertex3dv(&midpoints[(j+1)%3][0]);
				glVertex3dv(&center[0]);
				glVertex3dv(&midpoints[(j+2)%3][0]);
			glEnd();
		}
	}
}

void MeshViewerWidget::draw_graph() const 
{
	assert(pdgraph_ != NULL);

	const std::vector<DeformationGraph::Link>& edges = pdgraph_->getEdges();
	const int size = edges.size();

	glColor3f(0.0, 1.0, 0.0);
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

	const int BUFFERSIZE = 2048;
	char buffer[BUFFERSIZE];
	char *pBuf = buffer;
	char *pBufEnd = buffer + BUFFERSIZE;
	buffer[0] = '\0';

	if (_event->buttons() == RightButton) {
		for (std::set<InterMesh::VertexHandle>::iterator it = pdmesh_->getSVSet().begin(); 
			it != pdmesh_->getSVSet().end(); ++it) {
				sprintf(pBuf, "%d %lf %lf %lf\n%", it->idx(), mesh_.point(*it)[0],
					mesh_.point(*it)[1], mesh_.point(*it)[2]);
				int len = strlen(pBuf);
				pBuf += len;
				if (pBuf + len >= pBufEnd) {
					break;
				}
		}
		QMessageBox::information(this, "Selected Vertices", buffer);
	}
	

	// rotate, translate, scale, pick, or deform.
	if (_event->buttons() == LeftButton) {
		int mm = mouse_mode();
		if (mm < N_MOUSE_MODES) { // rotate, translate, scale.
			QGLViewerWidget::mousePressEvent(_event);
		} else {                  // pick, deform.
			switch (mm) {
			case MOUSE_PICK:
				processMousePickPress(_event);
				break;
			case MOUSE_DEFORM:
				processMouseDeformPress(_event);
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
				processMousePickMove(_event);
				break;
			case MOUSE_DEFORM:
				processMouseDeformMove(_event);
				break;
			}
		} // end of else
	} // end of if
}

void MeshViewerWidget::mouseReleaseEvent(QMouseEvent* _event)
{
	if (mesh_.n_vertices() == 0) return;
	if (_event->button() == LeftButton) {
		// rotate, translate, scale, pick, or deform.
		int mm = mouse_mode();
		if (mm < N_MOUSE_MODES) { // rotate, translate, scale.
			QGLViewerWidget::mouseReleaseEvent(_event);
		} else {                  // pick, deform.
			switch (mm) {
			case MOUSE_PICK:
				processMousePickRelease(_event);
				break;
			case MOUSE_DEFORM:
				processMouseDeformRelease(_event);
				break;
			}
		} // end of else
	}
}

void MeshViewerWidget::processPickHits(GLint hits, GLuint* buffer, bool controled, bool singlePick)
{
	int n_vertices = mesh_.n_vertices();
	if (!n_vertices) return; // no mesh.

	if (!controled) {
		pdmesh_->getSVSet().clear();
	}

	GLdouble winx, winy, winz;
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	GLdouble minZ = 10e8;
	int minZindex = n_vertices;

	GLuint* ptr = buffer;
	for (int i = 0; i < hits; ++i) {
		int names = *ptr;
		assert(names == 1);
		ptr += 3;
		int index = *ptr;
		assert(index <= n_vertices);
		if (index < n_vertices) {
			if (singlePick) {
				InterMesh::Point& p = mesh_.point(InterMesh::VertexHandle(index));
				gluProject(p[0], p[1], p[2], modelview_matrix(), projection_matrix(), 
					viewport, &winx, &winy, &winz);
				if (winz < minZ) {
					minZ = winz;
					minZindex = index;
				}	
			} else {
				pdmesh_->getSVSet().insert(InterMesh::VertexHandle(index));
			}
		}
		ptr++;
	}
	if (singlePick && minZindex != n_vertices) {
		pdmesh_->getSVSet().insert(InterMesh::VertexHandle(minZindex));
	}
	
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
	int i = 0;
	for (CITER cit = selectedVertices.begin(); cit != selectedVertices.end(); ++cit, ++i) {
		Vec3d center = mesh_.point(*cit);
		if (mouse_mode() == MOUSE_DEFORM) {
			glLoadName(i);
		}
		drawSelectBox(center, box_radius_);
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


void MeshViewerWidget::processMousePickPress(QMouseEvent* _event)
{
	// record the press pos.
	pick_press_point_.rx() = _event->x();
	pick_press_point_.ry() = _event->y();
}

void MeshViewerWidget::processMousePickMove(QMouseEvent* _event)
{
	// do nothing.
}

void MeshViewerWidget::processMousePickRelease(QMouseEvent* _event)
{
	std::vector<GLuint> selectBuf(mesh_.n_vertices()); // assert enough space.
	GLint hits;
	GLint viewport[4];

	GLdouble x1, y1, x2, y2, detX, detY, cenX, cenY;

	bool singlePick = false;

	glGetIntegerv(GL_VIEWPORT, viewport);
	glSelectBuffer(mesh_.n_vertices(), &selectBuf[0]);
	glRenderMode(GL_SELECT);
	
	glInitNames();
	glPushName(mesh_.n_vertices());

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	if (_event->x() == pick_press_point_.x() && _event->y() == pick_press_point_.y()) { // point select.
		singlePick = true;
		gluPickMatrix((GLdouble)(_event->x()), (GLdouble)(viewport[3]-_event->y()), 5.0, 5.0, viewport);
	} else { // rect select.
		// compute the pick rect.
		x1 = (GLdouble)(_event->x());
		y1 = (GLdouble)(viewport[3]-(_event->y()));
		x2 = (GLdouble)(pick_press_point_.x());
		y2 = (GLdouble)(viewport[3]-(pick_press_point_.y()));
		detX = fabs(x1-x2);
		if (detX == 0.0) { detX = 1.0; } // we should not call gluPickMatrix() with a zero-like detX, also detY.
		detY = fabs(y1-y2);
		if (detY == 0.0) { detY = 1.0; }
		cenX = (x1+x2)/2;
		cenY = (y1+y2)/2;
		gluPickMatrix(cenX, cenY, detX, detY, viewport);
	}
	glMultMatrixd(projection_matrix());
	draw_scene(draw_mode());

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glFlush();

	hits = glRenderMode(GL_RENDER);
	
	processPickHits(hits, &selectBuf[0], _event->modifiers() == ControlModifier, singlePick);
	updateGL();
}

void MeshViewerWidget::processMouseDeformPress(QMouseEvent* _event)
{
	std::vector<GLuint> selectBuf(mesh_.n_vertices()); // assert enough space.
	GLint hits;
	GLint viewport[4];

	glGetIntegerv(GL_VIEWPORT, viewport);
	glSelectBuffer(mesh_.n_vertices(), &selectBuf[0]);
	glRenderMode(GL_SELECT);
	
	glInitNames();
	glPushName(mesh_.n_vertices());

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	gluPickMatrix((GLdouble)(_event->x()), (GLdouble)(viewport[3]-_event->y()), 5.0, 5.0, viewport);
	glMultMatrixd(projection_matrix());
	draw_scene(draw_mode());

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glFlush();

	hits = glRenderMode(GL_RENDER);
	processDeformHits(hits, &selectBuf[0]);
	//updateGL();

}

void MeshViewerWidget::processMouseDeformMove(QMouseEvent* _event)
{
	if (selectedHandles.empty()) { return; } // no handle is selected.

	// get one handle's pos
	int handleIndex = selectedHandles[0];
	std::vector<InterMesh::VertexHandle>& handleIDs = pdmesh_->getHandleIDs();
	assert(handleIndex < handleIDs.size());
	InterMesh::VertexHandle vh = handleIDs[handleIndex];
	InterMesh::Point p = mesh_.point(vh);

	// get z-buffer of this handle.
	GLdouble winx, winy, winz;
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	gluProject(p[0], p[1], p[2], modelview_matrix(), projection_matrix(), 
			viewport, &winx, &winy, &winz);

	// get the translation vector.
	GLdouble objx, objy, objz;
	gluUnProject((GLdouble)(_event->x()), (GLdouble)(viewport[3]-_event->y()), 
			winz, modelview_matrix(), projection_matrix(), viewport, 
			&objx, &objy, &objz);
	typedef DeformableMesh3d::V3d V3d;
	V3d v(objx - p[0], objy - p[1], objz - p[2]);

	pdmesh_->translate(selectedHandles, v);
	updateGL();
}

void MeshViewerWidget::processMouseDeformRelease(QMouseEvent* _event)
{
	selectedHandles.clear();
}

void MeshViewerWidget::processDeformHits(GLint hits, GLuint* buffer)
{
	int n_vertices = mesh_.n_vertices();
	if (!n_vertices) return; // no mesh.

	GLuint* ptr = buffer;
	selectedHandles.clear();
	for (int i = 0; i < hits; ++i) {
		int names = *ptr;
		assert(names == 1);
		ptr += 3;
		int index = *ptr;
		assert(index <= n_vertices);
		if (index < n_vertices) {
			selectedHandles.push_back(index);
		}
		ptr++;
	}
	sort(selectedHandles.begin(), selectedHandles.end());
}
