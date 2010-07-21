#ifndef MESHTALENT_MESHVIEWERWIDGET_H
#define MESHTALENT_MESHVIEWERWIDGET_H

#include <QString>
#include <QMessageBox>
#include <QFileDialog>

#include <GL/glut.h>

#include "QGLViewerWidget.h"
#include "DeformableMesh3d.h"

using namespace meshtalent;

class MeshViewerWidget : public QGLViewerWidget {
	Q_OBJECT
public:
	typedef DeformableMesh3d::InterMesh InterMesh;
public:
	MeshViewerWidget(QWidget* parent = 0);
	~MeshViewerWidget();
public:
	void open_mesh_gui(QString fname);
	bool openMesh(const char* filename);
	InterMesh& mesh() { return mesh_; };
	const InterMesh& mesh() const { return mesh_; };
public slots:
	void open_mesh_query() {
        QString fileName = QFileDialog::getOpenFileName(this,
            tr("Open mesh file"),
            tr(""),
            tr("OBJ Files (*.obj);;"
            "OFF Files (*.off);;"
            "STL Files (*.stl);;"
            "All Files (*)"));
        if (!fileName.isEmpty())
            open_mesh_gui(fileName);
	}
protected:
	virtual void draw_scene(int drawmode);
private:
	void init_showlist();
	void draw_mesh_wireframe() const;
	void draw_mesh_solidflat() const;
	void draw_mesh_solidsmooth() const;
	void draw_mesh_pointset() const;
private:
	InterMesh mesh_;
	DeformableMesh3d* pdmesh_;
private:
	GLuint showlist_start_;
};

#endif // MESHTALENT_MESHVIEWERWIDGET_H
