#ifndef MESHTALENT_QGLVIEWERWIDGET_HH
#define MESHTALENT_QGLVIEWERWIDGET_HH

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <QGLWidget>

class QGLViewerWidget : public QGLWidget {
	Q_OBJECT
public:
	// Constructors.
	QGLViewerWidget( QWidget* _parent=0 );
	QGLViewerWidget( QGLFormat& _fmt, QWidget* _parent=0 );
	// Destructor.
	virtual ~QGLViewerWidget();
private:
	void init(void);
public:
	/* Sets the center and size of the whole scene. 
	   The _center is used as fixpoint for rotations and for adjusting
	   the camera/viewer (see view_all()). */
	void set_scene_pos(const OpenMesh::Vec3d& _center, float _radius);  

	/* view the whole scene: the eye point is moved far enough from the
	   center so that the whole scene is visible. */
	void view_all();

	float radius() const { return radius_; }
	const OpenMesh::Vec3d& center() const { return center_; }

	const GLdouble* modelview_matrix() const  { return modelview_matrix_;  }
	const GLdouble* projection_matrix() const { return projection_matrix_; }

	float fovy() const { return 45.0f; }
public:
	// draw modes.
	enum { WIRE_FRAME=0, SOLID_FLAT, SOLID_SMOOTH, POINT_SET, VORONOI_DIAGRAM, N_DRAW_MODES };
	void setDrawMode(int dm) { draw_mode_ = dm; updateGL(); }
	int draw_mode() const { return draw_mode_; }
public:
	// mouse modes.
	enum { MOUSE_ROTATE=0, MOUSE_TRANSLATE, MOUSE_SCALE, N_MOUSE_MODES };
	void setMouseMode(int mm) { mouse_mode_ = mm; }
	int mouse_mode() const { return mouse_mode_; }
protected:
	// draw the scene: will be called by the painGL() method.
	virtual void draw_scene(int drawmode);
	
	void setDefaultMaterial(void);
	void setDefaultLight(void);
	
private: // inherited
	// initialize OpenGL states (triggered by Qt)
	void initializeGL();
	// draw the scene (triggered by Qt)
	void paintGL();
	// handle resize events (triggered by Qt)
	void resizeGL(int w, int h);
protected:
	// Qt mouse events
	virtual void mousePressEvent(QMouseEvent*);
	virtual void mouseReleaseEvent(QMouseEvent*);
	virtual void mouseMoveEvent(QMouseEvent*);
	virtual void wheelEvent(QWheelEvent*);
	virtual void keyPressEvent(QKeyEvent*);
private:
	// updates projection matrix
	void update_projection_matrix();
	// translate the scene and update modelview matrix
	void translate(const OpenMesh::Vec3d& _trans);
	// rotate the scene (around its center) and update modelview matrix
	void rotate(const OpenMesh::Vec3d& _axis, float _angle);
private:
	int draw_mode_;
	int mouse_mode_;
	OpenMesh::Vec3d  center_;
	float            radius_;
	GLdouble projection_matrix_[16];
	GLdouble modelview_matrix_[16];

	// virtual trackball: map 2D screen point to unit sphere
	bool map_to_sphere(const QPoint& _point, OpenMesh::Vec3d& _result);
	
	QPoint           last_point_2D_;
	OpenMesh::Vec3d  last_point_3D_;
	bool             last_point_ok_;
};

#endif // MESHTALENT_QGLVIEWERWIDGET_HH
