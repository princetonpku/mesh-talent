#include <iomanip>
#include <sstream>
#include <algorithm>

#include <qapplication.h>
#include <QMouseEvent>

#include "QGLViewerWidget.h"

//#include <GL/glut.h>

#if !defined(M_PI)
#  define M_PI 3.1415926535897932
#endif

const double TRACKBALL_RADIUS = 0.6;

using namespace Qt;
using namespace OpenMesh;

QGLViewerWidget::QGLViewerWidget(QWidget* _parent)
	: QGLWidget(_parent)
{    
	init();
}

QGLViewerWidget::
QGLViewerWidget(QGLFormat& _fmt, QWidget* _parent)
	: QGLWidget(_fmt, _parent)
{
	init();
}

void QGLViewerWidget::init(void)
{
	// qt stuff
	setAttribute(Qt::WA_NoSystemBackground, true);
	setFocusPolicy(Qt::StrongFocus);
	//setAcceptDrops( true );  
	//setCursor(PointingHandCursor);

	// draw mode
	draw_mode_ = 0;
}

QGLViewerWidget::~QGLViewerWidget()
{
}

void QGLViewerWidget::setDefaultMaterial(void)
{
	GLfloat mat_a[] = {0.1, 0.1, 0.1, 1.0};
	GLfloat mat_d[] = {0.7, 0.7, 0.5, 1.0};
	GLfloat mat_s[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat shine[] = {120.0};
	
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
}

void QGLViewerWidget::setDefaultLight(void)
{
	GLfloat pos1[] = { 0.1,  0.1, -0.02, 0.0};
	GLfloat pos2[] = {-0.1,  0.1, -0.02, 0.0};
	GLfloat pos3[] = { 0.0,  0.0,  0.1,  0.0};
	GLfloat col1[] = { 0.7,  0.7,  0.8,  1.0};
	GLfloat col2[] = { 0.8,  0.7,  0.7,  1.0};
	GLfloat col3[] = { 1.0,  1.0,  1.0,  1.0};
	
	glEnable(GL_LIGHT0);    
	glLightfv(GL_LIGHT0,GL_POSITION, pos1);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,  col1);
	glLightfv(GL_LIGHT0,GL_SPECULAR, col1);
	
	glEnable(GL_LIGHT1);  
	glLightfv(GL_LIGHT1,GL_POSITION, pos2);
	glLightfv(GL_LIGHT1,GL_DIFFUSE,  col2);
	glLightfv(GL_LIGHT1,GL_SPECULAR, col2);
	
	glEnable(GL_LIGHT2);  
	glLightfv(GL_LIGHT2,GL_POSITION, pos3);
	glLightfv(GL_LIGHT2,GL_DIFFUSE,  col3);
	glLightfv(GL_LIGHT2,GL_SPECULAR, col3);
}

void QGLViewerWidget::initializeGL()
{  
	// OpenGL state
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glDisable( GL_DITHER );
	glEnable( GL_DEPTH_TEST );

	// Material
	setDefaultMaterial();
	// Lighting
	glLoadIdentity();
	setDefaultLight();  
	
	// scene pos and size
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
	set_scene_pos(Vec3d(0.0, 0.0, 0.0), 1.0);
}

void QGLViewerWidget::resizeGL( int _w, int _h )
{
	glViewport(0, 0, _w, _h);
	update_projection_matrix();
	updateGL();
}

void QGLViewerWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( projection_matrix_ );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( modelview_matrix_ );
	
	draw_scene(draw_mode_);
}

void QGLViewerWidget::draw_scene(int drawmode)
{  
	assert(drawmode < N_DRAW_MODES);
	switch (drawmode) {
	case WIRE_FRAME:
		glDisable(GL_LIGHTING);
		//glutWireTeapot(0.5);
		break;
	case SOLID_FLAT:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		//glutSolidTeapot(0.5);
		break;
	case SOLID_SMOOTH:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		//glutSolidTeapot(0.5);
		break;
	default:
		break;
	}
}

void QGLViewerWidget::mousePressEvent(QMouseEvent* _event)
{
	assert(mouse_mode_ < N_MOUSE_MODES);
	last_point_ok_ = map_to_sphere( last_point_2D_=_event->pos(),
				    last_point_3D_ );
}

void QGLViewerWidget::mouseMoveEvent(QMouseEvent* _event)
{  
	//assert(mouse_mode_ < N_MOUSE_MODES);

	QPoint newPoint2D = _event->pos(); 
  
	float  value_y;
	Vec3d  newPoint3D;
	bool   newPoint_hitSphere = map_to_sphere(newPoint2D, newPoint3D);
	
	float dx = newPoint2D.x() - last_point_2D_.x();
	float dy = newPoint2D.y() - last_point_2D_.y();
	
	float w  = width();
	float h  = height();

	// enable GL context
	makeCurrent();
  
	if (last_point_ok_) {
		switch (mouse_mode_) {
		case MOUSE_ROTATE:
			if (newPoint_hitSphere) {
				Vec3d axis = last_point_3D_ % newPoint3D;
				if (axis.sqrnorm() < 1e-7) {
					axis = Vec3d(1,0,0);
				} else {
					axis.normalize();
				}
				// find the amount of rotation
				Vec3d d = last_point_3D_ - newPoint3D;
				float t = 0.5*d.norm()/TRACKBALL_RADIUS;
				if (t<-1.0) t=-1.0;
				else if (t>1.0) t=1.0;
				float phi = 2.0 * asin(t);
				float  angle =  phi * 180.0 / M_PI;
				rotate( axis, angle );
			}
			break;
		case MOUSE_TRANSLATE:
			{
			float z = - (modelview_matrix_[ 2]*center_[0] + 
					modelview_matrix_[ 6]*center_[1] + 
					modelview_matrix_[10]*center_[2] + 
					modelview_matrix_[14]) /
						(modelview_matrix_[ 3]*center_[0] + 
						 modelview_matrix_[ 7]*center_[1] + 
						 modelview_matrix_[11]*center_[2] + 
						 modelview_matrix_[15]);

			float aspect     = w / h;
			float near_plane = 0.01 * radius_;
			float top        = tan(fovy()/2.0f*M_PI/180.0f) * near_plane;
			float right      = aspect*top;

			translate(Vec3d( 2.0*dx/w*right/near_plane*z, 
						-2.0*dy/h*top/near_plane*z, 
						0.0f));
			}
			break;
		case MOUSE_SCALE:
			value_y = radius_ * dy * 3.0 / h; // why use this formula?
			translate(Vec3d(0.0, 0.0, value_y));
			break;
		default:
			break;
		}
	} // end of if

	// remember this point
	last_point_2D_ = newPoint2D;
	last_point_3D_ = newPoint3D;
	last_point_ok_ = newPoint_hitSphere;
	
	// trigger redraw
	updateGL();
}

void QGLViewerWidget::mouseReleaseEvent(QMouseEvent* /* _event */ )
{  
	//assert(mouse_mode_ < N_MOUSE_MODES);
	last_point_ok_ = false;
}

void QGLViewerWidget::wheelEvent(QWheelEvent* _event)
{
	// Use the mouse wheel to zoom in/out
	float d = -(float)_event->delta() / 120.0 * 0.2 * radius_;
	translate(Vec3d(0.0, 0.0, d));
	updateGL();
	_event->accept();
}

void QGLViewerWidget::keyPressEvent( QKeyEvent* _event)
{
	switch(_event->key())
	{
    case Key_I:
		std::cout << "Scene radius: " << radius_ << std::endl;
		std::cout << "Scene center: " << center_ << std::endl;
		break;
    case Key_Q:
    case Key_Escape:
      qApp->quit();      
	}
	_event->ignore();
}

void QGLViewerWidget::translate(const Vec3d& _trans)
{
	// Translate the object by _trans
	// Update modelview_matrix_
	makeCurrent();
	glLoadIdentity();
	glTranslated( _trans[0], _trans[1], _trans[2] );
	glMultMatrixd( modelview_matrix_ );
	glGetDoublev( GL_MODELVIEW_MATRIX, modelview_matrix_);
}

void QGLViewerWidget::rotate(const Vec3d& _axis, float _angle)
{
	// Rotate around center center_, axis _axis, by angle _angle
	// Update modelview_matrix_

	Vec3d t( modelview_matrix_[0]*center_[0] + 
		modelview_matrix_[4]*center_[1] +
		modelview_matrix_[8]*center_[2] + 
		modelview_matrix_[12],
		modelview_matrix_[1]*center_[0] + 
		modelview_matrix_[5]*center_[1] +
		modelview_matrix_[9]*center_[2] + 
		modelview_matrix_[13],
		modelview_matrix_[2]*center_[0] + 
		modelview_matrix_[6]*center_[1] +
		modelview_matrix_[10]*center_[2] + 
		modelview_matrix_[14] );
	
	makeCurrent();
	glLoadIdentity();
	glTranslatef(t[0], t[1], t[2]);
	glRotated( _angle, _axis[0], _axis[1], _axis[2]);
	glTranslatef(-t[0], -t[1], -t[2]); 
	glMultMatrixd(modelview_matrix_);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
}

bool QGLViewerWidget::map_to_sphere(const QPoint& _v2D, Vec3d& _v3D)
{
	// This is actually doing the Sphere/Hyperbolic sheet hybrid thing,
    // based on Ken Shoemake's ArcBall in Graphics Gems IV, 1993.
    double x =  (2.0*_v2D.x() - width())/width();
    double y = -(2.0*_v2D.y() - height())/height();
    double xval = x;
    double yval = y;
    double x2y2 = xval*xval + yval*yval;

    const double rsqr = TRACKBALL_RADIUS*TRACKBALL_RADIUS;
    _v3D[0] = xval;
    _v3D[1] = yval;
    if (x2y2 < 0.5*rsqr) {
        _v3D[2] = sqrt(rsqr - x2y2);
    } else {
        _v3D[2] = 0.5*rsqr/sqrt(x2y2);
    }
    
    return true;
}

void QGLViewerWidget::update_projection_matrix()
{
	makeCurrent();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (GLfloat) width() / (GLfloat) height(),
			0.01*radius_, 100.0*radius_);
	glGetDoublev(GL_PROJECTION_MATRIX, projection_matrix_);
	glMatrixMode(GL_MODELVIEW);
}

void QGLViewerWidget::view_all()
{  
	translate( Vec3d( -(modelview_matrix_[0]*center_[0] + 
				modelview_matrix_[4]*center_[1] +
				modelview_matrix_[8]*center_[2] + 
				modelview_matrix_[12]),
				-(modelview_matrix_[1]*center_[0] + 
				modelview_matrix_[5]*center_[1] +
				modelview_matrix_[9]*center_[2] + 
				modelview_matrix_[13]),
				-(modelview_matrix_[2]*center_[0] + 
				modelview_matrix_[6]*center_[1] +
				modelview_matrix_[10]*center_[2] + 
				modelview_matrix_[14] +
				3.0*radius_) ) );
}

void QGLViewerWidget::set_scene_pos(const Vec3d& _cog, float _radius)
{
	center_ = _cog;
	radius_ = _radius;

	update_projection_matrix();
	view_all();
}
