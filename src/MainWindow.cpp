#include <QtGui>

#include "MainWindow.h"
#include "MeshViewerWidget.h"

MainWindow::MainWindow()
{
	viewer = new MeshViewerWidget(this);
	setCentralWidget(viewer);
	createActions();
	createMenus();
	createToolBars();
}

void MainWindow::createActions()
{
	openAction = new QAction(tr("&Open"), this);
    openAction->setIcon(QIcon("../images/open.png"));
    openAction->setShortcut(QKeySequence::Open);
    openAction->setStatusTip(tr("Open a mesh file"));
    connect(openAction, SIGNAL(triggered()), viewer, SLOT(open_mesh_query()));

    saveAction = new QAction(tr("&Save"), this);
    saveAction->setIcon(QIcon("../images/save.png"));
    saveAction->setShortcut(QKeySequence::Save);
    saveAction->setStatusTip(tr("Save the mesh to file"));
    connect(saveAction, SIGNAL(triggered()), this, SLOT(save()));

    saveAsAction = new QAction(tr("Save &As..."), this);
    saveAsAction->setStatusTip(tr("Save the mesh under a new name"));
    connect(saveAsAction, SIGNAL(triggered()), this, SLOT(saveAs()));

    exitAction = new QAction(tr("E&xit"), this);
    exitAction->setShortcut(tr("Ctrl+Q"));
    exitAction->setStatusTip(tr("Exit the application"));
    connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));

	wireFrameAction = new QAction(tr("&WireFrame"), this);
	wireFrameAction->setIcon(QIcon("../images/wireframe.png"));
	wireFrameAction->setStatusTip(tr("Using wireFrame showing method"));
	connect(wireFrameAction, SIGNAL(triggered()), this, SLOT(wireFrameShow()));

	solidFlatAction = new QAction(tr("Solid&Flat"), this);
	solidFlatAction->setIcon(QIcon("../images/solidflat.png"));
	solidFlatAction->setStatusTip(tr("Using solidflat showing method"));
	connect(solidFlatAction, SIGNAL(triggered()), this, SLOT(solidFlatShow()));

	solidSmoothAction = new QAction(tr("Solid&Smooth"), this);
	solidSmoothAction->setIcon(QIcon("../images/solidsmooth.png"));
	solidSmoothAction->setStatusTip(tr("Using solidsmooth showing method"));
	connect(solidSmoothAction, SIGNAL(triggered()), this, SLOT(solidSmoothShow()));

	pointSetAction = new QAction(tr("&PointSet"), this);
	pointSetAction->setIcon(QIcon("../images/pointset.png"));
	pointSetAction->setStatusTip(tr("Using pointset showing method"));
	connect(pointSetAction, SIGNAL(triggered()), this, SLOT(pointSetShow()));

	showGraphAction = new QAction(tr("Show &Graph"), this);
	showGraphAction->setIcon(QIcon("../images/showgraph.png"));
	showGraphAction->setStatusTip(tr("Show the deformation graph"));
	connect(showGraphAction, SIGNAL(triggered()), this, SLOT(graphShow()));

	viewAllAction = new QAction(tr("View &All"), this);
	viewAllAction->setIcon(QIcon("../images/viewall.png"));
	viewAllAction->setStatusTip(tr("View the whole object"));
	connect(viewAllAction, SIGNAL(triggered()), this, SLOT(viewAll()));

	mouseRotateAction = new QAction(tr("Mouse &Rotate"), this);
	mouseRotateAction->setIcon(QIcon("../images/mouserotate.png"));
	mouseRotateAction->setStatusTip(tr("Use mouse to rotate the object"));
	mouseRotateAction->setCheckable(true);
	mouseRotateAction->setChecked(true);
	connect(mouseRotateAction, SIGNAL(triggered()), this, SLOT(mouseRotate()));

	mouseTranslateAction = new QAction(tr("Mouse &Translate"), this);
	mouseTranslateAction->setIcon(QIcon("../images/mousetranslate.png"));
	mouseTranslateAction->setStatusTip(tr("Use mouse to translate the object"));
	mouseTranslateAction->setCheckable(true);
	connect(mouseTranslateAction, SIGNAL(triggered()), this, SLOT(mouseTranslate()));

	mouseScaleAction = new QAction(tr("Mouse &Scale"), this);
	mouseScaleAction->setIcon(QIcon("../images/mousescale.png"));
	mouseScaleAction->setStatusTip(tr("Use mouse to scale the object"));
	mouseScaleAction->setCheckable(true);
	connect(mouseScaleAction, SIGNAL(triggered()), this, SLOT(mouseScale()));

	mousePickAction = new QAction(tr("Mouse &Pick"), this);
	mousePickAction->setIcon(QIcon("../images/mousepick.png"));
	mousePickAction->setStatusTip(tr("Use mouse to pick vertices"));
	mousePickAction->setCheckable(true);
	connect(mousePickAction, SIGNAL(triggered()), this, SLOT(mousePick()));

	mouseDeformAction = new QAction(tr("Mouse &Deform"), this);
	mouseDeformAction->setIcon(QIcon("../images/mousedeform.png"));
	mouseDeformAction->setStatusTip(tr("Use mouse to deform the object"));
	mouseDeformAction->setCheckable(true);
	connect(mouseDeformAction, SIGNAL(triggered()), this, SLOT(mouseDeform()));

	genGraphAction = new QAction(tr("Generate Deformation&Graph"), this);
	genGraphAction->setIcon(QIcon("../image/deformationgraph.png"));
	genGraphAction->setStatusTip(tr("Generate the deformation graph"));
	connect(genGraphAction, SIGNAL(triggered()), viewer, SLOT(gen_graph_query()));
}

void MainWindow::createMenus()
{
	fileMenu = menuBar()->addMenu(tr("&File"));
	fileMenu->addAction(openAction);
	fileMenu->addAction(saveAction);
	fileMenu->addAction(saveAsAction);
	fileMenu->addSeparator();
	fileMenu->addAction(exitAction);

	editMenu = menuBar()->addMenu(tr("&Edit"));
	editMenu->addAction(genGraphAction);

	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->addAction(wireFrameAction);
	viewMenu->addAction(solidFlatAction);
	viewMenu->addAction(solidSmoothAction);
	viewMenu->addAction(pointSetAction);
	viewMenu->addAction(showGraphAction);
	viewMenu->addSeparator();
	viewMenu->addAction(viewAllAction);

	helpMenu = menuBar()->addMenu(tr("&Help"));
}

void MainWindow::createToolBars()
{
	fileToolBar = addToolBar(tr("&File"));
    fileToolBar->addAction(openAction);
    fileToolBar->addAction(saveAction);

	viewToolBar = addToolBar(tr("&View"));
	viewToolBar->addAction(wireFrameAction);
	viewToolBar->addAction(solidFlatAction);
	viewToolBar->addAction(solidSmoothAction);
	viewToolBar->addAction(pointSetAction);
	viewToolBar->addAction(showGraphAction);
	viewToolBar->addAction(viewAllAction);

	mouseToolBar = addToolBar(tr("&Mouse"));
	mouseToolBar->addAction(mouseRotateAction);
	mouseToolBar->addAction(mouseTranslateAction);
	mouseToolBar->addAction(mouseScaleAction);
	mouseToolBar->addAction(mousePickAction);
	mouseToolBar->addAction(mouseDeformAction);
}

void MainWindow::open()
{

}

bool MainWindow::save()
{
	return true;
}

bool MainWindow::saveAs()
{
	return true;
}

void MainWindow::wireFrameShow()
{
	viewer->setDrawMode(QGLViewerWidget::WIRE_FRAME);
}

void MainWindow::solidFlatShow()
{
	viewer->setDrawMode(QGLViewerWidget::SOLID_FLAT);
}

void MainWindow::solidSmoothShow()
{
	viewer->setDrawMode(QGLViewerWidget::SOLID_SMOOTH);
}

void MainWindow::pointSetShow()
{
	viewer->setDrawMode(QGLViewerWidget::POINT_SET);
}

void MainWindow::graphShow()
{
	viewer->setDrawMode(MeshViewerWidget::DRAW_GRAPH);
}

void MainWindow::viewAll()
{
	viewer->view_all();
	viewer->updateGL();
}

void MainWindow::mouseRotate()
{
	setAllMouseActionchecked(false);
	mouseRotateAction->setChecked(true);
	viewer->setMouseMode(QGLViewerWidget::MOUSE_ROTATE);
}

void MainWindow::mouseTranslate()
{
	setAllMouseActionchecked(false);
	mouseTranslateAction->setChecked(true);
	viewer->setMouseMode(QGLViewerWidget::MOUSE_TRANSLATE);
}

void MainWindow::mouseScale()
{
	setAllMouseActionchecked(false);
	mouseScaleAction->setChecked(true);
	viewer->setMouseMode(QGLViewerWidget::MOUSE_SCALE);
}

void MainWindow::mousePick()
{
	setAllMouseActionchecked(false);
	mousePickAction->setChecked(true);
	viewer->setMouseMode(MeshViewerWidget::MOUSE_PICK);
}

void MainWindow::mouseDeform()
{
	setAllMouseActionchecked(false);
	mouseDeformAction->setChecked(true);
	viewer->setMouseMode(MeshViewerWidget::MOUSE_DEFORM);
}

void MainWindow::setAllMouseActionchecked(bool b)
{
	mouseRotateAction->setChecked(b);
	mouseTranslateAction->setChecked(b);
	mouseScaleAction->setChecked(b);
	mousePickAction->setChecked(b);
	mouseDeformAction->setChecked(b);
}
