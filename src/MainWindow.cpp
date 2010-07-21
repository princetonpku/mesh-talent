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

	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->addAction(wireFrameAction);
	viewMenu->addAction(solidFlatAction);
	viewMenu->addAction(solidSmoothAction);
	viewMenu->addAction(pointSetAction);

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
