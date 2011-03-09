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
	createStatusBar();
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
    connect(saveAction, SIGNAL(triggered()), viewer, SLOT(save_mesh_query()));

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
	wireFrameAction->setCheckable(true);
	wireFrameAction->setChecked(true);
	connect(wireFrameAction, SIGNAL(triggered()), this, SLOT(wireFrameShow()));

	solidFlatAction = new QAction(tr("Solid&Flat"), this);
	solidFlatAction->setIcon(QIcon("../images/solidflat.png"));
	solidFlatAction->setStatusTip(tr("Using solidflat showing method"));
	solidFlatAction->setCheckable(true);
	solidFlatAction->setChecked(false);
	connect(solidFlatAction, SIGNAL(triggered()), this, SLOT(solidFlatShow()));

	solidSmoothAction = new QAction(tr("Solid&Smooth"), this);
	solidSmoothAction->setIcon(QIcon("../images/solidsmooth.png"));
	solidSmoothAction->setStatusTip(tr("Using solidsmooth showing method"));
	solidSmoothAction->setCheckable(true);
	solidSmoothAction->setChecked(false);
	connect(solidSmoothAction, SIGNAL(triggered()), this, SLOT(solidSmoothShow()));

	pointSetAction = new QAction(tr("&PointSet"), this);
	pointSetAction->setIcon(QIcon("../images/pointset.png"));
	pointSetAction->setStatusTip(tr("Using pointset showing method"));
	pointSetAction->setCheckable(true);
	pointSetAction->setChecked(false);
	connect(pointSetAction, SIGNAL(triggered()), this, SLOT(pointSetShow()));

	voronoiDiagramAction = new QAction(tr("&VoronoiDiagram"), this);
	voronoiDiagramAction->setIcon(QIcon("../images/voronoidiagram.png"));
	voronoiDiagramAction->setStatusTip(tr("Using voronoidiagram showing method"));
	voronoiDiagramAction->setCheckable(true);
	voronoiDiagramAction->setChecked(false);
	connect(voronoiDiagramAction, SIGNAL(triggered()), this, SLOT(voronoiDiagramShow()));

	showGraphAction = new QAction(tr("Show &Graph"), this);
	showGraphAction->setIcon(QIcon("../images/showgraph.png"));
	showGraphAction->setStatusTip(tr("Show the deformation graph"));
	showGraphAction->setCheckable(true);
	showGraphAction->setChecked(false);
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

	updateMeshAction = new QAction(tr("&Update Mesh"), this);
	updateMeshAction->setIcon(QIcon("../image/updatemesh.png"));
	updateMeshAction->setStatusTip(tr("Update the normals and the center of the mesh"));
	connect(updateMeshAction, SIGNAL(triggered()), viewer, SLOT(update_mesh()));

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
	editMenu->addSeparator();
	editMenu->addAction(updateMeshAction);

	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->addAction(wireFrameAction);
	viewMenu->addAction(solidFlatAction);
	viewMenu->addAction(solidSmoothAction);
	viewMenu->addAction(pointSetAction);
	viewMenu->addAction(voronoiDiagramAction);
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
	viewToolBar->addAction(voronoiDiagramAction);
	viewToolBar->addAction(showGraphAction);
	viewToolBar->addAction(viewAllAction);

	mouseToolBar = addToolBar(tr("&Mouse"));
	mouseToolBar->addAction(mouseRotateAction);
	mouseToolBar->addAction(mouseTranslateAction);
	mouseToolBar->addAction(mouseScaleAction);
	mouseToolBar->addAction(mousePickAction);
	mouseToolBar->addAction(mouseDeformAction);
}

void MainWindow::createStatusBar()
{
	QLabel* label = new QLabel(tr("no mesh"));
	label->setAlignment(Qt::AlignHCenter);

	statusBar()->addWidget(label);

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
	setAllViewActionChecked(false);
	wireFrameAction->setChecked(true);
	viewer->setDrawMode(QGLViewerWidget::WIRE_FRAME);
}

void MainWindow::solidFlatShow()
{
	setAllViewActionChecked(false);
	solidFlatAction->setChecked(true);
	viewer->setDrawMode(QGLViewerWidget::SOLID_FLAT);
}

void MainWindow::solidSmoothShow()
{
	setAllViewActionChecked(false);
	solidSmoothAction->setChecked(true);
	viewer->setDrawMode(QGLViewerWidget::SOLID_SMOOTH);
}

void MainWindow::pointSetShow()
{
	setAllViewActionChecked(false);
	pointSetAction->setChecked(true);
	viewer->setDrawMode(QGLViewerWidget::POINT_SET);
}

void MainWindow::voronoiDiagramShow()
{
	setAllViewActionChecked(false);
	voronoiDiagramAction->setChecked(true);
	viewer->setDrawMode(QGLViewerWidget::VORONOI_DIAGRAM);
}

void MainWindow::graphShow()
{
	if (!viewer->graphGened()) {
		int ret = QMessageBox::warning(this, tr("show the graph"), 
				tr("The graph has not been generated\nWould you like to generate it now?"), 
				QMessageBox::Ok | QMessageBox::Cancel);
		if (ret) {
			viewer->gen_graph_query();
		} else {
			return;
		}
	}
	setAllViewActionChecked(false);
	showGraphAction->setChecked(true);
	viewer->setDrawMode(MeshViewerWidget::DRAW_GRAPH);
}

void MainWindow::viewAll()
{
	viewer->view_all();
	viewer->updateGL();
}

void MainWindow::mouseRotate()
{
	setAllMouseActionChecked(false);
	mouseRotateAction->setChecked(true);
	viewer->setMouseMode(QGLViewerWidget::MOUSE_ROTATE);
}

void MainWindow::mouseTranslate()
{
	setAllMouseActionChecked(false);
	mouseTranslateAction->setChecked(true);
	viewer->setMouseMode(QGLViewerWidget::MOUSE_TRANSLATE);
}

void MainWindow::mouseScale()
{
	setAllMouseActionChecked(false);
	mouseScaleAction->setChecked(true);
	viewer->setMouseMode(QGLViewerWidget::MOUSE_SCALE);
}

void MainWindow::mousePick()
{
	setAllMouseActionChecked(false);
	mousePickAction->setChecked(true);
	viewer->setMouseMode(MeshViewerWidget::MOUSE_PICK);
}

void MainWindow::mouseDeform()
{
	setAllMouseActionChecked(false);
	mouseDeformAction->setChecked(true);
	viewer->setMouseMode(MeshViewerWidget::MOUSE_DEFORM);
	// let DeformableMesh3d get handles.
	viewer->getHandles();
}

void MainWindow::setAllMouseActionChecked(bool b)
{
	mouseRotateAction->setChecked(b);
	mouseTranslateAction->setChecked(b);
	mouseScaleAction->setChecked(b);
	mousePickAction->setChecked(b);
	mouseDeformAction->setChecked(b);
}

void MainWindow::setAllViewActionChecked(bool b)
{
	wireFrameAction->setChecked(b);
	solidFlatAction->setChecked(b);
	solidSmoothAction->setChecked(b);
	pointSetAction->setChecked(b);
	voronoiDiagramAction->setChecked(b);
	showGraphAction->setChecked(b);
}
