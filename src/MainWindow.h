#ifndef MESHTALENT_MAINWINDOW_H
#define MESHTALENT_MAINWINDOW_H

#include <QMainWindow>

class QAction;
class MeshViewerWidget;

class MainWindow : public QMainWindow {
	Q_OBJECT
public:
	MainWindow();
private slots:
	void open();
	bool save();
	bool saveAs();
private:
	void createActions();
	void createMenus();
	void createToolBars();
private slots:
	void wireFrameShow();
	void solidFlatShow();
	void solidSmoothShow();
	void pointSetShow();
	void graphShow();
	void viewAll();

	void mouseRotate();
	void mouseScale();
	void mouseTranslate();
	void mousePick();
	void mouseDeform();
private:
	// File Actions.
	QAction* openAction;
	QAction* saveAction;
	QAction* saveAsAction;
	QAction* exitAction;

	// Edit Actions.
	QAction* genGraphAction;

	// View Actions.
	QAction* wireFrameAction;
	QAction* solidFlatAction;
	QAction* solidSmoothAction;
	QAction* pointSetAction;
	QAction* showGraphAction;
	QAction* viewAllAction;

	// Mouse Actions.
	QAction* mouseRotateAction;
	QAction* mouseTranslateAction;
	QAction* mouseScaleAction;
	QAction* mousePickAction;
	QAction* mouseDeformAction;

	// Help Actions.

	// Menus.
	QMenu* fileMenu;
	QMenu* editMenu;
	QMenu* viewMenu;
	QMenu* helpMenu;
	// ToolBars.
	QToolBar* fileToolBar;
	QToolBar* viewToolBar;
	QToolBar* mouseToolBar;
private:
	void setAllMouseActionchecked(bool b);

private:
	MeshViewerWidget* viewer;
};
#endif // MESHTALENT_MAINWINDOW_H
