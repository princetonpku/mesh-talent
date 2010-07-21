#ifndef MESHTALENT_MAINWINDOW_H
#define MESHTALENT_MAINWINDOW_H

#include <QMainWindow>

class QAction;
class QGLViewerWidget;

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
private:
	// File Actions.
	QAction* openAction;
	QAction* saveAction;
	QAction* saveAsAction;
	QAction* exitAction;
	// Edit Actions.

	// View Actions.
	QAction* wireFrameAction;
	QAction* solidFlatAction;
	QAction* solidSmoothAction;
	QAction* pointSetAction;

	// Help Actions.

	// Menus.
	QMenu* fileMenu;
	QMenu* editMenu;
	QMenu* viewMenu;
	QMenu* helpMenu;
	// ToolBars.
	QToolBar* fileToolBar;
	QToolBar* viewToolBar;

private:
	QGLViewerWidget* viewer;
};
#endif // MESHTALENT_MAINWINDOW_H
