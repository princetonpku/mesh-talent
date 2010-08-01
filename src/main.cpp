#include <QApplication>
//#include <GL/glut.h>

#include "MainWindow.h"

int main(int argc, char** argv)
{
	QApplication app(argc, argv);
	//glutInit(&argc, argv);
	MainWindow mainWin;
	mainWin.resize(640, 480);
	mainWin.showMaximized();

	return app.exec();
}
