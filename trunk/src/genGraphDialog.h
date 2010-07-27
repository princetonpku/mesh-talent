#ifndef MESHTALENT_GENGRAPHDIALOG_H
#define MESHTALENT_GENGRAPHDIALOG_H

#include <QDialog>

class QLabel;
class QLineEdit;

class GenGraphDialog : public QDialog {
	Q_OBJECT
public:
	GenGraphDialog(QWidget* parent = 0);
public:
	// about graph.
	QLabel* nodeNumLabel;
	QLabel* relateNumLabel;
	QLabel* sampleScaleLabel;
	QLabel* deleteRadiusLabel;

	QLineEdit* nodeNumLineEdit;
	QLineEdit* relateNumLineEdit;
	QLineEdit* sampleScaleLineEdit;
	QLineEdit* deleteRadiusLineEdit;

	// about numeric.
	QLabel* ErotLabel;
	QLabel* EregLabel;
	QLabel* EconLabel;

	QLineEdit* ErotLineEdit;
	QLineEdit* EregLineEdit;
	QLineEdit* EconLineEdit;
};

#endif // MESHTALENT_GENGRAPHDIALOG_H
