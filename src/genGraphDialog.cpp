#include <QtGui>

#include "genGraphDialog.h"

GenGraphDialog::GenGraphDialog(QWidget* parent)
 : QDialog(parent)
{
	// about graph.
	nodeNumLabel = new QLabel(tr("&Node number:"));
	nodeNumLineEdit = new QLineEdit("300"); // init value should be read from config file.
	nodeNumLabel->setBuddy(nodeNumLineEdit);

	relateNumLabel = new QLabel(tr("&Relate number"));
	relateNumLineEdit = new QLineEdit("3");
	relateNumLabel->setBuddy(relateNumLineEdit);

	sampleScaleLabel = new QLabel(tr("&Sample scale"));
	sampleScaleLineEdit = new QLineEdit("20");
	sampleScaleLabel->setBuddy(sampleScaleLineEdit);

	deleteRadiusLabel = new QLabel(tr("&Delete radius"));
	deleteRadiusLineEdit = new QLineEdit("150");
	deleteRadiusLabel->setBuddy(deleteRadiusLineEdit);
	
	// about numeric.
	ErotLabel = new QLabel(tr("R&otation Energy"));
	ErotLineEdit = new QLineEdit("1");
	ErotLabel->setBuddy(ErotLineEdit);

	EregLabel = new QLabel(tr("R&egularization Energy"));
	EregLineEdit = new QLineEdit("10");
	EregLabel->setBuddy(EregLineEdit);

	EconLabel = new QLabel(tr("&Constraint Energy"));
	EconLineEdit = new QLineEdit("100");
	EconLabel->setBuddy(EconLineEdit);

	// layout
	// graph
	QHBoxLayout* nodeNumLayout = new QHBoxLayout;
	nodeNumLayout->addWidget(nodeNumLabel);
	nodeNumLayout->addWidget(nodeNumLineEdit);

	QHBoxLayout* relateNumLayout = new QHBoxLayout;
	relateNumLayout->addWidget(relateNumLabel);
	relateNumLayout->addWidget(relateNumLineEdit);

	QHBoxLayout* sampleScaleLayout = new QHBoxLayout;
	sampleScaleLayout->addWidget(sampleScaleLabel);
	sampleScaleLayout->addWidget(sampleScaleLineEdit);

	QHBoxLayout* deleteRadiusLayout = new QHBoxLayout;
	deleteRadiusLayout->addWidget(deleteRadiusLabel);
	deleteRadiusLayout->addWidget(deleteRadiusLineEdit);

	QVBoxLayout* graphLayout = new QVBoxLayout;
	graphLayout->addLayout(nodeNumLayout);
	graphLayout->addLayout(relateNumLayout);
	graphLayout->addLayout(sampleScaleLayout);
	graphLayout->addLayout(deleteRadiusLayout);

	QGroupBox* graphGroupBox = new QGroupBox(tr("Deformation Graph"));
	graphGroupBox->setLayout(graphLayout);

	// energy
	QHBoxLayout* ErotLayout = new QHBoxLayout;
	ErotLayout->addWidget(ErotLabel);
	ErotLayout->addWidget(ErotLineEdit);

	QHBoxLayout* EregLayout = new QHBoxLayout;
	EregLayout->addWidget(EregLabel);
	EregLayout->addWidget(EregLineEdit);

	QHBoxLayout* EconLayout = new QHBoxLayout;
	EconLayout->addWidget(EconLabel);
	EconLayout->addWidget(EconLineEdit);

	QVBoxLayout* energyLayout = new QVBoxLayout;
	energyLayout->addLayout(ErotLayout);
	energyLayout->addLayout(EregLayout);
	energyLayout->addLayout(EconLayout);

	QGroupBox* energyGroupBox = new QGroupBox(tr("Energy"));
	energyGroupBox->setLayout(energyLayout);

	QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | 
													   QDialogButtonBox::Cancel);
	connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

	QHBoxLayout* mainLayout = new QHBoxLayout;
	mainLayout->addWidget(graphGroupBox);
	mainLayout->addWidget(energyGroupBox);
	mainLayout->addWidget(buttonBox);
	setLayout(mainLayout);
}
