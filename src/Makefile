#############################################################################
# Makefile for building: mesh-talent
# Generated by qmake (2.01a) (Qt 4.6.2) on: ?? 7? 26 21:46:22 2010
# Project:  src.pro
# Template: app
# Command: /usr/bin/qmake -unix -o Makefile src.pro
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_NO_DEBUG -DQT_OPENGL_LIB -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
CFLAGS        = -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
INCPATH       = -I/usr/share/qt4/mkspecs/linux-g++ -I. -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I/usr/include/qt4/QtOpenGL -I/usr/include/qt4 -I. -Imath -I/usr/X11R6/include -I. -I/usr/include/suitesparse
LINK          = g++
LFLAGS        = -Wl,-O1
LIBS          = $(SUBLIBS)  -L/usr/lib -L/usr/X11R6/lib -lQtOpenGL -lQtGui -lQtCore -lGLU -lGL -lpthread -lOpenMeshCore -lglut -lcholmod
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = DeformableMesh3d.cpp \
		DeformationGraph.cpp \
		GaussNewtonSolver.cpp \
		genGraphDialog.cpp \
		main.cpp \
		MainWindow.cpp \
		MeshViewerWidget.cpp \
		QGLViewerWidget.cpp moc_genGraphDialog.cpp \
		moc_MainWindow.cpp \
		moc_MeshViewerWidget.cpp \
		moc_QGLViewerWidget.cpp
OBJECTS       = DeformableMesh3d.o \
		DeformationGraph.o \
		GaussNewtonSolver.o \
		genGraphDialog.o \
		main.o \
		MainWindow.o \
		MeshViewerWidget.o \
		QGLViewerWidget.o \
		moc_genGraphDialog.o \
		moc_MainWindow.o \
		moc_MeshViewerWidget.o \
		moc_QGLViewerWidget.o
DIST          = /usr/share/qt4/mkspecs/common/g++.conf \
		/usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/opengl.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		src.pro
QMAKE_TARGET  = mesh-talent
DESTDIR       = 
TARGET        = mesh-talent

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: src.pro  /usr/share/qt4/mkspecs/linux-g++/qmake.conf /usr/share/qt4/mkspecs/common/g++.conf \
		/usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/opengl.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		/usr/lib/libQtOpenGL.prl \
		/usr/lib/libQtGui.prl \
		/usr/lib/libQtCore.prl
	$(QMAKE) -unix -o Makefile src.pro
/usr/share/qt4/mkspecs/common/g++.conf:
/usr/share/qt4/mkspecs/common/unix.conf:
/usr/share/qt4/mkspecs/common/linux.conf:
/usr/share/qt4/mkspecs/qconfig.pri:
/usr/share/qt4/mkspecs/features/qt_functions.prf:
/usr/share/qt4/mkspecs/features/qt_config.prf:
/usr/share/qt4/mkspecs/features/exclusive_builds.prf:
/usr/share/qt4/mkspecs/features/default_pre.prf:
/usr/share/qt4/mkspecs/features/release.prf:
/usr/share/qt4/mkspecs/features/default_post.prf:
/usr/share/qt4/mkspecs/features/warn_on.prf:
/usr/share/qt4/mkspecs/features/qt.prf:
/usr/share/qt4/mkspecs/features/unix/opengl.prf:
/usr/share/qt4/mkspecs/features/unix/thread.prf:
/usr/share/qt4/mkspecs/features/moc.prf:
/usr/share/qt4/mkspecs/features/resources.prf:
/usr/share/qt4/mkspecs/features/uic.prf:
/usr/share/qt4/mkspecs/features/yacc.prf:
/usr/share/qt4/mkspecs/features/lex.prf:
/usr/share/qt4/mkspecs/features/include_source_dir.prf:
/usr/lib/libQtOpenGL.prl:
/usr/lib/libQtGui.prl:
/usr/lib/libQtCore.prl:
qmake:  FORCE
	@$(QMAKE) -unix -o Makefile src.pro

dist: 
	@$(CHK_DIR_EXISTS) .tmp/mesh-talent1.0.0 || $(MKDIR) .tmp/mesh-talent1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/mesh-talent1.0.0/ && $(COPY_FILE) --parents debugtrace.h DeformableMesh3d.h DeformationGraph.h GaussNewtonSolver.h genGraphDialog.h MainWindow.h MeshViewerWidget.h QGLViewerWidget.h math/matrix3d.h math/point3d.h math/SparseMatrix.h math/vector3d.h .tmp/mesh-talent1.0.0/ && $(COPY_FILE) --parents DeformableMesh3d.cpp DeformationGraph.cpp GaussNewtonSolver.cpp genGraphDialog.cpp main.cpp MainWindow.cpp MeshViewerWidget.cpp QGLViewerWidget.cpp .tmp/mesh-talent1.0.0/ && (cd `dirname .tmp/mesh-talent1.0.0` && $(TAR) mesh-talent1.0.0.tar mesh-talent1.0.0 && $(COMPRESS) mesh-talent1.0.0.tar) && $(MOVE) `dirname .tmp/mesh-talent1.0.0`/mesh-talent1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/mesh-talent1.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_moc_header_make_all: moc_genGraphDialog.cpp moc_MainWindow.cpp moc_MeshViewerWidget.cpp moc_QGLViewerWidget.cpp
compiler_moc_header_clean:
	-$(DEL_FILE) moc_genGraphDialog.cpp moc_MainWindow.cpp moc_MeshViewerWidget.cpp moc_QGLViewerWidget.cpp
moc_genGraphDialog.cpp: genGraphDialog.h
	/usr/bin/moc-qt4 $(DEFINES) $(INCPATH) genGraphDialog.h -o moc_genGraphDialog.cpp

moc_MainWindow.cpp: MainWindow.h
	/usr/bin/moc-qt4 $(DEFINES) $(INCPATH) MainWindow.h -o moc_MainWindow.cpp

moc_MeshViewerWidget.cpp: QGLViewerWidget.h \
		DeformableMesh3d.h \
		math/vector3d.h \
		math/matrix3d.h \
		math/point3d.h \
		MeshViewerWidget.h
	/usr/bin/moc-qt4 $(DEFINES) $(INCPATH) MeshViewerWidget.h -o moc_MeshViewerWidget.cpp

moc_QGLViewerWidget.cpp: QGLViewerWidget.h
	/usr/bin/moc-qt4 $(DEFINES) $(INCPATH) QGLViewerWidget.h -o moc_QGLViewerWidget.cpp

compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: compiler_moc_header_clean 

####### Compile

DeformableMesh3d.o: DeformableMesh3d.cpp DeformableMesh3d.h \
		math/vector3d.h \
		math/matrix3d.h \
		math/point3d.h \
		DeformationGraph.h \
		GaussNewtonSolver.h \
		math/SparseMatrix.h \
		debugtrace.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o DeformableMesh3d.o DeformableMesh3d.cpp

DeformationGraph.o: DeformationGraph.cpp DeformationGraph.h \
		math/vector3d.h \
		math/matrix3d.h \
		math/point3d.h \
		DeformableMesh3d.h \
		debugtrace.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o DeformationGraph.o DeformationGraph.cpp

GaussNewtonSolver.o: GaussNewtonSolver.cpp GaussNewtonSolver.h \
		math/vector3d.h \
		math/point3d.h \
		math/matrix3d.h \
		math/SparseMatrix.h \
		DeformableMesh3d.h \
		DeformationGraph.h \
		debugtrace.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o GaussNewtonSolver.o GaussNewtonSolver.cpp

genGraphDialog.o: genGraphDialog.cpp genGraphDialog.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o genGraphDialog.o genGraphDialog.cpp

main.o: main.cpp MainWindow.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

MainWindow.o: MainWindow.cpp MainWindow.h \
		MeshViewerWidget.h \
		QGLViewerWidget.h \
		DeformableMesh3d.h \
		math/vector3d.h \
		math/matrix3d.h \
		math/point3d.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o MainWindow.o MainWindow.cpp

MeshViewerWidget.o: MeshViewerWidget.cpp MeshViewerWidget.h \
		QGLViewerWidget.h \
		DeformableMesh3d.h \
		math/vector3d.h \
		math/matrix3d.h \
		math/point3d.h \
		genGraphDialog.h \
		DeformationGraph.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o MeshViewerWidget.o MeshViewerWidget.cpp

QGLViewerWidget.o: QGLViewerWidget.cpp QGLViewerWidget.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o QGLViewerWidget.o QGLViewerWidget.cpp

moc_genGraphDialog.o: moc_genGraphDialog.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_genGraphDialog.o moc_genGraphDialog.cpp

moc_MainWindow.o: moc_MainWindow.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_MainWindow.o moc_MainWindow.cpp

moc_MeshViewerWidget.o: moc_MeshViewerWidget.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_MeshViewerWidget.o moc_MeshViewerWidget.cpp

moc_QGLViewerWidget.o: moc_QGLViewerWidget.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_QGLViewerWidget.o moc_QGLViewerWidget.cpp

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

