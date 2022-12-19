QT += core gui multimedia opengl widgets printsupport

CONFIG += c++17

RC_ICONS = res/aci.ico

unix: QMAKE_CXXFLAGS += -Wno-deprecated-copy -Wno-implicit-fallthrough

unix:!macx: LIBS += -L$$PWD/lib/linux -lGL -lglut -lopenblas -lgfortran -lquadmath -lgomp

win32{
LIBS += -lopengl32
msvc:LIBS += -L$$PWD/lib/win-msvc -llibopenblas
gcc:LIBS += -L$$PWD/lib/win-gcc -lopenblas -lgfortran -lquadmath
}

DEFINES += ARMA_DONT_USE_ATLAS ARMA_USE_SUPERLU
DEFINES += QCUSTOMPLOT_USE_OPENGL

INCLUDEPATH += include \
    include/QCustomPlot \
    src

SOURCES += \
    src/CustomFilter.cpp \
    include/QCustomPlot/qcustomplot.cpp \
    src/CubicSpline.cpp \
    src/IntegrationScheme.cpp \
    src/alcongrator.cpp \
    src/MainWindow.cpp \
    src/resampling.cpp

HEADERS += \
    src/CustomFilter.h \
    include/QCustomPlot/qcustomplot.h \
    src/CubicSpline.h \
    src/IntegrationScheme.h \
    src/MainWindow.h \
    src/resampling.h

FORMS += \
    form/CustomFilter.ui \
    form/MainWindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    res/res.qrc

include(superlu.pro)
