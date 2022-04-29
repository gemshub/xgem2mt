TEMPLATE = app
LANGUAGE = C++
TARGET = demo

CONFIG -= qt
CONFIG -= warn_on
CONFIG += debug
CONFIG += console
CONFIG += c++17

XGEM2MT_CPP = ..
XGEM2MT_H   = $$XGEM2MT_CPP
DEPENDPATH += $$XGEM2MT_H
INCLUDEPATH += $$XGEM2MT_H

DEFINES += IPMGEMPLUGIN
DEFINES += NODEARRAYLEVEL
#DEFINES += NOPARTICLEARRAY
#DEFINES += OVERFLOW_EXCEPT  #compile with nan inf exceptions

#GEMS3K_CPP = ../../standalone/GEMS3K
#GEMS3K_H   = $$GEMS3K_CPP
#DEPENDPATH += $$GEMS3K_H
#INCLUDEPATH += $$GEMS3K_H
#GEM2MT_CPP = ../../standalone/nodearray-gem
#GEM2MT_H   = $$GEM2MT_CPP
#DEPENDPATH += $$GEM2MT_H
#INCLUDEPATH += $$GEM2MT_H

#include($$GEMS3K_CPP/gems3k.pri)
#include($$GEM2MT_CPP/gem2mt.pri)

OBJECTS_DIR = obj

LIBS += -lGEMS3K-static -lGEM2MT-static

HEADERS	 += $$XGEM2MT_H/xGEM2MT/ChemicalEngineArray.hpp \
            $$XGEM2MT_H/xGEM2MT/MassTransport.hpp \
            $$XGEM2MT_H/xGEM2MT/Eigen.hpp \
            $$XGEM2MT_H/xGEM2MT/Index.hpp \
            $$XGEM2MT_H/xGEM2MT/MassTransport_p.h

SOURCES	+=  $$XGEM2MT_CPP/xGEM2MT/ChemicalEngineArray.cpp \
            $$XGEM2MT_CPP/xGEM2MT/MassTransport.cpp \
            demo_gem2mt.cpp
            #demo_gem2mt_string.cpp
            #demo_nodearray_files.cpp
            #demo_nodearray_strings.cpp
