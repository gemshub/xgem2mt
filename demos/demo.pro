TEMPLATE = app
LANGUAGE = C++
TARGET = demo

CONFIG -= qt
CONFIG -= warn_on
CONFIG += debug
CONFIG += console
CONFIG += c++20

XGEM2MT_CPP = ..
XGEM2MT_H   = $$XGEM2MT_CPP
DEPENDPATH += $$XGEM2MT_H
INCLUDEPATH += $$XGEM2MT_H

#DEFINES += NO_NODEARRAYLEVEL
#DEFINES += USE_NLOHMANNJSON
#DEFINES += OVERFLOW_EXCEPT  #compile with nan inf exceptions

#GEMS3K_CPP = ../../GEMS3K/GEMS3K
#GEMS3K_H   = $$GEMS3K_CPP
#DEPENDPATH += $$GEMS3K_H
#INCLUDEPATH += $$GEMS3K_H
#include($$GEMS3K_CPP/gems3k.pri)

OBJECTS_DIR = obj

LIBS += -lGEMS3K

HEADERS	 += $$XGEM2MT_H/xGEM2MT/ChemicalEngineArray.hpp \
            $$XGEM2MT_H/xGEM2MT/MassTransport.hpp \
            $$XGEM2MT_H/xGEM2MT/Eigen.hpp \
            $$XGEM2MT_H/xGEM2MT/Index.hpp \
            $$XGEM2MT_H/xGEM2MT/MassTransport_p.h \
            $$XGEM2MT_H/xGEM2MT/m_gem2mt.h \
            $$XGEM2MT_H/xGEM2MT/particlearray.h

SOURCES	+=  $$XGEM2MT_CPP/xGEM2MT/ChemicalEngineArray.cpp \
            $$XGEM2MT_CPP/xGEM2MT/MassTransport.cpp \
            $$XGEM2MT_CPP/xGEM2MT/m_gem2mtbox.cpp \
            $$XGEM2MT_CPP/xGEM2MT/m_gem2mtfor.cpp \
            $$XGEM2MT_CPP/xGEM2MT/m_gem2mtsep.cpp \
            $$XGEM2MT_CPP/xGEM2MT/m_gem2mtt.cpp \
            $$XGEM2MT_CPP/xGEM2MT/particlearray.cpp \
            ../docs/develop/test_main.cpp \
            #demo_nodearray_files.cpp
            #demo_nodearray_strings.cpp
            #demo_gem2mt.cpp
            #demo_gem2mt_string.cpp
