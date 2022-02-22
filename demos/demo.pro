TEMPLATE = app
LANGUAGE = C++
TARGET = demo

CONFIG -= qt
CONFIG -= warn_on
CONFIG += debug
CONFIG += console
CONFIG += c++17

GEM2MT_CPP = ..
GEM2MT_H   = $$GEM2MT_CPP
DEPENDPATH += $$GEM2MT_H
INCLUDEPATH += $$GEM2MT_H

OBJECTS_DIR = obj

LIBS += -lGEMS3K-static

HEADERS	 += $$GEM2MT_H/xGEM2MT/ChemicalEngineArray.hpp \
            $$GEM2MT_H/xGEM2MT/Eigen.hpp \
            $$GEM2MT_H/xGEM2MT/Index.hpp


SOURCES	+=  $$GEM2MT_CPP/xGEM2MT/ChemicalEngineArray.cpp \
            demo_nodearray.cpp
