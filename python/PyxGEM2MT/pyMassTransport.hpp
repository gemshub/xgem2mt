// xGEM2MT is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
//
// Copyright (C) 2022 S.Dmytriieva, Dmitrii Kulik
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;
//PYBIND11_DECLARE_HOLDER_TYPE(T, SmartPtr<T>);

// xGEM2MT includes
#include <xGEM2MT/MassTransport.hpp>
using namespace xGEM2MT;


void exportMassTransport(py::module& m)
{
    py::class_<MassTransport, std::shared_ptr<MassTransport>>(m, "MassTransport")
            .def(py::init())
            .def("initializeGEM2MT", &MassTransport::initializeGEM2MT)
            .def("initializeGEM2MTFromJsonString", &MassTransport::initializeGEM2MTFromJsonString)
            .def("initializeNodes", &MassTransport::initializeNodes)
            .def("initializeNodesFromJsonStrings", &MassTransport::initializeNodesFromJsonStrings)
            .def("runGEM2MT", &MassTransport::runGEM2MT)
            .def("getChemicalEngineArray", &MassTransport::getChemicalEngineArray)
            .def("dump", &MassTransport::dump)
            .def("dumpGEM2MT", &MassTransport::dumpGEM2MT)
            .def("__repr__", [](const MassTransport& self) { std::stringstream ss; ss << self; return ss.str(); })
            ;

}
