// xGEM2MT is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
//
// Copyright (C) 2018-2026 Allan Leal, Dmtrii Kulik, G.D. Miron, S.Dmytriieva
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

#include "pyChemicalEngineArray.hpp"
#include "pyMassTransport.hpp"

PYBIND11_MODULE(PyxGEM2MT, m)
{
    m.doc() = "xGEMS2MT library for mass transport modeling with Gibbs energy minimization";
    xGEM2MT::update_loggers(false, "xGEM2MT.log", 3);
    exportChemicalEngineArray(m);
    exportMassTransport(m);
}

