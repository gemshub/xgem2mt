// xGEMS2MT is a C++ and Python library for mass transport modeling with Gibbs energy minimization
//
// Copyright (C) 2022-2026 S.Dmytriieva, Dmitrii Kulik
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

#pragma once

#include "MassTransport.hpp"
#include "m_gem2mt.h"

namespace xGEM2MT {

struct MassTransport::Impl
{
    /// Construct a default Impl instance
    Impl(char mode, size_t nnodes)
    {
        init_gem3k = false;
        gem2mt_task = std::unique_ptr<TGEM2MT>(new TGEM2MT(mode, nnodes));
    }

    bool init_gem3k = false;

    /// The TNodeArray instance from ChemicalEngineArray
    std::unique_ptr<TGEM2MT> gem2mt_task;
    //TGEM2MT gem2mt_task;

    /// The instance of ChemicalEngineArray linked to gem2mt_task internal data
    std::shared_ptr<ChemicalEngineArray> chem_array;

    /// Write IPM, DCH and DBR files (in binary, txt or json mode)
    GEMS3KGenerator::IOModes io_mode = GEMS3KGenerator::f_json;

    // To be done?
    // (1) Used in `Trans1D` for different transport function, the array is allocated but not used in mas transport calculations
    // Advective mass transport over DC gradients in Aq phase or Advective mass transport over IC gradients in Aq phase
    // PvDDc: Use diffusion coefficients for DC - DDc vector (+ -) (default -)
    //void usePvDDc(bool enable);
    // PvDIc: Use diffusion coefficients for IC - DIc vector (+ -) (default -)
    //void usePvDIc(bool enable);

    // (2) PvGrid: Use array of grid point locations (+ -) (default -)
    //void useArrayGridPoints(bool enable);
    // Set grid point location, size is nC [grid]
    // @param pndx: index in array
    // @param x: Array of initial mean particle type numbers per node
    // @param y: Minimum average total number of particles of each type per one node
    // @param z: Maximum average total number of particles of each type per one node
    //void setGrid(long int pndx, double x,  double y,  double z);

};

} // namespace xGEM2MT
