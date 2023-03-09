// xGEMS2MT is a C++ and Python library for mass transport modeling with Gibbs energy minimization
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

#pragma once

#include "MassTransport.hpp"
// GEMS3K includes
#include <nodearray-gem/m_gem2mt.h>

namespace xGEM2MT {

struct MassTransport::Impl
{
    /// Construct a default Impl instance
    Impl()
    {}

    /// The TNodeArray instance from ChemicalEngineArray
    std::unique_ptr<TGEM2MT> gem2mt_task;

    /// The instance of ChemicalEngineArray linked to gem2mt_task internal data
    std::shared_ptr<ChemicalEngineArray> chem_array;

    /// Write IPM, DCH and DBR files (in binary, txt or json mode)
    GEMS3KGenerator::IOModes io_mode = GEMS3KGenerator::f_json;

};

} // namespace xGEM2MT
