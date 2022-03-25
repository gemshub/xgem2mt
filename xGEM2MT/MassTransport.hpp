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

// C++ includes
#include <iostream>
#include <memory>
#include <string>

// xGEMSMT includes
#include <xGEM2MT/ChemicalEngineArray.hpp>

namespace xGEM2MT {

/// A wrapper class for TGEM2MT code
class MassTransport
{

public:

    /// Construct a MassTransport object
    MassTransport();

    /// Construct a copy of a MassTransport object.
    MassTransport(const MassTransport& other) = delete;

    /// Assign another MassTransport object to this.
    auto operator=(MassTransport other) -> MassTransport& = delete;

    /// Destroy this ChemicalEngineArray instance.
    virtual ~MassTransport();

    /// Initialize the MassTransport object from a GEM-Selektor generated file.
    /// @param filename The name of the file containing the definition of the gem2mt task
    auto initializeGEM2MT(std::string gem2mt_filename) -> void;

    /// Initialize internal nodes from a GEM-Selektor project files.
    /// @param filename The name of the file containing the definition of the chemical system
    auto initializeNodes(std::string ipm_filename, std::string dbr_filename) -> void;

    /// Call the mass-transport finite-difference coupled routine
    auto runGEM2MT() -> void;

    /// Dump the MassTransport object to file.
    /// @param filepath - IPM work structure file path&name
    auto dump() const -> std::string;

    /// Dump the MassTransport object to file.
    /// @param filepath - IPM work structure file path&name
    auto dumpGEM2MT(const std::string& gem2mt_filename, const std::string& ipm_filename =  {} ) const -> void;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

};

/// Output the state of the MassTransport object.
auto operator<<(std::ostream& out, const MassTransport& engine) -> std::ostream&;

} // namespace xGEM2MT
