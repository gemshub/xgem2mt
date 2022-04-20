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

#include "MassTransport.hpp"

// C++ includes
#include <chrono>
#include <iomanip>
#include <memory>

// GEMS3K includes
#include <GEMS3K/gems3k_impex.h>
#include <GEMS3K/io_simdjson.h>
#include <nodearray-gem/m_gem2mt.h>

namespace xGEM2MT {

struct MassTransport::Impl
{
    /// Construct a default Impl instance
    Impl()
    {}

    /// The TNodeArray instance from ChemicalEngineArray
    std::unique_ptr<TGEM2MT> gem2mt_task;

    /// Write IPM, DCH and DBR files (in binary, txt or json mode)
    GEMS3KGenerator::IOModes io_mode = GEMS3KGenerator::f_json;

};

MassTransport::MassTransport()
    : pimpl(new Impl())
{}

MassTransport::~MassTransport()
{}

void MassTransport::initializeGEM2MT(std::string gem2mt_filename)
{
    pimpl->gem2mt_task = std::unique_ptr<TGEM2MT>(new TGEM2MT(0));

    // Initialize the GEMS `nodes` member
    long int res = pimpl->gem2mt_task->ReadTask(gem2mt_filename.c_str(), "");

    // Check if there was a file read error during task initialization
    if(res == 1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the MassTransport object.\n"
                                 "Make sure the provided file exists relative to the working directory.");
}

void MassTransport::initializeGEM2MTString(std::string gem2mt_json)
{
    pimpl->gem2mt_task = std::unique_ptr<TGEM2MT>(new TGEM2MT(0));

    // Initialize the GEMS `nodes` member
    long int res = pimpl->gem2mt_task->ReadTaskString(gem2mt_json);

    // Check if there was a file read error during task initialization
    if(res == 1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the MassTransport object.");
}

void MassTransport::dumpGEM2MT(const std::string &gem2mt_filename, const std::string &ipm_filename) const
{
    pimpl->gem2mt_task->WriteTask(gem2mt_filename.c_str());
}

void MassTransport::initializeNodes(std::string ipm_filename, std::string dbr_filename)
{
    try{
        // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
        long int res = pimpl->gem2mt_task->MassTransInit( ipm_filename.c_str(), dbr_filename.c_str());

        // Check if there was a file read error during node initialization
        if(res == 1)
            throw std::runtime_error("\n*** ERROR ***\n"
                                     "Could not initialize the MassTransport TNodeArray data.\n"
                                     "Make sure the provided file exists relative to the working directory.");
    }
    catch(TError& err)
    {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the MassTransport TNodeArray data.\n"+err.mess);
    }
}

void MassTransport::initializeNodesFromJsonStrings(std::string dch_json, std::string ipm_json,
                                                  std::vector<std::string> dbr_json)
{
    try{
        // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
        long int res = pimpl->gem2mt_task->MassTransStringInit(dch_json, ipm_json, dbr_json);

        // Check if there was a file read error during node initialization
        if(res == 1)
            throw std::runtime_error("\n*** ERROR ***\n"
                                     "Could not initialize the MassTransport object from JSON strings.\n"
                                     " There was a problem during memory allocation.");
    }
    catch(TError& err)
    {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the MassTransport TNodeArray object from JSON strings.\n"+err.mess);
    }
}

void MassTransport::runGEM2MT()
{
    try{
        pimpl->gem2mt_task->RecCalc();
    }
    catch(TError& err)
    {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Error run the MassTransport .\n"+err.mess);
    }
}

string MassTransport::dump() const
{
    std::stringstream ss;
    io_formats::SimdJsonWrite out_format( ss, "", true );
    pimpl->gem2mt_task->to_text_file( out_format, true, false );
    return ss.str();
}

auto operator<<(std::ostream& out, const MassTransport& state) -> std::ostream&
{
    auto flags = out.flags();

    out << state.dump() << std::endl;

    // Recover the previous state of `out`
    out.flags(flags);
    return out;
}


} // namespace xGEM2MT
