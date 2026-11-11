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

#include "MassTransport.hpp"
#include "MassTransport_p.h"


// GEMS3K includes
#include <GEMS3K/io_simdjson.h>
#include <GEMS3K/gems3k_impex.h>
#include <chrono>
#include <sstream>
#include <cmath>
#include <filesystem>
#include <fstream>

namespace xGEM2MT {

MassTransport::MassTransport(char mode, size_t nnodes)
    : pimpl(new Impl(mode, nnodes))
{}

MassTransport::~MassTransport()
{}

void MassTransport::readTask(const std::string& gem2mt_filename)
{
    // Initialize the GEMS `nodes` member
    long int res = pimpl->gem2mt_task->ReadTask(gem2mt_filename, "");

    // Check if there was a file read error during task initialization
    if(res == 1) {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the MassTransport object.\n"
                                 "Make sure the provided file exists relative to the working directory.");
    }
}

void MassTransport::readTaskString(const std::string&  gem2mt_json)
{
    // Initialize the GEMS `nodes` member
    long int res = pimpl->gem2mt_task->ReadTaskString(gem2mt_json);

    // Check if there was a error during task initialization
    if(res == 1) {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the MassTransport object.");
    }
}

void MassTransport::dumpTask(const std::string &gem2mt_filename, const std::string &ipm_filename) const
{
    pimpl->gem2mt_task->WriteTask(gem2mt_filename);
    if(!ipm_filename.empty()){
        // to do write ipm files
    }
}

void MassTransport::initNodes(const std::string& ipm_filename, const std::string& dbr_filename)
{
    try{
        // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
        long int res = pimpl->gem2mt_task->MassTransInit( ipm_filename.c_str(), dbr_filename.c_str());

        // Check if there was a file read error during node initialization
        if(res == 1) {
            throw std::runtime_error("\n*** ERROR ***\n"
                                     "Could not initialize the MassTransport TNodeArray data.\n"
                                     "Make sure the provided file exists relative to the working directory.");
        }
    }
    catch(TError& err) {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the MassTransport TNodeArray data.\n"+err.mess);
    }
    pimpl->chem_array = std::shared_ptr<ChemicalEngineArray>(new ChemicalEngineArray(*this));
}

void MassTransport::initNodes(const std::string& ipm_filename)
{
    namespace fs = std::filesystem;
    std::vector<std::string> dbr_names;
    try {
        GEMS3KGenerator generator(ipm_filename);   // reads the file names in the *-dat.lst file
        dbr_names = generator.dbr_names();
    }
    catch(TError& err) {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not read "+ipm_filename+"\n"+err.mess);
    }
    if(dbr_names.empty())
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "No dbr files are listed in "+ipm_filename+
                                 " (list them after the dch and ipm files, separated by spaces).");

    // GEMS3K reads the node array only from a list file, with names relative to the folder of that file.
    // Write the names to a temporary hidden list file in the folder of the *-dat.lst file.
    const fs::path dat_path(ipm_filename);
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const fs::path lst_path = dat_path.parent_path() /
            ("." + dat_path.stem().string() + "-dbr-" + std::to_string(stamp) + ".lst");
    {
        std::ofstream lst(lst_path);
        if(!lst.good())
            throw std::runtime_error("\n*** ERROR ***\n"
                                     "Could not write the temporary file "+lst_path.string()+
                                     "\nUse initNodes(dat_lst, dbr_lst) with a *-dbr.lst file.");
        for(size_t i=0; i<dbr_names.size(); ++i)
            lst << (i ? ", " : "") << '"' << dbr_names[i] << '"';
    }
    struct RemoveFile {
        fs::path path;
        ~RemoveFile() { std::error_code ec; fs::remove(path, ec); }
    } remove_lst{lst_path};

    initNodes(ipm_filename, lst_path.string());
}

void MassTransport::initNodesFromStrings(const std::string& dch_json,
                                         const std::string& ipm_json,
                                         const std::vector<std::string>& dbr_json)
{
    try{
        // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
        long int res = pimpl->gem2mt_task->MassTransStringInit(dch_json, ipm_json, dbr_json);

        // Check if there was a file read error during node initialization
        if(res == 1) {
            throw std::runtime_error("\n*** ERROR ***\n"
                                     "Could not initialize the MassTransport object from JSON strings.\n"
                                     " There was a problem during memory allocation.");
        }
    }
    catch(TError& err) {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the MassTransport TNodeArray object from JSON strings.\n"+err.mess);
    }
    pimpl->chem_array = std::shared_ptr<ChemicalEngineArray>(new ChemicalEngineArray(*this));
}

void MassTransport::runMassTransport()
{
    try {
        pimpl->gem2mt_task->RecCalc();
    }
    catch(TError& err)  {
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Error run the MassTransport .\n"+err.mess);
    }
}

size_t MassTransport::nNodes() const
{
    return pimpl->gem2mt_task->nNodes();
}

std::string MassTransport::name() const
{
    return pimpl->gem2mt_task->name();
}

void MassTransport::setName(const std::string &taskname)
{
    pimpl->gem2mt_task->setName(taskname);
}

std::string MassTransport::comment() const
{
    return pimpl->gem2mt_task->comment();
}

void MassTransport::setComment(const std::string &comment)
{
    pimpl->gem2mt_task->setComment(comment);
}

void MassTransport::setDistributing(size_t node_ndx, size_t system_ndx)
{
    pimpl->gem2mt_task->setDistributing(node_ndx, system_ndx);
}

void MassTransport::setNodeType(size_t node_ndx, long type)
{
    pimpl->gem2mt_task->setNodeType(node_ndx, type);
}

void MassTransport::setSIA(char flag)
{
    pimpl->gem2mt_task->setSIA(flag);
}

void MassTransport::setOutput(bool enable)
{
    pimpl->gem2mt_task->setOutput(enable);
}

void MassTransport::setOutVTK(bool enable)
{
    pimpl->gem2mt_task->setOutVTK(enable);
}

void MassTransport::setTypeFluxPhase(char type)
{
    pimpl->gem2mt_task->setTypeFluxPhase(type);
}

void MassTransport::setTau(double start, double stop, double step)
{
    pimpl->gem2mt_task->setTau(start, stop, step);
}

void MassTransport::setSpatialDimensions(double x, double y, double z)
{
    pimpl->gem2mt_task->setSpatialDimensions(x, y, x);
}

void MassTransport::setTimeStepReductionFactor(double val)
{
    pimpl->gem2mt_task->setTimeStepReductionFactor(val);
}

void MassTransport::setInitialTotalNodeVolume(double val)
{
    pimpl->gem2mt_task->setInitialTotalNodeVolume(val);
}

void MassTransport::setFluidAdvectionVelocity(double val)
{
    pimpl->gem2mt_task->setFluidAdvectionVelocity(val);
}

void MassTransport::setInitialNodeEffectivePorosity(double val)
{
    pimpl->gem2mt_task->setInitialNodeEffectivePorosity(val);
}

void MassTransport::setInitialEffectivePermeability(double val)
{
    pimpl->gem2mt_task->setInitialEffectivePermeability(val);
}

void MassTransport::setInitialDispersivity(double val)
{
    pimpl->gem2mt_task->setInitialDispersivity(val);
}

void MassTransport::setInitialDiffusivity(double val)
{
    pimpl->gem2mt_task->setInitialDiffusivity(val);
}

void MassTransport::setInitialTortuosityFactor(double val)
{
    pimpl->gem2mt_task->setInitialTortuosityFactor(val);
}

void MassTransport::setCutoffICamount(double val)
{
    pimpl->gem2mt_task->setCutoffICamount(val);
}

void MassTransport::setCutoffMinimalAmountsIC(double val)
{
    pimpl->gem2mt_task->setCutoffMinimalAmountsIC(val);
}

size_t MassTransport::nPhases() const
{
    return pimpl->gem2mt_task->nPhases();
}

size_t MassTransport::nElements() const
{
    return pimpl->gem2mt_task->nElements();
}

size_t MassTransport::nComponents() const
{
    return pimpl->gem2mt_task->nComponents();
}

void MassTransport::setHydraulicParameters(size_t node_ndx, double Vt, double vp, double eps, double Km, double al, double Dif, double nto)
{
    pimpl->gem2mt_task->setHydraulicParameters(node_ndx, Vt, vp, eps, Km, al, Dif, nto);
}

std::string MassTransport::summary() const
{
    return pimpl->gem2mt_task->summary();
}

long int MassTransport::nodeSystem(size_t node_ndx) const
{
    return pimpl->gem2mt_task->nodeSystem(node_ndx);
}

long int MassTransport::nodeType(size_t node_ndx) const
{
    return pimpl->gem2mt_task->nodeType(node_ndx);
}

std::vector<size_t> MassTransport::nodesOfSystem(size_t system_ndx) const
{
    std::vector<size_t> nodes;
    for(size_t q=0; q<nNodes(); ++q)
        if(nodeSystem(q) == static_cast<long int>(system_ndx))
            nodes.push_back(q);
    return nodes;
}

std::vector<size_t> MassTransport::nodesOfType(long int type) const
{
    std::vector<size_t> nodes;
    for(size_t q=0; q<nNodes(); ++q)
        if(nodeType(q) == type)
            nodes.push_back(q);
    return nodes;
}

TransportProperties MassTransport::transportProperties(size_t node_ndx) const
{
    try {
        auto h = pimpl->gem2mt_task->hydraulicParameters(node_ndx);
        TransportProperties tp;
        tp.volume = h[0];  tp.velocity = h[1];  tp.porosity = h[2];  tp.permeability = h[3];
        tp.dispersivity = h[4];  tp.diffusivity = h[5];  tp.tortuosity = h[6];
        return tp;
    }
    catch(TError& err) {
        throw std::runtime_error(err.mess);
    }
}

void MassTransport::setNodes(const std::vector<size_t>& nodes, const NodeSetup& setup)
{
    const auto& t = setup.transport;
    const bool set_transport = t.volume || t.velocity || t.porosity || t.permeability ||
                               t.dispersivity || t.diffusivity || t.tortuosity;
    const bool set_chemistry = !std::holds_alternative<std::monostate>(setup.composition) ||
                               setup.temperature || setup.pressure;

    // Check everything before a node is changed
    for(auto q: nodes)
        if(q >= nNodes())
            throw std::runtime_error("Node index "+std::to_string(q)+" out of range (nodes: "+
                                     std::to_string(nNodes())+")");
    // Permitted ranges: porosity in (0, 1], tortuosity > 0, the others >= 0
    auto check = [](const std::optional<double>& v, const std::string& name, bool positive, double hi) {
        if(!v) return;
        const bool ok = (positive ? *v > 0. : *v >= 0.) && *v <= hi;
        if(!ok) {
            std::ostringstream msg;
            msg << "TransportProperties." << name << " = " << *v << " is not permitted: it must be "
                << (positive ? "> 0" : ">= 0");
            if(std::isfinite(hi)) msg << " and <= " << hi;
            throw std::runtime_error(msg.str());
        }
    };
    const double inf = std::numeric_limits<double>::infinity();
    check(t.volume, "volume", false, inf);
    check(t.porosity, "porosity", true, 1.);
    check(t.permeability, "permeability", false, inf);
    check(t.dispersivity, "dispersivity", false, inf);
    check(t.diffusivity, "diffusivity", false, inf);
    check(t.tortuosity, "tortuosity", true, inf);
    if(set_chemistry && !pimpl->chem_array)
        throw std::runtime_error("Nodes are not initialized: call initNodes first");

    // New element amounts in chemical-system order (no value: keep the composition of each node)
    std::optional<std::vector<double>> b;
    if(const auto* v = std::get_if<std::vector<double>>(&setup.composition)) {
        b = *v;    // the size is checked for each node
    }
    else if(const auto* m = std::get_if<std::map<std::string, double>>(&setup.composition)) {
        const auto& engine = *pimpl->chem_array;
        b = std::vector<double>(engine.numElements(), 0.);
        for(const auto& [name, amount]: *m) {
            auto i = engine.indexElement(name);
            if(i >= engine.numElements())
                throw std::runtime_error("Element "+name+" is not in the chemical system");
            (*b)[i] = amount;
        }
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    try {
        for(auto q: nodes) {
            if(set_transport) {
                auto h = pimpl->gem2mt_task->hydraulicParameters(q);
                pimpl->gem2mt_task->setHydraulicParameters(q,
                    t.volume.value_or(h[0]), t.velocity.value_or(h[1]), t.porosity.value_or(h[2]),
                    t.permeability.value_or(h[3]), t.dispersivity.value_or(h[4]),
                    t.diffusivity.value_or(h[5]), t.tortuosity.value_or(h[6]));
            }
            if(set_chemistry)
                pimpl->gem2mt_task->setNodeBulkComposition(q, b ? &*b : nullptr,
                    setup.temperature.value_or(nan), setup.pressure.value_or(nan));
        }
    }
    catch(TError& err) {
        throw std::runtime_error(err.mess);
    }
}

void MassTransport::setSystem(size_t system_ndx, const NodeSetup& setup)
{
    setNodes(nodesOfSystem(system_ndx), setup);
}

size_t MassTransport::nMGPfluxes() const
{
    return pimpl->gem2mt_task->nMGPfluxes();
}

void MassTransport::setNumberMGPfluxes(long nFD)
{
    pimpl->gem2mt_task->setNumberMGPfluxes(nFD);
}

void MassTransport::setFluxSourceReceive(size_t node_ndx, long source, long receive)
{
    pimpl->gem2mt_task->setFluxSourceReceive(node_ndx, source, receive);
}

void MassTransport::setFluxDefinition(size_t node_ndx, double order, double rate, double quantity, double val)
{
    pimpl->gem2mt_task->setFluxDefinition(node_ndx, order, rate, quantity, val);
}

void MassTransport::setFluxMGPid(size_t node_ndx, const std::string &id)
{
    pimpl->gem2mt_task->setFluxMGPid(node_ndx, id);
}

void MassTransport::setFluxID(size_t node_ndx, const std::string &id)
{
    pimpl->gem2mt_task->setFluxID(node_ndx, id);
}

size_t MassTransport::nPhaseGroups() const
{
    return pimpl->gem2mt_task->nPhaseGroups();
}

void MassTransport::setNumberPhaseGroups(long nMGP)
{
    pimpl->gem2mt_task->setNumberPhaseGroups(nMGP);
}

void MassTransport::setPhaseGroupID(size_t gndx, const std::string &id)
{
    pimpl->gem2mt_task->setPhaseGroupsID(gndx, id);
}

void MassTransport::setUnitsPhaseQuantities(size_t pndx, char units)
{
    pimpl->gem2mt_task->setUnitsPhaseQuantities(pndx, units);
}

void MassTransport::setPhaseGroupQuantities(size_t gndx, size_t pndx, double quantity)
{
    pimpl->gem2mt_task->setPhaseGroupsQuantities(gndx, pndx, quantity);
}

size_t MassTransport::nICsourceFluxes() const
{
    return pimpl->gem2mt_task->nICsourceFluxes();
}

void MassTransport::setNumberICsourceFluxes(long nSFD)
{
    pimpl->gem2mt_task->setNumberICsourceFluxes(nSFD);
}

void MassTransport::setICsourceQuantities(size_t gndx, size_t indx, double quantity)
{
    pimpl->gem2mt_task->setICsourceQuantities(gndx, indx, quantity);
}

size_t MassTransport::nParticleTypes() const
{
    return pimpl->gem2mt_task->nParticleTypes();
}

void MassTransport::setNumberParticleTypes(long num)
{
    pimpl->gem2mt_task->setNumberParticleTypes(num);
}

void MassTransport::setVTKfields(const std::vector<std::pair<int, int> > &vtk_fields)
{
    pimpl->gem2mt_task->setVTKfields(vtk_fields);
}

void MassTransport::setMonitoring(LogFunction_f new_funct) {

    ScriptFunction_f log_funct =  [this, new_funct](TNodeArray* node_data, long int step, double time) {
        if(pimpl->chem_array) {
            new_funct(*pimpl->chem_array, step, time);
        }
    };
    pimpl->gem2mt_task->setScriptFunc(log_funct);
}

void MassTransport::setParticle(size_t pndx, long pmean, long pmin, long pmax, const std::array<long, 6> &pparam)
{
    pimpl->gem2mt_task->setParticle(pndx, pmean, pmin, pmax, pparam);
}

std::shared_ptr<ChemicalEngineArray> MassTransport::getChemicalEngineArray()
{
     return pimpl->chem_array;
}

std::string MassTransport::dump() const
{
    std::stringstream ss;
    io_formats::SimdJsonWrite out_format(ss, "", true);
    pimpl->gem2mt_task->to_text_file( out_format, true, false );
    return ss.str();
}

std::ostream& operator<<(std::ostream& out, const MassTransport& engine)
{
    auto flags = out.flags();

    out << engine.summary();

    // Recover the previous state of `out`
    out.flags(flags);
    return out;
}

} // namespace xGEM2MT
