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

GEMFailures MassTransport::gemFailures() const
{
    const auto& f = pimpl->gem2mt_task->gem_failures;
    GEMFailures r;
    r.count = f.count;
    r.nodes.assign(f.nodes.begin(), f.nodes.end());
    r.first_step = f.first_step;
    r.first_time = f.first_time;
    r.first_node = f.first_node;
    r.first_message = f.first_message;
    return r;
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

std::string MassTransport::summary(size_t width, bool phases) const
{
    return pimpl->gem2mt_task->summary(width, phases);
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
    const bool has_composition = !std::holds_alternative<std::monostate>(setup.composition);
    const bool has_pore_water = !std::holds_alternative<std::monostate>(setup.pore_water);
    const bool has_solids = !std::holds_alternative<std::monostate>(setup.solids);
    const bool scaled = has_pore_water || has_solids || !setup.parts.empty();
    const bool set_chemistry = has_composition || scaled || setup.temperature || setup.pressure;

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

    // Element amounts in chemical-system order (no value: not given)
    auto to_vector = [&](const Composition& c, const std::string& what) -> std::optional<std::vector<double>> {
        std::optional<std::vector<double>> v;
        if(const auto* vec = std::get_if<std::vector<double>>(&c))
            v = *vec;
        else if(const auto* m = std::get_if<std::map<std::string, double>>(&c)) {
            const auto& engine = *pimpl->chem_array;
            v = std::vector<double>(engine.numElements(), 0.);
            for(const auto& [name, amount]: *m) {
                auto i = engine.indexElement(name);
                if(i >= engine.numElements())
                    throw std::runtime_error("Element "+name+" ("+what+") is not in the chemical system");
                (*v)[i] = amount;
            }
        }
        if(v && v->size() != pimpl->chem_array->numElements())
            throw std::runtime_error("The "+what+" has "+std::to_string(v->size())+" values, the chemical system has "+
                                     std::to_string(pimpl->chem_array->numElements())+" elements");
        return v;
    };
    if(has_composition && scaled)
        throw std::runtime_error("Give either composition, or pore_water and parts (solids), not both");
    if(scaled && !has_pore_water)
        throw std::runtime_error("pore_water is needed: it is the reference of the node (water_mass kg of water)");
    if(setup.water_mass && !(*setup.water_mass > 0.))
        throw std::runtime_error("water_mass must be > 0");
    const double water_saturation = setup.water_saturation.value_or(1.);
    if(!(water_saturation > 0. && water_saturation <= 1.))
        throw std::runtime_error("water_saturation must be > 0 and <= 1");
    const auto b = to_vector(setup.composition, "composition");             // as given (no scaling)
    const auto b_water = to_vector(setup.pore_water, "pore water");         // per kg of water

    // The parts, each normalised to 1 kg (element molar masses in kg/mol)
    struct PartIn { std::vector<double> b; bool solid; bool by_mass; double value; std::optional<double> density; };
    std::vector<PartIn> parts;
    const auto mm = pimpl->chem_array->elementMolarMasses();
    auto per_kg = [&](const Composition& c, const std::string& what) {
        auto v = to_vector(c, what);
        double mass = 0.;
        for(size_t i=0; i<v->size(); ++i) mass += (*v)[i] * mm[i];
        if(!(mass > 0.))
            throw std::runtime_error("The "+what+" has no mass");
        for(auto& x: *v) x /= mass;
        return *v;
    };
    if(has_solids) {
        if(setup.grain_density && !(*setup.grain_density > 0.))
            throw std::runtime_error("grain_density must be > 0");
        parts.push_back({per_kg(setup.solids, "solids"), true, false, 1., setup.grain_density});
    }
    for(size_t k=0; k<setup.parts.size(); ++k) {
        const auto& pt = setup.parts[k];
        const std::string what = "part " + std::to_string(k);
        if(int(bool(pt.fraction)) + int(bool(pt.mass)) != 1)
            throw std::runtime_error("Give either fraction or mass for "+what);
        if(pt.of != "pore" && pt.of != "solid")
            throw std::runtime_error("Give of=\"pore\" or of=\"solid\" for "+what);
        if(std::holds_alternative<std::monostate>(pt.material))
            throw std::runtime_error("No material for "+what);
        if(pt.density && !(*pt.density > 0.))
            throw std::runtime_error("density must be > 0 ("+what+")");
        if(pt.fraction && !(*pt.fraction >= 0. && *pt.fraction <= 1.))
            throw std::runtime_error("fraction must be 0 to 1 ("+what+")");
        if(pt.mass && !(*pt.mass >= 0.))
            throw std::runtime_error("mass must be >= 0 ("+what+")");
        parts.push_back({per_kg(pt.material, what), pt.of == "solid", bool(pt.mass),
                         pt.mass ? *pt.mass : *pt.fraction, pt.density});
    }
    double pore_sum = water_saturation, solid_sum = 0.;
    bool any_solid_fraction = false;
    for(const auto& pt: parts) {
        if(pt.by_mass) continue;
        if(pt.solid) { solid_sum += pt.value; any_solid_fraction = true; }
        else pore_sum += pt.value;
    }
    if(pore_sum > 1. + 1e-9)
        throw std::runtime_error("water_saturation and the pore fractions add up to more than 1");
    if(any_solid_fraction && std::fabs(solid_sum - 1.) > 1e-6)
        throw std::runtime_error("The solid fractions must add up to 1");
    const double water_mass = setup.water_mass.value_or(1.);

    // Each node is changed completely or not at all: if the chemistry fails, the previous
    // transport properties are restored (setNodeBulkComposition restores the chemistry itself)
    const double nan = std::numeric_limits<double>::quiet_NaN();
    try {
        for(auto q: nodes) {
            std::array<double, 7> old_h{};
            if(set_transport) {
                old_h = pimpl->gem2mt_task->hydraulicParameters(q);
                const auto& h = old_h;
                pimpl->gem2mt_task->setHydraulicParameters(q,
                    t.volume.value_or(h[0]), t.velocity.value_or(h[1]), t.porosity.value_or(h[2]),
                    t.permeability.value_or(h[3]), t.dispersivity.value_or(h[4]),
                    t.diffusivity.value_or(h[5]), t.tortuosity.value_or(h[6]));
            }
            if(set_chemistry) {
                try {
                    // Scaled composition from the porosity of this node (after its transport values are set)
                    std::optional<std::vector<double>> b_node = b;
                    if(scaled) {
                        // Volumes from GEM at the T and P of the node (L per kg; the pore water per kg of water)
                        const double T = setup.temperature.value_or(pimpl->chem_array->temperature(q));
                        const double P = setup.pressure.value_or(pimpl->chem_array->pressure(q));
                        // Volume per kg of a part: the volume that 1 kg of it adds to the pore water (1 kg of
                        // water), at the T and P of the node. A solid or a gas does not converge alone.
                        const double pore_water_volume = pimpl->gem2mt_task->materialVolume(*b_water, T, P);   // L per kg of water
                        auto volume_per_kg = [&](const std::vector<double>& bk) {
                            std::vector<double> mix(*b_water);
                            for(size_t i=0; i<mix.size(); ++i) mix[i] += bk[i];
                            const double v = pimpl->gem2mt_task->materialVolume(mix, T, P) - pore_water_volume;
                            if(!(v > 0.))
                                throw std::runtime_error("A part adds no volume to the pore water (it dissolves?): give its density");
                            return v;
                        };
                        const double porosity = pimpl->gem2mt_task->hydraulicParameters(q)[2];
                        const double water_volume = water_mass * pore_water_volume;                // L
                        const double pore_volume = water_volume / water_saturation;
                        const double solid_volume = pore_volume * (1. - porosity) / porosity;
                        // Mass (kg) of each part. Parts given by mass take their volume from their `of` volume;
                        // the solid fractions are shares of the solid volume that is left.
                        std::vector<double> kg(parts.size(), 0.), target(parts.size(), 0.), vpk(parts.size(), 0.);
                        double mass_pore = 0., mass_solid = 0.;
                        for(size_t k=0; k<parts.size(); ++k) {
                            const auto& pt = parts[k];
                            vpk[k] = pt.density ? 1. / *pt.density : volume_per_kg(pt.b);         // L/kg
                            if(pt.by_mass) {
                                kg[k] = pt.value;
                                (pt.solid ? mass_solid : mass_pore) += kg[k] * vpk[k];
                            }
                        }
                        if(any_solid_fraction && mass_solid > solid_volume * (1. + 1e-9))
                            throw std::runtime_error("The solid parts given by mass need more than the solid volume of node "+
                                                     std::to_string(q));
                        if(pore_sum + mass_pore / pore_volume > 1. + 1e-9)
                            throw std::runtime_error("The pore parts need more than the pore volume of node "+std::to_string(q));
                        for(size_t k=0; k<parts.size(); ++k) {
                            const auto& pt = parts[k];
                            if(pt.by_mass) continue;
                            target[k] = pt.solid ? (solid_volume - mass_solid) * pt.value : pore_volume * pt.value;  // L
                            kg[k] = target[k] / vpk[k];
                        }
                        auto compose = [&]() {
                            std::vector<double> bn(pimpl->chem_array->numElements(), 0.);
                            for(size_t i=0; i<bn.size(); ++i)
                                bn[i] = (*b_water)[i] * water_mass;
                            for(size_t k=0; k<parts.size(); ++k)
                                for(size_t i=0; i<bn.size(); ++i)
                                    bn[i] += parts[k].b[i] * kg[k];
                            return bn;
                        };
                        b_node = compose();

                        // A gas (or another fluid) can dissolve in part in the pore water, so that its volume in
                        // the node is smaller than planned. Then add the missing volume to the fluid parts and
                        // calculate again, until the node volume is pore volume + solid volume.
                        double fluid_target = 0.;
                        for(size_t k=0; k<parts.size(); ++k)
                            if(!parts[k].solid && !parts[k].by_mass && !parts[k].density) fluid_target += target[k];
                        if(fluid_target > 0.) {
                            const double node_target = pore_volume + solid_volume;     // L
                            bool converged = false;
                            for(int it=0; it<30 && !converged; ++it) {
                                pimpl->gem2mt_task->setNodeBulkComposition(q, &*b_node,
                                    setup.temperature.value_or(nan), setup.pressure.value_or(nan));
                                const double missing = node_target - pimpl->chem_array->systemVolume(q) * 1e3;
                                if(std::fabs(missing) <= 1e-5 * node_target) { converged = true; break; }
                                for(size_t k=0; k<parts.size(); ++k)
                                    if(!parts[k].solid && !parts[k].by_mass && !parts[k].density)
                                        kg[k] = std::max(0., kg[k] + missing * (target[k] / fluid_target) / vpk[k]);
                                b_node = compose();
                            }
                            if(!converged)
                                throw std::runtime_error("The volume of the fluid parts of node "+std::to_string(q)+
                                                         " did not converge (does the gas dissolve completely?)");
                        }
                    }
                    pimpl->gem2mt_task->setNodeBulkComposition(q, b_node ? &*b_node : nullptr,
                        setup.temperature.value_or(nan), setup.pressure.value_or(nan));
                }
                catch(...) {
                    if(set_transport)
                        pimpl->gem2mt_task->setHydraulicParameters(q, old_h[0], old_h[1], old_h[2],
                                                                   old_h[3], old_h[4], old_h[5], old_h[6]);
                    throw;
                }
            }
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
