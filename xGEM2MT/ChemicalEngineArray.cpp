// xGEMS2MT is a C++ and Python library for mass transport modeling with Gibbs energy minimization
//
// Copyright (C) 2022-2026 S.Dmytriieva, Dmitrii Kulik, G.D. Miron
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

#include "ChemicalEngineArray.hpp"
#include "MassTransport_p.h"

// C++ includes
//#include <chrono>
//#include <iomanip>
//#include <memory>

// GEMS3K includes
#include <GEMS3K/nodearray.h>
#include <GEMS3K/gems3k_impex.h>
#include <GEMS3K/jsonconfig.h>

namespace xGEM2MT {

void update_loggers( bool use_cout, const std::string& logfile_name, size_t log_level)
{
    gemsSettings().gems3k_update_loggers(use_cout, logfile_name, log_level);
}

struct ChemicalEngineArray::Impl
{
    /// Construct a default Impl instance
    Impl(size_t nodes): n_nodes(nodes)
    {}

    /// The options for ChemicalEngineArray
    ChemicalEngineArrayOptions options;

    /// Number of allocated nodes
    size_t n_nodes;

    /// The TNodeArray instance from ChemicalEngineArray
    std::shared_ptr<TNodeArray> node_array;

    /// Write IPM, DCH and DBR files (in binary, txt or json mode)
    GEMS3KGenerator::IOModes io_mode = GEMS3KGenerator::f_key_value;

    /// The formula matrix of the species
    Matrix formula_matrix;

    /// Molar masses of elements
    Vector elemMolMasses;

    /// Molar masses of species
    Vector molMasses;

};

ChemicalEngineArray::ChemicalEngineArray(size_t nodes)
    : pimpl(new Impl(nodes))
{}

ChemicalEngineArray::ChemicalEngineArray(const MassTransport& mt_class)
    : pimpl(new Impl(mt_class.pimpl->gem2mt_task->nodeArray()->nNodes()))
{
    pimpl->io_mode = GEMS3KGenerator::f_json;
    pimpl->node_array = mt_class.pimpl->gem2mt_task->nodeArray();
    reallocateEngineArrays();
}

ChemicalEngineArray::~ChemicalEngineArray()
{}

auto ChemicalEngineArray::initialize(std::string ipm_filename, std::string dbr_filename, std::vector<long int> nodeTypes) -> void
{
    // Allocate memory for the GEMS `node` member
    pimpl->node_array = TNodeArray::create(pimpl->n_nodes);

    // Initialize the GEMS `nodes` member
    size_t res = pimpl->node_array->GEM_init(ipm_filename.c_str(),
                                             dbr_filename.empty() ? nullptr : dbr_filename.c_str(),
                                             nodeTypes.empty() ? nullptr : nodeTypes.data(),
                                             false );

    // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the ChemicalEngineArray object.\n"
                                 "There was a problem during memory allocation.");

    // Check if there was a file read error during node initialization
    if(res == 1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the ChemicalEngineArray object.\n"
                                 "Make sure the provided file exists relative to the working directory.");

    // Initialize the GEMS input files mode
    GEMS3KGenerator generator_data(ipm_filename.c_str());
    pimpl->io_mode =   generator_data.files_mode();
    reallocateEngineArrays();
}

auto ChemicalEngineArray::initializeFromJsonStrings( std::string dch_json, std::string ipm_json,
                                                    std::vector<std::string> dbr_json, std::vector<long int> nodeTypes) -> void
{
    pimpl->io_mode = GEMS3KGenerator::f_json;
    // Allocate memory for the GEMS `node` member
    pimpl->node_array = TNodeArray::create(pimpl->n_nodes);

    // Initialize the GEMS `node` member
    size_t res = pimpl->node_array->GEM_init( dch_json, ipm_json, dbr_json, nodeTypes.empty() ? nullptr : nodeTypes.data() );

    // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the ChemicalEngineArray object from JSON strings. \n"
                                 " There was a problem during memory allocation.");

    // Check if there was a JSON parsing error during node initialization
    if(res == 1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not initialize the ChemicalEngineArray object from JSON strings.\n"
                                 "Make sure that these JSON strings are not empty and are in correct format.");

    reallocateEngineArrays();
}

long ChemicalEngineArray::reequilibrate(size_t ndx, bool warmstart)
{
    size_t  Mode = warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    TNode workNode(*pimpl->node_array->getCalcNode());
    return pimpl->node_array->RunGEM( &workNode, ndx, Mode, pimpl->node_array->pNodT1() );
}

void ChemicalEngineArray::reequilibrate(bool warmstart)
{
    size_t  Mode = warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    std::vector<long int> nodeFlags, retCodes;
    for( size_t ii=0; ii<pimpl->node_array->nNodes(); ++ii)
    {
        nodeFlags.push_back(1);
        retCodes.push_back(0);
    }
    pimpl->node_array->RunGEM( Mode, pimpl->node_array->nNodes(), pimpl->node_array->pNodT1(), nodeFlags.data(), retCodes.data() );
}

bool ChemicalEngineArray::equilibrate(size_t ndx)
{
    TestModeGEMParam  mode_param(pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA,
                                pimpl->options.useSIA,
                                pimpl->options.step,
                                pimpl->options.cdv,
                                pimpl->options.cez);
    return pimpl->node_array->CalcIPM_One( mode_param, ndx, nullptr );
}

bool ChemicalEngineArray::equilibrate()
{
    TestModeGEMParam  mode_param(pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA,
                                pimpl->options.useSIA,
                                pimpl->options.step,
                                pimpl->options.cdv,
                                pimpl->options.cez);
    return pimpl->node_array->CalcIPM_List( mode_param, 0, pimpl->node_array->nNodes()-1, nullptr );
}

void ChemicalEngineArray::dumpGEMS3K(const std::string &filepath, bool put_all) const
{
    try
    {
        ProcessProgressFunction messageF = [](const std::string& , long ){
            //std::cout << "TProcess GEM3k output" <<  message.c_str() << point << std::endl;
            return false;
        };
        pimpl->node_array->genGEMS3KInputFiles(filepath, messageF, pimpl->node_array->nNodes(),
                                               pimpl->io_mode, false, true, put_all, false );
    }
    catch(TError& err)
    {
        std::cout  << err.title << err.mess << std::endl;
    }
    catch(std::exception& e)
    {
        std::cout  << "std::exception: " << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout  << "unknown exception" << std::endl;
    }
}

void ChemicalEngineArray::nextTimePoint(bool always)
{
    DATACH* CH = pimpl->node_array->pCSD();       // DataCH structure
    DATABRPTR* C0 = pimpl->node_array->pNodT0();  // nodes at current time point
    DATABRPTR* C1 = pimpl->node_array->pNodT1();  // nodes at previous time point
    auto nnodes = pimpl->node_array->nNodes();


    for (size_t ii=0; ii<nnodes; ii++)
    {
        bool NeedCopy = always;
        if(!NeedCopy) {
            for(size_t ic=0; ic < CH->nICb-1; ic++) // Important!!!! use pm.E not 1 for systems without aqua( do not check charge)
            {
                auto dc = C0[ii]->bIC[ic] - C1[ii]->bIC[ic];
                if( fabs( dc ) > std::min( pimpl->options.cdv, (C1[ii]->bIC[ic] * 1e-3)))
                    NeedCopy = true;
            }
        }
        if( NeedCopy )
            pimpl->node_array->CopyNodeFromTo( ii, nnodes, C1, C0 );
    }
}

auto ChemicalEngineArray::reallocateEngineArrays() -> void
{
    // Assemble the formula matrix of the species
    const size_t E = numElements();
    const size_t N = numSpecies();
    pimpl->formula_matrix = Matrix::Zero(E, N);
    for(size_t i = 0; i < N; ++i)
        for(size_t j = 0; j < E; ++j)
            pimpl->formula_matrix(j, i) = pimpl->node_array->DCaJI(i, j);

    // Allocate memory for vector members
    pimpl->elemMolMasses.resize(E);
    pimpl->molMasses.resize(N);
}

auto ChemicalEngineArray::options() const -> const ChemicalEngineArrayOptions&
{
    return pimpl->options;
}

auto ChemicalEngineArray::setOptions(const ChemicalEngineArrayOptions& options) -> void
{
    pimpl->options = options;
}

auto ChemicalEngineArray::setWarmStart() -> void
{
    pimpl->options.warmstart = true;
}

auto ChemicalEngineArray::setColdStart() -> void
{
    pimpl->options.warmstart = false;
}

size_t ChemicalEngineArray::nNodes() const
{
    return pimpl->n_nodes;
}

auto ChemicalEngineArray::readDbrFromFile(size_t  ndx, std::string filename) -> void
{
    // Reads another dbr file with input system composition
    size_t res = pimpl->node_array->GEMS3k_read_dbr(ndx, filename.c_str(), pimpl->io_mode);

    // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not find the provided dbr file.\n"
                                 "There was a problem during memory allocation.");

    // Check if there was a file read error during node initialization
    if(res >= 1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not find the provided dbr file.\n"
                                 "Make sure the provided file path exists relative to the working directory.");
}

auto ChemicalEngineArray::writeDbrToFile(size_t  ndx, std::string filename) -> void
{
    // Write current node into a dbr file into path filename
    pimpl->node_array->GEMS3k_write_dbr( ndx, filename.c_str(), pimpl->io_mode);
}

auto ChemicalEngineArray::readDbrFromJsonString(size_t  ndx, std::string dbr_json) -> void
{
    // Reads another dbr file with input system composition
    size_t res = pimpl->node_array->GEMS3k_read_dbr( ndx, dbr_json, true );
    // res = 1;    // Temporary plug

    // Check if there was a system error during node initialization
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not process the provided DBR JSON string.\n"
                                 "There was a problem during memory allocation.");

    // Check if there was a JSON string read error during node initialization
    if(res >= 1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not process the provided DBR JSON string.\n"
                                 "Make sure that this JSON string is not empty and is in correct format and parameter position.");
}

auto ChemicalEngineArray::writeDbrToJsonString(size_t  ndx) -> const std::string
{
    // Writes another node dbr with input system composition into json string dbr_json
    std::string dbr_json = "";
    size_t res = pimpl->node_array->GEMS3k_write_dbr( ndx, dbr_json );

    // Check if there was a json write error
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not process the provided DBR into JSON string.\n"
                                 "Check the code.");

    return dbr_json;
}

auto ChemicalEngineArray::numElements() const -> size_t
{
    return pimpl->node_array->pCSD()->nIC;
}

auto ChemicalEngineArray::numSpecies() const -> size_t
{
    return pimpl->node_array->pCSD()->nDC;
}

auto ChemicalEngineArray::numPhases() const -> size_t
{
    return pimpl->node_array->pCSD()->nPH;
}

auto ChemicalEngineArray::numSpeciesInPhase_i(size_t iphase) const -> size_t
{
    return pimpl->node_array->pCSD()->nDCinPH[iphase];
}

auto ChemicalEngineArray::numSpeciesInPhase(size_t iphase) const -> size_t
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return pimpl->node_array->pCSD()->nDCinPH[iphase];
    }
    return 0;
}

auto ChemicalEngineArray::numSpeciesInPhase(std::string phase) const -> size_t
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return pimpl->node_array->pCSD()->nDCinPH[iphase];
    }
    return 0;
}

auto ChemicalEngineArray::elementName_i(size_t ielement) const -> std::string
{
    return pimpl->node_array->pCSD()->ICNL[ielement];
}

auto ChemicalEngineArray::elementName(size_t ielement) const -> std::string
{
    if(pimpl->node_array->getCalcNode()->check_IC_xCH(ielement)) {
        return pimpl->node_array->pCSD()->ICNL[ielement];
    }
    return {};
}

auto ChemicalEngineArray::speciesName_i(size_t ispecies) const -> std::string
{
    return pimpl->node_array->pCSD()->DCNL[ispecies];
}

auto ChemicalEngineArray::speciesName(size_t ispecies) const -> std::string
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        return pimpl->node_array->pCSD()->DCNL[ispecies];
    }
    return {};
}

auto ChemicalEngineArray::speciesCharge_i(size_t ispecies) const -> double
{
    MatrixConstRef W = formulaMatrix();
    auto chargeRow = numElements();
    return W(chargeRow - 1, ispecies);
}

auto ChemicalEngineArray::speciesCharge(size_t ispecies) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        return speciesCharge_i(ispecies);
    }
    return 0;
}

auto ChemicalEngineArray::speciesCharge(std::string species, std::optional<std::string> phase) const -> double
{
    auto ispecies = indexSpecies(species, phase);
    if(ispecies<numSpecies()) {
        return speciesCharge_i(ispecies);
    }
    return 0;
}

auto ChemicalEngineArray::phaseName_i(size_t iphase) const -> std::string
{
    return pimpl->node_array->pCSD()->PHNL[iphase];
}

auto ChemicalEngineArray::phaseName(size_t iphase) const -> std::string
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return pimpl->node_array->pCSD()->PHNL[iphase];
    }
    return {};
}

auto ChemicalEngineArray::setSpeciesUpperLimit_i(size_t ndx, size_t ispecies, double amount) -> void
{
    double bound = (amount < 0.0 ? 1e6 : amount);
    pimpl->node_array->pNodT1()[ndx]->dul[ispecies] = bound;
}

auto ChemicalEngineArray::setSpeciesLowerLimit_i(size_t ndx, size_t ispecies, double amount) -> void
{
    double bound = (amount < 0.0 ? 0.0 : amount);
    pimpl->node_array->pNodT1()[ndx]->dll[ispecies] = bound;
}

// These methods may be ambiguous, as in GEMS3K, species with the same name
//   may occur in more than one condensed phase!
// An overload including phase name or index is needed!
//
auto ChemicalEngineArray::setSpeciesUpperLimit(size_t ndx, std::string name, double amount, std::optional<std::string> phase) -> void
{
    auto ispecies = indexSpecies(name, phase);
    if(ispecies<numSpecies()) {
        setSpeciesUpperLimit_i(ndx, ispecies, amount);
    }
}

auto ChemicalEngineArray::setSpeciesLowerLimit(size_t ndx, std::string name, double amount, std::optional<std::string> phase) -> void
{
    auto ispecies = indexSpecies(name, phase);
    if(ispecies<numSpecies()) {
        setSpeciesLowerLimit_i(ndx, ispecies, amount);
    }
}

auto ChemicalEngineArray::setSpeciesUpperLimit(size_t ndx, size_t ispecies, double amount) -> void
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        setSpeciesUpperLimit_i(ndx, ispecies, amount);
    }
}

auto ChemicalEngineArray::setSpeciesLowerLimit(size_t ndx, size_t ispecies, double amount) -> void
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        setSpeciesLowerLimit_i(ndx, ispecies, amount);
    }
}

auto ChemicalEngineArray::speciesUpperLimits(size_t ndx) const -> Vector
{
    return Vector::Map(pimpl->node_array->pNodT1()[ndx]->dul, numSpecies());
}
auto ChemicalEngineArray::speciesLowerLimits(size_t ndx) const -> Vector
{
    return Vector::Map(pimpl->node_array->pNodT1()[ndx]->dll, numSpecies());
}

auto ChemicalEngineArray::setSpeciesUpperLimits(size_t ndx, VectorConstRef n) -> void
{
    double bound;
    for(size_t jj = 0; jj <numSpecies(); jj++) {
        bound = (n[jj] < 0.0? 1e6: n[jj]);
        pimpl->node_array->pNodT1()[ndx]->dul[jj] = bound;
    }
}

auto ChemicalEngineArray::setSpeciesLowerLimits(size_t ndx, VectorConstRef n) -> void
{
    double bound;
    for(size_t jj = 0; jj <numSpecies(); jj++) {
        bound = (n[jj] < 0.0? 0.0: n[jj]);
        pimpl->node_array->pNodT1()[ndx]->dll[jj] = bound;;
    }
}

auto ChemicalEngineArray::setStandardMolarGibbsEnergy_i(size_t ispecies, double TK, double P, double value) -> void
{
    auto idc = pimpl->node_array->getCalcNode()->DC_xDB_to_xCH(ispecies);
    pimpl->node_array->Set_DC_G0(idc, P, TK, value);
}

auto ChemicalEngineArray::setStandardMolarGibbsEnergy(size_t ispecies, double TK, double P, double value) -> void
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        setStandardMolarGibbsEnergy_i(ispecies, TK, P, value);
    }
}

auto ChemicalEngineArray::setStandardMolarGibbsEnergy(std::string name, double TK, double P, double value, std::optional<std::string> phase) -> void
{
    auto ispecies = indexSpecies(name, phase);
    if(ispecies<numSpecies()) {
        setStandardMolarGibbsEnergy_i(ispecies, TK, P, value);
    }
}

auto ChemicalEngineArray::indexElement(std::string element) const -> size_t
{
    auto ielement = pimpl->node_array->getCalcNode()->IC_name_to_xCH(element);
    return (ielement>=0 ? ielement : numElements()); // in case that element name was not found
}

auto ChemicalEngineArray::indexSpecies(std::string species, std::optional<std::string> phase) const -> size_t
{
    int ispecies = -1;
    if(phase.has_value()) {
        auto species_map = pimpl->node_array->getCalcNode()->DC_name_to_xCH_map(species);
        if(species_map.find(phase.value()) != species_map.end()) {
            return  species_map[phase.value()];
        }
    }
    else {
        ispecies =  pimpl->node_array->getCalcNode()->DC_name_to_xCH(species);
    }
    return (ispecies>=0 ? ispecies : numSpecies()); // in case that element name was not found
}

std::map<std::string, long> ChemicalEngineArray::indexSpeciesMap(std::string species) const
{
    return pimpl->node_array->getCalcNode()->DC_name_to_xCH_map(species);
}

auto ChemicalEngineArray::indexSpeciesAll(std::string species) const -> VectorXi
{
    std::vector<Index> occur;
    occur.resize(0);
    const Index size = numSpecies();
    for (Index i = 0; i < size; ++i)
        if (speciesName_i(i) == species)
            occur.push_back(i);
    Index n_species_found = occur.size();
    VectorXi result;
    result.resize(0);
    if (n_species_found > 0)
    {
        result.resize(n_species_found);
        for (Index is = 0; is < n_species_found; is++)
            result(is) = occur.at(is);
    }
    return result;
}

// Index of phase searched by name
auto ChemicalEngineArray::indexPhase(std::string phase) const -> size_t
{
    auto iphase = pimpl->node_array->getCalcNode()->Ph_name_to_xCH(phase);
    return (iphase>=0 ? iphase : numPhases()); // in case that element name was not found
}

auto ChemicalEngineArray::indexPhaseAll(std::string phase) const -> VectorXi
{
    std::vector<Index> occur;
    occur.resize(0);
    const Index size = numPhases();
    for (Index i = 0; i < size; ++i)
        if (phaseName_i(i) == phase)
            occur.push_back(i);
    Index n_phases_found = occur.size();
    VectorXi result;
    result.resize(0);
    if (n_phases_found > 0)
    {
        result.resize(n_phases_found);
        for (Index is = 0; is < n_phases_found; is++)
            result(is) = occur.at(is);
    }
    return result;
}

auto ChemicalEngineArray::indexPhaseWithSpecies(size_t ispecies) const -> size_t
{
    const size_t size = numPhases();
    size_t counter = 0;
    for(size_t i = 0; i < size; ++i) {
        counter += numSpeciesInPhase(i);
        if(counter > ispecies) {
            return i;
        }
    }
    return size; // in case that species name was not found
}

auto ChemicalEngineArray::indexFirstSpeciesInPhase_i(size_t iphase) const -> size_t
{
    size_t counter = 0;
    for(size_t i = 0; i < iphase; ++i)
        counter += numSpeciesInPhase_i(i);
    return counter;
}

auto ChemicalEngineArray::indexFirstSpeciesInPhase(size_t iphase) const -> size_t
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return indexFirstSpeciesInPhase_i(iphase);
    }
    return numSpecies();
}

auto ChemicalEngineArray::indexFirstSpeciesInPhase(std::string phase) const -> size_t
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return indexFirstSpeciesInPhase_i(iphase);
    }
    return numSpecies();
}

auto ChemicalEngineArray::elementMolarMasses() const -> VectorConstRef
{
    for(size_t i = 0; i < numElements(); ++i)
        pimpl->elemMolMasses[i] = pimpl->node_array->pCSD()->ICmm[i];
    return pimpl->elemMolMasses;
}

auto ChemicalEngineArray::speciesMolarMasses() const -> VectorConstRef
{
    for(size_t i = 0; i < numSpecies(); ++i)
        pimpl->molMasses[i] = pimpl->node_array->pCSD()->DCmm[i];
    return pimpl->molMasses;
}

auto ChemicalEngineArray::formulaMatrix() const -> MatrixConstRef
{
    return pimpl->formula_matrix;
}

auto ChemicalEngineArray::setPT(size_t ndx, double P, double T) const -> bool
{
    // Check if P, T are feasible
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    if(pimpl->node_array->getCalcNode()->check_TP( T, P )==false)
        return false;
    // Set temperature and pressure
    pimpl->node_array->pNodT1()[ndx]->TK = T;
    pimpl->node_array->pNodT1()[ndx]->P = P;
    return true;
}

auto ChemicalEngineArray::setB(size_t ndx, VectorConstRef b ) -> void
{
    // Set the mole amounts of elements
    for(size_t i = 0; i < numElements(); ++i)
        pimpl->node_array->pNodT1()[ndx]->bIC[i] = b[i];
}


auto ChemicalEngineArray::speciesAmount_i(size_t ndx, size_t ispecies) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->xDC[ispecies];
}

auto ChemicalEngineArray::speciesAmount(size_t ndx, size_t ispecies) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        return speciesAmount_i(ndx, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::speciesAmount(size_t ndx, std::string name, std::optional<std::string> phase) const -> double
{
    auto ispecies = indexSpecies(name, phase);
    if(ispecies<numSpecies()) {
        return speciesAmount_i(ndx, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::speciesAmounts(size_t ndx) const -> Vector
{
    return Vector::Map(pimpl->node_array->pNodT1()[ndx]->xDC, numSpecies());
}

auto ChemicalEngineArray::setCorrectionbIC(size_t ndx) -> void
{
    // Correction of the bIC vector
    size_t ii, jj;
    for( ii=0; ii<numElements(); ii++ )
        pimpl->node_array->pNodT1()[ndx]->bIC[ii] = 0.0;
    MatrixConstRef W = formulaMatrix();
    auto SA = speciesAmounts(ndx);
    for(jj = 0; jj <numSpecies(); jj++)
        for( ii=0; ii<numElements(); ii++ )
            pimpl->node_array->pNodT1()[ndx]->bIC[ii] += SA[jj] * W(ii , jj);
}

auto ChemicalEngineArray::setSpeciesAmounts(size_t ndx, VectorConstRef n) -> void
{
    // Updates all the species amount
    size_t jj;
    for(jj = 0; jj <numSpecies(); jj++)
        pimpl->node_array->pNodT1()[ndx]->xDC[jj] = n[jj];
    setCorrectionbIC(ndx);
}

auto ChemicalEngineArray::setSpeciesAmount_i(size_t ndx, size_t ispecies, double amount) -> void
{
    // Updates the species amount
    pimpl->node_array->pNodT1()[ndx]->xDC[ispecies] = amount;
    setCorrectionbIC(ndx);
}

auto ChemicalEngineArray::setSpeciesAmount(size_t ndx, std::string name, double amount, std::optional<std::string> phase) -> void
{
    auto ispecies = indexSpecies(name, phase);
    if(ispecies<numSpecies()) {
        setSpeciesAmount_i(ndx, ispecies, amount);
    }
}

auto ChemicalEngineArray::setSpeciesAmount(size_t ndx, size_t ispecies, double amount) -> void
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        setSpeciesAmount_i(ndx, ispecies, amount);
    }
}

auto ChemicalEngineArray::converged(size_t ndx) const -> bool
{
    const auto status = pimpl->node_array->pNodT1()[ndx]->NodeStatusCH;
    return status == OK_GEM_AIA || status == OK_GEM_SIA;
}

auto ChemicalEngineArray::numIterations(size_t ndx) const -> size_t
{
    return pimpl->node_array->pNodT1()[ndx]->IterDone;
}

auto ChemicalEngineArray::temperature(size_t ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->TK;
}

auto ChemicalEngineArray::pressure(size_t ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->P;
}

auto ChemicalEngineArray::elementAmounts(size_t ndx) const -> Vector
{
    return Vector::Map(pimpl->node_array->pNodT1()[ndx]->bIC, numElements());
}

auto ChemicalEngineArray::elementAmountsInPhase_i(size_t ndx, size_t iphase) const -> Vector
{
    MatrixConstRef W = formulaMatrix();
    VectorConstRef n = speciesAmounts(ndx);
    const size_t first = indexFirstSpeciesInPhase_i(iphase);
    const size_t size = numSpeciesInPhase_i(iphase);
    MatrixConstRef Wp = W.middleCols(first, size);
    VectorConstRef np = n.segment(first, size);
    Vector res = Wp * np;
    return res;
}

auto ChemicalEngineArray::elementAmountsInPhase(size_t ndx, size_t iphase) const -> Vector
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return elementAmountsInPhase_i(ndx, iphase);
    }
    return {};
}

auto ChemicalEngineArray::elementAmountsInPhase(size_t ndx, std::string phase) const -> Vector
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return elementAmountsInPhase_i(ndx, iphase);
    }
    return {};
}

auto ChemicalEngineArray::elementAmountsInSpecies(size_t ndx, VectorXiConstRef ispecies) const -> Vector
{
    MatrixConstRef W = formulaMatrix();
    VectorConstRef n = speciesAmounts(ndx);
    MatrixConstRef Wp = W(all, ispecies);
    VectorConstRef np = n(ispecies);
    Vector res = Wp * np;
    return res;
}

auto ChemicalEngineArray::systemVolume(size_t ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->Vs;
    // return pimpl->node_array->pNodT1()[ndx]->cVs();
}

auto ChemicalEngineArray::systemMass(size_t ndx) const -> double
{
    //return pimpl->node_array->pNodT1()[ndx]->cMs();
    return pimpl->node_array->pNodT1()[ndx]->Ms;
}

auto ChemicalEngineArray::ionicStrength(size_t ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->IC;
}

auto ChemicalEngineArray::pH(size_t ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->pH;
}

auto ChemicalEngineArray::pe(size_t ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->pe;
}

auto ChemicalEngineArray::Eh(size_t ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->Eh;
}

auto ChemicalEngineArray::systemGibbsEnergy(size_t ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->Gs;   // *RT to get it in J
}

auto ChemicalEngineArray::phaseDensities(size_t ndx) const -> Vector
{
    const size_t P = numPhases();
    Vector phDensities(P);

    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(size_t i = 0; i < P; ++i)
    {
        if( pimpl->node_array->getCalcNode()->Ph_Volume(i) != 0.0)
            phDensities[i] = pimpl->node_array->getCalcNode()->Ph_Mass(i)
                             /  pimpl->node_array->getCalcNode()->Ph_Volume(i);
        else
            phDensities[i] = 0.0;
    }
    return phDensities;
}

auto ChemicalEngineArray::phaseMasses(size_t ndx) const -> Vector
{
    const size_t P = numPhases();
    Vector phMasses(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(size_t i = 0; i < P; ++i)
        phMasses[i] = pimpl->node_array->getCalcNode()->Ph_Mass(i);
    return phMasses;
}

auto ChemicalEngineArray::phaseAmounts(size_t ndx) const -> Vector
{
    const size_t P = numPhases();
    Vector phAmounts(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(size_t i = 0; i < P; ++i)
        phAmounts[i] = pimpl->node_array->getCalcNode()->Ph_Moles(i);
    return phAmounts;
}

auto ChemicalEngineArray::phaseVolumes(size_t ndx) const -> Vector
{
    const size_t P = numPhases();
    Vector phVolumes(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(size_t i = 0; i < P; ++i)
        phVolumes[i] = pimpl->node_array->getCalcNode()->Ph_Volume(i);
    return phVolumes;
}

auto ChemicalEngineArray::phaseEnthalpies(size_t ndx) const -> Vector
{
    const size_t P = numPhases();
    Vector phEnthalpies(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(size_t i = 0; i < P; ++i)
        phEnthalpies[i] = pimpl->node_array->getCalcNode()->Ph_Enthalpy(i);
    return phEnthalpies;
}

auto ChemicalEngineArray::phaseEntropies(size_t ndx) const -> Vector
{
    const size_t P = numPhases();
    Vector phEntropies(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(size_t i = 0; i < P; ++i)
        phEntropies[i] = pimpl->node_array->getCalcNode()->Ph_Entropy(i);
    return phEntropies;
}

auto ChemicalEngineArray::phaseHeatCapacitiesConstP(size_t ndx) const -> Vector
{
    const size_t P = numPhases();
    Vector phHeatCapacitiesConstP(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(size_t i = 0; i < P; ++i)
        phHeatCapacitiesConstP[i] = pimpl->node_array->getCalcNode()->Ph_HeatCapacityCp(i);
    return phHeatCapacitiesConstP;
}

auto ChemicalEngineArray::phaseDensity_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    if(pimpl->node_array->getCalcNode()->Ph_Volume(iphase) > 1e-15) {
        return pimpl->node_array->getCalcNode()->Ph_Mass(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Volume(iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseDensity(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseDensity_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseDensity(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseDensity_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMass_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode()->Ph_Mass(iphase);
}

auto ChemicalEngineArray::phaseMass(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseMass_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMass(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseMass_i(ndx, iphase);
    }
    return 0.0;
}
auto ChemicalEngineArray::phaseAmount_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode()->Ph_Moles(iphase);
}

auto ChemicalEngineArray::phaseAmount(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseAmount_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseAmount(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseAmount_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseVolume_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode()->Ph_Volume(iphase);
}

auto ChemicalEngineArray::phaseVolume(size_t ndx, size_t  iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseVolume_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseVolume(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseVolume_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseEnthalpy_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode()->Ph_Enthalpy(iphase);
}

auto ChemicalEngineArray::phaseEnthalpy(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseEnthalpy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseEnthalpy(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseEnthalpy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseEntropy_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode()->Ph_Entropy(iphase);
}

auto ChemicalEngineArray::phaseEntropy(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseEntropy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseEntropy(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseEntropy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseHeatCapacityConstP_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode()->Ph_HeatCapacityCp(iphase);
}

auto ChemicalEngineArray::phaseHeatCapacityConstP(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseHeatCapacityConstP_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseHeatCapacityConstP(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseHeatCapacityConstP_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSatIndices(size_t ndx) const -> Vector
{
    const size_t P = numPhases();
    Vector phSatIndices(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(size_t i = 0; i < P; ++i)
        phSatIndices[i] = pimpl->node_array->getCalcNode()->Ph_SatInd(i);
    return phSatIndices;
}

auto ChemicalEngineArray::phaseSatIndex_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode()->Ph_SatInd(iphase);
}

auto ChemicalEngineArray::phaseSatIndex(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseSatIndex_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSatIndex(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseSatIndex_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMolarGibbsEnergy_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarGibbsEnergy = 0.0;
    if( pimpl->node_array->getCalcNode()->Ph_Mole(iphase) > 1e-15 ) {
        phMolarGibbsEnergy = pimpl->node_array->getCalcNode()->Ph_GibbsEnergy(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mole(iphase);
    }
    return phMolarGibbsEnergy;
}

auto ChemicalEngineArray::phaseMolarGibbsEnergy(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseMolarGibbsEnergy_i(ndx, iphase);
    }
    return 0;
}

auto ChemicalEngineArray::phaseMolarGibbsEnergy(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseMolarGibbsEnergy_i(ndx, iphase);
    }
    return 0;
}
auto ChemicalEngineArray::phaseMolarEnthalpy_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarEnthalpy = 0.0;
    if( pimpl->node_array->getCalcNode()->Ph_Mole(iphase) > 1e-15 )  {
        phMolarEnthalpy = pimpl->node_array->getCalcNode()->Ph_Enthalpy(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mole(iphase);
    }
    return phMolarEnthalpy;      
}

auto ChemicalEngineArray::phaseMolarEnthalpy(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseMolarEnthalpy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMolarEnthalpy(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseMolarEnthalpy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMolarVolume_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarVolume = 0.0;
    if(pimpl->node_array->getCalcNode()->Ph_Mole(iphase) > 1e-15) {
        phMolarVolume = pimpl->node_array->getCalcNode()->Ph_Volume(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mole(iphase);
    }
    return phMolarVolume;
}

auto ChemicalEngineArray::phaseMolarVolume(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseMolarVolume_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMolarVolume(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseMolarVolume_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMolarEntropy_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarEntropy = 0.0;
    if(pimpl->node_array->getCalcNode()->Ph_Mole(iphase) > 1e-15) {
        phMolarEntropy = pimpl->node_array->getCalcNode()->Ph_Entropy(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mole(iphase);
    }
    return phMolarEntropy;
}

auto ChemicalEngineArray::phaseMolarEntropy(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseMolarEntropy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMolarEntropy(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseMolarEntropy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::moleFractions(size_t ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const size_t N = numSpecies();
    Vector amountSp = Vector::Map(pimpl->node_array->getCalcNode()->pCNode()->xDC, N);
    Vector molFractions(N);

    size_t counter = 0;
    for(size_t k = 0; k < numPhases(); ++k)
    {
        size_t numSpInPh = numSpeciesInPhase(k);
        double amountPh = pimpl->node_array->getCalcNode()->Ph_Moles(k);
        if(numSpInPh == 1 )
        {
            if(amountPh > 0.0)
                molFractions(counter) = 1.0;
            else
                molFractions(counter) = 0.0;
            counter++;
            continue;
        }
        for(size_t i = 0; i < numSpInPh; ++i)
        {
            if(amountPh != 0.0)
                molFractions(counter) = amountSp(counter)/amountPh;
            else
                molFractions(counter) = 0.0;
            counter++;
        }
    }
    return molFractions;
}

auto ChemicalEngineArray::lnActivityCoefficients(size_t ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const size_t N = numSpecies();
    Vector ln_activity_coefficients(N);
    for(size_t i = 0; i < N; ++i)
        ln_activity_coefficients[i] = std::log(pimpl->node_array->getCalcNode()->Get_gDC(i));
    return ln_activity_coefficients;
}

auto ChemicalEngineArray::lnActivities(size_t ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const size_t N = numSpecies();
    Vector ln_activities(N);
    for(size_t i = 0; i < N; ++i)
        ln_activities[i] = std::log(pimpl->node_array->getCalcNode()->Get_aDC(i, true));
    return ln_activities;
}

auto ChemicalEngineArray::lnConcentrations(size_t ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const size_t N = numSpecies();
    Vector ln_concentrations(N);
    for(size_t i = 0; i < N; ++i)
        ln_concentrations[i] = std::log(pimpl->node_array->getCalcNode()->Get_cDC(i));
    return ln_concentrations;
}

auto ChemicalEngineArray::chemicalPotentials(size_t ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const size_t N = numSpecies();
    Vector chemPotentials(N);
    for(size_t i = 0; i < numSpecies(); ++i)
        chemPotentials[i] = pimpl->node_array->getCalcNode()->Get_muDC(i, false);
    return chemPotentials;
}

auto ChemicalEngineArray::phaseMolarHeatCapacityConstP_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarCp = 0.0;
    if(pimpl->node_array->getCalcNode()->Ph_Mole(iphase) > 1e-15) {
        phMolarCp = pimpl->node_array->getCalcNode()->Ph_HeatCapacityCp(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mole(iphase);
    }
    return phMolarCp;
}

auto ChemicalEngineArray::phaseMolarHeatCapacityConstP(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseMolarHeatCapacityConstP_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseMolarHeatCapacityConstP(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseMolarHeatCapacityConstP_i(ndx, iphase);
    }
    return 0.0;
}


auto ChemicalEngineArray::phaseSpecificHeatCapacityConstP_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecCp = 0.0;
    if(pimpl->node_array->getCalcNode()->Ph_Mass(iphase) > 1e-15) {
        phSpecCp = pimpl->node_array->getCalcNode()->Ph_HeatCapacityCp(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mass(iphase);
    }
    return phSpecCp;
}

auto ChemicalEngineArray::phaseSpecificHeatCapacityConstP(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseSpecificHeatCapacityConstP_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificHeatCapacityConstP(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseSpecificHeatCapacityConstP_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificGibbsEnergy_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecGibbsEnergy = 0.0;
    if(pimpl->node_array->getCalcNode()->Ph_Mass(iphase) > 1e-15) {
        phSpecGibbsEnergy = pimpl->node_array->getCalcNode()->Ph_GibbsEnergy(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mass(iphase);
    }
    return phSpecGibbsEnergy;
}

auto ChemicalEngineArray::phaseSpecificGibbsEnergy(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseSpecificGibbsEnergy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificGibbsEnergy(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseSpecificGibbsEnergy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificEnthalpy_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecEnthalpy = 0.0;
    if(pimpl->node_array->getCalcNode()->Ph_Mass(iphase) > 1e-15) {
        phSpecEnthalpy = pimpl->node_array->getCalcNode()->Ph_Enthalpy(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mass(iphase);
    }
    return phSpecEnthalpy;
}

auto ChemicalEngineArray::phaseSpecificEnthalpy(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseSpecificEnthalpy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificEnthalpy(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseSpecificEnthalpy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificVolume_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecVolume = 0.0;
    if(pimpl->node_array->getCalcNode()->Ph_Mass(iphase) > 1e-15) {
        phSpecVolume = pimpl->node_array->getCalcNode()->Ph_Volume(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mass(iphase);
    }
    return phSpecVolume;
}

auto ChemicalEngineArray::phaseSpecificVolume(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseSpecificVolume_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificVolume(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseSpecificVolume_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificEntropy_i(size_t ndx, size_t iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecEntropy = 0.0;
    if(pimpl->node_array->getCalcNode()->Ph_Mass(iphase) > 1e-15) {
        phSpecEntropy = pimpl->node_array->getCalcNode()->Ph_Entropy(iphase)
        /pimpl->node_array->getCalcNode()->Ph_Mass(iphase);
    }
    return phSpecEntropy;
}

auto ChemicalEngineArray::phaseSpecificEntropy(size_t ndx, size_t iphase) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_Phase_xCH(iphase)) {
        return phaseSpecificEntropy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::phaseSpecificEntropy(size_t ndx, std::string phase) const -> double
{
    auto iphase = indexPhase(phase);
    if(iphase<numPhases()) {
        return phaseSpecificEntropy_i(ndx, iphase);
    }
    return 0.0;
}

auto ChemicalEngineArray::systemEnthalpy(size_t ndx) const -> double   // in J
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    // enth = pimpl->node->phEnthalpies.sum(); // needs to be computed first anyway
    auto enth = 0.0;
    for(size_t i = 0; i < numPhases(); ++i)
        enth += pimpl->node_array->getCalcNode()->Ph_Enthalpy(i);
    return enth;  // in J
    // enth /= pimpl->node->pCNode()->Ms; in  kg
    // return pimpl->node->pCNode()->Hs;  // Not computed in GEMS3K so far!
}

auto ChemicalEngineArray::systemEntropy(size_t ndx) const -> double
{
    return 0.0;
}

auto ChemicalEngineArray::systemHeatCapacityConstP(size_t ndx) const -> double
{
    return 0.0;
}

// Uses curent node TK and P
auto ChemicalEngineArray::standardMolarGibbsEnergy_i(double TK, double P, size_t ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode()->DC_xDB_to_xCH(ispecies);
    auto stMolarG = pimpl->node_array->getCalcNode()->DC_G0(idc, P, TK,  false);
    return stMolarG;
}


auto ChemicalEngineArray::standardMolarGibbsEnergy(double TK, double P, size_t ispecies) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        return standardMolarGibbsEnergy_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::standardMolarGibbsEnergy(double TK, double P, std::string species, std::optional<std::string> phase) const -> double
{
    auto ispecies = indexSpecies(species, phase);
    if(ispecies<numSpecies()) {
        return standardMolarGibbsEnergy_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::standardMolarEnthalpy_i(double TK, double P, size_t ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode()->DC_xDB_to_xCH(ispecies);
    auto stMolarH = pimpl->node_array->getCalcNode()->DC_H0(idc, P, TK);
    return stMolarH;
}

auto ChemicalEngineArray::standardMolarEnthalpy(double TK, double P, size_t ispecies) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        return standardMolarEnthalpy_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::standardMolarEnthalpy(double TK, double P, std::string species, std::optional<std::string> phase) const -> double
{
    auto ispecies = indexSpecies(species, phase);
    if(ispecies<numSpecies()) {
        return standardMolarEnthalpy_i(TK, P, ispecies);
    }
    return 0.;
}


auto ChemicalEngineArray::standardMolarVolume_i(double TK, double P, size_t ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode()->DC_xDB_to_xCH(ispecies);
    auto stMolarV = pimpl->node_array->getCalcNode()->DC_V0(idc, P, TK);
    return stMolarV;
}

auto ChemicalEngineArray::standardMolarVolume(double TK, double P, size_t ispecies) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        return standardMolarVolume_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::standardMolarVolume(double TK, double P, std::string species, std::optional<std::string> phase) const -> double
{
    auto ispecies = indexSpecies(species, phase);
    if(ispecies<numSpecies()) {
        return standardMolarVolume_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::standardMolarEntropy_i(double TK, double P, size_t ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode()->DC_xDB_to_xCH(ispecies);
    auto stMolarS = pimpl->node_array->getCalcNode()->DC_S0(idc, P, TK);
    return stMolarS;
}

auto ChemicalEngineArray::standardMolarEntropy(double TK, double P, size_t ispecies) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        return standardMolarEntropy_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::standardMolarEntropy(double TK, double P, std::string species, std::optional<std::string> phase) const -> double
{
    auto ispecies = indexSpecies(species, phase);
    if(ispecies<numSpecies()) {
        return standardMolarEntropy_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::standardMolarHeatCapacityConstP_i(double TK, double P, size_t ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode()->DC_xDB_to_xCH( ispecies );
    auto stMolarCp = pimpl->node_array->getCalcNode()->DC_Cp0(idc, P, TK);
    return stMolarCp;
}

auto ChemicalEngineArray::standardMolarHeatCapacityConstP(double TK, double P, size_t ispecies) const -> double
{
    if(pimpl->node_array->getCalcNode()->check_DC_xCH(ispecies)) {
        return standardMolarHeatCapacityConstP_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::standardMolarHeatCapacityConstP(double TK, double P, std::string species, std::optional<std::string> phase) const -> double
{
    auto ispecies = indexSpecies(species, phase);
    if(ispecies<numSpecies()) {
        return standardMolarHeatCapacityConstP_i(TK, P, ispecies);
    }
    return 0.;
}

auto ChemicalEngineArray::aqueousPhaseName() const -> std::string
{
    for(Index ii = 0; ii < pimpl->node_array->pCSD()->nPS; ++ii) {
        if( pimpl->node_array->pCSD()->ccPH[ii] == 'a' ) {
            return phaseName_i(ii);
        }
    }
    return {};
}

auto ChemicalEngineArray::gasPhaseName() const -> std::string
{
    for(Index ii = 0; ii < pimpl->node_array->pCSD()->nPS; ++ii) {
        if( pimpl->node_array->pCSD()->ccPH[ii] == 'g' ||
            pimpl->node_array->pCSD()->ccPH[ii] == 'p' ||
            pimpl->node_array->pCSD()->ccPH[ii] == 'f') {
            return  phaseName_i(ii);
        }
    }
    return {};
}

// ----------------------------------------------

// Aquatic systems only (assuming aqueous phase is the first one and H2O-solvent
// is the last species in it)
//auto ChemicalEngineArray::speciesMolalities() const -> VectorConstRef
//{
//    for(size_t i = 0; i < numSpecies(); ++i)
//        pimpl->molalities[i] = pimpl->node->Get_cDC(i);
//    if(numSpecies() > numSpeciesInPhase(0))
//    {
//        Vector amountSp = Vector::Map(pimpl->node->pCNode()->xDC, numSpecies());
//        size_t H2Oindex = numSpeciesInPhase(0)-1;
//        double H2Omass = pimpl->node->Get_nDC(H2Oindex)*pimpl->node->DCmm(H2Oindex); // in kg
//        if( H2Omass != 0.0 )
//            for(size_t i = H2Oindex+1; i < numSpecies(); ++i)
//                pimpl->molalities[i] = speciesAmounts()(i)/H2Omass;
//    }
//    return pimpl->molalities;
//}


// Implemented on Oct 1, 2020
// auto ChemicalEngineArray::standardMolarGibbsEnergies() const -> VectorConstRef
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     for(size_t i = 0; i < numSpecies(); ++i)
//     {
//         auto idc = pimpl->node->DC_xDB_to_xCH( i );
//         pimpl->stMolarGibbsEnergies[i] = pimpl->node->DC_G0(idc, P, TK,  false);
//     }    
//     return pimpl->stMolarGibbsEnergies;
// }

// TBD
// auto ChemicalEngineArray::standardMolarInternalEnergy(size_t ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarU = pimpl->node->DC_U0(idc, P, TK); 
//     return stMolarU;
// }

// TBD
// auto ChemicalEngineArray::standardMolarHelmholtzEnergy(size_t ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarA = pimpl->node->DC_A0(idc, P, TK); 
//     return stMolarA;
// }


// TBD
// auto ChemicalEngineArray::standardMolarHeatCapacityConstV(size_t ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarCv = pimpl->node->DC_Cp0(idc, P, TK) - R_CONSTANT; 
//     return stMolarCv;
// }


// TBD
// auto ChemicalEngineArray::phaseMolarInternalEnergy(size_t iphase) const -> double
// {
//     return 0.0;
// }

// TBD
// auto ChemicalEngineArray::phaseMolarHelmholtzEnergy(size_t iphase) const -> double
// {
//     return 0.0;
// }


// TBD
// auto ChemicalEngineArray::phaseMolarHeatCapacityConstV(size_t iphase) const -> double
// {
//     return 0.0;
// }


// TBD
// auto ChemicalEngineArray::phaseSpecificInternalEnergy(size_t iphase) const -> double
// {
//     return 0.0;
// }

// TBD
// auto ChemicalEngineArray::phaseSpecificHelmholtzEnergy(size_t iphase) const -> double
// {
//     return 0.0;
// }

// auto ChemicalEngineArray::phaseSpecificHeatCapacityConstV(size_t iphase) const -> double
// {
//     return 0.0;
// }

auto operator<<(std::ostream& out, const ChemicalEngineArray& state) -> std::ostream&
{
    auto flags = out.flags();

    for(size_t ii=0; ii<state.nNodes(); ii++)
    {
        const double T = state.temperature(ii);
        const double P = state.pressure(ii);
        VectorConstRef n = state.speciesAmounts(ii);
        const Vector activity_coeffs = state.lnActivityCoefficients(ii).array().exp();
        const Vector activities = state.lnActivities(ii).array() * 0.4343; // .exp();
        const Vector chemical_potentials = state.chemicalPotentials(ii).array() / 1000.;
        const Vector concentrations = state.lnConcentrations(ii).array().exp();
        //    const Vector molalities = state.speciesMolalities().array();
        const Vector molfractions = state.moleFractions(ii).array();
        //    const Vector satindices = state.phaseSatIndices().array();

        const size_t num_phases = state.numPhases();
        const size_t bar_size = std::max(size_t(9), num_phases + 2) * 25;
        const std::string bar1(bar_size, '=');
        const std::string bar2(bar_size, '-');

        out << bar1 << std::endl;
        out << std::left << std::setw(25) << "Temperature[K]";
        out << std::left << std::setw(25) << "Temperature[C]";
        out << std::left << std::setw(25) << "Pressure[MPa]";
        out << std::endl << bar2 << std::endl;

        out << std::left << std::setw(25) << T;
        out << std::left << std::setw(25) << T - 273.15;
        out << std::left << std::setw(25) << P * 1e-6;
        out << std::endl;

        // Set output in scientific notation
        out << std::setprecision(6);

        // Output the table of the element-related state
        out << bar1 << std::endl;
        out << std::left << std::setw(25) << "Element";
        out << std::left << std::setw(25) << "Amount[mol]";
        for(size_t i = 0; i < state.numPhases(); ++i)
            out << std::left << std::setw(25) << state.phaseName(i) + "[mol]";
        out << std::endl;
        out << bar2 << std::endl;
        for(size_t i = 0; i < state.numElements(); ++i)
        {
            out << std::left << std::setw(25) << state.elementName(i);
            out << std::left << std::setw(25) << state.elementAmounts(ii)[i];
            for(size_t j = 0; j < state.numPhases(); ++j)
                out << std::left << std::setw(25) << state.elementAmountsInPhase(ii, j)[i];
            out << std::endl;
        }
        out << bar2 << std::endl;
        out << std::left << std::setw(25) << "PhaseAmount[mol]";
        out << std::left << std::setw(25) << 0.0;
        for(size_t j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.phaseAmounts(ii)[j];
        out << std::endl;
        out << std::left << std::setw(25) << "PhaseMass[kg]";
        out << std::left << std::setw(25) << 0.0;
        for(size_t j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.phaseMasses(ii)[j];
        out << std::endl;
        out << std::left << std::setw(25) << "PhaseVolume[m^3]";
        out << std::left << std::setw(25) << 0.0;
        for(size_t j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.phaseVolumes(ii)[j];
        out << std::endl;
        out << std::left << std::setw(25) << "PhaseDensity[kg/m^3]";
        out << std::left << std::setw(25) << 0.0;
        for(size_t j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.phaseDensities(ii)[j];
        out << std::endl;
        out << std::left << std::setw(25) << "PhaseSatIndex[lg]";
        out << std::left << std::setw(25) << 0.0;
        for(size_t j = 0; j < state.numPhases(); ++j)
            //        out << std::left << std::setw(25) << satindices[j];
            out << std::left << std::setw(25) << state.phaseSatIndices(ii)[j];
        out << std::endl;
        out << std::endl;

        // Output the table of the species-related state
        out << bar1 << std::endl;
        out << std::left << std::setw(25) << "Species";
        out << std::left << std::setw(25) << "Amount[mol]";
        out << std::left << std::setw(25) << "Concentration[-]";
        out << std::left << std::setw(25) << "ActivityCoeff[-]";
        //   out << std::left << std::setw(25) << "Molality[mol/kgw]";
        out << std::left << std::setw(25) << "MoleFraction[-]";
        out << std::left << std::setw(25) << "ChemPotential[kJ/mol]";
        out << std::left << std::setw(25) << "Activity[lg]";
        out << std::endl;
        out << bar2 << std::endl;
        for(size_t i = 0; i < state.numSpecies(); ++i)
        {
            out << std::left << std::setw(25) << state.speciesName(i);
            out << std::left << std::setw(25) << n[i];
            out << std::left << std::setw(25) << concentrations[i];
            out << std::left << std::setw(25) << activity_coeffs[i];
            //       out << std::left << std::setw(25) << molalities[i];
            out << std::left << std::setw(25) << molfractions[i];
            out << std::left << std::setw(25) << chemical_potentials[i];
            out << std::left << std::setw(25) << activities[i];
            out << std::endl;
        }

        // Output the table of the aqueous phase related state
        out << bar1 << std::endl;
        out << std::left << std::setw(25) << "Mass[kg]";
        out << std::left << std::setw(25) << "Volume[m^3]";
        out << std::left << std::setw(25) << "G[norm]";
        out << std::left << std::setw(25) << "H[kJ]";
        out << std::left << std::setw(25) << "IonicStrength[molal]";
        out << std::left << std::setw(25) << "pH";
        out << std::left << std::setw(25) << "pE";
        out << std::left << std::setw(25) << "Eh[V]";
        out << std::endl << bar2 << std::endl;
        out << std::left << std::setw(25) << state.systemMass(ii);
        out << std::left << std::setw(25) << state.systemVolume(ii);
        out << std::left << std::setw(25) << state.systemGibbsEnergy(ii);
        out << std::left << std::setw(25) << state.systemEnthalpy(ii)/1000.0;
        out << std::left << std::setw(25) << state.ionicStrength(ii);
        out << std::left << std::setw(25) << state.pH(ii);
        out << std::left << std::setw(25) << state.pe(ii);
        out << std::left << std::setw(25) << state.Eh(ii);
        out << std::endl << bar1 << std::endl;
    }
    // Recover the previous state of `out`
    out.flags(flags);

    return out;
}

} // namespace xGEM2MT
