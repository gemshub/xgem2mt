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

#include "ChemicalEngineArray.hpp"
#include "MassTransport_p.h"

// C++ includes
#include <chrono>
#include <iomanip>
#include <memory>

// GEMS3K includes
#include <GEMS3K/nodearray.h>
#include <GEMS3K/gems3k_impex.h>

namespace xGEM2MT {

struct ChemicalEngineArray::Impl
{
    /// Construct a default Impl instance
    Impl(long int nodes): n_nodes(nodes)
    {}

    /// The options for ChemicalEngineArray
    ChemicalEngineArrayOptions options;

    /// Number of allocated nodes
    long int n_nodes;

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

ChemicalEngineArray::ChemicalEngineArray(long int nodes)
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
    long int res = pimpl->node_array->GEM_init(ipm_filename.c_str(),
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
    long int res = pimpl->node_array->GEM_init( dch_json, ipm_json, dbr_json, nodeTypes.empty() ? nullptr : nodeTypes.data() );

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

long ChemicalEngineArray::reequilibrate(long ndx, bool warmstart)
{
    long int  Mode = pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    TNode workNode(pimpl->node_array->getCalcNode());
    return pimpl->node_array->RunGEM( workNode, ndx, Mode, pimpl->node_array->pNodT1() );
}

void ChemicalEngineArray::reequilibrate(bool warmstart)
{
    long int  Mode = pimpl->options.warmstart ? NEED_GEM_SIA : NEED_GEM_AIA;
    vector<long int> nodeFlags, retCodes;
    for( long int ii=0; ii<pimpl->node_array->nNodes(); ++ii)
    {
        nodeFlags.push_back(1);
        retCodes.push_back(0);
    }
    pimpl->node_array->RunGEM( Mode, pimpl->node_array->nNodes(), pimpl->node_array->pNodT1(), nodeFlags.data(), retCodes.data() );
}

bool ChemicalEngineArray::equilibrate(long ndx)
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

    for (long int ii=0; ii<nnodes; ii++)
    {
        bool NeedCopy = always;
        if(!NeedCopy) {
            for(long int ic=0; ic < CH->nICb-1; ic++) // do not check charge
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
    const Index E = numElements();
    const Index N = numSpecies();
    pimpl->formula_matrix = Matrix::Zero(E, N);
    for(Index i = 0; i < N; ++i)
        for(Index j = 0; j < E; ++j)
            pimpl->formula_matrix(j, i) = pimpl->node_array->DCaJI(i, j);

    // Allocate memory for vector members
    pimpl->elemMolMasses.resize(E);
    pimpl->molMasses.resize(N);
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

long ChemicalEngineArray::nNodes() const
{
    return pimpl->n_nodes;
}

auto ChemicalEngineArray::readDbrFromFile(long int  ndx, std::string filename) -> void
{
    // Reads another dbr file with input system composition
    long int res = pimpl->node_array->GEMS3k_read_dbr(ndx, filename.c_str(), pimpl->io_mode);

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

auto ChemicalEngineArray::writeDbrToFile(long int  ndx, std::string filename) -> void
{
    // Write current node into a dbr file into path filename
    pimpl->node_array->GEMS3k_write_dbr( ndx, filename.c_str(), pimpl->io_mode);
}

auto ChemicalEngineArray::readDbrFromJsonString(long int  ndx, std::string dbr_json) -> void
{
    // Reads another dbr file with input system composition
    long int res = pimpl->node_array->GEMS3k_read_dbr( ndx, dbr_json, true );
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

auto ChemicalEngineArray::writeDbrToJsonString(long int  ndx) -> const std::string
{
    // Writes another node dbr with input system composition into json string dbr_json
    std::string dbr_json = "";
    long int res = pimpl->node_array->GEMS3k_write_dbr( ndx, dbr_json );

    // Check if there was a json write error
    if(res == -1)
        throw std::runtime_error("\n*** ERROR ***\n"
                                 "Could not process the provided DBR into JSON string.\n"
                                 "Check the code.");

    return dbr_json;
}

auto ChemicalEngineArray::numElements() const -> Index
{
    return pimpl->node_array->pCSD()->nIC;
}

auto ChemicalEngineArray::numSpecies() const -> Index
{
    return pimpl->node_array->pCSD()->nDC;
}

auto ChemicalEngineArray::numPhases() const -> Index
{
    return pimpl->node_array->pCSD()->nPH;
}

auto ChemicalEngineArray::numSpeciesInPhase(Index iphase) const -> Index
{
    return pimpl->node_array->pCSD()->nDCinPH[iphase];
}

auto ChemicalEngineArray::elementName(Index ielement) const -> std::string
{
    std::string name = pimpl->node_array->pCSD()->ICNL[ielement];
    if( name.length() > MaxICnameLength )
        name.resize( MaxICnameLength );
    return name;
}

auto ChemicalEngineArray::speciesName(Index ispecies) const -> std::string
{
    std::string name = pimpl->node_array->pCSD()->DCNL[ispecies];
    if( name.length() > MaxDCnameLength )
        name.resize( MaxDCnameLength );
    return name;
}

auto ChemicalEngineArray::phaseName(Index iphase) const -> std::string
{
    std::string name = pimpl->node_array->pCSD()->PHNL[iphase];
    if( name.length() > MaxPHnameLength )
        name.resize( MaxPHnameLength );
    return name;
}

auto ChemicalEngineArray::indexElement(std::string element) const -> Index
{
    const Index size = numElements();
    for(Index i = 0; i < size; ++i)
        if(elementName(i) == element)
            return i;
    return size; // in case that element name was not found
}

auto ChemicalEngineArray::indexSpecies(std::string species) const -> Index
{
    const Index size = numSpecies();
    for(Index i = 0; i < size; ++i)
        if(speciesName(i) == species)
            return i;
    return size; // in case that species name was not found
}

auto ChemicalEngineArray::indexSpeciesAll(std::string species) const -> VectorXi
{
    std::vector<Index> occur;
    occur.resize(0);
    const Index size = numSpecies();
    for(Index i = 0; i < size; ++i)
        if(speciesName(i) == species)
            occur.push_back(i);
    Index n_species_found = occur.size();
    VectorXi result;
    result.resize(0);
    if(n_species_found > 0)
    {
        result.resize(n_species_found);
        for(Index is=0; is<n_species_found; is++)
            result(is) = occur.at(is);
    }
    return result;
}

// Index of phase searched by name
auto ChemicalEngineArray::indexPhase(std::string phase) const -> Index
{
    const Index size = numPhases();
    for(Index i = 0; i < size; ++i)
        if(phaseName(i) == phase)
            return i;
    return size; // in case that phase name was not found
}

auto ChemicalEngineArray::indexPhaseAll(std::string phase) const -> VectorXi
{
    std::vector<Index> occur;
    occur.resize(0);
    const Index size = numPhases();
    for(Index i = 0; i < size; ++i)
        if(phaseName(i) == phase)
            occur.push_back(i);
    Index n_phases_found = occur.size();
    VectorXi result;
    result.resize(0);
    if(n_phases_found > 0)
    {
        result.resize(n_phases_found);
        for(Index is=0; is<n_phases_found; is++)
            result(is) = occur.at(is);
    }
    return result;
}

auto ChemicalEngineArray::indexPhaseWithSpecies(Index ispecies) const -> Index
{
    const Index size = numPhases();
    Index counter = 0;
    for(Index i = 0; i < size; ++i)
    {
        counter += numSpeciesInPhase(i);
        if(counter > ispecies) return i;
    }
    return size; // in case that species name was not found
}

auto ChemicalEngineArray::indexFirstSpeciesInPhase(Index iphase) const -> Index
{
    Index counter = 0;
    for(Index i = 0; i < iphase; ++i)
        counter += numSpeciesInPhase(i);
    return counter;
}

auto ChemicalEngineArray::formulaMatrix() const -> MatrixConstRef
{
    return pimpl->formula_matrix;
}

auto ChemicalEngineArray::elementMolarMasses() const -> VectorConstRef
{
    for(Index i = 0; i < numElements(); ++i)
        pimpl->elemMolMasses[i] = pimpl->node_array->pCSD()->ICmm[i];
    return pimpl->elemMolMasses;
}

auto ChemicalEngineArray::speciesMolarMasses() const -> VectorConstRef
{
    for(Index i = 0; i < numSpecies(); ++i)
        pimpl->molMasses[i] = pimpl->node_array->pCSD()->DCmm[i];
    return pimpl->molMasses;
}

auto ChemicalEngineArray::speciesCharge(Index ispecies) const -> double
{
    MatrixConstRef W = formulaMatrix();
    auto chargeRow = numElements();
    return W(chargeRow -1 , ispecies);
}

auto ChemicalEngineArray::setStandardMolarGibbsEnergy(std::string name, double TK, double P, double value) -> void
{
    auto ispecies = indexSpecies(name);
    auto idc = pimpl->node_array->getCalcNode().DC_xDB_to_xCH(ispecies);
    pimpl->node_array->Set_DC_G0(idc, P, TK,  value);
}

// Uses curent node TK and P
auto ChemicalEngineArray::standardMolarGibbsEnergy(double TK, double P, Index ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode().DC_xDB_to_xCH(ispecies);
    auto stMolarG = pimpl->node_array->getCalcNode().DC_G0(idc, P, TK,  false);
    return stMolarG;
}


auto ChemicalEngineArray::standardMolarEnthalpy(double TK, double P, Index ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode().DC_xDB_to_xCH(ispecies);
    auto stMolarH = pimpl->node_array->getCalcNode().DC_H0(idc, P, TK);
    return stMolarH;
}

auto ChemicalEngineArray::standardMolarVolume(double TK, double P, Index ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode().DC_xDB_to_xCH(ispecies);
    auto stMolarV = pimpl->node_array->getCalcNode().DC_V0(idc, P, TK);
    return stMolarV;
}

auto ChemicalEngineArray::standardMolarEntropy(double TK, double P, Index ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode().DC_xDB_to_xCH(ispecies);
    auto stMolarS = pimpl->node_array->getCalcNode().DC_S0(idc, P, TK);
    return stMolarS;
}

auto ChemicalEngineArray::standardMolarHeatCapacityConstP(double TK, double P, Index ispecies) const -> double
{
    auto idc = pimpl->node_array->getCalcNode().DC_xDB_to_xCH( ispecies );
    auto stMolarCp = pimpl->node_array->getCalcNode().DC_Cp0(idc, P, TK);
    return stMolarCp;
}

auto ChemicalEngineArray::setPT(long int ndx, double P, double T) const -> bool
{
    // Check if P, T are feasible
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    if(pimpl->node_array->getCalcNode().check_TP( T, P )==false)
        return true;
    // Set temperature and pressure
    pimpl->node_array->pNodT1()[ndx]->TK = T;
    pimpl->node_array->pNodT1()[ndx]->P = P;
    return false;
}

auto ChemicalEngineArray::setB(long int ndx, VectorConstRef b ) -> void
{
    // Set the mole amounts of elements
    for(Index i = 0; i < numElements(); ++i)
        pimpl->node_array->pNodT1()[ndx]->bIC[i] = b[i];
}

auto ChemicalEngineArray::speciesAmount(long int ndx, Index ispecies) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->xDC[ispecies];
}

auto ChemicalEngineArray::speciesAmount(long int ndx, std::string name) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->xDC[indexSpecies(name)];
}

auto ChemicalEngineArray::speciesAmounts(long int ndx) const -> Vector
{
    return Vector::Map(pimpl->node_array->pNodT1()[ndx]->xDC, numSpecies());
}

auto ChemicalEngineArray::setCorrectionbIC(long int ndx) -> void
{
    // Correction of the bIC vector
    long int ii, jj;
    for( ii=0; ii<numElements(); ii++ )
        pimpl->node_array->pNodT1()[ndx]->bIC[ii] = 0.0;
    MatrixConstRef W = formulaMatrix();
    auto SA = speciesAmounts(ndx);
    for(jj = 0; jj <numSpecies(); jj++)
        for( ii=0; ii<numElements(); ii++ )
            pimpl->node_array->pNodT1()[ndx]->bIC[ii] += SA[jj] * W(ii , jj);
}

auto ChemicalEngineArray::setSpeciesAmounts(long int ndx, VectorConstRef n) -> void
{
    // Updates all the species amount
    long int jj;
    for(jj = 0; jj <numSpecies(); jj++)
        pimpl->node_array->pNodT1()[ndx]->xDC[jj] = n[jj];
    setCorrectionbIC(ndx);
}

// Set the species of name (string) with the new amount
// Caution: this may be ambiguous as in GEMS3K, species with the same name
//   may occur in more than one condensed phase!
// An overload including phase name is needed!
auto ChemicalEngineArray::setSpeciesAmount(long int ndx, std::string name, double amount) -> void
{
    auto ispecies = indexSpecies(name);
    // Updates the species amount
    pimpl->node_array->pNodT1()[ndx]->xDC[ispecies] = amount;
    setCorrectionbIC(ndx);
}

auto ChemicalEngineArray::setSpeciesAmount(long int ndx, Index ispecies, double amount) -> void
{
    // Updates the species amount
    pimpl->node_array->pNodT1()[ndx]->xDC[ispecies] = amount;
    setCorrectionbIC(ndx);
}

auto ChemicalEngineArray::speciesUpperLimits(long int ndx) const -> Vector
{
    return Vector::Map(pimpl->node_array->pNodT1()[ndx]->dul, numSpecies());
}
auto ChemicalEngineArray::speciesLowerLimits(long int ndx) const -> Vector
{
    return Vector::Map(pimpl->node_array->pNodT1()[ndx]->dll, numSpecies());
}

// These methods may be ambiguous, as in GEMS3K, species with the same name
//   may occur in more than one condensed phase!
// An overload including phase name or index is needed!
//
auto ChemicalEngineArray::setSpeciesUpperLimit(long int ndx, std::string name, double amount) -> void
{
    auto ispecies = indexSpecies(name);
    setSpeciesUpperLimit(ndx, ispecies, amount);
}

auto ChemicalEngineArray::setSpeciesUpperLimit(long int ndx, Index ispecies, double amount) -> void
{
    double bound = (amount < 0.0? 1e6: amount);
    pimpl->node_array->pNodT1()[ndx]->dul[ispecies] = bound;
}

auto ChemicalEngineArray::setSpeciesUpperLimits(long int ndx, VectorConstRef n) -> void
{
    double bound;
    for(Index jj = 0; jj <numSpecies(); jj++) {
        bound = (n[jj] < 0.0? 1e6: n[jj]);
        pimpl->node_array->pNodT1()[ndx]->dul[jj] = bound;
    }
}

auto ChemicalEngineArray::setSpeciesLowerLimit(long int ndx, std::string name, double amount) -> void
{
    auto ispecies = indexSpecies(name);
    setSpeciesLowerLimit(ndx, ispecies, amount);
}

auto ChemicalEngineArray::setSpeciesLowerLimit(long int ndx, Index ispecies, double amount) -> void
{
    double bound = (amount < 0.0? 0.0: amount);
    pimpl->node_array->pNodT1()[ndx]->dll[ispecies] = bound;
}

auto ChemicalEngineArray::setSpeciesLowerLimits(long int ndx, VectorConstRef n) -> void
{
    double bound;
    for(Index jj = 0; jj <numSpecies(); jj++) {
        bound = (n[jj] < 0.0? 0.0: n[jj]);
        pimpl->node_array->pNodT1()[ndx]->dll[jj] = bound;;
    }
}

auto ChemicalEngineArray::converged(long int ndx) const -> bool
{
    const auto status = pimpl->node_array->pNodT1()[ndx]->NodeStatusCH;
    return status == OK_GEM_AIA || status == OK_GEM_SIA;
}

auto ChemicalEngineArray::numIterations(long int ndx) const -> Index
{
    return pimpl->node_array->pNodT1()[ndx]->IterDone;
}

auto ChemicalEngineArray::temperature(long int ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->TK;
}

auto ChemicalEngineArray::pressure(long int ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->P;
}

auto ChemicalEngineArray::elementAmounts(long int ndx) const -> Vector
{
    return Vector::Map(pimpl->node_array->pNodT1()[ndx]->bIC, numElements());
}

auto ChemicalEngineArray::elementAmountsInPhase(long int ndx, Index iphase) const -> Vector
{
    MatrixConstRef W = formulaMatrix();
    VectorConstRef n = speciesAmounts(ndx);
    const Index first = indexFirstSpeciesInPhase(iphase);
    const Index size = numSpeciesInPhase(iphase);
    MatrixConstRef Wp = W.middleCols(first, size);
    VectorConstRef np = n.segment(first, size);
    Vector res = Wp * np;
    return res;
}

auto ChemicalEngineArray::elementAmountsInSpecies(long int ndx, VectorXiConstRef ispecies) const -> Vector
{
    MatrixConstRef W = formulaMatrix();
    VectorConstRef n = speciesAmounts(ndx);
    MatrixConstRef Wp = W(all, ispecies);
    VectorConstRef np = n(ispecies);
    Vector res = Wp * np;
    return res;
}

auto ChemicalEngineArray::systemVolume(long int ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->Vs;
    // return pimpl->node_array->pNodT1()[ndx]->cVs();
}

auto ChemicalEngineArray::systemMass(long int ndx) const -> double
{
    //return pimpl->node_array->pNodT1()[ndx]->cMs();
    return pimpl->node_array->pNodT1()[ndx]->Ms;
}

auto ChemicalEngineArray::ionicStrength(long int ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->IC;
}

auto ChemicalEngineArray::pH(long int ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->pH;
}

auto ChemicalEngineArray::pe(long int ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->pe;
}

auto ChemicalEngineArray::Eh(long int ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->Eh;
}

auto ChemicalEngineArray::systemGibbsEnergy(long int ndx) const -> double
{
    return pimpl->node_array->pNodT1()[ndx]->Gs;   // *RT to get it in J
}

auto ChemicalEngineArray::phaseDensities(long int ndx) const -> Vector
{
    const Index P = numPhases();
    Vector phDensities(P);

    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(Index i = 0; i < P; ++i)
    {
        if( pimpl->node_array->getCalcNode().Ph_Volume(i) != 0.0)
            phDensities[i] = pimpl->node_array->getCalcNode().Ph_Mass(i)
                    /  pimpl->node_array->getCalcNode().Ph_Volume(i);
        else
            phDensities[i] = 0.0;
    }
    return phDensities;
}

auto ChemicalEngineArray::phaseMasses(long int ndx) const -> Vector
{
    const Index P = numPhases();
    Vector phMasses(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(Index i = 0; i < P; ++i)
        phMasses[i] = pimpl->node_array->getCalcNode().Ph_Mass(i);
    return phMasses;
}

auto ChemicalEngineArray::phaseAmounts(long int ndx) const -> Vector
{
    const Index P = numPhases();
    Vector phAmounts(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(Index i = 0; i < P; ++i)
        phAmounts[i] = pimpl->node_array->getCalcNode().Ph_Moles(i);
    return phAmounts;
}

auto ChemicalEngineArray::phaseVolumes(long int ndx) const -> Vector
{
    const Index P = numPhases();
    Vector phVolumes(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(Index i = 0; i < P; ++i)
        phVolumes[i] = pimpl->node_array->getCalcNode().Ph_Volume(i);
    return phVolumes;
}

auto ChemicalEngineArray::phaseEnthalpies(long int ndx) const -> Vector
{
    const Index P = numPhases();
    Vector phEnthalpies(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(Index i = 0; i < P; ++i)
        phEnthalpies[i] = pimpl->node_array->getCalcNode().Ph_Enthalpy(i);
    return phEnthalpies;
}

auto ChemicalEngineArray::phaseEntropies(long int ndx) const -> Vector
{
    const Index P = numPhases();
    Vector phEntropies(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(Index i = 0; i < P; ++i)
        phEntropies[i] = pimpl->node_array->getCalcNode().Ph_Entropy(i);
    return phEntropies;
}

auto ChemicalEngineArray::phaseHeatCapacitiesConstP(long int ndx) const -> Vector
{
    const Index P = numPhases();
    Vector phHeatCapacitiesConstP(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(Index i = 0; i < P; ++i)
        phHeatCapacitiesConstP[i] = pimpl->node_array->getCalcNode().Ph_HeatCapacityCp(i);
    return phHeatCapacitiesConstP;
}

auto ChemicalEngineArray::phaseDensity(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    if(pimpl->node_array->getCalcNode().Ph_Volume(iphase) > 1e-15)
        return pimpl->node_array->getCalcNode().Ph_Mass(iphase)
                /pimpl->node_array->getCalcNode().Ph_Volume(iphase);
    else return 0.0;
}

auto ChemicalEngineArray::phaseMass(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode().Ph_Mass(iphase);
}

auto ChemicalEngineArray::phaseAmount(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode().Ph_Moles(iphase);
}

auto ChemicalEngineArray::phaseVolume(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode().Ph_Volume(iphase);
}

auto ChemicalEngineArray::phaseEnthalpy(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode().Ph_Enthalpy(iphase);
}

auto ChemicalEngineArray::phaseEntropy(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode().Ph_Entropy(iphase);
}

auto ChemicalEngineArray::phaseHeatCapacityConstP(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode().Ph_HeatCapacityCp(iphase);
}

auto ChemicalEngineArray::phaseSatIndices(long int ndx) const -> Vector
{
    const Index P = numPhases();
    Vector phSatIndices(P);
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    for(Index i = 0; i < P; ++i)
        phSatIndices[i] = pimpl->node_array->getCalcNode().Ph_SatInd(i);
    return phSatIndices;
}

auto ChemicalEngineArray::phaseSatIndex(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    return pimpl->node_array->getCalcNode().Ph_SatInd(iphase);
}

auto ChemicalEngineArray::phaseMolarGibbsEnergy(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarGibbsEnergy = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mole(iphase) > 1e-15 )
        phMolarGibbsEnergy = pimpl->node_array->getCalcNode().Ph_GibbsEnergy(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mole(iphase);
    return phMolarGibbsEnergy;
}

auto ChemicalEngineArray::phaseMolarEnthalpy(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarEnthalpy = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mole(iphase) > 1e-15 )
        phMolarEnthalpy = pimpl->node_array->getCalcNode().Ph_Enthalpy(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mole(iphase);
    return phMolarEnthalpy;
}

auto ChemicalEngineArray::phaseMolarVolume(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarVolume = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mole(iphase) > 1e-15 )
        phMolarVolume = pimpl->node_array->getCalcNode().Ph_Volume(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mole(iphase);
    return phMolarVolume;
}

auto ChemicalEngineArray::phaseMolarEntropy(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarEntropy = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mole(iphase) > 1e-15 )
        phMolarEntropy = pimpl->node_array->getCalcNode().Ph_Entropy(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mole(iphase);
    return phMolarEntropy;
}

auto ChemicalEngineArray::moleFractions(long int ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const Index N = numSpecies();
    Vector amountSp = Vector::Map(pimpl->node_array->getCalcNode().pCNode()->xDC, N);
    Vector molFractions(N);

    Index counter = 0;
    for(Index k = 0; k < numPhases(); ++k)
    {
        Index numSpInPh = numSpeciesInPhase(k);
        double amountPh = pimpl->node_array->getCalcNode().Ph_Moles(k);
        if(numSpInPh == 1 )
        {
            if(amountPh > 0.0)
                molFractions(counter) = 1.0;
            else
                molFractions(counter) = 0.0;
            counter++;
            continue;
        }
        for(Index i = 0; i < numSpInPh; ++i)
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

auto ChemicalEngineArray::lnActivityCoefficients(long int ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const Index N = numSpecies();
    Vector ln_activity_coefficients(N);
    for(Index i = 0; i < N; ++i)
        ln_activity_coefficients[i] = std::log(pimpl->node_array->getCalcNode().Get_gDC(i));
    return ln_activity_coefficients;
}

auto ChemicalEngineArray::lnActivities(long int ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const Index N = numSpecies();
    Vector ln_activities(N);
    for(Index i = 0; i < N; ++i)
        ln_activities[i] = std::log(pimpl->node_array->getCalcNode().Get_aDC(i, true));
    return ln_activities;
}

auto ChemicalEngineArray::lnConcentrations(long int ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const Index N = numSpecies();
    Vector ln_concentrations(N);
    for(Index i = 0; i < N; ++i)
        ln_concentrations[i] = std::log(pimpl->node_array->getCalcNode().Get_cDC(i));
    return ln_concentrations;
}

auto ChemicalEngineArray::chemicalPotentials(long int ndx) const -> Vector
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    const Index N = numSpecies();
    Vector chemPotentials(N);
    for(Index i = 0; i < numSpecies(); ++i)
        chemPotentials[i] = pimpl->node_array->getCalcNode().Get_muDC(i, false);
    return chemPotentials;
}


auto ChemicalEngineArray::phaseMolarHeatCapacityConstP(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phMolarCp = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mole(iphase) > 1e-15 )
        phMolarCp = pimpl->node_array->getCalcNode().Ph_HeatCapacityCp(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mole(iphase);
    return phMolarCp;
}

auto ChemicalEngineArray::phaseSpecificHeatCapacityConstP(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecCp = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mass(iphase) > 1e-15 )
        phSpecCp = pimpl->node_array->getCalcNode().Ph_HeatCapacityCp(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mass(iphase);
    return phSpecCp;
}

auto ChemicalEngineArray::phaseSpecificGibbsEnergy(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecGibbsEnergy = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mass(iphase) > 1e-15 )
        phSpecGibbsEnergy = pimpl->node_array->getCalcNode().Ph_GibbsEnergy(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mass(iphase);
    return phSpecGibbsEnergy;
}

auto ChemicalEngineArray::phaseSpecificEnthalpy(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecEnthalpy = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mass(iphase) > 1e-15 )
        phSpecEnthalpy = pimpl->node_array->getCalcNode().Ph_Enthalpy(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mass(iphase);
    return phSpecEnthalpy;
}

auto ChemicalEngineArray::phaseSpecificVolume(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecVolume = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mass(iphase) > 1e-15 )
        phSpecVolume = pimpl->node_array->getCalcNode().Ph_Volume(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mass(iphase);
    return phSpecVolume;
}

auto ChemicalEngineArray::phaseSpecificEntropy(long int ndx, Index iphase) const -> double
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    double phSpecEntropy = 0.0;
    if( pimpl->node_array->getCalcNode().Ph_Mass(iphase) > 1e-15 )
        phSpecEntropy = pimpl->node_array->getCalcNode().Ph_Entropy(iphase)
                /pimpl->node_array->getCalcNode().Ph_Mass(iphase);
    return phSpecEntropy;
}


auto ChemicalEngineArray::systemEnthalpy(long int ndx) const -> double   // in J
{
    pimpl->node_array->LinkToNode(ndx, pimpl->n_nodes, pimpl->node_array->pNodT1());
    // enth = pimpl->node->phEnthalpies.sum(); // needs to be computed first anyway
    auto enth = 0.0;
    for(Index i = 0; i < numPhases(); ++i)
        enth += pimpl->node_array->getCalcNode().Ph_Enthalpy(i);
    return enth;  // in J
    // enth /= pimpl->node->pCNode()->Ms; in  kg
    // return pimpl->node->pCNode()->Hs;  // Not computed in GEMS3K so far!
}

auto ChemicalEngineArray::systemEntropy(long int ndx) const -> double
{
    return 0.0;
}

auto ChemicalEngineArray::systemHeatCapacityConstP(long int ndx) const -> double
{
    return 0.0;
}

// ----------------------------------------------

// Aquatic systems only (assuming aqueous phase is the first one and H2O-solvent
// is the last species in it)
//auto ChemicalEngineArray::speciesMolalities() const -> VectorConstRef
//{
//    for(Index i = 0; i < numSpecies(); ++i)
//        pimpl->molalities[i] = pimpl->node->Get_cDC(i);
//    if(numSpecies() > numSpeciesInPhase(0))
//    {
//        Vector amountSp = Vector::Map(pimpl->node->pCNode()->xDC, numSpecies());
//        Index H2Oindex = numSpeciesInPhase(0)-1;
//        double H2Omass = pimpl->node->Get_nDC(H2Oindex)*pimpl->node->DCmm(H2Oindex); // in kg
//        if( H2Omass != 0.0 )
//            for(Index i = H2Oindex+1; i < numSpecies(); ++i)
//                pimpl->molalities[i] = speciesAmounts()(i)/H2Omass;
//    }
//    return pimpl->molalities;
//}


// Implemented on Oct 1, 2020
// auto ChemicalEngineArray::standardMolarGibbsEnergies() const -> VectorConstRef
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     for(Index i = 0; i < numSpecies(); ++i)
//     {
//         auto idc = pimpl->node->DC_xDB_to_xCH( i );
//         pimpl->stMolarGibbsEnergies[i] = pimpl->node->DC_G0(idc, P, TK,  false);
//     }    
//     return pimpl->stMolarGibbsEnergies;
// }

// TBD
// auto ChemicalEngineArray::standardMolarInternalEnergy(Index ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarU = pimpl->node->DC_U0(idc, P, TK); 
//     return stMolarU;
// }

// TBD
// auto ChemicalEngineArray::standardMolarHelmholtzEnergy(Index ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarA = pimpl->node->DC_A0(idc, P, TK); 
//     return stMolarA;
// }


// TBD
// auto ChemicalEngineArray::standardMolarHeatCapacityConstV(Index ispecies) const -> double
// {
//     double TK = pimpl->node->cTK();
//     double P = pimpl->node->cP();
//     auto idc = pimpl->node->DC_xDB_to_xCH( ispecies );
//     auto stMolarCv = pimpl->node->DC_Cp0(idc, P, TK) - R_CONSTANT; 
//     return stMolarCv;
// }


// TBD
// auto ChemicalEngineArray::phaseMolarInternalEnergy(Index iphase) const -> double
// {
//     return 0.0;
// }

// TBD
// auto ChemicalEngineArray::phaseMolarHelmholtzEnergy(Index iphase) const -> double
// {
//     return 0.0;
// }


// TBD
// auto ChemicalEngineArray::phaseMolarHeatCapacityConstV(Index iphase) const -> double
// {
//     return 0.0;
// }


// TBD
// auto ChemicalEngineArray::phaseSpecificInternalEnergy(Index iphase) const -> double
// {
//     return 0.0;
// }

// TBD
// auto ChemicalEngineArray::phaseSpecificHelmholtzEnergy(Index iphase) const -> double
// {
//     return 0.0;
// }

// auto ChemicalEngineArray::phaseSpecificHeatCapacityConstV(Index iphase) const -> double
// {
//     return 0.0;
// }

auto operator<<(std::ostream& out, const ChemicalEngineArray& state) -> std::ostream&
{
    auto flags = out.flags();

    for(long int ii=0; ii<state.nNodes(); ii++)
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

        const Index num_phases = state.numPhases();
        const Index bar_size = std::max(Index(9), num_phases + 2) * 25;
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
        for(Index i = 0; i < state.numPhases(); ++i)
            out << std::left << std::setw(25) << state.phaseName(i) + "[mol]";
        out << std::endl;
        out << bar2 << std::endl;
        for(Index i = 0; i < state.numElements(); ++i)
        {
            out << std::left << std::setw(25) << state.elementName(i);
            out << std::left << std::setw(25) << state.elementAmounts(ii)[i];
            for(Index j = 0; j < state.numPhases(); ++j)
                out << std::left << std::setw(25) << state.elementAmountsInPhase(ii, j)[i];
            out << std::endl;
        }
        out << bar2 << std::endl;
        out << std::left << std::setw(25) << "PhaseAmount[mol]";
        out << std::left << std::setw(25) << 0.0;
        for(Index j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.phaseAmounts(ii)[j];
        out << std::endl;
        out << std::left << std::setw(25) << "PhaseMass[kg]";
        out << std::left << std::setw(25) << 0.0;
        for(Index j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.phaseMasses(ii)[j];
        out << std::endl;
        out << std::left << std::setw(25) << "PhaseVolume[m^3]";
        out << std::left << std::setw(25) << 0.0;
        for(Index j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.phaseVolumes(ii)[j];
        out << std::endl;
        out << std::left << std::setw(25) << "PhaseDensity[kg/m^3]";
        out << std::left << std::setw(25) << 0.0;
        for(Index j = 0; j < state.numPhases(); ++j)
            out << std::left << std::setw(25) << state.phaseDensities(ii)[j];
        out << std::endl;
        out << std::left << std::setw(25) << "PhaseSatIndex[lg]";
        out << std::left << std::setw(25) << 0.0;
        for(Index j = 0; j < state.numPhases(); ++j)
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
        for(Index i = 0; i < state.numSpecies(); ++i)
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
