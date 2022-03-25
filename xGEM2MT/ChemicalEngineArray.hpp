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
#include <xGEM2MT/Index.hpp>
#include <xGEM2MT/Eigen.hpp>

namespace xGEM2MT {

/// A type that describes the options for ChemicalEngineArray
struct ChemicalEngineArrayOptions
{
    /// The flag that indicates if smart start initial approximation is used
    bool warmstart = false;
    char useSIA='*';      ///<   Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    long int step=0;    ///<   extern actual time iterator
    double cdv=1e-9;       ///<   cutoff for IC amount differences in the node between time steps (mol, 1e-9)
    double cez=1e-12;       ///<   cutoff for minimal amounts of IC in node bulk compositions (mol, 1e-12)

    std::string to_json_string(bool dense = false) const {return "";}; //TBD
};

/// A wrapper class for ChemicalEngineArray code
class ChemicalEngineArray
{

public:

    /// Construct a ChemicalEngineArray object for the 1D arrangement of nodes
    ChemicalEngineArray(long int n_nodes);

    /// Construct a copy of a ChemicalEngineArray object.
    ChemicalEngineArray(const ChemicalEngineArray& other) = delete;

    /// Assign another ChemicalEngineArray object to this.
    auto operator=(ChemicalEngineArray other) -> ChemicalEngineArray& = delete;

    /// Destroy this ChemicalEngineArray instance.
    virtual ~ChemicalEngineArray();

    /// Set the options of the ChemicalEngineArray instance
    auto setOptions(const ChemicalEngineArrayOptions& options) -> void;
    /// Set smart start initial approximation (faster convergence, can be less accurate)
    auto setWarmStart() -> void;
    /// Set cold start automatic initial approximation (slower convergence, more accurate)
    auto setColdStart() -> void;

    /// Get total number of nodes in the node array
    auto nNodes() const -> long int;

    /// Initialize the ChemicalEngineArray object from a GEM-Selektor project file.
    /// @param filename The name of the file containing the definition of the chemical system
    auto initialize(std::string ipm_filename, std::string dbr_filename, long int* nodeTypes) -> void;

    /// Initialize the ChemicalEngineArray object from JSON strings for dch, ipm and dbr exported from GEM-Selektor \ .
    /// Normally to be used after ChemicalEngineArray();
    /// @param dch_json The json string containing the definition of the chemical system
    /// @param ipm_json The json string containing the parameters and settings for GEMS3K IPM-3 algorithm
    /// @param dbr_json The vector of json string containing the input node composition of the chemical system
    /// @param nodeTypes The initial node contents from DATABR files will be distributed among
    ///                   nodes in array according to the distribution index list nodeTypes
    auto initializeFromJsonStrings(std::string dch_json, std::string ipm_json,
                                   std::vector<std::string> dbr_json, long int* nodeTypes) -> void;

    /// Recalculate the equilibrium state with preset P,T,B and other inputs for node ndx
    /// If @param warmstart is true then GEMS3K uses contents of node as initial
    ///   guess, otherwise (false) it gets the simplex LP initial guess.
    /// Returns the return code of GEMS3K GEM_run call
    auto reequilibrate(long int  ndx, bool warmstart) -> long int;

    /// Calls GEM IPM calculation for for all nodes of TNodeArray
    /// in a loop with an optional openmp parallelization
    /// If @param warmstart is true then GEMS3K uses contents of node as initial
    ///   guess, otherwise (false) it gets the simplex LP initial guess.
    auto reequilibrate(bool warmstart) -> void;

    /// Recalculate the equilibrium state with preset P,T,B and other inputs for node ndx
    /// Here we do a GEM calculation in box ndx
    /// Here we compare this node for current time and for previous time
    /// to test need to recalculate equilibrium in this node because its vector b has changed
    /// Returns true if GEM Run finished with OK_GEM_AIA or OK_GEM_SIA
    auto equilibrate(long int  ndx) -> bool;

    /// Recalculate the equilibrium state with preset P,T,B and other inputs for all nodes
    /// Here we compare this node for current time and for previous time
    /// to test need to recalculate equilibrium in this node because its vector b has changed
    /// Returns true if all GEM Run finished with OK_GEM_AIA or OK_GEM_SIA
    auto equilibrate() -> bool;

    /// Dump MULTI, DATACH and DATABR files structure.
    /// @param filepath - IPM work structure file path&name
    /// @param put_all -  if true print array of nodes for previous time point
    ///                           and for current time point
    auto dumpGEMS3K( const std::string& filepath, bool put_all = false ) const -> void;

    /// Copies array of nodes for current time point to the same place in the
    /// array of nodes for previous time point.
    auto nextTimePoint(bool always) -> void;

    /// Reads another DBR file (with input system composition, T,P etc.) to a node with absolute index ndx.
    /// The DBR file must be compatible with the currently loaded IPM and DCH files
    /// (see description of initialize() function call).
    /// @param Null-terminated (C) string containing a full path to the input DBR disk file.
    auto readDbrFromFile(long int ndx, std::string filename) -> void;

    /// Reads another DBR object (with input system composition, T,P etc.) from JSON string
    /// to a node with absolute index ndx. The DBR file must be compatible with the currently
    /// loaded IPM and DCH objects (see description of initialize() function call).
    /// @param Null-terminated (C) string containing a full path to the input DBR disk file.
    auto readDbrFromJsonString(long int  ndx, std::string dbr_json) -> void;

    /// Returns a DBR object (normally after some changes via API and GEM calculation) as JSON string
    /// from a node with absolute index ndx
    /// In case of error raises an exception or returns empty string.
    auto writeDbrToJsonString(long int  ndx) -> const std::string;

    /// Writes a DBR file (normally after some changes via API and GEM calculation) from a node
    /// with absolute index ndx.
    /// @param Null-terminated (C) string containing a full path to the output DBR disk file.
    auto writeDbrToFile(long int ndx, std::string filename) -> void;

    /// Return the number of elements.
    auto numElements() const -> Index;

    /// Return the number of species.
    auto numSpecies() const -> Index;

    /// Return the number of phases.
    auto numPhases() const -> Index;

    /// Return the number of species in a phase or 0, if the phase was not found
    /// @param iphase The index of the phase.
    auto numSpeciesInPhase(Index iphase) const -> Index;

    /// Return the name of an element.
    /// @param ielement The index of the element.
    /// In case of error raises an exception or returns empty string.
    auto elementName(Index ielement) const -> std::string;

    /// Return the name of a species.
    /// @param ispecies The index of the species.
    /// In case of error raises an exception or returns empty string.
    auto speciesName(Index ispecies) const -> std::string;

    /// Return the name of a phase
    /// @param iphase The index of the phase.
    /// In case of error returns null string.
    auto phaseName(Index iphase) const -> std::string;

    /// Return the index of an @param element name, or number of elements, if not found.
    /// @param Index The index of the element.
    auto indexElement(std::string element) const -> Index;

    /// Return the index of a @param species name, or number of species, if not found.
    /// @param Index The index of the species.
    auto indexSpecies(std::string species) const -> Index;

    /// Return the vector of all indices of a @param species name, number of species found is VectorXi.size().
    /// @param VectorXi contains all indices of the species found in the species name list in this system.
    auto indexSpeciesAll(std::string species) const -> VectorXi;

    /// Return the index of a @param phase name, or number of phases, if not found.
    /// @param Index The index of the phase.
    auto indexPhase(std::string phase) const -> Index;

    /// Return the vector of all indices of a @param phase name, number of indices found is VectorXi.size().
    /// @param VectorXi contains all indices of the phases found in the phase name list in this system.
    auto indexPhaseAll(std::string phase) const -> VectorXi;

    /// Return the index of the phase containing a @param species name (ambiguous).
    /// @param Index The index of the species.
    auto indexPhaseWithSpecies(Index ispecies) const -> Index;

    /// Return the @param Index of the first species in a phase.
    /// @param iphase The index of the phase.
    auto indexFirstSpeciesInPhase(Index iphase) const -> Index;

    /// Return charge of a species.
    auto speciesCharge(Index ispecies) const -> double;

    /// Return the formula matrix of elements in species (in moles per mole species)
    /// (rows: elements, cols: species).
    auto formulaMatrix() const -> MatrixConstRef;

    /// Return the molar masses of the elements (in units of kg/mol)
    auto elementMolarMasses() const -> VectorConstRef;

    /// Return the molar masses of the species (in units of kg/mol)
    auto speciesMolarMasses() const -> VectorConstRef;

    /// Set the species @param name StandardMolarGibbsEnergy @param value
    auto setStandardMolarGibbsEnergy(std::string name, double TK, double P, double value) -> void;

    /// Return the standard  molar Gibbs energy of the species (in units of J/mol).
    auto standardMolarGibbsEnergy(double TK, double P, Index ispecies) const -> double;

    /// Return the standard  molar enthalpy of the species (in units of J/mol).
    auto standardMolarEnthalpy(double TK, double P, Index ispecies) const -> double;

    /// Return the standard  molar volume of the species (in units of m3/mol).
    auto standardMolarVolume(double TK, double P, Index ispecies) const -> double;

    /// Return the standard  molar entropies of the species (in units of J/(mol*K)).
    auto standardMolarEntropy(double TK, double P, Index ispecies) const -> double;

    /// Return the standard  molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardMolarHeatCapacityConstP(double TK, double P, Index ispecies) const -> double;

    /// Sets pressure to @param P (in Pa) and temperature to @param T (in K)
    /// in ndx instance without computing the equilibrium state
    /// Returns false if P.T was reset o.k., or true if P or T is out of
    /// range for the lookup for thermodynamic data in the node (with old P,T retained).
    auto setPT(long int ndx, double P, double T) const -> bool;

    /// Sets the amounts of elements (vector b) in moles into ndx
    /// instance without computing the equilibrium state
    auto setB(long int ndx, VectorConstRef b) -> void;

    /// Sets the whole speciation vector xDC to amounts provided in n and resets bIC bector accordingly
    ///  (often used for input initial guess of GEM solution in reactive transport sims)
    auto setSpeciesAmounts(long int ndx, VectorConstRef n) -> void;
    /// Set the species name to the new amount in xDC vector and adds this to bIC vector
    auto setSpeciesAmount(long int ndx, std::string name, double amount) -> void;
    /// Set the species with index ispecies to the new amount in xDC vector and adds this to bIC vector
    auto setSpeciesAmount(long int ndx, Index ispecies, double amount) -> void;

    // Caution: this may be ambiguous as in GEMS3K, species with the same name
    //   may occur in more than one condensed phase!
    // An overload including phase name is needed!
    /// Set the species @param name @param amount upper bound in dul vector
    /// if amount < 0 then resets to default 1e6
    auto setSpeciesUpperLimit(long int ndx, std::string name, double amount) -> void;
    /// Set species with index @param ispecies upper bound to @param amount in dul vector
    /// if @param amount < 0 then resets to default 1e6
    auto setSpeciesUpperLimit(long int ndx, Index ispecies, double amount) -> void;
    /// set all the upper limits for the species dul (in units of mol)
    auto setSpeciesUpperLimits(long int ndx, VectorConstRef n) -> void;

    /// Set the species @param name @param amount lower bound in dll vector
    /// if @param amount < 0 then resets to default 0
    auto setSpeciesLowerLimit(long int ndx, std::string name, double amount) -> void;
    /// Set species with index @param ispecies upper bound to @param amount in dul vector
    /// if @param amount < 0 then resets to default 0
    auto setSpeciesLowerLimit(long int ndx, Index ispecies, double amount) -> void;
    /// set all the lower limits for the species dll (in units of mol)
    auto setSpeciesLowerLimits(long int ndx, VectorConstRef n) -> void;

    /// Return the convergence result of the equilibrium calculation
    auto converged(long int ndx) const -> bool;

    /// Return the number of iterations of the equilibrium calculation
    auto numIterations(long int ndx) const -> Index;

    /// Return the temperature (in units of K)
    auto temperature(long int ndx) const -> double;

    /// Return the pressure (in units of Pa)
    auto pressure(long int ndx) const -> double;

    /// Return the amounts of the elements (in units of mol)
    auto elementAmounts(long int ndx) const -> Vector;

    /// Return the amounts of the elements in a phase.
    /// @param iphase The index of the phase.
    auto elementAmountsInPhase(long int ndx, Index iphase) const -> Vector;

    /// Return the amounts of the elements in a group of species.
    /// @param ispecies The vector of indices of the species.
    auto elementAmountsInSpecies(long int ndx, VectorXiConstRef ispecies) const -> Vector;

    /// Return the amounts of the species (in units of mol)
    auto speciesAmounts(long int ndx) const -> Vector;

    /// Return the amount of the species with index @param ispecies (in units of mol)
    auto speciesAmount(long int ndx, Index ispecies) const -> double;

    /// Return the amount of the species with name @param name (in units of mol)
    auto speciesAmount(long int ndx, std::string name) const -> double;

    /// Return the upper limits for the species dul (in units of mol)
    auto speciesUpperLimits(long int ndx) const -> Vector;
    /// Return the lower limits for the species dll (in units of mol)
    auto speciesLowerLimits(long int ndx) const -> Vector;

    /// Return the mass of the system (in units of kg).
    auto systemMass(long int ndx) const -> double;

    /// Return the volume of the system (in units of m3).
    auto systemVolume(long int ndx) const -> double;

    /// Return the ionic strength of the aqueous phase (in units of molal).
    auto ionicStrength(long int ndx) const -> double;

    /// Return the pH of the aqueous phase.
    auto pH(long int ndx) const -> double;

    /// Return the pe of the aqueous phase.
    auto pe(long int ndx) const -> double;

    /// Return the Eh of the aqueous phase.
    auto Eh(long int ndx) const -> double;

    /// Return the total Gibbs energy of the system.
    auto systemGibbsEnergy(long int ndx) const -> double;

    /// Return the densities of the phases (in units of kg/m3).
    auto phaseDensities(long int ndx) const -> Vector;

    /// Return the density of the phase @param iphase (in units of kg/m3).
    auto phaseDensity(long int ndx, Index iphase) const -> double;

    /// Return the masses of the phases (in units of kg).
    auto phaseMasses(long int ndx) const -> Vector;

    /// Return the mass of the phase @param iphase (in units of kg).
    auto phaseMass(long int ndx, Index iphase) const -> double;

    /// Return the molar amounts of the phases (in units of mol).
    auto phaseAmounts(long int ndx) const -> Vector;

    /// Return the molar amount of the phase @param iphase (in units of mol).
    auto phaseAmount(long int ndx, Index iphase) const -> double;

    /// Return the volumes of the phases (in units of m3).
    auto phaseVolumes(long int ndx) const -> Vector;

    /// Return the volume of the phase @param iphase (in units of m3).
    auto phaseVolume(long int ndx, Index iphase) const -> double;

    /// Return the enthalpies of the phases (in units of J).
    auto phaseEnthalpies(long int ndx) const -> Vector;

    /// Returns the enthalpy of the phase @param iphase (in units of J).
    auto phaseEnthalpy(long int ndx, Index iphase) const -> double;

    /// Return the entropies of the phases (in units of J/K).
    auto phaseEntropies(long int ndx) const -> Vector;

    /// Returns the entropy of the phase @param iphase (in units of J/K).
    auto phaseEntropy(long int ndx, Index iphase) const -> double;

    /// Return the heat capacities Cp of the phases (in units of J/k).
    auto phaseHeatCapacitiesConstP(long int ndx) const -> Vector;

    /// Returns the heat capacity Cp of the phase @param iphase (in units of J/K).
    auto phaseHeatCapacityConstP(long int ndx, Index iphase) const -> double;

    /// Return the saturation (stability) indices  of the phases (in log10 units).
    auto phaseSatIndices(long int ndx) const -> Vector;

    /// Returns saturation (stability) index (in log10 units) of phase @param iphase.
    auto phaseSatIndex(long int ndx, Index iphase) const -> double;

    /// Return the molar Gibbs energies of the phases (in units of J/mol).
    auto phaseMolarGibbsEnergy(long int ndx, Index iphase) const -> double;

    /// Return the molar enthalpy of the phase @param iphase (in units of J/mol).
    auto phaseMolarEnthalpy(long int ndx, Index iphase) const -> double;

    /// Return the molar volume of the phase @param iphase (in units of m3/mol).
    auto phaseMolarVolume(long int ndx, Index iphase) const -> double;

    /// Return the molar entropy of the phase @param iphase (in units of J/(mol*K)).
    auto phaseMolarEntropy(long int ndx, Index iphase) const -> double;

    /// Return the mole fractions of the species.
    auto moleFractions(long int ndx) const -> Vector;

    /// Return the ln activity coefficients of the species (in mole fraction scale).
    auto lnActivityCoefficients(long int ndx) const -> Vector;

    /// Return the ln activities of the species.
    auto lnActivities(long int ndx) const -> Vector;

    /// Return the ln of concentrations of the species in their respective phases.
    auto lnConcentrations(long int ndx) const -> Vector;

    /// Return the chemical potentials of species (in units of J/mol).
    auto chemicalPotentials(long int ndx) const -> Vector;

    /// Return the molar isobaric heat capacity of phase @param iphase (in units of J/(mol*K)).
    auto phaseMolarHeatCapacityConstP(long int ndx, Index iphase) const -> double;

    /// Return the specific isobaric heat capacity of phase @param iphase (in units of J/(kg*K)).
    auto phaseSpecificHeatCapacityConstP(long int ndx, Index iphase) const -> double;

    /// Return the specific Gibbs energy of the phase @param iphase (in units of J/kg).
    auto phaseSpecificGibbsEnergy(long int ndx, Index iphase) const -> double;

    /// Return the specific enthalpy of the phase @param iphase (in units of J/kg).
    auto phaseSpecificEnthalpy(long int ndx, Index iphase) const -> double;

    /// Return the specific volume of the phase @param iphase (in units of m3/kg).
    auto phaseSpecificVolume(long int ndx, Index iphase) const -> double;

    /// Return the specific entropy of the phase @param iphase (in units of J/(kg*K)).
    auto phaseSpecificEntropy(long int ndx, Index iphase) const -> double;

    /// Return the total enthalpy of the system.
    auto systemEnthalpy(long int ndx) const -> double;

    /// Return the total entropy of the system.
    auto systemEntropy(long int ndx) const -> double;

    /// Return the total heat capacity Cp of the system.
    auto systemHeatCapacityConstP(long int ndx) const -> double;

    //---------------------------------------------------------------

    /// TDB Return the molalities of the species. Aquatic systems only (assuming
    /// the aqueous phase is the first one and H2O-solvent is the last species in it)
    //auto speciesMolalities() const -> VectorConstRef;

    /// TBD Return the standard  molar internal energies of the species (in units of J/mol).
    // auto standardMolarInternalEnergy(Index ispecies) const -> double;

    /// TBD Return the standard  molar Helmholtz energies of the species (in units of J/mol).
    // auto standardMolarHelmholtzEnergy(Index ispecies) const -> double;

    /// TBD Return the standard  molar isochoric heat capacities of the species (in units of J/(mol*K)).
    // auto standardMolarHeatCapacityConstV(Index ispecies) const -> double;

    /// TBD Return the molar internal energies of the phases (in units of J/mol).
    // auto phaseMolarInternalEnergy(Index iphase) const -> double;

    /// TBD Return the molar Helmholtz energies of the phases (in units of J/mol).
    // auto phaseMolarHelmholtzEnergy(Index iphase) const -> double;

    /// TBD Return the molar isochoric heat capacities of the phases (in units of J/(mol*K)).
    // auto phaseMolarHeatCapacityConstV(Index ispecies) const -> double;

    /// TBD Return the specific internal energies of the phases (in units of J/kg).
    // auto phaseSpecificInternalEnergy(Index iphase) const -> double;

    /// TBD Return the specific Helmholtz energies of the phases (in units of J/kg).
    // auto phaseSpecificHelmholtzEnergy(Index iphase) const -> double;

    /// TBD Return the specific isochoric heat capacities of the phases (in units of J/(kg*K)).
    // auto phaseSpecificHeatCapacityConstV(Index iphase) const -> double;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;

    void setCorrectionbIC(long ndx);
    /// Reallocates chemical engine vectors and matrix upon node initialization from files or json
    /// Used only inside of initialize() or initializeFromJsonStrings()
    auto reallocateEngineArrays() -> void;

};

/// Output the state of the ChemicalEngineArray object.
auto operator<<(std::ostream& out, const ChemicalEngineArray& engine) -> std::ostream&;

} // namespace xGEM2MT
