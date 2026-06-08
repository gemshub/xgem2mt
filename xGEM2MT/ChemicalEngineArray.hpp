/**
 * @file ChemicalEngineArray.hpp
 * @brief Header file for the ChemicalEngineArray class.
 *
 * xGEMS2MT is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
 * This file defines the classes and functions for Gibbs-energy minimization
 * thermodynamic modeling with xGEMS2MT for the list of systems. Units used throughout the API are:
 * - Temperature: Kelvin (K)
 * - Pressure: Pascals (Pa)
 * - Amounts: moles (mol)
 * - Mass: kilograms (kg)
 * - Volume: cubic meters (m³)
 * - Energies: Joules (J) or Joules per mole (J/mol) as appropriate.
 *
 * Example applications of these functions can be found in the demos at:
 * /demos/
 *
 * @note The examples below are indicative; please ensure that your vectors
 *       and matrices (such as VectorConstRef) are properly defined (e.g., using Eigen).
 *
 * @author S.Dmytriieva, Dmtrii Kulik, G.D. Miron
 * @date 2026
 *
 * license GNU General Public License v3 or later
 */

#pragma once

// C++ includes
#include <iostream>
#include <memory>
#include <string>
#include <optional>
#include <map>

// xGEMSMT includes
#include <xGEM2MT/Index.hpp>
#include <xGEM2MT/Eigen.hpp>

namespace xGEM2MT {

class MassTransport;

/**
   * @brief Updates logger settings.
   *
   * Configures logging behavior for the xGEMS2MT.
   *
   * @param use_cout If true, prints log messages to stdout.
   * @param logfile_name The filename for a rotating log file. If empty, file logging is disabled.
   * @param log_level Verbosity level for logging (trace = 0, debug = 1, info = 2, warn = 3, err = 4, critical = 5, off = 6 ).
   *
   * @code
   * // Example: Enable stdout logging and log to "chemical_log.txt" with verbosity level 2.
   * xGEMS2MT::update_loggers(true, "chemical_log.txt", 2);
   * @endcode
   */
void update_loggers(bool use_cout, const std::string &logfile_name, size_t log_level);

/**
   * @brief Options for configuring a ChemicalEngineArray instance.
   *
   * This structure contains settings (e.g., initial approximation method) that affect
   * the performance and accuracy of the equilibrium computations.
   *
   * @note The member `warmstart` indicates whether a "smart" (fast but possibly less
   *       accurate) initial approximation is used.
   *
   * @code
   * // Example: Disable warm start when a more robust (cold start) solution is needed.
   * xGEMS2MT::ChemicalEngineOptions options;
   * options.warmstart = false;
   * @endcode
   */
struct ChemicalEngineArrayOptions
{
    /// The flag that indicates if smart start initial approximation is used
    bool warmstart = false;
    /// Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    char useSIA='*';
    /// Extern actual time iterator
    long int step=0;
    /// Cutoff for IC amount differences in the node between time steps (mol, 1e-9)
    double cdv=1e-9;
    /// Cutoff for minimal amounts of IC in node bulk compositions (mol, 1e-12)
    double cez=1e-12;

    std::string to_json_string(bool dence=false) const
    { return {};}
};


/**
   * @class ChemicalEngineArray
   * @brief Class for equilibrium computations and thermodynamic analysis for the list of the systems.
   *
   * The ChemicalEngineArray performs Gibbs energy minimization to compute the equilibrium
   * state of a chemical systems. Its API provides methods to load system data, update
   * component amounts, set options, and query resulting thermodynamic properties.
   */
class ChemicalEngineArray : public std::enable_shared_from_this<ChemicalEngineArray>
{

public:
    // Initialization and Configuration
    /**
     * @brief Default constructor.
     *
     * Creates an empty ChemicalEngineArray instance for the 1D arrangement of nodes.
     * Use one of the initialization routines (such as initialize() or
     * initializeFromJsonStrings()) to configure the system.
     *
     * @param size (size_t) Number of nodes (systems) in the list.
     *
     * @code
     * // Example: Create an engine and then initialize it with a system file.
     * xGEMS2MT::ChemicalEngineArray engine(5);
     * std::vector<long int> nodeTypes={1,0,0,0,0};
     * std::string gems_init_lst = "nodearray-data/CalcColumn-dat.lst";
     * std::string gems_init_dbr_lst = "nodearray-data/CalcColumn-dbr.lst";
     * engine.initialize(gems_init_lst, gems_init_dbr_lst, nodeTypes); // Temperature in K, Pressure in Pa.
     * @endcode
     */
    ChemicalEngineArray(size_t size);

    /**
     * @brief Destructor.
     *
     * Cleans up resources. No explicit example is necessary.
     */
    virtual ~ChemicalEngineArray();

    /// @brief Deleted copy constructor.
    ChemicalEngineArray(const ChemicalEngineArray& other) = delete;

    /// @brief Deleted assignment operator.
    auto operator=(ChemicalEngineArray other) -> ChemicalEngineArray& = delete;

    /**
     * @brief Sets the options for the ChemicalEngine.
     *
     * @param options (const ChemicalEngineOptions&) Configuration options.
     *
     * @code
     * // Example: Disable warm start for higher accuracy.
     * xGEMS2MT::ChemicalEngineArrayOptions opt;
     * opt.warmstart = false;
     * engine.setOptions(opt);
     * @endcode
     */
    auto setOptions(const ChemicalEngineArrayOptions& options) -> void;

    /**
     * @brief Get the current chemical engine options.
     *
     * Returns a constant reference to the internal ChemicalEngineOptions,
     * which control solver behavior, convergence criteria, and output flags.
     * These options can be inspected or copied for modification.
     *
     * @return const ChemicalEngineArrayOptions& The current engine options.
     *
     * @code
     * // Example: Disable warm start for higher accuracy.
     * xGEMS2MT::ChemicalEngineOptions opt = engine.options();
     * opt.warmstart = false;
     * engine.setOptions(opt);
     * @endcode
     */
    auto options() const -> const ChemicalEngineArrayOptions&;

    /**
     * @brief Configures the engine to use a warm (smart) start.
     *
     * Uses previous equilibrium as initial guess (faster convergence).
     *
     * @code
     * // Example: Set smart initial approximation.
     * engine.setWarmStart();
     * @endcode
     */
    auto setWarmStart() -> void;

    /**
     * @brief Configures the engine to use a cold start.
     *
     * Uses a simplex LP initial guess (slower but may yield more accurate results).
     *
     * @code
     * // Example: Force a cold start.
     * engine.setColdStart();
     * @endcode
     */
    auto setColdStart() -> void;

    /**
     * @brief Get the total number of nodes in the node array
     *
     * @return size_t The size of the list.
     *
     * @code
     * auto size = engine.nNodes();
     * @endcode
     */
    auto nNodes() const -> size_t;

    /**
     * @brief Initializes the ChemicalEngineArray from a GEM-Selektor project files.
     *
     * Reads the files and configures the chemical systems defined inside it.
     *
     * @param ipm_filename (string) Path to the system definition file (e.g., "my-system-dat.lst").
     * @param dbr_filename (string) Path to the text file containing the list of DBR input files.
     * @param node_types (std::vector<long int>) The initial node contents from DATABR files
     *                           will be distributed among nodes in array.
     *
     * @code
     * xGEMS2MT::ChemicalEngineArray engine(5);
     * std::vector<long int> node_types={1,0,0,0,0};
     * std::string gems_init_lst = "nodearray-data/CalcColumn-dat.lst";
     * std::string gems_init_dbr_lst = "nodearray-data/CalcColumn-dbr.lst";
     * engine.initialize(gems_init_lst, gems_init_dbr_lst, node_types);
     * @endcode
     */
    auto initialize(std::string ipm_filename, std::string dbr_filename, std::vector<long int> node_types = {}) -> void;

    /**
     * @brief Initializes the ChemicalEngineArray using JSON strings.
     *
     * Configures the engine from three JSON strings:
     * - dch_json: Defines the chemical system.
     * - ipm_json: Contains calculation parameters and settings.
     * - dbr_json: Specifies the input node composition.
     * These strings can be obtained from the GEM-Selektor exported GEMS3K system files as JSON (*-imp.json, *-dbc.json, *-dbr.json).
     *
     * @param dch_json (string) JSON string describing the chemical system.
     * @param ipm_json (string) JSON string with parameter and algorithm settings.
     * @param dbr_json (std::vector<std::string>) JSON strings containing input nodes composition details.
     * @param node_types (std::vector<long int>) The initial node contents from DATABR files will be distributed among
     *                   nodes in array according to the distribution index list.
     *
     * @code
     * // Example: Initialize using JSON-based input.
     * std::string dch = "{}";
     * std::string ipm = "{}";
     * std::vector<std::string> dbr = {"{}"};
     * engine.initializeFromJsonStrings(dch, ipm, dbr);
     * @endcode
     */
    auto initializeFromJsonStrings(std::string dch_json, std::string ipm_json,
                                   std::vector<std::string> dbr_json, std::vector<long int> node_types = {}) -> void;

    /**
     * @brief Dump the MULTI, DATACH and  all DATABR files.
     *
     * @param filepath (string) IPM work structure file path&name
     * @param put_all (bool) If true print array of nodes for previous time point
     *                          and for current time point
    */
    auto dumpGEMS3K(const std::string& filepath, bool put_all = false) const -> void;

    /**
     * @brief Copies array of nodes for current time point to the same place in the
     * array of nodes for previous time point.
     *
     * @param always (bool) If false, copy only the nodes with essential changes
     *                      into the element quantities.
     */
    auto nextTimePoint(bool always) -> void;

    /**
     * @brief Reads a DBR file from disk to the node.
     *
     * Loads the structure of node-dependent data (composition, T, P and other).
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param filename (string) Path to the input DBR file (e.g., "*-dbr.json/.dat").
     *
     * @code
     * // Example: Read system composition from a DBR file to node 0.
     * engine.readDbrFromFile(0, "my-system-dbr-0-0000.json");
     * @endcode
     */
    auto readDbrFromFile(size_t ndx, std::string filename) -> void;

    /**
     * @brief Reads System DBR composition from a JSON string to the node.
     *
     * Loads the structure of node-dependent data (composition, T, P and other).
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param dbr_json (string) JSON string containing DBR composition data.
     *
     * @code
     * // Example: Update composition dynamically from a JSON string into node 0.
     * engine.readDbrFromJsonString(0, "{}");
     * @endcode
     */
    auto readDbrFromJsonString(size_t ndx, std::string dbr_json) -> void;

    /**
     * @brief Writes the current DBR for the node to a file (key-value format).
     *
     * Saves the state of the system (after equilibrium calculations) into a DBR file.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param filename (string) Path to the output DBR file (e.g., "my-system-dbr-0-0001.dat").
     *
     * @code
     * // Example: Save the current state to"my-system-dbr-0-0001.dat".
     * engine.writeDbrToFile(0, "my-system-dbr-0-0001.dat");
     * @endcode
     */
    auto writeDbrToFile(size_t ndx, std::string filename) -> void;


    /**
     * @brief Returns the current DBR as a JSON string for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return A JSON string representing the system composition.
     *
     * @code
     * // Example: Retrieve the system state as JSON for the node 0.
     * std::string json_state = engine.writeDbrToJsonString(0);
     * std::cout << json_state << std::endl;
     * @endcode
     */
    auto writeDbrToJsonString(size_t ndx) -> const std::string;

    /**
     * @brief Returns the number of elements in the system.
     *
     * @return (size_t) Total number of elements.
     *
     * @code
     * // Example: Get the number of elements.
     * size_t nElements = engine.numElements();
     * std::cout << "Elements: " << nElements << std::endl;
     * @endcode
     */
    auto numElements() const -> size_t;

    /**
     * @brief Returns the number of species in the system.
     *
     * @return (size_t) Total number of chemical species.
     *
     * @code
     * // Example: Display the number of species.
     * size_t speciesCount = engine.numSpecies();
     * std::cout << "Species: " << speciesCount << std::endl;
     * @endcode
     */
    auto numSpecies() const -> size_t;

    /**
     * @brief Returns the number of phases in the system.
     *
     * @return (size_t) Total number of phases.
     *
     * @code
     * // Example: Print the number of phases.
     * size_t phases = engine.numPhases();
     * std::cout << "Phases: " << phases << std::endl;
     * @endcode
     */
    auto numPhases() const -> size_t;

    /**
     * @brief Returns the number of species in a given phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (size_t) Index of the phase.
     * @return (size_t) Number of species in that phase.
     *
     * @code
     * // Example: For phase with index 0.
     * size_t sppInPhase = engine.numSpeciesInPhase(0);
     * std::cout << "Species in phase 0: " << sppInPhase << std::endl;
     * @endcode
     */
    auto numSpeciesInPhase(size_t iphase) const -> size_t;

    /**
     * @brief Returns the number of species in a phase identified by name.
     *
     * @param phase (std::string) Phases name.
     * @return (size_t) Number of species in that phase.
     *
     * @code
     * // Example: For phase with index 0.
     * size_t sppInPhase = engine.numSpeciesInPhase("aq_gen");
     * std::cout << "Species in phase aq_gen: " << sppInPhase << std::endl;
     * @endcode
     */
    auto numSpeciesInPhase(std::string phase) const -> size_t;

    /**
     * @brief Returns the number of species in a given phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (size_t) Index of the phase.
     * @return (size_t) Number of species in that phase.
     *
     * @code
     * // Example: For phase with index 0.
     * size_t sppInPhase = engine.numSpeciesInPhase_i(0);
     * std::cout << "Species in phase 0: " << sppInPhase << std::endl;
     * @endcode
     */
    auto numSpeciesInPhase_i(size_t iphase) const -> size_t;

    /**
     * @brief Returns the name of an element.
     * Access specified element with bounds checking.
     *
     * @param ielement (size_t) Index of the element.
     * @return (std::string) Element name.
     *
     * @code
     * // Example: Get the name of the first element.
     * std::string element = engine.elementName(0);
     * std::cout << "Element 0: " << element << std::endl;
     * @endcode
     */
    auto elementName(size_t ielement) const -> std::string;

    /**
     * @brief Returns the name of an element.
     * Access specified element without bounds checking.
     *
     * @param ielement (size_t) Index of the element.
     * @return (std::string) Element name.
     *
     * @code
     * // Example: Get the name of the first element.
     * std::string element = engine.elementName_i(0);
     * std::cout << "Element 0: " << element << std::endl;
     * @endcode
     */
    auto elementName_i(size_t ielement) const -> std::string;

    /**
     * @brief Returns the name of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (size_t) Index of the species.
     * @return (std::string) Species name.
     *
     * @code
     * // Example: Retrieve the name of species 0.
     * std::string spName = engine.speciesName(0);
     * std::cout << "Species 0: " << spName << std::endl;
     * @endcode
     */
    auto speciesName(size_t ispecies) const -> std::string;

    /**
     * @brief Returns the name of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (size_t) Index of the species.
     * @return (std::string) Species name.
     *
     * @code
     * // Example: Retrieve the name of species 0.
     * std::string spName = engine.speciesName_i(0);
     * std::cout << "Species 0: " << spName << std::endl;
     * @endcode
     */
    auto speciesName_i(size_t ispecies) const -> std::string;

    /**
     * @brief Returns the electrical charge of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (size_t) Index of the species.
     * @return (double) Charge of the species.
     *
     * @code
     * // Example: Get the charge of species with index 0.
     * double charge = engine.speciesCharge(0);
     * std::cout << "Charge: " << charge << std::endl;
     * @endcode
     */
    auto speciesCharge(size_t ispecies) const -> double;

    /**
     * @brief Returns the electrical charge of a species identified by name.
     *
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Charge of the species.
     *
     * @code
     * // Example: Get the charge of species with index 0.
     * double charge = engine.speciesCharge("SiO2");
     * std::cout << "Charge: " << charge << std::endl;
     * @endcode
     */
    auto speciesCharge(std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the electrical charge of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (size_t) Index of the species.
     * @return (double) Charge of the species.
     *
     * @code
     * // Example: Get the charge of species with index 0.
     * double charge = engine.speciesCharge_i(0);
     * std::cout << "Charge: " << charge << std::endl;
     * @endcode
     */
    auto speciesCharge_i(size_t ispecies) const -> double;

    /**
     * @brief Returns the name of the phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (size_t) Index of the phase.
     * @return (std::string) Phase name.
     *
     * @code
     * // Example: Retrieve the name of phase 0.
     * std::string phase = engine.phaseName(0);
     * std::cout << "Phase 0: " << phase << std::endl;
     * @endcode
     */
    auto phaseName(size_t iphase) const -> std::string;

    /**
     * @brief Returns the name of the phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (size_t) Index of the phase.
     * @return (std::string) Phase name.
     *
     * @code
     * // Example: Retrieve the name of phase 0.
     * std::string phase = engine.phaseName_i(0);
     * std::cout << "Phase 0: " << phase << std::endl;
     * @endcode
     */
    auto phaseName_i(size_t iphase) const -> std::string;

    /**
     * @brief Returns the index of an element by name.
     *
     * @param element (std::string) Element symbol or name.
     * @return (size_t) Index of the element (or total elements if not found).
     *
     * @code
     * // Example: Get the index of Oxygen ("O").
     * size_t idx = engine.indexElement("O");
     * std::cout << "size_t of O: " << idx << std::endl;
     * @endcode
     */
    auto indexElement(std::string element) const -> size_t;

    /**
     * @brief Returns the index of a species by name.
     *
     * @param species (std::string) Name of the species.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (size_t) Index of the species (or total species if not found).
     *
     * @code
     * // Example: Get the index of water ("H2O@").
     * size_t idx = engine.indexSpecies("H2O@");
     * std::cout << "size_t of H2O@: " << idx << std::endl;
     * @endcode
     */
    auto indexSpecies(std::string species, std::optional<std::string> phase = std::nullopt) const -> size_t;

    /**
     * @brief Returns all species indexes and the corresponding phase names matching the specified name.
     *
     * @param species (std::string) Name of the species.
     * @return (Map) Dictionary of phases and corresponding species index.
     *
     * @code
     * // Example: Retrieve all phase and indices where the species "Anorthite" appears.
     * std::map<std::string, long int> allIndices = engine.indexSpeciesAll("Anorthite");
     * @endcode
     */
    auto indexSpeciesMap(std::string species) const -> std::map<std::string, long int>;

    /**
     * @brief Returns all indices of species matching the specified name.
     *
     * @param species (std::string) Name of the species.
     * @return (VectorXi) Vector of indices.
     *
     * @code
     * // Example: Retrieve all indices where the species "H2O@" appears.
     * VectorXi allIndices = engine.indexSpeciesAll("H2O@");
     * @endcode
     */
    auto indexSpeciesAll(std::string species) const -> VectorXi;

    /**
     * @brief Returns the index of a phase by name.
     *
     * @param phase (std::string) Name of the phase.
     * @return (size_t) Index of the phase (or total phases if not found).
     *
     * @code
     * // Example: Get index for the "aq_gen" phase.
     * size_t aqueousIndex = engine.indexPhase("aq_gen");
     * @endcode
     */
    auto indexPhase(std::string phase) const -> size_t;

    /**
     * @brief Returns all indices for phases that match the given name.
     *
     * @param phase (std::string) Name of the phase.
     * @return (VectorXi) Vector of indices for matching phases.
     *
     * @code
     * // Example: Retrieve all indices for phase "SiO2".
     * VectorXi phaseIndices = engine.indexPhaseAll("SiO2");
     * @endcode
     */
    auto indexPhaseAll(std::string phase) const -> VectorXi;

    /**
     * @brief Returns the index of the phase containing a given species.
     *
     * @param ispecies (size_t) Index of the species.
     * @return (size_t) Index of the phase that contains the species.
     *
     * @code
     * // Example: Find the phase of the species with index 0.
     * size_t phaseIdx = engine.indexPhaseWithSpecies(0);
     * @endcode
     */
    auto indexPhaseWithSpecies(size_t ispecies) const -> size_t;

    /**
     * @brief Returns the index of the first species in a specified phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (size_t) Index of the phase.
     * @return (size_t) Index of the first species in that phase or numSpecies if illegal iphase.
     *
     * @code
     * // Example: Get the first species index in phase 0.
     * size_t firstSpec = engine.indexFirstSpeciesInPhase(0);
     * @endcode
     */
    auto indexFirstSpeciesInPhase(size_t iphase) const -> size_t;

    /**
     * @brief Returns the index of the first species in a specified by name phase.
     *
     * @param phase (std::string) Phase name.
     * @return (size_t) Index of the first species in that phase or numSpecies if illegal iphase.
     *
     * @code
     * // Example: Get the first species index in phase aq_gen.
     * size_t firstSpec = engine.indexFirstSpeciesInPhase("aq_gen");
     * @endcode
     */
    auto indexFirstSpeciesInPhase(std::string phase) const -> size_t;

    /**
     * @brief Returns the index of the first species in a specified phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (size_t) Index of the phase.
     * @return (size_t) Index of the first species in that phase or numSpecies if illegal iphase.
     *
     * @code
     * // Example: Get the first species index in phase 0.
     * size_t firstSpec = engine.indexFirstSpeciesInPhase_i(0);
     * @endcode
     */
    auto indexFirstSpeciesInPhase_i(size_t iphase) const -> size_t;

    /**
     * @brief Returns molar masses of elements.
     *
     * @return (VectorConstRef) Vector of molar masses (kg/mol) for each element.
     *
     * @code
     * // Example: Retrieve the element molar masses vector.
     * auto molarMasses = engine.elementMolarMasses();
     * @endcode
     */
    auto elementMolarMasses() const -> VectorConstRef;

    /**
     * @brief Returns molar masses of species.
     *
     * @return (VectorConstRef) Vector of species molar masses (kg/mol).
     *
     * @code
     * // Example: Get the species molar masses.
     * auto spMolarMasses = engine.speciesMolarMasses();
     * @endcode
     */
    auto speciesMolarMasses() const -> VectorConstRef;

    /**
     * @brief Returns the formula matrix (elements x species).
     *
     * @return (MatrixConstRef) Formula matrix where rows correspond to elements and
     * columns correspond to species.
     *
     * @code
     * // Example: Retrieve the formula matrix.
     * auto fMatrix = engine.formulaMatrix();
     * @endcode
     */
    auto formulaMatrix() const -> MatrixConstRef;

    /**
     * @brief Sets the complete species amounts vector for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param n (VectorConstRef) Vector of species amounts (mol).
     *
     * @code
     * // Example: Set all species amounts to 1.0 mol.
     * Eigen::VectorXd n(engine.numSpecies());
     * n.setOnes();
     * engine.setSpeciesAmounts(0, n);
     * @endcode
     */
    auto setSpeciesAmounts(size_t ndx, VectorConstRef n) -> void;

    /**
     * @brief Sets the amount for a species identified by name (as in the GEMS system) for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param name (std::string) Species name.
     * @param amount (double) New amount in mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Set the amount of "CaSO4@" to 0.01 mol.
     * engine.setSpeciesAmount(0, "CaSO4@", 0.01);
     * @endcode
     */
    auto setSpeciesAmount(size_t ndx, std::string name, double amount, std::optional<std::string> phase = std::nullopt) -> void;

    /**
     * @brief Sets the amount for a species identified by its index.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies (size_t) Index of the species.
     * @param amount (double) New amount in mol.
     *
     * @code
     * // Example: Set species at index 0 to 0.01 mol.
     * engine.setSpeciesAmount(0, 0, 0.01);
     * @endcode
     */
    auto setSpeciesAmount(size_t ndx, size_t ispecies, double amount) -> void;

    /**
     * @brief Sets the amount for a species identified by its index for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies (size_t) Index of the species.
     * @param amount (double) New amount in mol.
     *
     * @code
     * // Example: Set species at index 0 to 0.01 mol.
     * engine.setSpeciesAmount_i(0, 0, 0.01);
     * @endcode
     */
    auto setSpeciesAmount_i(size_t ndx, size_t ispecies, double amount) -> void;

    /**
     * @brief Sets an upper bound for a species identified by name for the node.
     *
     * If amount < 0, resets to default 1e6.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param name (std::string) Species name.
     * @param amount (double) Upper limit in mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Set the upper limit for "SiO2" to 0.1 mol.
     * engine.setSpeciesUpperLimit(0, "SiO2", 0.1);
     * @endcode
     */
    auto setSpeciesUpperLimit(size_t ndx, std::string name, double amount, std::optional<std::string> phase = std::nullopt) -> void;

    /**
     * @brief Sets a lower bound for a species identified by name for the node.
     *
     * If amount < 0, resets to default 0.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param name (std::string) Species name.
     * @param amount (double) Lower limit in mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Set the lower limit for "SiO2" to 0.05 mol.
     * engine.setSpeciesLowerLimit(0, "SiO2", 0.05);
     * @endcode
     */
    auto setSpeciesLowerLimit(size_t ndx, std::string name, double amount, std::optional<std::string> phase = std::nullopt) -> void;

    /**
     * @brief Sets an upper bound (maximum amount allowed to form) for a species identified by its index for the node.
     * Access specified element with bounds checking.
     *
     * If amount < 0, resets to default 1e6.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies (size_t) Species index.
     * @param amount (double) Upper limit in mol.
     *
     * @code
     * // Example: Set upper limit for species index 0.
     * engine.setSpeciesUpperLimit(0, 0, 0.1);
     * @endcode
     */
    auto setSpeciesUpperLimit(size_t ndx, size_t ispecies, double amount) -> void;

    /**
     * @brief Sets a lower bound (minimum amount allowed to form) for a species identified by its index for the node.
     * Access specified element with bounds checking.
     *
     * If amount < 0, resets to default 0.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies (size_t) Species index.
     * @param amount (double) Lower limit in mol.
     *
     * @code
     * // Example: Set lower limit for species index 0.
     * engine.setSpeciesLowerLimit(0, 0, 0.05);
     * @endcode
     */
    auto setSpeciesLowerLimit(size_t ndx, size_t ispecies, double amount) -> void;

    /**
     * @brief Sets an upper bound (maximum amount allowed to form) for a species identified by its index for the node.
     * Access specified element without bounds checking.
     *
     * If amount < 0, resets to default 1e6.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies (size_t) Species index.
     * @param amount (double) Upper limit in mol.
     *
     * @code
     * // Example: Set upper limit for species index 0.
     * engine.setSpeciesUpperLimit_i(0, 0, 0.1);
     * @endcode
     */
    auto setSpeciesUpperLimit_i(size_t ndx, size_t ispecies, double amount) -> void;

    /**
     * @brief Sets a lower bound (minimum amount allowed to form) for a species identified by its index for the node.
     * Access specified element without bounds checking.
     *
     * If amount < 0, resets to default 0.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies (size_t) Species index.
     * @param amount (double) Lower limit in mol.
     *
     * @code
     * // Example: Set lower limit for species index 0.
     * engine.setSpeciesLowerLimit_i(0, 0, 0.05);
     * @endcode
     */
    auto setSpeciesLowerLimit_i(size_t ndx, size_t ispecies, double amount) -> void;

    /**
     * @brief Sets the standard molar Gibbs energy for a species (@ T, P of the system).
     *
     * @param name (std::string) Species name.
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param value (double) Standard molar Gibbs energy in J/mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Set the standard molar Gibbs energy of H2O.
     * engine.setStandardMolarGibbsEnergy("H2O", 298.15, 101325, -237140); // Value in J/mol.
     * @endcode
     */
    auto setStandardMolarGibbsEnergy(std::string name, double TK, double P, double value, std::optional<std::string> phase = std::nullopt) -> void;

    /**
     * @brief Sets the standard molar Gibbs energy for a species (@ T, P of the system).
     * Access specified element with bounds checking.
     *
     * @param ispecies (size_t) Index of the species.
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param value (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: Set the standard molar Gibbs energy of H2O.
     * engine.setStandardMolarGibbsEnergy(0, 298.15, 101325, -237.13); // Value in J/mol.
     * @endcode
     */
    auto setStandardMolarGibbsEnergy(size_t ispecies, double TK, double P, double value) -> void;

    /**
     * @brief Sets the standard molar Gibbs energy for a species (@ T, P of the system).
     * Access specified element without bounds checking.
     *
     * @param ispecies (size_t) Index of the species.
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param value (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: Set the standard molar Gibbs energy of H2O.
     * engine.setStandardMolarGibbsEnergy_i(0, 298.15, 101325, -237.13); // Value in J/mol.
     * @endcode
     */
    auto setStandardMolarGibbsEnergy_i(size_t ispecies, double TK, double P, double value) -> void;

    /**
     * @brief Sets upper limits for all species for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param n (VectorConstRef) Vector of upper limits in mol.
     *
     * @code
     * // Example: Set all species' upper limits to values provided in vector "limits".
     * engine.setSpeciesUpperLimits(0, limits);
     * @endcode
     */
    auto setSpeciesUpperLimits(size_t ndx, VectorConstRef n) -> void;

    /**
     * @brief Sets lower limits for all species for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param n (VectorConstRef) Vector of lower limits in mol.
     *
     * @code
     * // Example: Set all species' lower limits.
     * engine.setSpeciesLowerLimits(0, lowerLimits);
     * @endcode
     */
    auto setSpeciesLowerLimits(size_t ndx, VectorConstRef n) -> void;

    /**
     * @brief Sets the pressure and temperature without computing equilibrium for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param P Pressure in Pascals (Pa).
     * @param T Temperature in Kelvin (K).
     * @return (bool) true if PT was set correctly, false if out of range.
     *
     * @code
     * // Example: Set PT to 101325 Pa and 298.15 K.
     * bool no_error = engine.setPT(0, 101325, 298.15);
     * std::cout << (no_error ? "PT set correctly" : "PT error") << std::endl;
     * @endcode
     */
    auto setPT(size_t ndx, double P, double T) const -> bool;

    /**
     * @brief Sets the amounts of elements (vector b) without computing equilibrium for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param b (VectorConstRef) Vector of element amounts in mol (order based on the index of elements in the chemical system).
     *
     * @code
     * // Example: Set elements amounts of a solution with 0.001 M MgCl2 0.001 M CO2
     * // elements order: C, Ca, Cl, H, Mg, O, Zz (charge)
     * Eigen::VectorXd b(7);
     * b << 0.001, 1e-09, 0.004, 110.6837, 0.002, 55.34405, 0;  // in mol
     * engine.setB(0, b);
     * @endcode
     */
    auto setB(size_t ndx, VectorConstRef b) -> void;

    /**
     * @brief Re-equilibrates the system with (or without) a warm start for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param warmstart (bool) If true, uses previous data as an initial guess.
     * @return (int) Return code of the re-equilibration.
     *
     * @code
     * // Example: Re-equilibrate with a warm start.
     * int status = engine.reequilibrate(true);
     * @endcode
     *
     * - 0: No GEM re-calculation needed
     * - 1: Need GEM calculation with LPP (automatic) initial approximation (AIA)
     * - 2: OK after GEM calculation with LPP AIA
     * - 3: Bad (not fully trustful) result after GEM calculation with LPP AIA
     * - 4: Failure (no result) in GEM calculation with LPP AIA
     * - 5: Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
     * - 6: OK after GEM calculation with SIA
     * - 7: Bad (not fully trustful) result after GEM calculation with SIA
     * - 8: Failure (no result) in GEM calculation with SIA
     * - 9: Terminal error in GEMS3K (e.g., memory corruption). Restart required.
     */
    auto reequilibrate(size_t ndx, bool warmstart) -> long int;

    /**
     * @brief Re-equilibrates the all systems with (or without) a warm start.
     *
     * @param warmstart (bool) If true, uses previous data as an initial guess.
     *
     * @code
     * // Example: Re-equilibrate with a warm start.
     * engine.reequilibrate(true);
     * @endcode
     *
     */
    auto reequilibrate(bool warmstart) -> void;

    /**
     * @brief Equilibrates the system in node ndx, according options and node settings.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (bool) Return true if OK after GEM calculation with LPP AIA or SIA.
     *
     * @code
     * // Example: Re-equilibrate.
     * auto status = engine.equilibrate(1);
     * @endcode
     *
     */
    auto equilibrate(size_t ndx) -> bool;

    /**
     * @brief Equilibrates all systems, according options and node settings.
     *
     * @return (bool) Return true if OK after GEM calculation with LPP AIA or SIA for all systems.
     *
     * @code
     * // Example: Re-equilibrate.
     * auto status = engine.equilibrate();
     * @endcode
     *
     */
    auto equilibrate() -> bool;

    /**
     * @brief Checks if the equilibrium calculation converged for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (bool) True if converged, false otherwise.
     *
     * @code
     * // Example: Check convergence status.
     * bool isConverged = engine.converged(0);
     * @endcode
     */
    auto converged(size_t ndx) const -> bool;

    /**
     * @brief Returns the number of iterations performed during equilibrium for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (size_t) The number of iterations.
     *
     * @code
     * // Example: Display the number of iterations.
     * size_t iters = engine.numIterations(0);
     * std::cout << "Iterations: " << iters << std::endl;
     * @endcode
     */
    auto numIterations(size_t ndx) const -> size_t;

    /**
     * @brief Returns the elapsed time of the equilibrium computation.
     *
     * @return (double) Elapsed time in seconds.
     *
     * @code
     * // Example: Get the computation time.
     * double timeSec = engine.elapsedTime();
     * std::cout << "Time: " << timeSec << " s" << std::endl;
     * @endcode
     */
    auto elapsedTime() const -> double;

    /**
     * @brief Returns the current temperature of the system for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Temperature in Kelvin (K).
     *
     * @code
     * // Example: Get the system temperature.
     * double temp = engine.temperature(0);
     * std::cout << "Temperature: " << temp << " K" << std::endl;
     * @endcode
     */
    auto temperature(size_t ndx) const -> double;

    /**
     * @brief Returns the current pressure of the system for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Pressure in Pascals (Pa).
     *
     * @code
     * // Example: Retrieve the system pressure.
     * double press = engine.pressure(0);
     * std::cout << "Pressure: " << press << " Pa" << std::endl;
     * @endcode
     */
    auto pressure(size_t ndx) const -> double;

    /**
     * @brief Returns the amounts of the elements for the node.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of element amounts in mol.
     *
     * @code
     * // Example: Get element amounts.
     * auto eAmounts = engine.elementAmounts(0);
     * @endcode
     */
    auto elementAmounts(size_t ndx) const -> Vector;

    /**
     * @brief Returns the amounts of elements in a specified phase for the node.
     * Access specified element with bounds checking.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Index of the phase.
     * @return A vector containing the amounts of each element in the phase (in mol).
     *
     * **Example usage:**
     * @code
     * // Retrieve the element amounts in phase index 0.
     * Eigen::VectorXd elementAmounts = engine.elementAmountsInPhase(0, 0);
     * std::cout << "Element amounts in phase 0: " << elementAmounts.transpose() << std::endl;
     * @endcode
     */
    auto elementAmountsInPhase(size_t ndx, size_t iphase) const -> Vector;

    /**
     * @brief Returns the amounts of elements in a specified phase for the node.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase Name of the phase.
     * @return A vector containing the amounts of each element in the phase (in mol).
     *
     * **Example usage:**
     * @code
     * // Retrieve the element amounts in phase aq_gen.
     * Eigen::VectorXd elementAmounts = engine.elementAmountsInPhase(0, "aq_gen");
     * std::cout << "Element amounts in phase aq_gen: " << elementAmounts.transpose() << std::endl;
     * @endcode
     */
    auto elementAmountsInPhase(size_t ndx, std::string phase) const -> Vector;

    /**
     * @brief Returns the amounts of elements in a specified phase for the node.
     * Access specified element without bounds checking.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Index of the phase.
     * @return A vector containing the amounts of each element in the phase (in mol).
     *
     * **Example usage:**
     * @code
     * // Retrieve the element amounts in phase index 0.
     * Eigen::VectorXd elementAmounts = engine.elementAmountsInPhase_i(0, 0);
     * std::cout << "Element amounts in phase 0: " << elementAmounts.transpose() << std::endl;
     * @endcode
     */
    auto elementAmountsInPhase_i(size_t ndx, size_t iphase) const -> Vector;

    /**
     * @brief Returns the amounts of elements in a specified group of species for the node.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies A vector containing the indices of the species whose elemental amounts should be retrieved.
     * @return A vector containing the amounts of each element corresponding to the selected species (in mol).
     *
     * **Example usage:**
     * @code
     * // Retrieve element amounts contributed by species indices {0, 2, 5}.
     * Eigen::VectorXi speciesIndices(3);
     * speciesIndices << 0, 2, 5;
     * Eigen::VectorXd elementAmounts = engine.elementAmountsInSpecies(0, speciesIndices);
     * std::cout << "Element amounts from selected species: " << elementAmounts.transpose() << std::endl;
     * @endcode
     */
    auto elementAmountsInSpecies(size_t ndx, VectorXiConstRef ispecies) const -> Vector;

    /**
     * @brief Returns the amounts of the species for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of species amounts in mol.
     *
     * @code
     * // Example: Retrieve the species amounts.
     * auto spAmounts = engine.speciesAmounts(0);
     * @endcode
     */
    auto speciesAmounts(size_t ndx) const -> Vector;

    /**
     * @brief Returns the amount of a species by index for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies (size_t) Index of the species.
     * @return (double) Amount in mol.
     *
     * @code
     * // Example: Get the species amount at index 0.
     * double amount = engine.speciesAmount(0, 0);
     * @endcode
     */
    auto speciesAmount(size_t ndx, size_t ispecies) const -> double;

    /**
     * @brief Returns the amount of a species identified by name for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param name (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Amount in mol.
     *
     * @code
     * // Example: Get the amount of "OH-".
     * double waterAmt = engine.speciesAmount(0, "OH-");
     * @endcode
     */
    auto speciesAmount(size_t ndx, std::string name, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the amount of a species by index for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param ispecies (size_t) Index of the species.
     * @return (double) Amount in mol.
     *
     * @code
     * // Example: Get the species amount at index 0.
     * double amount = engine.speciesAmount_i(0, 0);
     * @endcode
     */
    auto speciesAmount_i(size_t ndx, size_t ispecies) const -> double;

    /**
     * @brief Returns the upper limits for all species for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of upper limits in mol.
     *
     * @code
     * // Example: Fetch all species upper limits.
     * auto upperLimits = engine.speciesUpperLimits(1);
     * @endcode
     */
    auto speciesUpperLimits(size_t ndx) const -> Vector;

    /**
     * @brief Returns the lower limits for all species for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of lower limits in mol.
     *
     * @code
     * // Example: Get the species lower limits.
     * auto lowerLimits = engine.speciesLowerLimits(0);
     * @endcode
     */
    auto speciesLowerLimits(size_t ndx) const -> Vector;

    /**
     * @brief Returns the molalities of the species for the node.
     *
     * Assumes the aqueous phase is the first phase.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Species molalities (mol/kg H2O@).
     *
     * @code
     * // Example: Retrieve species molalities.
     * auto molalities = engine.speciesMolalities(0);
     * @endcode
     */
    auto speciesMolalities(size_t ndx) const -> VectorConstRef;

    /**
     * @brief Returns the mole fractions of the species for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Mole fractions.
     *
     * @code
     * // Example: Get species mole fractions.
     * auto moleFrac = engine.moleFractions(1);
     * @endcode
     */
    auto moleFractions(size_t ndx) const -> Vector;

    /**
     * @brief Returns the ln activity coefficients of the species (mole fraction scale) for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) ln Activity coefficients.
     *
     * @code
     * // Example: Retrieve ln activity coefficients.
     * auto lnActCoeff = engine.lnActivityCoefficients(1);
     * @endcode
     */
    auto lnActivityCoefficients(size_t ndx) const -> Vector;

    /**
     * @brief Returns the ln activities of the species for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) ln Activities.
     *
     * @code
     * // Example: Fetch ln activities.
     * auto lnActivities = engine.lnActivities(0);
     * @endcode
     */
    auto lnActivities(size_t ndx) const -> Vector;

    /**
     * @brief Returns the ln of species concentrations for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) ln concentrations.
     *
     * @code
     * // Example: Get ln concentrations.
     * auto lnConcs = engine.lnConcentrations(1);
     * @endcode
     */
    auto lnConcentrations(size_t ndx) const -> Vector;

    /**
     * @brief Returns the chemical potentials of the species for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Chemical potentials in J/mol.
     *
     * @code
     * // Example: Retrieve species chemical potentials.
     * auto chemPot = engine.chemicalPotentials(1);
     * @endcode
     */
    auto chemicalPotentials(size_t ndx) const -> Vector;

    /**
     * @brief Returns the standard molar Gibbs energy of a species.
     * Access specified element with bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: For species 0.
     * double stdGibbs = engine.standardMolarGibbsEnergy(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarGibbsEnergy(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar Gibbs energy of a species.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: For species 0.
     * double stdGibbs = engine.standardMolarGibbsEnergy(298.15, 101325, "SiO2");
     * @endcode
     */
    auto standardMolarGibbsEnergy(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar Gibbs energy of a species.
     * Access specified element without bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: For species 0.
     * double stdGibbs = engine.standardMolarGibbsEnergy_i(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarGibbsEnergy_i(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar enthalpy of a species.
     * Access specified element with bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Standard molar enthalpy in J/mol.
     *
     * @code
     * double stdEnthalpy = engine.standardMolarEnthalpy(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarEnthalpy(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar enthalpy of a species.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Standard molar enthalpy in J/mol.
     *
     * @code
     * double stdEnthalpy = engine.standardMolarEnthalpy(298.15, 101325, "SiO2");
     * @endcode
     */
    auto standardMolarEnthalpy(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar enthalpy of a species.
     * Access specified element without bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Standard molar enthalpy in J/mol.
     *
     * @code
     * double stdEnthalpy = engine.standardMolarEnthalpy_i(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarEnthalpy_i(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar volume of a species.
     * Access specified element with bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Standard molar volume in m³/mol.
     *
     * @code
     * double stdVolume = engine.standardMolarVolume(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarVolume(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar volume of a species.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Standard molar volume in m³/mol.
     *
     * @code
     * double stdVolume = engine.standardMolarVolume(298.15, 101325, "SiO2");
     * @endcode
     */
    auto standardMolarVolume(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar volume of a species.
     * Access specified element without bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Standard molar volume in m³/mol.
     *
     * @code
     * double stdVolume = engine.standardMolarVolume_i(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarVolume_i(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar entropy of a species.
     * Access specified element with bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Standard molar entropy in J/(mol·K).
     *
     * @code
     * double stdEntropy = engine.standardMolarEntropy(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarEntropy(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar entropy of a species.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Standard molar entropy in J/(mol·K).
     *
     * @code
     * double stdEntropy = engine.standardMolarEntropy(298.15, 101325, "SiO2");
     * @endcode
     */
    auto standardMolarEntropy(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar entropy of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (size_t) Species index.
     * @return (double) Standard molar entropy in J/(mol·K).
     *
     * @code
     * double stdEntropy = engine.standardMolarEntropy_i(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarEntropy_i(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar isobaric heat capacity of a species.
     * Access specified element with bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double heatCap = engine.standardMolarHeatCapacityConstP(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarHeatCapacityConstP(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the standard molar isobaric heat capacity of a species.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double heatCap = engine.standardMolarHeatCapacityConstP(298.15, 101325, "SiO2");
     * @endcode
     */
    auto standardMolarHeatCapacityConstP(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar isobaric heat capacity of a species.
     * Access specified element without bounds checking.
     *
     * @param T Temperature in Kelvin (K).
     * @param P Pressure in Pascals (Pa).
     * @param ispecies (size_t) Species index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double heatCap = engine.standardMolarHeatCapacityConstP_i(298.15, 101325, 0);
     * @endcode
     */
    auto standardMolarHeatCapacityConstP_i(double TK, double P, size_t ispecies) const -> double;

    /**
     * @brief Returns the molar Gibbs energy of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Molar Gibbs energy in J/mol.
     *
     * @code
     * double phaseGibbs = engine.phaseMolarGibbsEnergy(1, 0);
     * @endcode
     */
    auto phaseMolarGibbsEnergy(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar Gibbs energy of a phasefor the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Molar Gibbs energy in J/mol.
     *
     * @code
     * double phaseGibbs = engine.phaseMolarGibbsEnergy(1, "aq_gen");
     * @endcode
     */
    auto phaseMolarGibbsEnergy(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar Gibbs energy of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Molar Gibbs energy in J/mol.
     *
     * @code
     * double phaseGibbs = engine.phaseMolarGibbsEnergy_i(1, 0);
     * @endcode
     */
    auto phaseMolarGibbsEnergy_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar enthalpy of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Molar enthalpy in J/mol.
     *
     * @code
     * double phaseEnthalpy = engine.phaseMolarEnthalpy(1, 0);
     * @endcode
     */
    auto phaseMolarEnthalpy(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar enthalpy of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Molar enthalpy in J/mol.
     *
     * @code
     * double phaseEnthalpy = engine.phaseMolarEnthalpy(1, "aq_gen");
     * @endcode
     */
    auto phaseMolarEnthalpy(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar enthalpy of a phasefor the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Molar enthalpy in J/mol.
     *
     * @code
     * double phaseEnthalpy = engine.phaseMolarEnthalpy_i(1, 0);
     * @endcode
     */
    auto phaseMolarEnthalpy_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar volume of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Molar volume in m³/mol.
     *
     * @code
     * double phaseVol = engine.phaseMolarVolume(1, 0);
     * @endcode
     */
    auto phaseMolarVolume(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar volume of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Molar volume in m³/mol.
     *
     * @code
     * double phaseVol = engine.phaseMolarVolume(1, "aq_gen");
     * @endcode
     */
    auto phaseMolarVolume(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar volume of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Molar volume in m³/mol.
     *
     * @code
     * double phaseVol = engine.phaseMolarVolume_i(1, 0);
     * @endcode
     */
    auto phaseMolarVolume_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar entropy of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Molar entropy in J/(mol·K).
     *
     * @code
     * double phaseEnt = engine.phaseMolarEntropy(1, 0);
     * @endcode
     */
    auto phaseMolarEntropy(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar entropy of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Molar entropy in J/(mol·K).
     *
     * @code
     * double phaseEnt = engine.phaseMolarEntropy(1, "aq_gen");
     * @endcode
     */
    auto phaseMolarEntropy(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar entropy of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Molar entropy in J/(mol·K).
     *
     * @code
     * double phaseEnt = engine.phaseMolarEntropy_i(1, 0);
     * @endcode
     */
    auto phaseMolarEntropy_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar isobaric heat capacity of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double phaseCp = engine.phaseMolarHeatCapacityConstP(1, 0);
     * @endcode
     */
    auto phaseMolarHeatCapacityConstP(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar isobaric heat capacity of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double phaseCp = engine.phaseMolarHeatCapacityConstP(1, "aq_gen");
     * @endcode
     */
    auto phaseMolarHeatCapacityConstP(size_t ndx, std::string phase) const -> double;


    /**
     * @brief Returns the molar isobaric heat capacity of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double phaseCp = engine.phaseMolarHeatCapacityConstP_i(1, 0);
     * @endcode
     */
    auto phaseMolarHeatCapacityConstP_i(size_t ndx, size_t iphase) const -> double;


    /**
     * @brief Returns the specific Gibbs energy of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Specific Gibbs energy in J/kg.
     *
     * @code
     * double specGibbs = engine.phaseSpecificGibbsEnergy(1, 0);
     * @endcode
     */
    auto phaseSpecificGibbsEnergy(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific Gibbs energy of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Specific Gibbs energy in J/kg.
     *
     * @code
     * double specGibbs = engine.phaseSpecificGibbsEnergy(1, "aq_gen");
     * @endcode
     */
    auto phaseSpecificGibbsEnergy(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific Gibbs energy of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Specific Gibbs energy in J/kg.
     *
     * @code
     * double specGibbs = engine.phaseSpecificGibbsEnergy_i(1, 0);
     * @endcode
     */
    auto phaseSpecificGibbsEnergy_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific enthalpy of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Specific enthalpy in J/kg.
     *
     * @code
     * double specEnthalpy = engine.phaseSpecificEnthalpy(1, 0);
     * @endcode
     */
    auto phaseSpecificEnthalpy(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific enthalpy of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Specific enthalpy in J/kg.
     *
     * @code
     * double specEnthalpy = engine.phaseSpecificEnthalpy(1, "aq_gen");
     * @endcode
     */
    auto phaseSpecificEnthalpy(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific enthalpy of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Specific enthalpy in J/kg.
     *
     * @code
     * double specEnthalpy = engine.phaseSpecificEnthalpy_i(1, 0);
     * @endcode
     */
    auto phaseSpecificEnthalpy_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific volume of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Specific volume in m³/kg.
     *
     * @code
     * double specVolume = engine.phaseSpecificVolume(1, 0);
     * @endcode
     */
    auto phaseSpecificVolume(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific volume of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Specific volume in m³/kg.
     *
     * @code
     * double specVolume = engine.phaseSpecificVolume(1, "aq_gen");
     * @endcode
     */
    auto phaseSpecificVolume(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific volume of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Specific volume in m³/kg.
     *
     * @code
     * double specVolume = engine.phaseSpecificVolume_i(1, 0);
     * @endcode
     */
    auto phaseSpecificVolume_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific entropy of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Specific entropy in J/(kg·K).
     *
     * @code
     * double specEntropy = engine.phaseSpecificEntropy(1, 0);
     * @endcode
     */
    auto phaseSpecificEntropy(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific entropy of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Specific entropy in J/(kg·K).
     *
     * @code
     * double specEntropy = engine.phaseSpecificEntropy(1, "aq_gen");
     * @endcode
     */
    auto phaseSpecificEntropy(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific entropy of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Specific entropy in J/(kg·K).
     *
     * @code
     * double specEntropy = engine.phaseSpecificEntropy_i(1, 0);
     * @endcode
     */
    auto phaseSpecificEntropy_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific isobaric heat capacity of a phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Heat capacity in J/(kg·K).
     *
     * @code
     * double specCp = engine.phaseSpecificHeatCapacityConstP(1, 0);
     * @endcode
     */
    auto phaseSpecificHeatCapacityConstP(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the specific isobaric heat capacity of a phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Heat capacity in J/(kg·K).
     *
     * @code
     * double specCp = engine.phaseSpecificHeatCapacityConstP(1, "aq_gen");
     * @endcode
     */
    auto phaseSpecificHeatCapacityConstP(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific isobaric heat capacity of a phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Heat capacity in J/(kg·K).
     *
     * @code
     * double specCp = engine.phaseSpecificHeatCapacityConstP_i(1, 0);
     * @endcode
     */
    auto phaseSpecificHeatCapacityConstP_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the densities of all phases for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of densities in kg/m³.
     *
     * @code
     * auto densities = engine.phaseDensities(0);
     * @endcode
     */
    auto phaseDensities(size_t ndx) const -> Vector;

    /**
     * @brief Returns the density of a specified phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Density in kg/m³.
     *
     * @code
     * double density = engine.phaseDensity(1, 0);
     * @endcode
     */
    auto phaseDensity(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the density of a specified phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Density in kg/m³.
     *
     * @code
     * double density = engine.phaseDensity(1, "aq_gen");
     * @endcode
     */
    auto phaseDensity(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the density of a specified phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Density in kg/m³.
     *
     * @code
     * double density = engine.phaseDensity_i(1, 0);
     * @endcode
     */
    auto phaseDensity_i(size_t ndx, size_t iphase) const -> double;


    /**
     * @brief Returns the masses of all phases for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of masses in kg.
     *
     * @code
     * auto masses = engine.phaseMasses(1);
     * @endcode
     */
    auto phaseMasses(size_t ndx) const -> Vector;

    /**
     * @brief Returns the mass of a specified phase for the node.
     *  Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Mass in kg.
     *
     * @code
     * double mass = engine.phaseMass(1, 0);
     * @endcode
     */
    auto phaseMass(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the mass of a specified phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Mass in kg.
     *
     * @code
     * double mass = engine.phaseMass(1, "aq_gen");
     * @endcode
     */
    auto phaseMass(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the mass of a specified phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Mass in kg.
     *
     * @code
     * double mass = engine.phaseMass_i(1, 0);
     * @endcode
     */
    auto phaseMass_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar amounts of all phases for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of phase amounts in mol.
     *
     * @code
     * auto phaseMol = engine.phaseAmounts(1);
     * @endcode
     */
    auto phaseAmounts(size_t ndx) const -> Vector;

    /**
     * @brief Returns the molar amount of a specific phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Amount in mol.
     *
     * @code
     * double amt = engine.phaseAmount(1, 0);
     * @endcode
     */
    auto phaseAmount(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the molar amount of a specific phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Amount in mol.
     *
     * @code
     * double amt = engine.phaseAmount(1, "aq_gen");
     * @endcode
     */
    auto phaseAmount(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar amount of a specific phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Amount in mol.
     *
     * @code
     * double amt = engine.phaseAmount_i(1, 0);
     * @endcode
     */
    auto phaseAmount_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the volumes of all phases  for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of volumes in m³.
     *
     * @code
     * auto volumes = engine.phaseVolumes(1);
     * @endcode
     */
    auto phaseVolumes(size_t ndx) const -> Vector;

    /**
     * @brief Returns the volume of a specific phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Volume in m³.
     *
     * @code
     * double vol = engine.phaseVolume(1, 0);
     * @endcode
     */
    auto phaseVolume(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the volume of a specific phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Volume in m³.
     *
     * @code
     * double vol = engine.phaseVolume(1, "aq_gen");
     * @endcode
     */
    auto phaseVolume(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the volume of a specific phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Volume in m³.
     *
     * @code
     * double vol = engine.phaseVolume_i(1, 0);
     * @endcode
     */
    auto phaseVolume_i(size_t ndx, size_t iphase) const -> double;
    

    /**
     * @brief Returns the enthalpies of all phases for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of enthalpies in J.
     *
     * @code
     * auto enth = engine.phaseEnthalpies(1);
     * @endcode
     */
    auto phaseEnthalpies(size_t ndx) const -> Vector;

    /**
     * @brief Returns the enthalpy of a specific phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Enthalpy in J.
     *
     * @code
     * double enth0 = engine.phaseEnthalpy(1, 0);
     * @endcode
     */
    auto phaseEnthalpy(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the enthalpy of a specific phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Enthalpy in J.
     *
     * @code
     * double enth0 = engine.phaseEnthalpy(1, "aq_gen");
     * @endcode
     */
    auto phaseEnthalpy(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the enthalpy of a specific phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Enthalpy in J.
     *
     * @code
     * double enth0 = engine.phaseEnthalpy_i(1, 0);
     * @endcode
     */
    auto phaseEnthalpy_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the entropies of all phases for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of entropies in J/K.
     *
     * @code
     * auto ents = engine.phaseEntropies(1);
     * @endcode
     */
    auto phaseEntropies(size_t ndx) const -> Vector;

    /**
     * @brief Returns the entropy of a specified phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Entropy in J/(K).
     *
     * @code
     * double ent0 = engine.phaseEntropy(1, 0);
     * @endcode
     */
    auto phaseEntropy(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the entropy of a specified phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Entropy in J/(K).
     *
     * @code
     * double ent0 = engine.phaseEntropy(1, "aq_gen");
     * @endcode
     */
    auto phaseEntropy(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the entropy of a specified phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Entropy in J/(K).
     *
     * @code
     * double ent0 = engine.phaseEntropy_i(1, 0);
     * @endcode
     */
    auto phaseEntropy_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the heat capacities (Cp) of all phases for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of heat capacities in J/K.
     *
     * @code
     * auto cps = engine.phaseHeatCapacitiesConstP(1);
     * @endcode
     */
    auto phaseHeatCapacitiesConstP(size_t ndx) const -> Vector;

    /**
     * @brief Returns the heat capacity of a specified phase for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Heat capacity in J/K.
     *
     * @code
     * double cp0 = engine.phaseHeatCapacityConstP(1, 0);
     * @endcode
     */
    auto phaseHeatCapacityConstP(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the heat capacity of a specified phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Heat capacity in J/K.
     *
     * @code
     * double cp0 = engine.phaseHeatCapacityConstP(1, "aq_gen");
     * @endcode
     */
    auto phaseHeatCapacityConstP(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the heat capacity of a specified phase for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Heat capacity in J/K.
     *
     * @code
     * double cp0 = engine.phaseHeatCapacityConstP_i(1, 0);
     * @endcode
     */
    auto phaseHeatCapacityConstP_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the saturation indices of all phases (log₁₀ units) for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (VectorConstRef) Vector of saturation indices.
     *
     * @code
     * auto satIdx = engine.phaseSatIndices(1);
     * @endcode
     */
    auto phaseSatIndices(size_t ndx) const -> Vector;

    /**
     * @brief Returns the saturation index of a specific phase (log₁₀ units) for the node.
     * Access specified element with bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Saturation index.
     *
     * @code
     * double sat0 = engine.phaseSatIndex(1, 0);
     * @endcode
     */
    auto phaseSatIndex(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the saturation index of a specific phase (log₁₀ units) for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param phase (std::string) Phase name.
     * @return (double) Saturation index.
     *
     * @code
     * double sat0 = engine.phaseSatIndex(1, "aq_gen");
     * @endcode
     */
    auto phaseSatIndex(size_t ndx, std::string phase) const -> double;

    /**
     * @brief Returns the saturation index of a specific phase (log₁₀ units) for the node.
     * Access specified element without bounds checking.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param iphase (size_t) Phase index.
     * @return (double) Saturation index.
     *
     * @code
     * double sat0 = engine.phaseSatIndex_i(1, 0);
     * @endcode
     */
    auto phaseSatIndex_i(size_t ndx, size_t iphase) const -> double;

    /**
     * @brief Returns the total mass of the system for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) System mass in kg.
     *
     * @code
     * double sysMass = engine.systemMass(1);
     * std::cout << "System mass: " << sysMass << " kg" << std::endl;
     * @endcode
     */
    auto systemMass(size_t ndx) const -> double;

    /**
     * @brief Returns the total volume of the system for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Volume in m³.
     *
     * @code
     * double sysVol = engine.systemVolume(1);
     * std::cout << "System volume: " << sysVol << " m³" << std::endl;
     * @endcode
     */
    auto systemVolume(size_t ndx) const -> double;

    /**
     * @brief Returns the ionic strength of the aqueous phasefor the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Ionic strength in molal.
     *
     * @code
     * double ionicStr = engine.ionicStrength(1);
     * @endcode
     */
    auto ionicStrength(size_t ndx) const -> double;

    /**
     * @brief Returns the pH of the aqueous phasefor the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) pH (in the activity scale (-log10 molal)).
     *
     * @code
     * double ph = engine.pH(1);
     * @endcode
     */
    auto pH(size_t ndx) const -> double;

    /**
     * @brief Returns the pe of the aqueous phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) pe (in the activity scale (-log10 molal)).
     *
     * @code
     * double pe = engine.pe(1);
     * @endcode
     */
    auto pe(size_t ndx) const -> double;

    /**
     * @brief Returns the Eh of the aqueous phase for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Eh (V).
     *
     * @code
     * double eh = engine.Eh(1);
     * @endcode
     */
    auto Eh(size_t ndx) const -> double;

    /**
     * @brief Returns the total Gibbs energy of the system for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Gibbs energy in J/mol.
     *
     * @code
     * double sysGibbs = engine.systemGibbsEnergy(1);
     * @endcode
     */
    auto systemGibbsEnergy(size_t ndx) const -> double;

    /**
     * @brief Returns the total enthalpy of the system for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Enthalpy in J.
     *
     * @code
     * double sysEnthalpy = engine.systemEnthalpy(1);
     * @endcode
     */
    auto systemEnthalpy(size_t ndx) const -> double;

    /**
     * @brief Returns the total entropy of the system for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Entropy in J/K.
     *
     * @code
     * double sysEntropy = engine.systemEntropy(1);
     * @endcode
     */
    auto systemEntropy(size_t ndx) const -> double;

    /**
     * @brief Returns the total isobaric heat capacity of the system for the node.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @return (double) Heat capacity in J/K.
     *
     * @code
     * double sysCp = engine.systemHeatCapacityConstP(1);
     * @endcode
     */
    auto systemHeatCapacityConstP(size_t ndx) const -> double;

    /**
     * @brief Returns the aqueous phase name.
     *
     * @return (string) aqueous phase name. If empty, the aqueous phase is not in system.
     *
     * @code
     * auto aqueous_name = engine.aqueousPhaseName();
     * @endcode
    */
    auto aqueousPhaseName() const -> std::string;

    /**
     * @brief Returns the gaseous phase name.
     *
     * @return (string) gaseous phase name. If empty, the gaseous phase is not in system.
     *
     * @code
     * auto gaseous_name = engine.gasPhaseName();
     * @endcode
    */
    auto gasPhaseName() const -> std::string;

    ChemicalEngineArray(const MassTransport& mt_class);
    std::shared_ptr<ChemicalEngineArray> getptr() {
        return shared_from_this();
    }

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;

    void setCorrectionbIC(size_t ndx);
    /**
     * @brief Reallocates internal arrays.
     *
     * Should be invoked after the system composition changes to ensure proper memory
     * configuration.
    */
    auto reallocateEngineArrays() -> void;

};

/**
   * @brief Outputs the state of the ChemicalEngineArray.
   *
   * Overloaded stream operator for convenient printing of engine details.
   *
   * @param out Output stream.
   * @param engine ChemicalEngineArray instance.
   * @return (std::ostream&) Updated output stream.
   *
   * @code
   * // Example: Print the engine state.
   * std::cout << engine << std::endl;
   * @endcode
   */
auto operator<<(std::ostream& out, const ChemicalEngineArray& engine) -> std::ostream&;

} // namespace xGEM2MT
