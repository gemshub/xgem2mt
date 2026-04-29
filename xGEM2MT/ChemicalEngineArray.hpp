/**
 * @file ChemicalEngineArray.hpp
 * @brief Header file for the ChemicalEngineArray class.
 *
 * xGEMS2MT is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
 * This file defines the classes and functions for Gibbs-energy minimization
 * thermodynamic modeling with xGEMS2MT. Units used throughout the API are:
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
};


/**
   * @class ChemicalEngineArray
   * @brief Class for equilibrium computations and thermodynamic analysis for the list of  systems.
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
     * Creates an empty ChemicalEngineArray instance. Use one of the initialization routines
     * (such as initialize() or initializeFromJsonStrings()) to configure the system.
     *
     * @param size (size_t) Number of nodes (systems) in the list.
     *
     * @code
     * // Example: Create an engine and then initialize it with a system file.
     * xGEMS2MT::ChemicalEngineArray engine(5);
     * std::vector<long int> nodeTypes={1,0,0,0,0};
     * std::string gems_init_lst = "nodearray-data/CalcColumn-dat.lst";
     * std::string gems_init_dbr_lst = "nodearray-data/CalcColumn-dbr.lst";
     * engine.initialize(gems_init_lst, gems_init_dbr_lst, ); // Temperature in K, Pressure in Pa.
     * @endcode
     */
    /// Construct a ChemicalEngineArray object for the 1D arrangement of nodes
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
     * xGEMS::ChemicalEngineOptions opt;
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
     * @code
     * // Example: Disable warm start for higher accuracy.
     * xGEMS::ChemicalEngineOptions opt = engine.options();
     * opt.warmstart = false;
     * engine.setOptions(opt);
     * @endcode
     *
     * @return const ChemicalEngineArrayOptions& The current engine options.
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
    /// Set smart start initial approximation (faster convergence, can be less accurate)
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
    /// Set cold start automatic initial approximation (slower convergence, more accurate)
    auto setColdStart() -> void;

    /// Get total number of nodes in the node array
    auto nNodes() const -> long int;

    /**
     * @brief Initializes the ChemicalEngine from a GEM-Selektor project file.
     *
     * Reads the file and configures the chemical system defined inside it.
     *
     * @param filename Path to the system definition file (e.g., "my-system-dat.lst").
     *
     * @code
     * // Example: Initialize using the system definition file.
     * xGEMS::ChemicalEngine engine;
     * engine.initialize("my-system-dat.lst");
     * @endcode
     *//// Initialize the ChemicalEngineArray object from a GEM-Selektor project file.
    /// @param filename The name of the file containing the definition of the chemical system
    auto initialize(std::string ipm_filename, std::string dbr_filename, std::vector<long int> nodeTypes = {}) -> void;

    /**
     * @brief Initializes the ChemicalEngine using JSON strings.
     *
     * Configures the engine from three JSON strings:
     * - dch_json: Defines the chemical system.
     * - ipm_json: Contains calculation parameters and settings.
     * - dbr_json: Specifies the input node composition.
     * These strings can be obtained from the GEM-Selektor exported GEMS3K system files as JSON (*-imp.json, *-dbc.json, *-dbr.json).
     *
     * @param dch_json JSON string describing the chemical system.
     * @param ipm_json JSON string with parameter and algorithm settings.
     * @param dbr_json JSON string containing node composition details.
     *
     * @code
     * // Example: Initialize using JSON-based input.
     * std::string dch = "{}";
     * std::string ipm = "{}";
     * std::string dbr = "{}";
     * engine.initializeFromJsonStrings(dch, ipm, dbr);
     * @endcode
     */
    /// Initialize the ChemicalEngineArray object from JSON strings for dch, ipm and dbr exported from GEM-Selektor \ .
    /// Normally to be used after ChemicalEngineArray();
    /// @param dch_json The json string containing the definition of the chemical system
    /// @param ipm_json The json string containing the parameters and settings for GEMS3K IPM-3 algorithm
    /// @param dbr_json The vector of json string containing the input node composition of the chemical system
    /// @param nodeTypes The initial node contents from DATABR files will be distributed among
    ///                   nodes in array according to the distribution index list nodeTypes
    auto initializeFromJsonStrings(std::string dch_json, std::string ipm_json,
                                   std::vector<std::string> dbr_json, std::vector<long int> nodeTypes = {}) -> void;

    /// Dump MULTI, DATACH and DATABR files structure.
    /// @param filepath - IPM work structure file path&name
    /// @param put_all -  if true print array of nodes for previous time point
    ///                           and for current time point
    auto dumpGEMS3K( const std::string& filepath, bool put_all = false ) const -> void;

    /// Copies array of nodes for current time point to the same place in the
    /// array of nodes for previous time point.
    auto nextTimePoint(bool always) -> void;

    /**
     * @brief Reads a DBR file from disk.
     *
     * Loads system composition and parameter data from a DBR file.
     *
     * @param ndx (size_t) Index of the node (system) in the list.
     * @param filename Path to the input DBR file (e.g., "*-dbr.json/.dat").
     *
     * @code
     * // Example: Read system composition from a DBR file.
     * engine.readDbrFromFile("my-system-dbr-0-0000.json");
     * @endcode
     */
    /// Reads another DBR file (with input system composition, T,P etc.) to a node with absolute index ndx.
    /// The DBR file must be compatible with the currently loaded IPM and DCH files
    /// (see description of initialize() function call).
    /// @param Null-terminated (C) string containing a full path to the input DBR disk file.
    auto readDbrFromFile(long int ndx, std::string filename) -> void;

    /**
     * @brief Reads System DBR composition from a JSON string.
     *
     * Updates the system composition using JSON information.
     *
     * @param dbr_json JSON string containing DBR composition data.
     *
     * @code
     * // Example: Update composition dynamically from a JSON string.
     * engine.readDbrFromJsonString("{}");
     * @endcode
     */
    /// Reads another DBR object (with input system composition, T,P etc.) from JSON string
    /// to a node with absolute index ndx. The DBR file must be compatible with the currently
    /// loaded IPM and DCH objects (see description of initialize() function call).
    /// @param Null-terminated (C) string containing a full path to the input DBR disk file.
    auto readDbrFromJsonString(long int  ndx, std::string dbr_json) -> void;

    /**
     * @brief Writes the current DBR to a file (key-value format).
     *
     * Saves the state of the system (after equilibrium calculations) into a DBR file.
     *
     * @param filename Path to the output DBR file (e.g., "my-system-dbr-0-0001.dat").
     *
     * @code
     * // Example: Save the current state to"my-system-dbr-0-0001.dat".
     * engine.writeDbrToFile("my-system-dbr-0-0001.dat");
     * @endcode
     */
    /// Writes a DBR file (normally after some changes via API and GEM calculation) from a node
    /// with absolute index ndx.
    /// @param Null-terminated (C) string containing a full path to the output DBR disk file.
    auto writeDbrToFile(long int ndx, std::string filename) -> void;


    /**
     * @brief Returns the current DBR as a JSON string.
     *
     * @return A JSON string representing the system composition.
     *
     * @code
     * // Example: Retrieve the system state as JSON.
     * std::string json_state = engine.writeDbrToJsonString();
     * std::cout << json_state << std::endl;
     * @endcode
     */
    /// Returns a DBR object (normally after some changes via API and GEM calculation) as JSON string
    /// from a node with absolute index ndx
    /// In case of error raises an exception or returns empty string.
    auto writeDbrToJsonString(long int  ndx) -> const std::string;


    /**
     * @brief Returns the number of elements in the system.
     *
     * @return (Index) Total number of elements.
     *
     * @code
     * // Example: Get the number of elements.
     * Index nElements = engine.numElements();
     * std::cout << "Elements: " << nElements << std::endl;
     * @endcode
     */
    auto numElements() const -> Index;

    /**
     * @brief Returns the number of species in the system.
     *
     * @return (Index) Total number of chemical species.
     *
     * @code
     * // Example: Display the number of species.
     * Index speciesCount = engine.numSpecies();
     * std::cout << "Species: " << speciesCount << std::endl;
     * @endcode
     */
    auto numSpecies() const -> Index;

    /**
     * @brief Returns the number of phases in the system.
     *
     * @return (Index) Total number of phases.
     *
     * @code
     * // Example: Print the number of phases.
     * Index phases = engine.numPhases();
     * std::cout << "Phases: " << phases << std::endl;
     * @endcode
     */
    auto numPhases() const -> Index;

    /**
     * @brief Returns the number of species in a given phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Index of the phase.
     * @return (Index) Number of species in that phase.
     *
     * @code
     * // Example: For phase with index 0.
     * Index sppInPhase = engine.numSpeciesInPhase(0);
     * std::cout << "Species in phase 0: " << sppInPhase << std::endl;
     * @endcode
     */
    auto numSpeciesInPhase(Index iphase) const -> Index;

    /**
     * @brief Returns the number of species in a phase identified by name.
     *
     * @param phase (std::string) Phases name.
     * @return (Index) Number of species in that phase.
     *
     * @code
     * // Example: For phase with index 0.
     * Index sppInPhase = engine.numSpeciesInPhase("aq_gen");
     * std::cout << "Species in phase aq_gen: " << sppInPhase << std::endl;
     * @endcode
     */
    auto numSpeciesInPhase(std::string phase) const -> Index;

    /**
     * @brief Returns the number of species in a given phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Index of the phase.
     * @return (Index) Number of species in that phase.
     *
     * @code
     * // Example: For phase with index 0.
     * Index sppInPhase = engine.numSpeciesInPhase_i(0);
     * std::cout << "Species in phase 0: " << sppInPhase << std::endl;
     * @endcode
     */
    auto numSpeciesInPhase_i(Index iphase) const -> Index;

    /**
     * @brief Returns the name of an element.
     * Access specified element with bounds checking.
     *
     * @param ielement (Index) Index of the element.
     * @return (std::string) Element name.
     *
     * @code
     * // Example: Get the name of the first element.
     * std::string element = engine.elementName(0);
     * std::cout << "Element 0: " << element << std::endl;
     * @endcode
     */
    auto elementName(Index ielement) const -> std::string;

    /**
     * @brief Returns the name of an element.
     * Access specified element without bounds checking.
     *
     * @param ielement (Index) Index of the element.
     * @return (std::string) Element name.
     *
     * @code
     * // Example: Get the name of the first element.
     * std::string element = engine.elementName_i(0);
     * std::cout << "Element 0: " << element << std::endl;
     * @endcode
     */
    auto elementName_i(Index ielement) const -> std::string;

    /**
     * @brief Returns the name of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @return (std::string) Species name.
     *
     * @code
     * // Example: Retrieve the name of species 0.
     * std::string spName = engine.speciesName(0);
     * std::cout << "Species 0: " << spName << std::endl;
     * @endcode
     */
    auto speciesName(Index ispecies) const -> std::string;

    /**
     * @brief Returns the name of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @return (std::string) Species name.
     *
     * @code
     * // Example: Retrieve the name of species 0.
     * std::string spName = engine.speciesName_i(0);
     * std::cout << "Species 0: " << spName << std::endl;
     * @endcode
     */
    auto speciesName_i(Index ispecies) const -> std::string;

    /**
     * @brief Returns the electrical charge of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @return (double) Charge of the species.
     *
     * @code
     * // Example: Get the charge of species with index 0.
     * double charge = engine.speciesCharge(0);
     * std::cout << "Charge: " << charge << std::endl;
     * @endcode
     */
    auto speciesCharge(Index ispecies) const -> double;

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
     * @param ispecies (Index) Index of the species.
     * @return (double) Charge of the species.
     *
     * @code
     * // Example: Get the charge of species with index 0.
     * double charge = engine.speciesCharge_i(0);
     * std::cout << "Charge: " << charge << std::endl;
     * @endcode
     */
    auto speciesCharge_i(Index ispecies) const -> double;

    /**
     * @brief Returns the name of the phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Index of the phase.
     * @return (std::string) Phase name.
     *
     * @code
     * // Example: Retrieve the name of phase 0.
     * std::string phase = engine.phaseName(0);
     * std::cout << "Phase 0: " << phase << std::endl;
     * @endcode
     */
    auto phaseName(Index iphase) const -> std::string;

    /**
     * @brief Returns the name of the phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Index of the phase.
     * @return (std::string) Phase name.
     *
     * @code
     * // Example: Retrieve the name of phase 0.
     * std::string phase = engine.phaseName_i(0);
     * std::cout << "Phase 0: " << phase << std::endl;
     * @endcode
     */
    auto phaseName_i(Index iphase) const -> std::string;

    /**
     * @brief Returns the index of an element by name.
     *
     * @param element (std::string) Element symbol or name.
     * @return (Index) Index of the element (or total elements if not found).
     *
     * @code
     * // Example: Get the index of Oxygen ("O").
     * Index idx = engine.indexElement("O");
     * std::cout << "Index of O: " << idx << std::endl;
     * @endcode
     */
    auto indexElement(std::string element) const -> Index;

    /**
     * @brief Returns the index of a species by name.
     *
     * @param species (std::string) Name of the species.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (Index) Index of the species (or total species if not found).
     *
     * @code
     * // Example: Get the index of water ("H2O@").
     * Index idx = engine.indexSpecies("H2O@");
     * std::cout << "Index of H2O@: " << idx << std::endl;
     * @endcode
     */
    auto indexSpecies(std::string species, std::optional<std::string> phase = std::nullopt) const -> Index;

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
     * @return (Index) Index of the phase.
     *
     * @code
     * // Example: Get index for the "aq_gen" phase.
     * Index aqueousIndex = engine.indexPhase("aq_gen");
     * @endcode
     */
    auto indexPhase(std::string phase) const -> Index;

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
     * @param ispecies (Index) Index of the species.
     * @return (Index) Index of the phase that contains the species.
     *
     * @code
     * // Example: Find the phase of the species with index 0.
     * Index phaseIdx = engine.indexPhaseWithSpecies(0);
     * @endcode
     */
    auto indexPhaseWithSpecies(Index ispecies) const -> Index;

    /**
     * @brief Returns the index of the first species in a specified phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Index of the phase.
     * @return (Index) Index of the first species in that phase or numSpecies if illegal iphase.
     *
     * @code
     * // Example: Get the first species index in phase 0.
     * Index firstSpec = engine.indexFirstSpeciesInPhase(0);
     * @endcode
     */
    auto indexFirstSpeciesInPhase(Index iphase) const -> Index;

    /**
     * @brief Returns the index of the first species in a specified by name phase.
     *
     * @param phase (std::string) Phase name.
     * @return (Index) Index of the first species in that phase or numSpecies if illegal iphase.
     *
     * @code
     * // Example: Get the first species index in phase aq_gen.
     * Index firstSpec = engine.indexFirstSpeciesInPhase("aq_gen");
     * @endcode
     */
    auto indexFirstSpeciesInPhase(std::string phase) const -> Index;

    /**
     * @brief Returns the index of the first species in a specified phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Index of the phase.
     * @return (Index) Index of the first species in that phase or numSpecies if illegal iphase.
     *
     * @code
     * // Example: Get the first species index in phase 0.
     * Index firstSpec = engine.indexFirstSpeciesInPhase_i(0);
     * @endcode
     */
    auto indexFirstSpeciesInPhase_i(Index iphase) const -> Index;

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
     * @brief Sets the complete species amounts vector.
     *
     * @param n (VectorConstRef) Vector of species amounts (mol).
     *
     * @code
     * // Example: Set all species amounts to 1.0 mol.
     * Eigen::VectorXd n(engine.numSpecies());
     * n.setOnes();
     * engine.setSpeciesAmounts(n);
     * @endcode
     */
    auto setSpeciesAmounts(long int ndx, VectorConstRef n) -> void;

    /**
     * @brief Sets the amount for a species identified by name (as in the GEMS system).
     *
     * @param name (std::string) Species name.
     * @param amount (double) New amount in mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Set the amount of "CaSO4@" to 0.01 mol.
     * engine.setSpeciesAmount("CaSO4@", 0.01);
     * @endcode
     */
    auto setSpeciesAmount(long int ndx, std::string name, double amount, std::optional<std::string> phase = std::nullopt) -> void;

    /**
     * @brief Sets the amount for a species identified by its index.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @param amount (double) New amount in mol.
     *
     * @code
     * // Example: Set species at index 0 to 0.01 mol.
     * engine.setSpeciesAmount(0, 0.01);
     * @endcode
     */
    auto setSpeciesAmount(long int ndx, Index ispecies, double amount) -> void;

    /**
     * @brief Sets the amount for a species identified by its index.
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @param amount (double) New amount in mol.
     *
     * @code
     * // Example: Set species at index 0 to 0.01 mol.
     * engine.setSpeciesAmount_i(0, 0.01);
     * @endcode
     */
    auto setSpeciesAmount_i(long int ndx, Index ispecies, double amount) -> void;

    /**
     * @brief Sets an upper bound for a species identified by name.
     *
     * If amount < 0, resets to default 1e6.
     *
     * @param name (std::string) Species name.
     * @param amount (double) Upper limit in mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Set the upper limit for "SiO2" to 0.1 mol.
     * engine.setSpeciesUpperLimit("SiO2", 0.1);
     * @endcode
     */
    auto setSpeciesUpperLimit(long int ndx, std::string name, double amount, std::optional<std::string> phase = std::nullopt) -> void;

    /**
     * @brief Sets a lower bound for a species identified by name.
     *
     * If amount < 0, resets to default 0.
     *
     * @param name (std::string) Species name.
     * @param amount (double) Lower limit in mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Set the lower limit for "SiO2" to 0.05 mol.
     * engine.setSpeciesLowerLimit("SiO2", 0.05);
     * @endcode
     */
    auto setSpeciesLowerLimit(long int ndx, std::string name, double amount, std::optional<std::string> phase = std::nullopt) -> void;

    /**
     * @brief Sets an upper bound (maximum amount allowed to form) for a species identified by its index.
     * Access specified element with bounds checking.
     *
     * If amount < 0, resets to default 1e6.
     *
     * @param ispecies (Index) Species index.
     * @param amount (double) Upper limit in mol.
     *
     * @code
     * // Example: Set upper limit for species index 0.
     * engine.setSpeciesUpperLimit(0, 0.1);
     * @endcode
     */
    auto setSpeciesUpperLimit(long int ndx, Index ispecies, double amount) -> void;

    /**
     * @brief Sets a lower bound (minimum amount allowed to form) for a species identified by its index.
     * Access specified element with bounds checking.
     *
     * If amount < 0, resets to default 0.
     *
     * @param ispecies (Index) Species index.
     * @param amount (double) Lower limit in mol.
     *
     * @code
     * // Example: Set lower limit for species index 0.
     * engine.setSpeciesLowerLimit(0, 0.05);
     * @endcode
     */
    auto setSpeciesLowerLimit(long int ndx, Index ispecies, double amount) -> void;

    /**
     * @brief Sets an upper bound (maximum amount allowed to form) for a species identified by its index.
     * Access specified element without bounds checking.
     *
     * If amount < 0, resets to default 1e6.
     *
     * @param ispecies (Index) Species index.
     * @param amount (double) Upper limit in mol.
     *
     * @code
     * // Example: Set upper limit for species index 0.
     * engine.setSpeciesUpperLimit_i(0, 0.1);
     * @endcode
     */
    auto setSpeciesUpperLimit_i(long int ndx, Index ispecies, double amount) -> void;

    /**
     * @brief Sets a lower bound (minimum amount allowed to form) for a species identified by its index.
     * Access specified element without bounds checking.
     *
     * If amount < 0, resets to default 0.
     *
     * @param ispecies (Index) Species index.
     * @param amount (double) Lower limit in mol.
     *
     * @code
     * // Example: Set lower limit for species index 0.
     * engine.setSpeciesLowerLimit_i(0, 0.05);
     * @endcode
     */
    auto setSpeciesLowerLimit_i(long int ndx, Index ispecies, double amount) -> void;

    /**
     * @brief Sets the standard molar Gibbs energy for a species (@ T, P of the system).
     *
     * @param name (std::string) Species name.
     * @param value (double) Standard molar Gibbs energy in J/mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Set the standard molar Gibbs energy of H2O.
     * engine.setStandardMolarGibbsEnergy("H2O", -237140); // Value in J/mol.
     * @endcode
     */
    auto setStandardMolarGibbsEnergy(std::string name, double TK, double P, double value, std::optional<std::string> phase = std::nullopt) -> void;

    /**
     * @brief Sets the standard molar Gibbs energy for a species (@ T, P of the system).
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @param value (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: Set the standard molar Gibbs energy of H2O.
     * engine.setStandardMolarGibbsEnergy(0, -237.13); // Value in J/mol.
     * @endcode
     */
    auto setStandardMolarGibbsEnergy(Index ispecies, double TK, double P, double value) -> void;

    /**
     * @brief Sets the standard molar Gibbs energy for a species (@ T, P of the system).
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @param value (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: Set the standard molar Gibbs energy of H2O.
     * engine.setStandardMolarGibbsEnergy_i(0, -237.13); // Value in J/mol.
     * @endcode
     */
    auto setStandardMolarGibbsEnergy_i(Index ispecies, double TK, double P, double value) -> void;

    /**
     * @brief Sets upper limits for all species.
     *
     * @param n (VectorConstRef) Vector of upper limits in mol.
     *
     * @code
     * // Example: Set all species' upper limits to values provided in vector "limits".
     * engine.setSpeciesUpperLimits(limits);
     * @endcode
     */
    auto setSpeciesUpperLimits(long int ndx,VectorConstRef n) -> void;

    /**
     * @brief Sets lower limits for all species.
     *
     * @param n (VectorConstRef) Vector of lower limits in mol.
     *
     * @code
     * // Example: Set all species' lower limits.
     * engine.setSpeciesLowerLimits(lowerLimits);
     * @endcode
     */
    auto setSpeciesLowerLimits(long int ndx,VectorConstRef n) -> void;

    /**
     * @brief Sets the pressure and temperature without computing equilibrium.
     *
     * @param P Pressure in Pascals (Pa).
     * @param T Temperature in Kelvin (K).
     * @return (bool) true if PT was set correctly, false if out of range.
     *
     * @code
     * // Example: Set PT to 101325 Pa and 298.15 K.
     * bool no_error = engine.setPT(101325, 298.15);
     * std::cout << (no_error ? "PT set correctly" : "PT error") << std::endl;
     * @endcode
     */
    auto setPT(long int ndx, double P, double T) const -> bool;

    /**
     * @brief Sets the amounts of elements (vector b) without computing equilibrium.
     *
     * @param b (VectorConstRef) Vector of element amounts in mol (order based on the index of elements in the chemical system).
     *
     * @code
     * // Example: Set elements amounts of a solution with 0.001 M MgCl2 0.001 M CO2
     * // elements order: C, Ca, Cl, H, Mg, O, Zz (charge)
     * Eigen::VectorXd b(7);
     * b << 0.001, 1e-09, 0.004, 110.6837, 0.002, 55.34405, 0;  // in mol
     * engine.setB(b);
     * @endcode
     */
    auto setB(long int ndx, VectorConstRef b) -> void;

    /**
     * @brief Re-equilibrates the system with (or without) a warm start.
     *
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
    auto reequilibrate(long int  ndx, bool warmstart) -> long int;

    /**
     * @brief Re-equilibrates the system with (or without) a warm start.
     *
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
    /// Calls GEM IPM calculation for for all nodes of TNodeArray
    /// in a loop with an optional openmp parallelization
    /// If @param warmstart is true then GEMS3K uses contents of node as initial
    ///   guess, otherwise (false) it gets the simplex LP initial guess.
    auto reequilibrate(bool warmstart) -> void;

    /**
     * @brief Re-equilibrates the system with (or without) a warm start.
     *
     * @return (int) Return code of the re-equilibration.
     *
     * @code
     * // Example: Re-equilibrate.
     * int status = engine.reequilibrate();
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
    auto equilibrate(long int  ndx) -> bool;

    /**
     * @brief Re-equilibrates the system with (or without) a warm start.
     *
     * @return (int) Return code of the re-equilibration.
     *
     * @code
     * // Example: Re-equilibrate.
     * int status = engine.reequilibrate();
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
    /// Recalculate the equilibrium state with preset P,T,B and other inputs for all nodes
    /// Here we compare this node for current time and for previous time
    /// to test need to recalculate equilibrium in this node because its vector b has changed
    /// Returns true if all GEM Run finished with OK_GEM_AIA or OK_GEM_SIA
    auto equilibrate() -> bool;

    /**
     * @brief Checks if the equilibrium calculation converged.
     *
     * @return (bool) True if converged, false otherwise.
     *
     * @code
     * // Example: Check convergence status.
     * bool isConverged = engine.converged();
     * @endcode
     */
    auto converged(long int ndx) const -> bool;

    /**
     * @brief Returns the number of iterations performed during equilibrium.
     *
     * @return (Index) The number of iterations.
     *
     * @code
     * // Example: Display the number of iterations.
     * Index iters = engine.numIterations();
     * std::cout << "Iterations: " << iters << std::endl;
     * @endcode
     */
    auto numIterations(long int ndx) const -> Index;

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
     * @brief Returns the current temperature of the system.
     *
     * @return (double) Temperature in Kelvin (K).
     *
     * @code
     * // Example: Get the system temperature.
     * double temp = engine.temperature();
     * std::cout << "Temperature: " << temp << " K" << std::endl;
     * @endcode
     */
    auto temperature(long int ndx) const -> double;

    /**
     * @brief Returns the current pressure of the system.
     *
     * @return (double) Pressure in Pascals (Pa).
     *
     * @code
     * // Example: Retrieve the system pressure.
     * double press = engine.pressure();
     * std::cout << "Pressure: " << press << " Pa" << std::endl;
     * @endcode
     */
    auto pressure(long int ndx) const -> double;

    /**
     * @brief Returns the amounts of the elements.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @return (VectorConstRef) Vector of element amounts in mol.
     *
     * @code
     * // Example: Get element amounts.
     * auto eAmounts = engine.elementAmounts();
     * @endcode
     */
    auto elementAmounts(long int ndx) const -> Vector;

    /**
     * @brief Returns the amounts of elements in a specified phase.
     * Access specified element with bounds checking.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param iphase Index of the phase.
     * @return A vector containing the amounts of each element in the phase (in mol).
     *
     * **Example usage:**
     * @code
     * // Retrieve the element amounts in phase index 0.
     * Eigen::VectorXd elementAmounts = engine.elementAmountsInPhase(0);
     * std::cout << "Element amounts in phase 0: " << elementAmounts.transpose() << std::endl;
     * @endcode
     */
    auto elementAmountsInPhase(long int ndx, Index iphase) const -> Vector;

    /**
     * @brief Returns the amounts of elements in a specified phase.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param phase Name of the phase.
     * @return A vector containing the amounts of each element in the phase (in mol).
     *
     * **Example usage:**
     * @code
     * // Retrieve the element amounts in phase aq_gen.
     * Eigen::VectorXd elementAmounts = engine.elementAmountsInPhase("aq_gen");
     * std::cout << "Element amounts in phase aq_gen: " << elementAmounts.transpose() << std::endl;
     * @endcode
     */
    auto elementAmountsInPhase(long int ndx,std::string phase) const -> Vector;

    /**
     * @brief Returns the amounts of elements in a specified phase.
     * Access specified element without bounds checking.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param iphase Index of the phase.
     * @return A vector containing the amounts of each element in the phase (in mol).
     *
     * **Example usage:**
     * @code
     * // Retrieve the element amounts in phase index 0.
     * Eigen::VectorXd elementAmounts = engine.elementAmountsInPhase_i(0);
     * std::cout << "Element amounts in phase 0: " << elementAmounts.transpose() << std::endl;
     * @endcode
     */
    auto elementAmountsInPhase_i(long int ndx,Index iphase) const -> Vector;

    /**
     * @brief Returns the amounts of elements in a specified group of species.
     *
     * The order of the elements in the vector is based on the index of elements in the chemical system.
     *
     * @param ispecies A vector containing the indices of the species whose elemental amounts should be retrieved.
     * @return A vector containing the amounts of each element corresponding to the selected species (in mol).
     *
     * **Example usage:**
     * @code
     * // Retrieve element amounts contributed by species indices {0, 2, 5}.
     * Eigen::VectorXi speciesIndices(3);
     * speciesIndices << 0, 2, 5;
     * Eigen::VectorXd elementAmounts = engine.elementAmountsInSpecies(speciesIndices);
     * std::cout << "Element amounts from selected species: " << elementAmounts.transpose() << std::endl;
     * @endcode
     */
    auto elementAmountsInSpecies(long int ndx, VectorXiConstRef ispecies) const -> Vector;

    /**
     * @brief Returns the amounts of the species.
     *
     * @return (VectorConstRef) Vector of species amounts in mol.
     *
     * @code
     * // Example: Retrieve the species amounts.
     * auto spAmounts = engine.speciesAmounts();
     * @endcode
     */
    auto speciesAmounts(long int ndx) const -> Vector;

    /**
     * @brief Returns the amount of a species by index.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @return (double) Amount in mol.
     *
     * @code
     * // Example: Get the species amount at index 0.
     * double amount = engine.speciesAmount(0);
     * @endcode
     */
    auto speciesAmount(long int ndx, Index ispecies) const -> double;

    /**
     * @brief Returns the amount of a species identified by name.
     *
     * @param name (std::string) Species name.
     * @return (double) Amount in mol.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     *
     * @code
     * // Example: Get the amount of "OH-".
     * double waterAmt = engine.speciesAmount("OH-");
     * @endcode
     */
    auto speciesAmount(long int ndx, std::string name, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the amount of a species by index.
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Index of the species.
     * @return (double) Amount in mol.
     *
     * @code
     * // Example: Get the species amount at index 0.
     * double amount = engine.speciesAmount_i(0);
     * @endcode
     */
    auto speciesAmount_i(long int ndx, Index ispecies) const -> double;

    /**
     * @brief Returns the upper limits for all species.
     *
     * @return (VectorConstRef) Vector of upper limits in mol.
     *
     * @code
     * // Example: Fetch all species upper limits.
     * auto upperLimits = engine.speciesUpperLimits();
     * @endcode
     */
    auto speciesUpperLimits(long int ndx) const -> Vector;

    /**
     * @brief Returns the lower limits for all species.
     *
     * @return (VectorConstRef) Vector of lower limits in mol.
     *
     * @code
     * // Example: Get the species lower limits.
     * auto lowerLimits = engine.speciesLowerLimits();
     * @endcode
     */
    auto speciesLowerLimits(long int ndx) const -> Vector;

    /**
     * @brief Returns the molalities of the species.
     *
     * Assumes the aqueous phase is the first phase.
     *
     * @return (VectorConstRef) Species molalities (mol/kg H2O@).
     *
     * @code
     * // Example: Retrieve species molalities.
     * auto molalities = engine.speciesMolalities();
     * @endcode
     */
    auto speciesMolalities(long int ndx) const -> VectorConstRef;

    /**
     * @brief Returns the mole fractions of the species.
     *
     * @return (VectorConstRef) Mole fractions.
     *
     * @code
     * // Example: Get species mole fractions.
     * auto moleFrac = engine.moleFractions();
     * @endcode
     */
    auto moleFractions(long int ndx) const -> Vector;

    /**
     * @brief Returns the ln activity coefficients of the species (mole fraction scale).
     *
     * @return (VectorConstRef) ln Activity coefficients.
     *
     * @code
     * // Example: Retrieve ln activity coefficients.
     * auto lnActCoeff = engine.lnActivityCoefficients();
     * @endcode
     */
    auto lnActivityCoefficients(long int ndx) const -> Vector;

    /**
     * @brief Returns the ln activities of the species.
     *
     * @return (VectorConstRef) ln Activities.
     *
     * @code
     * // Example: Fetch ln activities.
     * auto lnActivities = engine.lnActivities();
     * @endcode
     */
    auto lnActivities(long int ndx) const -> Vector;

    /**
     * @brief Returns the ln of species concentrations.
     *
     * @return (VectorConstRef) ln concentrations.
     *
     * @code
     * // Example: Get ln concentrations.
     * auto lnConcs = engine.lnConcentrations();
     * @endcode
     */
    auto lnConcentrations(long int ndx) const -> Vector;

    /**
     * @brief Returns the chemical potentials of the species.
     *
     * @return (VectorConstRef) Chemical potentials in J/mol.
     *
     * @code
     * // Example: Retrieve species chemical potentials.
     * auto chemPot = engine.chemicalPotentials();
     * @endcode
     */
    auto chemicalPotentials(long int ndx) const -> Vector;

    /**
     * @brief Returns the standard molar Gibbs energy of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: For species 0.
     * double stdGibbs = engine.standardMolarGibbsEnergy(0);
     * @endcode
     */
    auto standardMolarGibbsEnergy(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar Gibbs energy of a species.
     *
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: For species 0.
     * double stdGibbs = engine.standardMolarGibbsEnergy("SiO2");
     * @endcode
     */
    auto standardMolarGibbsEnergy(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar Gibbs energy of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Standard molar Gibbs energy in J/mol.
     *
     * @code
     * // Example: For species 0.
     * double stdGibbs = engine.standardMolarGibbsEnergy_i(0);
     * @endcode
     */
    auto standardMolarGibbsEnergy_i(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar enthalpy of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Standard molar enthalpy in J/mol.
     *
     * @code
     * double stdEnthalpy = engine.standardMolarEnthalpy(0);
     * @endcode
     */
    auto standardMolarEnthalpy(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar enthalpy of a species.
     *
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Standard molar enthalpy in J/mol.
     *
     * @code
     * double stdEnthalpy = engine.standardMolarEnthalpy("SiO2");
     * @endcode
     */
    auto standardMolarEnthalpy(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar enthalpy of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Standard molar enthalpy in J/mol.
     *
     * @code
     * double stdEnthalpy = engine.standardMolarEnthalpy_i(0);
     * @endcode
     */
    auto standardMolarEnthalpy_i(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar volume of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Standard molar volume in m³/mol.
     *
     * @code
     * double stdVolume = engine.standardMolarVolume(0);
     * @endcode
     */
    auto standardMolarVolume(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar volume of a species.
     *
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Standard molar volume in m³/mol.
     *
     * @code
     * double stdVolume = engine.standardMolarVolume("SiO2");
     * @endcode
     */
    auto standardMolarVolume(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar volume of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Standard molar volume in m³/mol.
     *
     * @code
     * double stdVolume = engine.standardMolarVolume_i(0);
     * @endcode
     */
    auto standardMolarVolume_i(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar entropy of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Standard molar entropy in J/(mol·K).
     *
     * @code
     * double stdEntropy = engine.standardMolarEntropy(0);
     * @endcode
     */
    auto standardMolarEntropy(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar entropy of a species.
     *
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Standard molar entropy in J/(mol·K).
     *
     * @code
     * double stdEntropy = engine.standardMolarEntropy("SiO2");
     * @endcode
     */
    auto standardMolarEntropy(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar entropy of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Standard molar entropy in J/(mol·K).
     *
     * @code
     * double stdEntropy = engine.standardMolarEntropy_i(0);
     * @endcode
     */
    auto standardMolarEntropy_i(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar isobaric heat capacity of a species.
     * Access specified element with bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double heatCap = engine.standardMolarHeatCapacityConstP(0);
     * @endcode
     */
    auto standardMolarHeatCapacityConstP(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the standard molar isobaric heat capacity of a species.
     *
     * @param species (std::string) Species name.
     * @param optional phase (std::string) Name of the phase the species was included in. If default get the first index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double heatCap = engine.standardMolarHeatCapacityConstP("SiO2");
     * @endcode
     */
    auto standardMolarHeatCapacityConstP(double TK, double P, std::string species, std::optional<std::string> phase = std::nullopt) const -> double;

    /**
     * @brief Returns the standard molar isobaric heat capacity of a species.
     * Access specified element without bounds checking.
     *
     * @param ispecies (Index) Species index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double heatCap = engine.standardMolarHeatCapacityConstP_i(0);
     * @endcode
     */
    auto standardMolarHeatCapacityConstP_i(double TK, double P, Index ispecies) const -> double;

    /**
     * @brief Returns the molar Gibbs energy of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Molar Gibbs energy in J/mol.
     *
     * @code
     * double phaseGibbs = engine.phaseMolarGibbsEnergy(0);
     * @endcode
     */
    auto phaseMolarGibbsEnergy(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar Gibbs energy of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Molar Gibbs energy in J/mol.
     *
     * @code
     * double phaseGibbs = engine.phaseMolarGibbsEnergy("aq_gen");
     * @endcode
     */
    auto phaseMolarGibbsEnergy(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar Gibbs energy of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Molar Gibbs energy in J/mol.
     *
     * @code
     * double phaseGibbs = engine.phaseMolarGibbsEnergy_i(0);
     * @endcode
     */
    auto phaseMolarGibbsEnergy_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar enthalpy of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Molar enthalpy in J/mol.
     *
     * @code
     * double phaseEnthalpy = engine.phaseMolarEnthalpy(0);
     * @endcode
     */
    auto phaseMolarEnthalpy(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar enthalpy of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Molar enthalpy in J/mol.
     *
     * @code
     * double phaseEnthalpy = engine.phaseMolarEnthalpy("aq_gen");
     * @endcode
     */
    auto phaseMolarEnthalpy(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar enthalpy of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Molar enthalpy in J/mol.
     *
     * @code
     * double phaseEnthalpy = engine.phaseMolarEnthalpy_i(0);
     * @endcode
     */
    auto phaseMolarEnthalpy_i(Ilong int ndx, ndex iphase) const -> double;

    /**
     * @brief Returns the molar volume of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Molar volume in m³/mol.
     *
     * @code
     * double phaseVol = engine.phaseMolarVolume(0);
     * @endcode
     */
    auto phaseMolarVolume(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar volume of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Molar volume in m³/mol.
     *
     * @code
     * double phaseVol = engine.phaseMolarVolume("aq_gen");
     * @endcode
     */
    auto phaseMolarVolume(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar volume of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Molar volume in m³/mol.
     *
     * @code
     * double phaseVol = engine.phaseMolarVolume_i(0);
     * @endcode
     */
    auto phaseMolarVolume_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar entropy of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Molar entropy in J/(mol·K).
     *
     * @code
     * double phaseEnt = engine.phaseMolarEntropy(0);
     * @endcode
     */
    auto phaseMolarEntropy(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar entropy of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Molar entropy in J/(mol·K).
     *
     * @code
     * double phaseEnt = engine.phaseMolarEntropy("aq_gen");
     * @endcode
     */
    auto phaseMolarEntropy(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar entropy of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Molar entropy in J/(mol·K).
     *
     * @code
     * double phaseEnt = engine.phaseMolarEntropy_i(0);
     * @endcode
     */
    auto phaseMolarEntropy_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar isobaric heat capacity of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double phaseCp = engine.phaseMolarHeatCapacityConstP(0);
     * @endcode
     */
    auto phaseMolarHeatCapacityConstP(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar isobaric heat capacity of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double phaseCp = engine.phaseMolarHeatCapacityConstP("aq_gen");
     * @endcode
     */
    auto phaseMolarHeatCapacityConstP(long int ndx, std::string phase) const -> double;


    /**
     * @brief Returns the molar isobaric heat capacity of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Heat capacity in J/(mol·K).
     *
     * @code
     * double phaseCp = engine.phaseMolarHeatCapacityConstP_i(0);
     * @endcode
     */
    auto phaseMolarHeatCapacityConstP_i(long int ndx, Index iphase) const -> double;


    /**
     * @brief Returns the specific Gibbs energy of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Specific Gibbs energy in J/kg.
     *
     * @code
     * double specGibbs = engine.phaseSpecificGibbsEnergy(0);
     * @endcode
     */
    auto phaseSpecificGibbsEnergy(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific Gibbs energy of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Specific Gibbs energy in J/kg.
     *
     * @code
     * double specGibbs = engine.phaseSpecificGibbsEnergy("aq_gen");
     * @endcode
     */
    auto phaseSpecificGibbsEnergy(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific Gibbs energy of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Specific Gibbs energy in J/kg.
     *
     * @code
     * double specGibbs = engine.phaseSpecificGibbsEnergy_i(0);
     * @endcode
     */
    auto phaseSpecificGibbsEnergy_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific enthalpy of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Specific enthalpy in J/kg.
     *
     * @code
     * double specEnthalpy = engine.phaseSpecificEnthalpy(0);
     * @endcode
     */
    auto phaseSpecificEnthalpy(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific enthalpy of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Specific enthalpy in J/kg.
     *
     * @code
     * double specEnthalpy = engine.phaseSpecificEnthalpy("aq_gen");
     * @endcode
     */
    auto phaseSpecificEnthalpy(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific enthalpy of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Specific enthalpy in J/kg.
     *
     * @code
     * double specEnthalpy = engine.phaseSpecificEnthalpy_i(0);
     * @endcode
     */
    auto phaseSpecificEnthalpy_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific volume of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Specific volume in m³/kg.
     *
     * @code
     * double specVolume = engine.phaseSpecificVolume(0);
     * @endcode
     */
    auto phaseSpecificVolume(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific volume of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Specific volume in m³/kg.
     *
     * @code
     * double specVolume = engine.phaseSpecificVolume("aq_gen");
     * @endcode
     */
    auto phaseSpecificVolume(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific volume of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Specific volume in m³/kg.
     *
     * @code
     * double specVolume = engine.phaseSpecificVolume_i(0);
     * @endcode
     */
    auto phaseSpecificVolume_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific entropy of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Specific entropy in J/(kg·K).
     *
     * @code
     * double specEntropy = engine.phaseSpecificEntropy(0);
     * @endcode
     */
    auto phaseSpecificEntropy(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific entropy of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Specific entropy in J/(kg·K).
     *
     * @code
     * double specEntropy = engine.phaseSpecificEntropy("aq_gen");
     * @endcode
     */
    auto phaseSpecificEntropy(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific entropy of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Specific entropy in J/(kg·K).
     *
     * @code
     * double specEntropy = engine.phaseSpecificEntropy_i(0);
     * @endcode
     */
    auto phaseSpecificEntropy_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific isobaric heat capacity of a phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Heat capacity in J/(kg·K).
     *
     * @code
     * double specCp = engine.phaseSpecificHeatCapacityConstP(0);
     * @endcode
     */
    auto phaseSpecificHeatCapacityConstP(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the specific isobaric heat capacity of a phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Heat capacity in J/(kg·K).
     *
     * @code
     * double specCp = engine.phaseSpecificHeatCapacityConstP("aq_gen");
     * @endcode
     */
    auto phaseSpecificHeatCapacityConstP(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the specific isobaric heat capacity of a phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Heat capacity in J/(kg·K).
     *
     * @code
     * double specCp = engine.phaseSpecificHeatCapacityConstP_i(0);
     * @endcode
     */
    auto phaseSpecificHeatCapacityConstP_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the densities of all phases.
     *
     * @return (VectorConstRef) Vector of densities in kg/m³.
     *
     * @code
     * auto densities = engine.phaseDensities();
     * @endcode
     */
    auto phaseDensities(long int ndx) const -> Vector;

    /**
     * @brief Returns the density of a specified phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Density in kg/m³.
     *
     * @code
     * double density = engine.phaseDensity(0);
     * @endcode
     */
    auto phaseDensity(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the density of a specified phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Density in kg/m³.
     *
     * @code
     * double density = engine.phaseDensity("aq_gen");
     * @endcode
     */
    auto phaseDensity(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the density of a specified phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Density in kg/m³.
     *
     * @code
     * double density = engine.phaseDensity_i(0);
     * @endcode
     */
    auto phaseDensity_i(long int ndx, Index iphase) const -> double;


    /**
     * @brief Returns the masses of all phases.
     *
     * @return (VectorConstRef) Vector of masses in kg.
     *
     * @code
     * auto masses = engine.phaseMasses();
     * @endcode
     */
    auto phaseMasses(long int ndx) const -> Vector;

    /**
     * @brief Returns the mass of a specified phase.
     *  Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Mass in kg.
     *
     * @code
     * double mass = engine.phaseMass(0);
     * @endcode
     */
    auto phaseMass(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the mass of a specified phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Mass in kg.
     *
     * @code
     * double mass = engine.phaseMass("aq_gen");
     * @endcode
     */
    auto phaseMass(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the mass of a specified phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Mass in kg.
     *
     * @code
     * double mass = engine.phaseMass_i(0);
     * @endcode
     */
    auto phaseMass_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar amounts of all phases.
     *
     * @return (VectorConstRef) Vector of phase amounts in mol.
     *
     * @code
     * auto phaseMol = engine.phaseAmounts();
     * @endcode
     */
    auto phaseAmounts(long int ndx) const -> Vector;

    /**
     * @brief Returns the molar amount of a specific phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Amount in mol.
     *
     * @code
     * double amt = engine.phaseAmount(0);
     * @endcode
     */
    auto phaseAmount(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the molar amount of a specific phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Amount in mol.
     *
     * @code
     * double amt = engine.phaseAmount("aq_gen");
     * @endcode
     */
    auto phaseAmount(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the molar amount of a specific phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Amount in mol.
     *
     * @code
     * double amt = engine.phaseAmount_i(0);
     * @endcode
     */
    auto phaseAmount_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the volumes of all phases.
     *
     * @return (VectorConstRef) Vector of volumes in m³.
     *
     * @code
     * auto volumes = engine.phaseVolumes();
     * @endcode
     */
    auto phaseVolumes(long int ndx) const -> Vector;

    /**
     * @brief Returns the volume of a specific phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Volume in m³.
     *
     * @code
     * double vol = engine.phaseVolume(0);
     * @endcode
     */
    auto phaseVolume(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the volume of a specific phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Volume in m³.
     *
     * @code
     * double vol = engine.phaseVolume("aq_gen");
     * @endcode
     */
    auto phaseVolume(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the volume of a specific phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Volume in m³.
     *
     * @code
     * double vol = engine.phaseVolume_i(0);
     * @endcode
     */
    auto phaseVolume_i(long int ndx, Index iphase) const -> double;
    

    /**
     * @brief Returns the enthalpies of all phases.
     *
     * @return (VectorConstRef) Vector of enthalpies in J.
     *
     * @code
     * auto enth = engine.phaseEnthalpies();
     * @endcode
     */
    auto phaseEnthalpies(long int ndx) const -> Vector;

    /**
     * @brief Returns the enthalpy of a specific phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Enthalpy in J.
     *
     * @code
     * double enth0 = engine.phaseEnthalpy(0);
     * @endcode
     */
    auto phaseEnthalpy(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the enthalpy of a specific phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Enthalpy in J.
     *
     * @code
     * double enth0 = engine.phaseEnthalpy("aq_gen");
     * @endcode
     */
    auto phaseEnthalpy(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the enthalpy of a specific phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Enthalpy in J.
     *
     * @code
     * double enth0 = engine.phaseEnthalpy_i(0);
     * @endcode
     */
    auto phaseEnthalpy_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the entropies of all phases.
     *
     * @return (VectorConstRef) Vector of entropies in J/K.
     *
     * @code
     * auto ents = engine.phaseEntropies();
     * @endcode
     */
    auto phaseEntropies(long int ndx) const -> Vector;

    /**
     * @brief Returns the entropy of a specified phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Entropy in J/(K).
     *
     * @code
     * double ent0 = engine.phaseEntropy(0);
     * @endcode
     */
    auto phaseEntropy(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the entropy of a specified phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Entropy in J/(K).
     *
     * @code
     * double ent0 = engine.phaseEntropy("aq_gen");
     * @endcode
     */
    auto phaseEntropy(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the entropy of a specified phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Entropy in J/(K).
     *
     * @code
     * double ent0 = engine.phaseEntropy_i(0);
     * @endcode
     */
    auto phaseEntropy_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the heat capacities (Cp) of all phases.
     *
     * @return (VectorConstRef) Vector of heat capacities in J/K.
     *
     * @code
     * auto cps = engine.phaseHeatCapacitiesConstP();
     * @endcode
     */
    auto phaseHeatCapacitiesConstP(long int ndx) const -> Vector;

    /**
     * @brief Returns the heat capacity of a specified phase.
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Heat capacity in J/K.
     *
     * @code
     * double cp0 = engine.phaseHeatCapacityConstP(0);
     * @endcode
     */
    auto phaseHeatCapacityConstP(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the heat capacity of a specified phase.
     *
     * @param phase (std::string) Phase name.
     * @return (double) Heat capacity in J/K.
     *
     * @code
     * double cp0 = engine.phaseHeatCapacityConstP("aq_gen");
     * @endcode
     */
    auto phaseHeatCapacityConstP(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the heat capacity of a specified phase.
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Heat capacity in J/K.
     *
     * @code
     * double cp0 = engine.phaseHeatCapacityConstP_i(0);
     * @endcode
     */
    auto phaseHeatCapacityConstP_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the saturation indices of all phases (log₁₀ units).
     *
     * @return (VectorConstRef) Vector of saturation indices.
     *
     * @code
     * auto satIdx = engine.phaseSatIndices();
     * @endcode
     */
    auto phaseSatIndices(long int ndx) const -> Vector;

    /**
     * @brief Returns the saturation index of a specific phase (log₁₀ units).
     * Access specified element with bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Saturation index.
     *
     * @code
     * double sat0 = engine.phaseSatIndex(0);
     * @endcode
     */
    auto phaseSatIndex(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the saturation index of a specific phase (log₁₀ units).
     *
     * @param phase (std::string) Phase name.
     * @return (double) Saturation index.
     *
     * @code
     * double sat0 = engine.phaseSatIndex("aq_gen");
     * @endcode
     */
    auto phaseSatIndex(long int ndx, std::string phase) const -> double;

    /**
     * @brief Returns the saturation index of a specific phase (log₁₀ units).
     * Access specified element without bounds checking.
     *
     * @param iphase (Index) Phase index.
     * @return (double) Saturation index.
     *
     * @code
     * double sat0 = engine.phaseSatIndex_i(0);
     * @endcode
     */
    auto phaseSatIndex_i(long int ndx, Index iphase) const -> double;

    /**
     * @brief Returns the total mass of the system.
     *
     * @return (double) System mass in kg.
     *
     * @code
     * double sysMass = engine.systemMass();
     * std::cout << "System mass: " << sysMass << " kg" << std::endl;
     * @endcode
     */
    auto systemMass(long int ndx) const -> double;

    /**
     * @brief Returns the total volume of the system.
     *
     * @return (double) Volume in m³.
     *
     * @code
     * double sysVol = engine.systemVolume();
     * std::cout << "System volume: " << sysVol << " m³" << std::endl;
     * @endcode
     */
    auto systemVolume(long int ndx) const -> double;

    /**
     * @brief Returns the ionic strength of the aqueous phase.
     *
     * @return (double) Ionic strength in molal.
     *
     * @code
     * double ionicStr = engine.ionicStrength();
     * @endcode
     */
    auto ionicStrength(long int ndx) const -> double;

    /**
     * @brief Returns the pH of the aqueous phase.
     *
     * @return (double) pH (in the activity scale (-log10 molal)).
     *
     * @code
     * double ph = engine.pH();
     * @endcode
     */
    auto pH(long int ndx) const -> double;

    /**
     * @brief Returns the pe of the aqueous phase.
     *
     * @return (double) pe (in the activity scale (-log10 molal)).
     *
     * @code
     * double pe = engine.pe();
     * @endcode
     */
    auto pe(long int ndx) const -> double;

    /**
     * @brief Returns the Eh of the aqueous phase.
     *
     * @return (double) Eh (V).
     *
     * @code
     * double eh = engine.Eh();
     * @endcode
     */
    auto Eh(long int ndx) const -> double;

    /**
     * @brief Returns the total Gibbs energy of the system.
     *
     * @return (double) Gibbs energy in J/mol.
     *
     * @code
     * double sysGibbs = engine.systemGibbsEnergy();
     * @endcode
     */
    auto systemGibbsEnergy(long int ndx) const -> double;

    /**
     * @brief Returns the total enthalpy of the system.
     *
     * @return (double) Enthalpy in J.
     *
     * @code
     * double sysEnthalpy = engine.systemEnthalpy();
     * @endcode
     */
    auto systemEnthalpy(long int ndx) const -> double;

    /**
     * @brief Returns the total entropy of the system.
     *
     * @return (double) Entropy in J/K.
     *
     * @code
     * double sysEntropy = engine.systemEntropy();
     * @endcode
     */
    auto systemEntropy(long int ndx) const -> double;

    /**
     * @brief Returns the total isobaric heat capacity of the system.
     *
     * @return (double) Heat capacity in J/K.
     *
     * @code
     * double sysCp = engine.systemHeatCapacityConstP();
     * @endcode
     */
    auto systemHeatCapacityConstP(long int ndx) const -> double;
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
    auto gasPhaseName() const -> std::string;    /// TDB Return the molalities of the species. Aquatic systems only (assuming

    ChemicalEngineArray(const MassTransport& mt_class);
    std::shared_ptr<ChemicalEngineArray> getptr() {
        return shared_from_this();
    }

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;

    void setCorrectionbIC(long ndx);
    /**
     * @brief Reallocates internal arrays.
     *
     * Should be invoked after the system composition changes to ensure proper memory
     * configuration.
       */
    auto reallocateEngineArrays() -> void;

};

/**
   * @brief Outputs the state of the ChemicalEngine.
   *
   * Overloaded stream operator for convenient printing of engine details.
   *
   * @param out Output stream.
   * @param engine ChemicalEngine instance.
   * @return (std::ostream&) Updated output stream.
   *
   * @code
   * // Example: Print the engine state.
   * std::cout << engine << std::endl;
   * @endcode
   */
auto operator<<(std::ostream& out, const ChemicalEngineArray& engine) -> std::ostream&;

} // namespace xGEM2MT
