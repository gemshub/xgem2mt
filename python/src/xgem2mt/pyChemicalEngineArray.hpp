// xGEM2MT is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// xGEM2MT includes
#include <xGEM2MT/ChemicalEngineArray.hpp>
using namespace xGEM2MT;

/*
To do

1. auto initialize(std::string ipm_filename, std::string dbr_filename, std::vector<long int> node_types = {})
default parameter into python

2.    auto initializeFromJsonStrings(std::string dch_json, std::string ipm_json,
                                   std::vector<std::string> dbr_json, std::vector<long int> node_types = {})
default parameter into python

3. access function/property to
    /// Write IPM, DCH and DBR files (in binary, txt or json mode)
    GEMS3KGenerator::IOModes io_mode = GEMS3KGenerator::f_key_value;


*/


void exportChemicalEngineArray(py::module& m)
{

    m.def("update_loggers", &update_loggers,
          py::arg("use_cout"), py::arg("logfile_name"), py::arg("log_level"),
          R"doc(
  Updates logger settings for the ChemicalEngineArray.

  :param bool use_cout: If true, prints log messages to stdout.
  :param str logfile_name: The filename for a rotating log file. If empty, file logging is disabled.
  :param int log_level: Verbosity level for logging (trace=0, debug=1, info=2, warn=3, err=4, critical=5, off=6).

  **Example:**

  .. code-block:: python

      update_loggers(True, "chemical_log.txt", 2)
  )doc");



    py::class_<ChemicalEngineArrayOptions>(m, "ChemicalEngineArrayOptions",
        R"doc(
Configuration options for the ChemicalEngineArrayOptions.

Attributes:
    warmstart (bool): Use smart start for faster convergence.
    useSIA (char): Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-).
    step (int): Extern actual time iterator.
    cdv (float): Cutoff for IC amount differences in the node between time steps (mol, 1e-9).
    cez (float): Cutoff for minimal amounts of IC in node bulk compositions (mol, 1e-12).
)doc")
        .def(py::init<>())
        .def_readwrite("warmstart", &ChemicalEngineArrayOptions::warmstart,
            R"doc(
Whether to use smart start for faster convergence (may be less accurate).

:type: bool
)doc")
        .def_readwrite("useSIA", &ChemicalEngineArrayOptions::useSIA,
            R"doc(
Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-).

:type: char
)doc")
            .def_readwrite("step", &ChemicalEngineArrayOptions::step,
            R"doc(
Extern actual time iterator.

:type: int
)doc")
            .def_readwrite("cdv", &ChemicalEngineArrayOptions::cdv,
            R"doc(
Cutoff for IC amount differences in the node between time steps (mol, 1e-9).

:type: float
)doc")
            .def_readwrite("cez", &ChemicalEngineArrayOptions::cez,
            R"doc(
Cutoff for minimal amounts of IC in node bulk compositions (mol, 1e-12).
:type: float
)doc")
            .def("to_json_string", &ChemicalEngineArrayOptions::to_json_string, py::arg("dense")  = false)
            .def("__repr__", [](const ChemicalEngineArrayOptions& self) {
              std::stringstream ss; ss << self.to_json_string(); return ss.str();
         }, R"doc(
Returns a string representation of the ChemicalEngineArrayOptions instance.

**Example:**

.. code-block:: python

print(engine)
)doc");



    py::class_<ChemicalEngineArray, std::shared_ptr<ChemicalEngineArray>>(m, "ChemicalEngineArray",
                                                                          R"doc(
Class for equilibrium computations and thermodynamic analysis for the list of the systems.

 The ChemicalEngineArray performs Gibbs energy minimization to compute the equilibrium state of a chemical systems.
 Its API provides methods to load system data, update  component amounts, set options, and query resulting
 thermodynamic properties for all systems.
)doc")
            .def(py::init<size_t>(),
             R"doc(
  Default constructor for the ChemicalEngineArray.

  Creates an empty ChemicalEngineArray instance for the 1D arrangement of nodes.
  Use one of the initialization routines
  (such as ``initialize()`` or ``initializeFromJsonStrings()``) to configure the systems.

  :param int size:  Number of nodes (systems) in the list.

  **Example:**

  .. code-block:: python

        engine = ChemicalEngineArray(5)
        node_types=[1,0,0,0,0]
        gems_init_lst = "nodearray-data/CalcColumn-dat.lst"
        gems_init_dbr_lst = "nodearray-data/CalcColumn-dbr.lst"

        engine.initialize(gems_init_lst, gems_init_dbr_lst, node_types)
  )doc")
        .def_property("options",
                      &ChemicalEngineArray::options,
                      &ChemicalEngineArray::setOptions,
             R"doc(
  Gets or sets the options used by the ChemicalEngineArray.

:type: ChemicalEngineArrayOptions

**Example (get):**

.. code-block:: python

  opts = engine.options
  print("Warmstart:", opts.warmstart)

**Example (set):**

.. code-block:: python

  opts = engine.options
  opts.warmstart = False
  engine.options = opts
)doc")
            .def("setWarmStart", &ChemicalEngineArray::setWarmStart,
             R"doc(
Enables warm start for the ChemicalEngineArray.

Warm start uses the previous solution as an initial guess for faster convergence.

**Example:**

.. code-block:: python

  engine.setWarmStart()
)doc")
            .def("setColdStart", &ChemicalEngineArray::setColdStart,
             R"doc(
Enables cold start for the ChemicalEngineArray.

Cold start resets the initial guess to default values for robust convergence.

**Example:**

.. code-block:: python

  engine.setColdStart()
)doc")
            .def("nNodes", &ChemicalEngineArray::nNodes,
             R"doc(
Get the total number of nodes in the node array.

:return int: The size of the list.

**Example:**

.. code-block:: python

  size = engine.nNodes()
)doc")
            .def("initialize", &ChemicalEngineArray::initialize,
             R"doc(
 Initializes the ChemicalEngineArray from a GEM-Selektor project files.

 Reads the files and configures the chemical systems defined inside it.
 :param str  ipm_filename: Path to the system definition file (e.g., "my-system-dat.lst").
 :param str  dbr_filename: Path to the text file containing the list of DBR input files.
 :param list[int]  node_types: The initial node contents from DATABR files
                        will be distributed among nodes in array.

  **Example:**

  .. code-block:: python

        engine = ChemicalEngineArray(5)
        node_types=[1,0,0,0,0]
        gems_init_lst = "nodearray-data/CalcColumn-dat.lst"
        gems_init_dbr_lst = "nodearray-data/CalcColumn-dbr.lst"

        engine.initialize(gems_init_lst, gems_init_dbr_lst, node_types)
  )doc")
            .def("initializeFromJsonStrings", &ChemicalEngineArray::initializeFromJsonStrings,
             R"doc(
Initializes the ChemicalEngineArray using JSON strings.

Configures the engine from three JSON strings:
 - dch_json: Defines the chemical system.
 - ipm_json: Contains calculation parameters and settings.
 - dbr_json: Specifies the input node composition.
 These strings can be obtained from the GEM-Selektor exported GEMS3K system files as JSON (*-imp.json, *-dbc.json, *-dbr.json).

:param str dch_json: JSON string describing the chemical system.
:param str ipm_json: JSON string with parameter and algorithm settings.
:param list[str] dbr_json: JSON strings containing input nodes composition details.
:param list[int] node_types: The initial node contents from DATABR files will be distributed among
                              nodes in array according to the distribution index list.

  **Example:**

  .. code-block:: python

        engine = ChemicalEngineArray(5)
        node_types=[1,0,0,0,0]

        dch_json = Path("nodearray-data/CalcColumn-dch.json").read_text().replace('\n', ' ')
        ipm_json = Path("nodearray-data/CalcColumn-ipm.json").read_text().replace('\n', ' ')
        dbr_json1 = Path("nodearray-data/CalcColumn-dbr-0-0000.json").read_text().replace('\n', ' ')
        dbr_json0 = Path("nodearray-data/CalcColumn-dbr-0-0001.json").read_text().replace('\n', ' ')
        dbr_jsons = [dbr_json1, dbr_json0]

        engine.initializeFromJsonStrings(dch_json, ipm_json, dbr_jsons, node_types)
  )doc")
            .def("reequilibrate", py::overload_cast<size_t, bool>(&ChemicalEngineArray::reequilibrate),
             R"doc(
Re-equilibrates the system with (or without) a warm start for the node.

:param int ndx: size_t of the node (system) in the list.
:param bool warmstart: If true, uses previous data as an initial guess.
:return int: Return code of the re-equilibration.

**Example:**

.. code-block:: python

engine.reequilibrate(1, False)

**Return Codes**

   The function returns an integer code indicating the status:

   - 0: No GEM re-calculation needed
   - 1: Need GEM calculation with LPP (automatic) initial approximation (AIA)
   - 2: OK after GEM calculation with LPP AIA
   - 3: Bad (not fully trustful) result after GEM calculation with LPP AIA
   - 4: Failure (no result) in GEM calculation with LPP AIA
   - 5: Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation
   - 6: OK after GEM calculation with SIA
   - 7: Bad (not fully trustful) result after GEM calculation with SIA
   - 8: Failure (no result) in GEM calculation with SIA
   - 9: Terminal error in GEMS3K (e.g., memory corruption). Restart required.

  )doc")
            .def("reequilibrate", py::overload_cast<bool>(&ChemicalEngineArray::reequilibrate),
             R"doc(
Re-equilibrates the all system with (or without) a warm start for the node.

:param bool warmstart: If true, uses previous data as an initial guess.

**Example:**

.. code-block:: python

engine.reequilibrate(False)
    )doc")
            .def("equilibrate", py::overload_cast<size_t>(&ChemicalEngineArray::equilibrate),
             R"doc(
Equilibrates the system in node ndx, according options and node settings.

:param int ndx: size_t of the node (system) in the list.
:return bool: Return True if OK after GEM calculation with LPP AIA or SIA.

**Example:**

.. code-block:: python

   engine.equilibrate(1)
)doc")
            .def("equilibrate", py::overload_cast<>(&ChemicalEngineArray::equilibrate),
             R"doc(
Equilibrates all systems, according options and node settings.

:return bool: Return True if OK after GEM calculation with LPP AIA or SIA for all systems.

**Example:**

.. code-block:: python

   engine.equilibrate()
)doc")
            .def("dumpGEMS3K", &ChemicalEngineArray::dumpGEMS3K, py::arg("filepath"), py::arg("put_all")=false,
             R"doc(
Dump the MULTI, DATACH and  all DATABR files.

:param str filepath: IPM work structure file path&name
:param bool put_all: If True print array of nodes for previous time point
                       and for current time point

**Example:**

.. code-block:: python

engine.dumpGEMS3K("CalcColumn-dump-dat.lst", False)
)doc")
            .def("nextTimePoint", &ChemicalEngineArray::nextTimePoint,
             R"doc(
Copies array of nodes for current time point to the same place in the array of nodes for previous time point.

:param bool always: If False, copy only the nodes with essential changes  into the element quantities.

**Example:**

.. code-block:: python

engine.nextTimePoint(True)
)doc")
            .def("readDbrFromFile", &ChemicalEngineArray::readDbrFromFile,
             R"doc(
Reads a DBR file from disk to the node.

Loads the structure of node-dependent data (composition, T, P and other).

 :param int ndx: size_t of the node (system) in the list.
 :param str filename: Path to the input DBR file (e.g., "*-dbr.json/.dat").

**Example:**

.. code-block:: python

engine.readDbrFromFile(0, "my-system-dbr-0-0000.json")
)doc")
            .def("readDbrFromJsonString", &ChemicalEngineArray::readDbrFromJsonString,
             R"doc(
Reads System DBR composition from a JSON string to the node.

Loads the structure of node-dependent data (composition, T, P and other).

 :param int ndx: size_t of the node (system) in the list.
 :param str dbr_json: JSON string containing DBR composition data.

**Example:**

.. code-block:: python

dbr_json1 = Path("nodearray-data/CalcColumn-dbr-0-0000.json").read_text().replace('\n', ' ')
engine.readDbrFromJsonString(0, dbr_json1)
)doc")
            .def("writeDbrToFile", &ChemicalEngineArray::writeDbrToFile,
             R"doc(
Writes the current DBR for the node to a file according current io mode.

Saves the state of the system (after equilibrium calculations) into a DBR file.

:param int ndx: size_t of the node (system) in the list.
:param str filename: Path to the output DBR file (e.g., "my-system-dbr-0-0001.dat").

**Example:**

.. code-block:: python

engine.writeDbrToFile(1, "output.dbr")
)doc")
            .def("writeDbrToJsonString", &ChemicalEngineArray::writeDbrToJsonString,
             R"doc(
Returns the current DBR as a JSON string.

:param int ndx: size_t of the node (system) in the list.
:return: JSON string representing the system composition.
:rtype: str

**Example:**

.. code-block:: python

json_state = engine.writeDbrToJsonString(1)
print(json_state)
)doc")
        .def("numElements", &ChemicalEngineArray::numElements,
             R"doc(
  Returns the number of elements in the chemical system.

  **Example:**

  .. code-block:: python

      num_elements = engine.numElements()
      print(num_elements)
  )doc")
        .def("numSpecies", &ChemicalEngineArray::numSpecies,
             R"doc(
  Returns the number of species in the chemical system.

  **Example:**

  .. code-block:: python

      num_species = engine.numSpecies()
      print(num_species)
  )doc")
        .def("numPhases", &ChemicalEngineArray::numPhases,
             R"doc(
  Returns the number of phases in the chemical system.

  **Example:**

  .. code-block:: python

      num_phases = engine.numPhases()
      print(num_phases)
  )doc")
        .def("numSpeciesInPhase", static_cast<size_t (ChemicalEngineArray::*)(size_t) const>(&ChemicalEngineArray::numSpeciesInPhase),
             py::arg("iphase"),
             R"doc(
  Returns the number of species in a given phase.
  Access specified element with bounds checking.

  :param int iphase: Index of the phase.

  **Example:**

  .. code-block:: python

      num_species_in_phase = engine.numSpeciesInPhase(0)
      print(num_species_in_phase)
  )doc")
        .def("numSpeciesInPhase", static_cast<size_t (ChemicalEngineArray::*)(std::string) const>(&ChemicalEngineArray::numSpeciesInPhase),
             py::arg("phase"),
             R"doc(
  Returns the number of species in in a phase identified by name.

  :param str phase: Name of the phase.

  **Example:**

  .. code-block:: python

      num_species_in_phase = engine.numSpeciesInPhase("aq_gen")
      print(num_species_in_phase)
  )doc")
        .def("numSpeciesInPhase_i", &ChemicalEngineArray::numSpeciesInPhase_i,
             py::arg("iphase"),
             R"doc(
  Returns the number of species in a given phase.
  Access specified element without bounds checking.

  :param int iphase: Index of the phase.

  **Example:**

  .. code-block:: python

      num_species_in_phase = engine.numSpeciesInPhase_i(0)
      print(num_species_in_phase)
  )doc")
        .def("elementName", &ChemicalEngineArray::elementName, py::arg("index"),
             R"doc(
  Returns the name of an element by its index.
  Access specified element with bounds checking.

  :param int index: Index of the element.

  **Example:**

  .. code-block:: python

      element_name = engine.elementName(0)
      print(element_name)
  )doc")
        .def("elementName_i", &ChemicalEngineArray::elementName_i, py::arg("index"),
             R"doc(
  Returns the name of an element by its index.
  Access specified element without bounds checking.

  :param int index: Index of the element.

  **Example:**

  .. code-block:: python

      element_name = engine.elementName_i(0)
      print(element_name)
  )doc")
        .def("speciesName", &ChemicalEngineArray::speciesName, py::arg("index"),
             R"doc(
  Returns the name of a species by its index.
  Access specified element with bounds checking.

  :param int index: Index of the species.

  **Example:**

  .. code-block:: python

      species_name = engine.speciesName(0)
      print(species_name)
  )doc")
        .def("speciesName_i", &ChemicalEngineArray::speciesName_i, py::arg("index"),
             R"doc(
  Returns the name of a species by its index.
  Access specified element without bounds checking.

  :param int index: Index of the species.

  **Example:**

  .. code-block:: python

      species_name = engine.speciesName_i(0)
      print(species_name)
  )doc")
        .def("speciesCharge", static_cast<double (ChemicalEngineArray::*)(size_t) const>(&ChemicalEngineArray::speciesCharge), py::arg("index"),
             R"doc(
  Returns the electrical charge of a species by its index.
  Access specified element with bounds checking.

  :param int index: Index of the species.

  **Example:**

  .. code-block:: python

      charge = engine.speciesCharge(0)
      print(charge)
  )doc")
        .def("speciesCharge", static_cast<double (ChemicalEngineArray::*)(std::string, std::optional<std::string>) const>(&ChemicalEngineArray::speciesCharge),
             py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
  Returns the electrical charge of a species identified by name.

  :param str species: Name of the species.
  :param str phase: Name of the phase the species was included in. If None get the first index.

  **Example:**

  .. code-block:: python

      charge = engine.speciesCharge("SiO2")
      print(charge)
  )doc")
        .def("speciesCharge_i", &ChemicalEngineArray::speciesCharge_i, py::arg("index"),
             R"doc(
  Returns the electrical charge of a species by its index.
  Access specified element without bounds checking.

  :param int index: Index of the species.

  **Example:**

  .. code-block:: python

      charge = engine.speciesCharge_i(0)
      print(charge)
  )doc")
        .def("phaseName", &ChemicalEngineArray::phaseName, py::arg("index"),
             R"doc(
  Returns the name of a phase by its index.
  Access specified element with bounds checking.

  :param int index: Index of the phase.

  **Example:**

  .. code-block:: python

      phase_name = engine.phaseName(0)
      print(phase_name)
  )doc")
        .def("phaseName_i", &ChemicalEngineArray::phaseName_i, py::arg("index"),
             R"doc(
  Returns the name of a phase by its index.
  Access specified element without bounds checking.

  :param int index: Index of the phase.

  **Example:**

  .. code-block:: python

      phase_name = engine.phaseName_i(0)
      print(phase_name)
  )doc")
        .def("indexElement", &ChemicalEngineArray::indexElement, py::arg("name"),
             R"doc(
  Returns the index of an element by its name.

  :param str name: Name of the element.

  **Example:**

  .. code-block:: python

      index = engine.indexElement("O")
      print(index)
  )doc")
        .def("indexSpecies", &ChemicalEngineArray::indexSpecies, py::arg("name"), py::arg("phase") = py::none(),
             R"doc(
Returns the index of a species by its name.

:param str name: Name of the species.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  index = engine.indexSpecies("H2O@")
  print(index)
)doc")

        .def("indexSpeciesAll", &ChemicalEngineArray::indexSpeciesAll, py::arg("name"),
             R"doc(
Returns all indices of species matching the specified name.

:param str name: Name of the species.

**Example:**

.. code-block:: python

  indices = engine.indexSpeciesAll("H2O@")
  print(indices)
)doc")
        .def("indexSpeciesMap", &ChemicalEngineArray::indexSpeciesMap, py::arg("name"),
             R"doc(
Returns dictionary of all phase name and corresponding species index.

:param str name: Name of the species.

**Example:**

.. code-block:: python

  dict = engine.indexSpeciesMap("Anorthite")
  print(dict)
)doc")
        .def("indexPhase", &ChemicalEngineArray::indexPhase, py::arg("name"),
             R"doc(
Returns the index of a phase by its name.

:param str name: Name of the phase.

**Example:**

.. code-block:: python

  index = engine.indexPhase("aq_gen")
  print(index)
)doc")

        .def("indexPhaseAll", &ChemicalEngineArray::indexPhaseAll, py::arg("name"),
             R"doc(
Returns all indices of phases matching the specified name.

:param str name: Name of the phase.

**Example:**

.. code-block:: python

  indices = engine.indexPhaseAll("aq_gen")
  print(indices)
)doc")

        .def("indexPhaseWithSpecies", &ChemicalEngineArray::indexPhaseWithSpecies, py::arg("species_index"),
             R"doc(
Returns the index of the phase containing a given species.

:param int species_index: Index of the species.

**Example:**

.. code-block:: python

  phase_index = engine.indexPhaseWithSpecies(0)
  print(phase_index)
)doc")

        .def("indexFirstSpeciesInPhase", static_cast<size_t (ChemicalEngineArray::*)(size_t) const>(&ChemicalEngineArray::indexFirstSpeciesInPhase),
             py::arg("iphase"),
             R"doc(
Returns the index of the first species in a specified phase.
Access specified element with bounds checking.

:param int iphase: Index of the phase.

**Example:**

.. code-block:: python

  species_index = engine.indexFirstSpeciesInPhase(0)
  print(species_index)
)doc")
        .def("indexFirstSpeciesInPhase", static_cast<size_t (ChemicalEngineArray::*)(std::string) const>(&ChemicalEngineArray::indexFirstSpeciesInPhase),
             py::arg("phase"),
             R"doc(
Returns the index of the first species in a specified phase.

:param str phase: Name of the phase.

**Example:**

.. code-block:: python

  species_index = engine.indexFirstSpeciesInPhase("aq_gen")
  print(species_index)
)doc")
        .def("indexFirstSpeciesInPhase_i", &ChemicalEngineArray::indexFirstSpeciesInPhase_i, py::arg("iphase"),
             R"doc(
Returns the index of the first species in a specified phase.
Access specified element without bounds checking.

:param int iphase: Index of the phase.

**Example:**

.. code-block:: python

  species_index = engine.indexFirstSpeciesInPhase_i(0)
  print(species_index)
)doc")

        .def("elementMolarMasses", &ChemicalEngineArray::elementMolarMasses,
             py::return_value_policy::reference_internal,
             R"doc(
Returns the molar masses of all elements in the system (kg/mol).

**Example:**

.. code-block:: python

  molar_masses = engine.elementMolarMasses()
  print(molar_masses)
)doc")
        .def("speciesMolarMasses", &ChemicalEngineArray::speciesMolarMasses,
             py::return_value_policy::reference_internal,
             R"doc(
Returns the molar masses of all species in the system (kg/mol).

**Example:**

.. code-block:: python

  species_molar_masses = engine.speciesMolarMasses()
  print(species_molar_masses)
)doc")

        .def("formulaMatrix", &ChemicalEngineArray::formulaMatrix,
             py::return_value_policy::reference_internal,
             R"doc(
Returns the formula matrix of the system.

The formula matrix represents the stoichiometric coefficients of elements in species.

**Example:**

.. code-block:: python

  formula_matrix = engine.formulaMatrix()
  print(formula_matrix)
)doc")
        .def("setSpeciesUpperLimit", static_cast<void (ChemicalEngineArray::*)(size_t, std::string, double, std::optional<std::string>)>(&ChemicalEngineArray::setSpeciesUpperLimit),
             py::arg("ndx"), py::arg("name"), py::arg("limit"), py::arg("phase") = py::none(),
             R"doc(
Sets the upper limit for a species by its name for the node.

:param int ndx: Index of the node (system) in the list.
:param str name: Name of the species.
:param float limit: Upper limit for the species amount.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  engine.setSpeciesUpperLimit(1, "H2O@", 10.0)
)doc")
        .def("setSpeciesUpperLimit", static_cast<void (ChemicalEngineArray::*)(size_t, size_t, double)>(&ChemicalEngineArray::setSpeciesUpperLimit),
             py::arg("ndx"), py::arg("index"), py::arg("limit"),
             R"doc(
Sets the upper limit for a species by its index for the node.
Access specified element with bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int index: Index of the species.
:param float limit: Upper limit for the species amount.

**Example:**

.. code-block:: python

  engine.setSpeciesUpperLimit(1, 0, 10.0)
)doc")
        .def("setSpeciesUpperLimit_i", &ChemicalEngineArray::setSpeciesUpperLimit_i,
             py::arg("ndx"), py::arg("index"), py::arg("limit"),
             R"doc(
Sets the upper limit for a species by its index for the node.
Access specified element without bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int index: Index of the species.
:param float limit: Upper limit for the species amount.

**Example:**

.. code-block:: python

  engine.setSpeciesUpperLimit_i(1, 0, 10.0)
)doc")
        .def("setSpeciesLowerLimit", static_cast<void (ChemicalEngineArray::*)(size_t, std::string, double, std::optional<std::string>)>(&ChemicalEngineArray::setSpeciesLowerLimit),
             py::arg("ndx"), py::arg("name"), py::arg("limit"), py::arg("phase") = py::none(),
             R"doc(
Sets the lower limit for a species by its name for the node.

:param int ndx: Index of the node (system) in the list.
:param str name: Name of the species.
:param float limit: Lower limit for the species amount.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  engine.setSpeciesLowerLimit(1, "SiO2", 0.1)
)doc")
        .def("setSpeciesLowerLimit", static_cast<void (ChemicalEngineArray::*)(size_t, size_t, double)>(&ChemicalEngineArray::setSpeciesLowerLimit),
             py::arg("ndx"), py::arg("index"), py::arg("limit"),
             R"doc(
Sets the lower limit for a species by its index for the node.
Access specified element with bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int index: Index of the species.
:param float limit: Lower limit for the species amount.

**Example:**

.. code-block:: python

  engine.setSpeciesLowerLimit(1, 0, 0.1)
)doc")
        .def("setSpeciesLowerLimit_i", &ChemicalEngineArray::setSpeciesLowerLimit_i,
             py::arg("ndx"), py::arg("index"), py::arg("limit"),
             R"doc(
Sets the lower limit for a species by its index for the node.
Access specified element without bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int index: Index of the species.
:param float limit: Lower limit for the species amount.

**Example:**

.. code-block:: python

  engine.setSpeciesLowerLimit_i(1, 0, 0.1)
)doc")
        .def("setSpeciesAmount", static_cast<void (ChemicalEngineArray::*)(size_t, std::string, double, std::optional<std::string>)>(&ChemicalEngineArray::setSpeciesAmount),
             py::arg("ndx"), py::arg("name"), py::arg("amount"), py::arg("phase") = py::none(),
             R"doc(
Sets the amount of a species by its name for the node.

:param int ndx: Index of the node (system) in the list.
:param str name: Name of the species.
:param float amount: Amount of the species in moles.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  engine.setSpeciesAmount(1, "SiO2", 1.0)
)doc")
        .def("setSpeciesAmount", static_cast<void (ChemicalEngineArray::*)(size_t, size_t, double)>(&ChemicalEngineArray::setSpeciesAmount),
             py::arg("ndx"), py::arg("ispecies"), py::arg("amount"),
             R"doc(
Sets the amount of a species by its index for the node.
Access specified element with bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int ispecies: Index of the species.
:param float amount: Amount of the species in moles.

**Example:**

.. code-block:: python

  engine.setSpeciesAmount(1, 0, 1.0)
)doc")
        .def("setSpeciesAmount_i", &ChemicalEngineArray::setSpeciesAmount_i,
             py::arg("ndx"), py::arg("ispecies"), py::arg("amount"),
             R"doc(
Sets the amount of a species by its index for the node.
Access specified element without bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int ispecies: Index of the species.
:param float amount: Amount of the species in moles.

**Example:**

.. code-block:: python

  engine.setSpeciesAmount_i(1, 0, 1.0)
)doc")

        .def("setStandardMolarGibbsEnergy", static_cast<void (ChemicalEngineArray::*)(std::string, double, double, double, std::optional<std::string>)>(&ChemicalEngineArray::setStandardMolarGibbsEnergy),
             py::arg("name"), py::arg("TK"), py::arg("P"), py::arg("value"), py::arg("phase") = py::none(),
             R"doc(
Sets the standard molar Gibbs energy for a species (J/mol).

:param str name: Name of the species.
:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param float value: Standard molar Gibbs energy value (J/mol).
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

  engine.setStandardMolarGibbsEnergy("H2O", 298.15, 101325, -237.13)
)doc")
        .def("setStandardMolarGibbsEnergy", static_cast<void (ChemicalEngineArray::*)(size_t, double, double, double)>(&ChemicalEngineArray::setStandardMolarGibbsEnergy),
             py::arg("ispecies"), py::arg("TK"), py::arg("P"), py::arg("value"),
             R"doc(
Sets the standard molar Gibbs energy for a species (J/mol).
Access specified element with bounds checking.

:param int ispecies: Index of the species.
:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param float value: Standard molar Gibbs energy value (J/mol).

**Example:**

.. code-block:: python

  engine.setStandardMolarGibbsEnergy(0, 298.15, 101325, -237.13)
)doc")
        .def("setStandardMolarGibbsEnergy_i", &ChemicalEngineArray::setStandardMolarGibbsEnergy_i,
             py::arg("ispecies"), py::arg("TK"), py::arg("P"), py::arg("value"),
             R"doc(
Sets the standard molar Gibbs energy for a species (J/mol).
Access specified element with bounds checking.

:param int ispecies: Index of the species.
:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param float value: Standard molar Gibbs energy value (J/mol).

**Example:**

.. code-block:: python

  engine.setStandardMolarGibbsEnergy_i(0, 298.15, 101325, -237.13)
    )doc")
        .def("setPT", &ChemicalEngineArray::setPT,
             py::arg("ndx"), py::arg("P"), py::arg("T"),
             R"doc(
 Sets the temperature and pressure of the system for the node.

:param int ndx: Index of the node (system) in the list.
:param float P: Pressure in Pascals (Pa).
:param float T: Temperature in Kelvin (K).

:return bool: True if PT was set correctly, False if out of range.

**Example:**

.. code-block:: python

engine.setPT(101325, 298.15)
    )doc")

        .def("setB", &ChemicalEngineArray::setB,
             py::arg("ndx"), py::arg("b"),
             R"doc(
Sets the bulk composition of the system for the node.

:param int ndx: Index of the node (system) in the list.
:param list[float] composition: List of element amounts (mol).

Set element amounts of a solution with 0.001 M MgCl2 and 0.001 M CO2:
  - **elements order:** C, Ca, Cl, H, Mg, O, Zz (charge)
  - **Note:** The order of elements in the list should match the order in the chemical system.

**Example:**

.. code-block:: python

engine.setB(1, [0.001, 1e-09, 0.004, 110.6837, 0.002, 55.34405, 0])
    )doc")

    .def("converged", &ChemicalEngineArray::converged,
         R"doc(
Checks if the equilibrium computation has converged for the node.

:param int ndx: Index of the node (system) in the list.
:return bool: True if the computation has converged, False otherwise.

**Example:**

.. code-block:: python

is_converged = engine.converged(1)
print(is_converged)
)doc")

        .def("numIterations", &ChemicalEngineArray::numIterations,
             R"doc(
Returns the number of iterations performed during the last equilibrium computation for the node.

:param int ndx: Index of the node (system) in the list.
:return int: Number of iterations.

**Example:**

.. code-block:: python

iterations = engine.numIterations(1)
print(iterations)
   )doc")

        .def("temperature", &ChemicalEngineArray::temperature,
        R"doc(
Returns the current temperature of the system (K) for the node.

:param int ndx: Index of the node (system) in the list.
:return float: Temperature in Kelvin.

  **Example:**

  .. code-block:: python

 temp = engine.temperature(1)
 print(temp)
            )doc")

        .def("pressure", &ChemicalEngineArray::pressure,
            R"doc(
Returns the current pressure of the system (Pa) for the node.

:param int ndx: Index of the node (system) in the list.
:return float: Pressure in Pascals.

**Example:**

.. code-block:: python

pressure = engine.pressure(1)
print(pressure)
       )doc")

        .def("elementAmounts", &ChemicalEngineArray::elementAmounts,
         R"doc(
Returns the amounts of all elements in the system (mol) for the node.

:param int ndx: Index of the node (system) in the list.
:return list[float]: All elements amounts.

**Example:**

.. code-block:: python

 amounts = engine.elementAmounts(1)
 print(amounts)
   )doc")

        .def("elementAmountsInPhase", static_cast<Vector (ChemicalEngineArray::*)(size_t, size_t) const>(&ChemicalEngineArray::elementAmountsInPhase),
        py::arg("ndx"), py::arg("iphase"),
        R"doc(
Returns the amounts of all elements in a specific phase (mol) for the node.
Access specified element with bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int iphase: Index of the phase.
:return list[float]: All elements amounts in phase.

  **Example:**

  .. code-block:: python

 amounts = engine.elementAmountsInPhase(1, 0)
 print(amounts)
            )doc")
        .def("elementAmountsInPhase", static_cast<Vector (ChemicalEngineArray::*)(size_t, std::string) const>(&ChemicalEngineArray::elementAmountsInPhase),
        py::arg("ndx"), py::arg("phase"),
             R"doc(
Returns the amounts of all elements in a specific phase (mol) for the node.

:param int ndx: Index of the node (system) in the list.
:param str phase: Name of the phase.
:return list[float]: All elements amounts.

**Example:**

.. code-block:: python

amounts = engine.elementAmountsInPhase(1, "aq_gen"")
 print(amounts)
       )doc")
        .def("elementAmountsInPhase_i", &ChemicalEngineArray::elementAmountsInPhase_i,
            py::arg("ndx"), py::arg("iphase"),
            R"doc(
Returns the amounts of all elements in a specific phase (mol) for the node.
Access specified element without bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int iphase: Index of the phase.
:return list[float]: All elements amounts.

**Example:**

.. code-block:: python

amounts = engine.elementAmountsInPhase_i(1, 0)
print(amounts)
       )doc")

        .def("elementAmountsInSpecies", &ChemicalEngineArray::elementAmountsInSpecies,
            py::arg("ndx"), py::arg("ispecies"),
        R"doc(
Returns the amounts of all elements in a specific species (mol).

:param int ndx: Index of the node (system) in the list.
:param list[int] ispecies: Vector indexes of the species.
:return list[float]: All elements amounts.

    )doc")

        .def("speciesAmount", static_cast<double (ChemicalEngineArray::*)(size_t, size_t) const>(&ChemicalEngineArray::speciesAmount),
            py::arg("ndx"), py::arg("ispecies"),
             R"doc(
Returns the amount of a species by its index (mol) for the node.
Access specified element with bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int ispecies: Index of the species.
:return list[float]: All species amounts.

**Example:**

.. code-block:: python

amount = engine.speciesAmount(1, 0)
print(amount)
        )doc")

        .def("speciesAmount", static_cast<double (ChemicalEngineArray::*)(size_t, std::string, std::optional<std::string>) const>(&ChemicalEngineArray::speciesAmount),
        py::arg("ndx"), py::arg("name"), py::arg("phase") = py::none(),
             R"doc(
Returns the amount of a species by its name (mol) for the node.

:param int ndx: Index of the node (system) in the list.
:param str name: Name of the species.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

amount = engine.speciesAmount(1, "H2O@")
print(amount)
        )doc")

        .def("speciesAmount_i", &ChemicalEngineArray::speciesAmount_i,
            py::arg("ndx"), py::arg("ispecies"),
            R"doc(
Returns the amount of a species by its index (mol) for the node.
Access specified element without bounds checking.

:param int ndx: Index of the node (system) in the list.
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

amount = engine.speciesAmount_i(1, 0)
print(amount)
        )doc")

        .def("setSpeciesAmounts", &ChemicalEngineArray::setSpeciesAmounts,
         py::arg("ndx"), py::arg("amounts"),
          R"doc(
Sets the amounts of all species in the system (mol) for the node.

:param int ndx: Index of the node (system) in the list.
:param list[float] amounts: List of all species amounts.

**Example:**

.. code-block:: python

engine.setSpeciesAmounts(1, [1.0, 0.5, 0.2])
            )doc")

        .def("speciesAmounts", &ChemicalEngineArray::speciesAmounts,
             R"doc(
Returns the amounts of all species in the system (mol) for the node.

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

amounts = engine.speciesAmounts(1)
print(amounts)
  )doc")

        .def("setSpeciesUpperLimits", &ChemicalEngineArray::setSpeciesUpperLimits,
             R"doc(
Sets the upper limits for all species in the system (mol) for the node.

:param int ndx: Index of the node (system) in the list.
:param list[float] limits: List of upper limits for species amounts.

**Example:**

.. code-block:: python

engine.setSpeciesUpperLimits(1, [10.0, 5.0, 2.0])
    )doc")

        .def("setSpeciesLowerLimits", &ChemicalEngineArray::setSpeciesLowerLimits,
             R"doc(
Sets the lower limits for all species in the system (mol) for the node.

:param int ndx: Index of the node (system) in the list.
:param list[float] limits: List of lower limits for species amounts.

**Example:**

.. code-block:: python

engine.setSpeciesLowerLimits(1, [0.1, 0.05, 0.02])
            )doc")

        .def("speciesUpperLimits", &ChemicalEngineArray::speciesUpperLimits,
             R"doc(
Returns the upper limits for all species in the system (mol) for the node.

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

limits = engine.speciesUpperLimits(1)
print(limits)
            )doc")

        .def("speciesLowerLimits", &ChemicalEngineArray::speciesLowerLimits,
             R"doc(
Returns the lower limits for all species in the system (mol) for the node.

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

limits = engine.speciesLowerLimits(1)
print(limits)
    )doc")

        .def("moleFractions", &ChemicalEngineArray::moleFractions,
             R"doc(
Returns the mole fractions of all species in the system for the node.

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

fractions = engine.moleFractions(1)
print(fractions)
    )doc",
             py::return_value_policy::reference_internal)

        .def("speciesMolalities", &ChemicalEngineArray::speciesMolalities,
             R"doc(
Returns the molalities of all species in the system (mol/kg) for the node.

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

molalities = engine.speciesMolalities(1)
print(molalities)
    )doc")

        .def("lnActivityCoefficients", &ChemicalEngineArray::lnActivityCoefficients,
             R"doc(
Returns the natural logarithms of activity coefficients for all species for the node.

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

ln_coeffs = engine.lnActivityCoefficients(1)
print(ln_coeffs)
    )doc")

        .def("lnActivities", &ChemicalEngineArray::lnActivities,
             R"doc(
Returns the natural logarithms of activities for all species for the node.

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

ln_activities = engine.lnActivities(1)
print(ln_activities)
    )doc")

        .def("lnConcentrations", &ChemicalEngineArray::lnConcentrations,
             R"doc(
Returns the natural logarithms of concentrations for all species for the node.

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

ln_concentrations = engine.lnConcentrations(1)
print(ln_concentrations)
    )doc")

        .def("chemicalPotentials", &ChemicalEngineArray::chemicalPotentials,
             R"doc(
Returns the chemical potentials of all species in the node (J/mol).

:param int ndx: Index of the node (system) in the list.

**Example:**

.. code-block:: python

potentials = engine.chemicalPotentials(1)
print(potentials)
    )doc")

        .def("standardMolarGibbsEnergy", static_cast<double (ChemicalEngineArray::*)(double, double, size_t) const>(&ChemicalEngineArray::standardMolarGibbsEnergy),
             py::arg("TK"), py::arg("P"), py::arg("ispecies"),
             R"doc(
Returns the standard molar Gibbs energy of a specific species by its index (J/mol).
Access specified element with bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

gibbs_energy = engine.standardMolarGibbsEnergy(298.15, 101325, 0)
print(gibbs_energy)
        )doc")
        .def("standardMolarGibbsEnergy", static_cast<double (ChemicalEngineArray::*)(double, double, std::string, std::optional<std::string>) const>(&ChemicalEngineArray::standardMolarGibbsEnergy),
             py::arg("TK"), py::arg("P"), py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
Returns the standard molar Gibbs energy of a specific species by its index (J/mol).

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param str species: Name of the species.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

gibbs_energy = engine.standardMolarGibbsEnergy(298.15, 101325, "SiO2")
print(gibbs_energy)
                )doc")
        .def("standardMolarGibbsEnergy_i", &ChemicalEngineArray::standardMolarGibbsEnergy_i,
             py::arg("TK"), py::arg("P"), py::arg("ispecies"),
             R"doc(
Returns the standard molar Gibbs energy of a specific species by its index (J/mol).
Access specified element without bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

gibbs_energy = engine.standardMolarGibbsEnergy_i(298.15, 101325, 0)
print(gibbs_energy)
        )doc")

        .def("standardMolarEnthalpy", static_cast<double (ChemicalEngineArray::*)(double, double, size_t) const>(&ChemicalEngineArray::standardMolarEnthalpy),
             py::arg("TK"), py::arg("P"), py::arg("ispecies"),
             R"doc(
Returns the standard molar enthalpy of a specific species by its index (J/mol).
Access specified element with bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

enthalpy = engine.standardMolarEnthalpy(298.15, 101325, 0)
print(enthalpy)
    )doc")
        .def("standardMolarEnthalpy", static_cast<double (ChemicalEngineArray::*)(double, double, std::string, std::optional<std::string>) const>(&ChemicalEngineArray::standardMolarEnthalpy),
             py::arg("TK"), py::arg("P"), py::arg("species"), py::arg("phase") = py::none(),
             R"doc(
Returns the standard molar enthalpy of a specific species by its index (J/mol).

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param str species: Name of the species.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

enthalpy = engine.standardMolarEnthalpy(298.15, 101325, "SiO2")
print(enthalpy)
    )doc")
        .def("standardMolarEnthalpy_i", &ChemicalEngineArray::standardMolarEnthalpy_i,
             py::arg("TK"), py::arg("P"), py::arg("ispecies"),
             R"doc(
Returns the standard molar enthalpy of a specific species by its index (J/mol).
Access specified element without bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

enthalpy = engine.standardMolarEnthalpy_i(298.15, 101325, 0)
print(enthalpy)
    )doc")

        .def("standardMolarVolume", static_cast<double (ChemicalEngineArray::*)(double, double, size_t) const>(&ChemicalEngineArray::standardMolarVolume),
            py::arg("TK"), py::arg("P"), py::arg("ispecies"),
            R"doc(
Returns the standard molar volume of a species by its index (m³/mol).
Access specified element without bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

volume = engine.standardMolarVolume(298.15, 101325, 0)
print(volume)
    )doc")
        .def("standardMolarVolume", static_cast<double (ChemicalEngineArray::*)(double, double, std::string, std::optional<std::string>) const>(&ChemicalEngineArray::standardMolarVolume),
            py::arg("TK"), py::arg("P"), py::arg("species"), py::arg("phase") = py::none(),
      R"doc(
Returns the standard molar volume of a species by its index (m³/mol).

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param str species: Name of the species.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

volume = engine.standardMolarVolume(298.15, 101325, "SiO2")
print(volume)
    )doc")
        .def("standardMolarVolume_i", &ChemicalEngineArray::standardMolarVolume_i,
             py::arg("TK"), py::arg("P"), py::arg("ispecies"),
             R"doc(
Returns the standard molar volume of a species by its index (m³/mol).
Access specified element without bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

volume = engine.standardMolarVolume_i(298.15, 101325, 0)
print(volume)
    )doc")

        .def("standardMolarEntropy", static_cast<double (ChemicalEngineArray::*)(double, double, size_t) const>(&ChemicalEngineArray::standardMolarEntropy),
            py::arg("TK"), py::arg("P"), py::arg("ispecies"),
            R"doc(
Returns the standard molar entropy of a specific species by its index (J/K/mol).
Access specified element with bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

entropy = engine.standardMolarEntropy(298.15, 101325, 1)
print(entropy)
    )doc")
        .def("standardMolarEntropy", static_cast<double (ChemicalEngineArray::*)(double, double, std::string, std::optional<std::string>) const>(&ChemicalEngineArray::standardMolarEntropy),
            py::arg("TK"), py::arg("P"), py::arg("species"), py::arg("phase") = py::none(),
            R"doc(
Returns the standard molar entropy of a specific species by its index (J/K/mol).

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param str species: Name of the species.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

entropy = engine.standardMolarEntropy(298.15, 101325, "SiO2")
print(entropy)
    )doc")
        .def("standardMolarEntropy_i", &ChemicalEngineArray::standardMolarEntropy_i,
            py::arg("TK"), py::arg("P"), py::arg("ispecies"),
            R"doc(
Returns the standard molar entropy of a specific species by its index (J/K/mol).
Access specified element without bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

entropy = engine.standardMolarEntropy_i(298.15, 101325, 1)
print(entropy)
    )doc")

        .def("standardMolarHeatCapacityConstP",static_cast<double (ChemicalEngineArray::*)(double, double,  size_t) const>(&ChemicalEngineArray::standardMolarHeatCapacityConstP),
            py::arg("TK"), py::arg("P"), py::arg("ispecies"),
            R"doc(
Returns the standard molar heat capacity at constant pressure for species (J/K/mol).
Access specified element with bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

heat_capacity = engine.standardMolarHeatCapacityConstP(298.15, 101325, 0)
print(heat_capacity)
    )doc")
        .def("standardMolarHeatCapacityConstP", static_cast<double (ChemicalEngineArray::*)(double, double, std::string, std::optional<std::string>) const>(&ChemicalEngineArray::standardMolarHeatCapacityConstP),
            py::arg("TK"), py::arg("P"), py::arg("species"), py::arg("phase") = py::none(),
            R"doc(
Returns the standard molar heat capacity at constant pressure for species (J/K/mol).

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param str species: Name of the species.
:param str phase: Name of the phase the species was included in. If None get the first index.

**Example:**

.. code-block:: python

heat_capacity = engine.standardMolarHeatCapacityConstP(298.15, 101325, "SiO2")
print(heat_capacity)
    )doc")
        .def("standardMolarHeatCapacityConstP_i", &ChemicalEngineArray::standardMolarHeatCapacityConstP_i,
            py::arg("TK"), py::arg("P"), py::arg("ispecies"),
            R"doc(
Returns the standard molar heat capacity at constant pressure for species (J/K/mol).
Access specified element without bounds checking.

:param float T: Temperature in Kelvin (K).
:param float P: Pressure in Pascals (Pa).
:param int ispecies: Index of the species.

**Example:**

.. code-block:: python

heat_capacity = engine.standardMolarHeatCapacityConstP_i(298.15, 101325, 0)
print(heat_capacity)
    )doc");



            // .def("numElements", &ChemicalEngineArray::numElements)
            // .def("numSpecies", &ChemicalEngineArray::numSpecies)
            // .def("numPhases", &ChemicalEngineArray::numPhases)
            // .def("numSpeciesInPhase", &ChemicalEngineArray::numSpeciesInPhase)
            // .def("elementName", &ChemicalEngineArray::elementName)
            // .def("speciesName", &ChemicalEngineArray::speciesName)
            // .def("speciesCharge", &ChemicalEngineArray::speciesCharge)
            // .def("phaseName", &ChemicalEngineArray::phaseName)
            // .def("indexElement", &ChemicalEngineArray::indexElement)
            // .def("indexSpecies", &ChemicalEngineArray::indexSpecies)
            // .def("indexSpeciesAll", &ChemicalEngineArray::indexSpeciesAll)
            // .def("indexPhase", &ChemicalEngineArray::indexPhase)
            // .def("indexPhaseAll", &ChemicalEngineArray::indexPhaseAll)
            // .def("indexPhaseWithSpecies", &ChemicalEngineArray::indexPhaseWithSpecies)
            // .def("indexFirstSpeciesInPhase", &ChemicalEngineArray::indexFirstSpeciesInPhase)
            // .def("elementMolarMasses", &ChemicalEngineArray::elementMolarMasses, py::return_value_policy::reference_internal)
            // .def("speciesMolarMasses", &ChemicalEngineArray::speciesMolarMasses, py::return_value_policy::reference_internal)
            // .def("formulaMatrix", &ChemicalEngineArray::formulaMatrix, py::return_value_policy::reference_internal)
            // .def("setStandardMolarGibbsEnergy", &ChemicalEngineArray::setStandardMolarGibbsEnergy)
            // .def("standardMolarGibbsEnergy", &ChemicalEngineArray::standardMolarGibbsEnergy)
            // .def("standardMolarEnthalpy", &ChemicalEngineArray::standardMolarEnthalpy)
            // .def("standardMolarVolume", &ChemicalEngineArray::standardMolarVolume)
            // .def("standardMolarEntropy", &ChemicalEngineArray::standardMolarEntropy)
            // .def("standardMolarHeatCapacityConstP", &ChemicalEngineArray::standardMolarHeatCapacityConstP)
            // .def("setPT", &ChemicalEngineArray::setPT)
            // .def("setB", &ChemicalEngineArray::setB)
            // .def("setSpeciesAmounts", &ChemicalEngineArray::setSpeciesAmounts)
            // .def("setSpeciesAmount", py::overload_cast<long int, std::string, double>(&ChemicalEngineArray::setSpeciesAmount))
            // .def("setSpeciesAmount", py::overload_cast<long int, size_t, double>(&ChemicalEngineArray::setSpeciesAmount))
            // .def("setSpeciesUpperLimits", &ChemicalEngineArray::setSpeciesUpperLimits, py::return_value_policy::reference_internal)
            // .def("setSpeciesUpperLimit", py::overload_cast<long int, std::string, double>(&ChemicalEngineArray::setSpeciesUpperLimit))
            // .def("setSpeciesUpperLimit", py::overload_cast<long int, size_t, double>(&ChemicalEngineArray::setSpeciesUpperLimit))
            // .def("setSpeciesLowerLimits", &ChemicalEngineArray::setSpeciesLowerLimits, py::return_value_policy::reference_internal)
            // .def("setSpeciesLowerLimit", py::overload_cast<long int, std::string, double>(&ChemicalEngineArray::setSpeciesLowerLimit))
            // .def("setSpeciesLowerLimit", py::overload_cast<long int, size_t, double>(&ChemicalEngineArray::setSpeciesLowerLimit))
            // .def("converged", &ChemicalEngineArray::converged)
            // .def("numIterations", &ChemicalEngineArray::numIterations)
            // .def("temperature", &ChemicalEngineArray::temperature)
            // .def("pressure", &ChemicalEngineArray::pressure)
            // .def("elementAmounts", &ChemicalEngineArray::elementAmounts)
            // .def("elementAmountsInPhase", &ChemicalEngineArray::elementAmountsInPhase)
            // .def("elementAmountsInSpecies", &ChemicalEngineArray::elementAmountsInSpecies)
            // .def("speciesAmounts", &ChemicalEngineArray::speciesAmounts)
            // .def("speciesAmount", py::overload_cast<long int, size_t>(&ChemicalEngineArray::speciesAmount, py::const_))
            // .def("speciesAmount", py::overload_cast<long int, std::string>(&ChemicalEngineArray::speciesAmount, py::const_))
            // .def("speciesUpperLimits", &ChemicalEngineArray::speciesUpperLimits)
            // .def("speciesLowerLimits", &ChemicalEngineArray::speciesLowerLimits)
            // .def("systemMass", &ChemicalEngineArray::systemMass)
            // .def("systemVolume", &ChemicalEngineArray::systemVolume)
            // .def("ionicStrength", &ChemicalEngineArray::ionicStrength)
            // .def("pH", &ChemicalEngineArray::pH)
            // .def("pe", &ChemicalEngineArray::pe)
            // .def("Eh", &ChemicalEngineArray::Eh)
            // .def("systemGibbsEnergy", &ChemicalEngineArray::systemGibbsEnergy)
            // .def("phaseDensities", &ChemicalEngineArray::phaseDensities)
            // .def("phaseDensity", &ChemicalEngineArray::phaseDensity)
            // .def("phaseMasses", &ChemicalEngineArray::phaseMasses)
            // .def("phaseMass", &ChemicalEngineArray::phaseMass)
            // .def("phaseAmounts", &ChemicalEngineArray::phaseAmounts)
            // .def("phaseAmount", &ChemicalEngineArray::phaseAmount)
            // .def("phaseVolumes", &ChemicalEngineArray::phaseVolumes)
            // .def("phaseVolume", &ChemicalEngineArray::phaseVolume)
            // .def("phaseEnthalpies", &ChemicalEngineArray::phaseEnthalpies)
            // .def("phaseEnthalpy", &ChemicalEngineArray::phaseEnthalpy)
            // .def("phaseEntropies", &ChemicalEngineArray::phaseEntropies)
            // .def("phaseEntropy", &ChemicalEngineArray::phaseEntropy)
            // .def("phaseHeatCapacitiesConstP", &ChemicalEngineArray::phaseHeatCapacitiesConstP)
            // .def("phaseHeatCapacityConstP", &ChemicalEngineArray::phaseHeatCapacityConstP)
            // .def("phaseSatIndices", &ChemicalEngineArray::phaseSatIndices)
            // .def("phaseSatIndex", &ChemicalEngineArray::phaseSatIndex)
            // .def("phaseMolarGibbsEnergy", &ChemicalEngineArray::phaseMolarGibbsEnergy)
            // .def("phaseMolarEnthalpy", &ChemicalEngineArray::phaseMolarEnthalpy)
            // .def("phaseMolarVolume", &ChemicalEngineArray::phaseMolarVolume)
            // .def("phaseMolarEntropy", &ChemicalEngineArray::phaseMolarEntropy)
            // .def("moleFractions", &ChemicalEngineArray::moleFractions)
            // .def("lnActivityCoefficients", &ChemicalEngineArray::lnActivityCoefficients)
            // .def("lnActivities", &ChemicalEngineArray::lnActivities)
            // .def("lnConcentrations", &ChemicalEngineArray::lnConcentrations)
            // .def("chemicalPotentials", &ChemicalEngineArray::chemicalPotentials)
            // .def("phaseMolarHeatCapacityConstP", &ChemicalEngineArray::phaseMolarHeatCapacityConstP)
            // .def("phaseSpecificHeatCapacityConstP", &ChemicalEngineArray::phaseSpecificHeatCapacityConstP)
            // .def("phaseSpecificGibbsEnergy", &ChemicalEngineArray::phaseSpecificGibbsEnergy)
            // .def("phaseSpecificEnthalpy", &ChemicalEngineArray::phaseSpecificEnthalpy)
            // .def("phaseSpecificVolume", &ChemicalEngineArray::phaseSpecificVolume)
            // .def("phaseSpecificEntropy", &ChemicalEngineArray::phaseSpecificEntropy)
            // .def("systemEnthalpy", &ChemicalEngineArray::systemEnthalpy)
            // .def("systemEntropy", &ChemicalEngineArray::systemEntropy)
            // .def("systemHeatCapacityConstP", &ChemicalEngineArray::systemHeatCapacityConstP)
            // .def("__repr__", [](const ChemicalEngineArray& self) { std::stringstream ss; ss << self; return ss.str(); })
            // ;

    //   .def("speciesMolalities", &ChemicalEngineArray::speciesMolalities, py::return_value_policy::reference_internal)
    //   .def("standardMolarInternalEnergy", &ChemicalEngineArray::standardMolarInternalEnergy)
    //   .def("standardMolarHelmholtzEnergy", &ChemicalEngineArray::standardMolarHelmholtzEnergy)
    //   .def("standardMolarHeatCapacityConstV", &ChemicalEngineArray::standardMolarHeatCapacityConstV)
    //   .def("phaseMolarInternalEnergy", &ChemicalEngineArray::phaseMolarInternalEnergy)
    //   .def("phaseMolarHelmholtzEnergy", &ChemicalEngineArray::phaseMolarHelmholtzEnergy)
    //   .def("phaseMolarHeatCapacityConstV", &ChemicalEngineArray::phaseMolarHeatCapacitiesConstV)
    //   .def("phaseSpecificInternalEnergy", &ChemicalEngineArray::phaseSpecificInternalEnergy)
    //   .def("phaseSpecificHelmholtzEnergy", &ChemicalEngineArray::phaseSpecificHelmholtzEnergy)
    //   .def("phaseSpecificHeatCapacityConstV", &ChemicalEngineArray::phaseSpecificHeatCapacityConstVl)
}
