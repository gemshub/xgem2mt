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

:param int ndx: Index of the node (system) in the list.
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

:param int ndx: Index of the node (system) in the list.
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

 :param int ndx: Index of the node (system) in the list.
 :param str filename: Path to the input DBR file (e.g., "*-dbr.json/.dat").

**Example:**

.. code-block:: python

engine.readDbrFromFile(0, "my-system-dbr-0-0000.json")
)doc")
            .def("readDbrFromJsonString", &ChemicalEngineArray::readDbrFromJsonString,
             R"doc(
Reads System DBR composition from a JSON string to the node.

Loads the structure of node-dependent data (composition, T, P and other).

 :param int ndx: Index of the node (system) in the list.
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

:param int ndx: Index of the node (system) in the list.
:param str filename: Path to the output DBR file (e.g., "my-system-dbr-0-0001.dat").

**Example:**

.. code-block:: python

engine.writeDbrToFile(1, "output.dbr")
)doc")
            .def("writeDbrToJsonString", &ChemicalEngineArray::writeDbrToJsonString,
             R"doc(
Returns the current DBR as a JSON string.

:param int ndx: Index of the node (system) in the list.
:return: JSON string representing the system composition.
:rtype: str

**Example:**

.. code-block:: python

json_state = engine.writeDbrToJsonString(1)
print(json_state)
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
            // .def("setSpeciesAmount", py::overload_cast<long int, Index, double>(&ChemicalEngineArray::setSpeciesAmount))
            // .def("setSpeciesUpperLimits", &ChemicalEngineArray::setSpeciesUpperLimits, py::return_value_policy::reference_internal)
            // .def("setSpeciesUpperLimit", py::overload_cast<long int, std::string, double>(&ChemicalEngineArray::setSpeciesUpperLimit))
            // .def("setSpeciesUpperLimit", py::overload_cast<long int, Index, double>(&ChemicalEngineArray::setSpeciesUpperLimit))
            // .def("setSpeciesLowerLimits", &ChemicalEngineArray::setSpeciesLowerLimits, py::return_value_policy::reference_internal)
            // .def("setSpeciesLowerLimit", py::overload_cast<long int, std::string, double>(&ChemicalEngineArray::setSpeciesLowerLimit))
            // .def("setSpeciesLowerLimit", py::overload_cast<long int, Index, double>(&ChemicalEngineArray::setSpeciesLowerLimit))
            // .def("converged", &ChemicalEngineArray::converged)
            // .def("numIterations", &ChemicalEngineArray::numIterations)
            // .def("temperature", &ChemicalEngineArray::temperature)
            // .def("pressure", &ChemicalEngineArray::pressure)
            // .def("elementAmounts", &ChemicalEngineArray::elementAmounts)
            // .def("elementAmountsInPhase", &ChemicalEngineArray::elementAmountsInPhase)
            // .def("elementAmountsInSpecies", &ChemicalEngineArray::elementAmountsInSpecies)
            // .def("speciesAmounts", &ChemicalEngineArray::speciesAmounts)
            // .def("speciesAmount", py::overload_cast<long int, Index>(&ChemicalEngineArray::speciesAmount, py::const_))
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

    //   .def("speciesMolalities", &ChemicalEngine::speciesMolalities, py::return_value_policy::reference_internal)
    //   .def("standardMolarInternalEnergy", &ChemicalEngine::standardMolarInternalEnergy)
    //   .def("standardMolarHelmholtzEnergy", &ChemicalEngine::standardMolarHelmholtzEnergy)
    //   .def("standardMolarHeatCapacityConstV", &ChemicalEngine::standardMolarHeatCapacityConstV)
    //   .def("phaseMolarInternalEnergy", &ChemicalEngine::phaseMolarInternalEnergy)
    //   .def("phaseMolarHelmholtzEnergy", &ChemicalEngine::phaseMolarHelmholtzEnergy)
    //   .def("phaseMolarHeatCapacityConstV", &ChemicalEngine::phaseMolarHeatCapacitiesConstV)
    //   .def("phaseSpecificInternalEnergy", &ChemicalEngine::phaseSpecificInternalEnergy)
    //   .def("phaseSpecificHelmholtzEnergy", &ChemicalEngine::phaseSpecificHelmholtzEnergy)
    //   .def("phaseSpecificHeatCapacityConstV", &ChemicalEngine::phaseSpecificHeatCapacityConstVl)
}
