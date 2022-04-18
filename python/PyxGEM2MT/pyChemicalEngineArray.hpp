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

void exportChemicalEngineArray(py::module& m)
{
    py::class_<ChemicalEngineArrayOptions>(m, "ChemicalEngineArrayOptions")
            .def(py::init<>())
            .def_readwrite("warmstart", &ChemicalEngineArrayOptions::warmstart)
            .def_readwrite("useSIA", &ChemicalEngineArrayOptions::useSIA)
            .def_readwrite("step", &ChemicalEngineArrayOptions::step)
            .def_readwrite("cdv", &ChemicalEngineArrayOptions::cdv)
            .def_readwrite("cez", &ChemicalEngineArrayOptions::cez)
            .def("to_json_string", &ChemicalEngineArrayOptions::to_json_string, py::arg("dense")  = false)
            .def("__repr__", [](const ChemicalEngineArrayOptions& self)
    { std::stringstream ss; ss << self.to_json_string(); return ss.str(); })
            ;

    py::class_<ChemicalEngineArray>(m, "ChemicalEngineArray")
            .def(py::init<long int>())
            .def("setOptions", &ChemicalEngineArray::setOptions)
            .def("setWarmStart", &ChemicalEngineArray::setWarmStart)
            .def("setColdStart", &ChemicalEngineArray::setColdStart)
            .def("nNodes", &ChemicalEngineArray::nNodes)
            .def("initialize", &ChemicalEngineArray::initialize)
            .def("initializeFromJsonStrings", &ChemicalEngineArray::initializeFromJsonStrings)
            .def("reequilibrate", py::overload_cast<long int, bool>(&ChemicalEngineArray::reequilibrate))
            .def("reequilibrate", py::overload_cast<bool>(&ChemicalEngineArray::reequilibrate))
            .def("equilibrate", py::overload_cast<long int>(&ChemicalEngineArray::equilibrate))
            .def("equilibrate", py::overload_cast<>(&ChemicalEngineArray::equilibrate))
            .def("dumpGEMS3K", &ChemicalEngineArray::dumpGEMS3K, py::arg("filepath"), py::arg("put_all")=false)
            .def("nextTimePoint", &ChemicalEngineArray::nextTimePoint)
            .def("readDbrFromFile", &ChemicalEngineArray::readDbrFromFile)
            .def("readDbrFromJsonString", &ChemicalEngineArray::readDbrFromJsonString)
            .def("writeDbrToFile", &ChemicalEngineArray::writeDbrToFile)
            .def("writeDbrToJsonString", &ChemicalEngineArray::writeDbrToJsonString)
            .def("numElements", &ChemicalEngineArray::numElements)
            .def("numSpecies", &ChemicalEngineArray::numSpecies)
            .def("numPhases", &ChemicalEngineArray::numPhases)
            .def("numSpeciesInPhase", &ChemicalEngineArray::numSpeciesInPhase)
            .def("elementName", &ChemicalEngineArray::elementName)
            .def("speciesName", &ChemicalEngineArray::speciesName)
            .def("speciesCharge", &ChemicalEngineArray::speciesCharge)
            .def("phaseName", &ChemicalEngineArray::phaseName)
            .def("indexElement", &ChemicalEngineArray::indexElement)
            .def("indexSpecies", &ChemicalEngineArray::indexSpecies)
            .def("indexSpeciesAll", &ChemicalEngineArray::indexSpeciesAll)
            .def("indexPhase", &ChemicalEngineArray::indexPhase)
            .def("indexPhaseAll", &ChemicalEngineArray::indexPhaseAll)
            .def("indexPhaseWithSpecies", &ChemicalEngineArray::indexPhaseWithSpecies)
            .def("indexFirstSpeciesInPhase", &ChemicalEngineArray::indexFirstSpeciesInPhase)
            .def("elementMolarMasses", &ChemicalEngineArray::elementMolarMasses, py::return_value_policy::reference_internal)
            .def("speciesMolarMasses", &ChemicalEngineArray::speciesMolarMasses, py::return_value_policy::reference_internal)
            .def("formulaMatrix", &ChemicalEngineArray::formulaMatrix, py::return_value_policy::reference_internal)
            .def("setStandardMolarGibbsEnergy", &ChemicalEngineArray::setStandardMolarGibbsEnergy)
            .def("standardMolarGibbsEnergy", &ChemicalEngineArray::standardMolarGibbsEnergy)
            .def("standardMolarEnthalpy", &ChemicalEngineArray::standardMolarEnthalpy)
            .def("standardMolarVolume", &ChemicalEngineArray::standardMolarVolume)
            .def("standardMolarEntropy", &ChemicalEngineArray::standardMolarEntropy)
            .def("standardMolarHeatCapacityConstP", &ChemicalEngineArray::standardMolarHeatCapacityConstP)
            .def("setPT", &ChemicalEngineArray::setPT)
            .def("setB", &ChemicalEngineArray::setB)
            .def("setSpeciesAmounts", &ChemicalEngineArray::setSpeciesAmounts)
            .def("setSpeciesAmount", py::overload_cast<long int, std::string, double>(&ChemicalEngineArray::setSpeciesAmount))
            .def("setSpeciesAmount", py::overload_cast<long int, Index, double>(&ChemicalEngineArray::setSpeciesAmount))
            .def("setSpeciesUpperLimits", &ChemicalEngineArray::setSpeciesUpperLimits, py::return_value_policy::reference_internal)
            .def("setSpeciesUpperLimit", py::overload_cast<long int, std::string, double>(&ChemicalEngineArray::setSpeciesUpperLimit))
            .def("setSpeciesUpperLimit", py::overload_cast<long int, Index, double>(&ChemicalEngineArray::setSpeciesUpperLimit))
            .def("setSpeciesLowerLimits", &ChemicalEngineArray::setSpeciesLowerLimits, py::return_value_policy::reference_internal)
            .def("setSpeciesLowerLimit", py::overload_cast<long int, std::string, double>(&ChemicalEngineArray::setSpeciesLowerLimit))
            .def("setSpeciesLowerLimit", py::overload_cast<long int, Index, double>(&ChemicalEngineArray::setSpeciesLowerLimit))
            .def("converged", &ChemicalEngineArray::converged)
            .def("numIterations", &ChemicalEngineArray::numIterations)
            .def("temperature", &ChemicalEngineArray::temperature)
            .def("pressure", &ChemicalEngineArray::pressure)
            .def("elementAmounts", &ChemicalEngineArray::elementAmounts)
            .def("elementAmountsInPhase", &ChemicalEngineArray::elementAmountsInPhase)
            .def("elementAmountsInSpecies", &ChemicalEngineArray::elementAmountsInSpecies)
            .def("speciesAmounts", &ChemicalEngineArray::speciesAmounts)
            .def("speciesAmount", py::overload_cast<long int, Index>(&ChemicalEngineArray::speciesAmount, py::const_))
            .def("speciesAmount", py::overload_cast<long int, std::string>(&ChemicalEngineArray::speciesAmount, py::const_))
            .def("speciesUpperLimits", &ChemicalEngineArray::speciesUpperLimits)
            .def("speciesLowerLimits", &ChemicalEngineArray::speciesLowerLimits)
            .def("systemMass", &ChemicalEngineArray::systemMass)
            .def("systemVolume", &ChemicalEngineArray::systemVolume)
            .def("ionicStrength", &ChemicalEngineArray::ionicStrength)
            .def("pH", &ChemicalEngineArray::pH)
            .def("pe", &ChemicalEngineArray::pe)
            .def("Eh", &ChemicalEngineArray::Eh)
            .def("systemGibbsEnergy", &ChemicalEngineArray::systemGibbsEnergy)
            .def("phaseDensities", &ChemicalEngineArray::phaseDensities)
            .def("phaseDensity", &ChemicalEngineArray::phaseDensity)
            .def("phaseMasses", &ChemicalEngineArray::phaseMasses)
            .def("phaseMass", &ChemicalEngineArray::phaseMass)
            .def("phaseAmounts", &ChemicalEngineArray::phaseAmounts)
            .def("phaseAmount", &ChemicalEngineArray::phaseAmount)
            .def("phaseVolumes", &ChemicalEngineArray::phaseVolumes)
            .def("phaseVolume", &ChemicalEngineArray::phaseVolume)
            .def("phaseEnthalpies", &ChemicalEngineArray::phaseEnthalpies)
            .def("phaseEnthalpy", &ChemicalEngineArray::phaseEnthalpy)
            .def("phaseEntropies", &ChemicalEngineArray::phaseEntropies)
            .def("phaseEntropy", &ChemicalEngineArray::phaseEntropy)
            .def("phaseHeatCapacitiesConstP", &ChemicalEngineArray::phaseHeatCapacitiesConstP)
            .def("phaseHeatCapacityConstP", &ChemicalEngineArray::phaseHeatCapacityConstP)
            .def("phaseSatIndices", &ChemicalEngineArray::phaseSatIndices)
            .def("phaseSatIndex", &ChemicalEngineArray::phaseSatIndex)
            .def("phaseMolarGibbsEnergy", &ChemicalEngineArray::phaseMolarGibbsEnergy)
            .def("phaseMolarEnthalpy", &ChemicalEngineArray::phaseMolarEnthalpy)
            .def("phaseMolarVolume", &ChemicalEngineArray::phaseMolarVolume)
            .def("phaseMolarEntropy", &ChemicalEngineArray::phaseMolarEntropy)
            .def("moleFractions", &ChemicalEngineArray::moleFractions)
            .def("lnActivityCoefficients", &ChemicalEngineArray::lnActivityCoefficients)
            .def("lnActivities", &ChemicalEngineArray::lnActivities)
            .def("lnConcentrations", &ChemicalEngineArray::lnConcentrations)
            .def("chemicalPotentials", &ChemicalEngineArray::chemicalPotentials)
            .def("phaseMolarHeatCapacityConstP", &ChemicalEngineArray::phaseMolarHeatCapacityConstP)
            .def("phaseSpecificHeatCapacityConstP", &ChemicalEngineArray::phaseSpecificHeatCapacityConstP)
            .def("phaseSpecificGibbsEnergy", &ChemicalEngineArray::phaseSpecificGibbsEnergy)
            .def("phaseSpecificEnthalpy", &ChemicalEngineArray::phaseSpecificEnthalpy)
            .def("phaseSpecificVolume", &ChemicalEngineArray::phaseSpecificVolume)
            .def("phaseSpecificEntropy", &ChemicalEngineArray::phaseSpecificEntropy)
            .def("systemEnthalpy", &ChemicalEngineArray::systemEnthalpy)
            .def("systemEntropy", &ChemicalEngineArray::systemEntropy)
            .def("systemHeatCapacityConstP", &ChemicalEngineArray::systemHeatCapacityConstP)
            .def("__repr__", [](const ChemicalEngineArray& self) { std::stringstream ss; ss << self; return ss.str(); })
            ;

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
