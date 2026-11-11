// xGEM2MT is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
//
// Copyright (C) 2022-2026 S.Dmytriieva, Dmtrii Kulik, G.D. Miron
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

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <limits>
#include <map>
#include <optional>

#include "xGEM2MT/MassTransport.hpp"

namespace py = pybind11;
using namespace xGEM2MT;


// Node index or list of node indices
inline std::vector<size_t> toNodes(py::object nodes)
{
    if(py::isinstance<py::int_>(nodes))
        return {nodes.cast<size_t>()};
    return nodes.cast<std::vector<size_t>>();
}

// None, an xGEMS Material (anything with b_dict()), a dict or a sequence of element amounts
// None, an xGEMS Material (anything with b_dict()), a dict or a sequence of element amounts
inline Composition toComposition(py::object c)
{
    if(c.is_none())
        return std::monostate{};
    if(py::hasattr(c, "b_dict"))
        c = c.attr("b_dict")();
    if(py::isinstance<py::dict>(c))
        return c.cast<std::map<std::string, double>>();
    return c.cast<std::vector<double>>();
}

inline NodeSetup toSetup(py::object composition, std::optional<double> T, std::optional<double> P,
                         std::optional<TransportProperties> transport,
                         py::object pore_water = py::none(), py::object solids = py::none(),
                         std::optional<double> grain_density = std::nullopt,
                         std::optional<double> water_mass = std::nullopt,
                         std::optional<double> water_saturation = std::nullopt,
                         std::vector<Part> parts = {})
{
    NodeSetup setup;
    setup.composition = toComposition(composition);
    setup.temperature = T;
    setup.pressure = P;
    if(transport)
        setup.transport = *transport;
    setup.pore_water = toComposition(pore_water);
    setup.solids = toComposition(solids);
    setup.grain_density = grain_density;
    setup.water_mass = water_mass;
    setup.water_saturation = water_saturation;
    setup.parts = std::move(parts);
    return setup;
}

void exportTransportProperties(py::module& m)
{
    py::class_<Part>(m, "Part",
        R"doc(
        A part of a node for ``setNodes``/``setSystem(..., parts=[...])``: a phase or a solid, given as a
        material and the share of the node that it fills.

        ``of`` says which volume the part fills: ``"pore"`` (a gas, another liquid) or ``"solid"``. Give either
        ``fraction`` (a share of that volume) or ``mass``. The material can have any amount: it is scaled. The
        volume of the part comes from a GEM calculation at the T and P of the node (the volume that it adds to
        the pore water, with the models of the chemical system, e.g. a non-ideal gas), unless ``density`` is
        given.

        :param material: an xGEMS ``Material``, a dict ``{element: mol}`` or a list of element amounts.
        :param float fraction: share of the ``of`` volume, 0 to 1. Pore: water_saturation and the pore
            fractions add up to 1 or less. Solid: the solid fractions are shares of the solid volume that the
            parts given by mass leave, and they add up to 1.
        :param float mass: kg of the part in the node. Its volume is taken from the ``of`` volume.
        :param str of: ``"pore"`` or ``"solid"``.
        :param float density: kg/L; if not given, from GEM.

        **Example:**

        .. code-block:: python

            task.setSystem(1, pore_water=water, water_saturation=0.9,
                           parts=[Part(air, fraction=0.1, of="pore"),
                                  Part(quartz, fraction=0.99, of="solid"),
                                  Part(calcite, fraction=0.01, of="solid", density=2.71),
                                  Part(organic, mass=0.005, of="solid")])
        )doc")
        .def(py::init([](py::object material, std::optional<double> fraction, std::optional<double> mass,
                         std::optional<std::string> of, std::optional<double> density) {
                 Part pt;
                 pt.material = toComposition(material);
                 pt.fraction = fraction; pt.mass = mass; pt.of = of.value_or(""); pt.density = density;
                 return pt;
             }),
             py::arg("material"), py::kw_only(), py::arg("fraction") = py::none(), py::arg("mass") = py::none(),
             py::arg("of") = py::none(), py::arg("density") = py::none())
        .def_readwrite("fraction", &Part::fraction)
        .def_readwrite("mass", &Part::mass)
        .def_readwrite("of", &Part::of)
        .def_readwrite("density", &Part::density)
        .def("__repr__", [](const Part& p) {
            std::string r = "Part(";
            if(p.fraction) r += "fraction=" + py::repr(py::float_(*p.fraction)).cast<std::string>();
            if(p.mass) r += "mass=" + py::repr(py::float_(*p.mass)).cast<std::string>();
            r += ", of=\"" + p.of + "\"";
            if(p.density) r += ", density=" + py::repr(py::float_(*p.density)).cast<std::string>();
            return r + ")";
        });

    py::class_<GEMFailures>(m, "GEMFailures",
        "GEM calculations in the nodes that failed or were not trustful in the last run.")
        .def_readonly("ran", &GEMFailures::ran)
        .def_readonly("count", &GEMFailures::count)
        .def_readonly("nodes", &GEMFailures::nodes)
        .def_readonly("first_step", &GEMFailures::first_step)
        .def_readonly("first_time", &GEMFailures::first_time)
        .def_readonly("first_node", &GEMFailures::first_node)
        .def_readonly("first_message", &GEMFailures::first_message)
        .def("__bool__", [](const GEMFailures& f) { return f.count > 0; })
        .def("__repr__", [](const GEMFailures& f) {
            if(!f.ran) return std::string("GEMFailures(no run yet)");
            return "GEMFailures(count=" + std::to_string(f.count) + ", nodes=" + std::to_string(f.nodes.size()) +
                   (f.count ? ", first: step " + std::to_string(f.first_step) + " node " + std::to_string(f.first_node) +
                              ": " + f.first_message : std::string()) + ")";
        });

    py::class_<TransportProperties>(m, "TransportProperties",
        R"doc(
        Physical transport properties of a node. A property that is None is not changed.

        In mode 'A', ``velocity`` changes only the dispersion of the node (advection uses
        ``setFluidAdvectionVelocity``); ``volume`` and ``permeability`` are not used.
        Dispersion coefficient: D = dispersivity * velocity + diffusivity * porosity / tortuosity.

        :param float volume: total node volume, m³
        :param float velocity: fluid advection velocity, m/s
        :param float porosity: effective porosity, 0 to 1
        :param float permeability: effective permeability, m²
        :param float dispersivity: longitudinal dispersivity, m
        :param float diffusivity: diffusion coefficient in water, m²/s
        :param float tortuosity: tortuosity factor, 1 to 10
        )doc")
        .def(py::init([](std::optional<double> volume, std::optional<double> velocity,
                         std::optional<double> porosity, std::optional<double> permeability,
                         std::optional<double> dispersivity, std::optional<double> diffusivity,
                         std::optional<double> tortuosity) {
                 return TransportProperties{volume, velocity, porosity, permeability,
                                            dispersivity, diffusivity, tortuosity};
             }),
             py::kw_only(), py::arg("volume") = py::none(), py::arg("velocity") = py::none(),
             py::arg("porosity") = py::none(), py::arg("permeability") = py::none(),
             py::arg("dispersivity") = py::none(), py::arg("diffusivity") = py::none(),
             py::arg("tortuosity") = py::none())
        .def_readwrite("volume", &TransportProperties::volume)
        .def_readwrite("velocity", &TransportProperties::velocity)
        .def_readwrite("porosity", &TransportProperties::porosity)
        .def_readwrite("permeability", &TransportProperties::permeability)
        .def_readwrite("dispersivity", &TransportProperties::dispersivity)
        .def_readwrite("diffusivity", &TransportProperties::diffusivity)
        .def_readwrite("tortuosity", &TransportProperties::tortuosity)
        .def("__repr__", [](const TransportProperties& t) {
            std::string r = "TransportProperties(";
            const std::pair<const char*, const std::optional<double>*> f[] = {
                {"volume", &t.volume}, {"velocity", &t.velocity}, {"porosity", &t.porosity},
                {"permeability", &t.permeability}, {"dispersivity", &t.dispersivity},
                {"diffusivity", &t.diffusivity}, {"tortuosity", &t.tortuosity}};
            bool first = true;
            for(const auto& [name, v]: f)
                if(*v) { r += (first ? "" : ", ") + std::string(name) + "=" + py::repr(py::float_(**v)).cast<std::string>(); first = false; }
            return r + ")";
        });
}

void exportMassTransport(py::module& m)
{
    exportTransportProperties(m);
    py::class_<MassTransport>(m, "MassTransport",
                              R"doc(
        Mass Transport simulator for simple reactive mass transport models coupled with GEM
        calculation of equilibrium states in spatially distributed nodes (boxes, volumes) over time steps.

        The GEM2MT module can simulate reactive transport in three main modes:
        - **1. Sequential reactor chain and waves** (``'S'``)
        - **2. Flow-through and box-flux sequence** (``'F'``, ``'B'``)
        - **3. One-dimensional reactive transport** with advection/dispersion/diffusion (``'A'``, ``'C'``, ``'W'``)
        )doc")

        .def(py::init<char, size_t>(),
             py::arg("mode"),
             py::arg("nnodes"),
             R"doc(
             Constructs a MassTransport object.

             Creates the GEM2MT simulator for simple reactive mass transport models coupled with GEM
             calculation of equilibrium states in spatially distributed nodes (boxes, volumes) over time steps.
             The GEM2MT simulator can be configured by choosing its mode of operation and setting up the number
             of nodes (boxes, volumes).

             :param str mode: Code of GEM2MT mode of operation { S F B A C W }:
                              - ``S``: Sequential reactors chain with waves (moving fluids or solids between reactors)
                              - ``F``: Flow-through reactors chain with fluxes (simultaneous transport then equilibration)
                              - ``B``: Box-flux megasystem (with time-step integration of arbitrary MGP fluxes)
                              - ``A``: 1-D (advection-dispersion) reactive transport using classical finite-difference method
                              - ``C``: 1-D (advection-dispersion) reactive transport using implicit Crank-Nicolson scheme
                              - ``W``: 1-D (advection-diffusion) reactive transport using random-walk method
             :param int nnodes: Input number of local equilibrium cells (nodes). Typically, 50 to 250 nodes
                                provide relatively good spatial resolution and a reasonable computing time.

             **Example:**

             .. code-block:: python

                 task = MassTransport('F', 100)
             )doc")

        .def("getChemicalEngineArray", &MassTransport::getChemicalEngineArray,
             R"doc(
             Returns a shared pointer to the underlying ChemicalEngineArray instance.
             )doc")

        .def("readTask", &MassTransport::readTask,
             py::arg("gem2mt_filename"),
             R"doc(
             Initializes the MassTransport object from a GEM-Selektor generated file.

             :param str gem2mt_filename: The path to the file containing the definition of the gem2mt task.
             :raises RuntimeError: If there was a file read error during task initialization.
             )doc")

        .def("readTaskString", &MassTransport::readTaskString,
             py::arg("gem2mt_json"),
             R"doc(
             Initializes the MassTransport object from a JSON string.

             :param str gem2mt_json: A valid JSON string containing the definition of the gem2mt task.
             :raises RuntimeError: If there was an error during task initialization.
             )doc")

        .def("dump", &MassTransport::dump,
             R"doc(
             Dumps the MassTransport object state to a JSON string.

             :return str: A formatted JSON string representing the current state.
             )doc")

        .def("dumpTask", &MassTransport::dumpTask,
             py::arg("gem2mt_filename"),
             py::arg("ipm_filename") = std::string{},
             R"doc(
             Dumps the MassTransport object state and task configuration to files.

             :param str gem2mt_filename: The path where the main GEM2MT task file will be saved.
             :param str ipm_filename: The path where the IPM work structure file will be saved (optional).
             )doc")

        .def("initNodes", py::overload_cast<const std::string&>(&MassTransport::initNodes),
             py::arg("ipm_filename"),
             R"doc(
             Initializes all nodes from one GEMS3K ``*-dat.lst`` file.

             The initial systems are the dbr files listed in the ``*-dat.lst`` file (after the dch and ipm
             files, separated by spaces), in their order: the first file is system 0. A ``*-dbr.lst`` file is
             not necessary.

             .. note::
                 GEMS3K reads the initial systems only from a list file. This function writes a temporary
                 hidden list file in the folder of ``ipm_filename`` and removes it after reading. The folder
                 must be writable; if it is not, use the two-file form.

             :param str ipm_filename: Path to the ``*-dat.lst`` file.
             :raises RuntimeError: If a file cannot be read, or if no dbr files are listed.

             **Example:**

             .. code-block:: python

                 task.initNodes("my-system-dat.lst")
             )doc")

        .def("initNodes", py::overload_cast<const std::string&, const std::string&>(&MassTransport::initNodes),
             py::arg("ipm_filename"),
             py::arg("dbr_filename"),
             R"doc(
             Initializes all nodes (boxes) from GEM-Selektor IPM files.

             Reads the specified files and configures the chemical system defined inside them.

             :param str ipm_filename: Path to the system definition file (e.g., "my-system-dat.lst").
             :param str dbr_filename: Path to the file containing the list of system definition files (e.g., "my-system-dbr.lst").
             :raises RuntimeError: If there was a file read error during node initialization.

             **Example:**

             .. code-block:: python

                 task = MassTransport('F', 100)
                 task.initNodes("my-system-dat.lst", "my-system-dbr.lst")
             )doc")

        .def("initNodesFromStrings", &MassTransport::initNodesFromStrings,
             py::arg("dch_json"),
             py::arg("ipm_json"),
             py::arg("dbr_json"),
             R"doc(
             Initializes all nodes (boxes) using JSON strings.

             Configures the engine from three JSON inputs:
             - dch_json: Defines the chemical system.
             - ipm_json: Contains calculation parameters and settings for the GEMS3K IPM-3 algorithm.
             - dbr_json: Specifies the input node composition.
             These strings can be obtained from the GEM-Selektor exported GEMS3K system files as JSON (*-dch.json, *-ipm.json, *-dbr.json).

             :param str dch_json: JSON string describing the chemical system.
             :param str ipm_json: JSON string with parameter and algorithm settings.
             :param list[str] dbr_json: A vector of JSON strings containing node composition details.
             :raises RuntimeError: If there was an error during node initialization.

             **Example:**

             .. code-block:: python

                 dch = "{...}"
                 ipm = "{...}"
                 dbr = ["{...}", "{...}"]
                 engine = MassTransport('F', 51)
                 engine.initNodesFromStrings(dch, ipm, dbr)
             )doc")

        .def("runMassTransport",
             [](MassTransport& self) {
                 self.runMassTransport();
                 const auto f = self.gemFailures();
                 if(f.count > 0) {
                     std::string msg = std::to_string(f.count) + " GEM calculations failed or are not trustful in " +
                                       std::to_string(f.nodes.size()) + " nodes (first: step " + std::to_string(f.first_step) +
                                       ", node " + std::to_string(f.first_node) + ": " + f.first_message +
                                       "). See task.gemFailures() and xGEM2MT.log.";
                     py::module_::import("warnings").attr("warn")(msg, py::module_::import("builtins").attr("RuntimeWarning"), 2);
                 }
             },
             R"doc(
             Executes the coupled finite-difference mass-transport simulation loop.

             This method runs the main simulation routine, synchronizing the 1-D finite-difference
             mass transport with the GEM thermodynamic equilibrium calculations across all nodes.

             .. note::
                 Before calling this method, the simulator must be properly initialized.

             :raises RuntimeError: If the simulation fails to converge or hits a numerical instability.

             If GEM calculations in nodes fail or are not trustful, the run goes on (the node keeps its
             previous result) and Python gives a ``RuntimeWarning``; see ``gemFailures()``.
             )doc")

        .def("gemFailures", &MassTransport::gemFailures,
             R"doc(
             The GEM calculations that failed or were not trustful in the last run.

             :return GEMFailures: ``ran`` (False before the first run: then there is no result), ``count``
                 (0: all converged), ``nodes``, and ``first_step``, ``first_time``, ``first_node``,
                 ``first_message`` of the first failure. The calculations of ``initNodes`` and ``setNodes``
                 are not counted.
             )doc")

        .def_property_readonly("nNodes", &MassTransport::nNodes,
                               R"doc(
             Gets the total number of nodes (boxes).

             The number of nodes (boxes, reservoirs) is used for spatial discretization of
             the reactive transport problem (usually a 1-D column or a chain of boxes).

             :return int: The total number of nodes.

             **Example:**

             .. code-block:: python

                 total_nodes = task.nNodes()
             )doc")

        .def("name", &MassTransport::name,
             R"doc(
             Gets the full name of this GEM2MT task.

             :return str: The full name of the task.

             **Example:**

             .. code-block:: python

                 task_name = task.name()
             )doc")

        .def("setName", &MassTransport::setName,
             py::arg("taskname"),
             R"doc(
             Sets the full name of this GEM2MT task.

             :param str taskname: The name to assign to the GEM2MT task.

             **Example:**

             .. code-block:: python

                 task.setName("My_GEM2MT_Simulation")
             )doc")

        .def("comment", &MassTransport::comment,
             R"doc(
             Gets the optional enhancement comment for this GEM2MT task.

             :return str: The comment string.

             **Example:**

             .. code-block:: python

                 task_comment = task.comment()
             )doc")

        .def("setComment", &MassTransport::setComment,
             py::arg("comment"),
             R"doc(
             Sets the optional comment for this GEM2MT task.

             :param str comment: The comment text to assign to the task.

             **Example:**

             .. code-block:: python

                 task.setComment("Enable comment.")
             )doc")

        .def("setDistributing", &MassTransport::setDistributing,
             py::arg("node_ndx"), py::arg("system_ndx"),
             R"doc(
             Changes the initialization of the RT nodes (boxes).

             Each node (box) contains the index of its initial system composition variant.
             This method updates the index of the initial system variant for the specified node.

             .. note::
                 **Default distribution behavior:**
                 Assigns different fluid and rock composition indices to nodes (assuming fluid=0 and rock=1).
                 For ``A`` or ``C`` modes, nodes 0 and 1 contain the initial fluid. For other modes, node 0
                 contains initial 0, and all subsequent nodes contain 1.

             :param int node_ndx: The index of the target node.
             :param int system_ndx: The index of the initial system composition variant (from the GEMS3K .dbr file).

             **Example:**

             .. code-block:: python

                 # Set node 5 to use system composition variant 2
                 task.setDistributing(5, 2)
             )doc")

        .def("setNodeType", &MassTransport::setNodeType,
             py::arg("node_ndx"), py::arg("type"),
             R"doc(
             Changes the initialization type for nodes (boxes).

             This method defines the boundary condition type of the node:

             - **0**: Normal node
             - **1**: Dirichlet source (constant composition source)
             - **-1**: Dirichlet sink
             - **2**: Neumann source (constant gradient source)
             - **-2**: Neumann sink
             - **3**: Cauchy source (constant flux source)
             - **-3**: Cauchy sink
             - **4**: Input time-dependent function

             .. note::
                 **Default types:** 3 for node 0, -3 for the last node, and 0 for all other nodes
                 (except in ``F`` and ``B`` modes).

             :param int node_ndx: The index of the target node.
             :param int type: The boundary condition type to assign to the node.

             **Example:**

             .. code-block:: python

                 # Set node 0 as a Cauchy source (constant flux)
                 task.setNodeType(0, 3)
             )doc")

        .def("setSIA", &MassTransport::setSIA,
             py::arg("flag"),
             R"doc(
             Sets the initial approximation mode for the GEM IPM solver.

             Activating a smart initial approximation for GEM calculations at each node and
             time step can significantly reduce computational time. Deactivate this mode if
             error messages appear, or to verify results against the automatic initial approximation.

             The available modes are:

             - ``'+'``: Smart Initial Approximation (SIA)
             - ``'*'``: Internal Smart Initial Approximation (SIA internal)
             - ``'-'``: Automatic Initial Approximation (AIA) — takes longer but is more robust

             :param str flag: The character code indicating the desired mode ('+', '*', or '-').

             **Example:**

             .. code-block:: python

                 # Activate Smart Initial Approximation
                 task.setSIA('+')
             )doc")

        .def("setOutput", &MassTransport::setOutput,
             py::arg("enable"),
             R"doc(
             Enables or disables continuous debug output for nodes.

             .. note::
                 By default, this property is set to ``True``.

             :param bool enable: If ``True``, enables continuous debug logging; if ``False``, disables it.

             **Example:**

             .. code-block:: python

                 # Disable continuous debug output
                 task.setOutput(False)
             )doc")

        .def("setOutVTK", &MassTransport::setOutVTK,
             py::arg("enable"),
             R"doc(
             Enables or disables continuous debug output in VTK format.

             .. note::
                 By default, this property is set to ``False``.

             :param bool enable: If ``True``, enables continuous VTK output for nodes; if ``False``, disables it.

             **Example:**

             .. code-block:: python

                 # Enable continuous output to VTK format
                 task.setOutVTK(True)
             )doc")

        .def("setTypeFluxPhase", &MassTransport::setTypeFluxPhase,
             py::arg("type"),
             R"doc(
             Sets the type of MGP flux phase.

             Used in default initialization for MGP box-flux transport models.

             The available flux phase types are:

             - ``'0'``: Undefined
             - ``'1'``: Aqueous (aq)
             - ``'2'``: Gas
             - ``'3'``: Aqueous + Gas (aq+gas)
             - ``'4'``: Solids

             .. note::
                 By default, this property is set to ``1``.

             :param str type: The character code indicating the desired MGP flux phase type.

             **Example:**

             .. code-block:: python

                 # Set flux phase type to Gas ('2')
                 task.setTypeFluxPhase('2')
             )doc")

        .def("setTau", &MassTransport::setTau,
             py::arg("start"), py::arg("stop"), py::arg("step"),
             R"doc(
             Sets the physical time span and step size for RT models.

             The unit of time measurement is arbitrary (defined by the user), but all
             transport parameters (such as fluxes and velocities) must be consistent
             with the chosen time unit (e.g., seconds, hours, or days).

             The maximum number of time steps is automatically calculated from this time span.

             .. note::
                 **Default values are:**
                 - ``start`` = 0.0
                 - ``stop`` = 1000.0
                 - ``step`` = 1.0

             :param float start: The initial simulation time.
             :param float stop: The final simulation time.
             :param float step: The size of each time step.

             **Example:**

             .. code-block:: python

                 # Set simulation from t=0 to t=100 with a step size of 0.5
                 task.setTau(0.0, 100.0, 0.5)
             )doc")

        .def("setSpatialDimensions", &MassTransport::setSpatialDimensions,
             py::arg("x"), py::arg("y") = 0.0, py::arg("z") = 0.0,
             R"doc(
             Sets the spatial dimensions of the porous medium.

             The length ``x`` is used to calculate the spatial discretization interval:
             **dx = x / nC** (or node size), and to compute the linear pressure
             gradient (hydraulic head) as **(P_last - P_first) / x**.

             .. note::
                 The hydraulic head is not used directly to calculate advection velocity.
                 However, to maintain physical realism, a linear pressure gradient (from a few
                 millibars to a few bars) can be set in the initialization script. In this case,
                 each node will operate under a different pressure.

             .. warning::
                 Heat transport is not yet implemented in GEM2MT; all nodes are assumed
                 to be at the same temperature.

             Reasonable values for ``x`` range between 0.01 and 10.0 meters.
             By default, this property is set to ``x = 1.0``, ``y = 0.0``, ``z = 0.0``.

             :param float x: The spatial dimension along the X-axis (in meters).
             :param float y: [unused] Spatial dimension along the Y-axis. Reserved for future 2D extensions.
             :param float z: [unused] Spatial dimension along the Z-axis. Reserved for future 3D extensions.

             **Example:**

             .. code-block:: python

                 # Set a 1D column of 5.5 meters
                 task.setSpatialDimensions(5.5)
             )doc")

        .def("setTimeStepReductionFactor", &MassTransport::setTimeStepReductionFactor,
             py::arg("val"),
             R"doc(
             Sets the time-step reduction factor for advection, diffusion, and box-flux transport.

             .. note::
                 By default, this factor is typically set to ``1.0``.

             :param float val: The time-step reduction factor. Passing ``0.0`` resets the factor to ``1.0``.

             **Example:**

             .. code-block:: python

                 # Reset to default behavior
                 task.setTimeStepReductionFactor(0.0)
             )doc")

        .def("setInitialTotalNodeVolume", &MassTransport::setInitialTotalNodeVolume,
             py::arg("val"),
             R"doc(
             Sets the initial total node volume (in m³).

             This volume value is strictly used within the context of porosity-permeability relationships.

             .. note::
                 This value will be automatically copied to all nodes during the RT problem's
                 default initialization phase.

             :param float val: The initial total node volume in cubic meters (m³).

             **Example:**

             .. code-block:: python

                 # Set the initial volume of each node to 2.5 cubic meters
                 task.setInitialTotalNodeVolume(2.5)
             )doc")

        .def("setFluidAdvectionVelocity", &MassTransport::setFluidAdvectionVelocity,
             py::arg("val"),
             R"doc(
             Sets the fluid advection velocity (in m/s).

             Typically, the velocity value falls within the range: 10^-10 < v_p < 10^-3 m/s.

             .. note::
                 This value will be automatically copied to all nodes during the RT problem's
                 default initialization phase.

             :param float val: The fluid advection velocity value in meters per second (m/s).

             **Example:**

             .. code-block:: python

                 # Set advection velocity to 1e-5 m/s
                 task.setFluidAdvectionVelocity(1e-5)
             )doc")

        .def("setInitialNodeEffectivePorosity", &MassTransport::setInitialNodeEffectivePorosity,
             py::arg("val"),
             R"doc(
             Sets the initial effective porosity for nodes.

             Typically, the porosity value falls within the range: 0 < ε < 1.

             .. note::
                 - By default, this property is set to ``1.0``.
                 - This value will be automatically copied to all nodes during the RT problem's
                   default initialization phase.

             :param float val: The initial effective porosity value (dimensionless fraction).

             **Example:**

             .. code-block:: python

                 # Set the initial effective porosity to 25% (0.25)
                 task.setInitialNodeEffectivePorosity(0.25)
             )doc")

        .def("setInitialEffectivePermeability", &MassTransport::setInitialEffectivePermeability,
             py::arg("val"),
             R"doc(
             Sets the initial effective permeability (in m²).

             .. note::
                 - **Unit conversion:** 1 Darcy is approximately equal to 10^-12 m².
                 - This value will be automatically copied to all nodes during the RT problem's
                   default initialization phase.

             :param float val: The initial effective permeability value in square meters (m²).

             **Example:**

             .. code-block:: python

                 # Set the initial effective permeability to 10 Darcy (1e-11 m²)
                 task.setInitialEffectivePermeability(1e-11)
             )doc")

        .def("setInitialDispersivity", &MassTransport::setInitialDispersivity,
             py::arg("val"),
             R"doc(
             Sets the initial value of specific longitudinal dispersivity (in meters).

             .. note::
                 - **Typical values:** Range from 0.001 to 1.0 m.
                 - This value will be automatically copied to all nodes during the RT problem's
                   default initialization phase.

             :param float val: The initial specific longitudinal dispersivity value (in meters).

             **Example:**

             .. code-block:: python

                 # Set the initial longitudinal dispersivity to 0.1 meters
                 task.setInitialDispersivity(0.1)
             )doc")

        .def("setInitialDiffusivity", &MassTransport::setInitialDiffusivity,
             py::arg("val"),
             R"doc(
             Sets the initial general aqueous medium diffusivity (in m²/s).

             .. note::
                 - **Typical value:** Approximately 2 x 10^-9 m²/s.
                 - This value will be automatically copied to all nodes during the RT problem's
                   default initialization phase.

             :param float val: The initial general aqueous medium diffusivity value (in m²/s).

             **Example:**

             .. code-block:: python

                 # Set diffusivity to the typical value for aqueous ions
                 task.setInitialDiffusivity(2e-9)
             )doc")

        .def("setInitialTortuosityFactor", &MassTransport::setInitialTortuosityFactor,
             py::arg("val"),
             R"doc(
             Sets the initial tortuosity factor for nodes.

             Typically, the tortuosity factor value falls within the range: 1 <= τ_nto < 10.

             .. note::
                 - This value will be automatically copied to all nodes during the RT problem's
                   default initialization phase.

             :param float val: The initial tortuosity factor value (dimensionless factor).

             **Example:**

             .. code-block:: python

                 # Set the initial tortuosity factor to a typical value of 1.5
                 task.setInitialTortuosityFactor(1.5)
             )doc")

        .def("setCutoffICamount", &MassTransport::setCutoffICamount,
             py::arg("val"),
             R"doc(
             Sets the cutoff threshold for Independent Component (IC) amount differences between time steps.

             This cutoff is a performance optimization feature. If changes in the node bulk composition
             (at constant temperature and pressure, and in the absence of kinetics) after a time step
             are negligibly small, the chemical equilibrium state calculation is skipped.

             .. note::
                 - **Typical value:** Approximately 10^-9 moles.
                 - This cutoff can be adjusted by the user to balance computational speed and accuracy.

             :param float val: The cutoff threshold value in moles (mol).

             **Example:**

             .. code-block:: python

                 # Set the IC amount difference cutoff to 1e-9 moles
                 task.setCutoffICamount(1e-9)
             )doc")

        .def("setCutoffMinimalAmountsIC", &MassTransport::setCutoffMinimalAmountsIC,
             py::arg("val"),
             R"doc(
             Sets the cutoff threshold for minimal amounts of Independent Components (IC) in node bulk compositions.

             .. note::
                 **Typical value:** Approximately 10^-12 moles.

             :param float val: The cutoff threshold value in moles (mol).

             **Example:**

             .. code-block:: python

                 # Set the minimal IC amount cutoff to 1e-12 moles
                 task.setCutoffMinimalAmountsIC(1e-12)
             )doc")

        .def_property_readonly("nPhases", &MassTransport::nPhases,
                               R"doc(
             Gets the number of phases in the parent system (set automatically).

             :return int: The total number of phases in the system.

             **Example:**

             .. code-block:: python

                 phases_count = task.nPhases()
             )doc")

        .def_property_readonly("nElements", &MassTransport::nElements,
                               R"doc(
             Gets the number of independent components (IC) / elements (set automatically).

             :return int: The total number of IC (elements) in the system.

             **Example:**

             .. code-block:: python

                 elements_count = task.nElements()
             )doc")

        .def_property_readonly("nComponents", &MassTransport::nComponents,
                               R"doc(
             Gets the number of dependent components in multicomponent phases (set automatically).

             :return int: The total number of dependent components (species) in the system.

             **Example:**

             .. code-block:: python

                 components_count = task.nComponents()
             )doc")

        .def("setHydraulicParameters", &MassTransport::setHydraulicParameters,
             py::arg("node_ndx"), py::arg("Vt"), py::arg("vp"), py::arg("eps"),
             py::arg("Km"), py::arg("al"), py::arg("Dif"), py::arg("nto"),
             R"doc(
             Changes the initialization of hydraulic parameters for a specific RT node (box).

             This method updates the transport properties of an individual node within the simulation
             grid. Default initialization for all nodes is driven by 'setInitialTotalNodeVolume' and related methods.

             :param int node_ndx: The index of the target node.
             :param float Vt: The total volume of the node (Vt = dx * dy * dz in m³).
             :param float vp: The fluid advection velocity in the node (in m/s).
             :param float eps: The effective (actual) porosity in the node (0 < ε < 1).
             :param float Km: The effective permeability in the node (in m²).
             :param float al: The specific longitudinal dispersivity of the porous medium in the node (in meters).
             :param float Dif: The general diffusivity of dissolved matter in the node (in m²/s).
             :param float nto: The tortuosity factor in the node (typically 1.0 to 10.0).

             **Example:**

             .. code-block:: python

                 # Set hydraulic parameters for node 0
                 task.setHydraulicParameters(0, 1.0, 1e-4, 0.3, 1e-12, 0.1, 2e-9, 1.5)
             )doc")

        .def("summary", &MassTransport::summary, py::arg("width") = 100, py::arg("phases") = false,
             R"doc(
             Text summary of the task, also given by ``print(task)``.

             Shows the mode, the column, the time step that the program will use, the initial systems,
             and the nodes grouped by equal setup (system, role, node type, transport properties, and
             after initNodes: T, P, mass, volume). The last part tells which function changes
             each value.

             :param int width: Maximum line length (default 100). Long tables and texts continue on
                 the next lines. ``print(task)`` uses the default.
             :param bool phases: Also show a table of the phase amounts (mol) in the node groups, one
                 column for each phase present in a node (default False: only the phase names).
             )doc")
        .def("nodeSystem", &MassTransport::nodeSystem, py::arg("node_ndx"),
             "Index of the initial system (dbr file) the node was set up from.")
        .def("nodeType", &MassTransport::nodeType, py::arg("node_ndx"),
             "Node type of the node (see setNodeType).")
        .def("nodesOfSystem", &MassTransport::nodesOfSystem, py::arg("system_ndx"),
             R"doc(
             Indices of all nodes set up from initial system ``system_ndx``.

             In 'A' and 'C' modes nodes 0 and 1 get system 0 (inflowing fluid) and all other nodes system 1.
             )doc")
        .def("nodesOfType", &MassTransport::nodesOfType, py::arg("type"),
             "Indices of all nodes of the given node type (see setNodeType).")
        .def("setNodesType", &MassTransport::setNodesType, py::arg("nodes"), py::arg("type"),
             "As setNodeType(node_ndx, type), for a list of nodes. A node index out of range gives a RuntimeError "
             "before any node is changed.")

        .def("transportProperties", &MassTransport::transportProperties, py::arg("node_ndx"),
             "Transport properties of a node (all fields set). Available after initNodes().")

        .def("setNodes",
             [](MassTransport& self, py::object nodes, py::object composition,
                std::optional<double> T, std::optional<double> P, std::optional<TransportProperties> transport,
                py::object pore_water, py::object solids, std::optional<double> grain_density, std::optional<double> water_mass,
                std::optional<double> water_saturation, std::vector<Part> parts) {
                 self.setNodes(toNodes(nodes), toSetup(composition, T, P, transport, pore_water, solids, grain_density, water_mass, water_saturation, std::move(parts)));
             },
             py::arg("nodes"), py::kw_only(), py::arg("composition") = py::none(),
             py::arg("T") = py::none(), py::arg("P") = py::none(), py::arg("transport") = py::none(),
             py::arg("pore_water") = py::none(), py::arg("solids") = py::none(),
             py::arg("grain_density") = py::none(), py::arg("water_mass") = py::none(),
             py::arg("water_saturation") = py::none(), py::arg("parts") = std::vector<Part>{},
             R"doc(
             Sets up nodes: composition, temperature, pressure and transport properties.

             Call after initNodes() and before runMassTransport(). Only the arguments that are given are
             changed. The transport properties are set first. If ``composition``, ``T`` or ``P`` is given,
             each node is then equilibrated again from scratch.

             :param nodes: A node index or a list of node indices (see nodesOfSystem, nodesOfType).
             :param composition: Element amounts in mol: an xGEMS ``Material`` (its ``b_dict()`` is used),
                 a dict ``{element: amount}`` (elements not listed get 0), or a sequence in element order.
             :param T: Temperature in K.
             :param P: Pressure in Pa.
             :param TransportProperties transport: Physical transport properties; fields not set are not changed.
             :param pore_water: In place of ``composition``: the solution per kg of water (a ``Material`` with
                 1 kg of H2O, a dict or a sequence).
             :param solids: With ``pore_water``: the composition of the dry solid (scaled by its own mass, so
                 any amount can be given).
             :param float grain_density: kg of solid per L of solid grains for ``solids`` (quartz sand about
                 2.65). If not given, from a GEM calculation of the solid. From a measured bulk density:
                 grain density = bulk density / (1 - porosity).
             :param float water_mass: kg of pore water in each node (default 1).
             :param float water_saturation: share of the pore volume with pore water (default 1).
             :param parts: a list of ``Part`` objects: the other phases and solids of the node, e.g.
                 ``Part(air, fraction=0.1, of="pore")``, ``Part(quartz, fraction=1, of="solid")``.

             With ``pore_water`` and ``parts``/``solids`` the program calculates the node composition. For each
             node: pore volume = volume of the pore water / water_saturation, solid volume = pore volume ×
             (1 - porosity) / porosity (transport porosity of the node). Each part fills its share of these
             volumes. The volumes come from GEM calculations of each material alone at the T and P of the node.
             :raises RuntimeError: For a bad node index, a wrong size or an unknown element, or when the
                 GEM calculation fails.

             **Example:**

             .. code-block:: python

                 task.setNodes(range(40, 80), transport=TransportProperties(porosity=0.25))
             )doc")

        .def("setSystem",
             [](MassTransport& self, size_t system, py::object composition,
                std::optional<double> T, std::optional<double> P, std::optional<TransportProperties> transport,
                py::object pore_water, py::object solids, std::optional<double> grain_density, std::optional<double> water_mass,
                std::optional<double> water_saturation, std::vector<Part> parts) {
                 self.setSystem(system, toSetup(composition, T, P, transport, pore_water, solids, grain_density, water_mass, water_saturation, std::move(parts)));
             },
             py::arg("system_ndx"), py::kw_only(), py::arg("composition") = py::none(),
             py::arg("T") = py::none(), py::arg("P") = py::none(), py::arg("transport") = py::none(),
             py::arg("pore_water") = py::none(), py::arg("solids") = py::none(),
             py::arg("grain_density") = py::none(), py::arg("water_mass") = py::none(),
             py::arg("water_saturation") = py::none(), py::arg("parts") = std::vector<Part>{},
             R"doc(
             Sets up all nodes of initial system ``system_ndx``. The same as
             ``setNodes(task.nodesOfSystem(system_ndx), ...)``.

             **Example:**

             .. code-block:: python

                 task.setSystem(0, composition=fluid)     # inflowing fluid, nodes 0 and 1 in mode 'A'
                 task.setSystem(1, composition=column,
                                transport=TransportProperties(porosity=0.32, dispersivity=0.0067))
             )doc")

        .def_property_readonly("nMGPfluxes", &MassTransport::nMGPfluxes,
                               R"doc(
             Gets the number of MGP fluxes defined in the megasystem.

             :return int: The number of MGP fluxes.

             **Example:**

             .. code-block:: python

                 count = task.nMGPfluxes()
             )doc")

        .def("setNumberMGPfluxes", &MassTransport::setNumberMGPfluxes,
             py::arg("nFD"),
             R"doc(
             Sets the number of MGP fluxes defined in the megasystem (nFD >= 0).

             This counter is used in ``S``, ``F``, and ``B`` modes only, where it defines the total number of fluxes
             (including elemental source fluxes) to be defined in the box-flux RT model.

             If a unidirectional chain of boxes with fluxes of aqueous fluid between adjacent boxes is constructed,
             the minimum number of fluxes is **nFD = nNodes()**. This applies to ``S`` and ``F`` modes, where the fluxes
             are defined during default initialization. In ``B`` mode, it is possible to have more than one MGP,
             meaning **nFD > nNodes()**.

             :param int nFD: The total number of MGP fluxes to set.

             **Example:**

             .. code-block:: python

                 # Set the number of fluxes equal to the number of nodes
                 fluxes = task.nNodes()
                 task.setNumberMGPfluxes(fluxes)
             )doc")

        .def("setFluxSourceReceive", &MassTransport::setFluxSourceReceive,
             py::arg("pndx"), py::arg("source"), py::arg("receive"),
             R"doc(
             Changes the initialization of the source and receiver boxes in the flux definition.

             In a typical 1-D box-flux chain arrangement, the first flux connects box indices 0 (source)
             and 1 (receiver); the second flux connects 1 and 2; the third connects 2 and 3, and so on.
             The last cell should contain -1 to indicate that the flux acts as a sink.

             :param int pndx: The index of the target flux/node.
             :param int source: The source (outgoing) box index in the MGP flux definition.
             :param int receive: The receiver (incoming) box index in the MGP flux definition.

             **Example:**

             .. code-block:: python

                 # Set the first flux: flows from box 0 to box 1
                 task.setFluxSourceReceive(0, 0, 1)

                 # Define a sink for the last flux (e.g., at node 10)
                 task.setFluxSourceReceive(10, 10, -1)
             )doc")

        .def("setFluxDefinition", &MassTransport::setFluxDefinition,
             py::arg("pndx"), py::arg("order"), py::arg("rate"),
             py::arg("quantity") = 0.0, py::arg("val") = 0.0,
             R"doc(
             Changes the initialization of the flux definition.

             Configures the flux order, flux rate, and MGP quantity for a specific node (box).

             :param int pndx: The index of the target flux/node.
             :param float order: The MGP flux order {0, 1, 2, 3} (default: 0). Possible values:
                                 - **0**: Zero-order (absolute) flux, with velocity v in mass/time units.
                                 - **1**: First-order to source box mass flux, v in 1/time/kg units.
                                 - **2**: Second-order to source and sink MGP mass flux, v in 1/time/kg² units (proportional to the product of masses in both boxes).
                                 - **3**: First-order flux to MGP mass in receiving box, v in 1/time/kg units.
             :param float rate: The MGP flux rate v, depending on the MGP flux order (defaults: 1 for 'S'; 0.1 for 'B' and 'F').
             :param float quantity: [unused] The MGP quantity (default: 0).
             :param float val: Reserved for future extensions. Currently unused (default: 0).

             **Example:**

             .. code-block:: python

                 task.setFluxDefinition(0, 1.0, 0.5)
             )doc")

        .def("setFluxMGPid", &MassTransport::setFluxMGPid,
             py::arg("pndx"), py::arg("ids"),
             R"doc(
             Changes the initialization of the MGP ID to be moved in this flux.

             In ``S``, ``F``, and ``B`` modes, this parameter defines which MGP is moved
             by the flux (usually dubbed "aq").

             :param int pndx: The index of the target flux/node.
             :param str id: The MGP ID to be moved.

             **Example:**

             .. code-block:: python

                 task.setFluxMGPid(0, "aq")
             )doc")

        .def("setFluxID", &MassTransport::setFluxID,
             py::arg("pndx"), py::arg("id"),
             R"doc(
             Changes the initialization of the short string that identifies the flux.

             This parameter is not used by the internal program logic. It is provided solely
             as an option for the user to identify certain relevant fluxes more easily.

             :param int pndx: The index of the target flux/node.
             :param str id: The string ID used to identify the flux.

             **Example:**

             .. code-block:: python

                 task.setFluxID(0, "custom_flux_name")
             )doc")

        .def_property_readonly("nPhaseGroups", &MassTransport::nPhaseGroups,
             R"doc(
             Gets the number of mobile groups of phases.

             :return int: The number of mobile groups.

             **Example:**

             .. code-block:: python

                 groups = task.nPhaseGroups()
             )doc")

        .def("setNumberPhaseGroups", &MassTransport::setNumberPhaseGroups,
             py::arg("nMGP"),
             R"doc(
             Sets the number of mobile groups of phases (nMGP >= 0).

             This counter is used in ``S``, ``F``, and ``B`` modes only. It defines the number of MGPs
             (mobile phase groups) related to the bulk composition of phases in the outgoing box
             that can be moved by different fluxes.

             The most common case is to set **nMGP = 1**, i.e., one MGP composed of the aqueous
             electrolyte phase only.

             :param int nMGP: The number of mobile groups.

             **Example:**

             .. code-block:: python

                 task.setNumberPhaseGroups(1)
             )doc")

        .def("setPhaseGroupID", &MassTransport::setPhaseGroupID,
             py::arg("gndx"), py::arg("id"),
             R"doc(
             Changes the ID (a short alphanumeric string) that identifies the mobile phase group.

             :param int gndx: The index of the phase group.
             :param str id: The new ID for the phase group.

             **Example:**

             .. code-block:: python

                 task.setPhaseGroupID(0, "MGP_01")
             )doc")

        .def("setUnitsPhaseQuantities", &MassTransport::setUnitsPhaseQuantities,
             py::arg("pndx"), py::arg("units"),
             R"doc(
             Changes the units for setting phase quantities in the MGP.

             Specifies the unit code of measurement for phase quantities. The code is
             a single case-sensitive Latin letter or the '%' character.

             Default code is 'M' (moles).

             **Supported unit codes:**

             * **NUMBER OF MOLES:**
               - ``'Y'``: micromoles
               - ``'h'``: millimoles
               - ``'M'``: moles
             * **ABSOLUTE MASS:**
               - ``'y'``: milligrams
               - ``'g'``: grams
               - ``'G'``: kilograms
             * **VOLUME FRACTION:**
               - ``'v'``: cm³/cm³
               - ``'V'``: volume percent
               - ``'u'``: -lg(cm³/cm³)
             * **MOLE FRACTION:**
               - ``'n'``: mol/mol
               - ``'N'``: mol percent
               - ``'f'``: -lg(mol/mol)
             * **MASS FRACTION:**
               - ``'w'``: kg/kg
               - ``'%'``: mass percent
               - ``'P'``: -lg(kg/kg)
             * **VOLUME CONCENTRATION:**
               - ``'d'``: g/L solution
               - ``'e'``: mg/L solution
               - ``'b'``: µg/L solution
             * **MOLALITY:**
               - ``'m'``: mol/kg H₂O
               - ``'i'``: mmol/kg H₂O
               - ``'p'``: -lg(mol/kg H₂O)
             * **MOLARITY:**
               - ``'L'``: mol/L solution
               - ``'j'``: mmol/L
               - ``'q'``: -lg(mol/L)
             * **MASS CONCENTRATION:**
               - ``'C'``: mol/kg solution
               - ``'c'``: mass % in solution
               - ``'a'``: mg/kg solution

             :param int pndx: The index of the phase.
             :param str units: The unit code character from the supported list above.

             **Example:**

             .. code-block:: python

                 # Example usage: setting units to mass percent
                 task.setUnitsPhaseQuantities(0, '%')
             )doc")

        .def("setPhaseGroupQuantities", &MassTransport::setPhaseGroupQuantities,
             py::arg("gndx"), py::arg("pndx"), py::arg("quantity"),
             R"doc(
             Changes the quantities of phases defining the MGP composition.

             Typically, only the aqueous or gas phase is moved in an MGP; however,
             any combination of phases in the MGP is, in principle, possible.

             :param int gndx: The index of the phase group.
             :param int pndx: The index of the phase.
             :param float quantity: The quantity of the phase.

             **Example:**

             .. code-block:: python

                 task.setPhaseGroupQuantities(0, 2, 1.5)
             )doc")

        .def_property_readonly("nICsourceFluxes", &MassTransport::nICsourceFluxes,
             R"doc(
             Gets the number of elements (Independent Components / IC) source flux compositions.

             :return int: The number of IC source flux compositions.

             **Example:**

             .. code-block:: python

                 count = task.nICsourceFluxes()
             )doc")

        .def("setNumberICsourceFluxes", &MassTransport::setNumberICsourceFluxes,
             py::arg("nSFD"),
             R"doc(
             Sets the number of elemental constant (Independent Components / IC) source flux compositions (nSFD >= 0).

             This counter (used in ``F`` and ``B`` modes) defines the number of different
             elemental constant source compositions. The number of these special fluxes
             should be added to the number of normal MGP fluxes when defining the total
             number of rows (nFD) in the flux definition tables.

             :param int nSFD: The number of IC source flux compositions.

             **Example:**

             .. code-block:: python

                 task.setNumberICsourceFluxes(2)
             )doc")

        .def("setICsourceQuantities", &MassTransport::setICsourceQuantities,
             py::arg("gndx"), py::arg("indx"), py::arg("quantity"),
             R"doc(
             Changes a value in the table of bulk compositions for Independent Component (IC) source fluxes.

             In ``B`` and ``F`` modes, this optional table can be created by the user.
             The respective fluxes need to be specified separately for these nodes.

             :param int gndx: The index of the phase group.
             :param int indx: The index of the Independent Component (IC).
             :param float quantity: The quantity of the Independent Component (IC).

             **Example:**

             .. code-block:: python

                 task.setICsourceQuantities(0, 1, 10.5)
             )doc")

        .def_property_readonly("nParticleTypes", &MassTransport::nParticleTypes,
             R"doc(
             Gets the number of allocated particle types.

             :return int: The number of allocated particle types.

             **Example:**

             .. code-block:: python

                 types = task.nParticleTypes()
             )doc")

        .def("setNumberParticleTypes", &MassTransport::setNumberParticleTypes,
             py::arg("num"),
             R"doc(
             Sets the number of allocated particle types.

             This counter is used by the random walk or particle tracking algorithm.
             Usually, each MGP needs its own particle type.

             :param int num: The number of particle types (must be < 20).

             **Example:**

             .. code-block:: python

                 task.setNumberParticleTypes(5)
             )doc")

        .def("setParticle", &MassTransport::setParticle,
             py::arg("pndx"), py::arg("pmean"), py::arg("pmin"), py::arg("pmax"), py::arg("pparam"),
             R"doc(
             Changes the initialization parameters for particle statistical properties.

             :param int pndx: The index of the particle type.
             :param int pmean: The initial mean number of particles of this type per node.
             :param int pmin: The minimum average total number of particles of this type per node.
             :param int pmax: The maximum average total number of particles of this type per node.
             :param list[int] pparam: A list or tuple of exactly 6 integers defining particle types at t_0 or post-interruption.

             **Example:**

             .. code-block:: python

                 particle_defs = [0, 1, 0, 1, 5, 2]
                 task.setParticle(0, 100, 10, 500, particle_defs)
             )doc")

        .def("setVTKfields", &MassTransport::setVTKfields,
             py::arg("vtk_fields"),
             R"doc(
             Sets the list of selected fields and their indexes for VTK format.

             :param list[tuple[int,int]] vtk_fields: A vector of pairs mapping:
                                                     - ``first``: The object identifier.
                                                     - ``second``: The corresponding field index.

             **Example:**

             .. code-block:: python

                 # Object ID 33 mapped to field index 0, ID 38 mapped to field index 1
                 selected_fields = [(33, 0), (38, 1)]
                 task.setVTKfields(selected_fields)
             )doc")
        .def("setMonitoring", [](MassTransport& self, py::object python_callback) {
            LogFunction_f cxx_callback = [python_callback](ChemicalEngineArray& data, long int step, double time) {
                py::object py_data = py::cast(data, py::return_value_policy::reference);
                python_callback(py_data, step, time);
            };
            self.setMonitoring(cxx_callback);
        },
             R"doc(
             Sets the callback function for monitoring simulations and updating nodes/boxes over time.

             Registers a high-level monitoring function that will be executed at each simulation step.
             The internal decorator automatically forwards the engine context, step index, and physical time.

             Args:
                 python_callback (callable): The monitoring callback function.
                                    Must accept three arguments: (engine: ChemicalEngineArray, step: int, time: float)
             )doc")

        .def("__str__", [](const MassTransport& self) { return self.summary(); },
             "The same as summary().");

}
