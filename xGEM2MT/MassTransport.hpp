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

#include <memory>
#include <string>
#include <xGEM2MT/ChemicalEngineArray.hpp>

namespace xGEM2MT {

/**
 * @file MassTransport.hpp
 * @brief Header file for the Reactive Mass Transport simulator.
 *
 * Mass Transport simulator for simple reactive mass transport models coupled with GEM
 * calculation of equilibrium states in spatially distributed nodes (boxes, volumes) over time steps.
 *
 * The GEM2MT module can simulate reactive transport in three main modes:
 * - **1. Sequential reactor chain and waves** (`'S'`)
 * - **2. Flow-through and box-flux sequence** (`'F'`, `'B'`)
 * - **3. One-dimensional reactive transport** with advection/dispersion/diffusion (`'A'`, `'C'`, `'W'`)
 *
 * @author Dmtrii Kulik, G.D. Miron, S.Dmytriieva
 * @date 2026
 *
 * license GNU General Public License v3 or later
 */
class MassTransport
{

public:
    //
    std::shared_ptr<ChemicalEngineArray> getChemicalEngineArray();

    /**
     * @brief Constructs a MassTransport object.
     *
     * Creates the GEM2MT simulator for simple reactive mass transport models coupled with GEM
     * calculation of equilibrium states in spatially distributed nodes (boxes, volumes) over time steps.
     * The GEM2MT simulator can be configured by choosing its mode of operation and setting up the number
     * of nodes (boxes, volumes).
     *
     * @param mode   Code of GEM2MT mode of operation { S F B A C W }:
     *               - `S`: Sequential reactors chain with waves (moving fluids or solids between reactors)
     *               - `F`: Flow-through reactors chain with fluxes (simultaneous transport then equilibration)
     *               - `B`: Box-flux megasystem (with time-step integration of arbitrary MGP fluxes)
     *               - `A`: 1-D (advection-dispersion) reactive transport using classical finite-difference method
     *               - `C`: 1-D (advection-dispersion) reactive transport using implicit Crank-Nicolson scheme
     *               - `W`: 1-D (advection-diffusion) reactive transport using random-walk method
     * @param nnodes Input number of local equilibrium cells (nodes). Typically, 50 to 250 nodes
     *               provide relatively good spatial resolution and a reasonable computing time.
     *
     * @code
     * // Example: Create a flow-through reactor model with 100 nodes
     * MassTransport task('F', 100);
     * @endcode
     */
    MassTransport(char mode, size_t nnodes);

    /// Deleted copy constructor to prevent shallow copies of C-style pointers.
    MassTransport(const MassTransport& other) = delete;
    /// Deleted copy assignment operator to ensure resource exclusivity.
    MassTransport& operator=(const MassTransport& other) = delete;
    /// Allow moving resources from another instance.
    MassTransport(MassTransport&& other) noexcept = default;
    /// Allow moving resources from another instance.
    MassTransport& operator=(MassTransport&& other) noexcept = default;

    /**
     * @brief Destroys the MassTransport instance and releases allocated resources.
     */
    virtual ~MassTransport();

    /**
     * @brief Initializes the MassTransport object from a GEM-Selektor generated file.
     *
     * @param gem2mt_filename The path to the file containing the definition of the gem2mt task.
     *
     * @throws std::runtime_error If there was a file read error during task initialization.
     */
    void readTask(const std::string& gem2mt_filename);

    /**
     * @brief Initializes the MassTransport object from a JSON string.
     *
     * @param gem2mt_json A valid JSON string containing the definition of the gem2mt task.
     *
     * @throws std::runtime_error If there was an error during task initialization.
     */
    void readTaskString(const std::string& gem2mt_json);

    /**
     * @brief Dumps the MassTransport object state to a JSON string.
     *
     * @return std::string A formatted JSON string representing the current state.
     */
    std::string dump() const;

    /**
     * @brief Dumps the MassTransport object state and task configuration to files.
     *
     * @param gem2mt_filename The path where the main GEM2MT task file will be saved.
     * @param ipm_filename    The path where the IPM work structure file will be saved (optional).
     */
    void dumpTask(const std::string& gem2mt_filename, const std::string& ipm_filename = {}) const;

    /**
     * @brief Initializes all nodes (boxes) from GEM-Selektor IPM files.
     *
     * Reads the specified files and configures the chemical system defined inside them.
     *
     * @param ipm_filename Path to the system definition file (e.g., "my-system-dat.lst").
     * @param dbr_filename Path to the file containing the list of system definition files (e.g., "my-system-dbr.lst").
     *
     * @throws std::runtime_error If there was a file read error during node initialization.
     *
     * @code
     * MassTransport task('F', 100);
     * task.initNodes("my-system-dat.lst", "my-system-dbr.lst");
     * @endcode
     */
    void initNodes(const std::string& ipm_filename, const std::string& dbr_filename);

    /**
     * @brief Initializes all nodes (boxes) using JSON strings.
     *
     * Configures the engine from three JSON inputs:
     * - dch_json: Defines the chemical system.
     * - ipm_json: Contains calculation parameters and settings for the GEMS3K IPM-3 algorithm.
     * - dbr_json: Specifies the input node composition.
     * These strings can be obtained from the GEM-Selektor exported GEMS3K system files as JSON (*-dch.json, *-ipm.json, *-dbr.json).
     *
     * @param dch_json JSON string describing the chemical system.
     * @param ipm_json JSON string with parameter and algorithm settings.
     * @param dbr_json A vector of JSON strings containing node composition details.
     *
     * @throws std::runtime_error If there was an error during node initialization.
     *
     * @code
     * // Example: Initialize using JSON-based input.
     * std::string dch = "{...}";
     * std::string ipm = "{...}";
     * std::vector<std::string> dbr = {"{...}", "{...}"}; // Vector for nodes
     *
     * MassTransport engine('F', 2);
     * engine.initNodesFromStrings(dch, ipm, dbr);
     * @endcode
     */
    void initNodesFromStrings(const std::string& dch_json,
                              const std::string& ipm_json,
                              const std::vector<std::string>& dbr_json);

    /**
     * @brief Executes the coupled finite-difference mass-transport simulation loop.
     *
     * This method runs the main simulation routine, synchronizing the 1-D finite-difference
     * mass transport transport with the GEM thermodynamic equilibrium calculations across all nodes.
     *
     * @note Before calling this method, the simulator must be properly initialized.
     *
     * @throws std::runtime_error If the simulation fails to converge or hits a numerical instability.
     */
    void runMassTransport();

    /**
     * @brief Gets the total number of nodes (boxes).
     *
     * The number of nodes (boxes, reservoirs) is used for spatial discretization of
     * the reactive transport problem (usually a 1-D column or a chain of boxes).
     *
     * @return The total number of nodes.
     *
     * @code
     * size_t total_nodes = task.nNodes();
     * @endcode
     */
    size_t nNodes() const;

    /**
     * @brief Gets the full name of this GEM2MT task.
     *
     * @return The full name of the task.
     *
     * @code
     * std::string task_name = task.name();
     * @endcode
     */
    std::string name() const;
    /**
     * @brief Sets the full name of this GEM2MT task.
     *
     * @param taskname The name to assign to the GEM2MT task.
     *
     * @code
     * task.setName("My_GEM2MT_Simulation");
     * @endcode
     */
    void setName(const std::string& taskname);

    /**
     * @brief Gets the optional enhancement comment for this GEM2MT task.
     *
     * @return The comment string.
     *
     * @code
     * std::string task_comment = task.comment();
     * @endcode
     */
    std::string comment() const;
    /**
     * @brief Sets the optional comment for this GEM2MT task.
     *
     * @param comment The comment text to assign to the task.
     *
     * @code
     * task.setComment("Enable comment.");
     * @endcode
     */
    void setComment(const std::string& comment);

    /**
     * @brief Changes the initialization of the RT nodes (boxes).
     *
     * Each node (box) contains the index of its initial system composition variant.
     * This method updates the index of the initial system variant for the specified node.
     *
     * @note Default distribution behavior:
     * Assigns different fluid and rock composition indices to nodes (assuming fluid=0 and rock=1).
     * For `A` or `C` modes, nodes 0 and 1 contain the initial fluid. For other modes, node 0
     * contains initial 0, and all subsequent nodes contain 1.
     *
     * @param node_ndx   The index of the target node.
     * @param system_ndx The index of the initial system composition variant (from the GEMS3K .dbr file).
     *
     * @code
     * // Example: Set node 5 to use system composition variant 2
     * task.setDistributing(5, 2);
     * @endcode
     */
    void setDistributing(size_t node_ndx, size_t system_ndx);

    /**
     * @brief Changes the initialization type for nodes (boxes).
     *
     * This method defines the boundary condition type of the node:
     * - **0**: Normal node
     * - **1**: Dirichlet source (constant composition source)
     * - **-1**: Dirichlet sink
     * - **2**: Neumann source (constant gradient source)
     * - **-2**: Neumann sink
     * - **3**: Cauchy source (constant flux source)
     * - **-3**: Cauchy sink
     * - **4**: Input time-dependent function
     *
     * @note Default types: 3 for node 0, -3 for the last node, and 0 for all other nodes
     *       (except in 'F' and 'B' modes).
     *
     * @param node_ndx The index of the target node.
     * @param type     The boundary condition type to assign to the node.
     *
     * @code
     * // Set node 0 as a Cauchy source (constant flux)
     * task.setNodeType(0, 3);
     * @endcode
     */
    void setNodeType(size_t node_ndx, long int type);

    /**
     * @brief Sets the initial approximation mode for the GEM IPM solver.
     *
     * Activating a smart initial approximation for GEM calculations at each node and
     * time step can significantly reduce computational time. Deactivate this mode if
     * error messages appear, or to verify results against the automatic initial approximation.
     *
     * The available modes are:
     * - **'+'**: Smart Initial Approximation (SIA)
     * - **'*'**: Internal Smart Initial Approximation (SIA internal)
     * - **'-'**: Automatic Initial Approximation (AIA) — takes longer but is more robust
     *
     * @param flag The character code indicating the desired mode ('+', '*', or '-').
     *
     * @code
     * // Activate Smart Initial Approximation
     * task.setSIA('+');
     * @endcode
     */
    void setSIA(char flag);

    /**
     * @brief Enables or disables continuous debug output for nodes.
     *
     * @note By default, this property is set to `true`.
     *
     * @param enable If `true`, enables continuous debug logging; if `false`, disables it.
     *
     * @code
     * // Disable continuous debug output
     * task.setOutput(false);
     * @endcode
     */
    void setOutput(bool enable);

    /**
     * @brief Enables or disables continuous debug output in VTK format.
     *
     * @note By default, this property is set to `false`.
     *
     * @param enable If `true`, enables continuous VTK output for nodes; if `false`, disables it.
     *
     * @code
     * // Enable continuous output to VTK format
     * task.setOutVTK(true);
     * @endcode
     */
    void setOutVTK(bool enable);

    /**
     * @brief Sets the type of MGP flux phase.
     *
     * Used in default initialization for MGP box-flux transport models.
     *
     * The available flux phase types are:
     * - **0**: Undefined
     * - **1**: Aqueous (aq)
     * - **2**: Gas
     * - **3**: Aqueous + Gas (aq+gas)
     * - **4**: Solids
     *
     * @note By default, this property is set to `1`.
     *
     * @param type The character code indicating the desired MGP flux phase type.
     *
     * @code
     * // Set flux phase type to Gas (2)
     * task.setTypeFluxPhase('2');
     * @endcode
     */
    void setTypeFluxPhase(char type);

    /**
     * @brief Sets the physical time span and step size for RT models.
     *
     * The unit of time measurement is arbitrary (defined by the user), but all
     * transport parameters (such as fluxes and velocities) must be consistent
     * with the chosen time unit (e.g., seconds, hours, or days).
     *
     * The maximum number of time steps is automatically calculated from this time span.
     *
     * @note Default values are:
     * - `start` = 0.0
     * - `stop` = 1000.0
     * - `step` = 1.0
     *
     * @param start The initial simulation time.
     * @param stop  The final simulation time.
     * @param step  The size of each time step.
     *
     * @code
     * // Set simulation from t=0 to t=100 with a step size of 0.5
     * task.setTau(0.0, 100.0, 0.5);
     * @endcode
     */
    void setTau(double start, double stop, double step);

    /**
     * @brief Sets the spatial dimensions of the porous medium.
     *
     * The length `x` is used to calculate the spatial discretization interval:
     * \f\(dx = \frac{x}{nC} \f\) (or node size), and to compute the linear pressure
     * gradient (hydraulic head) as \f\(\frac{P_{last} - P_{first}}{x} \f\).
     *
     * @note The hydraulic head is not used directly to calculate advection velocity.
     * However, to maintain physical realism, a linear pressure gradient (from a few
     * millibars to a few bars) can be set in the initialization script. In this case,
     * each node will operate under a different pressure.
     *
     * @warning Heat transport is not yet implemented in GEM2MT; all nodes are assumed
     * to be at the same temperature.
     *
     * Reasonable values for `x` range between 0.01 and 10.0 meters.
     * By default, this property is set to `x = 1.0`, `y = 0.0`, `z = 0.0`.
     *
     * @param x The spatial dimension along the X-axis (in meters).
     * @param y [unused] Spatial dimension along the Y-axis. Reserved for future 2D extensions.
     * @param z [unused] Spatial dimension along the Z-axis. Reserved for future 3D extensions.
     *
     * @code
     * // Set a 1D column of 5.5 meters
     * task.setSpatialDimensions(5.5);
     * @endcode
     */
    void setSpatialDimensions(double x, double y=0, double z=0);

    /**
     * @brief Sets the time-step reduction factor for advection, diffusion, and box-flux transport.
     *
     * @note By default, this factor is typically set to `1.0`.
     *
     * @param val The time-step reduction factor. Passing `0.0` resets the factor to `1.0`.
     *
     * @code
     * // Reset to default behavior
     * task.setTimeStepReductionFactor(0.0);
     * @endcode
     */
    void setTimeStepReductionFactor(double val);

    /**
     * @brief Sets the initial total node volume (in m³).
     *
     * This volume value is strictly used within the context of porosity-permeability relationships.
     *
     * @note This value will be automatically copied to all nodes during the RT problem's
     *       default initialization phase.
     *
     * @param val The initial total node volume in cubic meters (m³).
     *
     * @code
     * // Set the initial volume of each node to 2.5 cubic meters
     * task.setInitialTotalNodeVolume(2.5);
     * @endcode
     */
    void setInitialTotalNodeVolume(double val);

    /**
     * @brief Sets the fluid advection velocity (in m/s).
     *
     * Typically, the velocity value falls within the range: \f$ 10^{-10} < v_p < 10^{-3} \f$ m/s.
     *
     * @note This value will be automatically copied to all nodes during the RT problem's
     *       default initialization phase.
     *
     * @param val The fluid advection velocity value in meters per second (m/s).
     *
     * @code
     * // Set advection velocity to 1e-5 m/s
     * task.setFluidAdvectionVelocity(1e-5);
     * @endcode
     */
    void setFluidAdvectionVelocity(double val);

    /**
     * @brief Sets the initial effective porosity for nodes.
     *
     * Typically, the porosity value falls within the range: \f\(0 < \varepsilon < 1 \f\).
     *
     * @note
     * - By default, this property is set to `1.0`.
     * - This value will be automatically copied to all nodes during the RT problem's
     *   default initialization phase.
     *
     * @param val The initial effective porosity value (dimensionless fraction).
     *
     * @code
     * // Set the initial effective porosity to 25% (0.25)
     * task.setInitialNodeEffectivePorosity(0.25);
     * @endcode
     */
    void setInitialNodeEffectivePorosity(double val);

    /**
     * @brief Sets the initial effective permeability (in m²).
     *
     * @note
     * - **Unit conversion**: 1 Darcy is approximately equal to \f$ 10^{-12} \f$ m².
     * - This value will be automatically copied to all nodes during the RT problem's
     *   default initialization phase.
     *
     * @param val The initial effective permeability value in square meters (m²).
     *
     * @code
     * // Set the initial effective permeability to 10 Darcy (1e-11 m²)
     * task.setInitialEffectivePermeability(1e-11);
     * @endcode
     */
    void setInitialEffectivePermeability(double val);

    /**
     * @brief Sets the initial value of specific longitudinal dispersivity (in meters).
     *
     * @note
     * - **Typical values**: Range from 0.001 to 1.0 m.
     * - This value will be automatically copied to all nodes during the RT problem's
     *   default initialization phase.
     *
     * @param val The initial specific longitudinal dispersivity value (in meters).
     *
     * @code
     * // Set the initial longitudinal dispersivity to 0.1 meters
     * task.setInitialDispersivity(0.1);
     * @endcode
     */
    void setInitialDispersivity(double val);

    /**
     * @brief Sets the initial general aqueous medium diffusivity (in m²/s).
     *
     * @note
     * - **Typical value**: Approximately \f$ 2 \times 10^{-9} \f$ m²/s.
     * - This value will be automatically copied to all nodes during the RT problem's
     *   default initialization phase.
     *
     * @param val The initial general aqueous medium diffusivity value (in m²/s).
     *
     * @code
     * // Set diffusivity to the typical value for aqueous ions
     * task.setInitialDiffusivity(2e-9);
     * @endcode
     */
    void setInitialDiffusivity(double val);

    /**
     * @brief Sets the initial tortuosity factor for nodes.
     *
     * Typically, the tortuosity factor value falls within the range: \f$ 1 \le \tau_{nto} < 10 \f$.
     *
     * @note
     * - This value will be automatically copied to all nodes during the RT problem's
     *   default initialization phase.
     *
     * @param val The initial tortuosity factor value (dimensionless factor).
     *
     * @code
     * // Set the initial tortuosity factor to a typical value of 1.5
     * task.setInitialTortuosityFactor(1.5);
     * @endcode
     */
    void setInitialTortuosityFactor(double val);

    /**
     * @brief Sets the cutoff threshold for Independent Component (IC) amount differences between time steps.
     *
     * This cutoff is a performance optimization feature. If changes in the node bulk composition
     * (at constant temperature and pressure, and in the absence of kinetics) after a time step
     * are negligibly small, the chemical equilibrium state calculation is skipped.
     *
     * @note
     * - **Typical value**: Approximately \f$ 10^{-9} \f$ moles.
     * - This cutoff can be adjusted by the user to balance computational speed and accuracy.
     *
     * @param val The cutoff threshold value in moles (mol).
     *
     * @code
     * // Set the IC amount difference cutoff to 1e-9 moles
     * task.setCutoffICamount(1e-9);
     * @endcode
     */
    void setCutoffICamount(double val);

    /**
     * @brief Sets the cutoff threshold for minimal amounts of Independent Components (IC) in node bulk compositions.
     *
     * @note **Typical value**: Approximately \f$ 10^{-12} \f$ moles.
     *
     * @param val The cutoff threshold value in moles (mol).
     *
     * @code
     * // Set the minimal IC amount cutoff to 1e-12 moles
     * task.setCutoffMinimalAmountsIC(1e-12);
     * @endcode
     */
    void setCutoffMinimalAmountsIC(double val);

    /**
     * @brief Gets the number of phases in the parent system (set automatically).
     *
     * @return The total number of phases in the system.
     *
     * @code
     * size_t phases_count = task.nPhases();
     * @endcode
     */
    size_t nPhases() const;

    /**
     * @brief Gets the number of independent components / elements (set automatically).
     *
     * @return The total number of independent components (elements) in the system.
     *
     * @code
     * size_t elements_count = task.nElements();
     * @endcode
     */
    size_t nElements() const;

    /**
     * @brief Gets the number of dependent components in multicomponent phases (set automatically).
     *
     * @return The total number of dependent components (species) in the system.
     *
     * @code
     * size_t components_count = task.nComponents();
     * @endcode
     */
    size_t nComponents() const;

    /**
     * @brief Changes the initialization of hydraulic parameters for a specific RT node (box).
     *
     * This method updates the transport properties of an individual node within the simulation
     * grid. Default initialization for all nodes is driven by 'setInitialTotalNodeVolume' and related methods.
     *
     * @param node_ndx   The index of the target node.
     * @param Vt         The total volume of the node (\f$ V_t = dx \times dy \times dz \f$ in m³).
     * @param vp         The fluid advection velocity in the node (in m/s).
     * @param eps        The effective (actual) porosity in the node (\f$ 0 < \varepsilon < 1 \f$).
     * @param Km         The effective permeability in the node (in m²).
     * @param al         The specific longitudinal dispersivity of the porous medium in the node (in meters).
     * @param Dif        The general diffusivity of dissolved matter in the node (in m²/s).
     * @param nto        The tortuosity factor in the node (typically 1.0 to 10.0).
     *
     * @code
     * // Set hydraulic parameters for node 0
     * task.setHydraulicParameters(0, 1.0, 1e-4, 0.3, 1e-12, 0.1, 2e-9, 1.5);
     * @endcode
     */
    void setHydraulicParameters(size_t node_ndx, double Vt, double vp, double eps, double Km, double al, double Dif, double nto);

    /**
     * @brief Gets the number of MGP fluxes defined in the megasystem.
     *
     * @return The number of MGP fluxes.
     *
     * @code
     * size_t count = task.nMGPfluxes();
     * @endcode
     */
    size_t nMGPfluxes() const;

    /**
     * @brief Sets the number of MGP fluxes defined in the megasystem (\f$ nFD \ge 0 \f$).
     *
     * This counter is used in 'S', 'F', and 'B' modes only, where it defines the total number of fluxes
     * (including elemental source fluxes) to be defined in the box-flux RT model.
     *
     * If a unidirectional chain of boxes with fluxes of aqueous fluid between adjacent boxes is constructed,
     * the minimum number of fluxes is \f$ nFD = nNodes() \f$. This applies to 'S' and 'F' modes, where the fluxes
     * are defined during default initialization. In 'B' mode, it is possible to have more than one MGP,
     * meaning \f$ nFD > nNodes() \f$.
     *
     * @param nFD The total number of MGP fluxes to set.
     *
     * @code
     * // Set the number of fluxes equal to the number of nodes
     * long int fluxes = static_cast<long int>(task.nNodes());
     * task.setNumberMGPfluxes(fluxes);
     * @endcode
     */
    void setNumberMGPfluxes(long int nFD);

    /**
     * @brief Changes the initialization of the source and receiver boxes in the flux definition.
     *
     * In a typical 1-D box-flux chain arrangement, the first flux connects box indices 0 (source)
     * and 1 (receiver); the second flux connects 1 and 2; the third connects 2 and 3, and so on.
     * The last cell should contain -1 to indicate that the flux acts as a sink.
     *
     * @param pndx     The index of the target flux/node.
     * @param source   The source (outgoing) box index in the MGP flux definition.
     * @param receive  The receiver (incoming) box index in the MGP flux definition.
     *
     * @code
     * // Set the first flux: flows from box 0 to box 1
     * task.setFluxSourceReceive(0, 0, 1);
     *
     * // Define a sink for the last flux (e.g., at node 10)
     * task.setFluxSourceReceive(10, 10, -1);
     * @endcode
     */
    void setFluxSourceReceive(size_t pndx, long int  source, long int  receive);

    /**
     * @brief Changes the initialization of the flux definition.
     *
     * Configures the flux order, flux rate, and MGP quantity for a specific node (box).
     *
     * @param pndx  The index of the target flux/node.
     * @param order The MGP flux order {0, 1, 2, 3} (default: 0). Possible values:
     *              - 0: Zero-order (absolute) flux, with velocity v in mass/time units.
     *              - 1: First-order to source box mass flux, v in 1/time/kg units.
     *              - 2: Second-order to source and sink MGP mass flux, v in 1/time/kg² units
     *                   (proportional to the product of masses in both boxes).
     *              - 3: First-order flux to MGP mass in receiving box, v in 1/time/kg units.
     * @param rate The MGP flux rate v, depending on the MGP flux order
     *             (defaults: 1 for 'S'; 0.1 for 'B' and 'F').
     * @param quantity [unused] The MGP quantity (default: 0).
     * @param val Reserved for future extensions. Currently unused (default: 0).
     *
     * @code
     * task.setFluxDefinition(0, 1, 0.5);
     * @endcode
     */
    void setFluxDefinition(size_t pndx, double order, double rate, double quantity=0, double val=0);

    /**
     * @brief Changes the initialization of the MGP ID to be moved in this flux.
     *
     * In 'S', 'F', and 'B' modes, this parameter defines which MGP is moved
     * by the flux (usually dubbed "aq").
     *
     * @param pndx  The index of the target flux/node.
     * @param id    The MGP ID to be moved.
     *
     * @code
     * task.setFluxMGPid(0, "aq");
     * @endcode
     */
    void setFluxMGPid(size_t pndx, const std::string& id);

    /**
     * @brief Changes the initialization of the short string that identifies the flux.
     *
     * This parameter is not used by the internal program logic. It is provided solely
     * as an option for the user to identify certain relevant fluxes more easily.
     *
     * @param pndx  The index of the target flux/node.
     * @param id    The string ID used to identify the flux.
     *
     * @code
     * task.setFluxID(0, "custom_flux_name");
     * @endcode
     */
    void setFluxID(size_t pndx, const std::string& id);

    /**
    * @brief Gets the number of mobile groups of phases.
    *
    * @return The number of mobile groups.
    *
    * @code
    * // Example usage:
    * size_t groups = task.nPhaseGroups();
    * @endcode
    */
    size_t nPhaseGroups() const;

    /**
     * @brief Sets the number of mobile groups of phases (nMGP >= 0).
     *
     * This counter is used in 'S', 'F', and 'B' modes only. It defines the number of MGPs
     * (mobile phase groups) related to the bulk composition of phases in the outgoing box
     * that can be moved by different fluxes.
     *
     * The most common case is to set nMGP = 1, i.e., one MGP composed of the aqueous
     * electrolyte phase only.
     *
     * @param nMGP The number of mobile groups.
     *
     * @code
     * // Example usage:
     * task.setNumberPhaseGroups(1);
     * @endcode
     */
    void setNumberPhaseGroups(long int nMGP);

    /**
     * @brief Changes the ID (a short alphanumeric string) that identifies the mobile phase group.
     *
     * @param gndx The index of the phase group.
     * @param id   The new ID for the phase group.
     *
     * @code
     * // Example usage:
     * task.setPhaseGroupID(0, "MGP_01");
     * @endcode
     */
    void setPhaseGroupID(size_t gndx, const std::string& id);

    /**
     * @brief Changes the units for setting phase quantities in the MGP.
     *
     * Specifies the unit code of measurement for phase quantities. The code is
     * a single case-sensitive Latin letter or the '%' character.
     *
     * Default code is 'M' (moles).
     *
     * **Supported unit codes:**
     * - **NUMBER OF MOLES:**
     *   - 'Y': micromoles
     *   - 'h': millimoles
     *   - 'M': moles
     * - **ABSOLUTE MASS:**
     *   - 'y': milligrams
     *   - 'g': grams
     *   - 'G': kilograms
     * - **VOLUME FRACTION:**
     *   - 'v': cm³/cm³
     *   - 'V': volume percent
     *   - 'u': -lg(cm³/cm³)
     * - **MOLE FRACTION:**
     *   - 'n': mol/mol
     *   - 'N': mol percent
     *   - 'f': -lg(mol/mol)
     * - **MASS FRACTION:**
     *   - 'w': kg/kg
     *   - '%': mass percent
     *   - 'P': -lg(kg/kg)
     * - **VOLUME CONCENTRATION:**
     *   - 'd': g/L solution
     *   - 'e': mg/L solution
     *   - 'b': µg/L solution
     * - **MOLALITY:**
     *   - 'm': mol/kg H₂O
     *   - 'i': mmol/kg H₂O
     *   - 'p': -lg(mol/kg H₂O)
     * - **MOLARITY:**
     *   - 'L': mol/L solution
     *   - 'j': mmol/L
     *   - 'q': -lg(mol/L)
     * - **MASS CONCENTRATION:**
     *   - 'C': mol/kg solution
     *   - 'c': mass % in solution
     *   - 'a': mg/kg solution
     *
     * @param pndx  The index of the phase.
     * @param units The unit code character from the supported list above.
     *
     * @code
     * // Example usage: setting units to mass percent
     * task.setUnitsPhaseQuantities(0, '%');
     * @endcode
     */
    void setUnitsPhaseQuantities(size_t pndx, char units);

    /**
     * @brief Changes the quantities of phases defining the MGP composition.
     *
     * Typically, only the aqueous or gas phase is moved in an MGP; however,
     * any combination of phases in the MGP is, in principle, possible.
     *
     * @param gndx     The index of the phase group.
     * @param pndx     The index of the phase.
     * @param quantity The quantity of the phase.
     *
     * @code
     * // Example usage:
     * task.setPhaseGroupQuantities(0, 2, 1.5);
     * @endcode
     */
    void setPhaseGroupQuantities(size_t gndx, size_t pndx, double quantity);

    /**
    * @brief Gets the number of elements (Independent Components / IC) source flux compositions.
    *
    * @return The number of IC source flux compositions.
    *
    * @code
    * // Example usage:
    * size_t count = task.nICsourceFluxes();
    * @endcode
    */
    size_t nICsourceFluxes() const;

    /**
     * @brief Sets the number of elemental constant (Independent Components / IC) source flux compositions (nSFD >= 0).
     *
     * This counter (used in 'F' and 'B' modes) defines the number of different
     * elemental constant source compositions. The number of these special fluxes
     * should be added to the number of normal MGP fluxes when defining the total
     * number of rows (nFD) in the flux definition tables.
     *
     * @param nSFD The number of IC source flux compositions.
     *
     * @code
     * // Example usage:
     * task.setNumberICsourceFluxes(2);
     * @endcode
     */
    void setNumberICsourceFluxes(long int nSFD);

    /**
     * @brief Changes a value in the table of bulk compositions for Independent Component (IC) source fluxes.
     *
     * In 'B' and 'F' modes, this optional table can be created by the user.
     * The respective fluxes need to be specified separately for these nodes.
     *
     * @param gndx     The index of the phase group.
     * @param indx     The index of the Independent Component (IC).
     * @param quantity The quantity of the Independent Component (IC).
     *
     * @code
     * // Example usage:
     * task.setICsourceQuantities(0, 1, 10.5);
     * @endcode
     */
    void setICsourceQuantities(size_t gndx, size_t indx, double quantity);

    /**
    * @brief Gets the number of allocated particle types.
    *
    * @return The number of allocated particle types.
    *
    * @code
    * // Example usage:
    * size_t types = task.nParticleTypes();
    * @endcode
    */
    size_t nParticleTypes() const;

    /**
      * @brief Sets the number of allocated particle types.
      *
      * This counter is used by the random walk or particle tracking algorithm.
      * Usually, each MGP needs its own particle type.
      *
      * @param num The number of particle types (must be < 20).
      *
      * @code
      * // Example usage:
      * task.setNumberParticleTypes(5);
      * @endcode
      */
    void setNumberParticleTypes(long int num);

    /**
     * @brief Changes the initialization parameters for particle statistical properties.
     *
     * @param pndx   The index of the particle type.
     * @param pmean  The initial mean number of particles of this type per node.
     * @param pmin   The minimum average total number of particles of this type per node.
     * @param pmax   The maximum average total number of particles of this type per node.
     * @param pparam A `std::array` of 6 elements defining particle types at t_0 or post-interruption.
     *
     * @code
     * std::array<long int, 6> particle_defs = {0, 1, 0, 1, 5, 2};
     * task.setParticle(0, 100, 10, 500, particle_defs);
     * @endcode
     */
    void setParticle(size_t pndx, long int pmean, long int pmin, long int pmax, const std::array<long int, 6>& pparam);

    /**
     * @brief Sets the list of selected fields and their indexes for VTK format.
     *
     * @param vtk_fields A vector of pairs mapping:
     *                   - `first`: The object identifier.
     *                   - `second`: The corresponding field index.
     *
     * @code
     * const std::vector<std::pair<int, int>> selected_fields = {
     *     {33, 0}, // Object ID 33 mapped to field index 0
     *     {38, 1}  // Object ID 38 mapped to field index 1
     * };
     * task.setVTKfields(selected_fields);
     * @endcode
     */
    void setVTKfields(const std::vector<std::pair<int, int>>& vtk_fields);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

    friend class ChemicalEngineArray;
};

/**
 * @brief Outputs the internal state of the MassTransport object to a stream.
 *
 * This operator formats the simulator's parameters, active mode, and current
 * node states into a human-readable format, suitable for logging or console output.
 *
 * @param out    The target output stream (e.g., `std::cout` or a file stream).
 * @param engine The MassTransport instance to print.
 * @return std::ostream& A reference to the destination output stream to allow chaining.
 *
 * @relates MassTransport
 *
 * @code
 * MassTransport task('F', 100);
 * std::cout << task << std::endl;
 * @endcode
 */
std::ostream& operator<<(std::ostream& out, const MassTransport& engine);

} // namespace xGEM2MT
