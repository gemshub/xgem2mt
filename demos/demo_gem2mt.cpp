// xGEMS2MT is a C++ and Python library for mass transport modeling with Gibbs energy minimization
//
// Copyright (C) 2026 S.Dmytriieva, Dmitrii Kulik
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License version 2.1 
// as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// xGEMS includes
#include <xGEM2MT/MassTransport.hpp>

struct DemoTask {
    /// Code of GEM2MT mode of operation { S F B A C W }
    char mode_operation = 'S';

    /// Input number of local equilibrium cells (nodes).
    size_t nnodes = 51;

    /// The path to the file containing the definition of the gem2mt task (optional).
    std::string gem2mt_filename;

    /// Path to the system definition file (e.g., "my-system-dat.lst").
    std::string ipmfiles_lst_name = "TestAD/CalcColumn-dat.lst";

    /// Path to the file containing the list of system definition files (e.g., "my-system-dbr.lst").
    std::string dbr_lst_name = "TestAD/CalcColumn-dbr.lst";
};

void show_usage(const std::string &name);
int extract_args(int argc, char* argv[], DemoTask& task_data);
static int task_from_file(const std::string& gem2mt_file, const std::string& ipm_lst, const std::string& dbr_lst);
static int task_A(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst);
static int task_C(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst);
static int task_W(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst);
static int task_F(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst);
static int task_S(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst);
static int task_B(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst);

int main(int argc, char **argv)
{
    xGEM2MT::update_loggers(true, "xGEM2MT.log", 2);

    DemoTask task_data;

    try  {

        if(extract_args( argc, argv,task_data )) {
            return 1;
        }

        switch (task_data.mode_operation) {
        case 'S':
            return  task_S(task_data.nnodes, task_data.ipmfiles_lst_name, task_data.dbr_lst_name);
        case 'F':
            return  task_F(task_data.nnodes, task_data.ipmfiles_lst_name, task_data.dbr_lst_name);
        case 'B':
            return  task_B(task_data.nnodes, task_data.ipmfiles_lst_name, task_data.dbr_lst_name);
        case 'A':
            return  task_A(task_data.nnodes, task_data.ipmfiles_lst_name, task_data.dbr_lst_name);
        case 'C':
            return  task_C(task_data.nnodes, task_data.ipmfiles_lst_name, task_data.dbr_lst_name);
        case 'W':
            return  task_W(task_data.nnodes, task_data.ipmfiles_lst_name, task_data.dbr_lst_name);
        case 'I':
            return  task_from_file(task_data.gem2mt_filename, task_data.ipmfiles_lst_name, task_data.dbr_lst_name);
        default:
            break;
        }
        return 0;
    }
    catch(std::runtime_error& xcpt)  {
        std::cout << xcpt.what() << std::endl;
    }
    return 1;
}

/**
'S'     Sequential reactors chain with waves (moving fluids or solids between reactors one-by-one)

This flexible scheme consists of many "reactors" (boxes) arranged into a 1-D chain, where at each step,
the whole fluid phase is moved (copied) from the previous into the next box. The constant composition "fresh"
fluid enters the very first box with a prescribed fluid-rock mass ratio and reacts there. The equilibrated fluid
is taken to the next box (at prescribed fluid-rock ratio, replacing the fluid left from the previous "wave")
and reacts there, and so on until the last box is reacted, thus completing the so-called "wave".
Any desired number of "waves" can be performed. At the beginning of second wave, the "fresh" fluid replaces
the old fluid phase in the first box. After equilibration, the fluid is taken into the second box to react, and so on.
From "wave" to "wave", the propagation of rock alteration zones can be followed by plotting amounts of phases in all boxes.
A variant of 'S' mode allows to move all solids instead of the fluid phase from box to box, like in a dialysis bag
(with an idea to follow the fluid evolution).

 See GUI example: CalDolCol2:G:CalcColumn:0:0:1:25:0:Test3:S
*/
int task_S(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst)
{
    // (1) Create the  MassTransport simulator
    xGEM2MT::MassTransport mt_task('S', nnodes);
    mt_task.setName("Test S mode (initial 2.5 mmol/L calcite)  1 to 96 C");

    // (2) Set up sizes, flags and values different from default
    // Set physical time iterator (start,end,step)
    mt_task.setTau(0, 30, 0.3);
    // Set number of mobile groups of phases
    mt_task.setNumberPhaseGroups(1);
    // Set number of MGP fluxes defined in the megasystem
    mt_task.setNumberMGPfluxes(nnodes);
    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task.setSIA('-');
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
    mt_task.setTypeFluxPhase('1');
    // Use non stop debug output for nodes
    mt_task.setOutput(true);

    // // ?? only to compare Set initial node effective porosity (0 < eps < 1), usually 1
    // mt_task.setInitialNodeEffectivePorosity(1e-9);
    // // ?? only to compare Set initial effective permeability, m2, usually 1
    // mt_task.setInitialEffectivePermeability(1e-12);
    // // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    // mt_task.setCutoffICamount(0);
    // // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    // mt_task.setCutoffMinimalAmountsIC(0);

    // (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    mt_task.initNodes(ipm_lst, dbr_lst);

    // (4) Change/define some other gem2mt arrays
    // Set units for aq phase quantities
    mt_task.setUnitsPhaseQuantities(0, 'n');
    for(size_t ii = 0; ii < mt_task.nMGPfluxes(); ++ii) {
        // Set the flux definition: flux order, flux rate, MGP quantity
        mt_task.setFluxDefinition(ii, 0., 1., 0, 0);
    }

    mt_task.dumpTask("gem2mt_out.dat");

    // (5) Here we call the mass-transport finite-difference coupled routine
    mt_task.runMassTransport();

    return 0;
}

/**
'F'     Flow-through reactors chain with fluxes (simultaneous transport then equilibration)

This scheme also consists of many "reactors" (boxes) arranged into a 1-D chain connected by 1st order fluxes.
At each (time) step, part of the fluid phase is moved from the previous into the next box simultaneously for all boxes.
The very first box is usually defined as a constant source of "fresh" fluid, and the last box in the chain - as having
the infinite sink for the fluid. 1st order flux means that the mass (of fluid) transported per (time) step is proportional
to the mass of fluid phase in the outgoing box. So, the flux rate constant (typically about 0.1 to 0.3) actually tells what
fraction of the fluid mass is transported per (time) step. The effect of this feedback is that, even if the boxes initially
contained different masses of fluid, with time the masses of fluid would become equal. This scheme simulates the propagation
of rock alteration zones, which can be followed by plotting amounts of phases in all boxes. The ODE integration (as in B mode)
is not performed. Compared with S mode, the F mode produces results more similar to those of classical advection RT algorithms.

 See GUI example: CalDolCol2:G:CalcColumn:0:0:1:25:0:Test2:F:
*/
int task_F(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst)
{
    // (1) Create the  MassTransport simulator
    xGEM2MT::MassTransport mt_task('F', nnodes);
    mt_task.setName("Test F mode (flow-through reactors)");

    // (2) Set up sizes, flags and values different from default
    // Set physical time iterator (start,end,step)
    mt_task.setTau(0, 1200, 1);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task.setSpatialDimensions(0, 0, 0);
    // Set number of mobile groups of phases
    mt_task.setNumberPhaseGroups(1);
    // Set number of MGP fluxes defined in the megasystem
    mt_task.setNumberMGPfluxes(nnodes);
    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task.setSIA('-');
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
    mt_task.setTypeFluxPhase('1');
    // Use non stop debug output for nodes
    mt_task.setOutput(true);

    // // ?? only to compare Set initial node effective porosity (0 < eps < 1), usually 1
    // mt_task.setInitialNodeEffectivePorosity(1e-9);
    // // ?? only to compare Set initial effective permeability, m2, usually 1
    // mt_task.setInitialEffectivePermeability(1e-12);
    // // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    // mt_task.setCutoffICamount(0);
    // // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    // mt_task.setCutoffMinimalAmountsIC(0);

    // (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    mt_task.initNodes(ipm_lst, dbr_lst);

    // (4) Change/define some other gem2mt arrays
    for(size_t ii=0; ii<mt_task.nPhaseGroups(); ++ii) {
        // Set ID of mobile phase group
        mt_task.setPhaseGroupID(ii, "phg"+std::to_string(ii+1));
        for(size_t k=0; k<mt_task.nPhases(); ++k) {
            // Set quantities of phases in MGP
            mt_task.setPhaseGroupQuantities(ii, k, (k==0 ? 1.: 0.));
        }
    }
    for(size_t k=0; k<mt_task.nPhases(); ++k) {
        // Set units for setting phase quantities in MGP (see PGT )
        mt_task.setUnitsPhaseQuantities(k, 'M');
    }
    for(size_t ii=0; ii<mt_task.nMGPfluxes(); ++ii) {
        // Set Source/Receive box index in the flux definition
        mt_task.setFluxSourceReceive(ii, ii, (ii<mt_task.nMGPfluxes()-1 ? ii+1: -1));
        // Set the flux defnition: flux order, flux rate, MGP quantity
        mt_task.setFluxDefinition(ii, 1., 0.5, 0, 0);
        //  Set the ID of MGP to move in this flux
        mt_task.setFluxMGPid(ii, "phg1");
        //  Set IDs of fluxes
        mt_task.setFluxID(ii, "qj");
    }

    // Set list of selected fields and indexes to VTK format
    mt_task.setVTKfields({{41,2},{41,6},{41,25},{41,26},{41,17},{41,13}});
    mt_task.setOutVTK(false);  // allocated but no write

    mt_task.dumpTask("gem2mt_out.dat");

    // (5) Here we call the mass-transport finite-difference coupled routine
    mt_task.runMassTransport();

    return 0;
}

/**
'B'     Box-flux megasystem (with time-step integration of complex arbitrary MGP fluxes)

This scheme is similar to the F scheme, but allows the most general setup of the box-flux dynamic transport
scheme, not only a unidirectionally connected chain of boxes. In the B scheme, shortcut fluxes and loops are allowed,
as well as zero-, first- and second-order fluxes, box (point) sinks and box (point) sources of given elemental composition.
The transport step is performed by the built-in integrator of a large  ODE system of time derivatives of masses of
each chemical element in each box (reservoir). After the transport step, chemical equilibria are re-calculated in all boxes.
Normally, with time (number of steps), such complex dynamic systems evolve to a steady state when masses and compositions of
reservoirs do not change anymore. This scheme may be good for simulating water-sediment columns or geochemical cycles.

 See GUI example: CalDolCol2:G:CalcColumn:0:0:1:25:0:Test4:B:
*/
int task_B(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst)
{
    // (1) Create the  MassTransport simulator
    xGEM2MT::MassTransport mt_task('B', nnodes);
    mt_task.setName("Test B mode (flow-through reactors), T = 1 to 99 C");

    // (2) Set up sizes, flags and values different from default
    // Set physical time iterator (start,end,step)
    mt_task.setTau(0, 10000, 1);
    // Set number of mobile groups of phases
    mt_task.setNumberPhaseGroups(1);
    // Set number of MGP fluxes defined in the megasystem
    mt_task.setNumberMGPfluxes(nnodes);
    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task.setSIA('-');
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
    mt_task.setTypeFluxPhase('1');
    // Use non stop debug output for nodes
    mt_task.setOutput(true);

    // Change from default: Other nodes normal, the last two set as constant source and sink
    mt_task.setNodeType(nnodes-2, 0);

    // // ?? only to compare Set initial node effective porosity (0 < eps < 1), usually 1
    // mt_task.setInitialNodeEffectivePorosity(1e-9);
    // // ?? only to compare Set initial effective permeability, m2, usually 1
    // mt_task.setInitialEffectivePermeability(1e-12);
    // // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    // mt_task.setCutoffICamount(0);
    // // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    // mt_task.setCutoffMinimalAmountsIC(0);

    // (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    mt_task.initNodes(ipm_lst, dbr_lst);

    // (4) Change/define some other gem2mt arrays
    // Set ID of the first mobile phase group
    mt_task.setPhaseGroupID(0, "Pg1");
    for(size_t ii=0; ii<mt_task.nMGPfluxes(); ++ii) {
        //  Set the ID of MGP to move in this flux
        mt_task.setFluxMGPid(ii, "Pg1");
    }
    // Change from default: flux order 1 (proportional to source MPG mass)
    mt_task.setFluxDefinition(0, 1., 0.1, 0, 0);

    mt_task.dumpTask("gem2mt_out.dat");

    // (5) Here we call the mass-transport finite-difference coupled routine
    mt_task.runMassTransport();

    return 0;
}

/**
'A'     1-D (advection-dispersion) reactive transport using classical finite-difference method

In this mode, the GEM2MT module performs a simple "textbook" finite-difference simulation of 1-D column reactive transport,
in which not "tracer" but the whole composition of aqueous phase including water is moved between nodes during the "transport loop",
followed by a "chemical equilibration loop" for all nodes at each time step (the so-called operator-splitting sequential
non-iterative coupling of reaction and transport).

The column of a given length l is thought to be divided into n = 100 to 1000 nodes of equal size (l/n) with implicit diameter
that defines the same characteristic volume (e.g. 1 cm3) of the node. Each node is characterized by several transport properties
(porosity, permeability, dispersivity, fluid advection velocity, temperature, pressure).

The "fresh" fluid is fed from the left into node 0 and 1 that are declared as infinite flux sources (type 3) and only contain fluid.
Other nodes are filled with "rock+fluid" of a given initial composition and porosity that may or may not correspond to real
fluid/rock volume ratio. The very last rock-fluid node is marked as an infinite constant-flux sink of aqueous fluid (type -3).

The time step duration is determined from the Courant number, and can be augmented by dividing by a prescribed time step factor (>1).
The advection velocity is input instead of the hydraulic head to simplify handling of chemical systems (that can be formally
processed at constant pressure, even though this is not quite physically correct).

 See GUI example: CalDolCol2:G:CalcColumn:0:0:1:25:0:Test4:B:
*/
int task_A(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst)
{
    // (1) Create the  MassTransport simulator
    xGEM2MT::MassTransport mt_task('A', nnodes);
    mt_task.setName("Test of 1D coupled advection problem (dissolved Ca, Mg)");

    // (2) Set up sizes, flags and values different from default
    // Set physical time iterator (start,end,step)
    mt_task.setTau(0, 500000, 5000);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task.setSpatialDimensions(0.2, 0, 0);
    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task.setSIA('-');
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
    mt_task.setTypeFluxPhase('1');
    // Use non stop debug output for nodes
    mt_task.setOutput(true);

    // Set advection/diffusion mass transport: time step reduction factor (usually 1)
    mt_task.setTimeStepReductionFactor(1.);
    // Set initial total node volume (m^3)
    mt_task.setInitialTotalNodeVolume(0.001);
    // Set fluid advection velocity (m/sec)
    mt_task.setFluidAdvectionVelocity(1e-7);
    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task.setInitialNodeEffectivePorosity(1.);
    // Set initial effective permeability, m2, usually 1
    mt_task.setInitialEffectivePermeability(1.);
    // Set initial value of specific longitudinal dispersivity (m), usually 1e-3
    mt_task.setInitialDispersivity(0.001);
    // Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
    mt_task.setInitialDiffusivity(0.);
    // Set initial tortuosity factor, usually 1
    mt_task.setInitialTortuosityFactor(1.);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task.setCutoffICamount(1e-9);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task.setCutoffMinimalAmountsIC(1e-11);

    // (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    mt_task.initNodes(ipm_lst, dbr_lst);
    mt_task.dumpTask("gem2mt_out.dat");

    // (4) Here we call the mass-transport finite-difference coupled routine
    mt_task.runMassTransport();

    return 0;
}

/**
'C'     1-D (advection-dispersion) reactive transport using implicit Crank-Nicolson scheme

In this mode, the problem setup is the same as in 'A' mode, with the only difference that the simulation
is performed using the Crank-Nicholson scheme, believed to be more robust regarding to numerical dispersion and/or
error accumulation. Both 'A' and 'C' methods are only applicable to advection-dominated cases.

 See GUI example: CalDolCol2:G:CalcColumn:0:0:1:25:0:1D-Cc-MgCl2-tmod:C:
*/
int task_C(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst)
{
    // (1) Create the  MassTransport simulator
    xGEM2MT::MassTransport mt_task('C', nnodes);
    mt_task.setName("Test of 1D coupled advection problem (dissolved Ca, Mg)");

    // (2) Set up sizes, flags and values different from default
    // Set physical time iterator (start,end,step)
    mt_task.setTau(0, 30000, 3);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task.setSpatialDimensions(0.2, 0, 0);
    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task.setSIA('-');
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
    mt_task.setTypeFluxPhase('1');
    // Use non stop debug output for nodes
    mt_task.setOutput(true);

    // Set advection/diffusion mass transport: time step reduction factor (usually 1)
    mt_task.setTimeStepReductionFactor(2.);
    // Set initial total node volume (m^3)
    mt_task.setInitialTotalNodeVolume(0.001);
    // Set fluid advection velocity (m/sec)
    mt_task.setFluidAdvectionVelocity(2e-6);
    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task.setInitialNodeEffectivePorosity(0.5);
    // Set initial effective permeability, m2, usually 1
    mt_task.setInitialEffectivePermeability(1e-12);
    // Set initial value of specific longitudinal dispersivity (m), usually 1e-3
    mt_task.setInitialDispersivity(0.001);
    // Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
    mt_task.setInitialDiffusivity(2e-09);
    // Set initial tortuosity factor, usually 1
    mt_task.setInitialTortuosityFactor(1.);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task.setCutoffICamount(1e-9);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task.setCutoffMinimalAmountsIC(1e-11);

    // (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    mt_task.initNodes(ipm_lst, dbr_lst);

    // (4) Change/define some other gem2mt arrays
    // Set list of selected fields and indexes to VTK format
    mt_task.setVTKfields({{33,0},{32,0},{38,1},{38,4}});

    mt_task.dumpTask("gem2mt_out.dat");

    // (5) Here we call the mass-transport finite-difference coupled routine
    mt_task.runMassTransport();

    return 0;
}

/**
'W'    1-D (advection-diffusion) reactive transport using the RW (random-walk) algorithm

This method is suitable for both advection- and diffusion-dominated cases of reactive transport.
Geometrical setup of 1-D column and its discretization into nodes is, in principle, the same as in the 'A' or 'C' modes
(with the difference that only one fluid source node is needed instead of two). Another difference is that solute diffusivity
(diffusion coefficient) Dif must be provided at input for each node. In the simplest case, the same Dif value is applied to
all aqueous solutes and water-solvent; as solutes, total dissolved elements are treated as default.

The random-walk method effectively simulates the Brownian movement of "particles" of different types (about 500 to 1500 particles per node).
At each time step, after the chemical equilibration loop, to one particle type, the elemental composition is assigned that corresponds
to the current composition of the aqueous fluid in the node in which this particle currently resides.
Other particle types can be assigned for different cations and anions, if their diffusion is to be considered individually,
or to solid phase compositions that are supposed to migrate as nanoparticles.

At the following transport loop, the velocity and the direction of each particle is determined by (pseudo)random numbers obtained
from a random number generator subroutine (Brownian movement) and augmented by a flow component in one direction
(if advection also takes place) or even by the electrical potentials on solids in nodes. The path (velocity times time step)
of each particle may move it within a given node, or from this to the next or to the previous node. The bookkeeping of particles
that have crossed the boundaries of each node from all the adjacent nodes results in re-calculation of augmented bulk composition of nodes,
and thus makes them ready for the next chemical equilibration loop. The time stepping should be selected such that no particle can cross
more than one node boundary in one time step.

 See GUI example: CalDolCol2:G:CalcColumn:0:0:1:25:0:1D-DifMgCl2-center:W:
*/
int task_W(size_t nnodes, const std::string& ipm_lst, const std::string& dbr_lst)
{
    // (1) Create the  MassTransport simulator
    xGEM2MT::MassTransport mt_task('W', nnodes);
    mt_task.setName("Test of 1D coupled advection problem (dissolved Ca, Mg)");

    // (2) Set up sizes, flags and values different from default
    // Set number of allocated particle types < 20
    mt_task.setNumberParticleTypes(1);
    // Set physical time iterator (start,end,step)
    mt_task.setTau(0, 1000000, 1000);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task.setSpatialDimensions(0.2, 0, 0);
    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task.setSIA('-');
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
    mt_task.setTypeFluxPhase('1');
    // Use non stop debug output for nodes
    mt_task.setOutput(true);

    // Set advection/diffusion mass transport: time step reduction factor (usually 1)
    mt_task.setTimeStepReductionFactor(5.);
    // Set initial total node volume (m^3)
    mt_task.setInitialTotalNodeVolume(0.001);
    // Set fluid advection velocity (m/sec)
    mt_task.setFluidAdvectionVelocity(2e-6);
    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task.setInitialNodeEffectivePorosity(0.5);
    // Set initial effective permeability, m2, usually 1
    mt_task.setInitialEffectivePermeability(1e-12);
    // Set initial value of specific longitudinal dispersivity (m), usually 1e-3
    mt_task.setInitialDispersivity(0.001);
    // Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
    mt_task.setInitialDiffusivity(2e-09);
    // Set initial tortuosity factor, usually 1
    mt_task.setInitialTortuosityFactor(1.);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task.setCutoffICamount(1e-9);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task.setCutoffMinimalAmountsIC(1e-11);

    // If need update input node distributing
    // Change index of initial system variant for the center node
    mt_task.setDistributing(nnodes/2, 0);
    // Change the type for for the center node
    mt_task.setNodeType(nnodes/2, 3);

    // (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
    // allocate and set default values for gem2mt arrays
    mt_task.initNodes(ipm_lst, dbr_lst);

    // (4) Change/define some other gem2mt arrays
    // Set list of selected fields and indexes to VTK format
    mt_task.setVTKfields({{33,0},{32,0},{38,1},{38,4}});
    // Set of particle statistic property
    mt_task.setParticle(0, 1000, 500, 1500, {0, 11, 21, 0, 0, 0});
    //If need update initial hydraulic parameters in node: Vt, vp, eps, Km, al, Dif,  nto
    //mt_task.setHydraulicParameters(1, 0.0015, 2.1e-06, 0.6, 1e-11, 0.0015, 3e-09, 1);

    mt_task.dumpTask("gem2mt_out.dat");

    // (5) Here we call the mass-transport finite-difference coupled routine
    mt_task.runMassTransport();

    return 0;
}

int task_from_file(const std::string& gem2mt_file, const std::string& ipm_lst, const std::string& dbr_lst)
{
    // Create the  MassTransport simulator
    xGEM2MT::MassTransport mt_task('S', 101);

    // Here we read the GEM2MT structure, prepared from GEMS or by hand
    mt_task.readTask(gem2mt_file);

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    mt_task.initNodes(ipm_lst, dbr_lst);

    //mt_task.dumpTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    mt_task.runMassTransport();

    return 0;
}

void show_usage(const std::string &name)
{
    std::cout << "Usage: " << name << " [ option(s) ] "
              << "\nThe GEM2MT simulator for simple reactive mass transport models coupled with GEM\n"
              << "Options:\n"
              << "\t-h,\t--help  \t\tshow this help message\n\n"
              // run mode
              << "\t-n,\t--number-of-nodes   NUM   \tInput number of local equilibrium cells (nodes)\n\n"
              << "\t-s,\t-S,\t--sequential         \tSequential reactors chain with waves (moving fluids or solids between reactors) \n"
              << "\t-f,\t-F,\t--flow-through       \tFlow-through reactors chain with fluxes (simultaneous transport then equilibration)\n"
              << "\t-b,\t-B,\t--box-flux           \tBox-flux megasystem (with time-step integration of arbitrary MGP fluxes)\n"
              << "\t-a,\t-A,\t--finite-difference  \t1-D (advection-dispersion) reactive transport using classical finite-difference method\n"
              << "\t-c,\t-C,\t--Crank-Nicolson     \t1-D (advection-dispersion) reactive transport using implicit Crank-Nicolson scheme\n"
              << "\t-w,\t-W,\t--random-walk        \t1-D (advection-diffusion) reactive transport using random-walk method\n"
              << "\t-i,\t--from-file    PATH  \tFrom the file containing the definition of the gem2mt task\n\n"
              // file path
              << "\t-sl,\t--system-lst  PATH   \tThe path to the file that describes the chemical system\n"
              << "\t-dl,\t--dbr-lst     PATH   \tThe path to the file containing the list of system definition files\n"
              << std::endl;
}


int extract_args(int argc, char* argv[], DemoTask& task_data)
{
    int i=0;
    for(i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if((arg == "-h") || (arg == "--help")) {
            show_usage( "demo_gem2mt" );
            return 1;
        }
        else if((arg == "-s") || (arg == "-S") || (arg == "--sequential")) {
            task_data.mode_operation = 'S';
            task_data.ipmfiles_lst_name = "TestS/CalcColumn-dat.lst";
            task_data.dbr_lst_name = "TestS/CalcColumn-dbr.lst";
        }
        else if((arg == "-f") || (arg == "-F") || (arg == "--flow-through")) {
            task_data.mode_operation = 'F';
            task_data.ipmfiles_lst_name = "TestF/CalcColumn-dat.lst";
            task_data.dbr_lst_name = "TestF/CalcColumn-dbr.lst";
        }
        else if((arg == "-b") || (arg == "-B") || (arg == "--box-flux")) {
            task_data.mode_operation = 'B';
            task_data.ipmfiles_lst_name = "TestB/CalcColumn-dat.lst";
            task_data.dbr_lst_name = "TestB/CalcColumn-dbr.lst";
        }
        else if((arg == "-a") || (arg == "-A") || (arg == "--finite-difference")) {
            task_data.mode_operation = 'A';
        }
        else if((arg == "-c") || (arg == "-C") || (arg == "--Crank-Nicolson")) {
            task_data.mode_operation = 'C';
        }
        else if((arg == "-w") || (arg == "-W") || (arg == "--random-walk")) {
            task_data.mode_operation = 'W';
        }
        else if((arg == "-i") || (arg == "--from-file")) {
            if(i + 1 < argc) {
                task_data.gem2mt_filename = argv[++i];
            } else {
                std::cerr << "--from-file option requires one argument." << std::endl;
                return 1;
            }
        }
        else if((arg == "-n") || (arg == "--number-of-nodes")) {
            if(i + 1 < argc) {
                task_data.nnodes = std::stoi(argv[++i]);
            } else {
                std::cerr << "--number-of-nodes option requires one argument." << std::endl;
                return 1;
            }
        }
        else if((arg == "-sl") || (arg == "--system-lst")) {
            if(i + 1 < argc) {
                task_data.ipmfiles_lst_name = argv[++i];
            } else {
                std::cerr << "--system-lst option requires one argument." << std::endl;
                return 1;
            }
        }
        else if((arg == "-dl") || (arg == "--dbr-lst")) {
            if(i + 1 < argc) {
                task_data.dbr_lst_name = argv[++i];
            } else {
                std::cerr << "--dbr-lst option requires one argument." << std::endl;
                return 1;
            }
        }
    }
    if(task_data.ipmfiles_lst_name.empty() || task_data.dbr_lst_name.empty()) {
        std::cerr << "--system-lst and --dbr-lst files required." << std::endl;
        return 1;
    }
    return 0;
}
