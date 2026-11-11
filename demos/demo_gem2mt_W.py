# xGEMS2MT is a C++ and Python library for mass transport modeling with Gibbs energy minimization
#
# Copyright (C) 2026 S.Dmytriieva, Dmitrii Kulik
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://gnu.org>.

from xgem2mt import *

"""
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
"""


def task_W(nnodes: int, ipm_lst: str, dbr_lst: str) -> int:
    try:
        # (1) Create the  MassTransport simulator
        mt_task = MassTransport("W", nnodes)
        mt_task.setName("Test of 1D coupled advection problem (dissolved Ca, Mg)")

        # (2) Set up sizes, flags and values different from default
        # Set number of allocated particle types < 20
        mt_task.setNumberParticleTypes(1)
        # Set physical time iterator (start,end,step)
        mt_task.setTau(0.0, 1000000.0, 1000.0)
        # Set spatial dimensions of the medium defines topology of nodes ( x y z )
        mt_task.setSpatialDimensions(0.2, 0.0, 0.0)
        # Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
        mt_task.setSIA("-")
        # Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
        mt_task.setTypeFluxPhase("1")
        # Use non stop debug output for nodes
        mt_task.setOutput(True)

        # Set advection/diffusion mass transport: time step reduction factor (usually 1)
        mt_task.setTimeStepReductionFactor(5.0)
        # Set initial total node volume (m^3)
        mt_task.setInitialTotalNodeVolume(0.001)
        # Set fluid advection velocity (m/sec)
        mt_task.setFluidAdvectionVelocity(2e-6)
        # Set initial node effective porosity (0 < eps < 1), usually 1
        mt_task.setInitialNodeEffectivePorosity(0.5)
        # Set initial effective permeability, m2, usually 1
        mt_task.setInitialEffectivePermeability(1e-12)
        # Set initial value of specific longitudinal dispersivity (m), usually 1e-3
        mt_task.setInitialDispersivity(0.001)
        # Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
        mt_task.setInitialDiffusivity(2e-09)
        # Set initial tortuosity factor, usually 1
        mt_task.setInitialTortuosityFactor(1.0)
        # Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
        mt_task.setCutoffICamount(1e-9)
        # Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
        mt_task.setCutoffMinimalAmountsIC(1e-11)

        # If need update input node distributing
        # Change index of initial system variant for the center node
        mt_task.setDistributing(int(nnodes / 2), 0)
        # Change the type for for the center node
        mt_task.setNodeType(int(nnodes / 2), 3)

        # (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
        # allocate and set default values for gem2mt arrays
        mt_task.initNodes(ipm_lst, dbr_lst)
        print(mt_task)  # setup summary: nodes, time step, properties, phases

        # (4) Change/define some other gem2mt arrays
        # Set list of selected fields and indexes to VTK format
        mt_task.setVTKfields([[33, 0], [32, 0], [38, 1], [38, 4]])
        # Set of particle statistic property
        mt_task.setParticle(0, 1000, 500, 1500, [0, 11, 21, 0, 0, 0])
        # If need update initial hydraulic parameters in node: Vt, vp, eps, Km, al, Dif,  nto
        # mt_task.setHydraulicParameters(1, 0.0015, 2.1e-06, 0.6, 1e-11, 0.0015, 3e-09, 1)

        mt_task.dumpTask("gem2mt_out.dat")

        # (5) Here we call the mass-transport finite-difference coupled routine
        mt_task.runMassTransport()

        return 0

    except RuntimeError as ex:
        print(ex)
        return 1


if __name__ == "__main__":
    task_W(101, "TestAD/CalcColumn-dat.lst", "TestAD/CalcColumn-dbr.lst")
