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
"""


def task_A(nnodes: int, ipm_lst: str, dbr_lst: str) -> int:
    try:
        # (1) Create the  MassTransport simulator
        mt_task = MassTransport("A", nnodes)
        mt_task.setName("Test of 1D coupled advection problem (dissolved Ca, Mg)")

        # (2) Set up sizes, flags and values different from default
        # Set physical time iterator (start,end,step)
        mt_task.setTau(0.0, 500000.0, 5000.0)
        # Set spatial dimensions of the medium defines topology of nodes ( x y z )
        mt_task.setSpatialDimensions(0.2, 0.0, 0.0)
        # Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
        mt_task.setSIA("-")
        # Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
        mt_task.setTypeFluxPhase("1")
        # Use non stop debug output for nodes
        mt_task.setOutput(True)

        # Set advection/diffusion mass transport: time step reduction factor (usually 1)
        mt_task.setTimeStepReductionFactor(1.0)
        # Set initial total node volume (m^3)
        mt_task.setInitialTotalNodeVolume(0.001)
        # Set fluid advection velocity (m/sec)
        mt_task.setFluidAdvectionVelocity(1e-7)
        # Set initial node effective porosity (0 < eps < 1), usually 1
        mt_task.setInitialNodeEffectivePorosity(1.0)
        # Set initial effective permeability, m2, usually 1
        mt_task.setInitialEffectivePermeability(1.0)
        # Set initial value of specific longitudinal dispersivity (m), usually 1e-3
        mt_task.setInitialDispersivity(0.001)
        # Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
        mt_task.setInitialDiffusivity(0.0)
        # Set initial tortuosity factor, usually 1
        mt_task.setInitialTortuosityFactor(1.0)
        # Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
        mt_task.setCutoffICamount(1e-9)
        # Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
        mt_task.setCutoffMinimalAmountsIC(1e-11)

        # (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
        # allocate and set default values for gem2mt arrays
        mt_task.initNodes(ipm_lst, dbr_lst)
        print(mt_task)  # setup summary: nodes, time step, properties, phases
        mt_task.dumpTask("gem2mt_out.dat")

        # (4) Here we call the mass-transport finite-difference coupled routine
        mt_task.runMassTransport()

        return 0

    except RuntimeError as ex:
        print(ex)
        return 1


if __name__ == "__main__":
    task_A(201, "TestAD/CalcColumn-dat.lst", "TestAD/CalcColumn-dbr.lst")
