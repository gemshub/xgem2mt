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
'C'     1-D (advection-dispersion) reactive transport using implicit Crank-Nicolson scheme

In this mode, the problem setup is the same as in 'A' mode, with the only difference that the simulation
is performed using the Crank-Nicholson scheme, believed to be more robust regarding to numerical dispersion and/or
error accumulation. Both 'A' and 'C' methods are only applicable to advection-dominated cases.

 See GUI example: CalDolCol2:G:CalcColumn:0:0:1:25:0:1D-Cc-MgCl2-tmod:C:
"""


def task_C(nnodes: int, ipm_lst: str, dbr_lst: str) -> int:
    try:
        # (1) Create the  MassTransport simulator
        mt_task = MassTransport("C", nnodes)
        mt_task.setName("Test of 1D coupled advection problem (dissolved Ca, Mg)")

        # (2) Set up sizes, flags and values different from default
        # Set physical time iterator (start,end,step)
        mt_task.setTau(0.0, 30000.0, 3.0)
        # Set spatial dimensions of the medium defines topology of nodes ( x y z )
        mt_task.setSpatialDimensions(0.2, 0.0, 0.0)
        # Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
        mt_task.setSIA("-")
        # Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
        mt_task.setTypeFluxPhase("1")
        # Use non stop debug output for nodes
        mt_task.setOutput(True)

        # Set advection/diffusion mass transport: time step reduction factor (usually 1)
        mt_task.setTimeStepReductionFactor(2.0)
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

        # (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
        # allocate and set default values for gem2mt arrays
        mt_task.initNodes(ipm_lst, dbr_lst)

        # (4) Change/define some other gem2mt arrays
        # Set list of selected fields and indexes to VTK format
        mt_task.setVTKfields([[33, 0], [32, 0], [38, 1], [38, 4]])

        mt_task.dumpTask("gem2mt_out.dat")

        # (5) Here we call the mass-transport finite-difference coupled routine
        mt_task.runMassTransport()

        return 0

    except RuntimeError as ex:
        print(ex)
        return 1


if __name__ == "__main__":
    task_C(101, "TestAD/CalcColumn-dat.lst", "TestAD/CalcColumn-dbr.lst")
