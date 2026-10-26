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
'B'     Box-flux megasystem (with time-step integration of complex arbitrary MGP fluxes)

This scheme is similar to the F scheme, but allows the most general setup of the box-flux dynamic transport
scheme, not only a unidirectionally connected chain of boxes. In the B scheme, shortcut fluxes and loops are allowed,
as well as zero-, first- and second-order fluxes, box (point) sinks and box (point) sources of given elemental composition.
The transport step is performed by the built-in integrator of a large  ODE system of time derivatives of masses of
each chemical element in each box (reservoir). After the transport step, chemical equilibria are re-calculated in all boxes.
Normally, with time (number of steps), such complex dynamic systems evolve to a steady state when masses and compositions of
reservoirs do not change anymore. This scheme may be good for simulating water-sediment columns or geochemical cycles.

 See GUI example: CalDolCol2:G:CalcColumn:0:0:1:25:0:Test4:B:
"""


def task_B(nnodes: int, ipm_lst: str, dbr_lst: str) -> int:
    try:
        # (1) Create the  MassTransport simulator
        mt_task = MassTransport("B", nnodes)
        mt_task.setName("Test B mode (flow-through reactors), T = 1 to 99 C")

        # (2) Set up sizes, flags and values different from default
        # Set physical time iterator (start,end,step)
        mt_task.setTau(0.0, 10000.0, 1.0)
        # Set number of mobile groups of phases
        mt_task.setNumberPhaseGroups(1)
        # Set number of MGP fluxes defined in the megasystem
        mt_task.setNumberMGPfluxes(nnodes)
        # Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
        mt_task.setSIA("-")
        # Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
        mt_task.setTypeFluxPhase("1")
        # Use non stop debug output for nodes
        mt_task.setOutput(True)

        # Change from default: Other nodes normal, the last two set as constant source and sink
        mt_task.setNodeType(nnodes - 2, 0)

        # (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
        # allocate and set default values for gem2mt arrays
        mt_task.initNodes(ipm_lst, dbr_lst)

        # (4) Change/define some other gem2mt arrays
        # Set ID of the first mobile phase group
        mt_task.setPhaseGroupID(0, "Pg1")
        for ii in range(mt_task.nMGPfluxes):
            #  Set the ID of MGP to move in this flux
            mt_task.setFluxMGPid(ii, "Pg1")

        # Change from default: flux order 1 (proportional to source MPG mass)
        mt_task.setFluxDefinition(0, 1.0, 0.1, 0, 0)

        mt_task.dumpTask("gem2mt_out.dat")

        # (5) Here we call the mass-transport finite-difference coupled routine
        mt_task.runMassTransport()

        return 0

    except RuntimeError as ex:
        print(ex)
        return 1


if __name__ == "__main__":
    task_B(51, "TestB/CalcColumn-dat.lst", "TestB/CalcColumn-dbr.lst")
