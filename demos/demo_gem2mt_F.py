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

def python_logger(engine: ChemicalEngineArray, step: int, time: float):
    print(f"[Python Logger] step={step} time={time} size={engine.nNodes()}.")
    bj = engine.elementAmounts(3)
    print(bj)

"""
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
"""


def task_F(nnodes: int, ipm_lst: str, dbr_lst: str) -> int:
    try:
        # (1) Create the  MassTransport simulator
        mt_task = MassTransport("F", nnodes)
        mt_task.setName("Test F mode (flow-through reactors)")

        # (2) Set up sizes, flags and values different from default
        # Set physical time iterator (start,end,step)
        mt_task.setTau(0.0, 1200.0, 1.0)
        # Set spatial dimensions of the medium defines topology of nodes ( x y z )
        mt_task.setSpatialDimensions(0, 0, 0)
        # Set number of mobile groups of phases
        mt_task.setNumberPhaseGroups(1)
        # Set number of MGP fluxes defined in the megasystem
        mt_task.setNumberMGPfluxes(nnodes)
        # Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
        mt_task.setSIA("-")
        # Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids )
        mt_task.setTypeFluxPhase("1")
        # Use non stop debug output for nodes
        mt_task.setOutput(False)

        # Set external monitoring
        mt_task.setMonitoring(python_logger)

        # (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
        # allocate and set default values for gem2mt arrays
        mt_task.initNodes(ipm_lst, dbr_lst)

        # (4) Change/define some other gem2mt arrays
        for ii in range(mt_task.nPhaseGroups):
            # Set ID of mobile phase group
            mt_task.setPhaseGroupID(ii, "phg" + str(ii + 1))
            for k in range(mt_task.nPhases):
                # Set quantities of phases in MGP
                mt_task.setPhaseGroupQuantities(ii, k, 1.0 if k == 0 else 0.0)

        for k in range(mt_task.nPhases):
            # Set units for setting phase quantities in MGP (see PGT )
            mt_task.setUnitsPhaseQuantities(k, "M")

        for ii in range(mt_task.nMGPfluxes):
            # Set Source/Receive box index in the flux definition
            mt_task.setFluxSourceReceive(
                ii, ii, ii + 1 if ii < mt_task.nMGPfluxes - 1 else -1
            )
            # Set the flux defnition: flux order, flux rate, MGP quantity
            mt_task.setFluxDefinition(ii, 1.0, 0.5, 0, 0)
            #  Set the ID of MGP to move in this flux
            mt_task.setFluxMGPid(ii, "phg1")
            #  Set IDs of fluxes
            mt_task.setFluxID(ii, "qj")

        # Set list of selected fields and indexes to VTK format
        mt_task.setVTKfields([[33, 0], [32, 0], [38, 1], [38, 4]])
        mt_task.setOutVTK(False)  # allocated but no write

        print(mt_task)  # setup summary after the last change: nodes, time step, properties, fluxes
        mt_task.dumpTask("gem2mt_out.dat")

        # (5) Here we call the mass-transport finite-difference coupled routine
        mt_task.runMassTransport()

        return 0

    except RuntimeError as ex:
        print(ex)
        return 1


if __name__ == "__main__":
    task_F(51, "TestF/CalcColumn-dat.lst", "TestF/CalcColumn-dbr.lst")
