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
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from xgem2mt import *

"""
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
"""


def task_S(nnodes: int, ipm_lst: str, dbr_lst: str) -> int:
    try:
        # (1) Create the  MassTransport simulator
        mt_task = MassTransport("S", nnodes)
        mt_task.setName("Test S mode (initial 2.5 mmol/L calcite)  1 to 96 C")

        # (2) Set up sizes, flags and values different from default
        # Set physical time iterator (start,end,step)
        mt_task.setTau(0.0, 30.0, 0.3)
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

        # (3) Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS,
        # allocate and set default values for gem2mt arrays
        mt_task.initNodes(ipm_lst, dbr_lst)
        print(mt_task)  # setup summary: nodes, time step, properties, phases

        # (4) Change/define some other gem2mt arrays
        # Set units for aq phase quantities
        mt_task.setUnitsPhaseQuantities(0, "n")
        for ii in range(mt_task.nMGPfluxes):
            # Set the flux definition: flux order, flux rate, MGP quantity
            mt_task.setFluxDefinition(ii, 0.0, 1.0, 0, 0)

        mt_task.dumpTask("gem2mt_out.dat")

        # (5) Here we call the mass-transport finite-difference coupled routine
        mt_task.runMassTransport()

        return 0

    except RuntimeError as ex:
        print(ex)
        return 1


if __name__ == "__main__":
    task_S(51, "TestS/CalcColumn-dat.lst", "TestS/CalcColumn-dbr.lst")
