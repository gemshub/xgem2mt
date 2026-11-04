
## xGEM2MT 
Linux, OSX, Windows

The **xGEM2MT** module is a tool for automation of one-dimensional (1-D) reactive mass transport simulations coupled with [GEMS3K](https://github.com/gemshub/GEMS3K) calculation of equilibrium states in spatially distributed nodes (boxes, volumes) over time steps. 

The **xGEM2MT** module can simulate reactive transport in three main modes:

* Sequential reactor chain and waves ('S');
* Flow-through and box-flux sequence ('F', 'B');
* One-dimensional reactive transport with advection/dispersion/diffusion ('A', 'C', 'W')

### License
* xGEM2MT is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* xGEM2MT is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
* You should have received a copy of the GNU Lesser General Public License along with xGEM2MT code. If not, see http://www.gnu.org/licenses/. 

### Installation and Build

#### Installation using Conda

xGEM2MT is not yet on conda-forge. Install it with [Conda](https://github.com/conda-forge/miniforge) from the xGEM2MT channel:

```bash
conda create -n xgem2mt -c https://gemshub.github.io/xgem2mt/stable -c conda-forge xgem2mt
conda activate xgem2mt
```

For the latest test build, use `https://gemshub.github.io/xgem2mt/prerelease` instead of `.../stable`.

Packages are available for Linux, macOS on Apple Silicon and Windows (not for Intel Macs).

#### How to build from source code

* In your home directory, make a folder named e.g. `~/xgem2mt`, open a terminal, and type the command line (do not forget a period) to download the actual "master" branch:
  
```sh
git clone https://github.com/gemshub/xgem2mt . 
```

* Install Dependencies by executing (if not using a Conda environment)

```sh
sudo ./install-dependencies.sh
```

* If using a Conda environment (recommended), create and activate the `xgem2mt` environment by executing:

```sh
conda devenv
conda activate xgem2mt
```

* Build the xgem2mt library and examples

```sh
cd ~/xgem2mt
mkdir build
cd build
cmake .. 
make
sudo make install
```

### Code Examples

For detailed information about how to use this project through its C++ and Python interface [View Documentation](https://xgems.readthedocs.io)

#### C++ API example.

The main functionality of xgem2mt is available through the `MassTransport` class. Its usage is demonstrated below using C++:

```cpp

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
    // Use non-stop debug output for nodes
    mt_task.setOutput(true);

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

    // (5) Here we call the mass-transport finite-difference coupled routine
    mt_task.runMassTransport();

    return 0;
}

```

> [!NOTE]
> You can generate GEMS system files using the print command in the GEM2MT module in the GUI, or prepare them in any other way.
>
> The file containing the list of system definition files must contain more than one file.

##### Python API example

The equivalent Python version is shown below:

```python3

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
        for ii in range(mt_task.nMGPfluxes):
            # Set the flux definition: flux order, flux rate, MGP quantity
            mt_task.setFluxDefinition(ii, 1.0, 0.5, 0, 0)
        
        # (5) Here we call the mass-transport finite-difference coupled routine
        mt_task.runMassTransport()

        return 0

    except RuntimeError as ex:
        print(ex)
        return 1


```

For more information on APIs and more examples, take a look at the `/demos` directory.

