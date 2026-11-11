# MassTransport in Python: functions and parameters

This page lists the functions that set up and run a `MassTransport` task. It gives the parameters,
units, permitted values and defaults. The text is for the 1-D modes 'A', 'C' and 'W', unless it
tells you differently.

## Sequence of calls

1. Create the task: `MassTransport(mode, nnodes)`.
2. Set the task options: time, length, solver start, output.
3. Set the transport parameters. These values go to all nodes.
4. If necessary, change which initial system each node gets: `setDistributing`.
5. Read the GEMS3K files: `initNodes`.
6. If necessary, change the compositions and properties of nodes.
7. Set the monitoring function: `setMonitoring`.
8. Run the task: `runMassTransport`.

Do the steps in this sequence. Some functions have no effect if you call them at a different time
(see the "When" column).

To see the setup at any step, use `print(task)` or `task.summary()`. It shows the mode, the column,
the time step that the program will use, the initial systems, and the nodes grouped by equal setup
(system, role, node type, transport properties; after `initNodes` also T, P, mass, volume and phases).
In modes 'S', 'F' and 'B' it also shows the mobile phase groups and the fluxes. Each value has the
function that changes it. If the `setTau` step is too large for the time step, it shows a warning.
`task.dump()` gives the full task as JSON.

## 1. Create the task

| Function | Parameters | Notes |
| --- | --- | --- |
| `MassTransport(mode, nnodes)` | `mode`: one character, see below. `nnodes`: number of nodes. | 50 to 250 nodes usually give a good resolution. |

| Mode | Model |
| --- | --- |
| `'S'` | Chain of reactors. Fluids or solids move from one reactor to the next as waves. |
| `'F'` | Chain of flow-through reactors. Transport, then equilibrium, at each step. |
| `'B'` | Box-flux system. Fluxes of mobile phase groups between boxes, integrated over time. |
| `'A'` | 1-D advection and dispersion. Explicit finite differences. |
| `'C'` | 1-D advection. At this time, it is the 'A' step with Courant number 1 (one node per step) and no dispersion. It is not a Crank-Nicolson scheme. |
| `'W'` | 1-D advection and diffusion with random-walk particles. |

## 2. Task options

| Function | Parameters | Default | Notes |
| --- | --- | --- | --- |
| `setName(text)`, `setComment(text)` | text | | Name and comment of the task. |
| `setTau(start, stop, step)` | times, in the time unit of the velocity (usually s) | 0, 1000, 1 | Modes 'A', 'C', 'W': see "Time step" below. |
| `setSpatialDimensions(x, y=0, z=0)` | `x`: length of the column, m | x = 1 | Node size: dx = x / nnodes. `y` and `z` are not used. |
| `setSIA(flag)` | `'+'`, `'*'` or `'-'` | | Start point of the GEM calculation in each node. `'+'`: the result of the previous step (smart initial approximation, fast). `'*'`: the internal state of the solver. `'-'`: a new start (automatic initial approximation, slow but robust). If the calculation fails, use `'-'`. |
| `setTypeFluxPhase(type)` | `'1'` aqueous, `'2'` gas, `'3'` aqueous + gas, `'4'` solids | `'1'` | The mobile phase. Other values give `'1'`. |
| `setOutput(enable)` | `True` or `False` | `True` | `True` writes the profile files `ICaq-log2.dat`, `ICdif-log2.dat` and `Ph-log2.dat` in the working folder. |
| `setOutVTK(enable)` | `True` or `False` | `False` | `True` writes the nodes in VTK format at each step. |
| `setCutoffICamount(val)` | mol | 1e-9 | If no element in a node changes by more than this value (or by more than 0.1 % of its amount, if that is smaller), the node is not equilibrated again in this step. Not used with `setSIA('-')`: then all nodes are equilibrated at each step. |
| `setCutoffMinimalAmountsIC(val)` | mol | 1e-12 | Smallest amount of an element in a node. Smaller amounts are set to this value. |

## 3. Transport parameters

Call these functions before `initNodes`. `initNodes` copies the values to all nodes.

| Function | Parameter | Unit | Used by mode 'A' |
| --- | --- | --- | --- |
| `setFluidAdvectionVelocity(v)` | pore velocity v | m/s | Advection of all nodes, the time step and dispersion. Must be 0 or more: the fluid flows from node 0 to the last node. |
| `setInitialDispersivity(al)` | longitudinal dispersivity | m | Dispersion. |
| `setInitialDiffusivity(Dif)` | diffusion coefficient in water | m²/s | Dispersion. Typical value: 2e-9. |
| `setInitialNodeEffectivePorosity(eps)` | porosity (0 to 1) | – | Diffusion part of dispersion only. Default 1. |
| `setInitialTortuosityFactor(nto)` | tortuosity factor (1 to 10) | – | Diffusion part of dispersion only. |
| `setInitialTotalNodeVolume(Vt)` | node volume | m³ | Not used. Kept in the node data. |
| `setInitialEffectivePermeability(Km)` | permeability | m² | Not used. Kept in the node data. 1 darcy = 1e-12 m². |
| `setTimeStepReductionFactor(tf)` | factor | – | Divides the time step by `tf`. Values less than 1 give 1. Default 1. |

Dispersion coefficient of a node: D = al · v + Dif · eps / nto.

## Time step (modes 'A' and 'C')

The program calculates the time step dt. It does not use the `step` of `setTau`.

* Mode 'A': dt = min(0.5 · dx / v, dx² / (4 · D)) / tf.
  The first term is the limit for advection (Courant number 0.5). The second term is the limit for
  dispersion. The program uses the largest D of all nodes.
* Mode 'C': dt = dx / v / tf.

The run stops at the first of these two conditions:

* The time is equal to or more than `stop`. The last step can go past `stop` by less than one dt.
* The number of steps is (stop − start) / step + 1.

Thus, set a `step` that is smaller than dt. If `step` is larger than dt, the run stops before `stop`.

Example (calcite–dolomite benchmark): L = 0.5 m, 80 nodes, v = 9.375e-6 m/s, al = 0.0067 m.
dx = 0.00625 m and D = 6.3e-8 m²/s. The advection limit is 333 s. The dispersion limit is 155 s.
Thus dt = 155 s. With `setTau(0, 21000, 10)`, the run does 136 steps and stops at 21144 s.

## 4. Initial systems and nodes

| Function | Parameters | When | Notes |
| --- | --- | --- | --- |
| `setDistributing(node, system)` | node index, initial system index | Before `initNodes` | The initial system that the node gets. |
| `initNodes(dat_lst)` | path of the `*-dat.lst` file | | Reads the GEMS3K files and equilibrates all nodes. The initial systems are the dbr files in the `*-dat.lst` file. |
| `initNodes(dat_lst, dbr_lst)` | paths of the `*-dat.lst` and `*-dbr.lst` files | | The same. The initial systems are the dbr files in the `*-dbr.lst` file. |
| `initNodesFromStrings(dch, ipm, dbr_list)` | JSON strings | | The same as `initNodes`, from strings. |
| `setNodeType(node, type)`, `setNodesType(nodes, type)` | node type, see below | Before `runMassTransport` | Mode 'A' does not use the node type. |

Default distribution in modes 'A' and 'C': nodes 0 and 1 get initial system 0 (the inflowing fluid).
All other nodes get initial system 1. Nodes 0 and 1 do not change during the run (constant source).
The last node is the outlet.

The GEMS3K files:

| File | Contents |
| --- | --- |
| `*-dch.dat` | The chemical system: elements, species, phases, thermodynamic data. |
| `*-ipm.dat` | The settings of the GEM solver. |
| `*-dbr-0-0000.dat`, ... | One composition each: bulk composition, T, P, equilibrium speciation. |
| `*-dat.lst` | `-t` (text files), the dch and ipm file names, then the dbr file names, separated by spaces. Put each name in double quotes; a comma directly after a quoted name is permitted (`"a.dat", "b.dat"`), a comma after a name without quotes is not. |
| `*-dbr.lst` | The dbr file names again, separated by commas. Only for `initNodes(dat_lst, dbr_lst)`. |

The position of a dbr file in the list gives the index of the initial system. The first file is
system 0. The numbers in the file names are not used.

If there is only one dbr file, all nodes start from it. The system numbers (0 for nodes 0 and 1,
1 for the other nodes in mode 'A') then only mark groups. Use `setSystem(0, ...)` and
`setSystem(1, ...)` to give the fluid and the column their compositions.

Mode 'A' and 'C' boundaries: nodes 0 and 1 (inflow) and the last node (outlet) keep their composition
during the run. To change the inflow, change the composition of nodes 0 and 1 (system 0). The node
type is not used in these modes.

Node types: 0 normal; 1 and −1 Dirichlet source and sink; 2 and −2 Neumann source and sink;
3 and −3 Cauchy source and sink; 4 time-dependent input. Default: 3 for nodes 0 and 1, −3 for the
last node, 0 for all other nodes.

## 5. Set up nodes after `initNodes`

Use these functions after `initNodes` and before `runMassTransport`. `nodes` is one index or a
list of indices. A node setup has two parts:

* The **composition** (chemical state): the element amounts, the temperature and the pressure.
  It changes during the run.
* The **transport properties** (physical parameters): the values of section 3 for each node.
  They do not change during the run.

| Function | Parameters | Notes |
| --- | --- | --- |
| `setSystem(system, composition=, T=, P=, transport=)` | see below | Sets up all nodes of an initial system. |
| `setNodes(nodes, composition=, T=, P=, transport=)` | see below | Sets up one node or a list of nodes. |
| `nodesOfSystem(system)` | initial system index | Returns the nodes that got this system. |
| `nodesOfType(type)` | node type | Returns the nodes of this type. |
| `nodeSystem(node)`, `nodeType(node)` | node index | Returns the system index or the node type. |
| `transportProperties(node)` | node index | Returns the `TransportProperties` of the node, with all fields set. |
| `setHydraulicParameters(node, Vt, vp, eps, Km, al, Dif, nto)` | all seven values | Old form for one node. No effect before `initNodes`. |

Parameters of `setSystem` and `setNodes`. Give only the parts that you want to change:

* `composition`: the element amounts in mol. Use one of these:
  * an xGEMS `Material` (the program uses its `b_dict()`),
  * a dict, for example `{"Ca": 3.4e-4, "C": 3.4e-4, "H": 111.0, "O": 55.51}`. Elements that are
    not in the dict get 0,
  * an array in the element order of the chemical system.
* `T`: temperature in K. `P`: pressure in Pa.
* `transport`: a `TransportProperties` object. Only its fields that are set are changed. The permitted values are: porosity more than 0 and not more than 1, tortuosity more than 0, the other fields 0 or more.

If you give `composition`, `T` or `P`, the program equilibrates each node again from the start.

`TransportProperties(volume=, velocity=, porosity=, permeability=, dispersivity=, diffusivity=, tortuosity=)`:

| Field | Section 3 name | Unit | Used by mode 'A' |
| --- | --- | --- | --- |
| `volume` | Vt | m³ | No |
| `velocity` | vp | m/s | Dispersion of the node only. Advection uses `setFluidAdvectionVelocity`. |
| `porosity` | eps | – | Diffusion part of dispersion |
| `permeability` | Km | m² | No |
| `dispersivity` | al | m | Dispersion |
| `diffusivity` | Dif | m²/s | Dispersion |
| `tortuosity` | nto | – | Diffusion part of dispersion |

Example:

```python
from xgem2mt import MassTransport, TransportProperties

task.initNodes("my-dat.lst", "my-dbr.lst")
task.setSystem(0, composition=fluid)                       # nodes 0 and 1 in mode 'A'
task.setSystem(1, composition=column,
               transport=TransportProperties(porosity=0.32, dispersivity=0.0067))
task.setNodes(range(60, 80), transport=TransportProperties(porosity=0.2))
```

An error (`RuntimeError`) occurs if a node index, an element name or the length of the
composition is not correct, or if the equilibrium calculation fails.

## 6. Monitor and run

| Function | Parameters | Notes |
| --- | --- | --- |
| `setMonitoring(f)` | function `f(engine, step, time)` | The program calls `f` after each time step. `engine` is the `ChemicalEngineArray` of the nodes. |
| `getChemicalEngineArray()` | | Returns the `ChemicalEngineArray`. Available after `initNodes`. |
| `runMassTransport()` | | Runs the task. An error in the transport gives a `RuntimeError`. |

Useful `ChemicalEngineArray` functions (the first argument is the node index):
`phaseAmount(n, "Calcite")`, `elementAmountsInPhase(n, "aq_gen")`, `phaseMass(n, "aq_gen")`,
`pH(n)`, `indexElement("Ca")`, `nNodes()`. The concentration in mol/kg is
`elementAmountsInPhase(n, "aq_gen") / phaseMass(n, "aq_gen")`.
