# MassTransport in Python: functions and parameters

This page lists the functions that set up and run a `MassTransport` task. It gives the parameters,
units, permitted values and defaults. The text is for the 1-D modes 'A', 'C' and 'W', unless it
tells you differently.

## Names used on this page

Some names look alike but are different values. This page uses the name in the first column. The
GEMS name is the short name in GEMS files, in `setHydraulicParameters` and in old documents.

| Name | Symbol | Unit | Meaning | GEMS name |
| --- | --- | --- | --- | --- |
| column length | L | m | length of the 1-D column | sizeLc[0] |
| node length | Δx | m | column length / number of nodes; the program calculates it | dx |
| velocity | v | m/s | pore velocity of the water | fVel, vp |
| porosity | φ | – | pore volume / bulk volume (0 to 1) | eps |
| **dispersivity** | α | m | a length of the porous medium: how much the flow spreads the solutes | al |
| **diffusivity** | D<sub>w</sub> | m²/s | diffusion coefficient of the solutes in water | Dif |
| tortuosity | τ | – | how much longer the diffusion path is than the column (1 to 10) | nto |
| **dispersion coefficient** | D | m²/s | D = α·v + D<sub>w</sub>·φ/τ; the program calculates it | hDl |
| time step | Δt | s | the program calculates it (see "Time step") | dTau |
| time step reduction factor | – | – | the time step is divided by this factor | tf |
| node volume | – | m³ | not used in mode 'A' | Vt |
| permeability | – | m² | not used in mode 'A' | Km |

Dispersivity, diffusivity and the dispersion coefficient are three different values. Porosity,
diffusivity and tortuosity change only the diffusion part of the dispersion coefficient.

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

To see the setup at any step, use `print(task)` or `task.summary(width=100, phases=False)`. It shows
the mode, the column, the time step that the program will use, the initial systems, and the nodes
grouped by equal setup (system, role, node type, transport properties; after `initNodes` also T, P,
mass and volume, and the names of the phases present). Lines are not longer than `width`.
`phases=True` adds a table of the phase amounts, with one column for each phase.
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
| `setTypeFluxPhase(type)` | `'1'` aqueous, `'2'` gas, `'3'` aqueous + gas, `'4'` solids | `'1'` | The mobile phase. Only modes 'S', 'F' and 'B' use it. Other values give `'1'`. |
| `setOutput(enable)` | `True` or `False` | `True` | `True` writes the profile files `ICaq-log2.dat`, `ICdif-log2.dat` and `Ph-log2.dat` in the working folder. |
| `setOutVTK(enable)` | `True` or `False` | `False` | `True` writes the nodes in VTK format at each step. |
| `setCutoffICamount(val)` | mol | 1e-9 | If no element in a node changes by more than this value (or by more than 0.1 % of its amount, if that is smaller), the node is not equilibrated again in this step. Not used with `setSIA('-')`: then all nodes are equilibrated at each step. |
| `setCutoffMinimalAmountsIC(val)` | mol | 1e-12 | Smallest amount of an element in a node. Smaller amounts are set to this value. |

## 3. Transport parameters

Call these functions before `initNodes`. `initNodes` copies the values to all nodes.

| Function | Parameter | Unit | Used by mode 'A' |
| --- | --- | --- | --- |
| `setFluidAdvectionVelocity(velocity)` | velocity (pore velocity) | m/s | Advection of all nodes, the time step and dispersion. Must be 0 or more: the fluid flows from node 0 to the last node. |
| `setInitialDispersivity(dispersivity)` | dispersivity | m | Dispersion coefficient. |
| `setInitialDiffusivity(diffusivity)` | diffusivity (diffusion coefficient in water) | m²/s | Diffusion part of the dispersion coefficient. Typical value: 2e-9. |
| `setInitialNodeEffectivePorosity(porosity)` | porosity (0 to 1) | – | Diffusion part of the dispersion coefficient only. Default 1. |
| `setInitialTortuosityFactor(tortuosity)` | tortuosity (1 to 10) | – | Diffusion part of the dispersion coefficient only. |
| `setInitialTotalNodeVolume(volume)` | node volume | m³ | Not used. Kept in the node data. |
| `setInitialEffectivePermeability(permeability)` | permeability | m² | Not used. Kept in the node data. 1 darcy = 1e-12 m². |
| `setTimeStepReductionFactor(factor)` | time step reduction factor | – | The time step is divided by this factor. Values less than 1 give 1. Default 1. |

Dispersion coefficient of a node = dispersivity · velocity + diffusivity · porosity / tortuosity.

## Time step (modes 'A' and 'C')

The program calculates the time step. It does not use the `step` of `setTau`.

* Mode 'A': time step = the smaller of 0.5 · node length / velocity (the advection limit) and
  node length² / (4 · dispersion coefficient) (the dispersion limit), divided by the time step reduction
  factor. The program uses the largest dispersion coefficient of all nodes.
* Mode 'C': time step = node length / velocity / time step reduction factor.

The run stops at the first of these two conditions:

* The time is equal to or more than `stop`. The last step can go past `stop` by less than one time step.
* The number of steps is (stop − start) / step + 1.

Thus, set a `step` that is smaller than the time step. If it is larger, the run stops before `stop`.

Example (calcite–dolomite benchmark): column length 0.5 m, 80 nodes, velocity 9.375e-6 m/s, dispersivity
0.0067 m, diffusivity 0. Node length = 0.00625 m, dispersion coefficient = 6.3e-8 m²/s. The advection limit
is 333 s. The dispersion limit is 155 s. Thus the time step is 155 s. With `setTau(0, 21000, 10)`, the run does 136 steps and stops at 21144 s.

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
  * an xGEMS `Material` (the program uses its `b_dict()`). An element of the system with a longer name, such as
    `Inert` (an inert solid) or `N_atm`, can be used in its formulas with xgems from 2026-11-12 (`1f9b319`) or later:
    `sand.add("Inert", 70, "g")`. With older xgems, give it as an element amount: `add({"Inert": 1.0}, 70, "g")`,
  * a dict, for example `{"Ca": 3.4e-4, "C": 3.4e-4, "H": 111.0, "O": 55.51}`. Elements that are
    not in the dict get 0,
  * an array in the element order of the chemical system.
* `T`: temperature in K. `P`: pressure in Pa.
* `transport`: a `TransportProperties` object. Only its fields that are set are changed. The permitted values are: porosity more than 0 and not more than 1, tortuosity more than 0, the other fields 0 or more.

In place of `composition`, you can give the node as parts, as measured in a column experiment. The program
then calculates the composition of each node (see section 8):

* `pore_solution`: the pore solution (the aqueous phase): a `Material`, a dict or an array, any amount. It is
  the reference of each node.
* `water_mass`: kg of liquid water (the H2O solvent) in each node. The program scales `pore_solution` to this
  water. Default: the water in `pore_solution` as given (from a GEM calculation at the T and P of the node).
  In the modes 'A', 'C' and 'W', all nodes need the same water mass. The first call with `pore_solution` sets
  it; later calls without `water_mass` scale their pore solution to it. `initNodes` resets it.
* `water_saturation`: the share of the pore volume that the pore solution fills. Default 1.
* `parts`: a list of `Part` objects, the other phases and solids of the node. `of` says which volume the part
  fills: `"pore"` (a gas, another liquid) or `"solid"`. Then give a share of that volume or a mass:
  * `Part(material, fraction=f, of="pore")`: the share f of the pore volume;
  * `Part(material, fraction=f, of="solid")`: the share f of the solid volume that the parts given by mass leave
    (the solid fractions add up to 1);
  * `Part(material, mass=m, of="solid")` (or `of="pore"`): m kg in the node. Its volume is taken from the
    `of` volume.

  The material can have any amount: the program scales it. Add `density=` (kg/L) to give the volume of a part;
  otherwise it comes from a GEM calculation (see section 8).
* `solids`, `grain_density`: a short form for one solid part: `parts=[Part(solids, fraction=1, of="solid",
  density=grain_density)]`. `grain_density` is optional.

Give either `composition`, or `pore_solution` and the parts, not both.

If you give a composition, `T` or `P`, the program equilibrates each node again from the start.

`TransportProperties(volume=, velocity=, porosity=, permeability=, dispersivity=, diffusivity=, tortuosity=)`:

| Field | GEMS name | Unit | Used by mode 'A' |
| --- | --- | --- | --- |
| `volume` | Vt | m³ | No |
| `velocity` | vp | m/s | Dispersion coefficient of the node only. Advection uses `setFluidAdvectionVelocity`. |
| `porosity` | eps | – | Diffusion part of the dispersion coefficient |
| `permeability` | Km | m² | No |
| `dispersivity` | al | m | Dispersion coefficient |
| `diffusivity` | Dif | m²/s | Diffusion part of the dispersion coefficient |
| `tortuosity` | nto | – | Diffusion part of the dispersion coefficient |

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
composition is not correct, a transport value is not permitted (it must be a finite number; the velocity 0 or
more, the porosity more than 0 and not more than 1, the tortuosity more than 0, the others 0 or more), or if
the equilibrium calculation fails.

A call changes all the nodes that you give or none: if one node fails, the nodes that the call already
changed get their earlier transport values and their exact earlier state back.

A call that changes the chemistry writes a short progress text. It has one line for each node with GEM
calculations, with the time, and one line for the nodes that get the state of an equal node. The GEM
calculations of one node can take minutes, so this text shows that the setup is not stuck.
`task.setProgressMessages(False)` hides the text.

```text
Setting up 28 nodes (2-29): composition and equilibrium (GEM)
  node 2 ... 2.3 s
  nodes 3-29: the same inputs as node 2, state copied (no GEM)
  done in 2.3 s
```

The transport properties of nodes exist only in the modes 'A', 'C' and 'W'. In the modes 'S', 'F' and 'B',
`transport=` gives an error, and `pore_solution` and `parts` use the initial porosity
(`setInitialNodeEffectivePorosity`) of the task.

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

## 7. Example: one exported GEMS3K system, a material for each zone

Use this procedure if you export only one system from GEMS (one dbr file) and give each zone of the
column its composition in Python. The repository has such an export: `demos/TestAD/CalcColumn-dat.lst`
(calcite-dolomite system, 298.15 K, 1 bar); you can use it in place of `MySystem-dat.lst`.

Before you export the system from the GEMS GUI, make sure that:

* The chemical system has all the elements, species and phases that your materials need. A material
  can use only the elements in the dch file.
* The temperatures and pressures of the export include the T and P of all nodes. The dch file has
  thermodynamic data only for these T and P points. At other values, the equilibrium calculation fails.

The composition in the dbr file is only the start value. `setSystem` and `setNodes` replace it.

```python
from xgems import ChemicalEngine, Material
from xgem2mt import MassTransport, TransportProperties

# (1) The task: mode 'A', 80 nodes
task = MassTransport("A", 80)
task.setTau(0.0, 21000.0, 10.0)            # s; the setTau step must be smaller than the time step
task.setSpatialDimensions(0.5)             # column length, m
task.setSIA("+")
task.setOutput(False)

# (2) Transport parameters for all nodes
task.setFluidAdvectionVelocity(9.375e-6)   # m/s, one value for the column
task.setInitialNodeEffectivePorosity(0.32)
task.setInitialDispersivity(0.0067)        # m
task.setInitialDiffusivity(2e-9)           # m2/s
task.setInitialTortuosityFactor(1.0)

# (3) Groups of nodes. Default in mode 'A': nodes 0 and 1 = group 0 (inflow), other nodes = group 1.
#     Give a third zone its own group number before initNodes.
for q in range(40, 60):
    task.setDistributing(q, 2)

# (4) Read the exported system. With one dbr file, all nodes start from it.
task.initNodes("MySystem-dat.lst")

# (5) The materials. The xGEMS engine gives the element list and the molar masses.
gems = ChemicalEngine("MySystem-dat.lst")
fluid = Material(gems, "MgCl2 solution")
fluid.add("H2O", 1.0, "kg")
fluid.add("MgCl2", 1e-3, "mol")
fluid.add("O2", 1e-5, "mol")              # a little O2: see the notes
rock = Material(gems, "calcite rock")
rock.add("H2O", 1.0, "kg")
rock.add("CaCO3", 3.4e-4, "mol")
rock.add("O2", 1e-4, "mol")
clay = Material(gems, "clay layer")
clay.add("H2O", 1.0, "kg")
clay.add("CaCO3", 1e-4, "mol")
clay.add("O2", 1e-4, "mol")

# (6) Set up the groups. A TransportProperties object changes only the values that you give;
#     the other values stay as in step (2).
task.setSystem(0, composition=fluid)                    # inflow: nodes 0 and 1
task.setSystem(1, composition=rock)                     # column
task.setSystem(2, composition=clay,                     # zone: nodes 40-59
               transport=TransportProperties(porosity=0.2, dispersivity=0.01))

# (7) Check the setup, then run
print(task)
task.runMassTransport()
```

Notes:

* Step (2) sets the transport parameters of all nodes. In step (6), `transport=` changes only the
  fields that you give for that group. In the example, nodes 40-59 have porosity 0.2 and dispersivity
  0.01 m. They keep the diffusivity 2e-9 m2/s and the tortuosity 1 of step (2).
* The advection velocity is one value for the column (`setFluidAdvectionVelocity`). A `velocity` in
  `TransportProperties` changes only the dispersion of these nodes.
* With one dbr file, the group numbers (0, 1, 2) are only labels for `setSystem`. You can also set any
  list of nodes without a group: `task.setNodes(range(70, 75), composition=rock, T=310.0)`.
* In mode 'A', nodes 0 and 1 are always the inflow and the last node is always the outlet. Give them the
  composition that you want at the boundaries.
* `print(task)` shows the groups, their roles and their values before the run. If the `setTau` step is
  too large, it shows a warning.
* Add a little O2 (or another redox buffer) to each material. Without it, the material has exactly the
  O and H of water, and the redox state is not defined. Then many GEM calculations fail at the first
  attempt: GEMS3K repeats them (from the previous result, with a small change of the composition, or
  with the other solver) and usually gets a result, but the run is slower (in the example 3 times) and
  `ipmlog.txt` and `xGEM2MT.log` show many errors.
* If a GEM calculation of a node gives a failure (ERR) or a result that is not trustful (BAD), xGEM2MT tries
  two rescues. First, it calculates the node again from the last converged result of this node (SIA, with
  the new bulk composition). If this also fails, and no element amount changed by more than the mass balance
  tolerance of GEM (`pa_DHB`, 1e-13 relative) since that result, the node keeps that result. `gemFailures().retried`
  and `gemFailures().kept` count these rescues; `xGEM2MT.log` has a line for each.
* After the run, `task.gemFailures()` gives the GEM calculations that failed also after these repeats and
  rescues (an ERR node keeps its speciation from before the step, a BAD node keeps the BAD result; the new
  bulk composition is used for the next time step). Python
  gives a `RuntimeWarning` if there are such failures, and `print(task)` shows them. `gemFailures().ran` is
  `False` before the first run: then there is no result, and `count` 0 does not mean that all calculations
  converged. The calculations of `initNodes` and `setNodes` are not counted. The messages in `ipmlog.txt`
  give the node number (`Error Node:47`).

## 8. From a column experiment to the model (mode 'A')

| Experiment | Model value | How |
| --- | --- | --- |
| column length L (m) | column length | `setSpatialDimensions(L)` |
| (resolution) | number of nodes n | `MassTransport("A", n)`. Rule of thumb: node length L/n not more than 2 × dispersivity. |
| flow rate Q (m³/s), cross-section area A (m²), porosity φ | velocity | v = Q / (A · φ); `setFluidAdvectionVelocity(v)` |
| porosity φ | porosity | pore volume / bulk volume, for example mass of water in the saturated column / (water density · A · L); `setInitialNodeEffectivePorosity(φ)` |
| tracer breakthrough curve (Cl⁻, Br⁻) | dispersivity | fit the curve; `setInitialDispersivity` |
| (literature) | diffusivity, tortuosity | about 2e-9 m²/s, 1 to 2 for sand; `setInitialDiffusivity`, `setInitialTortuosityFactor` |
| duration t (s) | time | `setTau(0, t, step)` with `step` smaller than the time step; pore volumes = v · t / L |
| inflowing solution (mol/kg water) | composition of nodes 0 and 1 | 1 kg of water and the dissolved amounts per kg of water |
| initial pore water (mol/kg water) | composition of the column nodes | 1 kg of water and the dissolved amounts per kg of water |
| solids packed in the column: minerals and their volume shares (or masses) | `parts` with `of="solid"` | the program scales them to each node (see below); a gas or another liquid: `of="pore"` |
| effluent | node n − 2 | the last node keeps its start composition (outlet); read the effluent from node n − 2 |

What each node needs:

* Nodes 0 and 1 (inflow): only the composition of the inflowing solution. Their transport values are not used:
  the transport does not change these nodes.
* Column nodes: the composition per kg of water (pore water and minerals) and the transport values
  (porosity, dispersivity, diffusivity, tortuosity).
* All nodes: the same mass of water, for example 1 kg. The absolute size of a node does not matter, only the
  amounts per kg of water. An inert solid is not necessary (it can make GEM faster, see the notebooks).

Two ways to give the compositions:

1. **Scaled by the program** (`pore_solution`, `water_mass`, `water_saturation`, `parts`). For each node:
   * pore volume = volume of the pore solution (water mass, with its solutes) / water saturation;
   * solid volume = pore volume · (1 − porosity) / porosity, with the transport porosity of the node;
   * the parts given by mass take their volume from their `of` volume; the parts given by a fraction fill their
     share of the pore volume, or of the solid volume that is left; mass = volume / (volume per kg).

   The volume per kg of a part comes from GEM at the T and P of the node, with the models of the chemical
   system (for example a non-ideal gas): the volume that 1 kg of the part adds to the pore solution. A part
   without water, such as a dry solid or a gas, does not converge alone. If you give a `density`, the program
   uses it. A gas that dissolves in part in the pore solution fills less volume than planned: the program then
   increases the fluid parts and calculates again, until the node volume is pore volume + solid volume.

   Thus a zone with a different porosity or saturation gets the correct amounts, and you give the porosity
   only once. The porosity of the chemical system (water and gas volume / node volume, which GEMS calculates)
   is then equal to the transport porosity, if the parts give the complete solid.

   ```python
   from xgem2mt import Part
   pore = Material(gems, "pore water"); pore.add("H2O", 1.0, "kg"); pore.add("O2", 1e-4, "mol")
   quartz = Material(gems, "quartz");   quartz.add("SiO2", 1.0, "mol")      # Si must be in the system
   calcite = Material(gems, "calcite"); calcite.add("CaCO3", 1.0, "mol")
   air = Material(gems, "air");         air.add("O2", 0.21, "mol"); air.add("N2", 0.79, "mol")  # N must be in the system
   task.setSystem(0, pore_solution=inflow)                                  # inflow: only the solution
   task.setSystem(1, pore_solution=pore, water_saturation=0.9,
                  parts=[Part(air, fraction=0.1, of="pore"),
                         Part(quartz, fraction=0.99, of="solid"),
                         Part(calcite, fraction=0.01, of="solid", density=2.71)])
   ```

   Mode 'A' moves only the aqueous phase. A gas or another liquid stays in its node: it reacts with the pore
   water there, but it does not move, and the saturation does not change the flow.

2. **Given by you** (`composition`): the program uses the amounts exactly as you give them, as the contents of
   one node. It does not scale them and it does not use the porosity for them. You must:
   * put the same mass of water in all nodes (for example 1 kg);
   * give the minerals per kg of that water: mol of mineral per kg of pore water =
     mol of mineral in the column / (φ · A · L · water density);
   * use the same porosity here as in the transport values.

   Example: a column of 1 L bulk volume with porosity 0.3 holds 0.3 kg of water and 0.01 mol of calcite.
   Per kg of water that is 0.01 / 0.3 = 0.033 mol of calcite. Thus each node gets 1 kg of water and 0.033 mol
   of calcite: `composition=` a `Material` with `add("H2O", 1, "kg")` and `add("CaCO3", 0.033, "mol")`.

If the parts give only the reactive minerals and not the complete solid, the porosity of the chemical system
is larger than the transport porosity. This does not change the transport.

The porosity is used twice: in the dispersion coefficient (diffusion part only), and to convert the mineral
content of the column to mol per kg of pore water. A wrong porosity gives a wrong mineral-to-water ratio.

Use only one `MassTransport` task at a time in a Python process. At this time, two tasks that exist at the
same time use the same internal data, and the results of the first task are wrong.
