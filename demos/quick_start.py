# Quick start of the README: a 1-D calcite-dolomite column ('A' mode).
# Run it from the repository root:  python demos/quick_start.py

from xgems import ChemicalEngine, Material
from xgem2mt import MassTransport, TransportProperties

system = "demos/TestAD/CalcColumn-dat.lst"        # one exported GEMS3K system

task = MassTransport("A", 80)                     # 1-D advection-dispersion, 80 nodes
task.setTau(0.0, 21000.0, 10.0)                   # time: start, stop, step (s)
task.setSpatialDimensions(0.5)                    # column length (m)
task.setFluidAdvectionVelocity(9.375e-6)          # m/s
task.setInitialNodeEffectivePorosity(0.32)        # transport values for all nodes
task.setInitialDispersivity(0.0067)               # m
task.setInitialDiffusivity(2e-9)                  # m2/s
task.setInitialTortuosityFactor(1.0)
task.initNodes(system)                            # one dbr file: all nodes start from it

gems = ChemicalEngine(system)                     # gives Material the elements and molar masses
fluid = Material(gems, "MgCl2 solution")
fluid.add("H2O", 1.0, "kg")
fluid.add("MgCl2", 1e-3, "mol")
fluid.add("O2", 1e-5, "mol")                      # a little O2 keeps the system oxidizing (faster)
rock = Material(gems, "calcite rock")
rock.add("H2O", 1.0, "kg")
rock.add("CaCO3", 1e-3, "mol")
rock.add("O2", 1e-4, "mol")

task.setSystem(0, composition=fluid)              # inflow: nodes 0 and 1
task.setSystem(1, composition=rock)               # the column
task.setNodes(range(40, 60),                      # a zone with other transport values
              transport=TransportProperties(porosity=0.2))

print(task)                                       # the setup at a glance
task.runMassTransport()
print("Run finished:", task.gemFailures())         # count=0: all GEM calculations converged
