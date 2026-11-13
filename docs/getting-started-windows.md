# Getting started on Windows

This guide takes a new user from a PC with nothing installed to a running xGEM2MT model.
It needs about 30 minutes and 3 GB of disk space. You do not build any C++ code: Conda installs
ready-made packages.

## What you need

- Windows 10 or 11, 64 bit, and an internet connection.
- No administrator rights are needed if you install Miniforge "for me only".
- Use a short folder path without spaces or special characters, for example `C:\work\xgem2mt`.

## 1. Install Miniforge (Conda)

1. Download the Windows installer from <https://github.com/conda-forge/miniforge#download>
   (`Miniforge3-Windows-x86_64.exe`).
2. Run it. Select "Just Me". Keep the other defaults.
3. Open **Miniforge Prompt** from the Start menu. Use this window for all commands below.

Check it:

```bat
conda --version
```

## 2. Create the environment

xGEM2MT is not yet on conda-forge. The packages come from two channels of the gemshub project
(the `prerelease` channels carry the newest build, with the node setup functions in this guide).

```bat
conda create -n xgem2mt -c https://gemshub.github.io/xgem2mt/prerelease -c https://gemshub.github.io/xgems/prerelease -c conda-forge python=3.12 xgem2mt xgems jupyter matplotlib pandas
conda activate xgem2mt
```

`xgems` is optional for xGEM2MT itself. The examples use its `Material` class to write a composition
from formulas and amounts (`1 kg H2O`, `1e-3 mol MgCl2`).
If Conda finds an old xGEM2MT build, clean the index cache and try again: `conda clean --index-cache`.

Check the installation:

```bat
python -c "import xgem2mt, xgems; print('ok')"
```

The line must print `ok`. An error here is the most important thing to report (see "Problems" below).

## 3. Get the examples

Install Git for Windows (<https://git-scm.com/download/win>) if you do not have it, then:

```bat
mkdir C:\work
cd C:\work
git clone https://github.com/gemshub/xgem2mt
git clone -b develop_xgem2mt https://github.com/gemshub/xgems-jupyter
```

- `xgem2mt\demos` has the GEMS3K data and the Python demo scripts for all modes.
- `xgems-jupyter\how-to-use-xgem2mt-examples` has the notebooks.

Without Git, use "Code > Download ZIP" on both GitHub pages and unpack the folders.

## 4. Run the first model

The script `demos/quick_start.py` is the Quick start of the [README](../README.md#quick-start-python).
Run it from the `xgem2mt` folder, because it reads `demos/TestAD/CalcColumn-dat.lst` by a relative path:

```bat
cd C:\work\xgem2mt
python demos\quick_start.py
```

It prints the model setup (`print(task)`), then runs the model. This takes about 5 seconds.
If the last line of the output starts with `Mode A:`, the installation works. Files `xGEM2MT.log` and `ipmlog.txt` appear in the current folder.
They list GEM calculations that failed or needed a second try. Delete them at any time.

## 5. Run the notebooks

```bat
cd C:\work\xgems-jupyter\how-to-use-xgem2mt-examples
jupyter lab
```

Open one of the folders:

| Folder | What it shows |
|---|---|
| `calcite-dolomite-1d-advection/calcite-dolomite-1d-advection.ipynb` | 'A' mode, the model set up from the GEMS3K files and the node functions |
| `calcite-dolomite-1d-advection/calcite-dolomite-1d-advection-material.ipynb` | the same, with a `Material` for each node group |
| `calcite-dolomite-1d-advection/calcite-dolomite-1d-advection-material-parts.ipynb` | the model from a column experiment: pore solution and solid parts, scaled by the program |
| `calcite-dolomite-column/calcite-dolomite-column.ipynb` | 'F' mode, flow-through reactors |

In each notebook choose the kernel **Python 3 (ipykernel)** and use "Run > Run All Cells".
Start with `calcite-dolomite-1d-advection-material.ipynb`.

## 6. Make your own model

1. In [GEM-Selektor](https://gems.web.psi.ch), open your project and export **one** chemical system as
   GEMS3K files (`...-dat.lst`, `...-dch.dat`, `...-ipm.dat`, `...-dbr-0-0000.dat`).
   Put them in a new folder.
2. Copy `demos/quick_start.py` to your own folder. Change the path of `system`, the column data and the
   compositions.
3. Call `print(task)` before `runMassTransport()`. It shows what the program will do, the time step, the
   role and the values of each node, and the function that changes each value.
4. All functions, parameters and units: [python-mass-transport.md](python-mass-transport.md).

## Problems

| Message or symptom | Cause and action |
|---|---|
| `conda` is not found | You are not in the Miniforge Prompt. Open it from the Start menu. |
| `PackagesNotFoundError` | Check the channel URLs for typing errors. Open one in a browser: it must show a folder list, not an error. |
| `ImportError: DLL load failed` | Run `conda activate xgem2mt` first. Install the Microsoft Visual C++ Redistributable (<https://aka.ms/vs/17/release/vc_redist.x64.exe>) and try again. |
| Antivirus stops `python.exe` or removes a `.dll` | Add the Miniforge folder to the exclusions. |
| `import xgems` and `import xgem2mt` abort with "logger ... already exists" | The installed GEMS3K is too old. Run `conda update -c https://gemshub.github.io/xgem2mt/prerelease -c https://gemshub.github.io/xgems/prerelease -c conda-forge xgem2mt xgems`. |
| "DBR_DAT fileopen error" | A path is wrong, or a file of the `-dat.lst` is missing from the folder. The details are in `xGEM2MT.log`. |
| The notebook cannot find `demos/...` or the `CalciteColumn` folder | Start Jupyter from the folder of the notebook, or give the full path. |
| Many error lines in `ipmlog.txt` but `0 failed` after the run | The GEM solver tried again and then succeeded. This is normal. Only the message after the run counts. |

Report problems to the maintainer with: the output of `conda list xgem2mt gems3k xgems`, the full error text,
and the `xGEM2MT.log` file.

## Build from source (not needed to use the program)

The CI workflow builds the Windows packages with Conda and Ninja
(`packaging/conda/xgem2mt/bld.bat`, maintainer notes in [../packaging/README.md](../packaging/README.md)).
To change the C++ code you need the Visual Studio 2022 Build Tools ("Desktop development with C++"),
the Conda packages `gems3k` and `optima` from the xGEM2MT channel, `eigen`, `spdlog`, `nlohmann_json`,
`pybind11`, `numpy`, `cmake` and `ninja`, and the CMake settings of `bld.bat`.
