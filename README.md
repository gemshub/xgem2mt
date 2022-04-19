
## xGEM2MT 
Linux, OSX, Windows

This is a numerical solver of chemical equilibria for thermodynamic modelling, extended (relative to GEMS3K code) with new C++ and Python APIs (compatible with [Reaktoro framework](http://reaktoro.org)), and supposed to replace GEMS3K in next-generation software.  


### Briefly about xGEM2MT

* Currently uses the [GEMS3K](https://bitbucket.org/gems4/gems3k) code that implements the improved GEM IPM-3 algorithm for Gibbs energy minization in very complex non-ideal chemical systems with two-sided metastability constraints.
* The xGEM2MT API extends the [xGEMS API](https://bitbucket.org/gems4/xgems/src/master/) to cover the TNodeArray class API of (standalone) GEM2MT as currently implemented in nodearray-gem example.
* Written in C/C++. The xGEM2MT code can be compiled as a standalone program (see /demos); or built as a static or dynamic library for coupling with a mass transport simulator or another code.

### License
* GEMS3K is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* GEMS3K is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
* You should have received a copy of the GNU Lesser General Public License along with GEMS3K code. If not, see http://www.gnu.org/licenses/. 

### Simple C++ API example

* Using the ChemicalEngineArray class API

```cpp
#include <xGEM2MT/ChemicalEngineArray.hpp>
using namespace xGEM2MT;


int main(int argc, char **argv)
{
    ChemicalEngineArray engine(5);
    std::vector<long int> nodeTypes={1,0,0,0,0};
    std::string gems_init_lst = "nodearray-data/CalcColumn-dat.lst";
    std::string gems_init_dbr_lst = "nodearray-data/CalcColumn-dbr.lst";

    try
    {
        engine.initialize(gems_init_lst, gems_init_dbr_lst, nodeTypes );

        Vector bj = engine.elementAmounts(0);
        std::cout << "Vector bj[0]: " << bj.transpose() << std::endl;
        bj = engine.elementAmounts(1);
        std::cout << "Vector bj[1]: " << bj.transpose() << std::endl;

        engine.equilibrate();
        engine.nextTimePoint(true); // copy results to pNodT0()

        engine.dumpGEMS3K("CalcColumn-dump-dat.lst");

        return 0;
    }
    catch (std::runtime_error& xcpt)
    {
        std::cout << xcpt.what() << std::endl;
    }
    return 1;
}
```

* Using  the Mass transport  API 

```cpp
#include <xGEM2MT/MassTransport.hpp>
using namespace xGEM2MT;

int main(int argc, char **argv)
{
    MassTransport mttask;
    std::string gem2mt_init = "gem2mt-data/CalcColumnAD.json";
    std::string gems_init_lst = "gem2mt-data/CalcColumn-dat.lst";
    std::string gems_init_dbr_lst = "gem2mt-data/CalcColumn-dbr.lst";

    try
    {
        mttask.initializeGEM2MT(gem2mt_init);
        mttask.initializeNodes(gems_init_lst, gems_init_dbr_lst);
        mttask.runGEM2MT();
        return 0;
    }
    catch (std::runtime_error& xcpt)
    {
        std::cout << xcpt.what() << std::endl;
    }
    return 1;
}
```

### Simple Python API example

* Using a json database file

```python3
from xgem2mt import *

def main():
    try:
        engine = ChemicalEngineArray(5)
        nodeTypes = [1,0,0,0,0]
        gems_init_lst = "nodearray-data/CalcColumn-dat.lst"
        gems_init_dbr_lst = "nodearray-data/CalcColumn-dbr.lst"

        engine.initialize(gems_init_lst, gems_init_dbr_lst, nodeTypes )
        bj = engine.elementAmounts(0)
        print(bj)

        engine.equilibrate();
        engine.nextTimePoint(True)
        print(engine)
    except RuntimeError as ex:
        print(ex)
```

For more information on APIs and more examples, take a look at the /demo directory.

### How to download xGEM2MT source code

* In your home directory, make a folder named e.g. ~/git/xGEM2MT.
* Change into ~/git/xGEM2MT and clone this repository from https://bitbucket.org/gems4/xgem2mt.git using a preinstalled free git client, e.g. SourceTree. 
* Alternatively on Mac OS X or linux, open a terminal and type in the command line:

```sh
cd ~/git/xGEM2MT
git clone https://bitbucket.org/gems4/xgem2mt.git .
```

### How to build xGEM2MT library and examples

xGEM2MT uses GEMS3K and Eigen3 as external projects. These are downloaded and compiled together with xGEM2MT during the make step. 
For compiling xGEM2MT python bindings you need a recent version of pybind11 installed (including development headers, 2.2.0 is the minimum version required) and Eigen library ([v3.3.9](https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz)). 

* For eigen3 download [v3.3.9](https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz) unpack the archive and install it as follows:

```sh
cd /eigen-3.3.9
mkdir build
cd build
cmake .. 
make install
```

* For pybind11. In case your distribution does not provide pybind11 development packages, or the version is too old, you can install it as follows:

```sh
cd ~/git/xGEM2MT/python
./install-dependencies.sh
```

* To install GEMS3K with the implementation of the TNodeArray and TGEM2MT classes execute script:

```sh
cd ~/git/xGEM2MT
./install-dependencies.sh
```

* To compile xGEM2MT with python bindings and demos:

`-DCMAKE_PREFIX_PATH=/your/local/codes/directory` is necessary if pybind11 is not installed in the default include path `/usr/include/` or `/usr/lib/cmake/`, provided that `/your/local/codes/directory/pybind11/pybind11Config.cmake` cmake config file is found at this path.

```sh
cd ~/git/xGEM2MT
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=/your/local/codes/directory
make 
```

* For building xGEM2MT for use from C++ codes only (i.e. no python in the system), replace the above cmake command with

```sh
cmake .. -DXGEMS_BUILD_PYTHON=OFF
```

* To install xGEM2MT into /usr/local as a library with /includes:

```sh
cd ~/git/xGEM2MT/build
sudo make install 
```

* To execute the c++ demo:

```sh
cd ~/git/xGEM2MT/build/bin
./demo_gem2mt
```

* To execute the python demo:

```sh
cd ~/git/xGEM2MT/demo
python3 demo_gem2mt.py
```

