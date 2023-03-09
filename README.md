
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

For more information on APIs and more examples, take a look at the `/demo` directory.

### How to download xGEM2MT source code

* In your home directory, make a folder named e.g. ~/git/xGEM2MT.
* Change into ~/git/xGEM2MT and clone this repository from https://bitbucket.org/gems4/xgem2mt.git using a preinstalled free git client, e.g. SourceTree. 
* Alternatively on Mac OS X or linux, open a terminal and type in the command line:

```sh
cd ~/git/xGEM2MT
git clone https://bitbucket.org/gems4/xgem2mt.git .
```

### How to build xGEM2MT library and examples

* Install Dependencies by executing in ```~/xgem2mt$```

```sh
#!bash
sudo ./install-dependencies.sh
```

* To build xGEM2MT and install it in your home directory or in the system directory (as in the example below), a typical sequence of commands can be executed in the terminal:

```sh
cd ~/xgem2mt
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
sudo make install
```

* To execute the c++ demo:

```sh
cd ~/build/bin
./demo_gem2mt
```

