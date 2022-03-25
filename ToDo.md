### Implemented

1. ChemicalEngineArray - created an xGEM2MT API that extends the xGEMS API to cover the TNodeArray class API.  

> All get/set functions works with  TNodeArray::na->pNodT1(), only dump use TNodeArray::na->pNodT0(). To put changes to TNodeArray::na->pNodT0() use function nextTimePoint(true);

2. demo_nodearray_files.cpp and demo_nodearray_strings.cpp are instructions how work with ChemicalEngineArray.

3. MassTransport is a wrapper class for TGEM2MT code.

4. demo_gem2mt.cpp is instructions how work with MassTransport class.

5. GEMS-198 Create Python API of xGEM2MT using pybind11 header

### To do

1. GEMS-198 Make demo for Python API of xGEM2MT and test API

2. Research and add to TGEM2MT code into GEMS3K and to MassTransport class functions for accsees to internal nodes and values. How user can set up and monitoring of mass transport process. Need suggestions with Dan and Dmitrii.

3. Implement test ?

### Problems&questions

1.  Into environment.devenv.yml we ask ```- gems3k=4.1.1```, but need installation of gems3k with GEM2MT implementation ```cmake .. -DBUILD_GEM2MT=ON```
Need suggestion about using install with nodearray-gem. 

