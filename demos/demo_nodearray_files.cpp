// xGEMS2MT is a C++ and Python library for mass transport modeling with Gibbs energy minimization
//
// Copyright (C) 2022 S.Dmytriieva, Dmitrii Kulik
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License version 2.1 
// as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

// xGEMS includes
#include <xGEM2MT/ChemicalEngineArray.hpp>
using namespace xGEM2MT;


int main(int argc, char **argv)
{
    // Testing= with input from IPM files
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
        //engine.equilibrate(0);
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
