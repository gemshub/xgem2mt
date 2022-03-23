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

// xGEMS includes
#include <xGEM2MT/MassTransport.hpp>
using namespace xGEM2MT;


int main(int argc, char **argv)
{
    // Testing= with input from IPM files
    MassTransport mttask;
    std::string gem2mt_init = "gem2mt-data/CalcColumnAD.json";
    std::string gems_init_lst = "gem2mt-data/CalcColumn-dat.lst";
    std::string gems_init_dbr_lst = "gem2mt-data/CalcColumn-dbr.lst";

    try
    {
        mttask.initializeGEM2MT(gem2mt_init);
        mttask.initializeNodes(gems_init_lst, gems_init_dbr_lst);

        //mttask.runGEM2MT();

        std::cout << "Readed: " << mttask;
        return 0;
    }
    catch (std::runtime_error& xcpt)
    {
        std::cout << xcpt.what() << std::endl;
    }
    return 1;
}
