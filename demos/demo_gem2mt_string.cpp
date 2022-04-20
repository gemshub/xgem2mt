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
#include <xGEM2MT/MassTransport.hpp>
using namespace xGEM2MT;

// Reading the whole file (given by path) into a string str
// and skipping the header e.g. {"dch": and the very last } at the end of file
int getJsonFile( const string path, string& str )
{
    ifstream f(path); //taking file as inputstream
    // string str;
    if(f) {
        ostringstream ss;
        ss << f.rdbuf(); // reading data until-from+1);
        str = ss.str();
        return 0;
    }
    cout << "Error reading file " << path << " into text string";
    return 1;
}

int main(int argc, char **argv)
{
    MassTransport mttask;
    std::string gem2mt_json;
    std::string dch_json;
    std::string ipm_json;
    std::string dbr_json;
    std::vector<std::string> dbr_jsons;

    try
    {
        getJsonFile( "gem2mt-data/CalcColumnAD.json", gem2mt_json );
        getJsonFile( "gem2mt-data/CalcColumn-dch.json", dch_json );
        getJsonFile( "gem2mt-data/CalcColumn-ipm.json", ipm_json );
        getJsonFile( "gem2mt-data/CalcColumn-dbr-0-0000.json", dbr_json );
        dbr_jsons.push_back(dbr_json);
        getJsonFile( "gem2mt-data/CalcColumn-dbr-0-0001.json", dbr_json );
        dbr_jsons.push_back(dbr_json);
        mttask.initializeGEM2MTString(gem2mt_json);
        mttask.initializeNodesFromJsonStrings(dch_json, ipm_json, dbr_jsons);

        mttask.runGEM2MT();
        std::cout << "Readed: " << mttask;
        return 0;
    }
    catch (std::runtime_error& xcpt)
    {
        std::cout << xcpt.what() << std::endl;
    }
    return 1;
}
