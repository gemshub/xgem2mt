# xGEMS2MT is a C++ and Python library for mass transport modeling with Gibbs energy minimization
#
# Copyright (C) 2022 S.Dmytriieva, Dmitrii Kulik
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from xgem2mt import *
from pathlib import Path

def main():
    try:
        mttask = MassTransport("S", 11)
        gem2mt_json = Path("gem2mt-data/CalcColumnAD.json").read_text().replace('\n', ' ')
        dch_json = Path("gem2mt-data/CalcColumn-dch.json").read_text().replace('\n', ' ')
        ipm_json = Path("gem2mt-data/CalcColumn-ipm.json").read_text().replace('\n', ' ')
        dbr_json1 = Path("gem2mt-data/CalcColumn-dbr-0-0000.json").read_text().replace('\n', ' ')
        dbr_json0 = Path("gem2mt-data/CalcColumn-dbr-0-0001.json").read_text().replace('\n', ' ')
        dbr_jsons = [dbr_json1, dbr_json0]

        mttask.readTaskString(gem2mt_json)
        mttask.initNodesFromStrings(dch_json, ipm_json, dbr_jsons)
        print(mttask)  # setup summary: nodes, time step, properties, phases

        mttask.runMassTransport()
        print(mttask.dump())  # the task as JSON

    except RuntimeError as ex:
        print(ex)

if __name__ == "__main__":
    main()

