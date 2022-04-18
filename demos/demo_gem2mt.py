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

def main():
    try:
        mttask = MassTransport()
        gem2mt_init = "gem2mt-data/CalcColumnAD.json"
        gems_init_lst = "gem2mt-data/CalcColumn-dat.lst"
        gems_init_dbr_lst = "gem2mt-data/CalcColumn-dbr.lst"

        mttask.initializeGEM2MT(gem2mt_init)
        mttask.initializeNodes(gems_init_lst, gems_init_dbr_lst)

        mttask.runGEM2MT()
        print(mttask)

    except RuntimeError as ex:
        print(ex)

if __name__ == "__main__":
    main()

