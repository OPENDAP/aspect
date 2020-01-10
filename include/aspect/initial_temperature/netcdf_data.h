//
// Created by Kodi Neumiller on 12/30/19.
//
/*
   Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

   This file is part of ASPECT.

   ASPECT is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   ASPECT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with ASPECT; see the file LICENSE.  If not see
   <http://www.gnu.org/licenses/>.
 */

#ifndef ASPECT_OPENDAP_NETCDF_DATA_H
#define ASPECT_OPENDAP_NETCDF_DATA_H

//#include <aspect/initial_temperature/interface.h>
//#include <aspect/simulator_access.h>

namespace aspect
{
    namespace InitialTemperature
    {
        using namespace std;
        /**
         * Used to store variable information inside an object.
         * The netcdf_data object can be used to look up values for variables in a netcdf file.
         * This will allow us to store variable names, id's, etc inside a prm file and then call the netcdf_data
         *  object when we need to lookup a specfic variable's data by reference.
         */
         class NetcdfData
         {
         public:
             /**
              * Constructor
              */
             NetcdfData() {}

             /**
              * Destructor
              */
             virtual ~NetcdfData() {}

             //** Getters **//
             string getVar1() {
                 return var1;
             }
             string getVar2() {
                 return var2;
             }
             string getVar3() {
                 return var3;
             }
             string getValuesVar() {
                 return valuesVar;
             }

             /** Setters **/
             void setVar1(string newVal) {
                 var1 = newVal;
             }
             void setVar2(string newVal) {
                 var2 = newVal;
             }
             void setVar3(string newVal) {
                 var3 = newVal;
             }
             void setValuesVar(string newVal) {
                 valuesVar = newVal;
             }

         private:
             string var1;
             string var2;
             string var3;
             string valuesVar;
         };
    }
}
#endif //ASPECT_OPENDAP_NETCDF_DATA_H
