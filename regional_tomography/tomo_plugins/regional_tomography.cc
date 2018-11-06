/*
  Copyright (C) 2016 - 2018 by the authors of the ASPECT code.

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


#include "regional_tomography.h"
#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {

   template <int dim>
   RegionalTomography<dim>::RegionalTomography()
   :
   surface_boundary_id(1)
   {}


   template <int dim>
   void
   RegionalTomography<dim>::initialize ()
   {
    Utilities::AsciiDataInitial<dim>::initialize(2);
    // Find the boundary indicator that represents the surface
    surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");

    std::set<types::boundary_id> surface_boundary_set;
    surface_boundary_set.insert(surface_boundary_id);

    // The input ascii table contains two components, the crust depth and the LAB depth
    ascii_data_lab.initialize(surface_boundary_set,2);
   }


   template <int dim>
   double
   RegionalTomography<dim>::
   ascii_grid_vs (const Point<dim> &position) const
   {
    const double vs_perturbation = Utilities::AsciiDataInitial<dim>::get_data_component(position,1);
    return (vs_perturbation)*0.01;
   }


    template <int dim>
    double
    RegionalTomography<dim>::initial_temperature (const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);

      double vs_perturbation;
      if (depth <= max_grid_depth - smoothing_length_scale)
      {
    	  vs_perturbation = ascii_grid_vs(position);
      }
      //add smoothing between the two models
      else if (depth > max_grid_depth - smoothing_length_scale && depth < max_grid_depth)
      {
    	   const double scale_factor = (depth-(max_grid_depth-smoothing_length_scale))/smoothing_length_scale;
    	   vs_perturbation = 0.0 *(scale_factor) + ascii_grid_vs(position)*(1.0-scale_factor);
      }
      else
      {
    	   vs_perturbation = 0.0;
      }
      const double density_perturbation = vs_to_density * vs_perturbation;

      const double isotherm_depth       =  ascii_data_lab.get_data_component(surface_boundary_id, position, 0);

      double temperature_perturbation;
      if (depth < max_grid_depth && depth > isotherm_depth)
      {
    	  // scale the density perturbation into a temperature perturbation
    	  temperature_perturbation =  -1./thermal_alpha* density_perturbation;
      }
      else
    	  // set heterogeneity to zero down to a specified depth
    	  temperature_perturbation = 0.0;



      double background_temperature = 0.0;

     if (depth > isotherm_depth)
            return  1673.0 + (depth - isotherm_depth) * 0.0005 + temperature_perturbation;
     else
    	    return 273.15 + (depth/isotherm_depth) * (1673.15 - 273.15);
    }

    template <int dim>
    void
    RegionalTomography<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Regional tomography");
        {
            prm.declare_entry ("Maximum grid depth", "350000.0",
          	                 Patterns::Double (0),
          	                 "The maximum depth of the Vs ascii grid. The model will read in  "
          	                 "Vs from golbal model or use temperature from adiabatic boundary below this depth.");
            prm.declare_entry ("Thermal expansion coefficient", "2e-5",
          	                 Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: $1/K$.");
            prm.declare_entry ("Smoothing length scale", "10000.0",
          	                 Patterns::Double (0),
          	                 "The depth range (above maximum grid depth) over which to smooth. "
          	                 "The boundary is smoothed using a depth weighted combination of Vs.");
            prm.declare_entry ("Vs to density", "0.25",
                      	      Patterns::Double (0),
                      	      "Vs to density scaling factor");
            Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
                                                              "regional_tomography_3d.txt");
        }
      prm.leave_subsection();
      }
      prm.leave_subsection();

      Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
         	    	     	    	         	        "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
         	    	     	    	         	        "litho.kenya.txt");
    }

    template <int dim>
    void
    RegionalTomography<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Regional tomography");
        {
          max_grid_depth           = prm.get_double ("Maximum grid depth");
          smoothing_length_scale   = prm.get_double ("Smoothing length scale");
          thermal_alpha            = prm.get_double ("Thermal expansion coefficient");
          vs_to_density            = prm.get_double ("Vs to density");

          Utilities::AsciiDataBase<dim>::parse_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      ascii_data_lab.initialize_simulator (this->get_simulator());

      // Note: parse_parameters will call initialize for us
      ascii_data_lab.parse_parameters(prm);

    }

  }

}


namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(RegionalTomography,
                                              "regional tomography",
                                              "An initial temperature condition that allows for discretizing "
                                              "is designed specifically for the ellipsoidal chunk geometry model.")
  }
}
