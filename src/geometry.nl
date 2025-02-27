/* file: geometry.nl
 * purpose: namelist for mesh setup
 *
 * Michael Borland, 1992
 */
#include "namelist.h"

#namelist define_geometry
    long nz = 0;
    long nr = 0;
    double zmin = 0;
    double zmax = 0;
    double rmax = 0;
    double zr_factor = 1;              /* z and r coordinates in boundary file are multiplied by this factor */
    STRING rootname = NULL;            /* rootname */
    STRING boundary = NULL;            /* SUPERFISH-format data defining the boundary           */
    STRING boundary_output = NULL;     /* mpl-format output showing the boundary                */
    STRING urmel_boundary_output = NULL; 
    STRING discrete_boundary_output = NULL;
    STRING interior_points = NULL;     /* mpl-format output showing the interior points         */
    STRING lower = "Dirichlet";        /* boundary conditions:                                  */
    STRING upper = "Neumann";          /* Dirichlet == Electric field parallel to boundary      */
    STRING right = "Neumann";          /* Neumann   == Electric field perpendicular to boundary */
    STRING left  = "Neumann";          /* defaults are the same as the SUPERFISH defaults       */
    long include_TE_fields = 0;        /* used to turn on TE-mode fields */
    long exclude_TM_fields = 0;        /* used to turn off TM-mode fields */
    long turn_off_Er = 0;
    long turn_off_Ez = 0;
    long turn_off_Ephi = 0;
    long turn_off_Br = 0;
    long turn_off_Bz = 0;
    long turn_off_Bphi = 0;
    long print_grids = 0;              /* logical--if non-zero, program prints grid picture to stdout */
    long radial_interpolation = 1;
    long longitudinal_interpolation = 1;
    long radial_smearing = 0;
    long longitudinal_smearing = 0;
#end

#namelist point
    int nt;
    double x;
    double y;
    double x0;
    double y0;
    double r;
    double theta;
    double a = 0;
    double b = 0;
    double aovrb = 1;
    double potential = 0;
    short material_id = 0;
    short ramp_potential = 0;
    short change_direction = 0;
#end

