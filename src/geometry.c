/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: geometry.c
 * purpose: process geometry definition for EM problem 
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"

#define N_BOUNDARY_TYPES 2
static char *boundary_type[N_BOUNDARY_TYPES] = {
  "Dirichlet", "Neumann"};

void process_geometry_definition(
                                 FIELDS *EM_problem,
                                 NAMELIST_TEXT *nl_text)
{
#include "geometry.h"
  FILE *fp_pt;
  char s[1024];
  long n_points, k, iz, ir;
  POINTSP *pnt = NULL;
  double xi, yi, xf, yf;
  FILE *fpUrmel = NULL;

  nz = nr = print_grids = include_TE_fields = exclude_TM_fields = 0;
  zmin = zmax = rmax = 0;
  boundary = boundary_output = interior_points = rootname =
    urmel_boundary_output = NULL;
  radial_interpolation = longitudinal_interpolation = 1;
  radial_smearing = longitudinal_smearing = 0;
  turn_off_Er = turn_off_Ez = turn_off_Ephi = 0;
  turn_off_Br = turn_off_Bz = turn_off_Bphi = 0;

  process_namelist(&define_geometry, nl_text);
  print_namelist(stdout, &define_geometry);
  if (rootname)
    {
      free(EM_problem->rootname);
      EM_problem->rootname = rootname;
    }
  if (!boundary)
    bomb("no boundary-points file named", NULL);
  if (nz < 2 || nr < 2)
    bomb("nz and nr must be greater than 2", NULL);
  if (zmax <= zmin)
    bomb("zmax must be greater than zmin", NULL);
  if (rmax <= 0)
    bomb("rmax must be greater than 0", NULL);
  if ((EM_problem->lower_bc = match_string(lower, boundary_type, N_BOUNDARY_TYPES, 0)) < 0)
    bomb("unrecognized boundary condition for lower boundary", NULL);
  if (EM_problem->lower_bc == NEUMANN)
    bomb("Neumann boundary conditions are meaningless on the lower boundary", NULL);
  if ((EM_problem->upper_bc = match_string(upper, boundary_type, N_BOUNDARY_TYPES, 0)) < 0)
    bomb("unrecognized boundary condition for upper boundary", NULL);
  if (EM_problem->upper_bc == DIRICHLET)
    bomb("Dirichlet boundary conditions are not implemented for the upper boundary", NULL);
  if ((EM_problem->right_bc = match_string(right, boundary_type, N_BOUNDARY_TYPES, 0)) < 0)
    bomb("unrecognized boundary condition for right boundary", NULL);
  if ((EM_problem->left_bc = match_string(left, boundary_type, N_BOUNDARY_TYPES, 0)) < 0)
    bomb("unrecognized boundary condition for left boundary", NULL);

  EM_problem->nz = nz;
  EM_problem->nr = nr;
  EM_problem->zmin = zmin;
  EM_problem->zmax = zmax;
  EM_problem->rmax = rmax;
  EM_problem->dz = (zmax - zmin) / (nz - 1);
  EM_problem->dr = rmax / (nr - 1);
  EM_problem->zt = 0;
  EM_problem->metal_flags = (short **)zarray_2d(sizeof(**EM_problem->metal_flags), nz + 1, nr + 1);
  EM_problem->modes = (radial_interpolation ? FL_RADIAL_INTERP : 0) + (longitudinal_interpolation ? FL_LONGIT_INTERP : 0) + (radial_smearing ? FL_RADIAL_SMEAR : 0) + (longitudinal_smearing ? FL_LONGIT_SMEAR : 0);
  EM_problem->area_c = tmalloc(sizeof(*EM_problem->area_c) * (nr + 1));
  EM_problem->area_o = tmalloc(sizeof(*EM_problem->area_o) * (nr + 1));
  EM_problem->PsiImposed = (double **)zarray_2d(sizeof(**EM_problem->PsiImposed), nz + 1, nr + 1);
  EM_problem->EzImposed = (double **)zarray_2d(sizeof(**EM_problem->EzImposed), nz + 1, nr + 1);
  EM_problem->ErImposed = (double **)zarray_2d(sizeof(**EM_problem->ErImposed), nz + 1, nr + 1);
  EM_problem->materialID = (short **)zarray_2d(sizeof(**EM_problem->materialID), nz + 1, nr + 1);
  EM_problem->Q = (double **)zarray_2d(sizeof(**EM_problem->Q), nz + 1, nr + 1);
  EM_problem->Jz = (double **)zarray_2d(sizeof(**EM_problem->Jz), nz + 1, nr + 1);
  if (!exclude_TM_fields)
    {
      EM_problem->Ez = (double **)zarray_2d(sizeof(**EM_problem->Ez), nz + 1, nr + 1);
      EM_problem->Er = (double **)zarray_2d(sizeof(**EM_problem->Er), nz + 1, nr + 1);
      EM_problem->Jr = (double **)zarray_2d(sizeof(**EM_problem->Jr), nz + 1, nr + 1);
      EM_problem->Bphi = (double **)zarray_2d(sizeof(**EM_problem->Bphi), nz + 1, nr + 1);
      EM_problem->Bphi2 = (double **)zarray_2d(sizeof(**EM_problem->Bphi2), nz + 1, nr + 1);
      EM_problem->Psi = (double **)zarray_2d(sizeof(**EM_problem->Psi), nz + 1, nr + 1);
      EM_problem->modes |= FL_TM_FIELDS;
    }
  if (include_TE_fields)
    {
      EM_problem->Ephi = (double **)zarray_2d(sizeof(**EM_problem->Ephi), nz + 1, nr + 1);
      EM_problem->Jphi = (double **)zarray_2d(sizeof(**EM_problem->Jphi), nz + 1, nr + 1);
      EM_problem->Bz = (double **)zarray_2d(sizeof(**EM_problem->Bz), nz + 1, nr + 1);
      EM_problem->Br = (double **)zarray_2d(sizeof(**EM_problem->Br), nz + 1, nr + 1);
      EM_problem->Bz2 = (double **)zarray_2d(sizeof(**EM_problem->Bz2), nz + 1, nr + 1);
      EM_problem->Br2 = (double **)zarray_2d(sizeof(**EM_problem->Br2), nz + 1, nr + 1);
      EM_problem->modes |= FL_TE_FIELDS;
    }
  if (!(EM_problem->modes & FL_TM_FIELDS) && !(EM_problem->modes & FL_TE_FIELDS))
    fputs("warning: all rf fields are disabled\n", stderr);
  EM_problem->turnOffEr = turn_off_Er;
  EM_problem->turnOffEz = turn_off_Ez;
  EM_problem->turnOffEphi = turn_off_Ephi;
  EM_problem->turnOffBr = turn_off_Br;
  EM_problem->turnOffBz = turn_off_Bz;
  EM_problem->turnOffBphi = turn_off_Bphi;

  EM_problem->constantEz = EM_problem->constantBz = 0;
  EM_problem->constantEr = EM_problem->constantBr = 0;
  EM_problem->constantEphi = EM_problem->constantBphi = 0;

  EM_problem->time_step = 0;

  if (urmel_boundary_output)
    {
      char *filename;
      filename = compose_filename(urmel_boundary_output, EM_problem->rootname);
      fpUrmel = fopen_e(filename, "w", 0);
      if (filename != urmel_boundary_output)
        free(filename);
    }

  /* process the POINTSP namelists to fill in the metal array */
  fp_pt = fopen_e(boundary, "r", 0);
  n_points = 0;
  while (get_namelist(s, 1024, fp_pt))
    {
      str_tolower(s);
      if (!scan_namelist(&namelist_text, s))
        bomb("namelist scanning failure", NULL);
      if (strcmp(namelist_text.group_name, "point") != 0 &&
          strcmp(namelist_text.group_name, "po") != 0)
        bomb("unrecognized namelist in boundary defintion", NULL);
      nt = 1;
      change_direction = 0;
      x0 = y0 = r = potential = material_id = ramp_potential = a = b = 0;
      x = y = theta = DBL_MAX;
      process_namelist(&point, &namelist_text);
      if (theta != DBL_MAX)
        theta *= PI / 180.0;
      free_namelist_text(&namelist_text);

      if (n_points == 0 && nt != 1)
        bomb("first point in geometry definition must have nt=1", NULL);

      /* store in array */
      pnt = trealloc(pnt, sizeof(*pnt) * (n_points + 1));
      pnt[n_points].type = nt;
      pnt[n_points].x = x * zr_factor;
      pnt[n_points].y = y * zr_factor;
      pnt[n_points].x0 = x0 * zr_factor;
      pnt[n_points].y0 = y0 * zr_factor;
      pnt[n_points].r = r * zr_factor;
      a *= zr_factor;
      b *= zr_factor;
      if (x != DBL_MAX)
        x *= zr_factor;
      if (y != DBL_MAX)
        y *= zr_factor;
      if (nt == 2)
        {
          if (theta == DBL_MAX && (x == DBL_MAX || y == DBL_MAX))
            bomb("If NT=2, must give theta or (x, y)", NULL);
          if (theta != DBL_MAX && (x != DBL_MAX || y != DBL_MAX))
            bomb("If NT=2, must give theta or (x, y), not both", NULL);
          if (r == 0)
            {
              /* ellipse */
              if (a == 0 && b == 0)
                bomb("Can't have NT=2, R=0, A=0, and B=0\n", NULL);
              if (a == 0 || b == 0)
                {
                  if (aovrb <= 0)
                    bomb("Can't have NT=2, R=0, AOVRB<=0 and A=0 or B=0\n", NULL);
                  if (a == 0)
                    a = b * aovrb;
                  else
                    b = a / aovrb;
                }
              if (theta == DBL_MAX)
                theta = atan2(y / b, x / a);
            }
          else
            {
              /* circle */
              if (theta == DBL_MAX)
                theta = atan2(y, x);
            }
        }
      if (change_direction)
        theta = PIx2 + theta;
      pnt[n_points].a = a;
      pnt[n_points].b = b;
      pnt[n_points].potential = potential;
      pnt[n_points].material = material_id;
      pnt[n_points].ramp_potential = ramp_potential;
      pnt[n_points].theta = theta;
      n_points++;
    }
  if (n_points == 0)
    bomb("too few boundary points in geometry file", NULL);

  /* set up to "draw" the boundary on the metal_flags array */
  set_up_drawing(nz, nr, zmin, zmax, 0.0L, rmax, 0.25, EM_problem->metal_flags, EM_problem->PsiImposed,
                 EM_problem->materialID);

  /* mark boundary points */
  xi = pnt[0].x + pnt[0].x0;
  yi = pnt[0].y + pnt[0].y0;
  if (fpUrmel)
    fprintf(fpUrmel, "%e, %e\n", yi, xi);
  starting_point(xi, yi, pnt[0].potential, pnt[0].material);
  for (k = 1; k < n_points; k++)
    {
      switch (pnt[k].type)
        {
        case 1: /* line */
          xf = pnt[k].x + pnt[k].x0;
          yf = pnt[k].y + pnt[k].y0;
          printf("line from (%e, %e) to (%e, %e)\n", xi, yi, xf, yf);
          if (fpUrmel)
            fprintf(fpUrmel, "%e, %e\n", yf, xf);
          draw_line(xi, yi, xf, yf, pnt[k].potential, pnt[k].material, pnt[k].ramp_potential ? pnt[k - 1].potential : pnt[k].potential);
          if (out_of_bounds())
            {
              printf("error: points on the boundary that were outside the simulation region for line:\n");
              printf("  xi = %e, yi=%e, xf=%e, yf=%e\n",
                     xi, yi, xf, yf);
              exit(1);
            }
          xi = xf;
          yi = yf;
          break;
        case 2: /* circle or ellipse */
          if (pnt[k].r != 0)
            {
              printf("circle of radius %e centered on (%e, %e),\n\tstarting at (%e, %e) going to angle %e\n",
                     pnt[k].r, pnt[k].x0, pnt[k].y0, xi, yi, pnt[k].theta);
              draw_circle(xi, yi, pnt[k].x0, pnt[k].y0,
                          pnt[k].r, pnt[k].theta, &xf, &yf, pnt[k].potential, pnt[k].material);
              if (fpUrmel)
                {
                  fprintf(fpUrmel, "-1, %e\n", pnt[k].r);
                  fprintf(fpUrmel, "%e, %e\n", yf, xf);
                }
            }
          else
            {
              printf("ellipse with semi-axes of (%e, %e) centered on (%e, %e),\n\tstarting at (%e, %e) going to angle %e\n",
                     pnt[k].a, pnt[k].b, pnt[k].x0, pnt[k].y0, xi, yi, pnt[k].theta);
              draw_ellipse(xi, yi, pnt[k].x0, pnt[k].y0,
                           pnt[k].a, pnt[k].b, pnt[k].theta, &xf, &yf, pnt[k].potential, pnt[k].material);
            }
          printf("  xf = %e, yf=%e\n", xf, yf);
          xi = xf;
          yi = yf;
          if (out_of_bounds())
            {
              printf("error: points on the boundary that were outside the simulation region for circle:/ellipse\n");
              printf("  x0 = %e, y0=%e, r=%e, a=%e, b=%e, theta=%e\n",
                     pnt[k].x0, pnt[k].y0, pnt[k].r, pnt[k].a, pnt[k].b, pnt[k].theta * 180 / PI);
              exit(1);
            }
          break;
        case 3: /* beginning of new shape */
          xi = pnt[k].x + pnt[k].x0;
          yi = pnt[k].y + pnt[k].y0;
          printf("new shape beginning with (%e, %e)\n", xi, yi);
          add_breakpoint(xi, yi, pnt[k].potential, pnt[k].material);
          break;
        case 4: /* definition of in-vacuum point */
          printf("marking (%e, %e) as interior point\n", pnt[k].x + pnt[k].x0, pnt[k].y + pnt[k].y0);
          set_vacuum_point(pnt[k].x + pnt[k].x0, pnt[k].y + pnt[k].y0);
          break;
        default: /* unknown/unimplemented */
          printf("warning: unknown/unimplemented point type code: nt=%d\n",
                 pnt[k].type);
          break;
        }
    }

  if (fpUrmel)
    {
      fclose(fpUrmel);
      fpUrmel = NULL;
    }

  if (boundary_output)
    dump_boundary_points(boundary_output, EM_problem->rootname);
  if (print_grids)
    type_flag_drawing(stdout, "boundary marked", EM_problem->metal_flags, nz, nr);

  /* Check all points (i*dx, j*dy) and find those that are outside
   * of the boundary.
   */
  fill_in_drawing(EM_problem->metal_flags, EM_problem->PsiImposed);
  if (print_grids)
    type_flag_drawing(stdout, "filled in", EM_problem->metal_flags, nz, nr);
  trim_drawing(EM_problem->metal_flags);
  if (print_grids)
    type_flag_drawing(stdout, "trimmed", EM_problem->metal_flags, nz, nr);

  if (discrete_boundary_output)
    dump_dboundary_points(discrete_boundary_output, EM_problem->rootname);
  if (interior_points)
    dump_interior_points(interior_points, EM_problem->rootname);

  /* Set flags to indicate effects of metal on fields
   */
  for (iz = 0; iz < nz + 1; iz++)
    {
      for (ir = 0; ir < nr + 1; ir++)
        {
          if (!(EM_problem->metal_flags[iz][ir] & FL_IS_METAL))
            continue;
          /* a ring of metal exists at this grid position, so Ephi must be zero */
          EM_problem->metal_flags[iz][ir] |= FL_EPHI_ZERO;
          if (iz && EM_problem->metal_flags[iz - 1][ir] & FL_IS_METAL)
            {
              /* metal exists from z=zmin+(iz-1)*dz to zmin+iz*dz at r=ir*dr, so Ez and Br
               * must be zero in this range, i.e., at the grid point at z=zmin+(iz-1/2)*dz, r=ir*dr.
               */
              EM_problem->metal_flags[iz][ir] |= FL_EZ_ZERO;
              EM_problem->metal_flags[iz][ir] |= FL_BR_ZERO;
            }
          if (ir == 0 || (ir && EM_problem->metal_flags[iz][ir - 1] & FL_IS_METAL))
            {
              /* metal exists from r=(ir-1)*dr to r=ir*dr at z=zmin+iz*dz, so Er and Bz 
               * must be zero in this range, i.e., at the grid point at r=(ir-1/2)*dr, z=zmin+iz*dz
               */
              EM_problem->metal_flags[iz][ir] |= FL_ER_ZERO;
              EM_problem->metal_flags[iz][ir] |= FL_BZ_ZERO;
            }
          if (iz && EM_problem->metal_flags[iz - 1][ir] & FL_IS_METAL &&
              (ir == 0 || (EM_problem->metal_flags[iz][ir - 1] & FL_IS_METAL && EM_problem->metal_flags[iz - 1][ir - 1] & FL_IS_METAL)))
            /* metal exists around the square [(ir-1)*dr, ir*dr] x [zmin+(iz-1)*dz, zmin+iz*dz], so Bphi must be
             * zero for this area, i.e., at the grid point at r=(ir-1/2)*dr, z=zmin+(iz-1/2)*dz
             */
            EM_problem->metal_flags[iz][ir] |= FL_BPHI_ZERO;
        }
    }

  /* set flags to indicate last active point in each row */
  for (ir = 0; ir < nr; ir++)
    {
      iz = nz - 1;
      while (iz >= 0 && EM_problem->metal_flags[iz][ir] & FL_IS_INSIDE)
        iz--;
      if (iz != 0 && (iz += 1) < nz - 1)
        {
          EM_problem->metal_flags[iz][ir] |= FL_LAST_ACTIVE_PT;
        }
    }
  /*
    if (print_grids)
    type_flag_drawing(stdout, "boundary conditions marked", EM_problem->metal_flags, nz, nr);
  */
  if (pnt)
    free(pnt);
}
