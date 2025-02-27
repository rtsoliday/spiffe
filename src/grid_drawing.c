/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file    : grid_drawing.c
 * purpose : routines for drawing on a grid--original version written for program maskgen
 *           marks boundary points on ``lattice'' and makes list of
 *           boundary points, then finds all exterior points.
 *
 * contents: set_up_drawing(), draw_line(), draw_circle(),
 *           fill_in_drawing(), type_drawing(), trim_drawing(),
 *           add_point(), dump_boundary_points()
 *
 * Michael Borland, 1988, 1991, 1992
 */
#include "mdb.h"
#include "spiffe.h"

#undef SIGN
static double sign_tmp;
#define SIGN(x) ((sign_tmp = (x)) > 0 ? 1 : (sign_tmp < 0 ? -1 : 0))

#define ARRAY_SIZE_INCREMENT 20

static int n_points_max, i_point;
static double *x_bnd, *y_bnd;
static double *x_dbnd, *y_dbnd;
static long *ix_dbnd, *iy_dbnd;
static short *material_bnd;
static long *index_brkpt;
static long n_brkpts = 0;

static int nx, ny;
static double xmin, xmax, ymin, ymax;
static double x_tolerance, y_tolerance;
static double dtheta, dx, dy;
static short **lattice;
static double **potential;
static short **material;
static long n_outside = 0;

static long vacuum_point_defined = 0;
static int ix_vacuum, iy_vacuum;

long getAllDiscreteBoundaryPoints(long **iz, long **ir)
{
  long i, j, ip1, ip2, nPts, nNew;

  ip1 = 0;
  ip2 = i_point;

  if (ip2 != ip1)
    {
      /* Record the (iz, ir) coordinates of the discretized boundary, skipping duplicates */
      *iz = tmalloc(10 * sizeof(**iz) * (ip2 - ip1));
      *ir = tmalloc(10 * sizeof(**ir) * (ip2 - ip1));
      j = 0;
      for (i = ip1; i < ip2; i++)
        {
          if (i > ip1 && ix_dbnd[i] == ix_dbnd[i - 1] && iy_dbnd[i] == iy_dbnd[i - 1])
            continue;
          (*iz)[j] = ix_dbnd[i];
          (*ir)[j] = iy_dbnd[i];
          j++;
        }
      nPts = j;

      do
        {
          nNew = 0;
          /* Insert new points to ensure that only one coordinate changes at a time */
          for (i = 1; i < nPts; i++)
            {
              long iz0, iz1, ir0, ir1, izn, irn, new;
              iz0 = (*iz)[i - 1];
              iz1 = (*iz)[i];
              ir0 = (*ir)[i - 1];
              ir1 = (*ir)[i];
              if (iz0 != iz1 && ir0 != ir1)
                {
                  new = 1;
                  if (lattice[iz0][ir1] & FL_IS_INSIDE)
                    {
                      izn = iz0;
                      irn = ir1;
                    }
                  else if (lattice[iz1][ir0] & FL_IS_INSIDE)
                    {
                      izn = iz1;
                      irn = ir0;
                    }
                  else if (lattice[iz0][ir1] & FL_IS_SURFACE)
                    {
                      izn = iz0;
                      irn = ir1;
                    }
                  else if (lattice[iz1][ir0] & FL_IS_SURFACE)
                    {
                      izn = iz1;
                      irn = ir0;
                    }
                  else if (iz0 == 0)
                    {
                      izn = iz0;
                      irn = ir1;
                    }
                  else
                    {
                      new = 0;
                      printf("boundary issue: (%ld, %ld) followed by (%ld, %ld)\n",
                             iz0, ir0, iz1, ir1);
                      printf("(%ld, %ld): %d\n", iz0, ir1, lattice[iz0][ir1]);
                      printf("(%ld, %ld): %d\n", iz1, ir0, lattice[iz1][ir0]);
                    }
                  if (new)
                    {
                      for (j = nPts; j > i; j--)
                        {
                          (*iz)[j] = (*iz)[j - 1];
                          (*ir)[j] = (*ir)[j - 1];
                        }
                      (*iz)[i] = izn;
                      (*ir)[i] = irn;
                      i++;
                      nPts++;
                      nNew++;
                    }
                }
            }
          printf("Added %ld boundary points\n", nNew);
        }
      while (nNew != 0);
      return nPts;
    }
  return 0;
}

void dumpBoundaryDebugData(char *filename, long *iz, long *ir, short *direction, long nPts)
{
  FILE *fp;
  long i;
  static char directionName[4] = {'r', 'u', 'l', 'd'};
  fp = fopen(filename, "w");
  fprintf(fp, "SDDS1\n&column name=iz type=long &end\n&column name=ir type=long &end\n");
  fprintf(fp, "&column name=z type=double &end\n&column name=r type=double &end\n");
  fprintf(fp, "&column name=direction type=long &end\n");
  fprintf(fp, "&data mode=ascii &end\n%ld\n", nPts);
  for (i = 0; i < nPts; i++)
    fprintf(fp, "%ld %ld %le %le %hd %c\n", iz[i], ir[i], iz[i] * dx + xmin, ir[i] * dy + ymin, direction ? direction[i] : (short)-1,
            (direction && direction[i] > 0) ? directionName[direction[i] - 1] : '?');
  fclose(fp);
}

void set_up_drawing(
                    int _nx,
                    int _ny,
                    double _xmin,
                    double _xmax,
                    double _ymin,
                    double _ymax,
                    double tolerance_fraction,
                    short **_lattice,
                    double **_potential,
                    short **_material)
{
  int i, j;

  nx = _nx;
  ny = _ny;
  xmin = _xmin;
  ymin = _ymin;
  xmax = _xmax;
  ymax = _ymax;
  x_tolerance = (dx = (xmax - xmin) / (nx - 1)) * tolerance_fraction;
  y_tolerance = (dy = (ymax - ymin) / (ny - 1)) * tolerance_fraction;
  lattice = _lattice;
  potential = _potential;
  material = _material;
  i_point = 0;
  n_points_max = ARRAY_SIZE_INCREMENT;
  x_bnd = array_1d(sizeof(double), 0, n_points_max - 1);
  y_bnd = array_1d(sizeof(double), 0, n_points_max - 1);
  x_dbnd = array_1d(sizeof(double), 0, n_points_max - 1);
  y_dbnd = array_1d(sizeof(double), 0, n_points_max - 1);
  ix_dbnd = array_1d(sizeof(long), 0, n_points_max - 1);
  iy_dbnd = array_1d(sizeof(long), 0, n_points_max - 1);
  material_bnd = array_1d(sizeof(short), 0, n_points_max - 1);

  index_brkpt = NULL;
  n_brkpts = 0;

  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      material[i][j] = lattice[i][j] = 0;
}

void starting_point(double x, double y, double potential, short material)
{
  add_point(x, y, potential, material);
}

void add_breakpoint(double x, double y, double potential, short material)
{
  add_point(x, y, potential, material);
  index_brkpt = trealloc(index_brkpt, (n_brkpts + 1) * sizeof(*index_brkpt));
  index_brkpt[n_brkpts] = i_point - 1;
  n_brkpts++;
}

void draw_line(
               double xi, double yi, /* initial coordinates */
               double xf, double yf, /* final coordinates */
               double potential,
               short material,
               double startPotential)
{
  double slope, intercept;
  double x, y, _dx, _dy;
  double potSlope;

  if (fabs(xi - xf) >= dx)
    {
      slope = (yf - yi) / (xf - xi);
      if (fabs(yf - yi) < dy / 1e6)
        slope = 0;
      intercept = yi - slope * xi;
      _dx = dx;
      if (fabs(slope * _dx) > fabs(dy))
        _dx = dy * .5 / fabs(slope);
      potSlope = (potential - startPotential) / (xf - xi);
      _dx = 0.5 * _dx;
      if (xi <= xf)
        {
          x = xi;
          do
            {
              y = slope * x + intercept;
              add_point(x, y, startPotential + potSlope * (x - xi), material);
            }
          while ((x += _dx) <= xf);
          add_point(xf, yf, potential, material);
        }
      else
        {
          x = xi;
          do
            {
              y = slope * x + intercept;
              add_point(x, y, startPotential + potSlope * (x - xi), material);
            }
          while ((x -= _dx) >= xf);
          add_point(xf, yf, potential, material);
        }
    }
  else if (fabs(yi - yf) >= dy)
    {
      slope = (xf - xi) / (yf - yi);
      if (fabs(xf - xi) < dx / 1e6)
        slope = 0;
      intercept = xi - slope * yi;
      potSlope = (potential - startPotential) / (yf - yi);
      _dy = dy * 0.5;
      if (yi <= yf)
        {
          y = yi;
          do
            {
              x = slope * y + intercept;
              add_point(x, y, startPotential + potSlope * (y - yi), material);
            }
          while ((y += _dy) <= yf);
          add_point(xf, yf, potential, material);
        }
      else
        {
          y = yi;
          do
            {
              x = slope * y + intercept;
              add_point(x, y, startPotential + potSlope * (y - yi), material);
            }
          while ((y -= _dy) >= yf);
          add_point(xf, yf, potential, material);
        }
    }
  else
    {
      add_point(xf, yf, potential, material);
    }
}

void draw_circle(
                 double xi, double yi,    /* initial point on arc */
                 double xc, double yc,    /* center of circle     */
                 double r, double thetaf, /* radius, final angle */
                 double *xf, double *yf,  /* for return of final point on arc */
                 double potential,
                 short material)
{
  double theta, x, y, _dtheta, theta0;
  int i, n, sense;

  theta = atan2(yi - yc, xi - xc);
  if (theta < thetaf)
    sense = 1;
  else
    sense = 0;

  if (r != 0)
    {
      /* set increment in theta for drawing arcs of circles */
      dtheta = MIN(dx, dy) / r / 2;
    }
  else
    {
      if (dx < dy)
        dtheta = dx / (xmax - xmin) / 2;
      else
        dtheta = dy / (ymax - ymin) / 2;
    }

  if (sense)
    {
      _dtheta = dtheta;
      while (theta > thetaf)
        theta -= PIx2;
    }
  else
    {
      while (theta < thetaf)
        theta += PIx2;
      _dtheta = -dtheta;
    }

  if ((n = fabs((thetaf - theta) / _dtheta)) < 15)
    _dtheta = SIGN(_dtheta) * fabs(thetaf - theta) / (n = 15);

  _dtheta *= 0.5;
  n *= 2;

  theta0 = theta;

  for (i = 1; i < n; i++)
    {
      x = xc + r * cos(theta);
      y = yc + r * sin(theta);
      add_point(x, y, potential, material);
      theta = theta0 + i * _dtheta;
    }

  x = xc + r * cos((thetaf + theta) / 2);
  y = yc + r * sin((thetaf + theta) / 2);
  add_point(x, y, potential, material);

  x = xc + r * cos(thetaf);
  y = yc + r * sin(thetaf);
  add_point(x, y, potential, material);

  *xf = x;
  *yf = y;
}

void draw_ellipse(
                  double xi, double yi,   /* initial point on arc */
                  double xc, double yc,   /* center of circle     */
                  double a, double b,     /* semi-axes */
                  double thetaf,          /* final angle */
                  double *xf, double *yf, /* for return of final point on arc */
                  double potential,
                  short material)
{
  double theta, x, y, _dtheta, theta0;
  int i, n, sense;

  theta = atan2((yi - yc) / b, (xi - xc) / a);
  if (theta < thetaf)
    sense = 1;
  else
    sense = 0;

  dtheta = MIN(dx / a, dy / b) / 2;

  if (sense)
    {
      _dtheta = dtheta;
      while (theta > thetaf)
        theta -= PIx2;
    }
  else
    {
      while (theta < thetaf)
        theta += PIx2;
      _dtheta = -dtheta;
    }

  if ((n = fabs((thetaf - theta) / _dtheta)) < 15)
    _dtheta = SIGN(_dtheta) * fabs(thetaf - theta) / (n = 15);

  _dtheta *= 0.5;
  n *= 2;

  theta0 = theta;

  for (i = 1; i < n; i++)
    {
      x = xc + a * cos(theta);
      y = yc + b * sin(theta);
      add_point(x, y, potential, material);
      theta = theta0 + i * _dtheta;
    }

  x = xc + a * cos((thetaf + theta) / 2);
  y = yc + b * sin((thetaf + theta) / 2);
  add_point(x, y, potential, material);

  x = xc + a * cos(thetaf);
  y = yc + b * sin(thetaf);
  add_point(x, y, potential, material);

  *xf = x;
  *yf = y;
}

void set_vacuum_point(double x, double y)
{
  ix_vacuum = (x - xmin) / dx + 0.5;
  iy_vacuum = (y - ymin) / dy + 0.5;
  vacuum_point_defined = 1;
  if (ix_vacuum < 0 || ix_vacuum >= nx || iy_vacuum < 0 || iy_vacuum >= ny)
    {
      printf("vacuum point out of bounds: z = %e, r = %e\n", x, y);
      n_outside++;
    }
}

/* routine: add_point()
 * purpose: add a point to arrays of boundary points 
 */

void add_point(double x, double y, double potential0, short material0)
{
  int i, j;

  if (i_point != 0 && x_bnd[i_point - 1] == (double)x && y_bnd[i_point - 1] == (double)y)
    return;

  if (i_point >= n_points_max)
    {
      n_points_max += ARRAY_SIZE_INCREMENT;
      x_bnd = (double *)trealloc(x_bnd, n_points_max * sizeof(double));
      y_bnd = (double *)trealloc(y_bnd, n_points_max * sizeof(double));
      x_dbnd = (double *)trealloc(x_dbnd, n_points_max * sizeof(double));
      y_dbnd = (double *)trealloc(y_dbnd, n_points_max * sizeof(double));
      ix_dbnd = (long *)trealloc(ix_dbnd, n_points_max * sizeof(long));
      iy_dbnd = (long *)trealloc(iy_dbnd, n_points_max * sizeof(long));
      material_bnd = (short *)trealloc(material_bnd, n_points_max * sizeof(short));
    }

  i = (x - xmin) / dx + 0.5;
  j = (y - ymin) / dy + 0.5;
  if (i >= 0 && j >= 0 && i < nx && j < ny)
    {
      lattice[i][j] = FL_IS_SURFACE;
      potential[i][j] = potential0;
      if (material[i][j] == 0)
        material[i][j] = material0;
    }
  else
    {
      printf("point out of bounds: z = %e, r = %e\n", x, y);
      n_outside++;
    }

  x_bnd[i_point] = x;
  y_bnd[i_point] = y;
  x_dbnd[i_point] = (ix_dbnd[i_point] = i) * dx + xmin;
  y_dbnd[i_point] = (iy_dbnd[i_point] = j) * dy + ymin;
  material_bnd[i_point] = material[i][j];
  i_point++;
}

long out_of_bounds(void)
{
  return (n_outside);
}

void dump_boundary_points(char *filename0, char *rootname)
{
  FILE *fp;
  int i, ib, ip1, ip2;
  char *filename;
  filename = compose_filename(filename0, rootname);
  fp = fopen_e(filename, "w", 0);
  if (filename != filename0)
    free(filename);
  fprintf(fp, "SDDS1\n&parameter name=spiffeOutputName, type=string, fixed_value=\"input boundary\" &end\n");
  fprintf(fp, "&column name=z, units=m, type=double &end\n");
  fprintf(fp, "&column name=r, units=m, type=double &end\n");
  fprintf(fp, "&column name=materialID, type=short &end\n&data mode=ascii &end\n");
  if (n_brkpts == 0)
    {
      fprintf(fp, "%d\n", i_point);
      /* dump x and y pairs, in meters */
      for (i = 0; i < i_point; i++)
        fprintf(fp, "%.6e %.6e %hd\n", x_bnd[i], y_bnd[i], material_bnd[i]);
    }
  else
    {
      for (ib = 0; ib <= n_brkpts; ib++)
        {
          if (ib == 0)
            {
              ip1 = 0;
              ip2 = index_brkpt[ib] - 1;
            }
          else if (ib == n_brkpts)
            {
              ip1 = index_brkpt[ib - 1];
              ip2 = i_point - 1;
            }
          else
            {
              ip1 = index_brkpt[ib - 1];
              ip2 = index_brkpt[ib] - 1;
            }
          fprintf(fp, "%d\n", ip2 - ip1 + 1);
          for (i = ip1; i <= ip2; i++)
            {
              fprintf(fp, "%.6e %.6e %hd\n", x_bnd[i], y_bnd[i], material_bnd[i]);
            }
        }
    }
  fclose(fp);
}

long getBoundaryPoints(long materialID, long **iz, long **ir, short **dir)
{
  long i, j, ip1, ip2, nPts;

  for (i = 0; i < i_point; i++)
    if (material_bnd[i] == materialID)
      break;
  if (i == i_point)
    return 0;
  ip1 = i;

  i += 1;
  for (; i < i_point; i++)
    if (material_bnd[i] != materialID)
      break;
  ip2 = i;

  if (ip2 != ip1)
    {
      long nNew;
      /* Record the (iz, ir) coordinates of the discretized boundary, skipping duplicates */
      *iz = tmalloc(10 * sizeof(**iz) * (ip2 - ip1));
      *ir = tmalloc(10 * sizeof(**ir) * (ip2 - ip1));
      j = 0;
      for (i = ip1; i < ip2; i++)
        {
          if (i > ip1 && ix_dbnd[i] == ix_dbnd[i - 1] && iy_dbnd[i] == iy_dbnd[i - 1])
            continue;
          (*iz)[j] = ix_dbnd[i];
          (*ir)[j] = iy_dbnd[i];
          j++;
        }
      nPts = j;
#ifdef DEBUG
      dumpBoundaryDebugData("dump1.sdds", *iz, *ir, NULL, nPts);
#endif

      /* Insert new points to ensure that only one coordinate changes at a time */
      do
        {
          nNew = 0;
          /* Insert new points to ensure that only one coordinate changes at a time */
          for (i = 1; i < nPts; i++)
            {
              long iz0, iz1, ir0, ir1, izn, irn, new;
              iz0 = (*iz)[i - 1];
              iz1 = (*iz)[i];
              ir0 = (*ir)[i - 1];
              ir1 = (*ir)[i];
              if (iz0 != iz1 && ir0 != ir1)
                {
                  new = 1;
                  if (lattice[iz0][ir1] & FL_IS_INSIDE)
                    {
                      izn = iz0;
                      irn = ir1;
                    }
                  else if (lattice[iz1][ir0] & FL_IS_INSIDE)
                    {
                      izn = iz1;
                      irn = ir0;
                    }
                  else if (lattice[iz0][ir1] & FL_IS_SURFACE)
                    {
                      izn = iz0;
                      irn = ir1;
                    }
                  else if (lattice[iz1][ir0] & FL_IS_SURFACE)
                    {
                      izn = iz1;
                      irn = ir0;
                    }
                  else if (iz0 == 0)
                    {
                      izn = iz0;
                      irn = ir1;
                    }
                  else
                    {
                      new = 0;
                      printf("boundary issue: (%ld, %ld) followed by (%ld, %ld)\n",
                             iz0, ir0, iz1, ir1);
                      printf("(%ld, %ld): %d\n", iz0, ir1, lattice[iz0][ir1]);
                      printf("(%ld, %ld): %d\n", iz1, ir0, lattice[iz1][ir0]);
                    }
                  if (new)
                    {
                      for (j = nPts; j > i; j--)
                        {
                          (*iz)[j] = (*iz)[j - 1];
                          (*ir)[j] = (*ir)[j - 1];
                        }
                      (*iz)[i] = izn;
                      (*ir)[i] = irn;
                      i++;
                      nPts++;
                      nNew++;
                    }
                }
            }
          printf("Added %ld boundary points\n", nNew);
        }
      while (nNew != 0);

#ifdef DEBUG
      dumpBoundaryDebugData("dump2.sdds", *iz, *ir, NULL, nPts);
#endif

      *dir = tmalloc(sizeof(**dir) * nPts);
      /* Determine which way the boundary faces */
      for (i = 1; i < nPts; i++)
        {
          long iz1, ir1, iz2, ir2;
          iz1 = (*iz)[i - 1];
          ir1 = (*ir)[i - 1];
          iz2 = (*iz)[i];
          ir2 = (*ir)[i];
          (*dir)[i - 1] = 0;
          if (iz1 == iz2)
            {
              /* z = constant segment */
              if (iz1 == 0 || lattice[iz1 - 1][ir2] & FL_IS_INSIDE || lattice[iz1 - 1][ir1] & FL_IS_INSIDE)
                {
                  (*dir)[i - 1] = EMITS_RIGHT;
                }
              else if ((*iz)[i] == (nx - 1) || lattice[iz1 + 1][ir2] & FL_IS_INSIDE || lattice[iz1 + 1][ir1] & FL_IS_INSIDE)
                {
                  (*dir)[i - 1] = EMITS_LEFT;
                }
              else
                {
                  if (iz1 > 1 && (lattice[iz1 - 2][ir2] & FL_IS_INSIDE || lattice[iz1 - 2][ir1] & FL_IS_INSIDE))
                    {
                      (*dir)[i - 1] = EMITS_RIGHT;
                    }
                  else if (iz1 < (nx - 2) && (lattice[iz1 + 2][ir2] & FL_IS_INSIDE || lattice[iz1 + 2][ir1] & FL_IS_INSIDE))
                    {
                      (*dir)[i - 1] = EMITS_LEFT;
                    }
                  else
                    {
                      if (lattice[iz1 - 1][ir2] & FL_IS_METAL || lattice[iz1 - 1][ir1] & FL_IS_METAL)
                        {
                          (*dir)[i - 1] = EMITS_RIGHT;
                        }
                      else if (lattice[iz1 + 1][ir2] & FL_IS_METAL || lattice[iz1 + 1][ir1] & FL_IS_METAL)
                        {
                          (*dir)[i - 1] = EMITS_LEFT;
                        }
                      else
                        {
                          printf("Can't determine which way constant-z surface at z=%le, r=%le, iz=%ld, ir=%ld emits.\n",
                                 (*iz)[i - 1] * dx + xmin, (*ir)[i - 1] * dy + ymin,
                                 (*iz)[i - 1], (*ir)[i - 1]);
                          bomb(NULL, NULL);
                        }
                    }
                }
            }
          else if ((*ir)[i] == (*ir)[i - 1])
            {
              /* r = constant segment */
              if (ir1 == 0 || lattice[iz1][ir1 - 1] & FL_IS_INSIDE || lattice[iz2][ir1 - 1] & FL_IS_INSIDE)
                {
                  (*dir)[i - 1] = EMITS_UP;
                }
              else if ((*ir)[i] == (ny - 1) || lattice[iz1][ir1 + 1] & FL_IS_INSIDE || lattice[iz2][ir1 + 1] & FL_IS_INSIDE)
                {
                  (*dir)[i - 1] = EMITS_DOWN;
                }
              else
                {
                  if (ir1 > 1 && (lattice[iz1][ir1 - 2] & FL_IS_INSIDE || lattice[iz2][ir1 - 2] & FL_IS_INSIDE))
                    {
                      (*dir)[i - 1] = EMITS_UP;
                    }
                  else if (ir1 < (ny - 2) && (lattice[iz1][ir1 + 2] & FL_IS_INSIDE || lattice[iz2][ir1 + 2] & FL_IS_INSIDE))
                    {
                      (*dir)[i - 1] = EMITS_DOWN;
                    }
                  else
                    {
                      if (lattice[iz1][ir1 - 1] & FL_IS_METAL || lattice[iz2][ir1 - 1] & FL_IS_METAL)
                        {
                          (*dir)[i - 1] = EMITS_UP;
                        }
                      else if (lattice[iz1][ir1 + 1] & FL_IS_METAL || lattice[iz2][ir1 + 1] & FL_IS_METAL)
                        {
                          (*dir)[i - 1] = EMITS_DOWN;
                        }
                      else
                        {
                          printf("Can't determine which way constant-r surface at z=%le, r=%le, iz=%ld, ir=%ld emits.\n",
                                 (*iz)[i - 1] * dx + xmin, (*ir)[i - 1] * dy + ymin,
                                 (*iz)[i - 1], (*ir)[i - 1]);
                          bomb(NULL, NULL);
                        }
                    }
                }
            }
          else if (i > 1)
            (*dir)[i - 1] = (*dir)[i - 2];
        }
#ifdef DEBUG
      dumpBoundaryDebugData("dump3.sdds", *iz, *ir, *dir, nPts);
#endif
      return nPts;
    }
  return 0;
}

void dump_dboundary_points(char *filename0, char *rootname)
{
  FILE *fp;
  char *filename;
  long *iz, *ir, i, nbp;
  if (!(nbp = getAllDiscreteBoundaryPoints(&iz, &ir)))
    return;

  filename = compose_filename(filename0, rootname);
  fp = fopen_e(filename, "w", 0);
  if (filename != filename0)
    free(filename);
  fprintf(fp, "SDDS1\n&parameter name=spiffeOutputName, type=string, fixed_value=\"input boundary\" &end\n");
  fprintf(fp, "&column name=iz, type=long &end\n");
  fprintf(fp, "&column name=ir, type=long &end\n");
  fprintf(fp, "&column name=z, units=m, type=double &end\n");
  fprintf(fp, "&column name=r, units=m, type=double &end\n");
  fprintf(fp, "&data mode=ascii &end\n");
  fprintf(fp, "%ld\n", nbp);
  for (i = 0; i < nbp; i++)
    fprintf(fp, "%ld %ld %le %le\n", iz[i], ir[i], xmin + iz[i] * dx, ymin + ir[i] * dy);
  fclose(fp);
  free(iz);
  free(ir);
}

void dumpSurfacePointData(char *rootname, FIELDS *EM)
{
  char filename[1000];
  long iz, ir;
  FILE *fp;

  sprintf(filename, "%s.surf", rootname);
  fp = fopen_e(filename, "w", 0);
  fprintf(fp, "SDDS1\n");
  fprintf(fp, "&column name=z type=double units=m &end\n");
  fprintf(fp, "&column name=r type=double units=m &end\n");
  fprintf(fp, "&column name=iz type=short &end\n");
  fprintf(fp, "&column name=ir type=short &end\n");
  fprintf(fp, "&column name=surface type=short &end\n");
  fprintf(fp, "&column name=potential type=double &end\n");
  fprintf(fp, "&column name=materialID type=short &end\n");
  fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  for (iz = 0; iz < EM->nz; iz++)
    {
      for (ir = 0; ir < EM->nr; ir++)
        {
          if (EM->metal_flags[iz][ir] & FL_IS_SURFACE || EM->materialID[iz][ir])
            fprintf(fp, "%e %e %ld %ld %hd %e %hd\n",
                    iz * EM->dz + EM->zmin,
                    ir * EM->dr, iz, ir,
                    EM->metal_flags[iz][ir] & FL_IS_SURFACE ? (short)1 : (short)0,
                    EM->PsiImposed[iz][ir],
                    EM->materialID[iz][ir]);
        }
    }
  fclose(fp);
}

void dumpGridPointData(char *rootname, FIELDS *EM)
{
  char filename[1000];
  long iz, ir;
  FILE *fp;

  sprintf(filename, "%s.grid", rootname);
  fp = fopen_e(filename, "w", 0);
  fprintf(fp, "SDDS1\n");
  fprintf(fp, "&column name=z type=double units=m &end\n");
  fprintf(fp, "&column name=r type=double units=m &end\n");
  fprintf(fp, "&column name=iz type=short &end\n");
  fprintf(fp, "&column name=ir type=short &end\n");
  fprintf(fp, "&column name=flags type=short &end\n");
  fprintf(fp, "&column name=surface type=short &end\n");
  fprintf(fp, "&column name=metal type=short &end\n");
  fprintf(fp, "&column name=inside type=short &end\n");
  fprintf(fp, "&column name=materialID type=short &end\n");
  fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  for (iz = 0; iz < EM->nz; iz++)
    {
      for (ir = 0; ir < EM->nr; ir++)
        {
          fprintf(fp, "%e %e %ld %ld %hd %hd %hd %hd %hd\n",
                  iz * EM->dz + EM->zmin,
                  ir * EM->dr, iz, ir,
                  EM->metal_flags[iz][ir],
                  EM->metal_flags[iz][ir] & FL_IS_SURFACE ? (short)1 : (short)0,
                  EM->metal_flags[iz][ir] & FL_IS_METAL ? (short)1 : (short)0,
                  EM->metal_flags[iz][ir] & FL_IS_INSIDE ? (short)1 : (short)0,
                  EM->materialID[iz][ir]);
        }
    }
  fclose(fp);
}

void dumpMetalPointData(char *rootname, FIELDS *EM)
{
  char filename[1000];
  long iz, ir;
  FILE *fp;

  sprintf(filename, "%s.metal", rootname);
  fp = fopen_e(filename, "w", 0);
  fprintf(fp, "SDDS1\n");
  fprintf(fp, "&column name=z type=double units=m &end\n");
  fprintf(fp, "&column name=r type=double units=m &end\n");
  fprintf(fp, "&column name=iz type=short &end\n");
  fprintf(fp, "&column name=ir type=short &end\n");
  fprintf(fp, "&column name=surface type=short &end\n");
  fprintf(fp, "&column name=metal type=short &end\n");
  fprintf(fp, "&column name=inside type=short &end\n");
  fprintf(fp, "&column name=materialID type=short &end\n");
  fprintf(fp, "&data mode=ascii no_row_counts=1 &end\n");
  for (iz = 0; iz < EM->nz; iz++)
    {
      for (ir = 0; ir < EM->nr; ir++)
        {
          if (EM->metal_flags[iz][ir] & FL_IS_INSIDE)
            fprintf(fp, "%e %e %ld %ld %hd %hd %hd %hd\n",
                    iz * EM->dz + EM->zmin,
                    ir * EM->dr, iz, ir,
                    EM->metal_flags[iz][ir] & FL_IS_SURFACE ? (short)1 : (short)0,
                    EM->metal_flags[iz][ir] & FL_IS_METAL ? (short)1 : (short)0,
                    EM->metal_flags[iz][ir] & FL_IS_INSIDE ? (short)1 : (short)0,
                    EM->materialID[iz][ir]);
        }
    }
  fclose(fp);
  dumpGridPointData(rootname, EM);
}

/* routine: fill_in_drawing()
 * purpose: fill in a drawing, defined by the closed list of boundary points
 *          (x_bnd, y_bnd)
 *
 */

void fill_in_drawing(short **_lattice, double **_potential)
{
  int i, j, i1, j1, changed, change;
  double potentialSum;
  long nSummed;

  if (!vacuum_point_defined)
    {
      /* move along the x axis and find a point that is outside the metal region */
      /* -- first see if there are points marked on the lower boundary */
      for (i = 0; i < nx; i++)
        if (_lattice[i][0] && !_lattice[i][1])
          break;
      if (i == nx)
        bomb("can't determine what points are in vacuum--use nt=4 point to mark a point that is in vacuum", NULL);
      if (i > 0)
        _lattice[i][1] = FL_IS_VACUUM;
    }
  else
    _lattice[ix_vacuum][iy_vacuum] = FL_IS_VACUUM;

  /* sweep over the mesh and fill in all exterior points (i.e., points exterior to metal) */
  changed = 1;
  while (changed)
    {
      changed = 0;
      for (i = 1; i < nx - 1; i++)
        {
          for (j = 1; j < ny - 1; j++)
            {
              if (_lattice[i][j] == FL_IS_VACUUM)
                {
                  if (_lattice[i][j + 1] == 0)
                    {
                      changed++;
                      _lattice[i][j + 1] = FL_IS_VACUUM;
                    }
                  if (_lattice[i][j - 1] == 0)
                    {
                      changed++;
                      _lattice[i][j - 1] = FL_IS_VACUUM;
                    }
                  if (_lattice[i + 1][j] == 0)
                    {
                      changed++;
                      _lattice[i + 1][j] = FL_IS_VACUUM;
                    }
                  if (_lattice[i - 1][j] == 0)
                    {
                      changed++;
                      _lattice[i - 1][j] = FL_IS_VACUUM;
                    }
                }
            }
          for (j = ny - 2; j > 0; j--)
            {
              if (_lattice[i][j] == FL_IS_VACUUM)
                {
                  if (_lattice[i][j + 1] == 0)
                    {
                      changed++;
                      _lattice[i][j + 1] = FL_IS_VACUUM;
                    }
                  if (_lattice[i][j - 1] == 0)
                    {
                      changed++;
                      _lattice[i][j - 1] = FL_IS_VACUUM;
                    }
                  if (_lattice[i + 1][j] == 0)
                    {
                      changed++;
                      _lattice[i + 1][j] = FL_IS_VACUUM;
                    }
                  if (_lattice[i - 1][j] == 0)
                    {
                      changed++;
                      _lattice[i - 1][j] = FL_IS_VACUUM;
                    }
                }
            }
        }
      for (i = nx - 2; i > 0; i--)
        {
          for (j = 1; j < ny - 1; j++)
            {
              if (_lattice[i][j] == FL_IS_VACUUM)
                {
                  if (_lattice[i][j + 1] == 0)
                    {
                      changed++;
                      _lattice[i][j + 1] = FL_IS_VACUUM;
                    }
                  if (_lattice[i][j - 1] == 0)
                    {
                      changed++;
                      _lattice[i][j - 1] = FL_IS_VACUUM;
                    }
                  if (_lattice[i + 1][j] == 0)
                    {
                      changed++;
                      _lattice[i + 1][j] = FL_IS_VACUUM;
                    }
                  if (_lattice[i - 1][j] == 0)
                    {
                      changed++;
                      _lattice[i - 1][j] = FL_IS_VACUUM;
                    }
                }
            }
          for (j = ny - 2; j > 0; j--)
            {
              if (_lattice[i][j] == FL_IS_VACUUM)
                {
                  if (_lattice[i][j + 1] == 0)
                    {
                      changed++;
                      _lattice[i][j + 1] = FL_IS_VACUUM;
                    }
                  if (_lattice[i][j - 1] == 0)
                    {
                      changed++;
                      _lattice[i][j - 1] = FL_IS_VACUUM;
                    }
                  if (_lattice[i + 1][j] == 0)
                    {
                      changed++;
                      _lattice[i + 1][j] = FL_IS_VACUUM;
                    }
                  if (_lattice[i - 1][j] == 0)
                    {
                      changed++;
                      _lattice[i - 1][j] = FL_IS_VACUUM;
                    }
                }
            }
        }
      printf("%d vacuum points found\n", changed);
    }

  /* unmark the vacuum points and mark the interior points */
  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
        {
          if (_lattice[i][j] == 0)
            _lattice[i][j] = FL_IS_INSIDE;
          else if (_lattice[i][j] == FL_IS_VACUUM)
            _lattice[i][j] = 0;
        }
    }

  /* make sure there are no points marked as inside that are really on the surface */
  for (i = 1; i < nx - 1; i++)
    {
      for (j = 1; j < ny - 1; j++)
        {
          if (_lattice[i][j] == FL_IS_INSIDE)
            {
              change = 0;
              potentialSum = 0;
              nSummed = 0;
              for (i1 = -1; i1 <= 1; i1++)
                {
                  for (j1 = -1; j1 <= 1; j1++)
                    {
                      if (i1 == 0 && j1 == 0)
                        continue;
                      if (!_lattice[i + i1][j + j1])
                        {
                          change++;
                        }
                      else
                        {
                          if (_lattice[i + i1][j + j1] == FL_IS_SURFACE)
                            {
                              potentialSum += _potential[i + i1][j + j1];
                              nSummed++;
                            }
                        }
                    }
                }
              if (change == 8)
                _lattice[i][j] = 0;
              else if (change)
                {
                  _lattice[i][j] = FL_IS_SURFACE;
                  if (nSummed)
                    _potential[i][j] = potentialSum / nSummed;
                  else
                    _potential[i][j] = 0;
                  printf("converted %d, %d  to surface with potential %e\n", i, j, _potential[i][j]);
                }
            }
        }
    }
}

void type_drawing(FILE *fp, char *label, short **_lattice, int nx, int ny)
{
  int i, j;

  fputc('\014', fp);
  if (label != NULL)
    fputs(label, fp);
  fputc('\n', fp);

  for (i = ny - 1; i >= 0; i--)
    {
      if (!(i % 5))
        fprintf(fp, "%3d-", i);
      else
        fputs("    ", fp);
      for (j = 0; j < nx; j++)
        {
          if (_lattice[j][i])
            fputc('#', fp);
          else
            fputc(' ', fp);
        }
      fputc('\n', fp);
    }
  fputs("    ", fp);
  for (j = 0; j < nx; j++)
    if (!(j % 5))
      fputc('|', fp);
    else
      fputc(' ', fp);
  fputs("\n  ", fp);
  for (j = 0; j < nx; j += 10)
    fprintf(fp, "%3d       ", j);
  fputc('\n', fp);
  fflush(fp);
}

void type_flag_drawing(FILE *fp, char *label, short **_lattice, int nx, int ny)
{
  int i, j;

  fputc('\014', fp);
  if (label != NULL)
    fputs(label, fp);
  fputc('\n', fp);

  for (i = ny; i >= 0; i--)
    {
      if (!(i % 5))
        fprintf(fp, "%3d-", i);
      else
        fputs("    ", fp);
      for (j = 0; j <= nx; j++)
        {
          if (_lattice[j][i])
            fprintf(fp, "%1x", _lattice[j][i]);
          else
            fputc(' ', fp);
        }
      fputc('\n', fp);
    }
  fputs("    ", fp);
  for (j = 0; j <= nx; j++)
    if (!(j % 5))
      fputc('|', fp);
    else
      fputc(' ', fp);
  fputs("\n  ", fp);
  for (j = 0; j <= nx; j += 10)
    fprintf(fp, "%3d       ", j);
  fputc('\n', fp);
  fflush(fp);
}

/* routine: trim_drawing()
 * purpose: eliminate any single-width filled areas that run along a
 *          boundary.
 */

void trim_drawing(short **_lattice)
{
  int i, j;

  /* left boundary */
  for (j = 0; j < ny; j++)
    if (_lattice[0][j] && !_lattice[1][j])
      _lattice[0][j] = 0;

  /* right boundary */
  i = nx - 1;
  for (j = 0; j < ny; j++)
    if (_lattice[i][j] && !_lattice[i - 1][j])
      _lattice[i][j] = 0;

  /* bottom boundary */
  for (i = 0; i < nx; i++)
    if (_lattice[i][0] && !_lattice[i][1])
      _lattice[i][0] = 0;

  /* top boundary */
  /**
     j = ny-1;
     for (i=0; i<nx; i++)
     if (_lattice[i][j] && !_lattice[i][j-1])
     _lattice[i][j] = 0;
  **/
}

void dump_interior_points(char *filename0, char *rootname)
{
  FILE *fp;
  int ix, iy, n_filled;
  char *filename;

  n_filled = 0;
  for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++)
      if (lattice[ix][iy])
        n_filled++;

  filename = compose_filename(filename0, rootname);
  fp = fopen_e(filename, "w", 0);
  fprintf(fp, "SDDS1\n&parameter name=spiffeOutputName, type=string, fixed_value=\"interior points\" &end\n");
  fprintf(fp, "&column name=iz, type=long &end\n");
  fprintf(fp, "&column name=ir, type=long &end\n");
  fprintf(fp, "&column name=z, units=m, type=double &end\n");
  fprintf(fp, "&column name=r, units=m, type=double &end\n");
  fprintf(fp, "&column name=lastActivePoint, type=short &end\n");
  fprintf(fp, "&column name=code, type=short &end\n&data mode=ascii &end\n");
  fprintf(fp, "%d\n", n_filled);
  if (filename != filename0)
    free(filename);

  /* dump x and y pairs of interior points, in meters */
  for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++)
      if (lattice[ix][iy])
        fprintf(fp, "%d %d %.6e %.6e %d %d\n", ix, iy, ix * dx + xmin, iy * dy + ymin,
                lattice[ix][iy] & FL_LAST_ACTIVE_PT ? 1 : 0, lattice[ix][iy]);

  fclose(fp);
}
