/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* routine: solve_poisson_cyl()
 * purpose: solve poisson's equation in cylindrical coordinates by
 *          simultaneous over-relaxation
 * Michael Borland, 1992, 1994, 1995
 */
#include "spiffe.h"

long solve_poisson_cyl(
                       double **Psi,
                       double **driving_term, /* -rho/epsilon_o for error charge */
                       short **flags,
                       long left_bc, long right_bc,
                       long nz, long nr,
                       double dz, double dr,
                       double accuracy,
                       long iteration_limit)
{
  long iz, ir, n;
  double residual, omega, r_jacobi;
  double abs_total_driving_term, abs_total_residual = 0;
  double a, b, c, d, e;

  if (accuracy <= 0)
    bomb("accuracy is <= 0 in solve_poisson_cyl", NULL);
  if (iteration_limit <= 1)
    bomb("iteration_limit is <= 1 in solve_poisson_cyl", NULL);

  abs_total_driving_term = 0.0;
  if (driving_term)
    for (iz = 0; iz < nz; iz++)
      for (ir = 0; ir < nr - 1; ir++)
        abs_total_driving_term += fabs(driving_term[iz][ir]);
  /*
    if (abs_total_driving_term==0)
    return(0);
  */

  r_jacobi = (cos(PI / nz) + sqr(dz / dr) * cos(PI / nr)) / (1 + sqr(dz / dr));

  omega = 1.0;
  a = b = 1 / sqr(dz);
  e = -2 * (1 / sqr(dz) + 1 / sqr(dr));
  for (n = 0; n < iteration_limit; n++)
    {
      abs_total_residual = 0;
      /* these BC implementations should be improved by use of 
       * extra cells to the left and right 
       */
      for (ir = 0; ir < nr - 1; ir++)
        {
          if (ir != 0)
            {
              c = 1 / sqr(dr) * (1 + 1. / (2 * ir));
              d = 1 / sqr(dr) * (1 - 1. / (2 * ir));
            }
          else
            c = d = 1 / sqr(dr);
          if (left_bc == DIRICHLET)
            {
              iz = 0;
              if ((ir + iz) % 2 == n % 2 && !(flags[iz][ir] & FL_IS_METAL))
                {
                  if (ir == 0)
                    {
                      residual = a * (Psi[iz + 1][ir] + Psi[iz + 1][ir]) + 4 * c * Psi[iz][ir + 1] + (e - 2 * c) * Psi[iz][ir] - (driving_term ? driving_term[iz][ir] : 0);
                      Psi[iz][ir] -= omega * residual / (e - 2 * c);
                    }
                  else
                    {
                      residual = a * (Psi[iz + 1][ir] + Psi[iz + 1][ir]) + c * Psi[iz][ir + 1] + d * Psi[iz][ir - 1] + e * Psi[iz][ir] - (driving_term ? driving_term[iz][ir] : 0);
                      Psi[iz][ir] -= omega * residual / e;
                    }
                  abs_total_residual += fabs(residual);
                }
            }
          for (iz = 1; iz < nz - 1; iz++)
            {
              if ((ir + iz) % 2 == n % 2 && !(flags[iz][ir] & FL_IS_METAL))
                {
                  if (ir == 0)
                    {
                      residual = a * (Psi[iz + 1][ir] + Psi[iz - 1][ir]) + 4 * c * Psi[iz][ir + 1] + (e - 2 * c) * Psi[iz][ir] - (driving_term ? driving_term[iz][ir] : 0);
                      Psi[iz][ir] -= omega * residual / (e - 2 * c);
                    }
                  else
                    {
                      residual = a * (Psi[iz + 1][ir] + Psi[iz - 1][ir]) + c * Psi[iz][ir + 1] + d * Psi[iz][ir - 1] + e * Psi[iz][ir] - (driving_term ? driving_term[iz][ir] : 0);
                      Psi[iz][ir] -= omega * residual / e;
                    }
                  abs_total_residual += fabs(residual);
                }
            }
          if (right_bc == DIRICHLET)
            {
              iz = nz - 1;
              if ((ir + iz) % 2 == n % 2 && !(flags[iz][ir] & FL_IS_METAL))
                {
                  if (ir == 0)
                    {
                      residual = a * (Psi[iz - 1][ir] + Psi[iz - 1][ir]) + 4 * c * Psi[iz][ir + 1] + (e - 2 * c) * Psi[iz][ir] - (driving_term ? driving_term[iz][ir] : 0);
                      Psi[iz][ir] -= omega * residual / (e - 2 * c);
                    }
                  else
                    {
                      residual = a * (Psi[iz - 1][ir] + Psi[iz - 1][ir]) + c * Psi[iz][ir + 1] + d * Psi[iz][ir - 1] + e * Psi[iz][ir] - (driving_term ? driving_term[iz][ir] : 0);
                      Psi[iz][ir] -= omega * residual / e;
                    }
                  abs_total_residual += fabs(residual);
                }
            }
        }
      if (n == 0)
        omega = 1.0 / (1 - sqr(r_jacobi) / 2);
      else
        omega = 1.0 / (1 - sqr(r_jacobi) * omega / 4);
      if (n > 0 && (abs_total_residual < accuracy * abs_total_driving_term || abs_total_residual == 0))
        return (n);
    }
  printf("warning: maximum number of iterations exceeded in solve_poisson_cyl--residual was %e\n", abs_total_residual);
  return (-1);
}
