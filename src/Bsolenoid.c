/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

#include "spiffe.h"

/* routine: B_solenoid() 
 * purpose: calculate the magnetic field due to an equispaced series of rings 
 *          carrying 1A current.
 */

void B_solenoid(double *B_rho, double *B_z, double rho, double z,
                double radius, double zStart, double zEnd, long coils, long symmetry)
{
  long iCoil;
  double dzCoil, Brho1, Bz1, zCoil;
  *B_rho = *B_z = 0;
  if (coils > 1)
    dzCoil = (zEnd - zStart) / (coils - 1);
  else
    dzCoil = 0;
  for (iCoil = 0; iCoil < coils; iCoil++)
    {
      zCoil = zStart + iCoil * dzCoil;
      B_ring(&Brho1, &Bz1, rho, z - zCoil, radius);
      *B_rho += Brho1;
      *B_z += Bz1;
      if (symmetry && zCoil)
        {
          B_ring(&Brho1, &Bz1, rho, z + zCoil, radius);
          switch (symmetry)
            {
            case 1:
              *B_rho += Brho1;
              *B_z += Bz1;
              break;
            case -1:
              *B_rho -= Brho1;
              *B_z -= Bz1;
              break;
            }
        }
    }
}

/* routine: B_ring()
 * purpose: calculate the magnetic field due to 1A current in a ring, the
 *          expressions being those found in Landau & Lifshitz, ELECTRODYNAMICS
 *          OF CONTINUOUS MEDIA, pg 112.  
 */

void B_ring(double *B_rho, double *B_z, double rho, double z, double radius)
{
  double radius2, z2, r2;
  double K, E, F;
  double s1, s2, s3, s4;
  double k;

  radius2 = radius * radius;
  z2 = z * z;
  r2 = rho * rho;
  s1 = radius2 + r2 + z2;
  s2 = radius2 - r2 - z2;

  radius2 = radius + rho;
  s3 = radius2 * radius2 + z2;

  radius2 = radius - rho;
  s4 = radius2 * radius2 + z2;

  k = sqrt(4 * radius * rho / s3);
  K = K_cei(k);
  E = E_cei(k);
  F = 2e-7 / sqrt(s3); /* 2e-7 is mu_o/(2 pi) */

  if (rho)
    *B_rho = F * z / rho * (-K + E * s1 / s4);
  else
    *B_rho = 0;
  *B_z = F * (K + E * s2 / s4);
}
