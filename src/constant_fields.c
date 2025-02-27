/*************************************************************************\
 * Copyright (c) 2002 The University of Chicago, as Operator of Argonne
 * National Laboratory.
 * Copyright (c) 2002 The Regents of the University of California, as
 * Operator of Los Alamos National Laboratory.
 * This file is distributed subject to a Software License Agreement found
 * in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: constant_fields.c
 * purpose: allow setting of constant fields
 *
 * Michael Borland, 1992 
 */
#include "spiffe.h"
#include "constant_fields.h"

void set_constant_fields(
                         FIELDS *EM_problem,
                         NAMELIST_TEXT *nl_text /* unparsed parameters in namelist text format */
                         )
{

  if (!EM_problem)
    bomb("can't set constant fields before geometry is defined", NULL);

  Ez = Er = Ephi = 0;
  Bz = Br = Bphi = 0;

  process_namelist(&constant_fields, nl_text);
  print_namelist(stdout, &constant_fields);

  EM_problem->constantEz = Ez;
  EM_problem->constantEr = Er;
  EM_problem->constantEphi = Ephi;
  EM_problem->constantBz = Bz;
  EM_problem->constantBr = Br;
  EM_problem->constantBphi = Bphi;
}
