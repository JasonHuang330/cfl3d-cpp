// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#pragma once
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <complex>
#include <cstdint>

namespace cgnstools_ns {

void dummycgns();
void dummy();

// --- Real CGNS read path (C API, libcgns), enough for the steady single-zone
//     delta case. Returns 0 on success, nonzero on error. ---
// Open CGNS file read-only; find a base with celldim==physdim==idimdesired;
// return file index iccg, base index ibase, and number of zones.
int ropencgns_c(const char* fname, int idimdesired,
                int& iccg, int& ibase, int& nzones);
// Read grid coords for zone igrid and transpose (i,j,k)->(j,k,i) into x,y,z.
// ialph==0: X->x, Y->y, Z->z. (ialph!=0 swaps y/z with sign, per getgrid.F.)
int getgrid_c(int iccg, int ibase, int igrid, int idim, int jdim, int kdim,
              int ialph, double* x, double* y, double* z);
// Close the CGNS file.
void closecgns_c(int iccg);

} // namespace cgnstools_ns
