// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "cgnstools.h"
#include <cgnslib.h>
#include <vector>
#include <cstdio>
#include <cstring>

namespace cgnstools_ns {

// Open CGNS file read-only, find a base of the right dimensionality,
// return iccg/ibase/nzones. Mirrors cgnstools.F:ropencgns.
int ropencgns_c(const char* fname, int idimdesired,
                int& iccg, int& ibase, int& nzones)
{
    int ier = cg_open(fname, CG_MODE_READ, &iccg);
    if (ier) { std::fprintf(stderr, "cgns: cannot open '%s'\n", fname); cg_error_print(); return ier; }
    int nbases = 0;
    cg_nbases(iccg, &nbases);
    ibase = 0;
    for (int ib = 1; ib <= nbases; ib++) {
        char bn[64]; int celldim = 0, physdim = 0;
        cg_base_read(iccg, ib, bn, &celldim, &physdim);
        if (celldim != physdim) { std::fprintf(stderr, "cgns: celldim!=physdim in base %d\n", ib); return 1; }
        if (celldim == idimdesired) ibase = ib;
    }
    if (ibase == 0) ibase = 1;               // fall back to the first base
    cg_nzones(iccg, ibase, &nzones);
    return 0;
}

// Read grid coordinates for zone igrid, transposing CGNS (i,j,k) ordering to
// CFL3D (j,k,i) ordering. Mirrors cgnstools.F:getgrid + reorderg.
int getgrid_c(int iccg, int ibase, int igrid, int idim, int jdim, int kdim,
              int ialph, double* x, double* y, double* z)
{
    CGNS_ENUMT(ZoneType_t) zt = CGNS_ENUMV(ZoneTypeNull);
    cg_zone_type(iccg, ibase, igrid, &zt);
    if (zt != CGNS_ENUMV(Structured)) { std::fprintf(stderr, "cgns: zone %d not Structured\n", igrid); return 1; }
    char zn[64]; cgsize_t isize[9] = {0};
    cg_zone_read(iccg, ibase, igrid, zn, isize);
    if (isize[0] != idim || isize[1] != jdim || isize[2] != kdim) {
        std::fprintf(stderr, "cgns: zone %d dim mismatch cgns=%lld,%lld,%lld deck=%d,%d,%d\n",
                     igrid, (long long)isize[0], (long long)isize[1], (long long)isize[2], idim, jdim, kdim);
        return 1;
    }
    const long n = (long)idim * jdim * kdim;
    std::vector<double> wk((size_t)n);
    cgsize_t rmin[3] = {1, 1, 1};
    cgsize_t rmax[3] = {idim, jdim, kdim};
    // xyz(j,k,i) = sign * wk(i,j,k)   [both 0-based flat below]
    auto transpose = [&](double* dst, double sign) {
        for (int i = 1; i <= idim; i++)
        for (int k = 1; k <= kdim; k++)
        for (int j = 1; j <= jdim; j++)
            dst[(j-1) + (long)jdim*((k-1) + (long)kdim*(i-1))] =
                sign * wk[(size_t)((i-1) + (long)idim*((j-1) + (long)jdim*(k-1)))];
    };
    int ier;
    ier = cg_coord_read(iccg, ibase, igrid, "CoordinateX", CGNS_ENUMV(RealDouble), rmin, rmax, wk.data());
    if (ier) { cg_error_print(); return ier; }
    transpose(x, 1.0);
    ier = cg_coord_read(iccg, ibase, igrid, "CoordinateY", CGNS_ENUMV(RealDouble), rmin, rmax, wk.data());
    if (ier) { cg_error_print(); return ier; }
    if (ialph == 0) transpose(y, 1.0); else transpose(z, 1.0);
    ier = cg_coord_read(iccg, ibase, igrid, "CoordinateZ", CGNS_ENUMV(RealDouble), rmin, rmax, wk.data());
    if (ier) { cg_error_print(); return ier; }
    if (ialph == 0) transpose(z, 1.0); else transpose(y, -1.0);
    return 0;
}

void closecgns_c(int iccg) { cg_close(iccg); }

void dummycgns() {}
void dummy() {}

} // namespace cgnstools_ns
