// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rcfl.h"
#include "runtime/fortran_io.h"
#include <cstdio>

namespace rcfl_ns {

// Read one coordinate array in Fortran implied-DO order
// (((a(j,k,i),j=1,jdim),k=1,kdim),i=1,idim) from the current record.
static bool read_coord_array(FortranRecordReader& rec, FortranArray3DRef<double>& a,
                             int jdim, int kdim, int idim) {
    for (int i = 1; i <= idim; i++) {
        for (int k = 1; k <= kdim; k++) {
            for (int j = 1; j <= jdim; j++) {
                double val;
                if (!rec.read_scalar(val)) return false;
                a(j, k, i) = val;
            }
        }
    }
    return true;
}

void rcfl(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
          int& jdim, int& kdim, int& idim, int& igrid, int& ialph, int& igeom_img, int& irr)
{
//***********************************************************************
//     Purpose:  Read grids in cfl3d format.
//     ialph     - flag to interpret angle of attack
//               = 0 alpha measured in x-z plane (cfl3d standard)
//               > 0 alpha measured in x-y plane (tlns3d standard)
//     igeom_img - flag for complex grid file
//               = 0 input grid is real valued
//               > 0 input grid is complex valued
//***********************************************************************

    FILE* fp = fortran_get_unit(1);
    if (!fp) { irr = 1; return; }

    // Each Fortran read(1) below is ONE unformatted sequential record
    // holding all three coordinate arrays back-to-back.
    FortranRecordReader rec(fp);

    if (ialph == 0) {
        // read(1) (((x...)),(((y...)),(((z...)))   -- x,y,z order
        if (!read_coord_array(rec, x, jdim, kdim, idim)) { irr = 1; return; }
        if (!read_coord_array(rec, y, jdim, kdim, idim)) { irr = 1; return; }
        if (!read_coord_array(rec, z, jdim, kdim, idim)) { irr = 1; return; }
        if (!rec.finish()) { irr = 1; return; }
    } else {
        // read(1) (((x...)),(((z...)),(((y...)))   -- x,z,y order (y last)
        if (!read_coord_array(rec, x, jdim, kdim, idim)) { irr = 1; return; }
        if (!read_coord_array(rec, z, jdim, kdim, idim)) { irr = 1; return; }
        if (!read_coord_array(rec, y, jdim, kdim, idim)) { irr = 1; return; }
        if (!rec.finish()) { irr = 1; return; }

        // Negate y after reading (tlns3d standard)
        for (int i = 1; i <= idim; i++) {
            for (int j = 1; j <= jdim; j++) {
                for (int k = 1; k <= kdim; k++) {
                    y(j, k, i) = -y(j, k, i);
                }
            }
        }
    }

    return;
}

} // namespace rcfl_ns
