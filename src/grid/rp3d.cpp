// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rp3d.h"
#include "runtime/fortran_io.h"
#include <cstdint>
#include <cstdio>

namespace rp3d_ns {

void rp3d(FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
          int& jdim, int& kdim, int& idim, int& igrid, int& ialph, int& igeom_img, int& irr)
{
    // Purpose: Read grids in plot3d format from unit 1 (unformatted binary).
    // ialph = 0: alpha measured in x-z plane (cfl3d standard), read x,y,z
    // ialph > 0: alpha measured in x-y plane (tlns3d standard), read x,z,y then negate y
    // igeom_img: flag for complex grid file (not used in read logic here)
    // igrid: not used in this subroutine body

    FILE* fp = fortran_get_unit(1);
    if (!fp) {
        irr = 1;
        return;
    }

    // Fortran unformatted sequential: one READ = one record = [int32 marker][data][int32 marker]
    // The entire multi-array read is ONE record.

    int32_t rec_start, rec_end;

    if (ialph == 0) {
        // read(1,end=999,err=999)
        //   (((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
        //   (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
        //   (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)

        if (fread(&rec_start, sizeof(int32_t), 1, fp) != 1) { irr = 1; return; }

        // Read x: k outer, j middle, i inner
        for (int k = 1; k <= kdim; k++)
            for (int j = 1; j <= jdim; j++)
                for (int i = 1; i <= idim; i++) {
                    double val;
                    if (fread(&val, sizeof(double), 1, fp) != 1) { irr = 1; return; }
                    x(j, k, i) = val;
                }

        // Read y: k outer, j middle, i inner
        for (int k = 1; k <= kdim; k++)
            for (int j = 1; j <= jdim; j++)
                for (int i = 1; i <= idim; i++) {
                    double val;
                    if (fread(&val, sizeof(double), 1, fp) != 1) { irr = 1; return; }
                    y(j, k, i) = val;
                }

        // Read z: k outer, j middle, i inner
        for (int k = 1; k <= kdim; k++)
            for (int j = 1; j <= jdim; j++)
                for (int i = 1; i <= idim; i++) {
                    double val;
                    if (fread(&val, sizeof(double), 1, fp) != 1) { irr = 1; return; }
                    z(j, k, i) = val;
                }

        if (fread(&rec_end, sizeof(int32_t), 1, fp) != 1) { irr = 1; return; }

    } else {
        // read(1,end=999,err=999)
        //   (((x(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
        //   (((z(j,k,i),i=1,idim),j=1,jdim),k=1,kdim),
        //   (((y(j,k,i),i=1,idim),j=1,jdim),k=1,kdim)
        // Note: file stores x, z, y (swapped order)

        if (fread(&rec_start, sizeof(int32_t), 1, fp) != 1) { irr = 1; return; }

        // Read x: k outer, j middle, i inner
        for (int k = 1; k <= kdim; k++)
            for (int j = 1; j <= jdim; j++)
                for (int i = 1; i <= idim; i++) {
                    double val;
                    if (fread(&val, sizeof(double), 1, fp) != 1) { irr = 1; return; }
                    x(j, k, i) = val;
                }

        // Read z (stored in file's second position): k outer, j middle, i inner
        for (int k = 1; k <= kdim; k++)
            for (int j = 1; j <= jdim; j++)
                for (int i = 1; i <= idim; i++) {
                    double val;
                    if (fread(&val, sizeof(double), 1, fp) != 1) { irr = 1; return; }
                    z(j, k, i) = val;
                }

        // Read y (stored in file's third position): k outer, j middle, i inner
        for (int k = 1; k <= kdim; k++)
            for (int j = 1; j <= jdim; j++)
                for (int i = 1; i <= idim; i++) {
                    double val;
                    if (fread(&val, sizeof(double), 1, fp) != 1) { irr = 1; return; }
                    y(j, k, i) = val;
                }

        if (fread(&rec_end, sizeof(int32_t), 1, fp) != 1) { irr = 1; return; }

        // Negate y after reading
        for (int i = 1; i <= idim; i++)
            for (int j = 1; j <= jdim; j++)
                for (int k = 1; k <= kdim; k++)
                    y(j, k, i) = -y(j, k, i);
    }

    return;

    // 999 continue
    // irr = 1
    // (handled inline via early returns above)
}

} // namespace rp3d_ns
