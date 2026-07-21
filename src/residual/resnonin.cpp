// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "resnonin.h"
#include "runtime/fortran_io.h"
#include <algorithm>
#include <cstdio>
#include <cstring>

namespace resnonin_ns {

void resnonin(int& nbl, int& jdim, int& kdim, int& idim,
              FortranArray2DRef<double> q,
              FortranArray1DRef<double> x,
              FortranArray1DRef<double> y,
              FortranArray1DRef<double> z,
              FortranArray2DRef<double> sj,
              FortranArray2DRef<double> sk,
              FortranArray2DRef<double> si,
              FortranArray1DRef<double> vol,
              FortranArray2DRef<double> res,
              FortranArray1DRef<int> nou,
              FortranArray2DRef<char[120]> bou,
              int& nbuf, int& ibufdim)
{
    // COMMON block aliases
    float& xcentrot  = cmn_noninertial.xcentrot;
    float& ycentrot  = cmn_noninertial.ycentrot;
    float& zcentrot  = cmn_noninertial.zcentrot;
    float& xrotrate  = cmn_noninertial.xrotrate;
    float& yrotrate  = cmn_noninertial.yrotrate;
    float& zrotrate  = cmn_noninertial.zrotrate;
    int32_t& isklton = cmn_sklton.isklton;

    // qiv(2), qiv(3), qiv(4) — 1-based Fortran → 0-based C array
    float* qiv = cmn_ivals.qiv;  // qiv[0..4], Fortran qiv(1..5)

    if (isklton > 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
                      "   adding NONINERTIAL source terms");
    }

    // rename xrotrate, yrotrate, zrotrate -> wx, wy, wz to save typing
    double wx = (double)xrotrate;
    double wy = (double)yrotrate;
    double wz = (double)zrotrate;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    // compute Omega x Uinf (constant with respect to solution or q)
    // qiv(2) = v0, qiv(3) = w0, qiv(4) = u0 (1-based Fortran indexing)
    double OmegaxUx = ( wy * (double)qiv[3] - wz * (double)qiv[2] );
    double OmegaxUy = ( wz * (double)qiv[1] - wx * (double)qiv[3] );
    double OmegaxUz = ( wx * (double)qiv[2] - wy * (double)qiv[1] );

    for (int i = 1; i <= idim1; i++) {
        for (int k = 1; k <= kdim1; k++) {
            for (int j = 1; j <= jdim1; j++) {

                // Find the indices corresponding to the 8 corners of a volume
                int ind000 = (i-1) * jdim * kdim + (k-1) * jdim + (j-1) + 1;
                int ind001 = (i-1) * jdim * kdim + (k-1) * jdim + (j  ) + 1;
                int ind010 = (i-1) * jdim * kdim + (k  ) * jdim + (j-1) + 1;
                int ind011 = (i-1) * jdim * kdim + (k  ) * jdim + (j  ) + 1;
                int ind100 = (i  ) * jdim * kdim + (k-1) * jdim + (j-1) + 1;
                int ind101 = (i  ) * jdim * kdim + (k-1) * jdim + (j  ) + 1;
                int ind110 = (i  ) * jdim * kdim + (k  ) * jdim + (j-1) + 1;
                int ind111 = (i  ) * jdim * kdim + (k  ) * jdim + (j  ) + 1;

                // find the radius vector (rx, ry, rz) from the center of rotation
                // to the center of the volume
                double rx = 0.125 * (
                    x(ind000) +
                    x(ind001) +
                    x(ind010) +
                    x(ind011) +
                    x(ind100) +
                    x(ind101) +
                    x(ind110) +
                    x(ind111) ) - (double)xcentrot;

                double ry = 0.125 * (
                    y(ind000) +
                    y(ind001) +
                    y(ind010) +
                    y(ind011) +
                    y(ind100) +
                    y(ind101) +
                    y(ind110) +
                    y(ind111) ) - (double)ycentrot;

                double rz = 0.125 * (
                    z(ind000) +
                    z(ind001) +
                    z(ind010) +
                    z(ind011) +
                    z(ind100) +
                    z(ind101) +
                    z(ind110) +
                    z(ind111) ) - (double)zcentrot;

                // compute centripetal pseudo-acceleration vector
                // (constant with respect to q) = -w x (w x r)
                double centaccx = wy * (wx*ry - wy*rx) - wz * (wz*rx - wx*rz);
                double centaccy = wz * (wy*rz - wz*ry) - wx * (wx*ry - wy*rx);
                double centaccz = wx * (wz*rx - wx*rz) - wy * (wy*rz - wz*ry);

                // Coriolis part of pseudo-acceleration
                double corx = 2. * ( wy * q(ind000,4) - wz * q(ind000,3) );
                double cory = 2. * ( wz * q(ind000,2) - wx * q(ind000,4) );
                double corz = 2. * ( wx * q(ind000,3) - wy * q(ind000,2) );

                double totx = -OmegaxUx + centaccx + corx;
                double toty = -OmegaxUy + centaccy + cory;
                double totz = -OmegaxUz + centaccz + corz;

                // increment residual (in conserved vars.)
                res(ind000,2) = res(ind000,2) + q(ind000,1) * totx * vol(ind000);
                res(ind000,3) = res(ind000,3) + q(ind000,1) * toty * vol(ind000);
                res(ind000,4) = res(ind000,4) + q(ind000,1) * totz * vol(ind000);
                res(ind000,5) = res(ind000,5) + q(ind000,1) * vol(ind000) * (
                    totx * q(ind000,2) +
                    toty * q(ind000,3) +
                    totz * q(ind000,4) );

            }
        }
    }

    return;
}

} // namespace resnonin_ns
