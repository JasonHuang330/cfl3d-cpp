// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "initnonin.h"
#include "init.h"

namespace initnonin_ns {

void init(int& nbl, int& jdim, int& kdim, int& idim,
          FortranArray4DRef<double> q,
          FortranArray4DRef<double> qj0,
          FortranArray4DRef<double> qk0,
          FortranArray4DRef<double> qi0,
          FortranArray4DRef<double> tj0,
          FortranArray4DRef<double> tk0,
          FortranArray4DRef<double> ti0,
          FortranArray3DRef<double> vol,
          FortranArray3DRef<double> volj0,
          FortranArray3DRef<double> volk0,
          FortranArray3DRef<double> voli0,
          int& nummem,
          FortranArray3DRef<double> x,
          FortranArray3DRef<double> y,
          FortranArray3DRef<double> z,
          FortranArray1DRef<int> nou,
          FortranArray2DRef<char[120]> bou,
          int& nbuf,
          int& ibufdim,
          int& iflagprnt)
{
    init_ns::init(nbl, jdim, kdim, idim, q, qj0, qk0, qi0,
                  tj0, tk0, ti0, vol, volj0, volk0, voli0,
                  nummem, x, y, z, nou, bou, nbuf, ibufdim, iflagprnt);
}

void initnonin(int& nbl, int& jdim, int& kdim, int& idim,
               FortranArray4DRef<double> q,
               FortranArray4DRef<double> qj0,
               FortranArray4DRef<double> qk0,
               FortranArray4DRef<double> qi0,
               FortranArray3DRef<double> vol,
               FortranArray3DRef<double> volj0,
               FortranArray3DRef<double> volk0,
               FortranArray3DRef<double> voli0,
               FortranArray3DRef<double> x,
               FortranArray3DRef<double> y,
               FortranArray3DRef<double> z)
{
    // local aliases for rotation rates from /noninertial/ common block
    double wx = (double)cmn_noninertial.xrotrate;
    double wy = (double)cmn_noninertial.yrotrate;
    double wz = (double)cmn_noninertial.zrotrate;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    for (int j = 1; j <= jdim1; j++) {
        for (int k = 1; k <= kdim1; k++) {
            for (int i = 1; i <= idim1; i++) {

                // compute cell centers
                double cx = 0.125 * (
                    x(j  , k  , i  ) + x(j  , k  , i+1) +
                    x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                    x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                    x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

                double cy = 0.125 * (
                    y(j  , k  , i  ) + y(j  , k  , i+1) +
                    y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                    y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                    y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

                double cz = 0.125 * (
                    z(j  , k  , i  ) + z(j  , k  , i+1) +
                    z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                    z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                    z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

                cx = cx - (double)cmn_noninertial.xcentrot;
                cy = cy - (double)cmn_noninertial.ycentrot;
                cz = cz - (double)cmn_noninertial.zcentrot;

                // increase velocity with rigid body rotation component:
                // Uinf + r x omega is the same as Uinf - omega x r
                q(j, k, i, 2) = q(j, k, i, 2) + ( cy * wz - cz * wy );
                q(j, k, i, 3) = q(j, k, i, 3) + ( cz * wx - cx * wz );
                q(j, k, i, 4) = q(j, k, i, 4) + ( cx * wy - cy * wx );
            }
        }
    }
}

} // namespace initnonin_ns
