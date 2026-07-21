// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "init_mast.h"
#include "init.h"
#include <cmath>

namespace init_mast_ns {

void init(int& nbl, int& jdim, int& kdim, int& idim,
          FortranArray4DRef<double> q,
          FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
          FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0,
          FortranArray3DRef<double> vol,
          FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0,
          int& nummem,
          FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
          FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
          int& nbuf, int& ibufdim, int& iflagprnt)
{
    (void)nbl; (void)jdim; (void)kdim; (void)idim;
    (void)q;
    (void)qj0; (void)qk0; (void)qi0;
    (void)tj0; (void)tk0; (void)ti0;
    (void)vol;
    (void)volj0; (void)volk0; (void)voli0;
    (void)nummem;
    (void)x; (void)y; (void)z;
    (void)nou; (void)bou;
    (void)nbuf; (void)ibufdim; (void)iflagprnt;
}

void init_mast()
{
    // COMMON /fluid/
    float& gamma  = cmn_fluid.gamma;
    float& gm1    = cmn_fluid.gm1;

    // COMMON /info/
    float& xmach  = cmn_info.xmach;
    float& alpha  = cmn_info.alpha;
    float& beta   = cmn_info.beta;

    // COMMON /ivals/
    float& p0     = cmn_ivals.p0;
    float& rho0   = cmn_ivals.rho0;
    float& c0     = cmn_ivals.c0;
    float& u0     = cmn_ivals.u0;
    float& v0     = cmn_ivals.v0;
    float& w0     = cmn_ivals.w0;
    float& et0    = cmn_ivals.et0;
    float& h0     = cmn_ivals.h0;
    float& pt0    = cmn_ivals.pt0;
    // qiv is float[5], 1-based in Fortran → 0-based in C
    float* qiv    = cmn_ivals.qiv;

    rho0   = 1.0f;
    c0     = 1.0f;
    p0     = rho0 * c0 * c0 / gamma;

    //   The wind axis system follows NASA SP-3070 (1972), with the exception that
    //   positive beta is in the opposite direction
    u0     = xmach * std::cos((double)alpha) * std::cos((double)beta);
    w0     = xmach * std::sin((double)alpha) * std::cos((double)beta);
    v0     = -xmach * std::sin((double)beta);

    float ei0 = p0 / ((gamma - 1.0f) * rho0);
    et0    = rho0 * (ei0 + 0.5f * (u0*u0 + v0*v0 + w0*w0));
    h0     = (et0 + p0) / rho0;
    pt0    = p0 * std::pow((double)(1.0f + 0.5f * gm1 * xmach * xmach),
                           (double)(gamma / gm1));

    qiv[0] = rho0;   // qiv(1)
    qiv[1] = u0;     // qiv(2)
    qiv[2] = v0;     // qiv(3)
    qiv[3] = w0;     // qiv(4)
    qiv[4] = p0;     // qiv(5)
}

} // namespace init_mast_ns
