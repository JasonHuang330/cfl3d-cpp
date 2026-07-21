// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rb_corr.h"

namespace rb_corr_ns {

void rb_corr(FortranArray2DRef<double> aesrfdat, double& cnt, double& cmt, int& maxaes)
{
    // COMMON /rbstmt1/
    // bmatrb(4,4) column-major flat: bmatrb[(j-1)*4+(n-1)]
    // gforcnrb(4), gforcnmrb(4), gforcsrb(4) — 0-based
    // stmrb(4,4), stmirb(4,4) column-major flat
    // xsrb(4), xxnrb(4), xnmrb(4), x0rb(4) — 0-based
    float* gforcnrb  = cmn_rbstmt1.gforcnrb;
    float* gforcsrb  = cmn_rbstmt1.gforcsrb;
    float* stmrb     = cmn_rbstmt1.stmrb;
    float* stmirb    = cmn_rbstmt1.stmirb;
    float* xsrb      = cmn_rbstmt1.xsrb;
    float* xxnrb     = cmn_rbstmt1.xxnrb;
    float* xnmrb     = cmn_rbstmt1.xnmrb;

    // COMMON /rbstmt2/
    float& tmass    = cmn_rbstmt2.tmass;
    float& qinfrb   = cmn_rbstmt2.qinfrb;
    float& greflrb  = cmn_rbstmt2.greflrb;
    float& gaccel   = cmn_rbstmt2.gaccel;
    float& crefrb   = cmn_rbstmt2.crefrb;
    float& xtmref   = cmn_rbstmt2.xtmref;
    float& areat    = cmn_rbstmt2.areat;

    // COMMON /trim/
    float& cnw = cmn_trim.cnw;
    float& cmy = cmn_trim.cmy;

    // rigid body displacement and velocity prediction via rigid body
    // equations of motion

    double grefl = aesrfdat(2, 1);
    cnw = (float)cnt;
    cmy = (float)cmt;

    // gforcsrb is 0-based (Fortran index n → C index n-1)
    gforcsrb[0] = 0.f;
    gforcsrb[1] = (float)((double)tmass * (double)gaccel
                          - (double)qinfrb * (double)areat
                            * grefl * grefl * cnt);
    gforcsrb[2] = 0.f;
    gforcsrb[3] = (float)((double)qinfrb * (double)areat
                          * grefl * grefl * (double)greflrb
                          * ((double)crefrb * cmt + (double)xtmref * cnt));

    // do n=1,4: xnmrb(n) = xsrb(n)
    for (int n = 1; n <= 4; n++) {
        xnmrb[n - 1] = xsrb[n - 1];
    }

    // do n=1,4: xsrb(n) = sum over j=1,4 of stmrb(n,j)*xxnrb(j)
    //                    + 0.5*stmirb(n,j)*(gforcsrb(j)+gforcnrb(j))
    // stmrb(n,j) column-major: stmrb[(j-1)*4 + (n-1)]
    for (int n = 1; n <= 4; n++) {
        double sum = 0.0;
        for (int j = 1; j <= 4; j++) {
            sum += (double)stmrb[(j - 1) * 4 + (n - 1)] * (double)xxnrb[j - 1]
                 + 0.5 * (double)stmirb[(j - 1) * 4 + (n - 1)]
                       * ((double)gforcsrb[j - 1] + (double)gforcnrb[j - 1]);
        }
        xsrb[n - 1] = (float)sum;
    }
}

} // namespace rb_corr_ns
