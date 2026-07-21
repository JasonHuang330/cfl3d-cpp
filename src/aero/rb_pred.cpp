// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rb_pred.h"
#include "runtime/fortran_io.h"
#include <cstdio>

namespace rb_pred_ns {

void rb_pred(int& nt)
{
    // COMMON /rbstmt1/
    // bmatrb(4,4) column-major flat: bmatrb[(j-1)*4+(n-1)]
    // stmrb(4,4)  column-major flat: stmrb[(j-1)*4+(n-1)]
    // stmirb(4,4) column-major flat: stmirb[(j-1)*4+(n-1)]
    float* gforcnrb  = cmn_rbstmt1.gforcnrb;   // [4], 0-based
    float* gforcnmrb = cmn_rbstmt1.gforcnmrb;  // [4], 0-based
    float* gforcsrb  = cmn_rbstmt1.gforcsrb;   // [4], 0-based
    float* stmrb     = cmn_rbstmt1.stmrb;      // [16], column-major
    float* stmirb    = cmn_rbstmt1.stmirb;     // [16], column-major
    float* xsrb      = cmn_rbstmt1.xsrb;       // [4], 0-based
    float* xxnrb     = cmn_rbstmt1.xxnrb;      // [4], 0-based

    // COMMON /rbstmt2/
    float& greflrb = cmn_rbstmt2.greflrb;

    // COMMON /trim/
    float& zrg0  = cmn_trim.zrg0;
    float& zrg1  = cmn_trim.zrg1;
    float& dzdt  = cmn_trim.dzdt;
    float& alf1  = cmn_trim.alf1;
    float& thtd1 = cmn_trim.thtd1;

    // COMMON /info/
    float& dt = cmn_info.dt;

    // rigid body displacement and velocity prediction via rigid body
    // equations of motion

    // do j = 1,4
    //    gforcnmrb(j) = gforcnrb(j)
    //    gforcnrb(j)  = gforcsrb(j)
    // end do
    for (int j = 1; j <= 4; j++) {
        gforcnmrb[j-1] = gforcnrb[j-1];
        gforcnrb[j-1]  = gforcsrb[j-1];
    }

    // do n=1,4
    //    xxnrb(n) = xsrb(n)
    // end do
    for (int n = 1; n <= 4; n++) {
        xxnrb[n-1] = xsrb[n-1];
    }

    // do n=1,4
    //    xsrb(n) = 0.
    //    do j=1,4
    //       xsrb(n) = xsrb(n) + stmrb(n,j)*xxnrb(j)
    //  .             + .5*stmirb(n,j)*(3.*gforcnrb(j) - gforcnmrb(j))
    //    end do
    // end do
    for (int n = 1; n <= 4; n++) {
        xsrb[n-1] = 0.f;
        for (int j = 1; j <= 4; j++) {
            // column-major: stmrb(n,j) → stmrb[(j-1)*4 + (n-1)]
            xsrb[n-1] = xsrb[n-1]
                + stmrb[(j-1)*4 + (n-1)] * xxnrb[j-1]
                + 0.5f * stmirb[(j-1)*4 + (n-1)]
                  * (3.f * gforcnrb[j-1] - gforcnmrb[j-1]);
        }
    }

    // zrg1 = zrg0 - xsrb(1)/greflrb
    zrg1 = zrg0 - xsrb[0] / greflrb;
    // dzdt = - xsrb(2)/greflrb
    dzdt = -xsrb[1] / greflrb;
    // alf1 =  xsrb(3)
    alf1 = xsrb[2];
    // thtd1=  xsrb(4)
    thtd1 = xsrb[3];

    // write(299,31029) dt*real(nt),zrg1,xsrb(3),gforcnrb(2),gforcnrb(4)
    // 31029 format(5(1x,e16.8))
    FILE* f299 = fortran_get_unit(299);
    fprintf(f299, " %16.8e %16.8e %16.8e %16.8e %16.8e\n",
            (double)(dt * (float)nt),
            (double)zrg1,
            (double)xsrb[2],
            (double)gforcnrb[1],
            (double)gforcnrb[3]);

    return;
}

} // namespace rb_pred_ns
