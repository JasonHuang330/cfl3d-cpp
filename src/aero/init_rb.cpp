// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "init_rb.h"
#include <cmath>
#include <cstdint>

namespace init_rb_ns {

// init: this function is defined in the 'init' module; the header declares it
// here as part of the init_rb module group. Body is empty (no-op delegation).
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
    (void)nbl; (void)jdim; (void)kdim; (void)idim;
    (void)q; (void)qj0; (void)qk0; (void)qi0;
    (void)tj0; (void)tk0; (void)ti0;
    (void)vol; (void)volj0; (void)volk0; (void)voli0;
    (void)nummem;
    (void)x; (void)y; (void)z;
    (void)nou; (void)bou;
    (void)nbuf; (void)ibufdim; (void)iflagprnt;
}

void init_rb(int& maxaes, int& maxbl,
             FortranArray1DRef<double> zorig,
             FortranArray2DRef<double> aesrfdat)
{
    (void)maxaes; (void)maxbl;

    // Local temp array (4x4), column-major: temp(i,j) -> temp_arr[(j-1)*4+(i-1)]
    float temp_arr[16];

    // Aliases for COMMON /rbstmt1/ flat arrays (column-major, 4x4 or 1D)
    // bmatrb(n,m) -> cmn_rbstmt1.bmatrb[(m-1)*4+(n-1)]
    // stmrb(n,m)  -> cmn_rbstmt1.stmrb[(m-1)*4+(n-1)]
    // stmirb(n,m) -> cmn_rbstmt1.stmirb[(m-1)*4+(n-1)]
    // gforcnrb(m) -> cmn_rbstmt1.gforcnrb[m-1]
    // etc.

    // Aliases for COMMON /rbstmt2/ scalars
    float& tmass   = cmn_rbstmt2.tmass;
    float& yinert  = cmn_rbstmt2.yinert;
    float& uinfrb  = cmn_rbstmt2.uinfrb;
    float& qinfrb  = cmn_rbstmt2.qinfrb;
    float& crefrb  = cmn_rbstmt2.crefrb;
    float& xtmref  = cmn_rbstmt2.xtmref;
    float& areat   = cmn_rbstmt2.areat;

    // Aliases for COMMON /fsum/
    float& sref = cmn_fsum.sref;
    float& cref = cmn_fsum.cref;
    float& xmc  = cmn_fsum.xmc;

    // Aliases for COMMON /info/
    float& dt    = cmn_info.dt;
    float& xmach = cmn_info.xmach;

    // Alias for COMMON /trim/
    float& zrg0 = cmn_trim.zrg0;

    //
    //     xcg      = xmc
    //     crefrb   = cref
    //     xtmref   = xcg - xmc
    //     areat    = sref
    //     grefl    = aesrfdat(2,1)
    //     uinfrb   = aesrfdat(3,1)
    //     qinfrb   = aesrfdat(4,1)
    //     ainf     = uinfrb/xmach
    //     timesc   = grefl/ainf
    //     dts      = dt*timesc
    //     zrg0     = zorig(1)
    //
    float xcg    = xmc;
    crefrb       = cref;
    xtmref       = xcg - xmc;
    areat        = sref;
    float grefl  = (float)aesrfdat(2, 1);
    uinfrb       = (float)aesrfdat(3, 1);
    qinfrb       = (float)aesrfdat(4, 1);
    float ainf   = (float)((double)uinfrb / (double)xmach);
    float timesc = (float)((double)grefl / (double)ainf);
    float dts    = (float)((double)dt * (double)timesc);
    zrg0         = (float)zorig(1);

    //
    //     do m=1,4
    //        do n = 1,4
    //           stmrb(n,m)   = 0.
    //           stmirb(n,m)  = 0.
    //           bmatrb(n,m)  = 0.
    //        end do
    //        gforcnrb(m)  = 0.
    //        gforcnmrb(m) = 0.
    //        gforcsrb(m)  = 0.
    //        x0rb(m)      = 0.
    //        xxnrb(m)     = x0rb(m)
    //        xsrb(m)      = xxnrb(m)
    //     end do
    //
    for (int m = 1; m <= 4; m++) {
        for (int n = 1; n <= 4; n++) {
            cmn_rbstmt1.stmrb[(m-1)*4 + (n-1)]  = 0.f;
            cmn_rbstmt1.stmirb[(m-1)*4 + (n-1)] = 0.f;
            cmn_rbstmt1.bmatrb[(m-1)*4 + (n-1)] = 0.f;
        }
        cmn_rbstmt1.gforcnrb[m-1]  = 0.f;
        cmn_rbstmt1.gforcnmrb[m-1] = 0.f;
        cmn_rbstmt1.gforcsrb[m-1]  = 0.f;
        cmn_rbstmt1.x0rb[m-1]      = 0.f;
        cmn_rbstmt1.xxnrb[m-1]     = cmn_rbstmt1.x0rb[m-1];
        cmn_rbstmt1.xsrb[m-1]      = cmn_rbstmt1.xxnrb[m-1];
    }

    //
    //     do i = 1,4
    //        stmrb(i,i)  = 1.
    //        stmirb(i,i) = dts
    //     enddo
    //
    for (int i = 1; i <= 4; i++) {
        cmn_rbstmt1.stmrb[(i-1)*4 + (i-1)]  = 1.f;
        cmn_rbstmt1.stmirb[(i-1)*4 + (i-1)] = dts;
    }

    //
    //     stmrb(1,2)   = dts
    //     stmrb(1,4)   = 0.
    //     stmrb(3,4)   = dts
    //     stmrb(2,4)   = 0.
    //     stmirb(1,2)  = 0.5*dts*dts
    //     stmirb(1,4)  = 0.
    //     stmirb(2,4)  = 0.
    //     stmirb(3,4)  = 0.5*dts*dts
    //     bmatrb(2,2)  = 1./tmass
    //     bmatrb(4,4)  = 1./yinert
    //
    // stmrb(n,m) -> [(m-1)*4+(n-1)]
    cmn_rbstmt1.stmrb[(2-1)*4 + (1-1)]  = dts;          // stmrb(1,2)
    cmn_rbstmt1.stmrb[(4-1)*4 + (1-1)]  = 0.f;          // stmrb(1,4)
    cmn_rbstmt1.stmrb[(4-1)*4 + (3-1)]  = dts;          // stmrb(3,4)
    cmn_rbstmt1.stmrb[(4-1)*4 + (2-1)]  = 0.f;          // stmrb(2,4)
    cmn_rbstmt1.stmirb[(2-1)*4 + (1-1)] = (float)(0.5 * (double)dts * (double)dts); // stmirb(1,2)
    cmn_rbstmt1.stmirb[(4-1)*4 + (1-1)] = 0.f;          // stmirb(1,4)
    cmn_rbstmt1.stmirb[(4-1)*4 + (2-1)] = 0.f;          // stmirb(2,4)
    cmn_rbstmt1.stmirb[(4-1)*4 + (3-1)] = (float)(0.5 * (double)dts * (double)dts); // stmirb(3,4)
    cmn_rbstmt1.bmatrb[(2-1)*4 + (2-1)] = (float)(1.0 / (double)tmass);   // bmatrb(2,2)
    cmn_rbstmt1.bmatrb[(4-1)*4 + (4-1)] = (float)(1.0 / (double)yinert);  // bmatrb(4,4)

    //
    //     overwrite stmi with matrix product stmi*bmat
    //     (theta*B in the references's notation)
    //
    //     do j=1,4
    //        do i=1,4
    //           temp(i,j) = 0.
    //           do k=1,4
    //              temp(i,j) = temp(i,j) + stmirb(i,k)*bmatrb(k,j)
    //           end do
    //        end do
    //     end do
    //     do j=1,4
    //        do i=1,4
    //           stmirb(i,j) = temp(i,j)
    //        end do
    //     end do
    //
    for (int j = 1; j <= 4; j++) {
        for (int i = 1; i <= 4; i++) {
            temp_arr[(j-1)*4 + (i-1)] = 0.f;
            for (int k = 1; k <= 4; k++) {
                temp_arr[(j-1)*4 + (i-1)] +=
                    cmn_rbstmt1.stmirb[(k-1)*4 + (i-1)] *
                    cmn_rbstmt1.bmatrb[(j-1)*4 + (k-1)];
            }
        }
    }
    for (int j = 1; j <= 4; j++) {
        for (int i = 1; i <= 4; i++) {
            cmn_rbstmt1.stmirb[(j-1)*4 + (i-1)] = temp_arr[(j-1)*4 + (i-1)];
        }
    }
}

} // namespace init_rb_ns
