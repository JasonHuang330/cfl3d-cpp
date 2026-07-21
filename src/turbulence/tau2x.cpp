// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "tau2x.h"

namespace tau2x_ns {

void tau(int& mgflag, int& nbl, int& jdim, int& kdim, int& idim,
         FortranArray2DRef<double> q, FortranArray2DRef<double> res,
         FortranArray2DRef<double> q1, FortranArray2DRef<double> qr,
         FortranArray2DRef<int> lw, FortranArray1DRef<double> w,
         FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
         int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& nblock,
         FortranArray1DRef<int> igridg, FortranArray1DRef<int> nblcg,
         FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg,
         FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg,
         FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
         FortranArray1DRef<int> iemg)
{
    for (int ibl = 1; ibl <= nbl; ibl++) {
        if (iemg(ibl) > 0) {
            int js_bl = jsg(ibl);
            int ks_bl = ksg(ibl);
            int is_bl = isg(ibl);
            int je_bl = jeg(ibl);
            int ke_bl = keg(ibl);
            int ie_bl = ieg(ibl);
            int kode  = mgflag;

            int jqr = je_bl - js_bl + 1;
            int kqr = ke_bl - ks_bl + 1;
            int iqr = ie_bl - is_bl;

            FortranArray4DRef<double> res_bl(&res(1, ibl),
                                             jdim, kdim, idim - 1, 5);
            FortranArray4DRef<double> qr_bl(&qr(1, ibl),
                                            jqr, kqr, iqr, 5);

            tau2x(jdim, kdim, idim, res_bl, qr_bl,
                  js_bl, ks_bl, is_bl, je_bl, ke_bl, ie_bl, kode);
        }
    }
}

void tau2x(int& jdim, int& kdim, int& idim,
           FortranArray4DRef<double> res, FortranArray4DRef<double> qr,
           int& js, int& ks, int& is, int& je, int& ke, int& ie, int& kode)
{
    int n, k, i, j, kk, ii, jj;

    //      put restricted  r from finer embedded mesh into coarser mesh
    //
    //      jdim,kdim,idim  coarser mesh indices
    //      js,ks,is        coarser mesh starting indices
    //      je,ke,ie        coarser mesh ending indices

    if (kode >= 2) {
        for (n = 1; n <= 5; n++) {
            kk = 0;
            for (k = ks; k <= ke - 1; k++) {
                kk = kk + 1;
                ii = 0;
                for (i = is; i <= ie - 1; i++) {
                    ii = ii + 1;
                    jj = 0;
                    for (j = js; j <= je - 1; j++) {
                        jj = jj + 1;
                        qr(jj, kk, ii, n) = qr(jj, kk, ii, n) - res(j, k, i, n);
                    }
                }
            }
        }
    }

    for (n = 1; n <= 5; n++) {
        kk = 0;
        for (k = ks; k <= ke - 1; k++) {
            kk = kk + 1;
            ii = 0;
            for (i = is; i <= ie - 1; i++) {
                ii = ii + 1;
                jj = 0;
                for (j = js; j <= je - 1; j++) {
                    jj = jj + 1;
                    res(j, k, i, n) = res(j, k, i, n) + qr(jj, kk, ii, n);
                }
            }
        }
    }
}

} // namespace tau2x_ns
