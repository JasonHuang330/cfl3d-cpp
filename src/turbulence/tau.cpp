// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "tau.h"
#include "tau2x.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <algorithm>

namespace tau_ns {

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
    // Access COMMON blocks
    int& levt    = cmn_mgrd.levt;
    int& kode    = cmn_mgrd.kode;
    int& mode    = cmn_mgrd.mode;
    int& ncyc    = cmn_mgrd.ncyc;
    int& mtt     = cmn_mgrd.mtt;
    int& icyc    = cmn_mgrd.icyc;
    int& level   = cmn_mgrd.level;
    int& lglobal = cmn_mgrd.lglobal;
    int& isklton = cmn_sklton.isklton;

    // Suppress unused variable warnings
    (void)levt; (void)ncyc; (void)mtt; (void)icyc;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;
    (void)jdim1; (void)kdim1;

    //
    // store q1 for later use in determining delta q
    // residual correction
    //

    if (kode >= 2 && mgflag >= 1) {
        if (level >= lglobal && mgflag < 2) goto label_1001;

        //
        // store q1=q   kode.ge.2
        //
        {
            int n    = jdim * kdim;
            int nplq = std::min(idim1, 999000 / n);
            int npl  = nplq;
            for (int i = 1; i <= idim1; i += nplq) {
                if (i + npl - 1 > idim1) npl = idim1 - i + 1;
                int nnpl = n * npl - jdim - 1;
                int ist  = (i - 1) * n;
                for (int l = 1; l <= 5; l++) {
                    for (int izz = 1; izz <= nnpl; izz++) {
                        q1(izz + ist, l) = q(izz + ist, l);
                    }
                }
            }
        }
    }
label_1001:;

    if (mode != 0) {
        //
        // qr=qr-res       mode.ne.0    kode.ge.2
        //
        // res=res+qr      mode.ne.0
        //

        if (level < lglobal) {
            //
            // residual correction from global grids
            //

            if (isklton > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " ");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " adding residual correction for block%6d from finer block%6d; kode=%3d",
                    nbl, nbl - 1, kode);
            }

            {
                int n    = jdim * kdim;
                int nplq = std::min(idim1, 999000 / n);
                int npl  = nplq;
                for (int i = 1; i <= idim1; i += nplq) {
                    if (i + npl - 1 > idim1) npl = idim1 - i + 1;
                    int nnpl = n * npl - jdim - 1;
                    int ist  = (i - 1) * n;

                    if (kode >= 2) {
                        for (int l = 1; l <= 5; l++) {
                            for (int izz = 1; izz <= nnpl; izz++) {
                                qr(izz + ist, l) = qr(izz + ist, l) - res(izz + ist, l);
                            }
                        }
                    }

                    for (int l = 1; l <= 5; l++) {
                        for (int izz = 1; izz <= nnpl; izz++) {
                            res(izz + ist, l) = res(izz + ist, l) + qr(izz + ist, l);
                        }
                    }
                }
            }

        } else {
            //
            // residual correction from embedded grids
            //

            for (int nblc = 1; nblc <= nblock; nblc++) {
                int igrid = igridg(nblc);
                if (nblc == nbl || iemg(igrid) == 0) continue;

                int nblcc = nblcg(nblc);
                if (nblcc == nbl) {
                    int jsc  = jsg(nblc);
                    int ksc  = ksg(nblc);
                    int isc  = isg(nblc);
                    int jec  = jeg(nblc);
                    int kec  = keg(nblc);
                    int iec  = ieg(nblc);
                    int lqrc = lw(17, nblc);

                    if (isklton > 0) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " ");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " adding residual correction for block%6d from embedded block%6d; kode=%3d",
                            nbl, nblc, kode);
                    }

                    // Build a 4D view of res for tau2x (sequence association)
                    // tau2x expects FortranArray4DRef<double> res(jdim,kdim,idim-1,5)
                    // and FortranArray4DRef<double> qr (the work array slice w(lqrc))
                    FortranArray4D<double> res4d(jdim, kdim, idim - 1, 5);
                    // Copy res into res4d (sequence association: same linear order)
                    {
                        int total = jdim * kdim * (idim - 1) * 5;
                        double* src = &res(1, 1);
                        double* dst = &res4d(1, 1, 1, 1);
                        for (int ii = 0; ii < total; ii++) dst[ii] = src[ii];
                    }

                    // Build 4D view of w(lqrc) for qr argument
                    // w(lqrc) is the start of the embedded qr data
                    // tau2x expects FortranArray4DRef<double> qr(jdim,kdim,idim-1,5)
                    FortranArray4D<double> qr4d(jdim, kdim, idim - 1, 5);
                    {
                        int total = jdim * kdim * (idim - 1) * 5;
                        double* src = &w(lqrc);
                        double* dst = &qr4d(1, 1, 1, 1);
                        for (int ii = 0; ii < total; ii++) dst[ii] = src[ii];
                    }

                    FortranArray4DRef<double> res4d_ref = res4d.ref();
                    FortranArray4DRef<double> qr4d_ref  = qr4d.ref();

                    tau2x_ns::tau2x(jdim, kdim, idim, res4d_ref, qr4d_ref,
                                    jsc, ksc, isc, jec, kec, iec, kode);

                    // Copy res4d back to res (sequence association)
                    {
                        int total = jdim * kdim * (idim - 1) * 5;
                        double* src = &res4d(1, 1, 1, 1);
                        double* dst = &res(1, 1);
                        for (int ii = 0; ii < total; ii++) dst[ii] = src[ii];
                    }
                }
            }
        }
    }

    return;
}

} // namespace tau_ns
