// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fa2xj.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <algorithm>

namespace fa2xj_ns {

void fa(int& nbl, FortranArray2DRef<int> lw, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<int> iw, int& ifamax, int& nwfa, int& nifa, FortranArray1DRef<int> nfajki, int& maxbl, int& maxseg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, int& nblock, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        // qjf is stored at w(lw(4,ibl)) with dimensions (kf,if-1,5,4)
        // jf=jdim, kf=kdim, if=idim
        // qjf(kf, if-1, 5, 4) = qjf(kdim, idim-1, 5, 4)
        FortranArray4DRef<double> qjf(
            reinterpret_cast<double*>(&w(lw(4, ibl))),
            kdim, idim - 1, 5, 4);

        int js_loc = jsg(ibl);
        int ks_loc = ksg(ibl);
        int is_loc = isg(ibl);
        int je_loc = jeg(ibl);
        int ke_loc = keg(ibl);
        int ie_loc = ieg(ibl);

        FortranArray1D<int> jbctyp_loc(2);
        jbctyp_loc(1) = 0;
        jbctyp_loc(2) = 0;

        int nsegj = nfajki(ibl);
        for (int iseg = 1; iseg <= nsegj; iseg++) {
            if (jbcinfo(ibl, iseg, 1, 1) == 21) {
                int nface = jbcinfo(ibl, iseg, 1, 8);
                if (nface == 1) {
                    jbctyp_loc(1) = 21;
                } else {
                    jbctyp_loc(2) = 21;
                }
            }
        }

        if (jbctyp_loc(1) == 21 || jbctyp_loc(2) == 21) {
            FortranArray1DRef<int> jbctyp_ref = jbctyp_loc.ref();
            fa2xj(jdim, kdim, idim, qjf,
                  js_loc, ks_loc, is_loc,
                  je_loc, ke_loc, ie_loc,
                  jbctyp_ref, w, iw, nwfa, nifa, ifamax,
                  nou, bou, nbuf, ibufdim);
        }
    }
}

void fa2xj(int& jf, int& kf, int& if_, FortranArray4DRef<double> qjf,
           int& js, int& ks, int& is, int& je, int& ke, int& ie,
           FortranArray1DRef<int> jbctyp, FortranArray1DRef<double> w,
           FortranArray1DRef<int> iw, int& nw, int& niw, int& neta,
           FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim)
{
    int jfm1 = jf - 1;
    int kfm1 = kf - 1;
    int ifm1 = if_ - 1;
    int nsi  = ifm1 / (ie - is);

    if (jbctyp(1) == 21) {
        // left boundary
        int mm = 2;
        iw(niw + 1) = js;
        iw(niw + 2) = ks;
        iw(niw + 3) = is;
        iw(niw + 4) = js;
        iw(niw + 5) = ke;
        iw(niw + 6) = ie;
        iw(niw + 7) = nw + 1;
        niw = niw + 7;
        neta = neta + 1;

        if (nsi == 2) {
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= ifm1; i += 2) {
                    for (int k = 1; k <= kfm1; k += 2) {
                        nw = nw + 1;
                        w(nw) = qjf(k,   i,   l, mm)
                              + qjf(k+1, i,   l, mm)
                              + qjf(k,   i+1, l, mm)
                              + qjf(k+1, i+1, l, mm);
                    }
                }
            }
        } else {
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= ifm1; i++) {
                    for (int k = 1; k <= kfm1; k += 2) {
                        nw = nw + 1;
                        w(nw) = qjf(k,   i, l, mm)
                              + qjf(k+1, i, l, mm);
                    }
                }
            }
        }

        if (cmn_sklton.isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "  summing flux at j=js  niw= %d  nw= %d", niw, nw);
        }
    }

    if (jbctyp(2) == 21) {
        // right boundary
        int mm = 4;
        iw(niw + 1) = je;
        iw(niw + 2) = ks;
        iw(niw + 3) = is;
        iw(niw + 4) = je;
        iw(niw + 5) = ke;
        iw(niw + 6) = ie;
        iw(niw + 7) = nw + 1;
        niw = niw + 7;
        neta = neta + 1;

        if (nsi == 2) {
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= ifm1; i += 2) {
                    for (int k = 1; k <= kfm1; k += 2) {
                        nw = nw + 1;
                        w(nw) = qjf(k,   i,   l, mm)
                              + qjf(k+1, i,   l, mm)
                              + qjf(k,   i+1, l, mm)
                              + qjf(k+1, i+1, l, mm);
                    }
                }
            }
        } else {
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= ifm1; i++) {
                    for (int k = 1; k <= kfm1; k += 2) {
                        nw = nw + 1;
                        w(nw) = qjf(k,   i, l, mm)
                              + qjf(k+1, i, l, mm);
                    }
                }
            }
        }

        if (cmn_sklton.isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "  summing flux at j=je  niw= %d  nw= %d", niw, nw);
        }
    }
}

} // namespace fa2xj_ns
