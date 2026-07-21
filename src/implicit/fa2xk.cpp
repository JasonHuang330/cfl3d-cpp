// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fa2xk.h"
#include "runtime/fortran_io.h"
#include <algorithm>
#include <cstdio>

namespace fa2xk_ns {

void fa(int& nbl, FortranArray2DRef<int> lw, FortranArray1DRef<double> w, int& mgwk,
        FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<int> iw, int& ifamax,
        int& nwfa, int& nifa, FortranArray1DRef<int> nfajki, int& maxbl, int& maxseg,
        FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg,
        FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg,
        FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
        FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
        FortranArray4DRef<int> ibcinfo, int& nblock, FortranArray1DRef<int> nblcg,
        FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        // qk0 is at lw(3,ibl), dimensions (jdim, idim-1, 5, 4)
        FortranArray4DRef<double> qk0(&w(lw(3, ibl)), jdim, idim - 1, 5, 4);

        int js_loc = jsg(ibl);
        int ks_loc = ksg(ibl);
        int is_loc = isg(ibl);
        int je_loc = jeg(ibl);
        int ke_loc = keg(ibl);
        int ie_loc = ieg(ibl);

        FortranArray1D<int> kbctyp_loc(2);
        kbctyp_loc(1) = 0;
        kbctyp_loc(2) = 0;

        int nsegk = nfajki(ibl);
        for (int iseg = 1; iseg <= nsegk; iseg++) {
            if (kbcinfo(ibl, iseg, 1, 1) == 21) {
                int nface = kbcinfo(ibl, iseg, 1, 8);
                if (nface == 5) {
                    kbctyp_loc(1) = 21;
                } else {
                    kbctyp_loc(2) = 21;
                }
            }
        }

        if (kbctyp_loc(1) == 21 || kbctyp_loc(2) == 21) {
            FortranArray1DRef<int> kbctyp_ref = kbctyp_loc.ref();
            fa2xk(jdim, kdim, idim, qk0, js_loc, ks_loc, is_loc,
                  je_loc, ke_loc, ie_loc, kbctyp_ref, w, iw,
                  nwfa, nifa, ifamax, nou, bou, nbuf, ibufdim);
        }
    }
}

void fa2xk(int& jf, int& kf, int& if_, FortranArray4DRef<double> qkf,
           int& js, int& ks, int& is, int& je, int& ke, int& ie,
           FortranArray1DRef<int> kbctyp, FortranArray1DRef<double> w,
           FortranArray1DRef<int> iw, int& nw, int& niw, int& neta,
           FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim)
{
    int& isklton = cmn_sklton.isklton;

    int jfm1 = jf - 1;
    int kfm1 = kf - 1;
    int ifm1 = if_ - 1;
    int nsi  = ifm1 / (ie - is);

    if (kbctyp(1) == 21) {
        // left boundary
        int mm        = 2;
        iw(niw + 1)   = js;
        iw(niw + 2)   = ks;
        iw(niw + 3)   = is;
        iw(niw + 4)   = je;
        iw(niw + 5)   = ks;
        iw(niw + 6)   = ie;
        iw(niw + 7)   = nw + 1;
        niw           = niw + 7;
        neta          = neta + 1;

        if (nsi == 2) {
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= ifm1; i += 2) {
                    for (int j = 1; j <= jfm1; j += 2) {
                        nw    = nw + 1;
                        w(nw) = qkf(j,   i,   l, mm)
                              + qkf(j+1, i,   l, mm)
                              + qkf(j,   i+1, l, mm)
                              + qkf(j+1, i+1, l, mm);
                    }
                }
            }
        } else {
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= ifm1; i++) {
                    for (int j = 1; j <= jfm1; j += 2) {
                        nw    = nw + 1;
                        w(nw) = qkf(j,   i, l, mm)
                              + qkf(j+1, i, l, mm);
                    }
                }
            }
        }

        if (isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                          "  summing flux at k=ks  niw= %d  nw= %d", niw, nw);
        }
    }

    if (kbctyp(2) == 21) {
        // right boundary
        int mm        = 4;
        iw(niw + 1)   = js;
        iw(niw + 2)   = ke;
        iw(niw + 3)   = is;
        iw(niw + 4)   = je;
        iw(niw + 5)   = ke;
        iw(niw + 6)   = ie;
        iw(niw + 7)   = nw + 1;
        niw           = niw + 7;
        neta          = neta + 1;

        if (nsi == 2) {
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= ifm1; i += 2) {
                    for (int j = 1; j <= jfm1; j += 2) {
                        nw    = nw + 1;
                        w(nw) = qkf(j,   i,   l, mm)
                              + qkf(j+1, i,   l, mm)
                              + qkf(j,   i+1, l, mm)
                              + qkf(j+1, i+1, l, mm);
                    }
                }
            }
        } else {
            for (int l = 1; l <= 5; l++) {
                for (int i = 1; i <= ifm1; i++) {
                    for (int j = 1; j <= jfm1; j += 2) {
                        nw    = nw + 1;
                        w(nw) = qkf(j,   i, l, mm)
                              + qkf(j+1, i, l, mm);
                    }
                }
            }
        }

        if (isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                          "  summing flux at k=ke  niw= %d  nw= %d", niw, nw);
        }
    }
}

} // namespace fa2xk_ns
