// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fa2xi.h"
#include "fa.h"
#include <algorithm>
#include <cstdio>

namespace fa2xi_ns {

void fa(int& nbl, FortranArray2DRef<int> lw, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, FortranArray1DRef<int> iw, int& ifamax, int& nwfa, int& nifa, FortranArray1DRef<int> nfajki, int& maxbl, int& maxseg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray4DRef<int> ibcinfo, int& nblock, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    fa_ns::fa(nbl, lw, w, mgwk, wk, nwork, iw, ifamax, nwfa, nifa, nfajki, maxbl, maxseg, jdimg, kdimg, idimg, jsg, ksg, isg, jeg, keg, ieg, jbcinfo, kbcinfo, ibcinfo, nblock, nblcg, nou, bou, nbuf, ibufdim);
}

void fa2xi(int& jf, int& kf, int& if_, FortranArray4DRef<double> qif, int& js, int& ks, int& is, int& je, int& ke, int& ie, FortranArray1DRef<int> ibctyp, FortranArray1DRef<double> w, FortranArray1DRef<int> iw, int& nw, int& niw, int& neta, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim)
{
    int& isklton = cmn_sklton.isklton;

    int jfm1 = jf - 1;
    int kfm1 = kf - 1;
    int ifm1 = if_ - 1;

    if (ibctyp(1) == 21) {
//
//     left boundary
//
        int mm        = 2;
        iw(niw+1) = js;
        iw(niw+2) = ks;
                        iw(niw+3) = is;
        iw(niw+4) = je;
        iw(niw+5) = ke;
                        iw(niw+6) = is;
        iw(niw+7) = nw+1;
        niw       = niw+7;
        neta      = neta+1;

        for (int l = 1; l <= 5; l++) {
        for (int k = 1; k <= kfm1; k += 2) {
        for (int j = 1; j <= jfm1; j += 2) {
            nw    = nw+1;
            w(nw) = qif(j,k,l,mm)
                   +qif(j+1,k,l,mm)
                   +qif(j,k+1,l,mm)
                   +qif(j+1,k+1,l,mm);
        }}}
        if (isklton == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "  summing flux at i=is  niw= %d  nw= %d", niw, nw);
        }
    }

    if (ibctyp(2) == 21) {
//
//     right boundary
//
        int mm        = 4;
        iw(niw+1) = js;
        iw(niw+2) = ks;
                        iw(niw+3) = ie;
        iw(niw+4) = je;
        iw(niw+5) = ke;
                        iw(niw+6) = ie;
        iw(niw+7) = nw+1;
        niw       = niw+7;
        neta      = neta+1;

        for (int l = 1; l <= 5; l++) {
        for (int k = 1; k <= kfm1; k += 2) {
        for (int j = 1; j <= jfm1; j += 2) {
            nw    = nw+1;
            w(nw) = qif(j,k,l,mm)
                   +qif(j+1,k,l,mm)
                   +qif(j,k+1,l,mm)
                   +qif(j+1,k+1,l,mm);
        }}}
        if (isklton == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "  summing flux at i=ie  niw= %d  nw= %d", niw, nw);
        }
    }
}

} // namespace fa2xi_ns
