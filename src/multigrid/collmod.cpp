// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "collmod.h"

namespace collmod_ns {

void collmod(FortranArray5DRef<double> xmdj, FortranArray5DRef<double> xmdk, FortranArray5DRef<double> xmdi, FortranArray5DRef<double> xmdjc, FortranArray5DRef<double> xmdkc, FortranArray5DRef<double> xmdic, int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2, int& nm, int& nmds, int& iaes, int& maxaes)
{
    int ii, iinc, kk, jj, ll;
    int i, j, k;

    //
    // j-surfaces
    //
    ii   = 0;
    iinc = 2;
    if (idim == 2) iinc = 1;
    for (i = 1; i <= idim; i += iinc) {
        ii = ii + 1;
        kk = 0;
        for (k = 1; k <= kdim; k += 2) {
            kk = kk + 1;
            for (ll = 1; ll <= 6; ll++) {
                xmdjc(kk,ii,ll,nm,iaes) = xmdj(k,i,ll,nm,iaes);
                xmdjc(kk,ii,ll,nm,iaes) = xmdj(k,i,ll,nm,iaes);
                xmdjc(kk,ii,ll,nm,iaes) = xmdj(k,i,ll,nm,iaes);
            }
        }
    }
    //
    // k-surfaces
    //
    ii   = 0;
    iinc = 2;
    if (idim == 2) iinc = 1;
    for (i = 1; i <= idim; i += iinc) {
        ii = ii + 1;
        jj = 0;
        for (j = 1; j <= jdim; j += 2) {
            jj = jj + 1;
            for (ll = 1; ll <= 6; ll++) {
                xmdkc(jj,ii,ll,nm,iaes) = xmdk(j,i,ll,nm,iaes);
                xmdkc(jj,ii,ll,nm,iaes) = xmdk(j,i,ll,nm,iaes);
                xmdkc(jj,ii,ll,nm,iaes) = xmdk(j,i,ll,nm,iaes);
            }
        }
    }
    //
    // i-surfaces
    //
    kk   = 0;
    for (k = 1; k <= kdim; k += 2) {
        kk = kk + 1;
        jj = 0;
        for (j = 1; j <= jdim; j += 2) {
            jj = jj + 1;
            for (ll = 1; ll <= 6; ll++) {
                xmdic(jj,kk,ll,nm,iaes) = xmdi(j,k,ll,nm,iaes);
                xmdic(jj,kk,ll,nm,iaes) = xmdi(j,k,ll,nm,iaes);
                xmdic(jj,kk,ll,nm,iaes) = xmdi(j,k,ll,nm,iaes);
            }
        }
    }
    //
    return;
}

} // namespace collmod_ns
