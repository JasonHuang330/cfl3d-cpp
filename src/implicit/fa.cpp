// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "fa.h"
#include "fa2xj.h"
#include "fa2xk.h"
#include "fa2xi.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace fa_ns {

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
    // COMMON /sklton/
    int& isklton = cmn_sklton.isklton;

    // Local arrays
    FortranArray1D<int> jbctyp(2);
    FortranArray1D<int> kbctyp(2);
    FortranArray1D<int> ibctyp(2);

    // Local scalars
    int nblc, nblcc, jf, kf, iff, js, ks, is, je, ke, ie;
    int neta;

    nfajki(1) = 0;
    nfajki(2) = 0;
    nfajki(3) = 0;

    nwfa = 0;
    nifa = 0;
    neta = 0;

    //
    // accumulate fluxes to ensure conservation - J direction
    //
    for (nblc = 1; nblc <= nblock; nblc++) {
        if (nbl == nblc) continue;  // go to 6500

        nblcc = nblcg(nblc);
        if (nblcc == nbl) {
            jf = jdimg(nblc);
            kf = kdimg(nblc);
            iff = idimg(nblc);
            js = jsg(nblc);
            ks = ksg(nblc);
            is = isg(nblc);
            je = jeg(nblc);
            ke = keg(nblc);
            ie = ieg(nblc);
            jbctyp(1) = jbcinfo(nblc, 1, 1, 1);
            jbctyp(2) = jbcinfo(nblc, 1, 1, 2);

            if ((isklton > 0 && jbctyp(1) == 21) ||
                (isklton > 0 && jbctyp(2) == 21)) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::memset(bou(nou(1), 1), ' ', 120);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 7: (1x,33hflux accumulation in j from block,i3,1x,16hfor use on block,i3)
                std::snprintf(bou(nou(1), 1), 120,
                    " flux accumulation in j from block%3d for use on block%3d",
                    nblc, nbl);
            }

            // w(lw(2,nblc)) passed as qjf — 4D array starting at that offset
            FortranArray4DRef<double> qjf_ref(&w(lw(2, nblc)), jf, kf, iff, 5);
            FortranArray1DRef<int> jbctyp_ref = jbctyp.ref();
            fa2xj_ns::fa2xj(jf, kf, iff, qjf_ref, js, ks, is, je, ke, ie,
                             jbctyp_ref, wk, iw, nwfa, nifa, neta,
                             nou, bou, nbuf, ibufdim);
        }
    }
    // 6500 continue

    nfajki(1) = neta;
    neta = 0;

    //
    // K-direction
    //
    for (nblc = 1; nblc <= nblock; nblc++) {
        if (nbl == nblc) continue;  // go to 6502

        nblcc = nblcg(nblc);
        if (nblcc == nbl) {
            jf = jdimg(nblc);
            kf = kdimg(nblc);
            iff = idimg(nblc);
            js = jsg(nblc);
            ks = ksg(nblc);
            is = isg(nblc);
            je = jeg(nblc);
            ke = keg(nblc);
            ie = ieg(nblc);
            kbctyp(1) = kbcinfo(nblc, 1, 1, 1);
            kbctyp(2) = kbcinfo(nblc, 1, 1, 2);

            if ((isklton > 0 && kbctyp(1) == 21) ||
                (isklton > 0 && kbctyp(2) == 21)) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::memset(bou(nou(1), 1), ' ', 120);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 17: (1x,33hflux accumulation in k from block,i3,1x,16hfor use on block,i3)
                std::snprintf(bou(nou(1), 1), 120,
                    " flux accumulation in k from block%3d for use on block%3d",
                    nblc, nbl);
            }

            // w(lw(3,nblc)) passed as qkf
            FortranArray4DRef<double> qkf_ref(&w(lw(3, nblc)), jf, kf, iff, 5);
            FortranArray1DRef<int> kbctyp_ref = kbctyp.ref();
            fa2xk_ns::fa2xk(jf, kf, iff, qkf_ref, js, ks, is, je, ke, ie,
                             kbctyp_ref, wk, iw, nwfa, nifa, neta,
                             nou, bou, nbuf, ibufdim);
        }
    }
    // 6502 continue

    nfajki(2) = neta;
    neta = 0;

    //
    // I-direction
    //
    for (nblc = 1; nblc <= nblock; nblc++) {
        if (nbl == nblc) continue;  // go to 6504

        nblcc = nblcg(nblc);
        if (nblcc == nbl) {
            jf = jdimg(nblc);
            kf = kdimg(nblc);
            iff = idimg(nblc);
            js = jsg(nblc);
            ks = ksg(nblc);
            is = isg(nblc);
            je = jeg(nblc);
            ke = keg(nblc);
            ie = ieg(nblc);
            ibctyp(1) = ibcinfo(nblc, 1, 1, 1);
            ibctyp(2) = ibcinfo(nblc, 1, 1, 2);

            if ((isklton > 0 && ibctyp(1) == 21) ||
                (isklton > 0 && ibctyp(2) == 21)) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::memset(bou(nou(1), 1), ' ', 120);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 71: (1x,33hflux accumulation in i from block,i3,1x,16hfor use on block,i3)
                std::snprintf(bou(nou(1), 1), 120,
                    " flux accumulation in i from block%3d for use on block%3d",
                    nblc, nbl);
            }

            // w(lw(4,nblc)) passed as qif
            FortranArray4DRef<double> qif_ref(&w(lw(4, nblc)), jf, kf, iff, 5);
            FortranArray1DRef<int> ibctyp_ref = ibctyp.ref();
            fa2xi_ns::fa2xi(jf, kf, iff, qif_ref, js, ks, is, je, ke, ie,
                             ibctyp_ref, wk, iw, nwfa, nifa, neta,
                             nou, bou, nbuf, ibufdim);
        }
    }
    // 6504 continue

    nfajki(3) = neta;
    return;
}

} // namespace fa_ns
