// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "avgint.h"
#include <algorithm>

namespace avgint_ns {

void avgint(FortranArray1DRef<double> xiec, FortranArray1DRef<double> etac,
            FortranArray1DRef<int> nblkc, int& nptc,
            FortranArray1DRef<double> xief, FortranArray1DRef<double> etaf,
            FortranArray1DRef<int> nblkf, int& nptf,
            int& j1c, int& j2c, int& k1c, int& k2c,
            int& j1f, int& j2f, int& k1f, int& k2f)
{
    int kc, kf, kfp, jc, jf, jfp;
    int llc, llf, llfjp, llfkp, llfjpkp;
    double xiefavg, etafavg;
    int jfhat, kfhat;
    double xiefhat, etafhat;
    int jchat, kchat;
    double xiechat, etachat;

    for (kc = k1c; kc <= k2c - 1; kc++) {
        kf  = (kc - 1) * 2 + 1;
        kfp = kf + 1;
        kfp = std::min(kfp, k2f - 1);
        for (jc = j1c; jc <= j2c - 1; jc++) {
            jf  = (jc - 1) * 2 + 1;
            jfp = jf + 1;
            jfp = std::min(jfp, j2f - 1);
            llc     = (j2c - j1c) * (kc - k1c) + (jc - j1c + 1);
            llf     = (j2f - j1f) * (kf  - k1f) + (jf  - j1f + 1);
            llfjp   = (j2f - j1f) * (kf  - k1f) + (jfp - j1f + 1);
            llfkp   = (j2f - j1f) * (kfp - k1f) + (jf  - j1f + 1);
            llfjpkp = (j2f - j1f) * (kfp - k1f) + (jfp - j1f + 1);

            // all finer points must lie in same block for averaging;
            // otherwise, use only data from lower corner point
            if (nblkf(llfjp)   != nblkf(llf)) llfjp   = llf;
            if (nblkf(llfkp)   != nblkf(llf)) llfkp   = llf;
            if (nblkf(llfjpkp) != nblkf(llf)) llfjpkp = llf;

            xiefavg = 0.25 * (xief(llf)   + xief(llfjp)
                            + xief(llfkp) + xief(llfjpkp));
            etafavg = 0.25 * (etaf(llf)   + etaf(llfjp)
                            + etaf(llfkp) + etaf(llfjpkp));
            jfhat   = (int)xiefavg;
            xiefhat = xiefavg - jfhat;
            kfhat   = (int)etafavg;
            etafhat = etafavg - kfhat;
            jchat   = (jfhat + 1) / 2;
            kchat   = (kfhat + 1) / 2;
            if (j1f == 1 && j2f == 2) {
                // 2d case
                xiechat = xiefhat;
            } else {
                if (jfhat / 2 * 2 == jfhat) {
                    xiechat = 0.5 * (1. + xiefhat);
                } else {
                    xiechat = 0.5 * xiefhat;
                }
            }
            if (k1f == 1 && k2f == 2) {
                // 2d case
                etachat = etafhat;
            } else {
                if (kfhat / 2 * 2 == kfhat) {
                    etachat = 0.5 * (1. + etafhat);
                } else {
                    etachat = 0.5 * etafhat;
                }
            }
            xiec(llc)  = jchat + xiechat;
            etac(llc)  = kchat + etachat;
            nblkc(llc) = nblkf(llf);
        }
    }

    return;
}

} // namespace avgint_ns
