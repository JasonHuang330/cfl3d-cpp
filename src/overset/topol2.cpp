// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "topol2.h"
#include "topol.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace topol2_ns {

void topol(int& jdim, int& kdim, int& nsub, FortranArray1DRef<int> jjmax, FortranArray1DRef<int> kkmax, int& lmax, int& l, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> xmid, FortranArray3DRef<double> ymid, FortranArray3DRef<double> zmid, FortranArray3DRef<double> xmide, FortranArray3DRef<double> ymide, FortranArray3DRef<double> zmide, int& limit, double& xc, double& yc, double& zc, double& xie, double& eta, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, int& ifit, int& itmax, int& igap, int& iok, FortranArray1DRef<int> lout, int& ic0, int& itoss0, int& jto, int& kto, int& iself, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid)
{
    topol_ns::topol(jdim, kdim, nsub, jjmax, kkmax, lmax, l, x, y, z,
                    xmid, ymid, zmid, xmide, ymide, zmide, limit, xc, yc, zc,
                    xie, eta, jimage, kimage, ifit, itmax, igap, iok, lout,
                    ic0, itoss0, jto, kto, iself, xif1, xif2, etf1, etf2,
                    nou, bou, nbuf, ibufdim, myid);
}

void topol2(int& jdim1, int& kdim1, int& msub1, FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1, int& lmax1, double& xie, double& eta, FortranArray3DRef<double> sarc, int& limit0, double& sc, int& jp, int& kp, int& lp, int& lsrch, int& itmax, double& xiet, double& etat, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    // COMMON /tol/
    float& epsc    = cmn_tol.epsc;
    // COMMON /filenam/
    char* grdmov   = cmn_filenam.grdmov;

    int jmax   = jjmax1(lp);
    int kmax   = kkmax1(lp);
    int js     = xif1(lp);
    int je     = xif2(lp);
    int ks     = etf1(lp);
    int ke     = etf2(lp);
    int intern = 0;
    int limit  = limit0;

    // label 6901 loop
    while (true) {
        intern = intern + 1;
        if (intern > itmax) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                " program terminated in dynamic patching routines - see file %-60.60s",
                grdmov);
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120, " stopping...intern > itmax ");
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            // termn8 does not return, but for safety:
            return;
        }

        // bi-linear fit in target cell jp,kp in "from" block lp
        double dy2 = sarc(jp+1, kp,   lp) - sarc(jp, kp, lp);
        double dy3 = sarc(jp+1, kp+1, lp) - sarc(jp, kp, lp);
        double dy4 = sarc(jp,   kp+1, lp) - sarc(jp, kp, lp);
        double a1  = sarc(jp, kp, lp);
        double a2  = dy2;
        double a3  = dy4;
        double a4  = dy3 - a2 - a3;

        if (lsrch == 1) {
            // for fixed eta
            xie = (sc - a1 - a3 * eta) / (a2 + a4 * eta);

            // check to ensure point is inside cell
            int imiss = 0;
            if ((float)xie < -(float)epsc || (float)xie > 1.0f + (float)epsc) {
                imiss = 1;
            }

            if (imiss == 0) {
                // go to 6902
                break;
            }

            // try new target cell
            int jpc = jp;
            int jinc;
            if ((float)xie >= 0.0f) jinc = (int)std::abs(xie);
            if ((float)xie <  0.0f) jinc = (int)std::abs(xie - 1.0);
            if (limit > 3 && intern > 5)  limit = 3;
            if (limit > 1 && intern > 10) limit = 1;
            jinc = std::min(jinc, limit);
            if ((float)xie > 1.0f) {
                jpc = jp + jinc;
            } else if ((float)xie < 0.0f) {
                jpc = jp - jinc;
            }
            jp = std::min(jpc, je - 1);
            jp = std::max(js, jp);
            // go to 6901
            continue;

        } else if (lsrch == 2) {
            // for fixed xie
            eta = (sc - a1 - a2 * xie) / (a3 + a4 * xie);

            int imiss = 0;
            if ((float)eta < -(float)epsc || (float)eta > 1.0f + (float)epsc) {
                imiss = 1;
            }
            if (imiss == 0) {
                // go to 6902
                break;
            }

            // try new target cell
            int kpc = kp;
            int kpl = kp;
            int kinc;
            if ((float)eta >= 0.0f) kinc = (int)std::abs(eta);
            if ((float)eta <  0.0f) kinc = (int)std::abs(eta - 1.0);
            if (limit > 3 && intern > 5)  limit = 3;
            if (limit > 1 && intern > 10) limit = 1;
            kinc = std::min(kinc, limit);
            if ((float)eta > 1.0f) {
                kpc = kp + kinc;
            } else if ((float)eta < 0.0f) {
                kpc = kp - kinc;
            }
            kp = std::min(kpc, ke - 1);
            kp = std::max(ks, kp);
            // go to 6901
            continue;
        }
        // if neither lsrch==1 nor lsrch==2, break out (shouldn't happen)
        break;
    }
    // label 6902
    xiet = jp + xie;
    etat = kp + eta;
}

} // namespace topol2_ns
