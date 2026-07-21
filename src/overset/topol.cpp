// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "topol.h"
#include "dsmin.h"
#include "xe.h"
#include "xe2.h"
#include <cmath>
#include <algorithm>

namespace topol_ns {

void topol(int& jdim, int& kdim, int& nsub, FortranArray1DRef<int> jjmax, FortranArray1DRef<int> kkmax, int& lmax, int& l, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> xmid, FortranArray3DRef<double> ymid, FortranArray3DRef<double> zmid, FortranArray3DRef<double> xmide, FortranArray3DRef<double> ymide, FortranArray3DRef<double> zmide, int& limit, double& xc, double& yc, double& zc, double& xie, double& eta, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, int& ifit, int& itmax, int& igap, int& iok, FortranArray1DRef<int> lout, int& ic0, int& itoss0, int& jto, int& kto, int& iself, FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2, FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid)
{
    // local arrays
    FortranArray1D<int> jfroz(itmax);
    FortranArray1D<int> kfroz(itmax);

    // COMMON /tol/
    float& epsc0 = cmn_tol.epsc0;

    int    idum1, idum2, idum3, idum4;
    double dum1, dum2, dum3;

    idum1 = 0;
    idum2 = 0;
    idum3 = 0;
    idum4 = 0;
    dum1  = 0.;
    dum2  = 0.;
    dum3  = 0.;

    for (int ll = 1; ll <= lmax; ll++)
        lout(ll) = 0;

    int iatmpt = 0;
    int ichk   = 0;
    int ifroze = 0;
    int ifroz  = 0;
    int ihuge  = 0;
    int jp     = (int)xie;
    int kp     = (int)eta;
    int lsav   = 1;
    double xiesav = xie;
    double etasav = eta;
    int jpsav  = jp;
    int kpsav  = kp;
    int jpc    = jp;
    int kpc    = kp;

    // to start search, use solution from last cell, if one exists; otherwise
    // start by searching for minimum distance point over all "from" blocks
    if ((float)xie < 1.f || (float)eta < 1.f) {
        int lflag_init = 999;
        dsmin_ns::dsmin(jdim, kdim, nsub, jjmax, kkmax, lmax, x, y, z,
                        xc, yc, zc, jp, kp, l, lout, lflag_init,
                        xif1, xif2, etf1, etf2);
    }

    // label 999
    label_999:

    iatmpt = iatmpt + 1;

    // all "from" blocks have been searched if iatmpt > lmax
    if (iatmpt > lmax) goto label_1000;

    {
        int jmax = jjmax(l);
        int kmax = kkmax(l);
        int js   = xif1(l);
        int je   = xif2(l);
        int ks   = etf1(l);
        int ke   = etf2(l);

        // avoid starting outside specified search range
        jp = std::min(jp, je - 1);
        kp = std::min(kp, ke - 1);
        jp = std::max(jp, js);
        kp = std::max(kp, ks);

        // do 5555 intern=1,itmax
        for (int intern = 1; intern <= itmax; intern++) {

            jfroz(intern) = jp;
            kfroz(intern) = kp;

            // find local xie, eta via Newton iteration in current target cell jp,kp
            int imiss = 0;

            if (itoss0 == 0) {
                // call general routine which first determines the best direction
                // for inversion, then call Newton iteration routine
                xe_ns::xe(jdim, kdim, nsub, l, x, y, z, xmid, ymid, zmid,
                           xmide, ymide, zmide, jp, kp, xc, yc, zc,
                           xie, eta, imiss, ifit, ic0,
                           nou, bou, nbuf, ibufdim, myid);
            } else {
                // best direction for inversion is known apriori
                double x1 = x(jp,   kp,   l);
                double y1 = y(jp,   kp,   l);
                double z1 = z(jp,   kp,   l);
                double x2 = x(jp+1, kp,   l);
                double y2 = y(jp+1, kp,   l);
                double z2 = z(jp+1, kp,   l);
                double x4 = x(jp,   kp+1, l);
                double y4 = y(jp,   kp+1, l);
                double z4 = z(jp,   kp+1, l);
                double x3 = x(jp+1, kp+1, l);
                double y3 = y(jp+1, kp+1, l);
                double z3 = z(jp+1, kp+1, l);
                double x5 = xmid(jp,   kp,   l);
                double y5 = ymid(jp,   kp,   l);
                double z5 = zmid(jp,   kp,   l);
                double x6 = xmid(jp,   kp+1, l);
                double y6 = ymid(jp,   kp+1, l);
                double z6 = zmid(jp,   kp+1, l);
                double x7 = xmide(jp,   kp,   l);
                double y7 = ymide(jp,   kp,   l);
                double z7 = zmide(jp,   kp,   l);
                double x8 = xmide(jp+1, kp,   l);
                double y8 = ymide(jp+1, kp,   l);
                double z8 = zmide(jp+1, kp,   l);

                if (itoss0 == 1) {
                    // use only y and z equations
                    xe2_ns::xe2(y1, y2, y3, y4, y5, y6, y7, y8, yc,
                                z1, z2, z3, z4, z5, z6, z7, z8, zc,
                                xie, eta, imiss, ifit,
                                nou, bou, nbuf, ibufdim, myid);
                }
                if (itoss0 == 2) {
                    // use only x and z equations
                    xe2_ns::xe2(x1, x2, x3, x4, x5, x6, x7, x8, xc,
                                z1, z2, z3, z4, z5, z6, z7, z8, zc,
                                xie, eta, imiss, ifit,
                                nou, bou, nbuf, ibufdim, myid);
                }
                if (itoss0 == 3) {
                    // use only x and y equations
                    xe2_ns::xe2(x1, x2, x3, x4, x5, x6, x7, x8, xc,
                                y1, y2, y3, y4, y5, y6, y7, y8, yc,
                                xie, eta, imiss, ifit,
                                nou, bou, nbuf, ibufdim, myid);
                }
            }

            // check to make sure that search did not find current "to" cell as the
            // target (a possibility for a grid which communicates with itself along
            // a branch cut). If so, move to a cell on the other side of the branch cut
            if (iself > 0) {
                if (imiss == 0 && (jto == jp - 1 && kto == kp - 1)) {
                    // determine in which direction the branch cut lies
                    int ibrdir = 0;
                    jpc = jpc - 1;
                    jpc = std::max(1, jpc);
                    if (jpc != jimage(l, jpc, kp)) ibrdir = 1;
                    jpc = jpc + 1;
                    jpc = std::min(jpc, jmax - 1);
                    if (jpc != jimage(l, jpc, kp)) ibrdir = 2;
                    kpc = kpc - 1;
                    kpc = std::max(1, kpc);
                    if (kpc != kimage(l, jp, kpc)) ibrdir = 3;
                    kpc = kpc + 1;
                    kpc = std::min(kpc, kmax - 1);
                    if (kpc != kimage(l, jp, kpc)) ibrdir = 4;

                    // set cell index to trigger a branch cut jump
                    if (ibrdir == 1) jp = jp - 1;
                    if (ibrdir == 2) jp = jp + 1;
                    if (ibrdir == 3) kp = kp - 1;
                    if (ibrdir == 4) kp = kp + 1;

                    jp = std::min(jp, jmax - 1);
                    kp = std::min(kp, kmax - 1);
                    jp = std::max(1, jp);
                    kp = std::max(1, kp);
                    jp = jimage(l, jp, kp);
                    kp = kimage(l, jp, kp);
                    imiss = 1;
                    continue; // go to 5555
                }
            }

            // current target cell correct if imiss = 0
            if (imiss == 0) goto label_5556;

            // update current guess for target cell based on result of Newton
            // iteration, with max allowable change set by limit
            int jinc = 0, kinc = 0;
            if ((float)xie >= 0.f) jinc = (int)std::abs(xie);
            if ((float)xie <  0.f) jinc = (int)std::abs(xie - 1.0);
            if ((float)eta >= 0.f) kinc = (int)std::abs(eta);
            if ((float)eta <  0.f) kinc = (int)std::abs(eta - 1.0);

            jinc = std::min(jinc, limit);
            kinc = std::min(kinc, limit);

            if ((float)xie > 1.0f) {
                jp = jp + jinc;
            } else if ((float)xie < 0.f) {
                jp = jp - jinc;
            }
            if ((float)eta > 1.0f) {
                kp = kp + kinc;
            } else if ((float)eta < 0.f) {
                kp = kp - kinc;
            }

            // keep within bounds of (expanded) "from" block
            jp = std::min(jp, jmax - 1);
            kp = std::min(kp, kmax - 1);
            jp = std::max(1, jp);
            kp = std::max(1, kp);

            jmax = jjmax(l);
            kmax = kkmax(l);
            js   = xif1(l);
            je   = xif2(l);
            ks   = etf1(l);
            ke   = etf2(l);

            // avoid cells outside specified search range
            jp = std::min(jp, je - 1);
            kp = std::min(kp, ke - 1);
            jp = std::max(jp, js);
            kp = std::max(kp, ks);

            // account for any branch cuts
            jpc = jimage(l, jp, kp);
            kpc = kimage(l, jp, kp);
            jp = jpc;
            kp = kpc;

            // search routine off track if local xie or eta become huge
            {
                float huge_val = 1.e+5f;
                if (std::abs((float)xie) > huge_val || std::abs((float)eta) > huge_val) {
                    ihuge = ihuge + 1;
                    if (ihuge > 1) goto label_1000;
                    int lflag_neg = -999;
                    dsmin_ns::dsmin(jdim, kdim, nsub, jjmax, kkmax, lmax, x, y, z,
                                    xc, yc, zc, jp, kp, l, lout, lflag_neg,
                                    xif1, xif2, etf1, etf2);
                    continue; // go to 5555
                }
            }

            // check for frozen convergence
            ifroz = 0;
            for (int ii = 1; ii <= intern; ii++) {
                int int_idx = intern - ii + 1;
                if (jp == jfroz(int_idx) && kp == kfroz(int_idx)) ifroz = 1;
            }
            if (ifroz == 1) {
                ifroze = ifroze + 1;
                if (ifroze > 1) goto label_1000;
                int lflag_neg = -999;
                dsmin_ns::dsmin(jdim, kdim, nsub, jjmax, kkmax, lmax, x, y, z,
                                xc, yc, zc, jp, kp, l, lout, lflag_neg,
                                xif1, xif2, etf1, etf2);
            }

            // 5555 continue (end of do loop body)
        } // end do 5555 intern=1,itmax
    } // end block after label_999

    label_1000:

    // search routine has been unsuccessful in the current "from" block
    if (iatmpt < lmax) {
        // find new starting cell by searching for minimum distance over
        // all "from" blocks not yet searched
        lout(l) = 1;
        int lflag_pos = 999;
        dsmin_ns::dsmin(jdim, kdim, nsub, jjmax, kkmax, lmax, x, y, z,
                        xc, yc, zc, jp, kp, l, lout, lflag_pos,
                        xif1, xif2, etf1, etf2);
        ifroze = 0;
        ihuge  = 0;
        goto label_999;
    } else {
        if (ichk == 0) {
            iok = 0;
            return;
        } else {
            l   = lsav;
            jp  = jpsav;
            kp  = kpsav;
            xie = xiesav;
            eta = etasav;
        }
        if (ichk == 1) {
            // point really does lie in the expanded cell of the "from" block
            goto label_5557;
        }
        if (ichk > 1) {
            // point lies in the expanded zone of two or more grids
            igap = 1;
            goto label_5557;
        }
    }

    label_5556:

    // search routine has been successful
    if (lmax > 1) {
        int jmax = jjmax(l);
        int kmax = kkmax(l);

        if (ic0 > 0) {
            int j = jp;
            int k = kp;
            double x1c = 0.25 * (x(j,   k,   l) + x(j+1, k,   l)
                               + x(j+1, k+1, l) + x(j,   k+1, l));
            double y1c = 0.25 * (y(j,   k,   l) + y(j+1, k,   l)
                               + y(j+1, k+1, l) + y(j,   k+1, l));
            double z1c = 0.25 * (z(j,   k,   l) + z(j+1, k,   l)
                               + z(j+1, k+1, l) + z(j,   k+1, l));
            if (std::abs((float)(x1c - xc)) > (float)epsc0 ||
                std::abs((float)(y1c - yc)) > (float)epsc0 ||
                std::abs((float)(z1c - zc)) > (float)epsc0) {
                ichk = ichk + 1;
                if (ichk == 1) {
                    lsav   = l;
                    jpsav  = jp;
                    kpsav  = kp;
                    xiesav = xie;
                    etasav = eta;
                }
                lout(l) = 1;
                int lflag_pos = 999;
                dsmin_ns::dsmin(jdim, kdim, nsub, jjmax, kkmax, lmax, x, y, z,
                                xc, yc, zc, jp, kp, l, lout, lflag_pos,
                                xif1, xif2, etf1, etf2);
                ifroze = 0;
                ihuge  = 0;
                goto label_999;
            }
        }

        if (jp == 1 || jp == jmax - 1 || kp == 1 || kp == kmax - 1) {
            ichk = ichk + 1;
            if (ichk == 1) {
                lsav   = l;
                jpsav  = jp;
                kpsav  = kp;
                xiesav = xie;
                etasav = eta;
            }
            lout(l) = 1;
            int lflag_pos = 999;
            dsmin_ns::dsmin(jdim, kdim, nsub, jjmax, kkmax, lmax, x, y, z,
                            xc, yc, zc, jp, kp, l, lout, lflag_pos,
                            xif1, xif2, etf1, etf2);
            ifroze = 0;
            ihuge  = 0;
            goto label_999;
        }
    }

    label_5557:

    // completed search for "to" cell center; convert to global xie, eta
    iok = 1;
    xie = xie + jp;
    eta = eta + kp;

    return;
} // end topol

} // namespace topol_ns
