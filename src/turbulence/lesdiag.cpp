// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "lesdiag.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace lesdiag_ns {

void lesdiag(int& myid, int& jdim, int& kdim, int& idim,
             FortranArray4DRef<double> q, FortranArray4DRef<double> ux,
             FortranArray3DRef<double> vist3d, FortranArray3DRef<double> vol,
             FortranArray4DRef<double> si, FortranArray4DRef<double> sj,
             FortranArray4DRef<double> sk,
             FortranArray3DRef<double> vor, FortranArray3DRef<double> smin,
             FortranArray3DRef<double> xjb, FortranArray4DRef<double> tursav,
             FortranArray3DRef<double> xkb, FortranArray3DRef<double> blnum,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& nbl, int& nummem,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y,
             FortranArray3DRef<double> z)
{
    // COMMON block aliases
    float& gamma         = cmn_fluid.gamma;
    float& gm1           = cmn_fluid.gm1;
    float& cbar          = cmn_fluid2.cbar;
    float& pr            = cmn_fluid2.pr;
    float& prt           = cmn_fluid2.prt;
    float& reue          = cmn_reyue.reue;
    float& tinf          = cmn_reyue.tinf;
    float& time          = cmn_unst.time;
    float& xmach         = cmn_info.xmach;
    int32_t& ntt         = cmn_info.ntt;
    int32_t& icyc        = cmn_mgrd.icyc;
    int32_t& ncyc        = cmn_mgrd.ncyc;
    int32_t& isklton     = cmn_sklton.isklton;
    int32_t& les_model   = cmn_lesinfo.les_model;
    int32_t& les_wallscale = cmn_lesinfo.les_wallscale;
    float& cs_smagorinsky = cmn_lesinfo.cs_smagorinsky;
    float& cs_wale        = cmn_lesinfo.cs_wale;
    float& cs_vreman      = cmn_lesinfo.cs_vreman;

    // Note: iexp and i2d are used but never declared in the original Fortran
    // (implicit integer, undefined). Translate faithfully.
    int iexp = 0;
    int i2d  = 0;

    // Note: (10.**(-iexp+1) is machine zero)
    double xminn = std::pow(10.0, (double)(-iexp + 1));

    double c2b  = (double)cbar / (double)tinf;
    double c2bp = c2b + 1.0;
    double re   = (double)reue / (double)xmach;

    if (isklton > 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "     Computing LES-type turbulent viscosity, block=%5d", nbl);

        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            "     WARNING: ivisc=25 still under development... use at your own risk!");

        if (les_model == 0) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "      no subgrid model");
        } else if (les_model == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "      Smag model, Cs=%10.5f", (double)cs_smagorinsky);
        } else if (les_model == 2) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "      WALE model, Cs=%10.5f", (double)cs_wale);
        } else if (les_model == 3) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "      Vreman model, Cs=%10.5f", (double)cs_vreman);
        }

        if (les_wallscale == 0 || les_model != 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "      no wall damping");
        } else {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "      van Driest wall damping on");
        }

        if (cs_smagorinsky != 0.0f && les_model != 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "     WARNING: cs_smagorinsky has been set but not used (les_model .ne. 1)");
        }
        if (cs_wale != 0.0f && les_model != 2) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "     WARNING: cs_wale has been set but not used (les_model .ne. 2)");
        }
        if (cs_vreman != 0.0f && les_model != 3) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "     WARNING: cs_vreman has been set but not used (les_model .ne. 3)");
        }
    }

    // acquire mu_t if employing an SGS model
    if (les_model == 0) {
        for (int i = 1; i <= idim - 1; i++) {
            for (int j = 1; j <= jdim - 1; j++) {
                for (int k = 1; k <= kdim - 1; k++) {
                    vist3d(j, k, i) = 0.0;
                }
            }
        }
    } else {
        for (int i = 1; i <= idim - 1; i++) {
            for (int j = 1; j <= jdim - 1; j++) {
                for (int k = 1; k <= kdim - 1; k++) {
                    // Set ypls = big number if nearest body is NOT in current block
                    double ypls;
                    int nblb = (int)(blnum(j, k, i) + 0.1);
                    if (nblb != nbl) {
                        ypls = 1000.0;
                    } else {
                        int jbb = (int)(xjb(j, k, i) + 0.1);
                        int kbb = (int)(xkb(j, k, i) + 0.1);
                        int ibb = (int)(tursav(j, k, i, 2) + 0.1);
                        if (jbb == jdim) jbb = jbb - 1;
                        if (kbb == kdim) kbb = kbb - 1;
                        if (ibb == idim) ibb = ibb - 1;
                        ibb = std::max(ibb, 1);
                        double tt  = (double)gamma * q(jbb, kbb, ibb, 5) / q(jbb, kbb, ibb, 1);
                        double wnu = c2bp * tt * std::sqrt(tt) / (c2b + tt) / q(jbb, kbb, ibb, 1);
                        double utau = std::sqrt(wnu * vor(jbb, kbb, ibb) / (q(jbb, kbb, ibb, 1) * re));
                        ypls = re * q(jbb, kbb, ibb, 1) * utau * std::abs(smin(j, k, i)) / wnu;
                    }

                    double deltaj = 2.0 * vol(j, k, i) / (sj(j, k, i, 4) + sj(j + 1, k, i, 4));
                    double deltak = 2.0 * vol(j, k, i) / (sk(j, k, i, 4) + sk(j, k + 1, i, 4));
                    double deltai = 2.0 * vol(j, k, i) / (si(j, k, i, 4) + si(j, k, i + 1, 4));
                    double delta;
                    if (i2d == 0) {
                        delta = std::pow(deltai * deltaj * deltak, 0.333333);
                    } else {
                        delta = std::sqrt(deltaj * deltak);
                    }

                    if (les_model == 1) {
                        // Standard Smagorinsky model:
                        // Scale delta by van Driest function
                        if (les_wallscale == 1) {
                            double vandriest = 1.0 - std::exp(-ypls / 25.0);
                            delta = delta * vandriest;
                        }
                        double s11 = ux(j, k, i, 1);
                        double s22 = ux(j, k, i, 5);
                        double s33 = ux(j, k, i, 9);
                        double s12 = 0.5 * (ux(j, k, i, 2) + ux(j, k, i, 4));
                        double s13 = 0.5 * (ux(j, k, i, 3) + ux(j, k, i, 7));
                        double s23 = 0.5 * (ux(j, k, i, 6) + ux(j, k, i, 8));
                        double xis = s11 * s11 + s22 * s22 + s33 * s33 +
                                     2.0 * s12 * s12 + 2.0 * s13 * s13 + 2.0 * s23 * s23;
                        vist3d(j, k, i) = q(j, k, i, 1) *
                            ((double)cs_smagorinsky * delta) * ((double)cs_smagorinsky * delta) *
                            std::sqrt(2.0 * xis) * (double)reue / (double)xmach;
                    } else if (les_model == 2) {
                        // WALE model (Flow, Turb, & Combust 62:183-200 1999)
                        double s11 = ux(j, k, i, 1);
                        double s22 = ux(j, k, i, 5);
                        double s33 = ux(j, k, i, 9);
                        double s12 = 0.5 * (ux(j, k, i, 2) + ux(j, k, i, 4));
                        double s13 = 0.5 * (ux(j, k, i, 3) + ux(j, k, i, 7));
                        double s23 = 0.5 * (ux(j, k, i, 6) + ux(j, k, i, 8));
                        double xis = s11 * s11 + s22 * s22 + s33 * s33 +
                                     2.0 * s12 * s12 + 2.0 * s13 * s13 + 2.0 * s23 * s23;
                        double g11 = ux(j, k, i, 1) * ux(j, k, i, 1) +
                                     ux(j, k, i, 2) * ux(j, k, i, 4) +
                                     ux(j, k, i, 3) * ux(j, k, i, 7);
                        double g22 = ux(j, k, i, 4) * ux(j, k, i, 2) +
                                     ux(j, k, i, 5) * ux(j, k, i, 5) +
                                     ux(j, k, i, 6) * ux(j, k, i, 8);
                        double g33 = ux(j, k, i, 7) * ux(j, k, i, 3) +
                                     ux(j, k, i, 8) * ux(j, k, i, 6) +
                                     ux(j, k, i, 9) * ux(j, k, i, 9);
                        double g12 = 0.5 * (ux(j, k, i, 1) * ux(j, k, i, 2) +
                                            ux(j, k, i, 2) * ux(j, k, i, 5) +
                                            ux(j, k, i, 3) * ux(j, k, i, 8) +
                                            ux(j, k, i, 4) * ux(j, k, i, 1) +
                                            ux(j, k, i, 5) * ux(j, k, i, 4) +
                                            ux(j, k, i, 6) * ux(j, k, i, 7));
                        double g13 = 0.5 * (ux(j, k, i, 1) * ux(j, k, i, 3) +
                                            ux(j, k, i, 2) * ux(j, k, i, 6) +
                                            ux(j, k, i, 3) * ux(j, k, i, 9) +
                                            ux(j, k, i, 7) * ux(j, k, i, 1) +
                                            ux(j, k, i, 8) * ux(j, k, i, 4) +
                                            ux(j, k, i, 9) * ux(j, k, i, 7));
                        double g23 = 0.5 * (ux(j, k, i, 4) * ux(j, k, i, 3) +
                                            ux(j, k, i, 5) * ux(j, k, i, 6) +
                                            ux(j, k, i, 6) * ux(j, k, i, 9) +
                                            ux(j, k, i, 7) * ux(j, k, i, 2) +
                                            ux(j, k, i, 8) * ux(j, k, i, 5) +
                                            ux(j, k, i, 9) * ux(j, k, i, 8));
                        // Note: Fortran updates g11/g22/g33 in sequence (each uses
                        // the ORIGINAL g11+g22+g33 sum, not the updated values)
                        // Fortran line 229: g11=g11-0.3333333*(g11+g22+g33)
                        // Fortran line 230: g22=g22-0.3333333*(g11+g22+g33)  <- uses UPDATED g11
                        // Fortran line 231: g33=g33-0.3333333*(g11+g22+g33)  <- uses UPDATED g11,g22
                        // Translate faithfully (sequential updates):
                        g11 = g11 - 0.3333333 * (g11 + g22 + g33);
                        g22 = g22 - 0.3333333 * (g11 + g22 + g33);
                        g33 = g33 - 0.3333333 * (g11 + g22 + g33);
                        double g_gamma = g11 * g11 + g22 * g22 + g33 * g33 +
                                         2.0 * g12 * g12 + 2.0 * g13 * g13 + 2.0 * g23 * g23;
                        double denom = std::pow(xis, 2.5) + std::pow(g_gamma, 1.25);
                        denom = std::max(denom, xminn);
                        vist3d(j, k, i) = q(j, k, i, 1) *
                            ((double)cs_wale * delta) * ((double)cs_wale * delta) *
                            std::pow(g_gamma, 1.5) / denom * (double)reue / (double)xmach;
                    } else if (les_model == 2) {
                        // Vreman model (Phys Fluids 16(10):3670-3681 2004)
                        // Note: original Fortran has "else if (les_model .eq. 2)" here too
                        // (dead code in original - faithfully translated)
                        double b11 = delta * delta * (ux(j, k, i, 1) * ux(j, k, i, 1) +
                                                      ux(j, k, i, 4) * ux(j, k, i, 4) +
                                                      ux(j, k, i, 7) * ux(j, k, i, 7));
                        double b22 = delta * delta * (ux(j, k, i, 2) * ux(j, k, i, 2) +
                                                      ux(j, k, i, 5) * ux(j, k, i, 5) +
                                                      ux(j, k, i, 8) * ux(j, k, i, 8));
                        double b33 = delta * delta * (ux(j, k, i, 3) * ux(j, k, i, 3) +
                                                      ux(j, k, i, 6) * ux(j, k, i, 6) +
                                                      ux(j, k, i, 9) * ux(j, k, i, 9));
                        double b12 = delta * delta * (ux(j, k, i, 1) * ux(j, k, i, 2) +
                                                      ux(j, k, i, 4) * ux(j, k, i, 5) +
                                                      ux(j, k, i, 7) * ux(j, k, i, 8));
                        double b13 = delta * delta * (ux(j, k, i, 1) * ux(j, k, i, 3) +
                                                      ux(j, k, i, 4) * ux(j, k, i, 6) +
                                                      ux(j, k, i, 7) * ux(j, k, i, 9));
                        double b23 = delta * delta * (ux(j, k, i, 2) * ux(j, k, i, 3) +
                                                      ux(j, k, i, 5) * ux(j, k, i, 6) +
                                                      ux(j, k, i, 8) * ux(j, k, i, 9));
                        double b_beta = b11 * b22 - b12 * b12 + b11 * b33 - b13 * b13 +
                                        b22 * b33 - b23 * b23;
                        double denom = ux(j, k, i, 1) * ux(j, k, i, 1) +
                                       ux(j, k, i, 5) * ux(j, k, i, 5) +
                                       ux(j, k, i, 9) * ux(j, k, i, 9) +
                                       2.0 * ux(j, k, i, 2) * ux(j, k, i, 4) +
                                       2.0 * ux(j, k, i, 3) * ux(j, k, i, 7) +
                                       2.0 * ux(j, k, i, 6) * ux(j, k, i, 8);
                        denom = std::max(denom, xminn);
                        vist3d(j, k, i) = q(j, k, i, 1) * (double)cs_vreman *
                            std::sqrt(b_beta / denom) * (double)reue / (double)xmach;
                    }
                }
            }
        }
    }

    // obtain LES diagnostic info
    if (icyc == ncyc) {
        c2b  = (double)cbar / (double)tinf;
        c2bp = c2b + 1.0;
        double eddytime = 1.0 / (double)xmach;
        double sumvel2   = 0.0;
        double sumux2    = 0.0;
        double sumvy2    = 0.0;
        double sumwz2    = 0.0;
        double sumux3    = 0.0;
        double sumvy3    = 0.0;
        double sumwz3    = 0.0;
        double sumu2     = 0.0;
        double sumv2     = 0.0;
        double sumw2     = 0.0;
        double sumuv     = 0.0;
        double sumuw     = 0.0;
        double sumvw     = 0.0;
        double sumenergy = 0.0;
        double fnu       = 0.0;
        double sumxis    = 0.0;
        double sumwis    = 0.0;
        double rho       = 0.0;
        int n = 0;

        for (int i = 1; i <= idim - 1; i++) {
            for (int j = 1; j <= jdim - 1; j++) {
                for (int k = 1; k <= kdim - 1; k++) {
                    n = n + 1;
                    sumvel2 = sumvel2 + q(j, k, i, 2) * q(j, k, i, 2) +
                                       q(j, k, i, 3) * q(j, k, i, 3) +
                                       q(j, k, i, 4) * q(j, k, i, 4);
                    sumux2 = sumux2 + (ux(j, k, i, 1) * ux(j, k, i, 1));
                    sumvy2 = sumvy2 + (ux(j, k, i, 5) * ux(j, k, i, 5));
                    sumwz2 = sumwz2 + (ux(j, k, i, 9) * ux(j, k, i, 9));
                    sumux3 = sumux3 + (ux(j, k, i, 1) * ux(j, k, i, 1) * ux(j, k, i, 1));
                    sumvy3 = sumvy3 + (ux(j, k, i, 5) * ux(j, k, i, 5) * ux(j, k, i, 5));
                    sumwz3 = sumwz3 + (ux(j, k, i, 9) * ux(j, k, i, 9) * ux(j, k, i, 9));
                    sumu2 = sumu2 + q(j, k, i, 2) * q(j, k, i, 2);
                    sumv2 = sumv2 + q(j, k, i, 3) * q(j, k, i, 3);
                    sumw2 = sumw2 + q(j, k, i, 4) * q(j, k, i, 4);
                    sumuv = sumuv + q(j, k, i, 2) * q(j, k, i, 3);
                    sumuw = sumuw + q(j, k, i, 2) * q(j, k, i, 4);
                    sumvw = sumvw + q(j, k, i, 3) * q(j, k, i, 4);
                    double energy = q(j, k, i, 5) / (double)gm1 +
                                    0.5 * q(j, k, i, 1) * (q(j, k, i, 2) * q(j, k, i, 2) +
                                                            q(j, k, i, 3) * q(j, k, i, 3) +
                                                            q(j, k, i, 4) * q(j, k, i, 4));
                    sumenergy = sumenergy + energy;
                    double s11 = ux(j, k, i, 1);
                    double s22 = ux(j, k, i, 5);
                    double s33 = ux(j, k, i, 9);
                    double s12 = 0.5 * (ux(j, k, i, 2) + ux(j, k, i, 4));
                    double s13 = 0.5 * (ux(j, k, i, 3) + ux(j, k, i, 7));
                    double s23 = 0.5 * (ux(j, k, i, 6) + ux(j, k, i, 8));
                    double w12 = 0.5 * (ux(j, k, i, 2) - ux(j, k, i, 4));
                    double w13 = 0.5 * (ux(j, k, i, 3) - ux(j, k, i, 7));
                    double w23 = 0.5 * (ux(j, k, i, 6) - ux(j, k, i, 8));
                    double xis = s11 * s11 + s22 * s22 + s33 * s33 +
                                 2.0 * s12 * s12 + 2.0 * s13 * s13 + 2.0 * s23 * s23;
                    double wis = 2.0 * w12 * w12 + 2.0 * w13 * w13 + 2.0 * w23 * w23;
                    double tt  = (double)gamma * q(j, k, i, 5) / q(j, k, i, 1);
                    fnu = fnu + (c2bp * tt * std::sqrt(tt) / (c2b + tt));
                    sumxis = sumxis + xis;
                    sumwis = sumwis + wis;
                    rho = rho + q(j, k, i, 1);
                }
            }
        }

        double avgux2  = sumux2  / (double)n;
        double avgvy2  = sumvy2  / (double)n;
        double avgwz2  = sumwz2  / (double)n;
        double avgux3  = sumux3  / (double)n;
        double avgvy3  = sumvy3  / (double)n;
        double avgwz3  = sumwz3  / (double)n;
        double avgvel2 = sumvel2 / (double)n;
        double avgu2   = sumu2   / (double)n;
        double avgv2   = sumv2   / (double)n;
        double avgw2   = sumw2   / (double)n;
        double avguv   = sumuv   / (double)n;
        double avguw   = sumuw   / (double)n;
        double avgvw   = sumvw   / (double)n;
        double avgfnu  = fnu     / (double)n;
        double avgrho  = rho     / (double)n;
        double avgwis  = sumwis  / (double)n;
        double xtime   = (double)time / eddytime;

        // for each processor, write out:
        // ntt,icyc,xtime,sumvel2,sumxis,avgux2,avgvy2,avgwz2,
        // avgux3,avgvy3,avgwz3,avgvel2,avgu2,avgv2,avgw2,avgfnu,avgrho,avgwis
        fortran_write_unit(500 + myid, "%5d%5d%15.5E%15.5E%15.5E"
            "%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E"
            "%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E"
            "%15.5E%15.5E%15.5E\n",
            (int)ntt, (int)icyc, xtime, sumvel2, sumxis,
            avgux2, avgvy2, avgwz2, avgux3, avgvy3, avgwz3, avgvel2,
            avgu2, avgv2, avgw2, avgfnu, avgrho, avgwis, sumenergy,
            avguv, avguw, avgvw);
    }

    return;
}

} // namespace lesdiag_ns
