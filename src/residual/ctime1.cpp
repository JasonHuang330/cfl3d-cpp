// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ctime1.h"
#include "ccomplex.h"
#include "q8sdot.h"
#include "q8smax.h"
#include "q8smin.h"
#include "xmukin.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cmath>

namespace ctime1_ns {

void ctime1(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q,
            FortranArray3DRef<double> vol, FortranArray3DRef<double> sj,
            FortranArray3DRef<double> sk, FortranArray3DRef<double> si,
            FortranArray3DRef<double> dtj, FortranArray2DRef<double> t,
            double& delt, FortranArray3DRef<double> vist3d, int& itur,
            double& dtmin, int& iout, int& ntime, FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& idef)
{
    // COMMON block aliases
    float& gamma    = cmn_fluid.gamma;
    float& pr       = cmn_fluid2.pr;
    float& xmach    = cmn_info.xmach;
    float& dt       = cmn_info.dt;
    float& reue     = cmn_reyue.reue;
    int32_t& i2d    = cmn_twod.i2d;
    float& cfltau   = cmn_unst.cfltau;

    int32_t& ita    = cmn_unst.ita;
    int32_t& iunst  = cmn_unst.iunst;
    float& cfltau0  = cmn_unst.cfltau0;
    float& cfltauMax = cmn_unst.cfltaumax;
    float& cprec    = cmn_precond.cprec;
    float& uref     = cmn_precond.uref;
    float& avn      = cmn_precond.avn;
    int32_t& iaxi2plane = cmn_axisym.iaxi2plane;
    int32_t& ncyc   = cmn_mgrd.ncyc;
    int32_t& icyc   = cmn_mgrd.icyc;

    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    double dttemp = (double)dt;

    // set flags for contributions from pseudo- and physical- time steps
    double fact1, fact2, tfacp1;
    if ((double)dt > 0.0) {
        if (ita < 0) {
            if ((double)cfltauMax > (double)cfltau && ncyc > 2) {
                dt = (float)(-(cfltau +
                    (cfltauMax - cfltau) * std::pow((icyc - 1.0) / (ncyc - 1.0), (double)cfltau0)));
            } else {
                dt = -cfltau;
            }
            fact1 = 1.0;
            fact2 = 1.0;
        } else {
            fact1 = 0.0;
            fact2 = 1.0;
        }
        double tfact;
        if (std::abs(ita) == 1) {
            tfact = 0.0e0;
        } else {
            tfact = 0.5e0;
        }
        tfacp1 = tfact + 1.0e0;
    } else {
        fact1 = 1.0;
        fact2 = 0.0;
        tfacp1 = 1.0e0;
    }

    if (iout > 0 && dttemp < 0.0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " ");
        nou(1) = std::min(nou(1) + 1, ibufdim);
        // format 7: (1x,42hcomputing time step distribution for block,i5)
        std::snprintf(bou(nou(1), 1), 120, " computing time step distribution for block%5d", nbl);
        nou(1) = std::min(nou(1) + 1, ibufdim);
        // format 17: (3x,25hsummary of time step data)
        std::snprintf(bou(nou(1), 1), 120, "   summary of time step data");
    }

    if (iout > 0 && dttemp > 0.0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120, " ");
        nou(1) = std::min(nou(1) + 1, ibufdim);
        // format 8: (1x,43hcomputing CFL number distribution for block,i5)
        std::snprintf(bou(nou(1), 1), 120, " computing CFL number distribution for block%5d", nbl);
        nou(1) = std::min(nou(1) + 1, ibufdim);
        // format 18: (3x,26hsummary of CFL number data)
        std::snprintf(bou(nou(1), 1), 120, "   summary of CFL number data");
    }

    // vterm = ccmaxrc(1.3333, gamma/pr) * xmach / reue
    double vterm_a = 1.3333;
    double vterm_b = (double)gamma / (double)pr;
    double vterm = ccomplex_ns::ccmaxrc(vterm_a, vterm_b) * (double)xmach / (double)reue;
    double term = 1.0;
    if (i2d == 1 || iaxi2plane == 1) term = 0.0;
    delt = 0.0;
    if (iunst > 0 && dttemp > 0.0) delt = dttemp;
    double dt_val = (double)dt;
    double cfl1_a = dt_val;
    double cfl1 = 1.0e0 / ccomplex_ns::ccabs(cfl1_a);
    double dt2 = dt_val * dt_val;
    int n = jdim * kdim - jdim - 1;

    if (iout > 0) {
        if (dttemp < 0.0) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            // format 1241: (3x,1hI,6x,3hcfl,10x,5hdtrms,9x,5hdtmin,9x,5hdtmax)
            std::snprintf(bou(nou(1), 1), 120, "   I      cfl          dtrms         dtmin         dtmax");
        } else {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            // format 1247: (3x,1hI,7x,2hdt,9x,6hcflrms,8x,6hcflmin,8x,6hcflmax)
            std::snprintf(bou(nou(1), 1), 120, "   I       dt         cflrms        cflmin        cflmax");
        }
    }

    for (int i = 1; i <= idim1; i++) {

        n = jdim * kdim - jdim - 1;

        if (itur == 0) {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 1) = 0.0;
            }
        } else {
            for (int k = 1; k <= kdim - 1; k++) {
                int izz = jdim * (k - 1);
                for (int j = 1; j <= jdim - 1; j++) {
                    t(izz + j, 2) = (double)gamma * q(j, k, i, 5) / q(j, k, i, 1);
                }
                // xmukin(jdim-1, t(izz+1,2), t(izz+1,1), tinf)
                {
                    int n_xmu = jdim - 1;
                    double tinf_val = (double)cmn_reyue.tinf;
                    FortranArray1DRef<double> temp_ref(&t(izz + 1, 2), jdim - 1);
                    FortranArray1DRef<double> visc_ref(&t(izz + 1, 1), jdim - 1);
                    xmukin_ns::xmukin(n_xmu, temp_ref, visc_ref, tinf_val);
                }
                t(izz + jdim, 1) = t(izz + jdim - 1, 1);
            }
        }

        if (itur >= 2) {
            for (int k = 1; k <= kdim - 1; k++) {
                int izz = jdim * (k - 1);
                for (int j = 1; j <= jdim - 1; j++) {
                    t(izz + j, 1) = t(izz + j, 1) + vist3d(j, k, i);
                }
            }
        }

        if ((double)cprec == 0.0) {
            // cdir$ ivdep
            for (int izz = 1; izz <= n; izz++) {
                double vrho = t(izz, 1) * vterm / (vol(izz, 1, i) * q(izz, 1, i, 1));

                double t11 = sj(izz, i, 1) + sj(izz + 1, i, 1);
                double t12 = sj(izz, i, 2) + sj(izz + 1, i, 2);
                double t13 = sj(izz, i, 3) + sj(izz + 1, i, 3);
                double t14 = sj(izz, i, 4) + sj(izz + 1, i, 4);
                double t1  = (t11 * q(izz, 1, i, 2) + t12 * q(izz, 1, i, 3)
                           +  t13 * q(izz, 1, i, 4)) * 0.5e0;
                // add cell face speed for moving grids
                double t20 = sj(izz, i, 5) + sj(izz + 1, i, 5);
                t1 = t1 + t20 * 0.5e0;
                {
                    double t1_tmp = t1;
                    t1 = ccomplex_ns::ccabs(t1_tmp) + t14 * vrho;
                }

                t11 = sk(izz, i, 1) + sk(izz + jdim, i, 1);
                t12 = sk(izz, i, 2) + sk(izz + jdim, i, 2);
                t13 = sk(izz, i, 3) + sk(izz + jdim, i, 3);
                double t15 = sk(izz, i, 4) + sk(izz + jdim, i, 4);
                double t2  = (t11 * q(izz, 1, i, 2) + t12 * q(izz, 1, i, 3)
                           +  t13 * q(izz, 1, i, 4)) * 0.5e0;
                // add cell face speed for moving grids
                t20 = sk(izz, i, 5) + sk(izz + jdim, i, 5);
                t2 = t2 + t20 * 0.5e0;
                {
                    double t2_tmp = t2;
                    t2 = ccomplex_ns::ccabs(t2_tmp) + t15 * vrho;
                }

                t11 = si(izz, i, 1) + si(izz, i + 1, 1);
                t12 = si(izz, i, 2) + si(izz, i + 1, 2);
                t13 = si(izz, i, 3) + si(izz, i + 1, 3);
                double t16 = si(izz, i, 4) + si(izz, i + 1, 4);
                double t3  = (t11 * q(izz, 1, i, 2) + t12 * q(izz, 1, i, 3)
                           +  t13 * q(izz, 1, i, 4)) * 0.5e0;
                // add cell face speed for moving grids
                t20 = si(izz, i, 5) + si(izz, i + 1, 5);
                t3 = t3 + t20 * 0.5e0;
                {
                    double t3_tmp = t3;
                    t3 = ccomplex_ns::ccabs(t3_tmp) + t16 * vrho;
                }

                t11 = (double)gamma * q(izz, 1, i, 5) / q(izz, 1, i, 1);
                t11 = std::sqrt(t11);

                dtj(izz, 1, i) = (t1 + t11) * t14 + (t2 + t11) * t15 + (t3 + t11) * t16 * term;
                dtj(izz, 1, i) = 0.5e0 * cfl1 * dtj(izz, 1, i);
                {
                    double dttemp_dt = dttemp / (double)dt;
                    double abs_val = dttemp_dt;
                    t(izz, 10) = vol(izz, 1, i) / dtj(izz, 1, i) * ccomplex_ns::ccabs(abs_val);
                }
            }
        } else {
            // cdir$ ivdep
            for (int izz = 1; izz <= n; izz++) {
                double vrho = t(izz, 1) * vterm / (vol(izz, 1, i) * q(izz, 1, i, 1));
                double c2 = (double)gamma * q(izz, 1, i, 5) / q(izz, 1, i, 1);
                double c = std::sqrt(c2);
                double vmag1 = q(izz, 1, i, 2) * q(izz, 1, i, 2)
                             + q(izz, 1, i, 3) * q(izz, 1, i, 3)
                             + q(izz, 1, i, 4) * q(izz, 1, i, 4);
                double avn_uref2 = (double)avn * (double)uref * (double)uref;
                double vel2 = ccomplex_ns::ccmax(vmag1, avn_uref2);
                double vel_tmp = ccomplex_ns::ccmin(c2, vel2);
                double vel = std::sqrt(vel_tmp);
                vel = (double)cprec * vel + (1.0 - (double)cprec) * c;
                double xm2 = (vel / c) * (vel / c);

                double t11 = sj(izz, i, 1) + sj(izz + 1, i, 1);
                double t12 = sj(izz, i, 2) + sj(izz + 1, i, 2);
                double t13 = sj(izz, i, 3) + sj(izz + 1, i, 3);
                double t14 = sj(izz, i, 4) + sj(izz + 1, i, 4);
                double t1  = (t11 * q(izz, 1, i, 2) + t12 * q(izz, 1, i, 3)
                           +  t13 * q(izz, 1, i, 4)) * 0.5e0;
                double xmave = t1 / c;
                double tt1j = 0.5 * (1.0 + xm2);
                double tt2j = 0.5 * std::sqrt(xmave * xmave * (1.0 - xm2) * (1.0 - xm2) + 4.0 * xm2);
                // add cell face speed for moving grids
                double t20 = sj(izz, i, 5) + sj(izz + 1, i, 5);
                t1 = tt1j * t1 + t20 * 0.5e0;
                {
                    double t1_tmp = t1;
                    t1 = ccomplex_ns::ccabs(t1_tmp) + t14 * vrho;
                }

                t11 = sk(izz, i, 1) + sk(izz + jdim, i, 1);
                t12 = sk(izz, i, 2) + sk(izz + jdim, i, 2);
                t13 = sk(izz, i, 3) + sk(izz + jdim, i, 3);
                double t15 = sk(izz, i, 4) + sk(izz + jdim, i, 4);
                double t2  = (t11 * q(izz, 1, i, 2) + t12 * q(izz, 1, i, 3)
                           +  t13 * q(izz, 1, i, 4)) * 0.5e0;
                xmave = t2 / c;
                double tt1k = 0.5 * (1.0 + xm2);
                double tt2k = 0.5 * std::sqrt(xmave * xmave * (1.0 - xm2) * (1.0 - xm2) + 4.0 * xm2);
                // add cell face speed for moving grids
                t20 = sk(izz, i, 5) + sk(izz + jdim, i, 5);
                t2 = tt1k * t2 + t20 * 0.5e0;
                {
                    double t2_tmp = t2;
                    t2 = ccomplex_ns::ccabs(t2_tmp) + t15 * vrho;
                }

                t11 = si(izz, i, 1) + si(izz, i + 1, 1);
                t12 = si(izz, i, 2) + si(izz, i + 1, 2);
                t13 = si(izz, i, 3) + si(izz, i + 1, 3);
                double t16 = si(izz, i, 4) + si(izz, i + 1, 4);
                double t3  = (t11 * q(izz, 1, i, 2) + t12 * q(izz, 1, i, 3)
                           +  t13 * q(izz, 1, i, 4)) * 0.5e0;
                xmave = t3 / c;
                double tt1i = 0.5 * (1.0 + xm2);
                double tt2i = 0.5 * std::sqrt(xmave * xmave * (1.0 - xm2) * (1.0 - xm2) + 4.0 * xm2);
                // add cell face speed for moving grids
                t20 = si(izz, i, 5) + si(izz, i + 1, 5);
                t3 = tt1i * t3 + t20 * 0.5e0;
                {
                    double t3_tmp = t3;
                    t3 = ccomplex_ns::ccabs(t3_tmp) + t16 * vrho;
                }

                double c1j = c * tt2j;
                double c1k = c * tt2k;
                double c1i = c * tt2i;

                dtj(izz, 1, i) = (t1 + c1j) * t14 + (t2 + c1k) * t15 + (t3 + c1i) * t16 * term;
                dtj(izz, 1, i) = 0.5e0 * cfl1 * dtj(izz, 1, i);
                {
                    double dttemp_dt = dttemp / (double)dt;
                    double abs_val = dttemp_dt;
                    t(izz, 10) = vol(izz, 1, i) / dtj(izz, 1, i) * ccomplex_ns::ccabs(abs_val);
                }
            }
        }

        if (dttemp > 0.0e0) {
            dt2 = dttemp * dttemp;
            // cdir$ ivdep
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 10) = dt2 / t(izz, 10);
                // code can only do 1st order temporal for pseudo-time term (with
                // subiterations). Therefore, pseudo-time term is divided here by
                // tfacp1. In several subroutines in af3f, dtj is multiplied by
                // tfacp1.
                dtj(izz, 1, i) = fact1 * dtj(izz, 1, i) / tfacp1
                               + fact2 * vol(izz, 1, i) / dttemp;
            }
            dt2 = (double)dt * (double)dt;

            // geometric conservation law terms for deforming grids
            if (idef > 0) {
                for (int izz = 1; izz <= n; izz++) {
                    t(izz, 17) = sj(izz + 1, i, 5) * sj(izz + 1, i, 4)
                               - sj(izz, i, 5) * sj(izz, i, 4)
                               + sk(izz + jdim, i, 5) * sk(izz + jdim, i, 4)
                               - sk(izz, i, 5) * sk(izz, i, 4)
                               + si(izz, i + 1, 5) * si(izz, i + 1, 4)
                               - si(izz, i, 5) * si(izz, i, 4);
                }
                for (int izz = jdim; izz <= n; izz += jdim) {
                    t(izz, 17) = sj(izz, i, 5) * sj(izz, i, 4)
                               - sj(izz - 1, i, 5) * sj(izz - 1, i, 4)
                               + sk(izz + jdim, i, 5) * sk(izz + jdim, i, 4)
                               - sk(izz, i, 5) * sk(izz, i, 4)
                               + si(izz, i + 1, 5) * si(izz, i + 1, 4)
                               - si(izz, i, 5) * si(izz, i, 4);
                }
                for (int izz = 1; izz <= n; izz++) {
                    dtj(izz, 1, i) = dtj(izz, 1, i) + t(izz, 17) / tfacp1;
                }
            }
        }

        for (int kk = 1; kk <= kdim1; kk++) {
            int jk = jdim * kk;
            t(jk, 10) = 0.0e0;
        }
        {
            int n_dot = n;
            FortranArray1DRef<double> t10_ref(&t(1, 10), n);
            double dtrms = q8sdot_ns::q8sdot(n_dot, t10_ref, n_dot, t10_ref);
            dtrms = std::sqrt(dtrms / (double)(jdim1 * kdim1));
            double dtmax = q8smax_ns::q8smax(n_dot, t10_ref);
            for (int kk = 1; kk <= kdim1; kk++) {
                int jk = jdim * kk;
                t(jk, 10) = dtmax;
            }
            FortranArray1DRef<double> t10_ref2(&t(1, 10), n);
            dtmin = q8smin_ns::q8smin(n_dot, t10_ref2);
            double dtpr_val = dttemp;
            double dtpr = ccomplex_ns::ccabs(dtpr_val);
            if (iout > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 1201: (i4,1x,e12.5,3(2x,e12.5))
                std::snprintf(bou(nou(1), 1), 120, "%4d %12.5E  %12.5E  %12.5E  %12.5E",
                              i, (float)dtpr, (float)dtrms, (float)dtmin, (float)dtmax);
            }
        }

        // cdir$ ivdep
        for (int izz = 1; izz <= jdim + 1; izz++) {
            dtj(izz + jdim - 1, kdim1, i) = 1.0e0;
        }

    } // end do 9000

    if ((double)dt != dttemp) {
        dt = (float)dttemp;
    }

    return;
}

} // namespace ctime1_ns
