// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "init.h"
#include "mms.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace init_ns {

void init(int& nbl, int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q,
          FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
          FortranArray4DRef<double> qi0, FortranArray4DRef<double> tj0,
          FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0,
          FortranArray3DRef<double> vol, FortranArray3DRef<double> volj0,
          FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0,
          int& nummem, FortranArray3DRef<double> x, FortranArray3DRef<double> y,
          FortranArray3DRef<double> z, FortranArray1DRef<int> nou,
          FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& iflagprnt)
{
    // COMMON block references
    float& gamma   = cmn_fluid.gamma;
    float& gm1     = cmn_fluid.gm1;
    float& xmach   = cmn_info.xmach;
    float& alpha   = cmn_info.alpha;
    float& beta    = cmn_info.beta;
    float& p0      = cmn_ivals.p0;
    float& rho0    = cmn_ivals.rho0;
    float& c0      = cmn_ivals.c0;
    float& u0      = cmn_ivals.u0;
    float& v0      = cmn_ivals.v0;
    float& w0      = cmn_ivals.w0;
    float& et0     = cmn_ivals.et0;
    float& h0      = cmn_ivals.h0;
    float& pt0     = cmn_ivals.pt0;
    float& rhot0   = cmn_ivals.rhot0;
    float* qiv     = cmn_ivals.qiv;    // 0-based: qiv[0]..qiv[4]
    float* tur10   = cmn_ivals.tur10;  // 0-based: tur10[0]..tur10[6]
    int32_t* ivisc = cmn_reyue.ivisc;  // 0-based: ivisc[0]..ivisc[2]
    int32_t& ivmx  = cmn_maxiv.ivmx;
    int32_t& iexact_trunc = cmn_mms.iexact_trunc;
    int32_t& iexact_disc  = cmn_mms.iexact_disc;
    int32_t& issglrrw2012 = cmn_reystressmodel.issglrrw2012;

    // c     write(15,904)nbl
    // 904 format(29h initial conditions for block,i3)

    rho0   = 1.0e0f;
    c0     = 1.0e0f;
    p0     = (float)((double)rho0 * (double)c0 * (double)c0 / (double)gamma);
    // The wind axis system follows NASA SP-3070 (1972), with the exception that
    // positive beta is in the opposite direction
    u0     = (float)((double)xmach * std::cos((double)alpha) * std::cos((double)beta));
    w0     = (float)((double)xmach * std::sin((double)alpha) * std::cos((double)beta));
    v0     = (float)(-(double)xmach * std::sin((double)beta));
    double ei0 = (double)p0 / ((double)(gamma - 1.0e0f) * (double)rho0);
    et0    = (float)((double)rho0 * (ei0 + 0.5e0 * ((double)u0*(double)u0 + (double)v0*(double)v0 + (double)w0*(double)w0)));
    h0     = (float)(((double)et0 + (double)p0) / (double)rho0);
    pt0    = (float)((double)p0 * std::pow(1.0e0 + 0.5e0*(double)gm1*(double)xmach*(double)xmach, (double)gamma/(double)gm1));
    rhot0  = (float)((double)rho0 * std::pow(1.0e0 + 0.5e0*(double)gm1*(double)xmach*(double)xmach, 1.0e0/(double)gm1));
    qiv[0] = rho0;   // qiv(1)
    qiv[1] = u0;     // qiv(2)
    qiv[2] = v0;     // qiv(3)
    qiv[3] = w0;     // qiv(4)
    qiv[4] = p0;     // qiv(5)

    // c     write(15,499) idim,jdim,kdim
    // c     write(15,500) ...

    double a1    = (double)u0;
    double a2    = (double)v0;
    double a3    = (double)w0;
    int jkpro = jdim * kdim;
    int idim1 = idim - 1;
    int npl   = 999000 / jkpro;
    int nplq  = npl;

    // Initialize q array
    for (int i = 1; i <= idim; i += nplq) {
        if (i + npl - 1 > idim) npl = idim - i + 1;
        int nv = npl * jkpro;
        // cdir$ ivdep
        // q(izz,1,i,comp) — flat indexing: izz=1..nv covers jdim*kdim*npl elements
        double* q_ptr1 = &q(1, 1, i, 1);
        double* q_ptr2 = &q(1, 1, i, 2);
        double* q_ptr3 = &q(1, 1, i, 3);
        double* q_ptr4 = &q(1, 1, i, 4);
        double* q_ptr5 = &q(1, 1, i, 5);
        for (int izz = 0; izz < nv; izz++) {
            q_ptr1[izz] = (double)rho0;
            q_ptr2[izz] = a1;
            q_ptr3[izz] = a2;
            q_ptr4[izz] = a3;
            q_ptr5[izz] = (double)p0;
        }
    }

    // Initialize qj0 array
    npl  = 999000 / kdim;
    nplq = npl;
    for (int i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        int nv = npl * kdim;
        for (int m = 1; m <= 4; m++) {
            // cdir$ ivdep
            // qj0(izz,i,comp,m) — flat indexing: izz=1..nv covers kdim*npl elements
            double* qj0_ptr1 = &qj0(1, i, 1, m);
            double* qj0_ptr2 = &qj0(1, i, 2, m);
            double* qj0_ptr3 = &qj0(1, i, 3, m);
            double* qj0_ptr4 = &qj0(1, i, 4, m);
            double* qj0_ptr5 = &qj0(1, i, 5, m);
            for (int izz = 0; izz < nv; izz++) {
                qj0_ptr1[izz] = (double)rho0;
                qj0_ptr2[izz] = a1;
                qj0_ptr3[izz] = a2;
                qj0_ptr4[izz] = a3;
                qj0_ptr5[izz] = (double)p0;
            }
        }
    }

    // Initialize qk0 array
    npl  = 999000 / kdim;
    nplq = npl;
    for (int i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        int nv = jdim * npl;
        for (int m = 1; m <= 4; m++) {
            // cdir$ ivdep
            // qk0(izz,i,comp,m) — flat indexing: izz=1..nv covers jdim*npl elements
            double* qk0_ptr1 = &qk0(1, i, 1, m);
            double* qk0_ptr2 = &qk0(1, i, 2, m);
            double* qk0_ptr3 = &qk0(1, i, 3, m);
            double* qk0_ptr4 = &qk0(1, i, 4, m);
            double* qk0_ptr5 = &qk0(1, i, 5, m);
            for (int izz = 0; izz < nv; izz++) {
                qk0_ptr1[izz] = (double)rho0;
                qk0_ptr2[izz] = a1;
                qk0_ptr3[izz] = a2;
                qk0_ptr4[izz] = a3;
                qk0_ptr5[izz] = (double)p0;
            }
        }
    }

    // Initialize qi0 array
    for (int m = 1; m <= 4; m++) {
        // cdir$ ivdep
        // qi0(izz,1,comp,m) — flat indexing: izz=1..jkpro covers jdim*kdim elements
        double* qi0_ptr1 = &qi0(1, 1, 1, m);
        double* qi0_ptr2 = &qi0(1, 1, 2, m);
        double* qi0_ptr3 = &qi0(1, 1, 3, m);
        double* qi0_ptr4 = &qi0(1, 1, 4, m);
        double* qi0_ptr5 = &qi0(1, 1, 5, m);
        for (int izz = 0; izz < jkpro; izz++) {
            qi0_ptr1[izz] = (double)rho0;
            qi0_ptr2[izz] = a1;
            qi0_ptr3[izz] = a2;
            qi0_ptr4[izz] = a3;
            qi0_ptr5[izz] = (double)p0;
        }
    }

    // For field eqn turbulence models:
    if (ivisc[0] >= 4 || ivisc[1] >= 4 || ivisc[2] >= 4) {
        if (ivisc[0] == 4 || ivisc[1] == 4 || ivisc[2] == 4) {
            if ((float)tur10[0] < 0.0f) tur10[0] = 0.1f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 0.0f;
        } else if (ivisc[0] == 5 || ivisc[1] == 5 || ivisc[2] == 5) {
            if ((float)tur10[0] < 0.0f) tur10[0] = 1.341946f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 0.0f;
        } else if (ivisc[0] == 11 || ivisc[1] == 11 || ivisc[2] == 11 ||
                   ivisc[0] == 10 || ivisc[1] == 10 || ivisc[2] == 10 ||
                   ivisc[0] ==  9 || ivisc[1] ==  9 || ivisc[2] ==  9 ||
                   ivisc[0] == 13 || ivisc[1] == 13 || ivisc[2] == 13 ||
                   ivisc[0] == 15 || ivisc[1] == 15 || ivisc[2] == 15) {
            if ((float)tur10[0] < 0.0f) tur10[0] = 1.0e-17f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 1.0e-9f;
        } else if (ivisc[0] ==  8 || ivisc[1] ==  8 || ivisc[2] ==  8 ||
                   ivisc[0] == 12 || ivisc[1] == 12 || ivisc[2] == 12 ||
                   ivisc[0] == 14 || ivisc[1] == 14 || ivisc[2] == 14) {
            if ((float)tur10[0] < 0.0f) tur10[0] = 9.0e-8f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 9.0e-9f;
        } else if (ivisc[0] == 16 || ivisc[1] == 16 || ivisc[2] == 16) {
            if ((float)tur10[0] < 0.0f) tur10[0] = 1.5589e-6f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 9.0e-9f;
        } else if (ivisc[0] == 25 || ivisc[1] == 25 || ivisc[2] == 25) {
            // tur10 and tur20 not used for model 25, but we set them anyway for safety:
            if ((float)tur10[0] < 0.0f) tur10[0] = 1.0e-6f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 9.0e-9f;
        } else if (ivisc[0] ==  6 || ivisc[1] ==  6 || ivisc[2] ==  6 ||
                   ivisc[0] ==  7 || ivisc[1] ==  7 || ivisc[2] ==  7) {
            if ((float)tur10[0] < 0.0f) tur10[0] = 1.0e-6f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 9.0e-9f;
        } else if (ivisc[0] == 30 || ivisc[1] == 30 || ivisc[2] == 30) {
            if ((float)tur10[0] < 0.0f) tur10[0] = 1.0e-6f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 9.0e-9f;
            if ((float)tur10[2] < 0.0f) tur10[2] = 1.0f;
        } else if (ivisc[0] == 40 || ivisc[1] == 40 || ivisc[2] == 40) {
            if ((float)tur10[0] < 0.0f) tur10[0] = 1.0e-6f;
            if ((float)tur10[1] < 0.0f) tur10[1] = 9.0e-9f;
            if ((float)tur10[2] < 0.0f) tur10[2] = 1.0f;
            if ((float)tur10[3] < 0.0f) {
                double turbintensity_inf_percent = 100.0 / (double)xmach *
                    std::sqrt((double)tur10[1] / 1.5);
                if ((float)turbintensity_inf_percent <= 1.3f) {
                    tur10[3] = (float)(1173.51 - 589.428 * turbintensity_inf_percent +
                        0.2196 / (turbintensity_inf_percent * turbintensity_inf_percent));
                } else {
                    tur10[3] = (float)(331.5 * std::pow(turbintensity_inf_percent - 0.5658, -0.671));
                }
            }
        } else if (ivisc[0] == 72 || ivisc[1] == 72 || ivisc[2] == 72) {
            if ((float)tur10[0] == -1.0f) tur10[0] = (float)(-2.0/3.0 * 9.0e-9);
            if ((float)tur10[1] == -1.0f) tur10[1] = (float)(-2.0/3.0 * 9.0e-9);
            if ((float)tur10[2] == -1.0f) tur10[2] = (float)(-2.0/3.0 * 9.0e-9);
            if ((float)tur10[3] == -1.0f) tur10[3] = 0.0f;
            if ((float)tur10[4] == -1.0f) tur10[4] = 0.0f;
            if ((float)tur10[5] == -1.0f) tur10[5] = 0.0f;
            if (issglrrw2012 == 6) {
                if ((float)tur10[6] == -1.0f) tur10[6] = (float)(1.0 / std::sqrt(1.0e-6));
            } else {
                if ((float)tur10[6] == -1.0f) tur10[6] = 1.0e-6f;
            }
            // Write out some diagnostic info:
            // CFL3D has default tur10 values; can be overridden via
            // turbintensity_inf_percent and/or eddy_visc_inf (see readkey.F)
            if (nbl == 1 && iflagprnt == 0) {
                double tke = -((double)tur10[0] + (double)tur10[1] + (double)tur10[2]) / 2.0;
                double tu  = std::sqrt(tke / 1.5) / (double)xmach;
                double xnu_ratio = tke / (double)tur10[6];
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " For Stress-Omega RSM, ivisc=72 (all blocks):");
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " initial Tu =%12.5E", tu);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " initial xnu_ratio =%12.5E", xnu_ratio);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " initial tke =%12.5E", tke);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " initial omega =%12.5E", (double)tur10[6]);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " tur10(1) =%12.5E", (double)tur10[0]);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " tur10(2) =%12.5E", (double)tur10[1]);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " tur10(3) =%12.5E", (double)tur10[2]);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " tur10(4) =%12.5E", (double)tur10[3]);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " tur10(5) =%12.5E", (double)tur10[4]);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " tur10(6) =%12.5E", (double)tur10[5]);
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    " tur10(7) =%12.5E", (double)tur10[6]);
            }
        }

        // Fill tj0, tk0, ti0 with tur10 values
        for (int nn = 1; nn <= nummem; nn++) {
            for (int m = 1; m <= 4; m++) {
                for (int i = 1; i <= idim - 1; i++) {
                    for (int j = 1; j <= jdim; j++) {
                        tk0(j, i, nn, m) = (double)tur10[nn - 1];
                    }
                    for (int k = 1; k <= kdim; k++) {
                        tj0(k, i, nn, m) = (double)tur10[nn - 1];
                    }
                }
                for (int k = 1; k <= kdim; k++) {
                    for (int j = 1; j <= jdim; j++) {
                        ti0(j, k, nn, m) = (double)tur10[nn - 1];
                    }
                }
            }
        }

    } // end if ivisc >= 4

    // boundary volumes (default to interior cells)
    if (ivmx > 0) {
        int kk = 0;
        for (int m = 1; m <= 4; m++) {
            if (m == 1) kk = 1;
            if (m == 2) kk = std::min(2, kdim - 1);
            if (m == 3) kk = kdim - 1;
            if (m == 4) kk = std::max(1, kdim - 2);
            for (int i = 1; i <= idim - 1; i++) {
                for (int j = 1; j <= jdim; j++) {
                    volk0(j, i, m) = vol(j, kk, i);
                }
            }
        }
        int jj = 0;
        for (int m = 1; m <= 4; m++) {
            if (m == 1) jj = 1;
            if (m == 2) jj = std::min(2, jdim - 1);
            if (m == 3) jj = jdim - 1;
            if (m == 4) jj = std::max(1, jdim - 2);
            // the sgi f77 compiler with -O3 won't do the following
            // assignment correctly if the k-loop is innermost
            for (int k = 1; k <= kdim; k++) {
                for (int i = 1; i <= idim - 1; i++) {
                    volj0(k, i, m) = vol(jj, k, i);
                }
            }
        }
        int ii = 0;
        for (int m = 1; m <= 4; m++) {
            if (m == 1) ii = 1;
            if (m == 2) ii = std::min(2, idim - 1);
            if (m == 3) ii = idim - 1;
            if (m == 4) ii = std::max(1, idim - 2);
            for (int k = 1; k <= kdim; k++) {
                for (int j = 1; j <= jdim; j++) {
                    voli0(j, k, m) = vol(j, k, ii);
                }
            }
        }
    }

    // Overwrite with exact soln if doing MMS
    if (iexact_trunc != 0 || iexact_disc != 0) {
        mms_ns::exact_flow_q(jdim, kdim, idim, x, y, z, q, iexact_trunc, iexact_disc);
    }

    return;
}

} // namespace init_ns
