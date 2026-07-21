// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "l2norm2.h"
#include "l2norm.h"
#include "resadd.h"
#include <cmath>
#include <algorithm>

namespace l2norm2_ns {

void l2norm(int& nbl, int& ntime, double& rmsl, int& irdq, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray3DRef<double> vol)
{
    l2norm_ns::l2norm(nbl, ntime, rmsl, irdq, jdim, kdim, idim, res, vol);
}

void l2norm2(int& nbl, int& ntime, double& rmsl, int& irdq, int& jdim, int& kdim, int& idim,
             FortranArray4DRef<double> res, FortranArray3DRef<double> vol,
             FortranArray4DRef<double> qc0, FortranArray4DRef<double> dqc0,
             FortranArray4DRef<double> q, FortranArray3DRef<double> blank)
{
    // Access COMMON block fields
    float& dt_f      = cmn_info.dt;
    int32_t& ita     = cmn_unst.ita;
    int32_t& ncyc    = cmn_mgrd.ncyc;

    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int nplq  = std::min(idim1, 999000 / (jdim * kdim));
    int npl   = nplq;
    int nt    = jdim * kdim;
    rmsl      = 0.0;

    double tfact;
    if (std::abs(ita) == 1) {
        tfact = 0.0;
    } else {
        tfact = 0.5 / (double)dt_f;
    }
    double tfacp1 = tfact + 1.0 / (double)dt_f;

    // Zero out boundary cells in res, set qc0=q at boundaries, zero dqc0
    for (int i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        for (int l = 1; l <= 5; l++) {
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                // j loop: zero boundary at k=kdim
                for (int j = 1; j <= jdim; j++) {
                    res(j, kdim, ii, l)  = 0.0;
                    qc0(j, kdim, ii, l)  = q(j, kdim, ii, l);
                    dqc0(j, kdim, ii, l) = 0.0;
                }
                // k loop: zero boundary at j=jdim
                for (int k = 1; k <= kdim1; k++) {
                    qc0(jdim, k, ii, l)  = q(jdim, k, ii, l);
                    dqc0(jdim, k, ii, l) = 0.0;
                    res(jdim, k, ii, l)  = 0.0;
                }
            }
        }
    }

    // Compute factdqc0 and factqc0
    double factdqc0, factqc0;
    if ((float)dt_f > 0.0f) {
        factdqc0 = std::min((double)(std::abs(ita) - 1), 1.0);
        if ((float)factdqc0 < 0.0f) factdqc0 = 0.0;
        factqc0 = 0.0;
        if (ncyc > 1) factqc0 = 1.0;
    } else {
        factqc0  = 0.0;
        factdqc0 = 0.0;
    }

    // Only density residual is monitored (l=1), using flat indexing
    npl = nplq;
    for (int i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        int n = npl * nt - jdim - 1;
        for (int izz = 1; izz <= n; izz++) {
            double resminus = res(izz, 1, i, 1)
                + factdqc0 * tfact * vol(izz, 1, i) * dqc0(izz, 1, i, 1)
                - factqc0 * tfacp1 * vol(izz, 1, i)
                  * (q(izz, 1, i, 1) - qc0(izz, 1, i, 1));
            rmsl = rmsl + resminus * resminus * blank(izz, 1, i);
        }
    }
}

void resadd(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qc0, FortranArray4DRef<double> dqc0,
            FortranArray4DRef<double> res, FortranArray3DRef<double> vol,
            int& iover, FortranArray3DRef<double> blank)
{
    resadd_ns::resadd(jdim, kdim, idim, q, qc0, dqc0, res, vol, iover, blank);
}

} // namespace l2norm2_ns
