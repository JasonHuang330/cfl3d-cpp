// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "coll2q.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <algorithm>

namespace coll2q_ns {

void coll2q(FortranArray4DRef<double> q, FortranArray4DRef<double> qc,
            FortranArray3DRef<double> vol, FortranArray3DRef<double> volc,
            int& jdim, int& kdim, int& idim,
            int& jc, int& kc, int& ic,
            FortranArray4DRef<double> res, FortranArray4DRef<double> qr,
            FortranArray4DRef<double> qw,
            int& js, int& ks, int& is, int& je, int& ke, int& ie,
            int& nbl, int& nblc,
            FortranArray3DRef<double> vistf, FortranArray3DRef<double> vistc,
            FortranArray4DRef<double> tursavf, FortranArray4DRef<double> tursavc,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& nummem)
{
    // COMMON blocks
    int* ivisc = cmn_reyue.ivisc;  // 0-based C array, ivisc[0]=ivisc(1), etc.
    int& isklton = cmn_sklton.isklton;

    // Diagnostic output
    if (isklton > 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::memset(bou(nou(1), 1), ' ', 120);
        // write blank line
        nou(1) = std::min(nou(1) + 1, ibufdim);
        // format 7: "restricting variables and residual from embedded blockNNN to coarser blockNNN"
        std::snprintf(bou(nou(1), 1), 120,
            " restricting variables and residual from embedded block%3d to coarser block%3d",
            nbl, nblc);
    }

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    // Compute qw = q * vol using flat indexing over j,k plane
    // nv = jdim*kdim; loop izz=1..nv
    int nv = jdim * kdim;
    for (int n = 1; n <= 5; n++) {
        for (int i = 1; i <= idim1; i++) {
            // flat pointer access: qw(izz,1,i,n) for izz=1..nv
            double* qw_ptr  = &qw(1, 1, i, n);
            double* q_ptr   = &q(1, 1, i, n);
            double* vol_ptr = &vol(1, 1, i);
            for (int izz = 0; izz < nv; izz++) {
                qw_ptr[izz] = q_ptr[izz] * vol_ptr[izz];
            }
        }
    }

    // Semi-coarsening / directional refinement check
    int nsi = (idim - 1) / (ie - is);

    if (nsi == 2) {
        // 8-point restriction (i steps by 2)
        for (int n = 1; n <= 5; n++) {
            int kk = 0;
            int kq = ks - 1;
            for (int k = 1; k <= kdim1; k += 2) {
                kk++;
                kq++;
                int ii = 0;
                int iq = is - 1;
                for (int i = 1; i <= idim1; i += 2) {
                    ii++;
                    iq++;
                    int jj = 0;
                    int jq = js - 1;
                    for (int j = 1; j <= jdim1; j += 2) {
                        jj++;
                        jq++;
                        qc(jq, kq, iq, n) = (qw(j,   k,   i,   n) + qw(j,   k,   i+1, n)
                                            + qw(j+1, k,   i,   n) + qw(j+1, k,   i+1, n)
                                            + qw(j,   k+1, i,   n) + qw(j,   k+1, i+1, n)
                                            + qw(j+1, k+1, i,   n) + qw(j+1, k+1, i+1, n))
                                           / (vol(j,   k,   i)   + vol(j,   k,   i+1)
                                            + vol(j+1, k,   i)   + vol(j+1, k,   i+1)
                                            + vol(j,   k+1, i)   + vol(j,   k+1, i+1)
                                            + vol(j+1, k+1, i)   + vol(j+1, k+1, i+1));
                        qr(jj, kk, ii, n) = (res(j,   k,   i,   n) + res(j,   k,   i+1, n)
                                           + res(j+1, k,   i,   n) + res(j+1, k,   i+1, n)
                                           + res(j,   k+1, i,   n) + res(j,   k+1, i+1, n)
                                           + res(j+1, k+1, i,   n) + res(j+1, k+1, i+1, n));
                        if ((ivisc[0] > 1 && n == 1) || (ivisc[1] > 1 && n == 1)
                            || (ivisc[2] > 1 && n == 1)) {
                            vistc(jq, kq, iq) = 0.125e0 * (vistf(j,   k,   i)
                                                           + vistf(j,   k,   i+1)
                                                           + vistf(j+1, k,   i)
                                                           + vistf(j+1, k,   i+1)
                                                           + vistf(j,   k+1, i)
                                                           + vistf(j,   k+1, i+1)
                                                           + vistf(j+1, k+1, i)
                                                           + vistf(j+1, k+1, i+1));
                        }
                        if ((ivisc[0] == 4 && n == 1) || (ivisc[1] == 4 && n == 1)
                            || (ivisc[2] == 4 && n == 1)
                            || (ivisc[0] == 5 && n == 1) || (ivisc[1] == 5 && n == 1)
                            || (ivisc[2] == 5 && n == 1)) {
                            tursavc(jq, kq, iq, n) = 0.125e0 * (tursavf(j,   k,   i,   n)
                                                                + tursavf(j,   k,   i+1, n)
                                                                + tursavf(j+1, k,   i,   n)
                                                                + tursavf(j+1, k,   i+1, n)
                                                                + tursavf(j,   k+1, i,   n)
                                                                + tursavf(j,   k+1, i+1, n)
                                                                + tursavf(j+1, k+1, i,   n)
                                                                + tursavf(j+1, k+1, i+1, n));
                        }
                    }
                }
            }
        }
        // ivisc > 5: restrict tursavc for all nummem
        if (ivisc[0] > 5 || ivisc[1] > 5 || ivisc[2] > 5) {
            for (int n = 1; n <= nummem; n++) {
                int kk = 0;
                int kq = ks - 1;
                for (int k = 1; k <= kdim1; k += 2) {
                    kk++;
                    kq++;
                    int ii = 0;
                    int iq = is - 1;
                    for (int i = 1; i <= idim1; i += 2) {
                        ii++;
                        iq++;
                        int jj = 0;
                        int jq = js - 1;
                        for (int j = 1; j <= jdim1; j += 2) {
                            jj++;
                            jq++;
                            tursavc(jq, kq, iq, n) = 0.125e0 * (tursavf(j,   k,   i,   n)
                                                                + tursavf(j,   k,   i+1, n)
                                                                + tursavf(j+1, k,   i,   n)
                                                                + tursavf(j+1, k,   i+1, n)
                                                                + tursavf(j,   k+1, i,   n)
                                                                + tursavf(j,   k+1, i+1, n)
                                                                + tursavf(j+1, k+1, i,   n)
                                                                + tursavf(j+1, k+1, i+1, n));
                        }
                    }
                }
            }
        }

    } else {
        // 4-point restriction (i steps by 1)
        for (int n = 1; n <= 5; n++) {
            int kk = 0;
            int kq = ks - 1;
            for (int k = 1; k <= kdim1; k += 2) {
                kk++;
                kq++;
                int ii = 0;
                int iq = is - 1;
                for (int i = 1; i <= idim1; i += 1) {
                    ii++;
                    iq++;
                    int jj = 0;
                    int jq = js - 1;
                    for (int j = 1; j <= jdim1; j += 2) {
                        jj++;
                        jq++;
                        qc(jq, kq, iq, n) = (qw(j,   k,   i, n) + qw(j+1, k,   i, n)
                                            + qw(j,   k+1, i, n) + qw(j+1, k+1, i, n))
                                           / (vol(j,   k,   i)   + vol(j+1, k,   i)
                                            + vol(j,   k+1, i)   + vol(j+1, k+1, i));
                        qr(jj, kk, ii, n) = res(j,   k,   i, n) + res(j+1, k,   i, n)
                                          + res(j,   k+1, i, n) + res(j+1, k+1, i, n);
                        if ((ivisc[0] > 1 && n == 1) || (ivisc[1] > 1 && n == 1)
                            || (ivisc[2] > 1 && n == 1)) {
                            vistc(jq, kq, iq) = 0.25e0 * (vistf(j,   k,   i)
                                                          + vistf(j+1, k,   i)
                                                          + vistf(j,   k+1, i)
                                                          + vistf(j+1, k+1, i));
                        }
                        if ((ivisc[0] == 4 && n == 1) || (ivisc[1] == 4 && n == 1)
                            || (ivisc[2] == 4 && n == 1)
                            || (ivisc[0] == 5 && n == 1) || (ivisc[1] == 5 && n == 1)
                            || (ivisc[2] == 5 && n == 1)) {
                            tursavc(jq, kq, iq, n) = 0.25e0 * (tursavf(j,   k,   i, n)
                                                               + tursavf(j+1, k,   i, n)
                                                               + tursavf(j,   k+1, i, n)
                                                               + tursavf(j+1, k+1, i, n));
                        }
                    }
                }
            }
        }
        // ivisc > 5: restrict tursavc for all nummem
        if (ivisc[0] > 5 || ivisc[1] > 5 || ivisc[2] > 5) {
            for (int n = 1; n <= nummem; n++) {
                int kk = 0;
                int kq = ks - 1;
                for (int k = 1; k <= kdim1; k += 2) {
                    kk++;
                    kq++;
                    int ii = 0;
                    int iq = is - 1;
                    for (int i = 1; i <= idim1; i += 1) {
                        ii++;
                        iq++;
                        int jj = 0;
                        int jq = js - 1;
                        for (int j = 1; j <= jdim1; j += 2) {
                            jj++;
                            jq++;
                            tursavc(jq, kq, iq, n) = 0.25e0 * (tursavf(j,   k,   i, n)
                                                               + tursavf(j+1, k,   i, n)
                                                               + tursavf(j,   k+1, i, n)
                                                               + tursavf(j+1, k+1, i, n));
                        }
                    }
                }
            }
        }
    }

    return;
}

} // namespace coll2q_ns
