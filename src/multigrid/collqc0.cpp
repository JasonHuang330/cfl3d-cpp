// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "collqc0.h"
#include "collq.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <algorithm>

namespace collqc0_ns {

void collq(FortranArray4DRef<double> q, FortranArray4DRef<double> qq,
           FortranArray3DRef<double> vol, FortranArray3DRef<double> volc,
           int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2,
           FortranArray4DRef<double> res, FortranArray4DRef<double> qr,
           FortranArray4DRef<double> qw, FortranArray3DRef<double> vistf,
           FortranArray3DRef<double> vistc, FortranArray4DRef<double> tursavf,
           FortranArray4DRef<double> tursavc, int& nbl,
           FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim, int& nummem)
{
    collq_ns::collq(q, qq, vol, volc, jdim, kdim, idim, jj2, kk2, ii2,
                    res, qr, qw, vistf, vistc, tursavf, tursavc, nbl,
                    nou, bou, nbuf, ibufdim, nummem);
}

void collqc0(FortranArray4DRef<double> qc0, FortranArray4DRef<double> qc0c,
             FortranArray3DRef<double> vol, FortranArray3DRef<double> volc,
             int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2,
             FortranArray4DRef<double> dqc0, FortranArray4DRef<double> dqc0c,
             FortranArray4DRef<double> qw, int& nbl,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    int& isklton = cmn_sklton.isklton;

    int nbl1 = nbl + 1;
    if (isklton > 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            " restricting qc0       and     dqc0 from finer block%4d to coarser block%4d",
            nbl, nbl1);
    }

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;
    int jjl   = jj2 - 1;
    int kkl   = kk2 - 1;
    int iil   = ii2 - 1;
    // n = jj2*kk2  (computed but unused in original)

    // restrict Q(n)
    for (int n = 1; n <= 5; n++) {
        if (idim > 2) {
            int kk = 0;
            for (int k = 1; k <= kdim1; k += 2) {
                kk = kk + 1;
                int jj = 0;
                for (int j = 1; j <= jdim1; j += 2) {
                    jj = jj + 1;
                    int ii = 0;
                    for (int i = 1; i <= idim1; i += 2) {
                        ii = ii + 1;
                        qc0c(jj, kk, ii, n) =
                            qc0(j,   k,   i,   n) + qc0(j,   k,   i+1, n)
                          + qc0(j+1, k,   i,   n) + qc0(j+1, k,   i+1, n)
                          + qc0(j,   k+1, i,   n) + qc0(j,   k+1, i+1, n)
                          + qc0(j+1, k+1, i,   n) + qc0(j+1, k+1, i+1, n);
                    }
                }
            }
        } else {
            int ii = 1;
            int i  = 1;
            int kk = 0;
            for (int k = 1; k <= kdim1; k += 2) {
                kk = kk + 1;
                int jj = 0;
                for (int j = 1; j <= jdim1; j += 2) {
                    jj = jj + 1;
                    qc0c(jj, kk, ii, n) =
                        qc0(j,   k,   i, n) + qc0(j+1, k,   i, n)
                      + qc0(j,   k+1, i, n) + qc0(j+1, k+1, i, n);
                }
            }
        }
    }

    // restrict Q(n)-Q(n-1)
    for (int n = 1; n <= 5; n++) {
        if (idim > 2) {
            int kk = 0;
            for (int k = 1; k <= kdim1; k += 2) {
                kk = kk + 1;
                int jj = 0;
                for (int j = 1; j <= jdim1; j += 2) {
                    jj = jj + 1;
                    int ii = 0;
                    for (int i = 1; i <= idim1; i += 2) {
                        ii = ii + 1;
                        dqc0c(jj, kk, ii, n) =
                            dqc0(j,   k,   i,   n) + dqc0(j,   k,   i+1, n)
                          + dqc0(j+1, k,   i,   n) + dqc0(j+1, k,   i+1, n)
                          + dqc0(j,   k+1, i,   n) + dqc0(j,   k+1, i+1, n)
                          + dqc0(j+1, k+1, i,   n) + dqc0(j+1, k+1, i+1, n);
                    }
                }
            }
        } else {
            int ii = 1;
            int i  = 1;
            int kk = 0;
            for (int k = 1; k <= kdim1; k += 2) {
                kk = kk + 1;
                int jj = 0;
                for (int j = 1; j <= jdim1; j += 2) {
                    jj = jj + 1;
                    dqc0c(jj, kk, ii, n) =
                        dqc0(j,   k,   i, n) + dqc0(j+1, k,   i, n)
                      + dqc0(j,   k+1, i, n) + dqc0(j+1, k+1, i, n);
                }
            }
        }
    }
}

} // namespace collqc0_ns
