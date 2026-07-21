// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rechk.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace rechk_ns {

void rechk(int& mdim1, int& ndim1, FortranArray3DRef<int> jimage, FortranArray3DRef<int> kimage, int& msub1, int& jmax1, int& kmax1, int& l, FortranArray2DRef<double> x1, FortranArray2DRef<double> y1, FortranArray2DRef<double> z1, int& xie1, int& xie2, int& eta1, int& eta2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    // COMMON /tol/ epsc,epsc0,eps,epscoll
    // eps is renamed to epsreen in cmn_tol
    float& epsc     = cmn_tol.epsc;
    float& epsc0    = cmn_tol.epsc0;
    float& eps      = cmn_tol.epsreen;
    float& epscoll  = cmn_tol.epscoll;

    // COMMON /sklt1/ isklt1
    int32_t& isklt1 = cmn_sklt1.isklt1;

    // local variables
    int js, jq, ks, kq;
    int kcheck, jcheck;
    int jcount, kcount;
    int iflag, iflago;
    int jpr, kpr;
    double delxyz;

    // default to no branch cuts
    for (int j = 1; j <= jmax1+1; j++) {
        for (int k = 1; k <= kmax1+1; k++) {
            jimage(l,j,k) = j;
            kimage(l,j,k) = k;
        }
    }

    // check for C-grid type branch cut on k=1 and k=kmax1 boundaries
    // branch cuts occur only on block edges, so if search range does not
    // include block edges, don't bother to check for cuts
    if (eta1 > 1 && eta2 < kmax1) goto label11;

    for (int k = 1; k <= kmax1; k += kmax1-1) {
        kcheck = 0;
        for (js = xie1; js <= xie2-1; js++) {
            for (jq = js+1; jq <= xie2; jq++) {
                delxyz = std::abs(x1(js,k)-x1(jq,k))
                       + std::abs(y1(js,k)-y1(jq,k))
                       + std::abs(z1(js,k)-z1(jq,k));
                if ((float)delxyz < (float)eps) {
                    kcheck = kcheck+1;
                    if (k == 1) {
                        jimage(l,js+1,1) = jq;
                        kimage(l,js+1,1) = 2;
                        jimage(l,jq,1)   = js+1;
                        kimage(l,jq,1)   = 2;
                    }
                    if (k == kmax1) {
                        jimage(l,js+1,kmax1+1) = jq;
                        kimage(l,js+1,kmax1+1) = kmax1;
                        jimage(l,jq,kmax1+1)   = js+1;
                        kimage(l,jq,kmax1+1)   = kmax1;
                    }
                }
            }
        }
        if (kcheck == 0) continue; // go to 10
        if (k == 1) {
            // correction for border points between start/end
            // of branch cuts on the boundary
            for (int j = xie1+1; j <= xie2-1; j++) {
                if (jimage(l,j+1,1) == j+1) {
                    jimage(l,jimage(l,j,1),1) = jimage(l,j,1);
                    kimage(l,jimage(l,j,1),1) = 1;
                    jimage(l,j,1)             = j;
                    kimage(l,j,1)             = 1;
                    goto label1809;
                }
            }
            label1809:;

            // image point at ends of boundary of extended grid
            if (jimage(l,xie1+1,1) == xie2) {
                jimage(l,xie1,1)   = xie2+1;
                kimage(l,xie1,1)   = 2;
                jimage(l,xie2+1,1) = xie1;
                kimage(l,xie2+1,1) = 2;
            }
        }

        if (k == kmax1) {
            // correction for border points between start/end
            // of branch cuts on the boundary
            for (int j = xie1+1; j <= xie2-1; j++) {
                if (jimage(l,j+1,kmax1+1) == j+1) {
                    jimage(l,jimage(l,j,kmax1+1),kmax1+1) = jimage(l,j,kmax1+1);
                    kimage(l,jimage(l,j,kmax1+1),kmax1+1) = kmax1+1;
                    jimage(l,j,kmax1+1)                   = j;
                    kimage(l,j,kmax1+1)                   = kmax1+1;
                    goto label1808;
                }
            }
            label1808:;

            // image point at left boundary of extended grid
            if (jimage(l,xie1+1,kmax1+1) == xie2) {
                jimage(l,xie1,kmax1+1)   = xie2+1;
                kimage(l,xie1,kmax1+1)   = kmax1;
                jimage(l,xie2+1,kmax1+1) = xie1;
                kimage(l,xie2+1,kmax1+1) = kmax1;
            }
        }
        // end of do 10 loop (k loop)
    }

    label11:;

    // check for C-grid type branch cut on j=1 and j=jmax1 boundaries
    // branch cuts occur only on block edges, so if search range does not
    // include block edges, don't bother to check for cuts
    if (xie1 > 1 && xie2 < jmax1) goto label21;

    for (int j = 1; j <= jmax1; j += jmax1-1) {
        jcheck = 0;
        for (ks = eta1; ks <= eta2-1; ks++) {
            for (kq = ks+1; kq <= eta2; kq++) {
                delxyz = std::abs(x1(j,ks)-x1(j,kq))
                       + std::abs(y1(j,ks)-y1(j,kq))
                       + std::abs(z1(j,ks)-z1(j,kq));
                if ((float)delxyz < (float)eps) {
                    jcheck = jcheck+1;
                    if (j == 1) {
                        jimage(l,1,ks+1) = 2;
                        kimage(l,1,ks+1) = kq;
                        jimage(l,1,kq)   = 2;
                        kimage(l,1,kq)   = ks+1;
                    }
                    if (j == jmax1) {
                        jimage(l,jmax1+1,ks+1) = jmax1;
                        kimage(l,jmax1+1,ks+1) = kq;
                        jimage(l,jmax1+1,kq)   = jmax1;
                        kimage(l,jmax1+1,kq)   = ks+1;
                    }
                    goto label2712; // go to 2712 (next ks iteration)
                }
            }
            label2712:;
        }
        if (jcheck == 0) continue; // go to 20
        if (j == 1) {
            // correction for border points between start/end
            // of branch cuts on the boundary
            for (int k = eta1+1; k <= eta2-1; k++) {
                if (kimage(l,1,k+1) == k+1) {
                    jimage(l,1,kimage(l,1,k)) = 1;
                    kimage(l,1,kimage(l,1,k)) = kimage(l,1,k);
                    jimage(l,1,k)             = 1;
                    kimage(l,1,k)             = k;
                    goto label2809;
                }
            }
            label2809:;

            // image points at ends of boundary of extended grid
            if (kimage(l,1,eta1+1) == eta2) {
                jimage(l,1,eta1)   = 2;
                kimage(l,1,eta1)   = eta2+1;
                jimage(l,1,eta2+1) = 2;
                kimage(l,1,eta2+1) = eta1;
            }
        }

        if (j == jmax1) {
            // correction for border points between start/end
            // of branch cuts on the boundary
            for (int k = eta1+1; k <= eta2-1; k++) {
                if (kimage(l,jmax1+1,k+1) == k+1) {
                    jimage(l,jmax1+1,kimage(l,jmax1+1,k)) = jmax1+1;
                    kimage(l,jmax1+1,kimage(l,jmax1+1,k)) = kimage(l,jmax1+1,k);
                    jimage(l,jmax1+1,k)                   = jmax1+1;
                    kimage(l,jmax1+1,k)                   = k;
                    goto label2808;
                }
            }
            label2808:;

            // image point at bottom boundary of extended grid
            if (kimage(l,jmax1+1,eta1+1) == eta2) {
                jimage(l,jmax1+1,eta1)   = jmax1;
                kimage(l,jmax1+1,eta1)   = eta2+1;
                jimage(l,jmax1+1,eta2+1) = jmax1;
                kimage(l,jmax1+1,eta2+1) = eta1;
            }
        }
        // end of do 20 loop (j loop)
    }

    label21:;

    // check for O-grid type branch cut on k=1 and k=kmax1 boundaries
    jcount = 0;

    // branch cuts occur only on block edges, so if search range does not
    // include block edges, don't bother to check for cuts
    if (eta1 > 1 && eta2 < kmax1) goto label31;

    for (int j = xie1; j <= xie2; j++) {
        delxyz = std::abs(x1(j,1)-x1(j,kmax1))
               + std::abs(y1(j,1)-y1(j,kmax1))
               + std::abs(z1(j,1)-z1(j,kmax1));
        if ((float)delxyz < (float)eps) jcount = jcount+1;
    }
    if (jcount > 1 && jcount != xie2-xie1+1) {
        if (isklt1 > 1) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                " WARNING...O-type branch cut does"
                " not extend over entire k=constant  boundary, as"
                " assumed");
        }
    }
    if (jcount > 1) {
        for (int j = xie1; j <= xie2+1; j++) {
            jimage(l,j,1)       = j;
            kimage(l,j,1)       = kmax1;
            jimage(l,j,kmax1+1) = j;
            kimage(l,j,kmax1+1) = 2;
        }
    }

    label31:;

    // check for O-grid type branch cut on j=1 and j=jmax1 boundaries
    kcount = 0;

    // branch cuts occur only on block edges, so if search range does not
    // include block edges, don't bother to check for cuts
    if (xie1 > 1 && xie2 < jmax1) goto label41;

    for (int k = eta1; k <= eta2; k++) {
        delxyz = std::abs(x1(1,k)-x1(jmax1,k))
               + std::abs(y1(1,k)-y1(jmax1,k))
               + std::abs(z1(1,k)-z1(jmax1,k));
        if ((float)delxyz < (float)eps) kcount = kcount+1;
    }
    if (kcount > 1 && kcount != eta2-eta1+1) {
        if (isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            std::snprintf(bou(nou(4),4), 120,
                " WARNING...O-type branch cut does"
                " not extend over entire k=constant  boundary, as"
                " assumed");
        }
    }
    if (kcount > 1) {
        for (int k = eta1; k <= eta2+1; k++) {
            jimage(l,1,k)       = jmax1;
            kimage(l,1,k)       = k;
            jimage(l,jmax1+1,k) = 2;
            kimage(l,jmax1+1,k) = k;
        }
    }

    label41:;

    // print out branch cut topology
    // do 50 j=xie1,xie2+1,xie2+1-xie1
    for (int j = xie1; j <= xie2+1; j += (xie2+1-xie1)) {
        iflag = 0;
        jpr   = j;
        if (j == jmax1+1) jpr = jmax1;
        for (int k = eta1+1; k <= eta2; k++) {
            if (jimage(l,j,k) != j || kimage(l,j,k) != k) iflag = 1;
        }
        iflago = 0;
        for (int k = eta1+1; k <= eta2; k++) {
            if (jimage(l,j,k) == j || kimage(l,j,k) != k) iflago = 1;
        }
        if (iflag != 0 && iflago == 0 && isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            // 101 format(' ','     j=',i3,' is an O-type (periodic)',' branch cut boundary')
            std::snprintf(bou(nou(4),4), 120,
                "      j=%3d is an O-type (periodic) branch cut boundary", jpr);
        }
        if (iflag != 0 && iflago != 0 && isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            // 102 format(' ','     j=',i3,' is a C-type branch cut boundary')
            std::snprintf(bou(nou(4),4), 120,
                "      j=%3d is a C-type branch cut boundary", jpr);
        }
    }

    // do 60 k=eta1,eta2+1,eta2+1-eta1
    for (int k = eta1; k <= eta2+1; k += (eta2+1-eta1)) {
        kpr = k;
        if (k == kmax1+1) kpr = kmax1;
        iflag = 0;
        for (int j = xie1+1; j <= xie2; j++) {
            if (jimage(l,j,k) != j || kimage(l,j,k) != k) iflag = 1;
        }
        iflago = 0;
        for (int j = xie1+1; j <= xie2; j++) {
            if (jimage(l,j,k) != j || kimage(l,j,k) == k) iflago = 1;
        }
        if (iflag != 0 && iflago == 0 && isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            // 201 format(' ','     k=',i3,' is an O-type (periodic)',' branch cut boundary')
            std::snprintf(bou(nou(4),4), 120,
                "      k=%3d is an O-type (periodic) branch cut boundary", kpr);
        }
        if (iflag != 0 && iflago != 0 && isklt1 > 0) {
            nou(4) = std::min(nou(4)+1, ibufdim);
            // 202 format(' ','     k=',i3,' is a C-type branch cut boundary')
            std::snprintf(bou(nou(4),4), 120,
                "      k=%3d is a C-type branch cut boundary", kpr);
        }
    }

    return;
}

} // namespace rechk_ns
