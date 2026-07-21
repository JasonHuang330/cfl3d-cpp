// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "arclen.h"
#include "arc.h"
#include <cmath>
#include <algorithm>

namespace arclen_ns {

void arc(int& jdim1, int& kdim1, int& msub1, int& msub2,
         FortranArray1DRef<int> jjmax1, FortranArray1DRef<int> kkmax1,
         int& lmax1,
         FortranArray3DRef<double> x1, FortranArray3DRef<double> y1, FortranArray3DRef<double> z1,
         int& limit0,
         FortranArray1DRef<int> jjmax2, FortranArray1DRef<int> kkmax2,
         FortranArray3DRef<double> x2, FortranArray3DRef<double> y2, FortranArray3DRef<double> z2,
         FortranArray1DRef<double> xie2, FortranArray1DRef<double> eta2,
         FortranArray1DRef<int> mblkpt,
         int& ifit, int& itmax, int& jcorr, int& kcorr,
         FortranArray3DRef<double> sxie, FortranArray3DRef<double> seta,
         FortranArray3DRef<double> sxie2, FortranArray3DRef<double> seta2,
         FortranArray2DRef<double> xie2s, FortranArray2DRef<double> eta2s,
         int& intmx, int& icheck,
         FortranArray1DRef<int> nblkj, FortranArray1DRef<int> nblkk,
         FortranArray1DRef<int> jmm, FortranArray1DRef<int> kmm,
         int& j21, int& j22, int& k21, int& k22, int& npt,
         FortranArray1DRef<int> xif1, FortranArray1DRef<int> xif2,
         FortranArray1DRef<int> etf1, FortranArray1DRef<int> etf2,
         FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
         int& nbuf, int& ibufdim,
         FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    arc_ns::arc(jdim1, kdim1, msub1, msub2, jjmax1, kkmax1, lmax1,
                x1, y1, z1, limit0, jjmax2, kkmax2, x2, y2, z2,
                xie2, eta2, mblkpt, ifit, itmax, jcorr, kcorr,
                sxie, seta, sxie2, seta2, xie2s, eta2s,
                intmx, icheck, nblkj, nblkk, jmm, kmm,
                j21, j22, k21, k22, npt,
                xif1, xif2, etf1, etf2,
                nou, bou, nbuf, ibufdim, mblk2nd, maxbl);
}

void arclen(int& idim, int& jdim, int& kdim,
            FortranArray3DRef<double> arci, FortranArray3DRef<double> arcj,
            FortranArray3DRef<double> arck,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y,
            FortranArray3DRef<double> z,
            int& nbl,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid)
{
    int& iexp    = cmn_zero.iexp;
    int& isklton = cmn_sklton.isklton;

    // tolerance for switch to linear blending function
    // (10.**(-iexp) is machine zero)
    float tol = std::max(1.e-07f, (float)std::pow(10.0, (double)(-iexp + 1)));

    // arc length measured from i=1 surface
    for (int j = 1; j <= jdim; j++) {
        for (int k = 1; k <= kdim; k++) {
            arci(j, k, 1) = 0.0;
        }
    }
    int nlinear = 0;
    for (int j = 1; j <= jdim; j++) {
        for (int k = 1; k <= kdim; k++) {
            for (int i = 2; i <= idim; i++) {
                double ds = std::sqrt(
                    (x(j,k,i) - x(j,k,i-1)) * (x(j,k,i) - x(j,k,i-1)) +
                    (y(j,k,i) - y(j,k,i-1)) * (y(j,k,i) - y(j,k,i-1)) +
                    (z(j,k,i) - z(j,k,i-1)) * (z(j,k,i) - z(j,k,i-1)));
                arci(j, k, i) = arci(j, k, i-1) + ds;
            }
            if ((float)arci(j, k, idim) < (float)tol) {
                nlinear = nlinear + 1;
                for (int i = 2; i <= idim; i++) {
                    arci(j, k, i) = (float)(i - 1);
                }
            }
        }
    }
    if (nlinear > 0 && isklton == 1) {
//         nou(1) = min(nou(1)+1,ibufdim)
//         write(bou(nou(1),1),'('' linear blending functions''
//    .      '' in i-direction, block '',i3)') nlinear,nbl
    }

    // arc length measured from j=1 surface
    for (int i = 1; i <= idim; i++) {
        for (int k = 1; k <= kdim; k++) {
            arcj(1, k, i) = 0.0;
        }
    }
    nlinear = 0;
    for (int i = 1; i <= idim; i++) {
        for (int k = 1; k <= kdim; k++) {
            for (int j = 2; j <= jdim; j++) {
                double ds = std::sqrt(
                    (x(j,k,i) - x(j-1,k,i)) * (x(j,k,i) - x(j-1,k,i)) +
                    (y(j,k,i) - y(j-1,k,i)) * (y(j,k,i) - y(j-1,k,i)) +
                    (z(j,k,i) - z(j-1,k,i)) * (z(j,k,i) - z(j-1,k,i)));
                arcj(j, k, i) = arcj(j-1, k, i) + ds;
            }
            if ((float)arcj(jdim, k, i) < (float)tol) {
                nlinear = nlinear + 1;
                for (int j = 2; j <= jdim; j++) {
                    arcj(j, k, i) = (float)(j - 1);
                }
            }
        }
    }
    if (nlinear > 0 && isklton == 1) {
//         nou(1) = min(nou(1)+1,ibufdim)
//         write(bou(nou(1),1),'('' linear blending functions''
//    .      '' in j-direction, block '',i3)') nlinear,nbl
    }

    // arc length measured from k=1 surface
    for (int j = 1; j <= jdim; j++) {
        for (int i = 1; i <= idim; i++) {
            arck(j, 1, i) = 0.0;
        }
    }
    nlinear = 0;
    for (int j = 1; j <= jdim; j++) {
        for (int i = 1; i <= idim; i++) {
            for (int k = 2; k <= kdim; k++) {
                double ds = std::sqrt(
                    (x(j,k,i) - x(j,k-1,i)) * (x(j,k,i) - x(j,k-1,i)) +
                    (y(j,k,i) - y(j,k-1,i)) * (y(j,k,i) - y(j,k-1,i)) +
                    (z(j,k,i) - z(j,k-1,i)) * (z(j,k,i) - z(j,k-1,i)));
                arck(j, k, i) = arck(j, k-1, i) + ds;
            }
            if ((float)arck(j, kdim, i) < (float)tol) {
                nlinear = nlinear + 1;
                for (int k = 2; k <= kdim; k++) {
                    arck(j, k, i) = (float)(k - 1);
                }
            }
        }
    }
    if (nlinear > 0 && isklton == 1) {
//         nou(1) = min(nou(1)+1,ibufdim)
//         write(bou(nou(1),1),'('' linear blending functions''
//    .      '' in k-direction, block '',i3)') nlinear,nbl
    }

    // normalize arc lengths
    for (int i = 1; i <= idim; i++) {
        for (int j = 1; j <= jdim; j++) {
            for (int k = 1; k <= kdim; k++) {
                arci(j, k, i) = arci(j, k, i) / arci(j, k, idim);
                arcj(j, k, i) = arcj(j, k, i) / arcj(jdim, k, i);
                arck(j, k, i) = arck(j, k, i) / arck(j, kdim, i);
            }
        }
    }
}

} // namespace arclen_ns
