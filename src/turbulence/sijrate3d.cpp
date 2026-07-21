// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "sijrate3d.h"

namespace sijrate3d_ns {

void sijrate3d(int& idim, int& jdim, int& kdim,
               FortranArray4DRef<double> q,
               FortranArray4DRef<double> ux,
               FortranArray3DRef<double> vol,
               FortranArray4DRef<double> si,
               FortranArray4DRef<double> sj,
               FortranArray4DRef<double> sk,
               FortranArray4DRef<double> vx)
{
    // Local variables
    int i, j, k, m;
    double xc, yc, zc, tc, uu, fac;
    double s11p, s22p, s33p, s12p, s13p, s23p;
    double s11m, s22m, s33m, s12m, s13m, s23m;

    // initialize
    for (m = 1; m <= 6; m++) {
        for (i = 1; i <= idim-1; i++) {
            for (k = 1; k <= kdim-1; k++) {
                for (j = 1; j <= jdim-1; j++) {
                    vx(j,k,i,m) = 0.;
                }
            }
        }
    }

    // Note:
    //   s11 = ux(j,k,i,1)
    //   s22 = ux(j,k,i,5)
    //   s33 = ux(j,k,i,9)
    //   s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
    //   s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
    //   s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))

    // j-direction:
    if (jdim > 2) {
        for (i = 1; i <= idim-1; i++) {
            for (k = 1; k <= kdim-1; k++) {
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(sj(j  ,k,i,1)*sj(j  ,k,i,4) +
                               sj(j+1,k,i,1)*sj(j+1,k,i,4)) /
                               vol(j,k,i);
                    yc = 0.5*(sj(j  ,k,i,2)*sj(j  ,k,i,4) +
                               sj(j+1,k,i,2)*sj(j+1,k,i,4)) /
                               vol(j,k,i);
                    zc = 0.5*(sj(j  ,k,i,3)*sj(j  ,k,i,4) +
                               sj(j+1,k,i,3)*sj(j+1,k,i,4)) /
                               vol(j,k,i);
                    tc = 0.5*(sj(j  ,k,i,5)*sj(j  ,k,i,4) +
                               sj(j+1,k,i,5)*sj(j+1,k,i,4)) /
                               vol(j,k,i);
                    uu = xc*q(j,k,i,2) + yc*q(j,k,i,3) + zc*q(j,k,i,4) + tc;
                    if (j >= 2 && j <= jdim-2) {
                        fac = 2.;
                        s11p = ux(j+1,k,i,1);
                        s22p = ux(j+1,k,i,5);
                        s33p = ux(j+1,k,i,9);
                        s12p = 0.5*(ux(j+1,k,i,2) + ux(j+1,k,i,4));
                        s13p = 0.5*(ux(j+1,k,i,3) + ux(j+1,k,i,7));
                        s23p = 0.5*(ux(j+1,k,i,6) + ux(j+1,k,i,8));
                        s11m = ux(j-1,k,i,1);
                        s22m = ux(j-1,k,i,5);
                        s33m = ux(j-1,k,i,9);
                        s12m = 0.5*(ux(j-1,k,i,2) + ux(j-1,k,i,4));
                        s13m = 0.5*(ux(j-1,k,i,3) + ux(j-1,k,i,7));
                        s23m = 0.5*(ux(j-1,k,i,6) + ux(j-1,k,i,8));
                    } else if (j == 1) {
                        fac = 1.;
                        s11p = ux(j+1,k,i,1);
                        s22p = ux(j+1,k,i,5);
                        s33p = ux(j+1,k,i,9);
                        s12p = 0.5*(ux(j+1,k,i,2) + ux(j+1,k,i,4));
                        s13p = 0.5*(ux(j+1,k,i,3) + ux(j+1,k,i,7));
                        s23p = 0.5*(ux(j+1,k,i,6) + ux(j+1,k,i,8));
                        s11m = ux(j  ,k,i,1);
                        s22m = ux(j  ,k,i,5);
                        s33m = ux(j  ,k,i,9);
                        s12m = 0.5*(ux(j  ,k,i,2) + ux(j  ,k,i,4));
                        s13m = 0.5*(ux(j  ,k,i,3) + ux(j  ,k,i,7));
                        s23m = 0.5*(ux(j  ,k,i,6) + ux(j  ,k,i,8));
                    } else if (j == jdim-1) {
                        fac = 1.;
                        s11p = ux(j  ,k,i,1);
                        s22p = ux(j  ,k,i,5);
                        s33p = ux(j  ,k,i,9);
                        s12p = 0.5*(ux(j  ,k,i,2) + ux(j  ,k,i,4));
                        s13p = 0.5*(ux(j  ,k,i,3) + ux(j  ,k,i,7));
                        s23p = 0.5*(ux(j  ,k,i,6) + ux(j  ,k,i,8));
                        s11m = ux(j-1,k,i,1);
                        s22m = ux(j-1,k,i,5);
                        s33m = ux(j-1,k,i,9);
                        s12m = 0.5*(ux(j-1,k,i,2) + ux(j-1,k,i,4));
                        s13m = 0.5*(ux(j-1,k,i,3) + ux(j-1,k,i,7));
                        s23m = 0.5*(ux(j-1,k,i,6) + ux(j-1,k,i,8));
                    }
                    vx(j,k,i,1) = vx(j,k,i,1) + uu*(s11p-s11m)/fac;
                    vx(j,k,i,2) = vx(j,k,i,2) + uu*(s12p-s12m)/fac;
                    vx(j,k,i,3) = vx(j,k,i,3) + uu*(s13p-s13m)/fac;
                    vx(j,k,i,4) = vx(j,k,i,4) + uu*(s22p-s22m)/fac;
                    vx(j,k,i,5) = vx(j,k,i,5) + uu*(s23p-s23m)/fac;
                    vx(j,k,i,6) = vx(j,k,i,6) + uu*(s33p-s33m)/fac;
                }
            }
        }
    }

    // k-direction:
    if (kdim > 2) {
        for (i = 1; i <= idim-1; i++) {
            for (j = 1; j <= jdim-1; j++) {
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(sk(j,k  ,i,1)*sk(j,k  ,i,4) +
                               sk(j,k+1,i,1)*sk(j,k+1,i,4)) /
                               vol(j,k,i);
                    yc = 0.5*(sk(j,k  ,i,2)*sk(j,k  ,i,4) +
                               sk(j,k+1,i,2)*sk(j,k+1,i,4)) /
                               vol(j,k,i);
                    zc = 0.5*(sk(j,k  ,i,3)*sk(j,k  ,i,4) +
                               sk(j,k+1,i,3)*sk(j,k+1,i,4)) /
                               vol(j,k,i);
                    tc = 0.5*(sk(j,k  ,i,5)*sk(j,k  ,i,4) +
                               sk(j,k+1,i,5)*sk(j,k+1,i,4)) /
                               vol(j,k,i);
                    uu = xc*q(j,k,i,2) + yc*q(j,k,i,3) + zc*q(j,k,i,4) + tc;
                    if (k >= 2 && k <= kdim-2) {
                        fac = 2.;
                        s11p = ux(j,k+1,i,1);
                        s22p = ux(j,k+1,i,5);
                        s33p = ux(j,k+1,i,9);
                        s12p = 0.5*(ux(j,k+1,i,2) + ux(j,k+1,i,4));
                        s13p = 0.5*(ux(j,k+1,i,3) + ux(j,k+1,i,7));
                        s23p = 0.5*(ux(j,k+1,i,6) + ux(j,k+1,i,8));
                        s11m = ux(j,k-1,i,1);
                        s22m = ux(j,k-1,i,5);
                        s33m = ux(j,k-1,i,9);
                        s12m = 0.5*(ux(j,k-1,i,2) + ux(j,k-1,i,4));
                        s13m = 0.5*(ux(j,k-1,i,3) + ux(j,k-1,i,7));
                        s23m = 0.5*(ux(j,k-1,i,6) + ux(j,k-1,i,8));
                    } else if (k == 1) {
                        fac = 1.;
                        s11p = ux(j,k+1,i,1);
                        s22p = ux(j,k+1,i,5);
                        s33p = ux(j,k+1,i,9);
                        s12p = 0.5*(ux(j,k+1,i,2) + ux(j,k+1,i,4));
                        s13p = 0.5*(ux(j,k+1,i,3) + ux(j,k+1,i,7));
                        s23p = 0.5*(ux(j,k+1,i,6) + ux(j,k+1,i,8));
                        s11m = ux(j,k  ,i,1);
                        s22m = ux(j,k  ,i,5);
                        s33m = ux(j,k  ,i,9);
                        s12m = 0.5*(ux(j,k  ,i,2) + ux(j,k  ,i,4));
                        s13m = 0.5*(ux(j,k  ,i,3) + ux(j,k  ,i,7));
                        s23m = 0.5*(ux(j,k  ,i,6) + ux(j,k  ,i,8));
                    } else if (k == kdim-1) {
                        fac = 1.;
                        s11p = ux(j,k  ,i,1);
                        s22p = ux(j,k  ,i,5);
                        s33p = ux(j,k  ,i,9);
                        s12p = 0.5*(ux(j,k  ,i,2) + ux(j,k  ,i,4));
                        s13p = 0.5*(ux(j,k  ,i,3) + ux(j,k  ,i,7));
                        s23p = 0.5*(ux(j,k  ,i,6) + ux(j,k  ,i,8));
                        s11m = ux(j,k-1,i,1);
                        s22m = ux(j,k-1,i,5);
                        s33m = ux(j,k-1,i,9);
                        s12m = 0.5*(ux(j,k-1,i,2) + ux(j,k-1,i,4));
                        s13m = 0.5*(ux(j,k-1,i,3) + ux(j,k-1,i,7));
                        s23m = 0.5*(ux(j,k-1,i,6) + ux(j,k-1,i,8));
                    }
                    vx(j,k,i,1) = vx(j,k,i,1) + uu*(s11p-s11m)/fac;
                    vx(j,k,i,2) = vx(j,k,i,2) + uu*(s12p-s12m)/fac;
                    vx(j,k,i,3) = vx(j,k,i,3) + uu*(s13p-s13m)/fac;
                    vx(j,k,i,4) = vx(j,k,i,4) + uu*(s22p-s22m)/fac;
                    vx(j,k,i,5) = vx(j,k,i,5) + uu*(s23p-s23m)/fac;
                    vx(j,k,i,6) = vx(j,k,i,6) + uu*(s33p-s33m)/fac;
                }
            }
        }
    }

    // i-direction:
    if (idim > 2) {
        for (k = 1; k <= kdim-1; k++) {
            for (j = 1; j <= jdim-1; j++) {
                for (i = 1; i <= idim-1; i++) {
                    xc = 0.5*(si(j,k,i  ,1)*si(j,k,i  ,4) +
                               si(j,k,i+1,1)*si(j,k,i+1,4)) /
                               vol(j,k,i);
                    yc = 0.5*(si(j,k,i  ,2)*si(j,k,i  ,4) +
                               si(j,k,i+1,2)*si(j,k,i+1,4)) /
                               vol(j,k,i);
                    zc = 0.5*(si(j,k,i  ,3)*si(j,k,i  ,4) +
                               si(j,k,i+1,3)*si(j,k,i+1,4)) /
                               vol(j,k,i);
                    tc = 0.5*(si(j,k,i  ,5)*si(j,k,i  ,4) +
                               si(j,k,i+1,5)*si(j,k,i+1,4)) /
                               vol(j,k,i);
                    uu = xc*q(j,k,i,2) + yc*q(j,k,i,3) + zc*q(j,k,i,4) + tc;
                    if (i >= 2 && i <= idim-2) {
                        fac = 2.;
                        s11p = ux(j,k,i+1,1);
                        s22p = ux(j,k,i+1,5);
                        s33p = ux(j,k,i+1,9);
                        s12p = 0.5*(ux(j,k,i+1,2) + ux(j,k,i+1,4));
                        s13p = 0.5*(ux(j,k,i+1,3) + ux(j,k,i+1,7));
                        s23p = 0.5*(ux(j,k,i+1,6) + ux(j,k,i+1,8));
                        s11m = ux(j,k,i-1,1);
                        s22m = ux(j,k,i-1,5);
                        s33m = ux(j,k,i-1,9);
                        s12m = 0.5*(ux(j,k,i-1,2) + ux(j,k,i-1,4));
                        s13m = 0.5*(ux(j,k,i-1,3) + ux(j,k,i-1,7));
                        s23m = 0.5*(ux(j,k,i-1,6) + ux(j,k,i-1,8));
                    } else if (i == 1) {
                        fac = 1.;
                        s11p = ux(j,k,i+1,1);
                        s22p = ux(j,k,i+1,5);
                        s33p = ux(j,k,i+1,9);
                        s12p = 0.5*(ux(j,k,i+1,2) + ux(j,k,i+1,4));
                        s13p = 0.5*(ux(j,k,i+1,3) + ux(j,k,i+1,7));
                        s23p = 0.5*(ux(j,k,i+1,6) + ux(j,k,i+1,8));
                        s11m = ux(j,k,i  ,1);
                        s22m = ux(j,k,i  ,5);
                        s33m = ux(j,k,i  ,9);
                        s12m = 0.5*(ux(j,k,i  ,2) + ux(j,k,i  ,4));
                        s13m = 0.5*(ux(j,k,i  ,3) + ux(j,k,i  ,7));
                        s23m = 0.5*(ux(j,k,i  ,6) + ux(j,k,i  ,8));
                    } else if (i == idim-1) {
                        fac = 1.;
                        s11p = ux(j,k,i  ,1);
                        s22p = ux(j,k,i  ,5);
                        s33p = ux(j,k,i  ,9);
                        s12p = 0.5*(ux(j,k,i  ,2) + ux(j,k,i  ,4));
                        s13p = 0.5*(ux(j,k,i  ,3) + ux(j,k,i  ,7));
                        s23p = 0.5*(ux(j,k,i  ,6) + ux(j,k,i  ,8));
                        s11m = ux(j,k,i-1,1);
                        s22m = ux(j,k,i-1,5);
                        s33m = ux(j,k,i-1,9);
                        s12m = 0.5*(ux(j,k,i-1,2) + ux(j,k,i-1,4));
                        s13m = 0.5*(ux(j,k,i-1,3) + ux(j,k,i-1,7));
                        s23m = 0.5*(ux(j,k,i-1,6) + ux(j,k,i-1,8));
                    }
                    vx(j,k,i,1) = vx(j,k,i,1) + uu*(s11p-s11m)/fac;
                    vx(j,k,i,2) = vx(j,k,i,2) + uu*(s12p-s12m)/fac;
                    vx(j,k,i,3) = vx(j,k,i,3) + uu*(s13p-s13m)/fac;
                    vx(j,k,i,4) = vx(j,k,i,4) + uu*(s22p-s22m)/fac;
                    vx(j,k,i,5) = vx(j,k,i,5) + uu*(s23p-s23m)/fac;
                    vx(j,k,i,6) = vx(j,k,i,6) + uu*(s33p-s33m)/fac;
                }
            }
        }
    }

    return;
}

} // namespace sijrate3d_ns
