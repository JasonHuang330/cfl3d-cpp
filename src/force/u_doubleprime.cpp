// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "u_doubleprime.h"
#include <cmath>

namespace u_doubleprime_ns {

void u_doubleprime(int& idim, int& jdim, int& kdim, FortranArray4DRef<double> q, FortranArray4DRef<double> ux, FortranArray3DRef<double> vol, FortranArray4DRef<double> si, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> vx)
{
    int m, i, j, k;
    double xc, yc, zc, fac;
    double dudxp, dudyp, dudzp;
    double dvdxp, dvdyp, dvdzp;
    double dwdxp, dwdyp, dwdzp;
    double dudxm, dudym, dudzm;
    double dvdxm, dvdym, dvdzm;
    double dwdxm, dwdym, dwdzm;

    // initialize
    for (m = 1; m <= 4; m++) {
        for (i = 1; i <= idim-1; i++) {
            for (k = 1; k <= kdim-1; k++) {
                for (j = 1; j <= jdim-1; j++) {
                    vx(j,k,i,m) = 0.;
                }
            }
        }
    }

    // j-direction:
    if (jdim > 2) {
        for (i = 1; i <= idim-1; i++) {
            for (k = 1; k <= kdim-1; k++) {
                for (j = 1; j <= jdim-1; j++) {
                    xc = 0.5*(sj(j  ,k  ,i  ,1)*sj(j  ,k  ,i  ,4)+
                              sj(j+1,k  ,i  ,1)*sj(j+1,k  ,i  ,4))/
                              vol(j,k,i);
                    yc = 0.5*(sj(j  ,k  ,i  ,2)*sj(j  ,k  ,i  ,4)+
                              sj(j+1,k  ,i  ,2)*sj(j+1,k  ,i  ,4))/
                              vol(j,k,i);
                    zc = 0.5*(sj(j  ,k  ,i  ,3)*sj(j  ,k  ,i  ,4)+
                              sj(j+1,k  ,i  ,3)*sj(j+1,k  ,i  ,4))/
                              vol(j,k,i);
                    if (j >= 2 && j <= jdim-2) {
                        fac = 2.;
                        dudxp = ux(j+1,k  ,i  ,1);
                        dudyp = ux(j+1,k  ,i  ,2);
                        dudzp = ux(j+1,k  ,i  ,3);
                        dvdxp = ux(j+1,k  ,i  ,4);
                        dvdyp = ux(j+1,k  ,i  ,5);
                        dvdzp = ux(j+1,k  ,i  ,6);
                        dwdxp = ux(j+1,k  ,i  ,7);
                        dwdyp = ux(j+1,k  ,i  ,8);
                        dwdzp = ux(j+1,k  ,i  ,9);
                        dudxm = ux(j-1,k  ,i  ,1);
                        dudym = ux(j-1,k  ,i  ,2);
                        dudzm = ux(j-1,k  ,i  ,3);
                        dvdxm = ux(j-1,k  ,i  ,4);
                        dvdym = ux(j-1,k  ,i  ,5);
                        dvdzm = ux(j-1,k  ,i  ,6);
                        dwdxm = ux(j-1,k  ,i  ,7);
                        dwdym = ux(j-1,k  ,i  ,8);
                        dwdzm = ux(j-1,k  ,i  ,9);
                    } else if (j == 1) {
                        fac = 1.;
                        dudxp = ux(j+1,k  ,i  ,1);
                        dudyp = ux(j+1,k  ,i  ,2);
                        dudzp = ux(j+1,k  ,i  ,3);
                        dvdxp = ux(j+1,k  ,i  ,4);
                        dvdyp = ux(j+1,k  ,i  ,5);
                        dvdzp = ux(j+1,k  ,i  ,6);
                        dwdxp = ux(j+1,k  ,i  ,7);
                        dwdyp = ux(j+1,k  ,i  ,8);
                        dwdzp = ux(j+1,k  ,i  ,9);
                        dudxm = ux(j  ,k  ,i  ,1);
                        dudym = ux(j  ,k  ,i  ,2);
                        dudzm = ux(j  ,k  ,i  ,3);
                        dvdxm = ux(j  ,k  ,i  ,4);
                        dvdym = ux(j  ,k  ,i  ,5);
                        dvdzm = ux(j  ,k  ,i  ,6);
                        dwdxm = ux(j  ,k  ,i  ,7);
                        dwdym = ux(j  ,k  ,i  ,8);
                        dwdzm = ux(j  ,k  ,i  ,9);
                    } else if (j == jdim-1) {
                        fac = 1.;
                        dudxp = ux(j  ,k  ,i  ,1);
                        dudyp = ux(j  ,k  ,i  ,2);
                        dudzp = ux(j  ,k  ,i  ,3);
                        dvdxp = ux(j  ,k  ,i  ,4);
                        dvdyp = ux(j  ,k  ,i  ,5);
                        dvdzp = ux(j  ,k  ,i  ,6);
                        dwdxp = ux(j  ,k  ,i  ,7);
                        dwdyp = ux(j  ,k  ,i  ,8);
                        dwdzp = ux(j  ,k  ,i  ,9);
                        dudxm = ux(j-1,k  ,i  ,1);
                        dudym = ux(j-1,k  ,i  ,2);
                        dudzm = ux(j-1,k  ,i  ,3);
                        dvdxm = ux(j-1,k  ,i  ,4);
                        dvdym = ux(j-1,k  ,i  ,5);
                        dvdzm = ux(j-1,k  ,i  ,6);
                        dwdxm = ux(j-1,k  ,i  ,7);
                        dwdym = ux(j-1,k  ,i  ,8);
                        dwdzm = ux(j-1,k  ,i  ,9);
                    }
                    vx(j,k,i,2) = vx(j,k,i,2) + (xc*(dudxp-dudxm)
                                                 + yc*(dudyp-dudym)
                                                 + zc*(dudzp-dudzm))/fac;
                    vx(j,k,i,3) = vx(j,k,i,3) + (xc*(dvdxp-dvdxm)
                                                 + yc*(dvdyp-dvdym)
                                                 + zc*(dvdzp-dvdzm))/fac;
                    vx(j,k,i,4) = vx(j,k,i,4) + (xc*(dwdxp-dwdxm)
                                                 + yc*(dwdyp-dwdym)
                                                 + zc*(dwdzp-dwdzm))/fac;
                }
            }
        }
    }

    // k-direction:
    if (kdim > 2) {
        for (i = 1; i <= idim-1; i++) {
            for (j = 1; j <= jdim-1; j++) {
                for (k = 1; k <= kdim-1; k++) {
                    xc = 0.5*(sk(j  ,k  ,i  ,1)*sk(j  ,k  ,i  ,4)+
                              sk(j  ,k+1,i  ,1)*sk(j  ,k+1,i  ,4))/
                              vol(j,k,i);
                    yc = 0.5*(sk(j  ,k  ,i  ,2)*sk(j  ,k  ,i  ,4)+
                              sk(j  ,k+1,i  ,2)*sk(j  ,k+1,i  ,4))/
                              vol(j,k,i);
                    zc = 0.5*(sk(j  ,k  ,i  ,3)*sk(j  ,k  ,i  ,4)+
                              sk(j  ,k+1,i  ,3)*sk(j  ,k+1,i  ,4))/
                              vol(j,k,i);
                    if (k >= 2 && k <= kdim-2) {
                        fac = 2.;
                        dudxp = ux(j  ,k+1,i  ,1);
                        dudyp = ux(j  ,k+1,i  ,2);
                        dudzp = ux(j  ,k+1,i  ,3);
                        dvdxp = ux(j  ,k+1,i  ,4);
                        dvdyp = ux(j  ,k+1,i  ,5);
                        dvdzp = ux(j  ,k+1,i  ,6);
                        dwdxp = ux(j  ,k+1,i  ,7);
                        dwdyp = ux(j  ,k+1,i  ,8);
                        dwdzp = ux(j  ,k+1,i  ,9);
                        dudxm = ux(j  ,k-1,i  ,1);
                        dudym = ux(j  ,k-1,i  ,2);
                        dudzm = ux(j  ,k-1,i  ,3);
                        dvdxm = ux(j  ,k-1,i  ,4);
                        dvdym = ux(j  ,k-1,i  ,5);
                        dvdzm = ux(j  ,k-1,i  ,6);
                        dwdxm = ux(j  ,k-1,i  ,7);
                        dwdym = ux(j  ,k-1,i  ,8);
                        dwdzm = ux(j  ,k-1,i  ,9);
                    } else if (k == 1) {
                        fac = 1.;
                        dudxp = ux(j  ,k+1,i  ,1);
                        dudyp = ux(j  ,k+1,i  ,2);
                        dudzp = ux(j  ,k+1,i  ,3);
                        dvdxp = ux(j  ,k+1,i  ,4);
                        dvdyp = ux(j  ,k+1,i  ,5);
                        dvdzp = ux(j  ,k+1,i  ,6);
                        dwdxp = ux(j  ,k+1,i  ,7);
                        dwdyp = ux(j  ,k+1,i  ,8);
                        dwdzp = ux(j  ,k+1,i  ,9);
                        dudxm = ux(j  ,k  ,i  ,1);
                        dudym = ux(j  ,k  ,i  ,2);
                        dudzm = ux(j  ,k  ,i  ,3);
                        dvdxm = ux(j  ,k  ,i  ,4);
                        dvdym = ux(j  ,k  ,i  ,5);
                        dvdzm = ux(j  ,k  ,i  ,6);
                        dwdxm = ux(j  ,k  ,i  ,7);
                        dwdym = ux(j  ,k  ,i  ,8);
                        dwdzm = ux(j  ,k  ,i  ,9);
                    } else if (k == kdim-1) {
                        fac = 1.;
                        dudxp = ux(j  ,k  ,i  ,1);
                        dudyp = ux(j  ,k  ,i  ,2);
                        dudzp = ux(j  ,k  ,i  ,3);
                        dvdxp = ux(j  ,k  ,i  ,4);
                        dvdyp = ux(j  ,k  ,i  ,5);
                        dvdzp = ux(j  ,k  ,i  ,6);
                        dwdxp = ux(j  ,k  ,i  ,7);
                        dwdyp = ux(j  ,k  ,i  ,8);
                        dwdzp = ux(j  ,k  ,i  ,9);
                        dudxm = ux(j  ,k-1,i  ,1);
                        dudym = ux(j  ,k-1,i  ,2);
                        dudzm = ux(j  ,k-1,i  ,3);
                        dvdxm = ux(j  ,k-1,i  ,4);
                        dvdym = ux(j  ,k-1,i  ,5);
                        dvdzm = ux(j  ,k-1,i  ,6);
                        dwdxm = ux(j  ,k-1,i  ,7);
                        dwdym = ux(j  ,k-1,i  ,8);
                        dwdzm = ux(j  ,k-1,i  ,9);
                    }
                    vx(j,k,i,2) = vx(j,k,i,2) + (xc*(dudxp-dudxm)
                                                 + yc*(dudyp-dudym)
                                                 + zc*(dudzp-dudzm))/fac;
                    vx(j,k,i,3) = vx(j,k,i,3) + (xc*(dvdxp-dvdxm)
                                                 + yc*(dvdyp-dvdym)
                                                 + zc*(dvdzp-dvdzm))/fac;
                    vx(j,k,i,4) = vx(j,k,i,4) + (xc*(dwdxp-dwdxm)
                                                 + yc*(dwdyp-dwdym)
                                                 + zc*(dwdzp-dwdzm))/fac;
                }
            }
        }
    }

    // i-direction:
    if (idim > 2) {
        for (k = 1; k <= kdim-1; k++) {
            for (j = 1; j <= jdim-1; j++) {
                for (i = 1; i <= idim-1; i++) {
                    xc = 0.5*(si(j  ,k  ,i  ,1)*si(j  ,k  ,i  ,4)+
                              si(j  ,k  ,i+1,1)*si(j  ,k  ,i+1,4))/
                              vol(j,k,i);
                    yc = 0.5*(si(j  ,k  ,i  ,2)*si(j  ,k  ,i  ,4)+
                              si(j  ,k  ,i+1,2)*si(j  ,k  ,i+1,4))/
                              vol(j,k,i);
                    zc = 0.5*(si(j  ,k  ,i  ,3)*si(j  ,k  ,i  ,4)+
                              si(j  ,k  ,i+1,3)*si(j  ,k  ,i+1,4))/
                              vol(j,k,i);
                    if (i >= 2 && i <= idim-2) {
                        fac = 2.;
                        dudxp = ux(j  ,k  ,i+1,1);
                        dudyp = ux(j  ,k  ,i+1,2);
                        dudzp = ux(j  ,k  ,i+1,3);
                        dvdxp = ux(j  ,k  ,i+1,4);
                        dvdyp = ux(j  ,k  ,i+1,5);
                        dvdzp = ux(j  ,k  ,i+1,6);
                        dwdxp = ux(j  ,k  ,i+1,7);
                        dwdyp = ux(j  ,k  ,i+1,8);
                        dwdzp = ux(j  ,k  ,i+1,9);
                        dudxm = ux(j  ,k  ,i-1,1);
                        dudym = ux(j  ,k  ,i-1,2);
                        dudzm = ux(j  ,k  ,i-1,3);
                        dvdxm = ux(j  ,k  ,i-1,4);
                        dvdym = ux(j  ,k  ,i-1,5);
                        dvdzm = ux(j  ,k  ,i-1,6);
                        dwdxm = ux(j  ,k  ,i-1,7);
                        dwdym = ux(j  ,k  ,i-1,8);
                        dwdzm = ux(j  ,k  ,i-1,9);
                    } else if (i == 1) {
                        fac = 1.;
                        dudxp = ux(j  ,k  ,i+1,1);
                        dudyp = ux(j  ,k  ,i+1,2);
                        dudzp = ux(j  ,k  ,i+1,3);
                        dvdxp = ux(j  ,k  ,i+1,4);
                        dvdyp = ux(j  ,k  ,i+1,5);
                        dvdzp = ux(j  ,k  ,i+1,6);
                        dwdxp = ux(j  ,k  ,i+1,7);
                        dwdyp = ux(j  ,k  ,i+1,8);
                        dwdzp = ux(j  ,k  ,i+1,9);
                        dudxm = ux(j  ,k  ,i  ,1);
                        dudym = ux(j  ,k  ,i  ,2);
                        dudzm = ux(j  ,k  ,i  ,3);
                        dvdxm = ux(j  ,k  ,i  ,4);
                        dvdym = ux(j  ,k  ,i  ,5);
                        dvdzm = ux(j  ,k  ,i  ,6);
                        dwdxm = ux(j  ,k  ,i  ,7);
                        dwdym = ux(j  ,k  ,i  ,8);
                        dwdzm = ux(j  ,k  ,i  ,9);
                    } else if (i == idim-1) {
                        fac = 1.;
                        dudxp = ux(j  ,k  ,i  ,1);
                        dudyp = ux(j  ,k  ,i  ,2);
                        dudzp = ux(j  ,k  ,i  ,3);
                        dvdxp = ux(j  ,k  ,i  ,4);
                        dvdyp = ux(j  ,k  ,i  ,5);
                        dvdzp = ux(j  ,k  ,i  ,6);
                        dwdxp = ux(j  ,k  ,i  ,7);
                        dwdyp = ux(j  ,k  ,i  ,8);
                        dwdzp = ux(j  ,k  ,i  ,9);
                        dudxm = ux(j  ,k  ,i-1,1);
                        dudym = ux(j  ,k  ,i-1,2);
                        dudzm = ux(j  ,k  ,i-1,3);
                        dvdxm = ux(j  ,k  ,i-1,4);
                        dvdym = ux(j  ,k  ,i-1,5);
                        dvdzm = ux(j  ,k  ,i-1,6);
                        dwdxm = ux(j  ,k  ,i-1,7);
                        dwdym = ux(j  ,k  ,i-1,8);
                        dwdzm = ux(j  ,k  ,i-1,9);
                    }
                    vx(j,k,i,2) = vx(j,k,i,2) + (xc*(dudxp-dudxm)
                                                 + yc*(dudyp-dudym)
                                                 + zc*(dudzp-dudzm))/fac;
                    vx(j,k,i,3) = vx(j,k,i,3) + (xc*(dvdxp-dvdxm)
                                                 + yc*(dvdyp-dvdym)
                                                 + zc*(dvdzp-dvdzm))/fac;
                    vx(j,k,i,4) = vx(j,k,i,4) + (xc*(dwdxp-dwdxm)
                                                 + yc*(dwdyp-dwdym)
                                                 + zc*(dwdzp-dwdzm))/fac;
                }
            }
        }
    }

    for (k = 1; k <= kdim-1; k++) {
        for (j = 1; j <= jdim-1; j++) {
            for (i = 1; i <= idim-1; i++) {
                vx(j,k,i,1) = std::sqrt(vx(j,k,i,2)*vx(j,k,i,2) +
                                         vx(j,k,i,3)*vx(j,k,i,3) +
                                         vx(j,k,i,4)*vx(j,k,i,4));
            }
        }
    }

    return;
}

} // namespace u_doubleprime_ns
