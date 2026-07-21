// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "sijrate2d.h"

namespace sijrate2d_ns {

void sijrate2d(int& idim, int& jdim, int& kdim,
               FortranArray4DRef<double> q,
               FortranArray4DRef<double> qj0,
               FortranArray4DRef<double> qk0,
               FortranArray3DRef<double> bcj,
               FortranArray3DRef<double> bck,
               FortranArray3DRef<double> vol,
               FortranArray4DRef<double> sj,
               FortranArray4DRef<double> sk,
               FortranArray4DRef<double> vx)
{
    // Local scalar temporaries
    int i, j, k, jm, jp, km, kp;
    double voll, volu;
    double skm1u, skm3u, skm4u, skm1x, skm3x, skm4x;
    double skp1u, skp3u, skp4u, skp1x, skp3x, skp4x;
    double sjm1u, sjm3u, sjm4u, sjm1x, sjm3x, sjm4x;
    double sjp1u, sjp3u, sjp4u, sjp1x, sjp3x, sjp4x;
    double xc, zc, tc, xp, zp, xm, zm, uu;
    double d2udn2, d2wdn2, d2udn3, d2wdn3;
    double dxunx, dxwnz, dxunz, dxwnx;

    // Put q (u and w) velocities in vx(1) and (2), and fill edges with
    // appropriate BC values:
    for (i = 1; i <= idim-1; i++) {
        for (j = 1; j <= jdim-1; j++) {
            for (k = 1; k <= kdim-1; k++) {
                vx(j,k,i,1) = q(j,k,i,2);
                vx(j,k,i,2) = q(j,k,i,4);
            }
        }
        // Get ghost cell values
        for (j = 1; j <= jdim-1; j++) {
            vx(j,0,i,1) = qk0(j,i,2,1)*(1.-bck(j,i,1)) +
                           (2.*qk0(j,i,2,1)-q(j,1,i,2))*bck(j,i,1);
            vx(j,0,i,2) = qk0(j,i,4,1)*(1.-bck(j,i,1)) +
                           (2.*qk0(j,i,4,1)-q(j,1,i,4))*bck(j,i,1);
            vx(j,kdim,i,1) = qk0(j,i,2,3)*(1.-bck(j,i,2)) +
                              (2.*qk0(j,i,2,3)-q(j,kdim-1,i,2))*bck(j,i,2);
            vx(j,kdim,i,2) = qk0(j,i,4,3)*(1.-bck(j,i,2)) +
                              (2.*qk0(j,i,4,3)-q(j,kdim-1,i,4))*bck(j,i,2);
        }
        for (k = 1; k <= kdim-1; k++) {
            vx(0,k,i,1) = qj0(k,i,2,1)*(1.-bcj(k,i,1)) +
                           (2.*qj0(k,i,2,1)-q(1,k,i,2))*bcj(k,i,1);
            vx(0,k,i,2) = qj0(k,i,4,1)*(1.-bcj(k,i,1)) +
                           (2.*qj0(k,i,4,1)-q(1,k,i,4))*bcj(k,i,1);
            vx(jdim,k,i,1) = qj0(k,i,2,3)*(1.-bcj(k,i,2)) +
                              (2.*qj0(k,i,2,3)-q(jdim-1,k,i,2))*bcj(k,i,2);
            vx(jdim,k,i,2) = qj0(k,i,4,3)*(1.-bcj(k,i,2)) +
                              (2.*qj0(k,i,4,3)-q(jdim-1,k,i,4))*bcj(k,i,2);
        }
    }

    // compute material derivative of principal strain directions (ignore time term)
    // the following assumes that the boundary terms are at ghost cells
    //   do for each component (2-D for now only):
    //   j-direction:
    for (i = 1; i <= idim-1; i++) {
        for (j = 1; j <= jdim-1; j++) {
            if (j == 1) {
                jm = j;
            } else {
                jm = j-1;
            }
            if (j == jdim-1) {
                jp = j;
            } else {
                jp = j+1;
            }
            for (k = 1; k <= kdim-1; k++) {
                voll = vol(jm,k,i);
                volu = vol(jp,k,i);
                skm1u = sk(jm,k+1,i,1);
                skm3u = sk(jm,k+1,i,3);
                skm4u = sk(jm,k+1,i,4);
                skm1x = sk(jm,k,i,1);
                skm3x = sk(jm,k,i,3);
                skm4x = sk(jm,k,i,4);
                skp1u = sk(jp,k+1,i,1);
                skp3u = sk(jp,k+1,i,3);
                skp4u = sk(jp,k+1,i,4);
                skp1x = sk(jp,k,i,1);
                skp3x = sk(jp,k,i,3);
                skp4x = sk(jp,k,i,4);
                xc = 0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4) +
                          sj(j,k,i  ,1)*sj(j,k,i  ,4))/vol(j,k,i);
                zc = 0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4) +
                          sj(j,k,i  ,3)*sj(j,k,i  ,4))/vol(j,k,i);
                tc = 0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4) +
                          sj(j,k,i  ,5)*sj(j,k,i  ,4))/vol(j,k,i);
                xp = sj(j+1,k,i,1)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volu));
                zp = sj(j+1,k,i,3)*sj(j+1,k,i,4)/(0.5*(vol(j,k,i)+volu));
                xm = sj(j  ,k,i,1)*sj(j  ,k,i,4)/(0.5*(vol(j,k,i)+voll));
                zm = sj(j  ,k,i,3)*sj(j  ,k,i,4)/(0.5*(vol(j,k,i)+voll));
                uu = xc*q(j,k,i,2) + zc*q(j,k,i,4) + tc;
                d2udn2 = xp*(vx(j+1,k,i,1)-vx(j  ,k,i,1)) -
                         xm*(vx(j  ,k,i,1)-vx(j-1,k,i,1));
                d2wdn2 = zp*(vx(j+1,k,i,2)-vx(j  ,k,i,2)) -
                         zm*(vx(j  ,k,i,2)-vx(j-1,k,i,2));
                d2udn3 = zp*(vx(j+1,k,i,1)-vx(j  ,k,i,1)) -
                         zm*(vx(j  ,k,i,1)-vx(j-1,k,i,1));
                d2wdn3 = xp*(vx(j+1,k,i,2)-vx(j  ,k,i,2)) -
                         xm*(vx(j  ,k,i,2)-vx(j-1,k,i,2));
                vx(j,k,i,3) = 0.5*uu*(d2udn2-d2wdn2);
                vx(j,k,i,4) = 0.5*uu*(d2udn3+d2wdn3);
                // x-deriv terms:
                xp = 0.25*(sk(j,k+1,i,1)*sk(j,k+1,i,4) +
                           skp1u*skp4u+sk(j,k,i,1)*sk(j,k,i,4) +
                           skp1x*skp4x)/(0.5*(vol(j,k,i)+volu));
                zp = 0.25*(sk(j,k+1,i,3)*sk(j,k+1,i,4) +
                           skp3u*skp4u+sk(j,k,i,3)*sk(j,k,i,4) +
                           skp3x*skp4x)/(0.5*(vol(j,k,i)+volu));
                xm = 0.25*(sk(j,k+1,i,1)*sk(j,k+1,i,4) +
                           skm1u*skm4u+sk(j,k,i,1)*sk(j,k,i,4) +
                           skm1x*skm4x)/(0.5*(vol(j,k,i)+voll));
                zm = 0.25*(sk(j,k+1,i,3)*sk(j,k+1,i,4) +
                           skm3u*skm4u+sk(j,k,i,3)*sk(j,k,i,4) +
                           skm3x*skm4x)/(0.5*(vol(j,k,i)+voll));
                dxunx = 0.25*(xp*(vx(jp ,k+1,i,1)-vx(jp ,k-1,i,1)-
                                  vx(j  ,k-1,i,1)+vx(j  ,k+1,i,1))
                             -xm*(vx(jm ,k+1,i,1)-vx(jm ,k-1,i,1)-
                                  vx(j  ,k-1,i,1)+vx(j  ,k+1,i,1)));
                dxwnz = 0.25*(zp*(vx(jp ,k+1,i,2)-vx(jp ,k-1,i,2)-
                                  vx(j  ,k-1,i,2)+vx(j  ,k+1,i,2))
                             -zm*(vx(jm ,k+1,i,2)-vx(jm ,k-1,i,2)-
                                  vx(j  ,k-1,i,2)+vx(j  ,k+1,i,2)));
                dxunz = 0.25*(zp*(vx(jp ,k+1,i,1)-vx(jp ,k-1,i,1)-
                                  vx(j  ,k-1,i,1)+vx(j  ,k+1,i,1))
                             -zm*(vx(jm ,k+1,i,1)-vx(jm ,k-1,i,1)-
                                  vx(j  ,k-1,i,1)+vx(j  ,k+1,i,1)));
                dxwnx = 0.25*(xp*(vx(jp ,k+1,i,2)-vx(jp ,k-1,i,2)-
                                  vx(j  ,k-1,i,2)+vx(j  ,k+1,i,2))
                             -xm*(vx(jm ,k+1,i,2)-vx(jm ,k-1,i,2)-
                                  vx(j  ,k-1,i,2)+vx(j  ,k+1,i,2)));
                vx(j,k,i,3) = vx(j,k,i,3) + 0.5*uu*(dxunx-dxwnz);
                vx(j,k,i,4) = vx(j,k,i,4) + 0.5*uu*(dxunz+dxwnx);
            }
        }
    }

    // k-direction:
    for (i = 1; i <= idim-1; i++) {
        for (k = 1; k <= kdim-1; k++) {
            if (k == 1) {
                km = k;
            } else {
                km = k-1;
            }
            if (k == kdim-1) {
                kp = k;
            } else {
                kp = k+1;
            }
            for (j = 1; j <= jdim-1; j++) {
                voll = vol(j,km,i);
                volu = vol(j,kp,i);
                sjm1u = sj(j+1,km,i,1);
                sjm3u = sj(j+1,km,i,3);
                sjm4u = sj(j+1,km,i,4);
                sjm1x = sj(j,km,i,1);
                sjm3x = sj(j,km,i,3);
                sjm4x = sj(j,km,i,4);
                sjp1u = sj(j+1,kp,i,1);
                sjp3u = sj(j+1,kp,i,3);
                sjp4u = sj(j+1,kp,i,4);
                sjp1x = sj(j,kp,i,1);
                sjp3x = sj(j,kp,i,3);
                sjp4x = sj(j,kp,i,4);
                xc = 0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4) +
                          sk(j,k,i  ,1)*sk(j,k,i  ,4))/vol(j,k,i);
                zc = 0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4) +
                          sk(j,k,i  ,3)*sk(j,k,i  ,4))/vol(j,k,i);
                tc = 0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4) +
                          sk(j,k,i  ,5)*sk(j,k,i  ,4))/vol(j,k,i);
                xp = sk(j,k+1,i,1)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volu));
                zp = sk(j,k+1,i,3)*sk(j,k+1,i,4)/(0.5*(vol(j,k,i)+volu));
                xm = sk(j,k  ,i,1)*sk(j,k  ,i,4)/(0.5*(vol(j,k,i)+voll));
                zm = sk(j,k  ,i,3)*sk(j,k  ,i,4)/(0.5*(vol(j,k,i)+voll));
                uu = xc*q(j,k,i,2) + zc*q(j,k,i,4) + tc;
                d2udn2 = xp*(vx(j,k+1,i,1)-vx(j,k  ,i,1)) -
                         xm*(vx(j,k  ,i,1)-vx(j,k-1,i,1));
                d2wdn2 = zp*(vx(j,k+1,i,2)-vx(j,k  ,i,2)) -
                         zm*(vx(j,k  ,i,2)-vx(j,k-1,i,2));
                d2udn3 = zp*(vx(j,k+1,i,1)-vx(j,k  ,i,1)) -
                         zm*(vx(j,k  ,i,1)-vx(j,k-1,i,1));
                d2wdn3 = xp*(vx(j,k+1,i,2)-vx(j,k  ,i,2)) -
                         xm*(vx(j,k  ,i,2)-vx(j,k-1,i,2));
                vx(j,k,i,3) = vx(j,k,i,3) + 0.5*uu*(d2udn2-d2wdn2);
                vx(j,k,i,4) = vx(j,k,i,4) + 0.5*uu*(d2udn3+d2wdn3);
                // x-deriv terms:
                xp = 0.25*(sj(j+1,k,i,1)*sj(j+1,k,i,4) +
                           sjp1u*sjp4u+sj(j,k,i,1)*sj(j,k,i,4) +
                           sjp1x*sjp4x)/(0.5*(vol(j,k,i)+volu));
                zp = 0.25*(sj(j+1,k,i,3)*sj(j+1,k,i,4) +
                           sjp3u*sjp4u+sj(j,k,i,3)*sj(j,k,i,4) +
                           sjp3x*sjp4x)/(0.5*(vol(j,k,i)+volu));
                xm = 0.25*(sj(j+1,k,i,1)*sj(j+1,k,i,4) +
                           sjm1u*sjm4u+sj(j,k,i,1)*sj(j,k,i,4) +
                           sjm1x*sjm4x)/(0.5*(vol(j,k,i)+voll));
                zm = 0.25*(sj(j+1,k,i,3)*sj(j+1,k,i,4) +
                           sjm3u*sjm4u+sj(j,k,i,3)*sj(j,k,i,4) +
                           sjm3x*sjm4x)/(0.5*(vol(j,k,i)+voll));
                dxunx = 0.25*(xp*(vx(j+1,kp ,i,1)-vx(j-1,kp ,i,1)-
                                  vx(j-1,k  ,i,1)+vx(j+1,k  ,i,1))
                             -xm*(vx(j+1,km ,i,1)-vx(j-1,km ,i,1)-
                                  vx(j-1,k  ,i,1)+vx(j+1,k  ,i,1)));
                dxwnz = 0.25*(zp*(vx(j+1,kp ,i,2)-vx(j-1,kp ,i,2)-
                                  vx(j-1,k  ,i,2)+vx(j+1,k  ,i,2))
                             -zm*(vx(j+1,km ,i,2)-vx(j-1,km ,i,2)-
                                  vx(j-1,k  ,i,2)+vx(j+1,k  ,i,2)));
                dxunz = 0.25*(zp*(vx(j+1,kp ,i,1)-vx(j-1,kp ,i,1)-
                                  vx(j-1,k  ,i,1)+vx(j+1,k  ,i,1))
                             -zm*(vx(j+1,km ,i,1)-vx(j-1,km ,i,1)-
                                  vx(j-1,k  ,i,1)+vx(j+1,k  ,i,1)));
                dxwnx = 0.25*(xp*(vx(j+1,kp ,i,2)-vx(j-1,kp ,i,2)-
                                  vx(j-1,k  ,i,2)+vx(j+1,k  ,i,2))
                             -xm*(vx(j+1,km ,i,2)-vx(j-1,km ,i,2)-
                                  vx(j-1,k  ,i,2)+vx(j+1,k  ,i,2)));
                vx(j,k,i,3) = vx(j,k,i,3) + 0.5*uu*(dxunx-dxwnz);
                vx(j,k,i,4) = vx(j,k,i,4) + 0.5*uu*(dxunz+dxwnx);
            }
        }
    }

    return;
}

} // namespace sijrate2d_ns
