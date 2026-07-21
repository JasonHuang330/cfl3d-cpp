// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bcnonin.h"
#include "bc.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <algorithm>

namespace bcnonin_ns {

// bc() is a thin dispatcher wrapper calling bc_ns::bc
void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork,
        double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0,
        FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
        FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat,
        int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil, int& nummem)
{
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf, ibufdim,
              maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg, nblg,
              nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
              ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek,
              lwdat, myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

void bcnonin(int& jdim, int& kdim, int& idim,
             FortranArray4DRef<double> q, FortranArray4DRef<double> qj0,
             FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
             FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
             FortranArray4DRef<double> si, FortranArray3DRef<double> bcj,
             FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
             FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk,
             FortranArray4DRef<double> xtbi, FortranArray4DRef<double> atbj,
             FortranArray4DRef<double> atbk, FortranArray4DRef<double> atbi,
             int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
             int& nface, int& iuns, FortranArray1DRef<int> nou,
             FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
             FortranArray3DRef<double> x, FortranArray3DRef<double> y,
             FortranArray3DRef<double> z, int& nbl)
{
    // COMMON block aliases
    int32_t& isklton   = cmn_sklton.isklton;
    float& xcentrot    = cmn_noninertial.xcentrot;
    float& ycentrot    = cmn_noninertial.ycentrot;
    float& zcentrot    = cmn_noninertial.zcentrot;
    float& xrotrate    = cmn_noninertial.xrotrate;
    float& yrotrate    = cmn_noninertial.yrotrate;
    float& zrotrate    = cmn_noninertial.zrotrate;

    // Local variables
    double dx[2], dy[2], dz[2];
    double wx, wy, wz;
    double cx, cy, cz, gx, gy, gz;
    int i, j, k, l;

    wx = (double)xrotrate;
    wy = (double)yrotrate;
    wz = (double)zrotrate;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;


    //**************************************************************************
    //     j=1 boundary NON-INERTIAL freestream                    bctype 1000
    //**************************************************************************
    if (nface == 3) {
        if (isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   j=   1  NONINERTIAL freestream         type 1000"
                "  i=%3d,%3d  k=%3d,%3d",
                ista, iend, ksta, kend);
        }
        for (i = ista; i <= iend1; i++) {
        for (k = ksta; k <= kend1; k++) {

            j = 1;

            cx = 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            cy = 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            cz = 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            dx[0] = 0.25 * (
                x(j  , k  , i  ) - x(j+1, k  , i  ) +
                x(j  , k  , i+1) - x(j+1, k  , i+1) +
                x(j  , k+1, i  ) - x(j+1, k+1, i  ) +
                x(j  , k+1, i+1) - x(j+1, k+1, i+1) );

            dy[0] = 0.25 * (
                y(j  , k  , i  ) - y(j+1, k  , i  ) +
                y(j  , k  , i+1) - y(j+1, k  , i+1) +
                y(j  , k+1, i  ) - y(j+1, k+1, i  ) +
                y(j  , k+1, i+1) - y(j+1, k+1, i+1) );

            dz[0] = 0.25 * (
                z(j  , k  , i  ) - z(j+1, k  , i  ) +
                z(j  , k  , i+1) - z(j+1, k  , i+1) +
                z(j  , k+1, i  ) - z(j+1, k+1, i  ) +
                z(j  , k+1, i+1) - z(j+1, k+1, i+1) );

            j = 2;

            dx[1] = dx[0] + cx - 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            dy[1] = dy[0] + cy - 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            dz[1] = dz[0] + cz - 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            cx = cx - (double)xcentrot;
            cy = cy - (double)ycentrot;
            cz = cz - (double)zcentrot;

            for (l = 1; l <= 2; l++) {
                gx = cx + dx[l-1];
                gy = cy + dy[l-1];
                gz = cz + dz[l-1];

                // Uinf + r x omega is the same as Uinf - omega x r
                qj0(k, i, 2, l) = qj0(k, i, 2, l) + ( gy * wz - gz * wy );
                qj0(k, i, 3, l) = qj0(k, i, 3, l) + ( gz * wx - gx * wz );
                qj0(k, i, 4, l) = qj0(k, i, 4, l) + ( gx * wy - gy * wx );

                bcj(k, i, 1) = 0.0;
            }
        } // k
        } // i
    } // nface==3


    //**************************************************************************
    //     j=jdim boundary NON-INERTIAL freestream                 bctype 1000
    //**************************************************************************
    if (nface == 4) {
        if (isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   j=jdim  NONINERTIAL freestream         type 1000"
                "  i=%3d,%3d  k=%3d,%3d",
                ista, iend, ksta, kend);
        }
        for (i = ista; i <= iend1; i++) {
        for (k = ksta; k <= kend1; k++) {

            j = jdim1;

            cx = 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            cy = 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            cz = 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            dx[0] = 0.25 * (
                x(j+1, k  , i  ) - x(j  , k  , i  ) +
                x(j+1, k  , i+1) - x(j  , k  , i+1) +
                x(j+1, k+1, i  ) - x(j  , k+1, i  ) +
                x(j+1, k+1, i+1) - x(j  , k+1, i+1) );

            dy[0] = 0.25 * (
                y(j+1, k  , i  ) - y(j  , k  , i  ) +
                y(j+1, k  , i+1) - y(j  , k  , i+1) +
                y(j+1, k+1, i  ) - y(j  , k+1, i  ) +
                y(j+1, k+1, i+1) - y(j  , k+1, i+1) );

            dz[0] = 0.25 * (
                z(j+1, k  , i  ) - z(j  , k  , i  ) +
                z(j+1, k  , i+1) - z(j  , k  , i+1) +
                z(j+1, k+1, i  ) - z(j  , k+1, i  ) +
                z(j+1, k+1, i+1) - z(j  , k+1, i+1) );

            j = jdim1 - 1;

            dx[1] = dx[0] + cx - 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            dy[1] = dy[0] + cy - 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            dz[1] = dz[0] + cz - 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            cx = cx - (double)xcentrot;
            cy = cy - (double)ycentrot;
            cz = cz - (double)zcentrot;

            // l=3,4 in Fortran; dx[l-3] maps to dx[0] and dx[1]
            for (l = 3; l <= 4; l++) {
                gx = cx + dx[l-3];
                gy = cy + dy[l-3];
                gz = cz + dz[l-3];

                // Uinf + r x omega is the same as Uinf - omega x r
                qj0(k, i, 2, l) = qj0(k, i, 2, l) + ( gy * wz - gz * wy );
                qj0(k, i, 3, l) = qj0(k, i, 3, l) + ( gz * wx - gx * wz );
                qj0(k, i, 4, l) = qj0(k, i, 4, l) + ( gx * wy - gy * wx );

                bcj(k, i, 2) = 0.0;
            }
        } // k
        } // i
    } // nface==4


    //**************************************************************************
    //     k=1 boundary NON-INERTIAL freestream                    bctype 1000
    //**************************************************************************
    if (nface == 5) {
        if (isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   k=   1  NONINERTIAL freestream         type 1000"
                "  i=%3d,%3d  j=%3d,%3d",
                ista, iend, ksta, kend);
        }
        for (i = ista; i <= iend1; i++) {
        for (j = jsta; j <= jend1; j++) {

            k = 1;

            cx = 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            cy = 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            cz = 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            dx[0] = 0.25 * (
                x(j  , k  , i  ) - x(j  , k+1, i  ) +
                x(j  , k  , i+1) - x(j  , k+1, i+1) +
                x(j+1, k  , i  ) - x(j+1, k+1, i  ) +
                x(j+1, k  , i+1) - x(j+1, k+1, i+1) );

            dy[0] = 0.25 * (
                y(j  , k  , i  ) - y(j  , k+1, i  ) +
                y(j  , k  , i+1) - y(j  , k+1, i+1) +
                y(j+1, k  , i  ) - y(j+1, k+1, i  ) +
                y(j+1, k  , i+1) - y(j+1, k+1, i+1) );

            dz[0] = 0.25 * (
                z(j  , k  , i  ) - z(j  , k+1, i  ) +
                z(j  , k  , i+1) - z(j  , k+1, i+1) +
                z(j+1, k  , i  ) - z(j+1, k+1, i  ) +
                z(j+1, k  , i+1) - z(j+1, k+1, i+1) );

            k = 2;

            dx[1] = dx[0] + cx - 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            dy[1] = dy[0] + cy - 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            dz[1] = dz[0] + cz - 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            cx = cx - (double)xcentrot;
            cy = cy - (double)ycentrot;
            cz = cz - (double)zcentrot;

            for (l = 1; l <= 2; l++) {
                gx = cx + dx[l-1];
                gy = cy + dy[l-1];
                gz = cz + dz[l-1];

                // Uinf + r x omega is the same as Uinf - omega x r
                qk0(j, i, 2, l) = qk0(j, i, 2, l) + ( gy * wz - gz * wy );
                qk0(j, i, 3, l) = qk0(j, i, 3, l) + ( gz * wx - gx * wz );
                qk0(j, i, 4, l) = qk0(j, i, 4, l) + ( gx * wy - gy * wx );
                bck(j, i, 1) = 0.0;
            }
        } // j
        } // i
    } // nface==5


    //**************************************************************************
    //     k=kdim boundary NON-INERTIAL freestream                 bctype 1000
    //**************************************************************************
    if (nface == 6) {
        if (isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   k=kdim  NONINERTIAL freestream         type 1000"
                "  i=%3d,%3d  j=%3d,%3d",
                ista, iend, ksta, kend);
        }
        for (i = ista; i <= iend1; i++) {
        for (j = jsta; j <= jend1; j++) {

            k = kdim1;

            cx = 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            cy = 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            cz = 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            dx[0] = 0.25 * (
                x(j  , k+1, i  ) - x(j  , k  , i  ) +
                x(j  , k+1, i+1) - x(j  , k  , i+1) +
                x(j+1, k+1, i  ) - x(j+1, k  , i  ) +
                x(j+1, k+1, i+1) - x(j+1, k  , i+1) );

            dy[0] = 0.25 * (
                y(j  , k+1, i  ) - y(j  , k  , i  ) +
                y(j  , k+1, i+1) - y(j  , k  , i+1) +
                y(j+1, k+1, i  ) - y(j+1, k  , i  ) +
                y(j+1, k+1, i+1) - y(j+1, k  , i+1) );

            dz[0] = 0.25 * (
                z(j  , k+1, i  ) - z(j  , k  , i  ) +
                z(j  , k+1, i+1) - z(j  , k  , i+1) +
                z(j+1, k+1, i  ) - z(j+1, k  , i  ) +
                z(j+1, k+1, i+1) - z(j+1, k  , i+1) );

            k = kdim1 - 1;

            dx[1] = dx[0] + cx - 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            dy[1] = dy[0] + cy - 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            dz[1] = dz[0] + cz - 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            cx = cx - (double)xcentrot;
            cy = cy - (double)ycentrot;
            cz = cz - (double)zcentrot;

            // l=3,4 in Fortran; dx[l-3] maps to dx[0] and dx[1]
            for (l = 3; l <= 4; l++) {
                gx = cx + dx[l-3];
                gy = cy + dy[l-3];
                gz = cz + dz[l-3];

                // Uinf + r x omega is the same as Uinf - omega x r
                qk0(j, i, 2, l) = qk0(j, i, 2, l) + ( gy * wz - gz * wy );
                qk0(j, i, 3, l) = qk0(j, i, 3, l) + ( gz * wx - gx * wz );
                qk0(j, i, 4, l) = qk0(j, i, 4, l) + ( gx * wy - gy * wx );

                bck(j, i, 2) = 0.0;
            }
        } // j
        } // i
    } // nface==6


    //**************************************************************************
    //     i=1 boundary NON-INERTIAL freestream                    bctype 1000
    //**************************************************************************
    if (nface == 1) {
        if (isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   i=   1  NONINERTIAL freestream         type 1000"
                "  j=%3d,%3d  k=%3d,%3d",
                ista, iend, ksta, kend);
        }
        for (k = ksta; k <= kend1; k++) {
        for (j = jsta; j <= jend1; j++) {

            i = 1;

            cx = 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            cy = 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            cz = 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            dx[0] = 0.25 * (
                x(j  , k  , i  ) - x(j  , k  , i+1) +
                x(j  , k+1, i  ) - x(j  , k+1, i+1) +
                x(j+1, k  , i  ) - x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) - x(j+1, k+1, i+1) );

            dy[0] = 0.25 * (
                y(j  , k  , i  ) - y(j  , k  , i+1) +
                y(j  , k+1, i  ) - y(j  , k+1, i+1) +
                y(j+1, k  , i  ) - y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) - y(j+1, k+1, i+1) );

            dz[0] = 0.25 * (
                z(j  , k  , i  ) - z(j  , k  , i+1) +
                z(j  , k+1, i  ) - z(j  , k+1, i+1) +
                z(j+1, k  , i  ) - z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) - z(j+1, k+1, i+1) );

            i = 2;

            dx[1] = dx[0] + cx - 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            dy[1] = dy[0] + cy - 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            dz[1] = dz[0] + cz - 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            cx = cx - (double)xcentrot;
            cy = cy - (double)ycentrot;
            cz = cz - (double)zcentrot;

            for (l = 1; l <= 2; l++) {
                gx = cx + dx[l-1];
                gy = cy + dy[l-1];
                gz = cz + dz[l-1];

                // Uinf + r x omega is the same as Uinf - omega x r
                qi0(j, k, 2, l) = qi0(j, k, 2, l) + ( gy * wz - gz * wy );
                qi0(j, k, 3, l) = qi0(j, k, 3, l) + ( gz * wx - gx * wz );
                qi0(j, k, 4, l) = qi0(j, k, 4, l) + ( gx * wy - gy * wx );
                bci(j, i, 1) = 0.0;
            }
        } // j
        } // k
    } // nface==1


    //**************************************************************************
    //     i=idim boundary NON-INERTIAL freestream                 bctype 1000
    //**************************************************************************
    if (nface == 2) {
        if (isklton == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   i=idim  NONINERTIAL freestream         type 1000"
                "  j=%3d,%3d  k=%3d,%3d",
                ista, iend, ksta, kend);
        }
        for (k = ksta; k <= kend1; k++) {
        for (j = jsta; j <= jend1; j++) {

            i = idim1;

            cx = 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            cy = 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            cz = 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            dx[0] = 0.25 * (
                x(j  , k  , i+1) - x(j  , k  , i  ) +
                x(j  , k+1, i+1) - x(j  , k+1, i  ) +
                x(j+1, k  , i+1) - x(j+1, k  , i  ) +
                x(j+1, k+1, i+1) - x(j+1, k+1, i  ) );

            dy[0] = 0.25 * (
                y(j  , k  , i+1) - y(j  , k  , i  ) +
                y(j  , k+1, i+1) - y(j  , k+1, i  ) +
                y(j+1, k  , i+1) - y(j+1, k  , i  ) +
                y(j+1, k+1, i+1) - y(j+1, k+1, i  ) );

            dz[0] = 0.25 * (
                z(j  , k  , i+1) - z(j  , k  , i  ) +
                z(j  , k+1, i+1) - z(j  , k+1, i  ) +
                z(j+1, k  , i+1) - z(j+1, k  , i  ) +
                z(j+1, k+1, i+1) - z(j+1, k+1, i  ) );

            i = idim1 - 1;

            dx[1] = dx[0] + cx - 0.125 * (
                x(j  , k  , i  ) + x(j  , k  , i+1) +
                x(j  , k+1, i  ) + x(j  , k+1, i+1) +
                x(j+1, k  , i  ) + x(j+1, k  , i+1) +
                x(j+1, k+1, i  ) + x(j+1, k+1, i+1) );

            dy[1] = dy[0] + cy - 0.125 * (
                y(j  , k  , i  ) + y(j  , k  , i+1) +
                y(j  , k+1, i  ) + y(j  , k+1, i+1) +
                y(j+1, k  , i  ) + y(j+1, k  , i+1) +
                y(j+1, k+1, i  ) + y(j+1, k+1, i+1) );

            dz[1] = dz[0] + cz - 0.125 * (
                z(j  , k  , i  ) + z(j  , k  , i+1) +
                z(j  , k+1, i  ) + z(j  , k+1, i+1) +
                z(j+1, k  , i  ) + z(j+1, k  , i+1) +
                z(j+1, k+1, i  ) + z(j+1, k+1, i+1) );

            cx = cx - (double)xcentrot;
            cy = cy - (double)ycentrot;
            cz = cz - (double)zcentrot;

            // l=3,4 in Fortran; dx[l-3] maps to dx[0] and dx[1]
            for (l = 3; l <= 4; l++) {
                gx = cx + dx[l-3];
                gy = cy + dy[l-3];
                gz = cz + dz[l-3];

                // Uinf + r x omega is the same as Uinf - omega x r
                qi0(j, k, 2, l) = qi0(j, k, 2, l) + ( gy * wz - gz * wy );
                qi0(j, k, 3, l) = qi0(j, k, 3, l) + ( gz * wx - gx * wz );
                qi0(j, k, 4, l) = qi0(j, k, 4, l) + ( gx * wy - gy * wx );
                bci(j, i, 2) = 0.0;
            }
        } // j
        } // k
    } // nface==2

} // bcnonin

} // namespace bcnonin_ns
