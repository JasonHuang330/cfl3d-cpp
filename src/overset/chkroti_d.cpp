// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "chkroti_d.h"
#include "grdmove.h"
#include "termn8.h"
#include "ccomplex.h"
#include "chkrotj_d.h"
#include "chkrotk_d.h"
#include <cstdio>
#include <cmath>

namespace chkroti_d_ns {

void chkroti_d(int& nbl, int& jdim, int& kdim, int& idim,
               FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
               int& nblp, int& jdimp, int& kdimp, int& idimp,
               int& nface, FortranArray4DRef<double> bcdata,
               FortranArray3DRef<double> xyzip,
               int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
               int& mdim, int& ndim, int& lcnt,
               FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
               int& maxbl, FortranArray1DRef<double> period_miss, int& lbcprd,
               FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
               int& nbuf, int& ibufdim, int& myid)
{
    // COMMON /conversion/ radtodeg
    float& radtodeg = cmn_conversion.radtodeg;

    int ip = 1;

    double dthtx = bcdata(1,1,ip,2) / (double)radtodeg;
    double dthty = bcdata(1,1,ip,3) / (double)radtodeg;
    double dthtz = bcdata(1,1,ip,4) / (double)radtodeg;

    // check that only one of dthtx,dthty,dthtz is non-zero
    if (dthtx != 0) {
        if (dthty != 0. || dthtz != 0.) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " 2 of the 3 dtht values currently must = 0");
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }
    if (dthty != 0) {
        if (dthtx != 0. || dthtz != 0.) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " 2 of the 3 dtht values currently must = 0");
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }
    if (dthtz != 0) {
        if (dthtx != 0. || dthty != 0.) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " 2 of the 3 dtht values currently must = 0");
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }
    }

    double eps = 0;

    // i=1 boundary
    if (nface == 1) {

        int i  = 1;
        ip = idimp;

        // check that jdim = jdimp and kdim = kdimp
        if (jdim != jdimp || kdim != kdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " periodic block face must be of the same 2 dimensions (and orientation) as the current block face");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " jdim,kdim,jdimp,kdimp = %5d%5d%5d%5d", jdim, kdim, jdimp, kdimp);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // rotate periodic block face (entire face)
        int jdum = jdimp;
        int kdum = kdimp;
        int idum = 1;
        // Build 3D array views into xyzip slices for grdmove
        FortranArray3DRef<double> xip(&xyzip(1,1,1), jdum, kdum, idum);
        FortranArray3DRef<double> yip(&xyzip(1,1,2), jdum, kdum, idum);
        FortranArray3DRef<double> zip(&xyzip(1,1,3), jdum, kdum, idum);
        grdmove_ns::grdmove(nblp, jdum, kdum, idum,
                            xip, yip, zip,
                            xorig(nblp), yorig(nblp), zorig(nblp),
                            xorig(nblp), yorig(nblp), zorig(nblp),
                            dthtx, dthty, dthtz);

        // check for point match
        for (int j = jsta; j <= jend; j++) {
            for (int k = ksta; k <= kend; k++) {
                double xt1 = x(j,k,i);
                double yt1 = y(j,k,i);
                double zt1 = z(j,k,i);
                double xr1 = xyzip(j,k,1);
                double yr1 = xyzip(j,k,2);
                double zr1 = xyzip(j,k,3);
                double dist = std::sqrt((xr1-xt1)*(xr1-xt1) + (yr1-yt1)*(yr1-yt1) + (zr1-zt1)*(zr1-zt1));
                eps = ccomplex_ns::ccmax(eps, dist);
            }
        }

        period_miss(lcnt) = eps;

    }

    // i=idim boundary
    if (nface == 2) {

        int i  = idim;
        ip = 1;

        // check that jdim = jdimp and kdim = kdimp
        if (jdim != jdimp || kdim != kdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " periodic block face must be of the same 2 dimensions (and orientation) as the current block face");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " jdim,kdim,jdimp,kdimp = %5d%5d%5d%5d", jdim, kdim, jdimp, kdimp);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // rotate periodic block face (entire face)
        int jdum = jdimp;
        int kdum = kdimp;
        int idum = 1;
        // Build 3D array views into xyzip slices for grdmove
        FortranArray3DRef<double> xip(&xyzip(1,1,1), jdum, kdum, idum);
        FortranArray3DRef<double> yip(&xyzip(1,1,2), jdum, kdum, idum);
        FortranArray3DRef<double> zip(&xyzip(1,1,3), jdum, kdum, idum);
        grdmove_ns::grdmove(nblp, jdum, kdum, idum,
                            xip, yip, zip,
                            xorig(nblp), yorig(nblp), zorig(nblp),
                            xorig(nblp), yorig(nblp), zorig(nblp),
                            dthtx, dthty, dthtz);

        // check for point match
        for (int j = jsta; j <= jend; j++) {
            for (int k = ksta; k <= kend; k++) {
                double xt1 = x(j,k,i);
                double yt1 = y(j,k,i);
                double zt1 = z(j,k,i);
                double xr1 = xyzip(j,k,1);
                double yr1 = xyzip(j,k,2);
                double zr1 = xyzip(j,k,3);
                double dist = std::sqrt((xr1-xt1)*(xr1-xt1) + (yr1-yt1)*(yr1-yt1) + (zr1-zt1)*(zr1-zt1));
                eps = ccomplex_ns::ccmax(eps, dist);
            }
        }

        period_miss(lcnt) = eps;

    }

    // formats:
    // 101: " 2 of the 3 dtht values currently must = 0"
    // 102: " periodic block face must be of the same 2 dimensions (and orientation) as the current block face"
    // 103: " jdim,kdim,jdimp,kdimp = %5d%5d%5d%5d"
}

void chkrot(int& nbl, int& jdim, int& kdim, int& idim,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
            int& nblp, int& jdimp, int& kdimp, int& idimp,
            FortranArray3DRef<double> xp, FortranArray3DRef<double> yp, FortranArray3DRef<double> zp,
            int& nface, FortranArray4DRef<double> bcdata,
            FortranArray3DRef<double> xyzjp, FortranArray3DRef<double> xyzkp, FortranArray3DRef<double> xyzip,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& mdim, int& ndim, int& lcnt,
            FortranArray1DRef<double> xorig, FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
            int& maxbl, FortranArray1DRef<double> period_miss, int& lbcprd,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid)
{
    chkrotj_d_ns::chkrotj_d(nbl, jdim, kdim, idim, x, y, z,
                              nblp, jdimp, kdimp, idimp,
                              nface, bcdata, xyzjp,
                              ista, iend, jsta, jend, ksta, kend,
                              mdim, ndim, lcnt,
                              xorig, yorig, zorig,
                              maxbl, period_miss, lbcprd,
                              nou, bou, nbuf, ibufdim, myid);

    chkrotk_d_ns::chkrotk_d(nbl, jdim, kdim, idim, x, y, z,
                              nblp, jdimp, kdimp, idimp,
                              nface, bcdata, xyzkp,
                              ista, iend, jsta, jend, ksta, kend,
                              mdim, ndim, lcnt,
                              xorig, yorig, zorig,
                              maxbl, period_miss, lbcprd,
                              nou, bou, nbuf, ibufdim, myid);

    chkroti_d(nbl, jdim, kdim, idim, x, y, z,
              nblp, jdimp, kdimp, idimp,
              nface, bcdata, xyzip,
              ista, iend, jsta, jend, ksta, kend,
              mdim, ndim, lcnt,
              xorig, yorig, zorig,
              maxbl, period_miss, lbcprd,
              nou, bou, nbuf, ibufdim, myid);
}

} // namespace chkroti_d_ns
