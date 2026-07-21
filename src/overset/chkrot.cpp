// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "chkrot.h"
#include "grdmove.h"
#include "ccomplex.h"
#include "termn8.h"
#include <cstdio>
#include <cmath>

namespace chkrot_ns {

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
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid)
{
    // local variables
    int ip, i, j, k, jp, kp;
    int jdum, kdum, idum;
    double dthtx, dthty, dthtz;
    double eps;
    double xt1, yt1, zt1, xr1, yr1, zr1;

    ip = 1;

    // The caller (bc_period) hands xyzjp/xyzkp/xyzip in as flat work-array
    // slices; chkrot's Fortran contract re-dimensions them per its own
    // `dimension xyzjp(kdim,idim,3),xyzkp(jdim,idim,3),xyzip(jdim,kdim,3)`.
    // Re-view them here so xyz*(j,k,mm) indexes with the correct strides
    // (otherwise the periodic-face coordinates are read with the wrong shape,
    // yielding a spurious "periodic mismatch" and zeroed ghost cells).
    xyzjp = FortranArray3DRef<double>(xyzjp.data(), kdim, idim, 3);
    xyzkp = FortranArray3DRef<double>(xyzkp.data(), jdim, idim, 3);
    xyzip = FortranArray3DRef<double>(xyzip.data(), jdim, kdim, 3);

    dthtx = bcdata(1,1,ip,2) / (double)cmn_conversion.radtodeg;
    dthty = bcdata(1,1,ip,3) / (double)cmn_conversion.radtodeg;
    dthtz = bcdata(1,1,ip,4) / (double)cmn_conversion.radtodeg;

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

    eps = 0;

    // i=1 boundary
    if (nface == 1) {

        i  = 1;
        ip = idimp;

        // check that jdim = jdimp and kdim = kdimp
        if (jdim != jdimp || kdim != kdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " periodic block face must be of the same 2 dimensions (and orientation) as the current block face");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " jdim,kdim,jdimp,kdimp = %5d%5d%5d%5d", jdim, kdim, jdimp, kdimp);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // load x,y,z values on periodic face into work array
        for (int jj = 1; jj <= jdimp; jj++) {
            for (int kk = 1; kk <= kdimp; kk++) {
                xyzip(jj,kk,1) = xp(jj,kk,ip);
                xyzip(jj,kk,2) = yp(jj,kk,ip);
                xyzip(jj,kk,3) = zp(jj,kk,ip);
            }
        }

        // rotate periodic block face (entire face)
        jdum = jdimp;
        kdum = kdimp;
        idum = 1;
        {
            FortranArray3DRef<double> xslice(&xyzip(1,1,1), jdum, kdum, idum);
            FortranArray3DRef<double> yslice(&xyzip(1,1,2), jdum, kdum, idum);
            FortranArray3DRef<double> zslice(&xyzip(1,1,3), jdum, kdum, idum);
            grdmove_ns::grdmove(nblp, jdum, kdum, idum,
                                xslice, yslice, zslice,
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                dthtx, dthty, dthtz);
        }

        // check for point match
        for (int jj = jsta; jj <= jend; jj++) {
            for (int kk = ksta; kk <= kend; kk++) {
                xt1 = x(jj,kk,i);
                yt1 = y(jj,kk,i);
                zt1 = z(jj,kk,i);
                xr1 = xyzip(jj,kk,1);
                yr1 = xyzip(jj,kk,2);
                zr1 = xyzip(jj,kk,3);
                double val = std::sqrt((xr1-xt1)*(xr1-xt1)+(yr1-yt1)*(yr1-yt1)+(zr1-zt1)*(zr1-zt1));
                eps = ccomplex_ns::ccmax(eps, val);
            }
        }

        period_miss(lcnt) = eps;

    }

    // i=idim boundary
    if (nface == 2) {

        i  = idim;
        ip = 1;

        // check that jdim = jdimp and kdim = kdimp
        if (jdim != jdimp || kdim != kdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " periodic block face must be of the same 2 dimensions (and orientation) as the current block face");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " jdim,kdim,jdimp,kdimp = %5d%5d%5d%5d", jdim, kdim, jdimp, kdimp);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // load x,y,z values on periodic face into work array
        for (int jj = 1; jj <= jdimp; jj++) {
            for (int kk = 1; kk <= kdimp; kk++) {
                xyzip(jj,kk,1) = xp(jj,kk,ip);
                xyzip(jj,kk,2) = yp(jj,kk,ip);
                xyzip(jj,kk,3) = zp(jj,kk,ip);
            }
        }

        // rotate periodic block face (entire face)
        jdum = jdimp;
        kdum = kdimp;
        idum = 1;
        {
            FortranArray3DRef<double> xslice(&xyzip(1,1,1), jdum, kdum, idum);
            FortranArray3DRef<double> yslice(&xyzip(1,1,2), jdum, kdum, idum);
            FortranArray3DRef<double> zslice(&xyzip(1,1,3), jdum, kdum, idum);
            grdmove_ns::grdmove(nblp, jdum, kdum, idum,
                                xslice, yslice, zslice,
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                dthtx, dthty, dthtz);
        }

        // check for point match
        for (int jj = jsta; jj <= jend; jj++) {
            for (int kk = ksta; kk <= kend; kk++) {
                xt1 = x(jj,kk,i);
                yt1 = y(jj,kk,i);
                zt1 = z(jj,kk,i);
                xr1 = xyzip(jj,kk,1);
                yr1 = xyzip(jj,kk,2);
                zr1 = xyzip(jj,kk,3);
                double val = std::sqrt((xr1-xt1)*(xr1-xt1)+(yr1-yt1)*(yr1-yt1)+(zr1-zt1)*(zr1-zt1));
                eps = ccomplex_ns::ccmax(eps, val);
            }
        }

        period_miss(lcnt) = eps;

    }


    // j=1 boundary
    if (nface == 3) {

        j  = 1;
        jp = jdimp;

        // check that idim = idimp and kdim = kdimp
        if (idim != idimp || kdim != kdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " periodic block face must be of the same 2 dimensions (and orientation) as the current block face");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " idim,kdim,idimp,kdimp = %5d%5d%5d%5d", idim, kdim, idimp, kdimp);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // load x,y,z values on periodic face into work array
        for (int ii = 1; ii <= idimp; ii++) {
            for (int kk = 1; kk <= kdimp; kk++) {
                xyzjp(kk,ii,1) = xp(jp,kk,ii);
                xyzjp(kk,ii,2) = yp(jp,kk,ii);
                xyzjp(kk,ii,3) = zp(jp,kk,ii);
            }
        }

        // rotate periodic block face (entire face)
        jdum = 1;
        kdum = kdimp;
        idum = idimp;
        {
            FortranArray3DRef<double> xslice(&xyzjp(1,1,1), jdum, kdum, idum);
            FortranArray3DRef<double> yslice(&xyzjp(1,1,2), jdum, kdum, idum);
            FortranArray3DRef<double> zslice(&xyzjp(1,1,3), jdum, kdum, idum);
            grdmove_ns::grdmove(nblp, jdum, kdum, idum,
                                xslice, yslice, zslice,
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                dthtx, dthty, dthtz);
        }

        // check for point match
        for (int ii = ista; ii <= iend; ii++) {
            for (int kk = ksta; kk <= kend; kk++) {
                xt1 = x(j,kk,ii);
                yt1 = y(j,kk,ii);
                zt1 = z(j,kk,ii);
                xr1 = xyzjp(kk,ii,1);
                yr1 = xyzjp(kk,ii,2);
                zr1 = xyzjp(kk,ii,3);
                double val = std::sqrt((xr1-xt1)*(xr1-xt1)+(yr1-yt1)*(yr1-yt1)+(zr1-zt1)*(zr1-zt1));
                eps = ccomplex_ns::ccmax(eps, val);
            }
        }

        period_miss(lcnt) = eps;

    }

    // j=jdim boundary
    if (nface == 4) {

        j  = jdim;
        jp = 1;

        // check that idim = idimp and kdim = kdimp
        if (idim != idimp || kdim != kdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " periodic block face must be of the same 2 dimensions (and orientation) as the current block face");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " idim,kdim,idimp,kdimp = %5d%5d%5d%5d", idim, kdim, idimp, kdimp);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // load x,y,z values on periodic face into work array
        for (int ii = 1; ii <= idimp; ii++) {
            for (int kk = 1; kk <= kdimp; kk++) {
                xyzjp(kk,ii,1) = xp(jp,kk,ii);
                xyzjp(kk,ii,2) = yp(jp,kk,ii);
                xyzjp(kk,ii,3) = zp(jp,kk,ii);
            }
        }

        // rotate periodic block face (entire face)
        jdum = 1;
        kdum = kdimp;
        idum = idimp;
        {
            FortranArray3DRef<double> xslice(&xyzjp(1,1,1), jdum, kdum, idum);
            FortranArray3DRef<double> yslice(&xyzjp(1,1,2), jdum, kdum, idum);
            FortranArray3DRef<double> zslice(&xyzjp(1,1,3), jdum, kdum, idum);
            grdmove_ns::grdmove(nblp, jdum, kdum, idum,
                                xslice, yslice, zslice,
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                dthtx, dthty, dthtz);
        }

        // check for point match
        for (int ii = ista; ii <= iend; ii++) {
            for (int kk = ksta; kk <= kend; kk++) {
                xt1 = x(j,kk,ii);
                yt1 = y(j,kk,ii);
                zt1 = z(j,kk,ii);
                xr1 = xyzjp(kk,ii,1);
                yr1 = xyzjp(kk,ii,2);
                zr1 = xyzjp(kk,ii,3);
                double val = std::sqrt((xr1-xt1)*(xr1-xt1)+(yr1-yt1)*(yr1-yt1)+(zr1-zt1)*(zr1-zt1));
                eps = ccomplex_ns::ccmax(eps, val);
            }
        }

        period_miss(lcnt) = eps;

    }


    // k=1 boundary
    if (nface == 5) {

        k  = 1;
        kp = kdimp;

        // check that idim = idimp and jdim = jdimp
        if (idim != idimp || jdim != jdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " periodic block face must be of the same 2 dimensions (and orientation) as the current block face");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " idim,jdim,idimp,jdimp = %5d%5d%5d%5d", idim, jdim, idimp, jdimp);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // load x,y,z values on periodic face into work array
        for (int ii = 1; ii <= idimp; ii++) {
            for (int jj = 1; jj <= jdimp; jj++) {
                xyzkp(jj,ii,1) = xp(jj,kp,ii);
                xyzkp(jj,ii,2) = yp(jj,kp,ii);
                xyzkp(jj,ii,3) = zp(jj,kp,ii);
            }
        }

        // rotate periodic block face (entire face)
        jdum = jdimp;
        kdum = 1;
        idum = idimp;
        {
            FortranArray3DRef<double> xslice(&xyzkp(1,1,1), jdum, kdum, idum);
            FortranArray3DRef<double> yslice(&xyzkp(1,1,2), jdum, kdum, idum);
            FortranArray3DRef<double> zslice(&xyzkp(1,1,3), jdum, kdum, idum);
            grdmove_ns::grdmove(nblp, jdum, kdum, idum,
                                xslice, yslice, zslice,
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                dthtx, dthty, dthtz);
        }

        // check for point match
        for (int ii = ista; ii <= iend; ii++) {
            for (int jj = jsta; jj <= jend; jj++) {
                xt1 = x(jj,k,ii);
                yt1 = y(jj,k,ii);
                zt1 = z(jj,k,ii);
                xr1 = xyzkp(jj,ii,1);
                yr1 = xyzkp(jj,ii,2);
                zr1 = xyzkp(jj,ii,3);
                double val = std::sqrt((xr1-xt1)*(xr1-xt1)+(yr1-yt1)*(yr1-yt1)+(zr1-zt1)*(zr1-zt1));
                eps = ccomplex_ns::ccmax(eps, val);
            }
        }

        period_miss(lcnt) = eps;

    }

    // k=kdim boundary
    if (nface == 6) {

        k  = kdim;
        kp = 1;

        // check that idim = idimp and jdim = jdimp
        if (idim != idimp || jdim != jdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " periodic block face must be of the same 2 dimensions (and orientation) as the current block face");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " idim,jdim,idimp,jdimp = %5d%5d%5d%5d", idim, jdim, idimp, jdimp);
            int ierrflg = -1;
            termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
        }

        // load x,y,z values on periodic face into work array
        for (int ii = 1; ii <= idimp; ii++) {
            for (int jj = 1; jj <= jdimp; jj++) {
                xyzkp(jj,ii,1) = xp(jj,kp,ii);
                xyzkp(jj,ii,2) = yp(jj,kp,ii);
                xyzkp(jj,ii,3) = zp(jj,kp,ii);
            }
        }

        // rotate periodic block face (entire face)
        jdum = jdimp;
        kdum = 1;
        idum = idimp;
        {
            FortranArray3DRef<double> xslice(&xyzkp(1,1,1), jdum, kdum, idum);
            FortranArray3DRef<double> yslice(&xyzkp(1,1,2), jdum, kdum, idum);
            FortranArray3DRef<double> zslice(&xyzkp(1,1,3), jdum, kdum, idum);
            grdmove_ns::grdmove(nblp, jdum, kdum, idum,
                                xslice, yslice, zslice,
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                xorig(nblp), yorig(nblp), zorig(nblp),
                                dthtx, dthty, dthtz);
        }

        // check for point match
        for (int ii = ista; ii <= iend; ii++) {
            for (int jj = jsta; jj <= jend; jj++) {
                xt1 = x(jj,k,ii);
                yt1 = y(jj,k,ii);
                zt1 = z(jj,k,ii);
                xr1 = xyzkp(jj,ii,1);
                yr1 = xyzkp(jj,ii,2);
                zr1 = xyzkp(jj,ii,3);
                double val = std::sqrt((xr1-xt1)*(xr1-xt1)+(yr1-yt1)*(yr1-yt1)+(zr1-zt1)*(zr1-zt1));
                eps = ccomplex_ns::ccmax(eps, val);
            }
        }

        period_miss(lcnt) = eps;

    }

    return;
}

} // namespace chkrot_ns
