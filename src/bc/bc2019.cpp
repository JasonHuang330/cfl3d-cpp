// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2019.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>

namespace bc2019_ns {

// Inline helpers for ccmincr/ccmaxcr (avoid reference-binding issues)
static inline double ccmincr(double a, double b) { return std::min(a, b); }
static inline double ccmaxcr(double a, double b) { return std::max(a, b); }

// ============================================================
// bc dispatcher
// ============================================================
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
    int iuns = 0;

    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        FortranArray4DRef<double> q(&w(lw(1,ibl)), jdim, kdim, idim, 5);
        FortranArray4DRef<double> qj0(&w(lw(2,ibl)), kdim, idim-1, 5, 4);
        FortranArray4DRef<double> qk0(&w(lw(3,ibl)), jdim, idim-1, 5, 4);
        FortranArray4DRef<double> qi0(&w(lw(4,ibl)), jdim, kdim, 5, 4);
        FortranArray4DRef<double> sj(&w(lw(5,ibl)), jdim, kdim, idim-1, 5);
        FortranArray4DRef<double> sk(&w(lw(6,ibl)), jdim, kdim, idim-1, 5);
        FortranArray4DRef<double> si(&w(lw(7,ibl)), jdim, kdim, idim, 5);
        FortranArray3DRef<double> bcj(&w(lw(8,ibl)), kdim, idim-1, 2);
        FortranArray3DRef<double> bck(&w(lw(9,ibl)), jdim, idim-1, 2);
        FortranArray3DRef<double> bci(&w(lw(10,ibl)), jdim, kdim, 2);

        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim, kdim, idim, nummem);
        FortranArray4DRef<double> tj0(&wk(lw2(2,ibl)), kdim, idim-1, nummem, 4);
        FortranArray4DRef<double> tk0(&wk(lw2(3,ibl)), jdim, idim-1, nummem, 4);
        FortranArray4DRef<double> ti0(&wk(lw2(4,ibl)), jdim, kdim, nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim, kdim, idim);
        FortranArray4DRef<double> vj0(&wk(lw2(6,ibl)), kdim, idim-1, 1, 4);
        FortranArray4DRef<double> vk0(&wk(lw2(7,ibl)), jdim, idim-1, 1, 4);
        FortranArray4DRef<double> vi0(&wk(lw2(8,ibl)), jdim, kdim, 1, 4);

        // j-face segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) == 2019) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                int ldata = lw2(9,ibl);
                FortranArray4DRef<double> bcdata(&wk(ldata), mdim, ndim, 2, 12);
                int ifilno = bcfilej(ibl,iseg,1);
                char (&filname)[80] = bcfiles(ifilno);
                bc2019(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // k-face segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 2019) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                int ldata = lw2(9,ibl);
                FortranArray4DRef<double> bcdata(&wk(ldata), mdim, ndim, 2, 12);
                int ifilno = bcfilek(ibl,iseg,1);
                char (&filname)[80] = bcfiles(ifilno);
                bc2019(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // i-face segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 2019) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                int ldata = lw2(9,ibl);
                FortranArray4DRef<double> bcdata(&wk(ldata), mdim, ndim, 2, 12);
                int ifilno = bcfilei(ibl,iseg,1);
                char (&filname)[80] = bcfiles(ifilno);
                bc2019(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }
    }
}

// ============================================================
// bc2019 subroutine
// ============================================================
void bc2019(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0,
            FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj,
            FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0,
            FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0,
            FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& mdim, int& ndim,
            FortranArray4DRef<double> bcdata,
            char (&filname)[80],
            int& iuns,
            FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim,
            int& myid, int& nummem)
{
    // COMMON block aliases
    float  gamma   = cmn_fluid.gamma;
    float  gm1     = cmn_fluid.gm1;
    float  p0      = cmn_ivals.p0;
    float  rho0    = cmn_ivals.rho0;
    float  xmach   = cmn_info.xmach;
    float* tur10   = cmn_ivals.tur10;   // 0-based
    int*   ivisc   = cmn_reyue.ivisc;   // 0-based
    int    level   = cmn_mgrd.level;
    int    lglobal = cmn_mgrd.lglobal;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // this bc makes use of only one plane of data
    int ip = 1;


    // *** standard boundary condition bctype=2019 ***

    // -----------------------------------------------------------------------
    // j=1 boundary   overflow-type-nozzle total BCs   type 2019
    // -----------------------------------------------------------------------
    if (nface == 3) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,3) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            int js = (i - ista) * (kend - ksta) + 1;
            (void)js;

            for (int k = ksta; k <= kend1; k++) {
                int kk = k - ksta + 1;

                // bcdata(kk,ii,ip,1) is Pt/Ptinf
                // bcdata(kk,ii,ip,2) is Tt/Ttinf

                // like overflow (vrgas=1): freestream
                double finf   = 1. + 0.5 * (gamma - 1.) * xmach * xmach;
                double p0inf  = p0 * std::pow(finf, (double)gamma / (double)gm1);
                double t0inf  = (double)p0 / (double)rho0 * finf;
                double p00 = 0., t00 = 0.;
                if (bcdata(kk,ii,ip,1) != 0.) p00 = p0inf * bcdata(kk,ii,ip,1);
                if (bcdata(kk,ii,ip,2) != 0.) t00 = t0inf * bcdata(kk,ii,ip,2);
                double gm1h  = 0.5 * (gamma - 1.);
                double gdgm1 = (double)gamma / (double)gm1;
                // external conditions
                double rho = qj0(k,i,1,1);
                double u   = qj0(k,i,2,1);
                double v   = qj0(k,i,3,1);
                double w   = qj0(k,i,4,1);
                double p   = qj0(k,i,5,1);
                double v2  = 0.5 * (u*u + v*v + w*w);
                double ei  = p / (rho * gm1);
                double c2  = (double)gamma * (double)gm1 * ei;
                double t   = (double)gm1 * ei;
                double fact = 1. + (double)gm1 * v2 / c2;
                if (bcdata(kk,ii,ip,2) == 0.) t00 = t * fact;
                if (bcdata(kk,ii,ip,1) == 0.) p00 = p * std::pow(fact, gdgm1);
                double t0gr = t00 * gamma;
                // need outward-pointing (- sign for ilo,jlo,klo):
                double uun = -(u*sj(1,k,i,1) + v*sj(1,k,i,2) + w*sj(1,k,i,3) + sj(1,k,i,5));
                double un  = -uun * sj(1,k,i,1);
                double vn  = -uun * sj(1,k,i,2);
                double vw  = -uun * sj(1,k,i,3);
                double wn;  // used but never assigned in Fortran (bug) - uninitialized
                double ut  = u - un;
                double vt  = v - vn;
                double wt  = w - wn;
                double uut2 = ut*ut + vt*vt + wt*wt;
                // internal conditions
                double rhox = q(1,k,i,1);
                double ux   = q(1,k,i,2);
                double vx   = q(1,k,i,3);
                double wx   = q(1,k,i,4);
                double px   = q(1,k,i,5);
                double v2x  = 0.5 * (ux*ux + vx*vx + wx*wx);
                double eix  = px / (rhox * gm1);
                double c2x  = (double)gamma * (double)gm1 * eix;
                double cx   = std::sqrt(c2x);
                // need outward-pointing (- sign for ilo,jlo,klo):
                double uunx = -(ux*sj(1,k,i,1) + vx*sj(1,k,i,2) + wx*sj(1,k,i,3) + sj(1,k,i,5));
                double r2   = uunx + cx / gm1h;
                // flow angle
                double tan2p1;
                if (v2 == 0.) {
                    tan2p1 = 1.;
                } else if (uun == 0.) {
                    tan2p1 = 1.e20;
                } else {
                    tan2p1 = uut2 / (uun*uun) + 1.;
                    tan2p1 = ccmincr(tan2p1, 1.e20);
                }
                // find speed of sound
                double aterm  = 1. + tan2p1 / gm1h;
                double hbtrmm = tan2p1 * r2;
                double cterm  = hbtrmm * gm1h * r2 - t0gr;
                double sqterm = hbtrmm*hbtrmm - aterm*cterm;
                double c      = (hbtrmm + std::sqrt(ccmaxcr(sqterm, 0.))) / aterm;
                c2 = c * c;
                // final values for setting BC
                double uunold = uun;
                uun    = -std::abs(r2 - c / gm1h);
                double v2old  = v2;
                v2     = 0.5 * uun*uun * tan2p1;
                fact   = 1. + (double)gm1 * v2 / c2;
                p      = p00 * std::pow(fact, -gdgm1);
                rho    = (double)gamma * p / c2;
                double vratio;
                if (v2old == 0.) {
                    vratio = uun;
                    u = 1.;
                    v = 0.;
                    w = 0.;
                } else if (uunold == 0. || tan2p1 == 1.e20) {
                    vratio = std::sqrt(v2 / v2old);
                } else {
                    vratio = uun / uunold;
                }

                qj0(k,i,1,1) = rho;
                qj0(k,i,1,2) = qj0(k,i,1,1);
                qj0(k,i,2,1) = u * vratio;
                qj0(k,i,2,2) = qj0(k,i,2,1);
                qj0(k,i,3,1) = v * vratio;
                qj0(k,i,3,2) = qj0(k,i,3,1);
                qj0(k,i,4,1) = w * vratio;
                qj0(k,i,4,2) = qj0(k,i,4,1);
                qj0(k,i,5,1) = p;
                qj0(k,i,5,2) = qj0(k,i,5,1);

                bcj(k,i,1) = 0.0;
            } // k loop
        } // i loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,1) = vist3d(1,k,i);
                    vj0(k,i,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            int kk = k - ksta + 1;
                            double ubar = -(qj0(k,i,2,1)*sj(1,k,i,1) +
                                            qj0(k,i,3,1)*sj(1,k,i,2) +
                                            qj0(k,i,4,1)*sj(1,k,i,3));
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(kk,ii,ip,2+l);
                            if ((float)ubar < 0.f) {
                                tj0(k,i,l,1) = t1;
                                tj0(k,i,l,2) = t1;
                            } else {
                                tj0(k,i,l,1) = tursav(1,k,i,l);
                                tj0(k,i,l,2) = tj0(k,i,l,1);
                            }
                        }
                    }
                }
            }
        }
    } // nface == 3


    // -----------------------------------------------------------------------
    // j=jdim boundary   overflow-type-nozzle total BCs   type 2019
    // -----------------------------------------------------------------------
    if (nface == 4) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,3) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            int js = (i - ista) * (kend - ksta) + 1;
            (void)js;

            for (int k = ksta; k <= kend1; k++) {
                int kk = k - ksta + 1;

                // bcdata(kk,ii,ip,1) is Pt/Ptinf
                // bcdata(kk,ii,ip,2) is Tt/Ttinf

                // like overflow (vrgas=1): freestream
                double finf   = 1. + 0.5 * (gamma - 1.) * xmach * xmach;
                double p0inf  = p0 * std::pow(finf, (double)gamma / (double)gm1);
                double t0inf  = (double)p0 / (double)rho0 * finf;
                double p00 = 0., t00 = 0.;
                if (bcdata(kk,ii,ip,1) != 0.) p00 = p0inf * bcdata(kk,ii,ip,1);
                if (bcdata(kk,ii,ip,2) != 0.) t00 = t0inf * bcdata(kk,ii,ip,2);
                double gm1h  = 0.5 * (gamma - 1.);
                double gdgm1 = (double)gamma / (double)gm1;
                // external conditions
                double rho = qj0(k,i,1,3);
                double u   = qj0(k,i,2,3);
                double v   = qj0(k,i,3,3);
                double w   = qj0(k,i,4,3);
                double p   = qj0(k,i,5,3);
                double v2  = 0.5 * (u*u + v*v + w*w);
                double ei  = p / (rho * gm1);
                double c2  = (double)gamma * (double)gm1 * ei;
                double t   = (double)gm1 * ei;
                double fact = 1. + (double)gm1 * v2 / c2;
                if (bcdata(kk,ii,ip,2) == 0.) t00 = t * fact;
                if (bcdata(kk,ii,ip,1) == 0.) p00 = p * std::pow(fact, gdgm1);
                double t0gr = t00 * gamma;
                // need outward-pointing (+ sign for jhi):
                double uun = +(u*sj(jdim,k,i,1) + v*sj(jdim,k,i,2) + w*sj(jdim,k,i,3) +
                               sj(jdim,k,i,5));
                double un  = +uun * sj(jdim,k,i,1);
                double vn  = +uun * sj(jdim,k,i,2);
                double vw  = +uun * sj(jdim,k,i,3);
                double wn;  // used but never assigned in Fortran (bug) - uninitialized
                double ut  = u - un;
                double vt  = v - vn;
                double wt  = w - wn;
                double uut2 = ut*ut + vt*vt + wt*wt;
                // internal conditions
                double rhox = q(jdim1,k,i,1);
                double ux   = q(jdim1,k,i,2);
                double vx   = q(jdim1,k,i,3);
                double wx   = q(jdim1,k,i,4);
                double px   = q(jdim1,k,i,5);
                double v2x  = 0.5 * (ux*ux + vx*vx + wx*wx);
                double eix  = px / (rhox * gm1);
                double c2x  = (double)gamma * (double)gm1 * eix;
                double cx   = std::sqrt(c2x);
                // need outward-pointing (+ sign for jhi):
                double uunx = +(ux*sj(jdim,k,i,1) + vx*sj(jdim,k,i,2) + wx*sj(jdim,k,i,3) +
                                sj(jdim,k,i,5));
                double r2   = uunx + cx / gm1h;
                // flow angle
                double tan2p1;
                if (v2 == 0.) {
                    tan2p1 = 1.;
                } else if (uun == 0.) {
                    tan2p1 = 1.e20;
                } else {
                    tan2p1 = uut2 / (uun*uun) + 1.;
                    tan2p1 = ccmincr(tan2p1, 1.e20);
                }
                // find speed of sound
                double aterm  = 1. + tan2p1 / gm1h;
                double hbtrmm = tan2p1 * r2;
                double cterm  = hbtrmm * gm1h * r2 - t0gr;
                double sqterm = hbtrmm*hbtrmm - aterm*cterm;
                double c      = (hbtrmm + std::sqrt(ccmaxcr(sqterm, 0.))) / aterm;
                c2 = c * c;
                // final values for setting BC
                double uunold = uun;
                uun    = -std::abs(r2 - c / gm1h);
                double v2old  = v2;
                v2     = 0.5 * uun*uun * tan2p1;
                fact   = 1. + (double)gm1 * v2 / c2;
                p      = p00 * std::pow(fact, -gdgm1);
                rho    = (double)gamma * p / c2;
                double vratio;
                if (v2old == 0.) {
                    vratio = uun;
                    u = 1.;
                    v = 0.;
                    w = 0.;
                } else if (uunold == 0. || tan2p1 == 1.e20) {
                    vratio = std::sqrt(v2 / v2old);
                } else {
                    vratio = uun / uunold;
                }

                qj0(k,i,1,3) = rho;
                qj0(k,i,1,4) = qj0(k,i,1,3);
                qj0(k,i,2,3) = u * vratio;
                qj0(k,i,2,4) = qj0(k,i,2,3);
                qj0(k,i,3,3) = v * vratio;
                qj0(k,i,3,4) = qj0(k,i,3,3);
                qj0(k,i,4,3) = w * vratio;
                qj0(k,i,4,4) = qj0(k,i,4,3);
                qj0(k,i,5,3) = p;
                qj0(k,i,5,4) = qj0(k,i,5,3);

                bcj(k,i,2) = 0.0;
            } // k loop
        } // i loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,3) = vist3d(jdim-1,k,i);
                    vj0(k,i,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            int kk = k - ksta + 1;
                            double ubar = qj0(k,i,2,3)*sj(jdim,k,i,1) +
                                          qj0(k,i,3,3)*sj(jdim,k,i,2) +
                                          qj0(k,i,4,3)*sj(jdim,k,i,3);
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(kk,ii,ip,2+l);
                            if ((float)ubar < 0.f) {
                                tj0(k,i,l,3) = t1;
                                tj0(k,i,l,4) = t1;
                            } else {
                                tj0(k,i,l,3) = tursav(jdim-1,k,i,l);
                                tj0(k,i,l,4) = tj0(k,i,l,3);
                            }
                        }
                    }
                }
            }
        }
    } // nface == 4


    // -----------------------------------------------------------------------
    // k=1 boundary   overflow-type-nozzle total BCs   type 2019
    // -----------------------------------------------------------------------
    if (nface == 5) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,3) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            int js = (i - ista) * (jend - jsta) + 1;
            (void)js;

            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;

                // bcdata(jj,ii,ip,1) is Pt/Ptinf
                // bcdata(jj,ii,ip,2) is Tt/Ttinf

                // like overflow (vrgas=1): freestream
                double finf   = 1. + 0.5 * (gamma - 1.) * xmach * xmach;
                double p0inf  = p0 * std::pow(finf, (double)gamma / (double)gm1);
                double t0inf  = (double)p0 / (double)rho0 * finf;
                double p00 = 0., t00 = 0.;
                if (bcdata(jj,ii,ip,1) != 0.) p00 = p0inf * bcdata(jj,ii,ip,1);
                if (bcdata(jj,ii,ip,2) != 0.) t00 = t0inf * bcdata(jj,ii,ip,2);
                double gm1h  = 0.5 * (gamma - 1.);
                double gdgm1 = (double)gamma / (double)gm1;
                // external conditions
                double rho = qk0(j,i,1,1);
                double u   = qk0(j,i,2,1);
                double v   = qk0(j,i,3,1);
                double w   = qk0(j,i,4,1);
                double p   = qk0(j,i,5,1);
                double v2  = 0.5 * (u*u + v*v + w*w);
                double ei  = p / (rho * gm1);
                double c2  = (double)gamma * (double)gm1 * ei;
                double t   = (double)gm1 * ei;
                double fact = 1. + (double)gm1 * v2 / c2;
                if (bcdata(jj,ii,ip,2) == 0.) t00 = t * fact;
                if (bcdata(jj,ii,ip,1) == 0.) p00 = p * std::pow(fact, gdgm1);
                double t0gr = t00 * gamma;
                // need outward-pointing (- sign for klo):
                double uun = -(u*sk(j,1,i,1) + v*sk(j,1,i,2) + w*sk(j,1,i,3) + sk(j,1,i,5));
                double un  = -uun * sk(j,1,i,1);
                double vn  = -uun * sk(j,1,i,2);
                double vw  = -uun * sk(j,1,i,3);
                double wn;  // used but never assigned in Fortran (bug) - uninitialized
                double ut  = u - un;
                double vt  = v - vn;
                double wt  = w - wn;
                double uut2 = ut*ut + vt*vt + wt*wt;
                // internal conditions
                double rhox = q(j,1,i,1);
                double ux   = q(j,1,i,2);
                double vx   = q(j,1,i,3);
                double wx   = q(j,1,i,4);
                double px   = q(j,1,i,5);
                double v2x  = 0.5 * (ux*ux + vx*vx + wx*wx);
                double eix  = px / (rhox * gm1);
                double c2x  = (double)gamma * (double)gm1 * eix;
                double cx   = std::sqrt(c2x);
                // need outward-pointing (- sign for klo):
                double uunx = -(ux*sk(j,1,i,1) + vx*sk(j,1,i,2) + wx*sk(j,1,i,3) + sk(j,1,i,5));
                double r2   = uunx + cx / gm1h;
                // flow angle
                double tan2p1;
                if (v2 == 0.) {
                    tan2p1 = 1.;
                } else if (uun == 0.) {
                    tan2p1 = 1.e20;
                } else {
                    tan2p1 = uut2 / (uun*uun) + 1.;
                    tan2p1 = ccmincr(tan2p1, 1.e20);
                }
                // find speed of sound
                double aterm  = 1. + tan2p1 / gm1h;
                double hbtrmm = tan2p1 * r2;
                double cterm  = hbtrmm * gm1h * r2 - t0gr;
                double sqterm = hbtrmm*hbtrmm - aterm*cterm;
                double c      = (hbtrmm + std::sqrt(ccmaxcr(sqterm, 0.))) / aterm;
                c2 = c * c;
                // final values for setting BC
                double uunold = uun;
                uun    = -std::abs(r2 - c / gm1h);
                double v2old  = v2;
                v2     = 0.5 * uun*uun * tan2p1;
                fact   = 1. + (double)gm1 * v2 / c2;
                p      = p00 * std::pow(fact, -gdgm1);
                rho    = (double)gamma * p / c2;
                double vratio;
                if (v2old == 0.) {
                    vratio = uun;
                    u = 1.;
                    v = 0.;
                    w = 0.;
                } else if (uunold == 0. || tan2p1 == 1.e20) {
                    vratio = std::sqrt(v2 / v2old);
                } else {
                    vratio = uun / uunold;
                }

                qk0(j,i,1,1) = rho;
                qk0(j,i,1,2) = qk0(j,i,1,1);
                qk0(j,i,2,1) = u * vratio;
                qk0(j,i,2,2) = qk0(j,i,2,1);
                qk0(j,i,3,1) = v * vratio;
                qk0(j,i,3,2) = qk0(j,i,3,1);
                qk0(j,i,4,1) = w * vratio;
                qk0(j,i,4,2) = qk0(j,i,4,1);
                qk0(j,i,5,1) = p;
                qk0(j,i,5,2) = qk0(j,i,5,1);

                bck(j,i,1) = 0.0;
            } // j loop
        } // i loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,1) = vist3d(j,1,i);
                    vk0(j,i,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            double ubar = -(qk0(j,i,2,1)*sk(j,1,i,1) +
                                            qk0(j,i,3,1)*sk(j,1,i,2) +
                                            qk0(j,i,4,1)*sk(j,1,i,3));
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(jj,ii,ip,2+l);
                            if ((float)ubar < 0.f) {
                                tk0(j,i,l,1) = t1;
                                tk0(j,i,l,2) = t1;
                            } else {
                                tk0(j,i,l,1) = tursav(j,1,i,l);
                                tk0(j,i,l,2) = tk0(j,i,l,1);
                            }
                        }
                    }
                }
            }
        }
    } // nface == 5


    // -----------------------------------------------------------------------
    // k=kdim boundary   overflow-type-nozzle total BCs   type 2019
    // -----------------------------------------------------------------------
    if (nface == 6) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,3) > -1.e10f) itrflg1 = 1;

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            int js = (i - ista) * (jend - jsta) + 1;
            (void)js;

            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;

                // bcdata(jj,ii,ip,1) is Pt/Ptinf
                // bcdata(jj,ii,ip,2) is Tt/Ttinf

                // like overflow (vrgas=1): freestream
                double finf   = 1. + 0.5 * (gamma - 1.) * xmach * xmach;
                double p0inf  = p0 * std::pow(finf, (double)gamma / (double)gm1);
                double t0inf  = (double)p0 / (double)rho0 * finf;
                double p00 = 0., t00 = 0.;
                if (bcdata(jj,ii,ip,1) != 0.) p00 = p0inf * bcdata(jj,ii,ip,1);
                if (bcdata(jj,ii,ip,2) != 0.) t00 = t0inf * bcdata(jj,ii,ip,2);
                double gm1h  = 0.5 * (gamma - 1.);
                double gdgm1 = (double)gamma / (double)gm1;
                // external conditions
                double rho = qk0(j,i,1,3);
                double u   = qk0(j,i,2,3);
                double v   = qk0(j,i,3,3);
                double w   = qk0(j,i,4,3);
                double p   = qk0(j,i,5,3);
                double v2  = 0.5 * (u*u + v*v + w*w);
                double ei  = p / (rho * gm1);
                double c2  = (double)gamma * (double)gm1 * ei;
                double t   = (double)gm1 * ei;
                double fact = 1. + (double)gm1 * v2 / c2;
                if (bcdata(jj,ii,ip,2) == 0.) t00 = t * fact;
                if (bcdata(jj,ii,ip,1) == 0.) p00 = p * std::pow(fact, gdgm1);
                double t0gr = t00 * gamma;
                // need outward-pointing (+ sign for khi):
                double uun = +(u*sk(j,kdim,i,1) + v*sk(j,kdim,i,2) + w*sk(j,kdim,i,3) +
                               sk(j,kdim,i,5));
                double un  = +uun * sk(j,kdim,i,1);
                double vn  = +uun * sk(j,kdim,i,2);
                double vw  = +uun * sk(j,kdim,i,3);
                double wn;  // used but never assigned in Fortran (bug) - uninitialized
                double ut  = u - un;
                double vt  = v - vn;
                double wt  = w - wn;
                double uut2 = ut*ut + vt*vt + wt*wt;
                // internal conditions
                double rhox = q(j,kdim1,i,1);
                double ux   = q(j,kdim1,i,2);
                double vx   = q(j,kdim1,i,3);
                double wx   = q(j,kdim1,i,4);
                double px   = q(j,kdim1,i,5);
                double v2x  = 0.5 * (ux*ux + vx*vx + wx*wx);
                double eix  = px / (rhox * gm1);
                double c2x  = (double)gamma * (double)gm1 * eix;
                double cx   = std::sqrt(c2x);
                // need outward-pointing (+ sign for khi):
                double uunx = +(ux*sk(j,kdim,i,1) + vx*sk(j,kdim,i,2) + wx*sk(j,kdim,i,3) +
                                sk(j,kdim,i,5));
                double r2   = uunx + cx / gm1h;
                // flow angle
                double tan2p1;
                if (v2 == 0.) {
                    tan2p1 = 1.;
                } else if (uun == 0.) {
                    tan2p1 = 1.e20;
                } else {
                    tan2p1 = uut2 / (uun*uun) + 1.;
                    tan2p1 = ccmincr(tan2p1, 1.e20);
                }
                // find speed of sound
                double aterm  = 1. + tan2p1 / gm1h;
                double hbtrmm = tan2p1 * r2;
                double cterm  = hbtrmm * gm1h * r2 - t0gr;
                double sqterm = hbtrmm*hbtrmm - aterm*cterm;
                double c      = (hbtrmm + std::sqrt(ccmaxcr(sqterm, 0.))) / aterm;
                c2 = c * c;
                // final values for setting BC
                double uunold = uun;
                uun    = -std::abs(r2 - c / gm1h);
                double v2old  = v2;
                v2     = 0.5 * uun*uun * tan2p1;
                fact   = 1. + (double)gm1 * v2 / c2;
                p      = p00 * std::pow(fact, -gdgm1);
                rho    = (double)gamma * p / c2;
                double vratio;
                if (v2old == 0.) {
                    vratio = uun;
                    u = 1.;
                    v = 0.;
                    w = 0.;
                } else if (uunold == 0. || tan2p1 == 1.e20) {
                    vratio = std::sqrt(v2 / v2old);
                } else {
                    vratio = uun / uunold;
                }

                qk0(j,i,1,3) = rho;
                qk0(j,i,1,4) = qk0(j,i,1,3);
                qk0(j,i,2,3) = u * vratio;
                qk0(j,i,2,4) = qk0(j,i,2,3);
                qk0(j,i,3,3) = v * vratio;
                qk0(j,i,3,4) = qk0(j,i,3,3);
                qk0(j,i,4,3) = w * vratio;
                qk0(j,i,4,4) = qk0(j,i,4,3);
                qk0(j,i,5,3) = p;
                qk0(j,i,5,4) = qk0(j,i,5,3);

                bck(j,i,2) = 0.0;
            } // j loop
        } // i loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,3) = vist3d(j,kdim-1,i);
                    vk0(j,i,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            double ubar = qk0(j,i,2,3)*sk(j,kdim,i,1) +
                                          qk0(j,i,3,3)*sk(j,kdim,i,2) +
                                          qk0(j,i,4,3)*sk(j,kdim,i,3);
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(jj,ii,ip,2+l);
                            if ((float)ubar < 0.f) {
                                tk0(j,i,l,3) = t1;
                                tk0(j,i,l,4) = t1;
                            } else {
                                tk0(j,i,l,3) = tursav(j,kdim-1,i,l);
                                tk0(j,i,l,4) = tk0(j,i,l,3);
                            }
                        }
                    }
                }
            }
        }
    } // nface == 6


    // -----------------------------------------------------------------------
    // i=1 boundary   overflow-type-nozzle total BCs   type 2019
    // -----------------------------------------------------------------------
    if (nface == 1) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,3) > -1.e10f) itrflg1 = 1;

        for (int k = ksta; k <= kend1; k++) {
            int kk = k - ksta + 1;
            int js = (k - ksta) * (jend - jsta) + 1;
            (void)js;

            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;

                // bcdata(jj,kk,ip,1) is Pt/Ptinf
                // bcdata(jj,kk,ip,2) is Tt/Ttinf

                // like overflow (vrgas=1): freestream
                double finf   = 1. + 0.5 * (gamma - 1.) * xmach * xmach;
                double p0inf  = p0 * std::pow(finf, (double)gamma / (double)gm1);
                double t0inf  = (double)p0 / (double)rho0 * finf;
                double p00 = 0., t00 = 0.;
                if (bcdata(jj,kk,ip,1) != 0.) p00 = p0inf * bcdata(jj,kk,ip,1);
                if (bcdata(jj,kk,ip,2) != 0.) t00 = t0inf * bcdata(jj,kk,ip,2);
                double gm1h  = 0.5 * (gamma - 1.);
                double gdgm1 = (double)gamma / (double)gm1;
                // external conditions
                double rho = qi0(j,k,1,1);
                double u   = qi0(j,k,2,1);
                double v   = qi0(j,k,3,1);
                double w   = qi0(j,k,4,1);
                double p   = qi0(j,k,5,1);
                double v2  = 0.5 * (u*u + v*v + w*w);
                double ei  = p / (rho * gm1);
                double c2  = (double)gamma * (double)gm1 * ei;
                double t   = (double)gm1 * ei;
                double fact = 1. + (double)gm1 * v2 / c2;
                if (bcdata(jj,kk,ip,2) == 0.) t00 = t * fact;
                if (bcdata(jj,kk,ip,1) == 0.) p00 = p * std::pow(fact, gdgm1);
                double t0gr = t00 * gamma;
                // need outward-pointing (- sign for ilo):
                double uun = -(u*si(j,k,1,1) + v*si(j,k,1,2) + w*si(j,k,1,3) + si(j,k,1,5));
                double un  = -uun * si(j,k,1,1);
                double vn  = -uun * si(j,k,1,2);
                double vw  = -uun * si(j,k,1,3);
                double wn;  // used but never assigned in Fortran (bug) - uninitialized
                double ut  = u - un;
                double vt  = v - vn;
                double wt  = w - wn;
                double uut2 = ut*ut + vt*vt + wt*wt;
                // internal conditions
                double rhox = q(j,k,1,1);
                double ux   = q(j,k,1,2);
                double vx   = q(j,k,1,3);
                double wx   = q(j,k,1,4);
                double px   = q(j,k,1,5);
                double v2x  = 0.5 * (ux*ux + vx*vx + wx*wx);
                double eix  = px / (rhox * gm1);
                double c2x  = (double)gamma * (double)gm1 * eix;
                double cx   = std::sqrt(c2x);
                // need outward-pointing (- sign for ilo):
                double uunx = -(ux*si(j,k,1,1) + vx*si(j,k,1,2) + wx*si(j,k,1,3) + si(j,k,1,5));
                double r2   = uunx + cx / gm1h;
                // flow angle
                double tan2p1;
                if (v2 == 0.) {
                    tan2p1 = 1.;
                } else if (uun == 0.) {
                    tan2p1 = 1.e20;
                } else {
                    tan2p1 = uut2 / (uun*uun) + 1.;
                    tan2p1 = ccmincr(tan2p1, 1.e20);
                }
                // find speed of sound
                double aterm  = 1. + tan2p1 / gm1h;
                double hbtrmm = tan2p1 * r2;
                double cterm  = hbtrmm * gm1h * r2 - t0gr;
                double sqterm = hbtrmm*hbtrmm - aterm*cterm;
                double c      = (hbtrmm + std::sqrt(ccmaxcr(sqterm, 0.))) / aterm;
                c2 = c * c;
                // final values for setting BC
                double uunold = uun;
                uun    = -std::abs(r2 - c / gm1h);
                double v2old  = v2;
                v2     = 0.5 * uun*uun * tan2p1;
                fact   = 1. + (double)gm1 * v2 / c2;
                p      = p00 * std::pow(fact, -gdgm1);
                rho    = (double)gamma * p / c2;
                double vratio;
                if (v2old == 0.) {
                    vratio = uun;
                    u = 1.;
                    v = 0.;
                    w = 0.;
                } else if (uunold == 0. || tan2p1 == 1.e20) {
                    vratio = std::sqrt(v2 / v2old);
                } else {
                    vratio = uun / uunold;
                }

                qi0(j,k,1,1) = rho;
                qi0(j,k,1,2) = qi0(j,k,1,1);
                qi0(j,k,2,1) = u * vratio;
                qi0(j,k,2,2) = qi0(j,k,2,1);
                qi0(j,k,3,1) = v * vratio;
                qi0(j,k,3,2) = qi0(j,k,3,1);
                qi0(j,k,4,1) = w * vratio;
                qi0(j,k,4,2) = qi0(j,k,4,1);
                qi0(j,k,5,1) = p;
                qi0(j,k,5,2) = qi0(j,k,5,1);

                bci(j,k,1) = 0.0;
            } // j loop
        } // k loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,1) = vist3d(j,k,1);
                    vi0(j,k,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            double ubar = -(qi0(j,k,2,1)*si(j,k,1,1) +
                                            qi0(j,k,3,1)*si(j,k,1,2) +
                                            qi0(j,k,4,1)*si(j,k,1,3));
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(jj,kk,ip,2+l);
                            if ((float)ubar < 0.f) {
                                ti0(j,k,l,1) = t1;
                                ti0(j,k,l,2) = t1;
                            } else {
                                ti0(j,k,l,1) = tursav(j,k,1,l);
                                ti0(j,k,l,2) = ti0(j,k,l,1);
                            }
                        }
                    }
                }
            }
        }
    } // nface == 1


    // -----------------------------------------------------------------------
    // i=idim boundary   overflow-type-nozzle total BCs   type 2019
    // -----------------------------------------------------------------------
    if (nface == 2) {
        // check turbulence data flag
        int itrflg1 = 0;
        if ((float)bcdata(1,1,ip,3) > -1.e10f) itrflg1 = 1;

        for (int k = ksta; k <= kend1; k++) {
            int kk = k - ksta + 1;
            int js = (k - ksta) * (jend - jsta) + 1;
            (void)js;

            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;

                // bcdata(jj,kk,ip,1) is Pt/Ptinf
                // bcdata(jj,kk,ip,2) is Tt/Ttinf

                // like overflow (vrgas=1): freestream
                double finf   = 1. + 0.5 * (gamma - 1.) * xmach * xmach;
                double p0inf  = p0 * std::pow(finf, (double)gamma / (double)gm1);
                double t0inf  = (double)p0 / (double)rho0 * finf;
                double p00 = 0., t00 = 0.;
                if (bcdata(jj,kk,ip,1) != 0.) p00 = p0inf * bcdata(jj,kk,ip,1);
                if (bcdata(jj,kk,ip,2) != 0.) t00 = t0inf * bcdata(jj,kk,ip,2);
                double gm1h  = 0.5 * (gamma - 1.);
                double gdgm1 = (double)gamma / (double)gm1;
                // external conditions
                double rho = qi0(j,k,1,3);
                double u   = qi0(j,k,2,3);
                double v   = qi0(j,k,3,3);
                double w   = qi0(j,k,4,3);
                double p   = qi0(j,k,5,3);
                double v2  = 0.5 * (u*u + v*v + w*w);
                double ei  = p / (rho * gm1);
                double c2  = (double)gamma * (double)gm1 * ei;
                double t   = (double)gm1 * ei;
                double fact = 1. + (double)gm1 * v2 / c2;
                if (bcdata(jj,kk,ip,2) == 0.) t00 = t * fact;
                if (bcdata(jj,kk,ip,1) == 0.) p00 = p * std::pow(fact, gdgm1);
                double t0gr = t00 * gamma;
                // need outward-pointing (+ sign for ihi):
                double uun = +(u*si(j,k,idim,1) + v*si(j,k,idim,2) + w*si(j,k,idim,3) +
                               si(j,k,idim,5));
                double un  = +uun * si(j,k,idim,1);
                double vn  = +uun * si(j,k,idim,2);
                double vw  = +uun * si(j,k,idim,3);
                double wn;  // used but never assigned in Fortran (bug) - uninitialized
                double ut  = u - un;
                double vt  = v - vn;
                double wt  = w - wn;
                double uut2 = ut*ut + vt*vt + wt*wt;
                // internal conditions
                double rhox = q(j,k,idim1,1);
                double ux   = q(j,k,idim1,2);
                double vx   = q(j,k,idim1,3);
                double wx   = q(j,k,idim1,4);
                double px   = q(j,k,idim1,5);
                double v2x  = 0.5 * (ux*ux + vx*vx + wx*wx);
                double eix  = px / (rhox * gm1);
                double c2x  = (double)gamma * (double)gm1 * eix;
                double cx   = std::sqrt(c2x);
                // need outward-pointing (+ sign for ihi):
                double uunx = +(ux*si(j,k,idim,1) + vx*si(j,k,idim,2) + wx*si(j,k,idim,3) +
                                si(j,k,idim,5));
                double r2   = uunx + cx / gm1h;
                // flow angle
                double tan2p1;
                if (v2 == 0.) {
                    tan2p1 = 1.;
                } else if (uun == 0.) {
                    tan2p1 = 1.e20;
                } else {
                    tan2p1 = uut2 / (uun*uun) + 1.;
                    tan2p1 = ccmincr(tan2p1, 1.e20);
                }
                // find speed of sound
                double aterm  = 1. + tan2p1 / gm1h;
                double hbtrmm = tan2p1 * r2;
                double cterm  = hbtrmm * gm1h * r2 - t0gr;
                double sqterm = hbtrmm*hbtrmm - aterm*cterm;
                double c      = (hbtrmm + std::sqrt(ccmaxcr(sqterm, 0.))) / aterm;
                c2 = c * c;
                // final values for setting BC
                double uunold = uun;
                uun    = -std::abs(r2 - c / gm1h);
                double v2old  = v2;
                v2     = 0.5 * uun*uun * tan2p1;
                fact   = 1. + (double)gm1 * v2 / c2;
                p      = p00 * std::pow(fact, -gdgm1);
                rho    = (double)gamma * p / c2;
                double vratio;
                if (v2old == 0.) {
                    vratio = uun;
                    u = 1.;
                    v = 0.;
                    w = 0.;
                } else if (uunold == 0. || tan2p1 == 1.e20) {
                    vratio = std::sqrt(v2 / v2old);
                } else {
                    vratio = uun / uunold;
                }

                qi0(j,k,1,3) = rho;
                qi0(j,k,1,4) = qi0(j,k,1,3);
                qi0(j,k,2,3) = u * vratio;
                qi0(j,k,2,4) = qi0(j,k,2,3);
                qi0(j,k,3,3) = v * vratio;
                qi0(j,k,3,4) = qi0(j,k,3,3);
                qi0(j,k,4,3) = w * vratio;
                qi0(j,k,4,4) = qi0(j,k,4,3);
                qi0(j,k,5,3) = p;
                qi0(j,k,5,4) = qi0(j,k,5,3);

                bci(j,k,2) = 0.0;
            } // j loop
        } // k loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,3) = vist3d(j,k,idim-1);
                    vi0(j,k,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            double ubar = qi0(j,k,2,3)*si(j,k,idim,1) +
                                          qi0(j,k,3,3)*si(j,k,idim,2) +
                                          qi0(j,k,4,3)*si(j,k,idim,3);
                            double t1 = (1 - itrflg1) * (double)tur10[l-1]
                                      + itrflg1 * bcdata(jj,kk,ip,2+l);
                            if ((float)ubar < 0.f) {
                                ti0(j,k,l,3) = t1;
                                ti0(j,k,l,4) = t1;
                            } else {
                                ti0(j,k,l,3) = tursav(j,k,idim-1,l);
                                ti0(j,k,l,4) = ti0(j,k,l,3);
                            }
                        }
                    }
                }
            }
        }
    } // nface == 2

} // bc2019

} // namespace bc2019_ns
