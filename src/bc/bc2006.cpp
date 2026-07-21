// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2006.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace bc2006_ns {

// -----------------------------------------------------------------------
// bc dispatcher for bctype=2006
// -----------------------------------------------------------------------
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

        FortranArray4DRef<double> q(   &w(lw(1,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray4DRef<double> qj0( &w(lw(2,ibl)),  kdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qk0( &w(lw(3,ibl)),  jdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qi0( &w(lw(4,ibl)),  jdim,   kdim,   5,      4);
        FortranArray4DRef<double> sj(  &w(lw(5,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> sk(  &w(lw(6,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> si(  &w(lw(7,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray3DRef<double> bcj( &w(lw(8,ibl)),  kdim,   idim-1, 2);
        FortranArray3DRef<double> bck( &w(lw(9,ibl)),  jdim,   idim-1, 2);
        FortranArray3DRef<double> bci( &w(lw(10,ibl)), jdim,   kdim,   2);
        FortranArray4DRef<double> xtbj(&w(lw(11,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbk(&w(lw(12,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbi(&w(lw(13,ibl)), jdim,   kdim,   3,      2);
        FortranArray4DRef<double> atbj(&w(lw(14,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbk(&w(lw(15,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbi(&w(lw(16,ibl)), jdim,   kdim,   3,      2);

        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0(   &wk(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0(   &wk(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0(   &wk(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0(   &wk(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0(   &wk(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0(   &wk(lw2(8,ibl)), jdim,   kdim,   1,      4);

        FortranArray3DRef<double> x(&w(lw(17,ibl)), jdim, kdim, idim);
        FortranArray3DRef<double> y(&w(lw(18,ibl)), jdim, kdim, idim);
        FortranArray3DRef<double> z(&w(lw(19,ibl)), jdim, kdim, idim);

        // j-face segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) == 2006) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                int lbcdata = lw2(9,ibl);
                FortranArray4DRef<double> bcdata(&wk(lbcdata), mdim, ndim, 2, 12);
                int ifilj = bcfilej(ibl,iseg,1);
                char filname[80];
                std::memcpy(filname, bcfiles(ifilj), 80);
                // get nblc from bcdata
                int nblc = (int)(float)bcdata(1,1,1,1);
                int jdimc = 1, kdimc = 1, idimc = 1;
                FortranArray4DRef<double> qj0c(nullptr, 1, 1, 5, 4);
                FortranArray4DRef<double> qk0c(nullptr, 1, 1, 5, 4);
                FortranArray4DRef<double> qi0c(nullptr, 1, 1, 5, 4);
                int ichk = 0;
                if (nblc != 0) {
                    jdimc = jdimg(nblc);
                    kdimc = kdimg(nblc);
                    idimc = idimg(nblc);
                    qj0c = FortranArray4DRef<double>(&wk(lw2(2,nblc)), kdimc, idimc-1, 5, 4);
                    qk0c = FortranArray4DRef<double>(&wk(lw2(3,nblc)), jdimc, idimc-1, 5, 4);
                    qi0c = FortranArray4DRef<double>(&wk(lw2(4,nblc)), jdimc, kdimc,   5, 4);
                    ichk = 1;
                }
                bc2006(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       x, y, z, nblc, jdimc, kdimc, idimc,
                       qj0c, qk0c, qi0c, ichk,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // k-face segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 2006) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                int lbcdata = lw2(9,ibl);
                FortranArray4DRef<double> bcdata(&wk(lbcdata), mdim, ndim, 2, 12);
                int ifilk = bcfilek(ibl,iseg,1);
                char filname[80];
                std::memcpy(filname, bcfiles(ifilk), 80);
                int nblc = (int)(float)bcdata(1,1,1,1);
                int jdimc = 1, kdimc = 1, idimc = 1;
                FortranArray4DRef<double> qj0c(nullptr, 1, 1, 5, 4);
                FortranArray4DRef<double> qk0c(nullptr, 1, 1, 5, 4);
                FortranArray4DRef<double> qi0c(nullptr, 1, 1, 5, 4);
                int ichk = 0;
                if (nblc != 0) {
                    jdimc = jdimg(nblc);
                    kdimc = kdimg(nblc);
                    idimc = idimg(nblc);
                    qj0c = FortranArray4DRef<double>(&wk(lw2(2,nblc)), kdimc, idimc-1, 5, 4);
                    qk0c = FortranArray4DRef<double>(&wk(lw2(3,nblc)), jdimc, idimc-1, 5, 4);
                    qi0c = FortranArray4DRef<double>(&wk(lw2(4,nblc)), jdimc, kdimc,   5, 4);
                    ichk = 1;
                }
                bc2006(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       x, y, z, nblc, jdimc, kdimc, idimc,
                       qj0c, qk0c, qi0c, ichk,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // i-face segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 2006) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                int lbcdata = lw2(9,ibl);
                FortranArray4DRef<double> bcdata(&wk(lbcdata), mdim, ndim, 2, 12);
                int ifili = bcfilei(ibl,iseg,1);
                char filname[80];
                std::memcpy(filname, bcfiles(ifili), 80);
                int nblc = (int)(float)bcdata(1,1,1,1);
                int jdimc = 1, kdimc = 1, idimc = 1;
                FortranArray4DRef<double> qj0c(nullptr, 1, 1, 5, 4);
                FortranArray4DRef<double> qk0c(nullptr, 1, 1, 5, 4);
                FortranArray4DRef<double> qi0c(nullptr, 1, 1, 5, 4);
                int ichk = 0;
                if (nblc != 0) {
                    jdimc = jdimg(nblc);
                    kdimc = kdimg(nblc);
                    idimc = idimg(nblc);
                    qj0c = FortranArray4DRef<double>(&wk(lw2(2,nblc)), kdimc, idimc-1, 5, 4);
                    qk0c = FortranArray4DRef<double>(&wk(lw2(3,nblc)), jdimc, idimc-1, 5, 4);
                    qi0c = FortranArray4DRef<double>(&wk(lw2(4,nblc)), jdimc, kdimc,   5, 4);
                    ichk = 1;
                }
                bc2006(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       x, y, z, nblc, jdimc, kdimc, idimc,
                       qj0c, qk0c, qi0c, ichk,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }
    } // end ibl loop
}

// -----------------------------------------------------------------------
// bc2006: radial equilibrium pressure BC
// -----------------------------------------------------------------------
void bc2006(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj, FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk,
            FortranArray4DRef<double> xtbi,
            FortranArray4DRef<double> atbj, FortranArray4DRef<double> atbk,
            FortranArray4DRef<double> atbi,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& mdim, int& ndim,
            FortranArray4DRef<double> bcdata,
            char (&filname)[80], int& iuns,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y,
            FortranArray3DRef<double> z,
            int& nblc, int& jdimc, int& kdimc, int& idimc,
            FortranArray4DRef<double> qj0c, FortranArray4DRef<double> qk0c,
            FortranArray4DRef<double> qi0c,
            int& ichk,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid, int& nummem)
{
    // COMMON block references
    float& p0      = cmn_ivals.p0;
    int&   level   = cmn_mgrd.level;
    int&   lglobal = cmn_mgrd.lglobal;
    int*   ivisc   = cmn_reyue.ivisc;  // 0-based

    // local dimension helpers
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // only one plane of data
    int ipp = 1;

    // integration direction and sense
    int  lijk = std::abs((int)(float)bcdata(1,1,1,3));
    int  ldir = (int)(float)bcdata(1,1,1,3) /
                std::abs((int)(float)bcdata(1,1,1,3));

    double xfact = 1.0, yfact = 1.0, zfact = 1.0;
    if (bcdata(1,1,1,4) == 1.0) {
        xfact = 0.0;
    } else if (bcdata(1,1,1,4) == 2.0) {
        yfact = 0.0;
    } else if (bcdata(1,1,1,4) == 3.0) {
        zfact = 0.0;
    } else {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120,
            " stopping...bcdata(4) = %g  should be +1, +2, or +3",
            bcdata(1,1,1,4));
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }



    // =========================================================
    // j=1 boundary   bctype 2006
    // =========================================================
    if (nface == 3) {
        int j = 1;

        // extrapolate all primitives except pressure
        for (int l = 1; l <= 4; l++)
            for (int i = ista; i <= iend1; i++)
                for (int k = ksta; k <= kend1; k++) {
                    qj0(k,i,l,1) = q(1,k,i,l);
                    qj0(k,i,l,2) = qj0(k,i,l,1);
                    bcj(k,i,1)   = 0.0;
                }

        if (lijk == 1) {
            // radial integration in i-direction (k is theta-direction)
            int i0, is, ie, ii0, iinc, iic;
            if (ldir > 0) {
                i0   = ista;
                is   = ista + 1;
                ie   = iend;
                ii0  = 1;
                iinc = 1;
                iic  = idimc - 1;
            } else {
                i0   = iend;
                is   = iend - 1;
                ie   = ista;
                ii0  = iend - ista;
                iinc = 0;
                iic  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        bcdata(kk,ii0,ipp,2) = qj0c(k,iic,5,1) / (double)p0;
                    }
                } else {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        bcdata(kk,ii0,ipp,2) = q(j,k,ii0,5) / (double)p0;
                    }
                }
            }

            for (int k = ksta; k <= kend1; k++) {
                int kk  = k - ksta + 1;
                double pm1 = bcdata(kk,ii0,ipp,2) * (double)p0;
                int i   = i0;
                int kp  = k  + 1;
                int jp  = j  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                  + x(jp,k,i) + x(jp,kp,i));
                double yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                  + y(jp,k,i) + y(jp,kp,i));
                double zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                  + z(jp,k,i) + z(jp,kp,i));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (i = is; ldir > 0 ? i <= ie : i >= ie; i += ldir) {
                    int im  = i - iinc;
                    kp  = k  + 1;
                    jp  = j  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                      + x(jp,k,i) + x(jp,kp,i));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                      + y(jp,k,i) + y(jp,kp,i));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                      + z(jp,k,i) + z(jp,kp,i));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(j,k,im,1);
                    double dnxav = 0.5*(sk(j,k,im,1) + sk(j,k+1,im,1));
                    double dnyav = 0.5*(sk(j,k,im,2) + sk(j,k+1,im,2));
                    double dnzav = 0.5*(sk(j,k,im,3) + sk(j,k+1,im,3));
                    double qthet = q(j,k,im,2)*dnxav + q(j,k,im,3)*dnyav
                                 + q(j,k,im,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qj0(k,im,5,1) = pm1 + dp;
                    qj0(k,im,5,2) = qj0(k,im,5,1);
                    r0  = r1;
                    pm1 = qj0(k,im,5,1);
                }
            }
        } // lijk==1

        if (lijk == 3) {
            // radial integration in k-direction (i is theta-direction)
            int k0, ks, ke, kk0, kinc, kkc;
            if (ldir > 0) {
                k0   = ksta;
                ks   = ksta + 1;
                ke   = kend;
                kk0  = 1;
                kinc = 1;
                kkc  = kdimc - 1;
            } else {
                k0   = kend;
                ks   = kend - 1;
                ke   = ksta;
                kk0  = kend - ksta;
                kinc = 0;
                kkc  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        bcdata(kk0,ii,ipp,2) = qj0c(kkc,i,5,1) / (double)p0;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        bcdata(kk0,ii,ipp,2) = q(j,kk0,i,5) / (double)p0;
                    }
                }
            }

            for (int i = ista; i <= iend1; i++) {
                int ii  = i - ista + 1;
                double pm1 = bcdata(kk0,ii,ipp,2) * (double)p0;
                int k   = k0;
                int ip  = i  + 1;
                int jp  = j  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,k,ip)
                                  + x(jp,k,i) + x(jp,k,ip));
                double yav0 = 0.25*(y(j,k,i)  + y(j,k,ip)
                                  + y(jp,k,i) + y(jp,k,ip));
                double zav0 = 0.25*(z(j,k,i)  + z(j,k,ip)
                                  + z(jp,k,i) + z(jp,k,ip));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (k = ks; ldir > 0 ? k <= ke : k >= ke; k += ldir) {
                    int km  = k - kinc;
                    ip  = i  + 1;
                    jp  = j  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,k,ip)
                                      + x(jp,k,i) + x(jp,k,ip));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,k,ip)
                                      + y(jp,k,i) + y(jp,k,ip));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,k,ip)
                                      + z(jp,k,i) + z(jp,k,ip));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(j,km,i,1);
                    double dnxav = 0.5*(si(j,km,i,1) + si(j,km,i+1,1));
                    double dnyav = 0.5*(si(j,km,i,2) + si(j,km,i+1,2));
                    double dnzav = 0.5*(si(j,km,i,3) + si(j,km,i+1,3));
                    double qthet = q(j,km,i,2)*dnxav + q(j,km,i,3)*dnyav
                                 + q(j,km,i,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qj0(km,i,5,1) = pm1 + dp;
                    qj0(km,i,5,2) = qj0(km,i,5,1);
                    r0  = r1;
                    pm1 = qj0(km,i,5,1);
                }
            }
        } // lijk==3

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++)
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,1) = vist3d(1,k,i);
                    vj0(k,i,1,2) = vist3d(1,k,i);
                }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++)
                    for (int i = ista; i <= iend1; i++)
                        for (int k = ksta; k <= kend1; k++) {
                            tj0(k,i,l,1) = tursav(1,k,i,l);
                            tj0(k,i,l,2) = tursav(1,k,i,l);
                        }
            }
        }
    } // nface==3



    // =========================================================
    // j=jdim boundary   bctype 2006
    // =========================================================
    if (nface == 4) {
        int j = jdim - 1;

        // extrapolate all primitives except pressure
        for (int l = 1; l <= 4; l++)
            for (int i = ista; i <= iend1; i++)
                for (int k = ksta; k <= kend1; k++) {
                    qj0(k,i,l,3) = q(jdim1,k,i,l);
                    qj0(k,i,l,4) = qj0(k,i,l,3);
                    bcj(k,i,2)   = 0.0;
                }

        if (lijk == 1) {
            // radial integration in i-direction (k is theta-direction)
            int i0, is, ie, ii0, iinc, iic;
            if (ldir > 0) {
                i0   = ista;
                is   = ista + 1;
                ie   = iend;
                ii0  = 1;
                iinc = 1;
                iic  = idimc - 1;
            } else {
                i0   = iend;
                is   = iend - 1;
                ie   = ista;
                ii0  = iend - ista;
                iinc = 0;
                iic  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        bcdata(kk,ii0,ipp,2) = qj0c(k,iic,5,3) / (double)p0;
                    }
                } else {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        bcdata(kk,ii0,ipp,2) = q(j,k,ii0,5) / (double)p0;
                    }
                }
            }

            for (int k = ksta; k <= kend1; k++) {
                int kk  = k - ksta + 1;
                double pm1 = bcdata(kk,ii0,ipp,2) * (double)p0;
                int i   = i0;
                int kp  = k  + 1;
                int jp  = j  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                  + x(jp,k,i) + x(jp,kp,i));
                double yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                  + y(jp,k,i) + y(jp,kp,i));
                double zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                  + z(jp,k,i) + z(jp,kp,i));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (i = is; ldir > 0 ? i <= ie : i >= ie; i += ldir) {
                    int im  = i - iinc;
                    kp  = k  + 1;
                    jp  = j  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                      + x(jp,k,i) + x(jp,kp,i));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                      + y(jp,k,i) + y(jp,kp,i));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                      + z(jp,k,i) + z(jp,kp,i));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(j,k,im,1);
                    double dnxav = 0.5*(sk(j,k,im,1) + sk(j,k+1,im,1));
                    double dnyav = 0.5*(sk(j,k,im,2) + sk(j,k+1,im,2));
                    double dnzav = 0.5*(sk(j,k,im,3) + sk(j,k+1,im,3));
                    double qthet = q(j,k,im,2)*dnxav + q(j,k,im,3)*dnyav
                                 + q(j,k,im,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qj0(k,im,5,3) = pm1 + dp;
                    qj0(k,im,5,4) = qj0(k,im,5,3);
                    r0  = r1;
                    pm1 = qj0(k,im,5,3);
                }
            }
        } // lijk==1

        if (lijk == 3) {
            // radial integration in k-direction (i is theta-direction)
            int k0, ks, ke, kk0, kinc, kkc;
            if (ldir > 0) {
                k0   = ksta;
                ks   = ksta + 1;
                ke   = kend;
                kk0  = 1;
                kinc = 1;
                kkc  = kdimc - 1;
            } else {
                k0   = kend;
                ks   = kend - 1;
                ke   = ksta;
                kk0  = kend - ksta;
                kinc = 0;
                kkc  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        bcdata(kk0,ii,ipp,2) = qj0c(kkc,i,5,3) / (double)p0;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        bcdata(kk0,ii,ipp,2) = q(j,kk0,i,5) / (double)p0;
                    }
                }
            }

            for (int i = ista; i <= iend1; i++) {
                int ii  = i - ista + 1;
                double pm1 = bcdata(kk0,ii,ipp,2) * (double)p0;
                int k   = k0;
                int ip  = i  + 1;
                int jp  = j  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,k,ip)
                                  + x(jp,k,i) + x(jp,k,ip));
                double yav0 = 0.25*(y(j,k,i)  + y(j,k,ip)
                                  + y(jp,k,i) + y(jp,k,ip));
                double zav0 = 0.25*(z(j,k,i)  + z(j,k,ip)
                                  + z(jp,k,i) + z(jp,k,ip));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                // NOTE: r0 not computed in Fortran for nface==4/lijk==3 outer loop
                // (Fortran do 260 loop starts with k=k0 but r0 is not set before inner loop)
                // Faithfully translate: r0 is set from xav0 etc. as in nface==3
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (k = ks; ldir > 0 ? k <= ke : k >= ke; k += ldir) {
                    int km  = k - kinc;
                    ip  = i  + 1;
                    jp  = j  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,k,ip)
                                      + x(jp,k,i) + x(jp,k,ip));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,k,ip)
                                      + y(jp,k,i) + y(jp,k,ip));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,k,ip)
                                      + z(jp,k,i) + z(jp,k,ip));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(j,km,i,1);
                    double dnxav = 0.5*(si(j,km,i,1) + si(j,km,i+1,1));
                    double dnyav = 0.5*(si(j,km,i,2) + si(j,km,i+1,2));
                    double dnzav = 0.5*(si(j,km,i,3) + si(j,km,i+1,3));
                    double qthet = q(j,km,i,2)*dnxav + q(j,km,i,3)*dnyav
                                 + q(j,km,i,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qj0(km,i,5,3) = pm1 + dp;
                    qj0(km,i,5,4) = qj0(km,i,5,3);
                    r0  = r1;
                    pm1 = qj0(km,i,5,3);
                }
            }
        } // lijk==3

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++)
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,3) = vist3d(jdim1,k,i);
                    vj0(k,i,1,4) = vist3d(jdim1,k,i);
                }
        }
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++)
                    for (int i = ista; i <= iend1; i++)
                        for (int k = ksta; k <= kend1; k++) {
                            tj0(k,i,l,3) = tursav(jdim1,k,i,l);
                            tj0(k,i,l,4) = tursav(jdim1,k,i,l);
                        }
            }
        }
    } // nface==4



    // =========================================================
    // k=1 boundary   bctype 2006
    // =========================================================
    if (nface == 5) {
        int k = 1;

        // extrapolate all primitives except pressure
        for (int l = 1; l <= 4; l++)
            for (int i = ista; i <= iend1; i++)
                for (int j = jsta; j <= jend1; j++) {
                    qk0(j,i,l,1) = q(j,1,i,l);
                    qk0(j,i,l,2) = qk0(j,i,l,1);
                    bck(j,i,1)   = 0.0;
                }

        if (lijk == 1) {
            // radial integration in i-direction (j is theta-direction)
            int i0, is, ie, ii0, iinc, iic;
            if (ldir > 0) {
                i0   = ista;
                is   = ista + 1;
                ie   = iend;
                ii0  = 1;
                iinc = 1;
                iic  = idimc - 1;
            } else {
                i0   = iend;
                is   = iend - 1;
                ie   = ista;
                ii0  = iend - ista;
                iinc = 0;
                iic  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        bcdata(jj,ii0,ipp,2) = qk0c(j,iic,5,1);
                    }
                } else {
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        bcdata(jj,ii0,ipp,2) = q(j,k,ii0,5) / (double)p0;
                    }
                }
            }

            for (int j = jsta; j <= jend1; j++) {
                int jj  = j - jsta + 1;
                double pm1 = bcdata(jj,ii0,ipp,2) * (double)p0;
                int i   = i0;
                int kp  = k  + 1;
                int jp  = j  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                  + x(jp,k,i) + x(jp,kp,i));
                double yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                  + y(jp,k,i) + y(jp,kp,i));
                double zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                  + z(jp,k,i) + z(jp,kp,i));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (i = is; ldir > 0 ? i <= ie : i >= ie; i += ldir) {
                    int im  = i - iinc;
                    kp  = k  + 1;
                    jp  = j  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                      + x(jp,k,i) + x(jp,kp,i));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                      + y(jp,k,i) + y(jp,kp,i));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                      + z(jp,k,i) + z(jp,kp,i));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(j,k,im,1);
                    double dnxav = 0.5*(sj(j,k,im,1) + sj(j+1,k,im,1));
                    double dnyav = 0.5*(sj(j,k,im,2) + sj(j+1,k,im,2));
                    double dnzav = 0.5*(sj(j,k,im,3) + sj(j+1,k,im,3));
                    double qthet = q(j,k,im,2)*dnxav + q(j,k,im,3)*dnyav
                                 + q(j,k,im,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qk0(j,im,5,1) = pm1 + dp;
                    qk0(j,im,5,2) = qk0(j,im,5,1);
                    r0  = r1;
                    pm1 = qk0(j,im,5,1);
                }
            }
        } // lijk==1

        if (lijk == 2) {
            // radial integration in j-direction (i is theta-direction)
            int j0, js, je, jj0, jinc, jjc;
            if (ldir > 0) {
                j0   = jsta;
                js   = jsta + 1;
                je   = jend;
                jj0  = 1;
                jinc = 1;
                jjc  = jdimc - 1;
            } else {
                j0   = jend;
                js   = jend - 1;
                je   = jsta;
                jj0  = jend - jsta;
                jinc = 0;
                jjc  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        bcdata(jj0,ii,ipp,2) = qk0c(jjc,i,5,1) / (double)p0;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        bcdata(jj0,ii,ipp,2) = q(jj0,k,i,5) / (double)p0;
                    }
                }
            }

            for (int i = ista; i <= iend1; i++) {
                int ii  = i - ista + 1;
                double pm1 = bcdata(jj0,ii,ipp,2) * (double)p0;
                int j   = j0;
                int ip  = i  + 1;
                int kp  = k  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                  + x(j,k,ip) + x(j,kp,ip));
                double yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                  + y(j,k,ip) + y(j,kp,ip));
                double zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                  + z(j,k,ip) + z(j,kp,ip));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (j = js; ldir > 0 ? j <= je : j >= je; j += ldir) {
                    int jm  = j - jinc;
                    kp  = k  + 1;
                    ip  = i  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                      + x(j,k,ip) + x(j,kp,ip));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                      + y(j,k,ip) + y(j,kp,ip));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                      + z(j,k,ip) + z(j,kp,ip));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(jm,k,i,1);
                    double dnxav = 0.5*(si(jm,k,i,1) + si(jm,k,i+1,1));
                    double dnyav = 0.5*(si(jm,k,i,2) + si(jm,k,i+1,2));
                    double dnzav = 0.5*(si(jm,k,i,3) + si(jm,k,i+1,3));
                    double qthet = q(jm,k,i,2)*dnxav + q(jm,k,i,3)*dnyav
                                 + q(jm,k,i,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qk0(jm,i,5,1) = pm1 + dp;
                    qk0(jm,i,5,2) = qk0(jm,i,5,1);
                    r0  = r1;
                    pm1 = qk0(jm,i,5,1);
                }
            }
        } // lijk==2

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++)
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,1) = vist3d(j,1,i);
                    vk0(j,i,1,2) = vist3d(j,1,i);
                }
        }
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++)
                    for (int i = ista; i <= iend1; i++)
                        for (int j = jsta; j <= jend1; j++) {
                            tk0(j,i,l,1) = tursav(j,1,i,l);
                            tk0(j,i,l,2) = tursav(j,1,i,l);
                        }
            }
        }
    } // nface==5



    // =========================================================
    // k=kdim boundary   bctype 2006
    // =========================================================
    if (nface == 6) {
        int k = kdim - 1;

        // extrapolate all primitives except pressure
        for (int l = 1; l <= 4; l++)
            for (int i = ista; i <= iend1; i++)
                for (int j = jsta; j <= jend1; j++) {
                    qk0(j,i,l,3) = q(j,kdim1,i,l);
                    qk0(j,i,l,4) = qk0(j,i,l,3);
                    bck(j,i,2)   = 0.0;
                }

        if (lijk == 1) {
            // radial integration in i-direction (j is theta-direction)
            int i0, is, ie, ii0, iinc, iic;
            if (ldir > 0) {
                i0   = ista;
                is   = ista + 1;
                ie   = iend;
                ii0  = 1;
                iinc = 1;
                iic  = idimc - 1;
            } else {
                i0   = iend;
                is   = iend - 1;
                ie   = ista;
                ii0  = iend - ista;
                iinc = 0;
                iic  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        bcdata(jj,ii0,ipp,2) = qk0c(j,iic,5,3) / (double)p0;
                    }
                } else {
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        bcdata(jj,ii0,ipp,2) = q(j,k,ii0,5) / (double)p0;
                    }
                }
            }

            for (int j = jsta; j <= jend1; j++) {
                int jj  = j - jsta + 1;
                double pm1 = bcdata(jj,ii0,ipp,2) * (double)p0;
                int i   = i0;
                int kp  = k  + 1;
                int jp  = j  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                  + x(jp,k,i) + x(jp,kp,i));
                double yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                  + y(jp,k,i) + y(jp,kp,i));
                double zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                  + z(jp,k,i) + z(jp,kp,i));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (i = is; ldir > 0 ? i <= ie : i >= ie; i += ldir) {
                    int im  = i - iinc;
                    kp  = k  + 1;
                    jp  = j  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                      + x(jp,k,i) + x(jp,kp,i));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                      + y(jp,k,i) + y(jp,kp,i));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                      + z(jp,k,i) + z(jp,kp,i));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(j,k,im,1);
                    double dnxav = 0.5*(sj(j,k,im,1) + sj(j+1,k,im,1));
                    double dnyav = 0.5*(sj(j,k,im,2) + sj(j+1,k,im,2));
                    double dnzav = 0.5*(sj(j,k,im,3) + sj(j+1,k,im,3));
                    double qthet = q(j,k,im,2)*dnxav + q(j,k,im,3)*dnyav
                                 + q(j,k,im,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qk0(j,im,5,3) = pm1 + dp;
                    qk0(j,im,5,4) = qk0(j,im,5,3);
                    r0  = r1;
                    pm1 = qk0(j,im,5,3);
                }
            }
        } // lijk==1

        if (lijk == 2) {
            // radial integration in j-direction (i is theta-direction)
            int j0, js, je, jj0, jinc, jjc;
            if (ldir > 0) {
                j0   = jsta;
                js   = jsta + 1;
                je   = jend;
                jj0  = 1;
                jinc = 1;
                jjc  = jdimc - 1;
            } else {
                j0   = jend;
                js   = jend - 1;
                je   = jsta;
                jj0  = jend - jsta;
                jinc = 0;
                jjc  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        bcdata(jj0,ii,ipp,2) = qk0c(jjc,i,5,3) / (double)p0;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        bcdata(jj0,ii,ipp,2) = q(jj0,k,i,5) / (double)p0;
                    }
                }
            }

            for (int i = ista; i <= iend1; i++) {
                int ii  = i - ista + 1;
                double pm1 = bcdata(jj0,ii,ipp,2) * (double)p0;
                int j   = j0;
                int ip  = i  + 1;
                int kp  = k  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                  + x(j,k,ip) + x(j,kp,ip));
                double yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                  + y(j,k,ip) + y(j,kp,ip));
                double zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                  + z(j,k,ip) + z(j,kp,ip));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (j = js; ldir > 0 ? j <= je : j >= je; j += ldir) {
                    int jm  = j - jinc;
                    kp  = k  + 1;
                    ip  = i  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                      + x(j,k,ip) + x(j,kp,ip));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                      + y(j,k,ip) + y(j,kp,ip));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                      + z(j,k,ip) + z(j,kp,ip));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(jm,k,i,1);
                    double dnxav = 0.5*(si(jm,k,i,1) + si(jm,k,i+1,1));
                    double dnyav = 0.5*(si(jm,k,i,2) + si(jm,k,i+1,2));
                    double dnzav = 0.5*(si(jm,k,i,3) + si(jm,k,i+1,3));
                    double qthet = q(jm,k,i,2)*dnxav + q(jm,k,i,3)*dnyav
                                 + q(jm,k,i,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qk0(jm,i,5,3) = pm1 + dp;
                    qk0(jm,i,5,4) = qk0(jm,i,5,3);
                    r0  = r1;
                    pm1 = qk0(jm,i,5,1);  // Fortran line 978: pm1 = qk0(jm,i,5,1) (note: slot 1, not 3)
                }
            }
        } // lijk==2

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++)
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,3) = vist3d(j,kdim1,i);
                    vk0(j,i,1,4) = vist3d(j,kdim1,i);
                }
        }
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++)
                    for (int i = ista; i <= iend1; i++)
                        for (int j = jsta; j <= jend1; j++) {
                            tk0(j,i,l,3) = tursav(j,kdim1,i,l);
                            tk0(j,i,l,4) = tursav(j,kdim1,i,l);
                        }
            }
        }
    } // nface==6



    // =========================================================
    // i=1 boundary   bctype 2006
    // =========================================================
    if (nface == 1) {
        int i = 1;

        // extrapolate all primitives except pressure
        for (int l = 1; l <= 4; l++)
            for (int k = ksta; k <= kend1; k++)
                for (int j = jsta; j <= jend1; j++) {
                    qi0(j,k,l,1) = q(j,k,1,l);
                    qi0(j,k,l,2) = qi0(j,k,l,1);
                    bci(j,k,1)   = 0.0;
                }

        if (lijk == 2) {
            // radial integration in j-direction (k is theta-direction)
            int j0, js, je, jj0, jinc, jjc;
            if (ldir > 0) {
                j0   = jsta;
                js   = jsta + 1;
                je   = jend;
                jj0  = 1;
                jinc = 1;
                jjc  = jdimc - 1;
            } else {
                j0   = jend;
                js   = jend - 1;
                je   = jsta;
                jj0  = jend - jsta;
                jinc = 0;
                jjc  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        bcdata(jj0,kk,ipp,2) = qi0c(jjc,k,5,1) / (double)p0;
                    }
                } else {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        bcdata(jj0,kk,ipp,2) = q(jj0,k,i,5) / (double)p0;
                    }
                }
            }

            for (int k = ksta; k <= kend1; k++) {
                int kk  = k - ksta + 1;
                double pm1 = bcdata(jj0,kk,ipp,2) * (double)p0;
                int j   = j0;
                int ip  = i  + 1;
                int kp  = k  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                  + x(j,k,ip) + x(j,kp,ip));
                double yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                  + y(j,k,ip) + y(j,kp,ip));
                double zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                  + z(j,k,ip) + z(j,kp,ip));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (j = js; ldir > 0 ? j <= je : j >= je; j += ldir) {
                    int jm  = j - jinc;
                    kp  = k  + 1;
                    ip  = i  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                      + x(j,k,ip) + x(j,kp,ip));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                      + y(j,k,ip) + y(j,kp,ip));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                      + z(j,k,ip) + z(j,kp,ip));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(jm,k,i,1);
                    double dnxav = 0.5*(sk(jm,k,i,1) + sk(jm,k+1,i,1));
                    double dnyav = 0.5*(sk(jm,k,i,2) + sk(jm,k+1,i,2));
                    double dnzav = 0.5*(sk(jm,k,i,3) + sk(jm,k+1,i,3));
                    double qthet = q(jm,k,i,2)*dnxav + q(jm,k,i,3)*dnyav
                                 + q(jm,k,i,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qi0(jm,k,5,1) = pm1 + dp;
                    qi0(jm,k,5,2) = qi0(jm,k,5,1);
                    r0  = r1;
                    pm1 = qi0(jm,k,5,3);  // Fortran line 1102: pm1 = qi0(jm,k,5,3)
                }
            }
        } // lijk==2

        if (lijk == 3) {
            // radial integration in k-direction (j is theta-direction)
            int k0, ks, ke, kk0, kinc, kkc;
            if (ldir > 0) {
                k0   = ksta;
                ks   = ksta + 1;
                ke   = kend;
                kk0  = 1;
                kinc = 1;
                kkc  = kdimc - 1;
            } else {
                k0   = kend;
                ks   = kend - 1;
                ke   = ksta;
                kk0  = kend - ksta;
                kinc = 0;
                kkc  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        bcdata(jj,kk0,ipp,2) = qi0c(j,kkc,5,1) / (double)p0;
                    }
                } else {
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        bcdata(jj,kk0,ipp,2) = q(j,kk0,i,5) / (double)p0;
                    }
                }
            }

            for (int j = jsta; j <= jend1; j++) {
                int jj  = j - jsta + 1;
                double pm1 = bcdata(jj,kk0,ipp,2) * (double)p0;
                int k   = k0;
                int ip  = i  + 1;
                int jp  = j  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,k,ip)
                                  + x(jp,k,i) + x(jp,k,ip));
                double yav0 = 0.25*(y(j,k,i)  + y(j,k,ip)
                                  + y(jp,k,i) + y(jp,k,ip));
                double zav0 = 0.25*(z(j,k,i)  + z(j,k,ip)
                                  + z(jp,k,i) + z(jp,k,ip));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (k = ks; ldir > 0 ? k <= ke : k >= ke; k += ldir) {
                    int km  = k - kinc;
                    ip  = i  + 1;
                    jp  = j  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,k,ip)
                                      + x(jp,k,i) + x(jp,k,ip));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,k,ip)
                                      + y(jp,k,i) + y(jp,k,ip));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,k,ip)
                                      + z(jp,k,i) + z(jp,k,ip));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(j,km,i,1);
                    double dnxav = 0.5*(sj(j,km,i,1) + sj(j+1,km,i,1));
                    double dnyav = 0.5*(sj(j,km,i,2) + sj(j+1,km,i,2));
                    double dnzav = 0.5*(sj(j,km,i,3) + sj(j+1,km,i,3));
                    double qthet = q(j,km,i,2)*dnxav + q(j,km,i,3)*dnyav
                                 + q(j,km,i,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qi0(j,km,5,1) = pm1 + dp;
                    qi0(j,km,5,2) = qi0(j,km,5,1);
                    r0  = r1;
                    pm1 = qi0(j,km,5,1);
                }
            }
        } // lijk==3

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++)
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,1) = vist3d(j,k,1);
                    vi0(j,k,1,2) = vist3d(j,k,1);
                }
        }
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++)
                    for (int k = ksta; k <= kend1; k++)
                        for (int j = jsta; j <= jend1; j++) {
                            ti0(j,k,l,1) = tursav(j,k,1,l);
                            ti0(j,k,l,2) = tursav(j,k,1,l);
                        }
            }
        }
    } // nface==1



    // =========================================================
    // i=idim boundary   bctype 2006
    // =========================================================
    if (nface == 2) {
        int i = idim - 1;

        // extrapolate all primitives except pressure
        for (int l = 1; l <= 4; l++)
            for (int k = ksta; k <= kend1; k++)
                for (int j = jsta; j <= jend1; j++) {
                    qi0(j,k,l,3) = q(j,k,idim1,l);
                    qi0(j,k,l,4) = qi0(j,k,l,3);
                    bci(j,k,2)   = 0.0;
                }

        if (lijk == 2) {
            // radial integration in j-direction (k is theta-direction)
            int j0, js, je, jj0, jinc, jjc;
            if (ldir > 0) {
                j0   = jsta;
                js   = jsta + 1;
                je   = jend;
                jj0  = 1;
                jinc = 1;
                jjc  = jdimc - 1;
            } else {
                j0   = jend;
                js   = jend - 1;
                je   = jsta;
                jj0  = jend - jsta;
                jinc = 0;
                jjc  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        bcdata(jj0,kk,ipp,2) = qi0c(jjc,k,5,3) / (double)p0;
                    }
                } else {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        bcdata(jj0,kk,ipp,2) = q(jj0,k,i,5) / (double)p0;
                    }
                }
            }

            for (int k = ksta; k <= kend1; k++) {
                int kk  = k - ksta + 1;
                double pm1 = bcdata(jj0,kk,ipp,2) * (double)p0;
                int j   = j0;
                int ip  = i  + 1;
                int kp  = k  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                  + x(j,k,ip) + x(j,kp,ip));
                double yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                  + y(j,k,ip) + y(j,kp,ip));
                double zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                  + z(j,k,ip) + z(j,kp,ip));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (j = js; ldir > 0 ? j <= je : j >= je; j += ldir) {
                    int jm  = j - jinc;
                    kp  = k  + 1;
                    ip  = i  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
                                      + x(j,k,ip) + x(j,kp,ip));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
                                      + y(j,k,ip) + y(j,kp,ip));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
                                      + z(j,k,ip) + z(j,kp,ip));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(jm,k,i,1);
                    double dnxav = 0.5*(sk(jm,k,i,1) + sk(jm,k+1,i,1));
                    double dnyav = 0.5*(sk(jm,k,i,2) + sk(jm,k+1,i,2));
                    double dnzav = 0.5*(sk(jm,k,i,3) + sk(jm,k+1,i,3));
                    double qthet = q(jm,k,i,2)*dnxav + q(jm,k,i,3)*dnyav
                                 + q(jm,k,i,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qi0(jm,k,5,3) = pm1 + dp;
                    qi0(jm,k,5,4) = qi0(jm,k,5,3);
                    r0  = r1;
                    pm1 = qi0(jm,k,5,3);
                }
            }
        } // lijk==2



        if (lijk == 3) {
            // radial integration in k-direction (j is theta-direction)
            int k0, ks, ke, kk0, kinc, kkc;
            if (ldir > 0) {
                k0   = ksta;
                ks   = ksta + 1;
                ke   = kend;
                kk0  = 1;
                kinc = 1;
                kkc  = kdimc - 1;
            } else {
                k0   = kend;
                ks   = kend - 1;
                ke   = ksta;
                kk0  = kend - ksta;
                kinc = 0;
                kkc  = 1;
            }

            if (nblc != 0) {
                if (ichk > 0) {
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        bcdata(jj,kk0,ipp,2) = qi0c(j,kkc,5,3) / (double)p0;
                    }
                } else {
                    for (int j = jsta; j <= jend1; j++) {
                        int jj = j - jsta + 1;
                        bcdata(jj,kk0,ipp,2) = q(j,kk0,i,5) / (double)p0;
                    }
                }
            }

            for (int j = jsta; j <= jend1; j++) {
                int jj  = j - jsta + 1;
                double pm1 = bcdata(jj,kk0,ipp,2) * (double)p0;
                int k   = k0;
                int ip  = i  + 1;
                int jp  = j  + 1;
                double xav0 = 0.25*(x(j,k,i)  + x(j,k,ip)
                                  + x(jp,k,i) + x(jp,k,ip));
                double yav0 = 0.25*(y(j,k,i)  + y(j,k,ip)
                                  + y(jp,k,i) + y(jp,k,ip));
                double zav0 = 0.25*(z(j,k,i)  + z(j,k,ip)
                                  + z(jp,k,i) + z(jp,k,ip));
                xav0 = xav0 * xfact;
                yav0 = yav0 * yfact;
                zav0 = zav0 * zfact;
                double r0 = std::sqrt(xav0*xav0 + yav0*yav0 + zav0*zav0);

                for (k = ks; ldir > 0 ? k <= ke : k >= ke; k += ldir) {
                    int km  = k - kinc;
                    ip  = i  + 1;
                    jp  = j  + 1;
                    double xav1 = 0.25*(x(j,k,i)  + x(j,k,ip)
                                      + x(jp,k,i) + x(jp,k,ip));
                    double yav1 = 0.25*(y(j,k,i)  + y(j,k,ip)
                                      + y(jp,k,i) + y(jp,k,ip));
                    double zav1 = 0.25*(z(j,k,i)  + z(j,k,ip)
                                      + z(jp,k,i) + z(jp,k,ip));
                    xav1 = xav1 * xfact;
                    yav1 = yav1 * yfact;
                    zav1 = zav1 * zfact;
                    double r1    = std::sqrt(xav1*xav1 + yav1*yav1 + zav1*zav1);
                    double rho   = q(j,km,i,1);
                    double dnxav = 0.5*(sj(j,km,i,1) + sj(j+1,km,i,1));
                    double dnyav = 0.5*(sj(j,km,i,2) + sj(j+1,km,i,2));
                    double dnzav = 0.5*(sj(j,km,i,3) + sj(j+1,km,i,3));
                    double qthet = q(j,km,i,2)*dnxav + q(j,km,i,3)*dnyav
                                 + q(j,km,i,4)*dnzav;
                    double rav   = 0.5*(r1 + r0);
                    double dr    = r1 - r0;
                    double dp    = dr * rho * qthet*qthet / rav;
                    qi0(j,km,5,3) = pm1 + dp;
                    qi0(j,km,5,4) = qi0(j,km,5,3);
                    r0  = r1;
                    pm1 = qi0(j,km,5,3);
                }
            }
        } // lijk==3

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++)
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,3) = vist3d(j,k,idim1);
                    vi0(j,k,1,4) = vist3d(j,k,idim1);
                }
        }
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++)
                    for (int k = ksta; k <= kend1; k++)
                        for (int j = jsta; j <= jend1; j++) {
                            ti0(j,k,l,3) = tursav(j,k,idim1,l);
                            ti0(j,k,l,4) = tursav(j,k,idim1,l);
                        }
            }
        }
    } // nface==2

    return;
} // bc2006

} // namespace bc2006_ns
