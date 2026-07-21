// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2004.h"
#include <cstdio>
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>

namespace bc2004_ns {

// -----------------------------------------------------------------------
// bc() dispatcher
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
    int iuns, irelv;
    char filname[80];

    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);
        int jdim1 = jdim - 1;
        int kdim1 = kdim - 1;
        int idim1 = idim - 1;

        // Build array views from w using lw offsets
        FortranArray4DRef<double> q   (&w(lw(1,ibl)), jdim,  kdim,  idim,  5);
        FortranArray4DRef<double> qj0 (&w(lw(2,ibl)), kdim,  idim1, 5,     4);
        FortranArray4DRef<double> qk0 (&w(lw(3,ibl)), jdim,  idim1, 5,     4);
        FortranArray4DRef<double> qi0 (&w(lw(4,ibl)), jdim,  kdim,  5,     4);
        FortranArray4DRef<double> sj  (&w(lw(5,ibl)), jdim,  kdim,  idim1, 5);
        FortranArray4DRef<double> sk  (&w(lw(6,ibl)), jdim,  kdim,  idim1, 5);
        FortranArray4DRef<double> si  (&w(lw(7,ibl)), jdim,  kdim,  idim,  5);
        FortranArray3DRef<double> bcj (&w(lw(8,ibl)), kdim,  idim1, 2);
        FortranArray3DRef<double> bck (&w(lw(9,ibl)), jdim,  idim1, 2);
        FortranArray3DRef<double> bci (&w(lw(10,ibl)),jdim,  kdim,  2);
        FortranArray4DRef<double> xtbj(&w(lw(11,ibl)),kdim,  idim1, 3,     2);
        FortranArray4DRef<double> xtbk(&w(lw(12,ibl)),jdim,  idim1, 3,     2);
        FortranArray4DRef<double> xtbi(&w(lw(13,ibl)),jdim,  kdim,  3,     2);
        FortranArray4DRef<double> atbj(&w(lw(14,ibl)),kdim,  idim1, 3,     2);
        FortranArray4DRef<double> atbk(&w(lw(15,ibl)),jdim,  idim1, 3,     2);
        FortranArray4DRef<double> atbi(&w(lw(16,ibl)),jdim,  kdim,  3,     2);

        // Build array views from wk using lw2 offsets
        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,  kdim,  idim,  nummem);
        FortranArray4DRef<double> tj0   (&wk(lw2(2,ibl)), kdim,  idim1, nummem,4);
        FortranArray4DRef<double> tk0   (&wk(lw2(3,ibl)), jdim,  idim1, nummem,4);
        FortranArray4DRef<double> ti0   (&wk(lw2(4,ibl)), jdim,  kdim,  nummem,4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,  kdim,  idim);
        FortranArray4DRef<double> vj0   (&wk(lw2(6,ibl)), kdim,  idim1, 1,     4);
        FortranArray4DRef<double> vk0   (&wk(lw2(7,ibl)), jdim,  idim1, 1,     4);
        FortranArray4DRef<double> vi0   (&wk(lw2(8,ibl)), jdim,  kdim,  1,     4);
        FortranArray3DRef<double> smin  (&wk(lw2(9,ibl)), jdim1, kdim1, idim1);
        FortranArray3DRef<double> snj0  (&wk(lw2(10,ibl)),jdim1, kdim1, idim1);
        FortranArray3DRef<double> snk0  (&wk(lw2(11,ibl)),jdim1, kdim1, idim1);
        FortranArray3DRef<double> sni0  (&wk(lw2(12,ibl)),jdim1, kdim1, idim1);
        FortranArray3DRef<double> snjm  (&wk(lw2(13,ibl)),jdim1, kdim1, idim1);
        FortranArray3DRef<double> snkm  (&wk(lw2(14,ibl)),jdim1, kdim1, idim1);
        FortranArray3DRef<double> snim  (&wk(lw2(15,ibl)),jdim1, kdim1, idim1);

        iuns  = irotat(ibl);
        irelv = idefrm(ibl);

        // Loop over j-face segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            int bctype = jbcinfo(ibl,iseg,1,1);
            int iuse3  = 0;
            if      (bctype == 2004) iuse3 = 0;
            else if (bctype == 2014) iuse3 = 0;
            else if (bctype == 2024) iuse3 = 2;
            else if (bctype == 2034) iuse3 = 1;
            else continue;

            int ista  = jbcinfo(ibl,iseg,1,2);
            int iend  = jbcinfo(ibl,iseg,1,3);
            int jsta  = jbcinfo(ibl,iseg,1,4);
            int jend  = jbcinfo(ibl,iseg,1,5);
            int ksta  = jbcinfo(ibl,iseg,1,6);
            int kend  = jbcinfo(ibl,iseg,1,7);
            int nface = jbcinfo(ibl,iseg,1,8);
            int mdim  = lwdat(ibl,iseg,3);
            int ndim  = lwdat(ibl,iseg,4);
            FortranArray4DRef<double> bcdata(&wk(lwdat(ibl,iseg,1)), mdim, ndim, 2, 12);
            int ifilj = bcfilej(ibl,iseg,1);
            std::memcpy(filname, bcfiles(ifilj), 80);

            bc2004(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, smin, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns, irelv,
                   snj0, snk0, sni0, ntime, snjm, snkm, snim,
                   nou, bou, nbuf, ibufdim, myid, nummem, iuse3);
        }

        // Loop over k-face segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            int bctype = kbcinfo(ibl,iseg,1,1);
            int iuse3  = 0;
            if      (bctype == 2004) iuse3 = 0;
            else if (bctype == 2014) iuse3 = 0;
            else if (bctype == 2024) iuse3 = 2;
            else if (bctype == 2034) iuse3 = 1;
            else continue;

            int ista  = kbcinfo(ibl,iseg,1,2);
            int iend  = kbcinfo(ibl,iseg,1,3);
            int jsta  = kbcinfo(ibl,iseg,1,4);
            int jend  = kbcinfo(ibl,iseg,1,5);
            int ksta  = kbcinfo(ibl,iseg,1,6);
            int kend  = kbcinfo(ibl,iseg,1,7);
            int nface = kbcinfo(ibl,iseg,1,8);
            int mdim  = lwdat(ibl,iseg,3);
            int ndim  = lwdat(ibl,iseg,4);
            FortranArray4DRef<double> bcdata(&wk(lwdat(ibl,iseg,1)), mdim, ndim, 2, 12);
            int ifilk = bcfilek(ibl,iseg,1);
            std::memcpy(filname, bcfiles(ifilk), 80);

            bc2004(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, smin, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns, irelv,
                   snj0, snk0, sni0, ntime, snjm, snkm, snim,
                   nou, bou, nbuf, ibufdim, myid, nummem, iuse3);
        }

        // Loop over i-face segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            int bctype = ibcinfo(ibl,iseg,1,1);
            int iuse3  = 0;
            if      (bctype == 2004) iuse3 = 0;
            else if (bctype == 2014) iuse3 = 0;
            else if (bctype == 2024) iuse3 = 2;
            else if (bctype == 2034) iuse3 = 1;
            else continue;

            int ista  = ibcinfo(ibl,iseg,1,2);
            int iend  = ibcinfo(ibl,iseg,1,3);
            int jsta  = ibcinfo(ibl,iseg,1,4);
            int jend  = ibcinfo(ibl,iseg,1,5);
            int ksta  = ibcinfo(ibl,iseg,1,6);
            int kend  = ibcinfo(ibl,iseg,1,7);
            int nface = ibcinfo(ibl,iseg,1,8);
            int mdim  = lwdat(ibl,iseg,3);
            int ndim  = lwdat(ibl,iseg,4);
            FortranArray4DRef<double> bcdata(&wk(lwdat(ibl,iseg,1)), mdim, ndim, 2, 12);
            int ifili = bcfilei(ibl,iseg,1);
            std::memcpy(filname, bcfiles(ifili), 80);

            bc2004(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                   bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                   ista, iend, jsta, jend, ksta, kend, nface,
                   tursav, tj0, tk0, ti0, smin, vist3d, vj0, vk0, vi0,
                   mdim, ndim, bcdata, filname, iuns, irelv,
                   snj0, snk0, sni0, ntime, snjm, snkm, snim,
                   nou, bou, nbuf, ibufdim, myid, nummem, iuse3);
        }
    }
}

// -----------------------------------------------------------------------
// bc2004() - solid wall viscous BC
// -----------------------------------------------------------------------
void bc2004(int& jdim, int& kdim, int& idim,
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
            FortranArray3DRef<double> smin, FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& mdim, int& ndim,
            FortranArray4DRef<double> bcdata,
            char (&filname)[80],
            int& iuns, int& irelv,
            FortranArray3DRef<double> snj0, FortranArray3DRef<double> snk0,
            FortranArray3DRef<double> sni0,
            int& ntime,
            FortranArray3DRef<double> snjm, FortranArray3DRef<double> snkm,
            FortranArray3DRef<double> snim,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid, int& nummem, int& iuse3)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    float& gm1     = cmn_fluid.gm1;
    float& gp1     = cmn_fluid.gp1;
    float& gm1g    = cmn_fluid.gm1g;
    float& gp1g    = cmn_fluid.gp1g;
    float& ggm1    = cmn_fluid.ggm1;
    float& pr      = cmn_fluid2.pr;
    float& prt     = cmn_fluid2.prt;
    float& cbar    = cmn_fluid2.cbar;
    float& xmach   = cmn_info.xmach;
    int32_t& ivmx  = cmn_maxiv.ivmx;
    int32_t& level = cmn_mgrd.level;
    int32_t& lglobal = cmn_mgrd.lglobal;
    float& reue    = cmn_reyue.reue;
    float& tinf    = cmn_reyue.tinf;
    int32_t* ivisc = cmn_reyue.ivisc;
    int32_t* iwf   = cmn_wallfun.iwf;
    int32_t& issglrrw2012 = cmn_reystressmodel.issglrrw2012;
    float& uub_2034 = cmn_wallvel_2034.uub_2034;
    float& vvb_2034 = cmn_wallvel_2034.vvb_2034;
    float& wwb_2034 = cmn_wallvel_2034.wwb_2034;
    float& vvb_xaxisrot2034 = cmn_wallvel_2034.vvb_xaxisrot2034;

    // Data statement arrays (wall function polynomial coefficients)
    static const double a0[7] = {2.354039, 0.1179840, -4.2899192e-04, 2.0404148e-06,
                                  -5.1775775e-09, 6.2687308e-12, -2.916958e-15};
    static const double a1[5] = {5.777191, 6.8756983e-02, -7.1582745e-06, 1.5594904e-09,
                                  -1.4865778e-13};
    static const double a2[3] = {31.08654, 5.0429072e-02, -2.0072314e-8};

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // this bc makes use of only one plane of data
    int ip = 1;


    // **************************************************************************
    // j=1 boundary   viscous wall with T & cq specified   bctype 2004
    // **************************************************************************
    if (nface == 3) {

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            for (int k = ksta; k <= kend1; k++) {
                int kk = k - ksta + 1;
                double cq  = bcdata(kk,ii,ip,2);
                double pb  = q(1,k,i,5);
                double dpb = q(2,k,i,5) - q(1,k,i,5);
                pb = pb - dpb/2.0;
                if (pb <= 0.0) pb = q(1,k,i,5);
                double c2 = (double)gamma * q(1,k,i,5) / q(1,k,i,1);
                if (bcdata(kk,ii,ip,1) > 0.) {
                    c2 = bcdata(kk,ii,ip,1);
                } else if (bcdata(kk,ii,ip,1) < 0.) {
                    c2 = 1.e0 + (double)gm1 * 0.5e0 * (double)xmach * (double)xmach;
                } else {
                    double xm2;
                    if (iuns > 0 && irelv > 0) {
                        xm2 = (q(1,k,i,2)-xtbj(k,i,1,1))*(q(1,k,i,2)-xtbj(k,i,1,1))
                            + (q(1,k,i,3)-xtbj(k,i,2,1))*(q(1,k,i,3)-xtbj(k,i,2,1))
                            + (q(1,k,i,4)-xtbj(k,i,3,1))*(q(1,k,i,4)-xtbj(k,i,3,1));
                    } else {
                        xm2 = q(1,k,i,2)*q(1,k,i,2)
                            + q(1,k,i,3)*q(1,k,i,3)
                            + q(1,k,i,4)*q(1,k,i,4);
                    }
                    xm2 = xm2 / c2;
                    c2  = c2 * (1. + 0.5*(double)gm1*xm2);
                }
                double uub, vvb, wwb;
                if (iuse3 == 1) {
                    uub = (double)uub_2034;
                    vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034 * sj(1,k,i,3);
                    wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034 * sj(1,k,i,2);
                } else {
                    uub = 0.; vvb = 0.; wwb = 0.;
                }
                if (iuns > 0 && irelv > 0) {
                    if (iuse3 == 1) {
                        uub = (double)uub_2034 + xtbj(k,i,1,1);
                        vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034*sj(1,k,i,3) + xtbj(k,i,2,1);
                        wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034*sj(1,k,i,2) + xtbj(k,i,3,1);
                    } else {
                        uub = xtbj(k,i,1,1); vvb = xtbj(k,i,2,1); wwb = xtbj(k,i,3,1);
                    }
                }
                uub = uub + (double)xmach*cq*sj(1,k,i,1)*c2/((double)gamma*pb);
                vvb = vvb + (double)xmach*cq*sj(1,k,i,2)*c2/((double)gamma*pb);
                wwb = wwb + (double)xmach*cq*sj(1,k,i,3)*c2/((double)gamma*pb);
                qj0(k,i,1,1) = (double)gamma*pb/c2;
                qj0(k,i,2,1) = uub; qj0(k,i,3,1) = vvb; qj0(k,i,4,1) = wwb;
                qj0(k,i,5,1) = pb;
                bcj(k,i,1)   = 1.0;
                double f23 = 0.0;
                int j2 = std::min(2, jdim1);
                if (j2 == 1) f23 = 0.0;
                double z1 =   2.0 + 1.5*f23;
                double z2 =       - 0.5*f23;
                double z3 = -(2.0 +     f23);
                qj0(k,i,1,2) = z1*q(1,k,i,1) + z2*q(j2,k,i,1) + z3*qj0(k,i,1,1);
                qj0(k,i,2,2) = z1*q(1,k,i,2) + z2*q(j2,k,i,2) + z3*qj0(k,i,2,1);
                qj0(k,i,3,2) = z1*q(1,k,i,3) + z2*q(j2,k,i,3) + z3*qj0(k,i,3,1);
                qj0(k,i,4,2) = z1*q(1,k,i,4) + z2*q(j2,k,i,4) + z3*qj0(k,i,4,1);
                qj0(k,i,5,2) = z1*q(1,k,i,5) + z2*q(j2,k,i,5) + z3*qj0(k,i,5,1);
            }
        }

        if (ivmx >= 2) {
            if (level >= lglobal && ntime != 0) {
                if (iwf[1] == 0) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        vj0(k,i,1,1) = 0.; vj0(k,i,1,2) = 0.;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        int j = 1;
                        double c2b  = (double)cbar / (double)tinf;
                        double c2bp = c2b + 1.0;
                        double uu = std::sqrt(
                            (q(j,k,i,2)-qj0(k,i,2,1))*(q(j,k,i,2)-qj0(k,i,2,1)) +
                            (q(j,k,i,3)-qj0(k,i,3,1))*(q(j,k,i,3)-qj0(k,i,3,1)) +
                            (q(j,k,i,4)-qj0(k,i,4,1))*(q(j,k,i,4)-qj0(k,i,4,1)) );
                        double dist = (ivmx == 2) ? snj0(j,k,i) : std::abs(smin(j,k,i));
                        double tt   = (double)gamma*qj0(k,i,5,1)/qj0(k,i,1,1);
                        double fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                        double xnplus;
                        if (rc <= 20.24) {
                            xnplus = std::sqrt(rc);
                        } else if (rc <= 435.) {
                            xnplus = a0[0]+a0[1]*rc+(a0[2]*rc)*rc+(a0[3]*rc)*rc*rc
                                   +(a0[4]*rc)*rc*rc*rc+(a0[5]*rc)*rc*rc*rc*rc
                                   +(a0[6]*rc)*rc*rc*rc*rc*rc;
                        } else if (rc <= 4000.) {
                            xnplus = a1[0]+a1[1]*rc+(a1[2]*rc)*rc+(a1[3]*rc)*rc*rc
                                   +(a1[4]*rc)*rc*rc*rc;
                        } else {
                            xnplus = a2[0]+a2[1]*rc+a2[2]*rc*rc;
                        }
                        double xnplussav = xnplus;
                        if (xnplus > 10.) {
                            bool conv = false;
                            for (int num = 1; num <= 10; num++) {
                                double f_val = rc/xnplus - 2.44*std::log(xnplus) - 5.2;
                                double dfdn  = -rc/(xnplus*xnplus) - 2.44/xnplus;
                                double delta = -f_val/dfdn;
                                xnplus = std::abs(xnplus + delta);
                                if (std::abs(delta) < 1.e-3) { conv = true; break; }
                            }
                            if (!conv) xnplus = xnplussav;
                        }
                        double dudy = uu / dist;
                        double xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/
                                      (q(j,k,i,1)*dudy)*(double)xmach/(double)reue - 1.);
                        vj0(k,i,1,1) = xmut; vj0(k,i,1,2) = 0.;
                    }
                }
            }
        }

        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal && ntime != 0) {
            if (ivmx == 4 || ivmx == 5) {
                for (int i = ista; i <= iend1; i++)
                for (int k = ksta; k <= kend1; k++) {
                    tj0(k,i,1,1) = -tursav(1,k,i,1);
                    tj0(k,i,2,1) = -tursav(1,k,i,2);
                    tj0(k,i,1,2) = 2.*tj0(k,i,1,1) - tursav(1,k,i,1);
                    tj0(k,i,2,2) = 2.*tj0(k,i,2,1) - tursav(1,k,i,2);
                }
            }
            if (ivmx >= 6) {
                double c2b  = (double)cbar / (double)tinf;
                double c2bp = c2b + 1.0;
                double re   = (double)reue / (double)xmach;
                double beta1;
                if      (ivmx == 6)                              beta1 = .075;
                else if (ivmx == 7 || ivmx == 30 || ivmx == 40) beta1 = .075;
                else if (ivmx == 8 || ivmx == 12 || ivmx == 14) beta1 = .83;
                else if (ivmx == 72)                             beta1 = .0708;
                else                                             beta1 = .075;
                int j = 1;
                if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        double tt   = (double)gamma*qj0(k,i,5,1)/qj0(k,i,1,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                        tj0(k,i,1,1) = 2.*(2.*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)
                                     - tursav(j,k,i,1);
                        tj0(k,i,2,1) = -tursav(j,k,i,2);
                        tj0(k,i,1,2) = 2.*tj0(k,i,1,1) - tursav(1,k,i,1);
                        tj0(k,i,2,2) = 2.*tj0(k,i,2,1) - tursav(1,k,i,2);
                    }
                } else if (ivmx == 15) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = 0.5*(tursav(j,k,i,2)+tursav(j+1,k,i,2))/3.;
                        tj0(k,i,1,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,1);
                        tj0(k,i,2,1) = -tursav(j,k,i,2);
                        tj0(k,i,1,2) = 0.; tj0(k,i,2,2) = 0.;
                    }
                } else if (ivmx == 16) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        tj0(k,i,1,1) = -tursav(j,k,i,1);
                        tj0(k,i,2,1) = -tursav(j,k,i,2);
                        tj0(k,i,1,2) = 0.; tj0(k,i,2,2) = 0.;
                    }
                } else if (ivmx == 70) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        double tke1, tke2;
                        for (int iv = 1; iv <= 6; iv++) {
                            tj0(k,i,iv,1) = -tursav(j,k,i,iv);
                            tj0(k,i,iv,2) = std::numeric_limits<double>::max();
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        tke1 = -0.5*(tursav(j,k,i,1)+tursav(j,k,i,2)+tursav(j,k,i,3));
                        tke2 = -0.5*(tursav(j+1,k,i,1)+tursav(j+1,k,i,2)+tursav(j+1,k,i,3));
                        double dkdy = 0.5*(tke1+tke2)/3.*0.5;
                        tj0(k,i,7,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,7);
                        tj0(k,i,7,2) = 0.0;
                    }
                } else if (ivmx == 72 && issglrrw2012 != 6) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        for (int iv = 1; iv <= 6; iv++) {
                            tj0(k,i,iv,1) = -tursav(j,k,i,iv);
                            tj0(k,i,iv,2) = 2*tj0(k,i,iv,1) - tursav(j,k,i,iv);
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        tj0(k,i,7,1) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,7);
                        tj0(k,i,7,2) = 2*tj0(k,i,7,1) - tursav(j,k,i,7);
                    }
                } else if (ivmx == 72 && issglrrw2012 == 6) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        for (int iv = 1; iv <= 7; iv++) {
                            tj0(k,i,iv,1) = -tursav(j,k,i,iv);
                            tj0(k,i,iv,2) = 2*tj0(k,i,iv,1) - tursav(j,k,i,iv);
                        }
                    }
                } else if (ivmx == 25) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        tj0(k,i,1,1) = 0.; tj0(k,i,2,1) = 0.;
                        tj0(k,i,1,2) = 0.; tj0(k,i,2,2) = 0.;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        double tt   = (double)gamma*qj0(k,i,5,1)/qj0(k,i,1,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        tj0(k,i,1,1) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,1);
                        tj0(k,i,2,1) = -tursav(j,k,i,2);
                        tj0(k,i,1,2) = 2.*tj0(k,i,1,1) - tursav(1,k,i,1);
                        tj0(k,i,2,2) = 2.*tj0(k,i,2,1) - tursav(1,k,i,2);
                    }
                }
                if (ivmx == 30) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            int kk = k - ksta + 1;
                            if (iuse3 != 2) {
                                tj0(k,i,3,1) = tursav(j,k,i,3);
                            } else {
                                double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                                double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                                double dist = std::abs(smin(j,k,i));
                                double hee  = (1.+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                                double rhoee= std::pow((double)gamma*q(j,k,i,5), 1./(double)gamma);
                                double uee2 = ((double)gm1*rhoee*hee-(double)gamma*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                                double uee  = (uee2 > 0.) ? std::sqrt(uee2) : (double)xmach;
                                tj0(k,i,3,1) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
                                               uee*bcdata(kk,ii,ip,3)*(double)reue/(double)xmach;
                            }
                            tj0(k,i,3,2) = 2.*tj0(k,i,3,1) - tursav(1,k,i,3);
                        }
                    }
                } else if (ivmx == 40) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            int kk = k - ksta + 1;
                            tj0(k,i,3,1) = tursav(j,k,i,3);
                            tj0(k,i,4,1) = tursav(j,k,i,4);
                        }
                    }
                }
            }
        }

    } // end nface==3


    // **************************************************************************
    // j=jdim boundary   viscous wall with T & cq specified   bctype 2004
    // **************************************************************************
    if (nface == 4) {

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            for (int k = ksta; k <= kend1; k++) {
                int kk = k - ksta + 1;
                double cq  = bcdata(kk,ii,ip,2);
                double pb  = q(jdim1,k,i,5);
                double dpb = q(jdim1,k,i,5) - q(jdim1-1,k,i,5);
                pb = pb + dpb/2.0;
                if (pb <= 0.0) pb = q(jdim1,k,i,5);
                double c2 = (double)gamma * q(jdim1,k,i,5) / q(jdim1,k,i,1);
                if (bcdata(kk,ii,ip,1) > 0.) {
                    c2 = bcdata(kk,ii,ip,1);
                } else if (bcdata(kk,ii,ip,1) < 0.) {
                    c2 = 1.e0 + (double)gm1 * 0.5e0 * (double)xmach * (double)xmach;
                } else {
                    double xm2;
                    if (iuns > 0 && irelv > 0) {
                        xm2 = (q(jdim1,k,i,2)-xtbj(k,i,1,2))*(q(jdim1,k,i,2)-xtbj(k,i,1,2))
                            + (q(jdim1,k,i,3)-xtbj(k,i,2,2))*(q(jdim1,k,i,3)-xtbj(k,i,2,2))
                            + (q(jdim1,k,i,4)-xtbj(k,i,3,2))*(q(jdim1,k,i,4)-xtbj(k,i,3,2));
                    } else {
                        xm2 = q(jdim1,k,i,2)*q(jdim1,k,i,2)
                            + q(jdim1,k,i,3)*q(jdim1,k,i,3)
                            + q(jdim1,k,i,4)*q(jdim1,k,i,4);
                    }
                    xm2 = xm2 / c2;
                    c2  = c2 * (1. + 0.5*(double)gm1*xm2);
                }
                double uub, vvb, wwb;
                if (iuse3 == 1) {
                    uub = (double)uub_2034;
                    vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034 * sj(jdim,k,i,3);
                    wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034 * sj(jdim,k,i,2);
                } else {
                    uub = 0.; vvb = 0.; wwb = 0.;
                }
                if (iuns > 0 && irelv > 0) {
                    if (iuse3 == 1) {
                        uub = (double)uub_2034 + xtbj(k,i,1,2);
                        vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034*sj(jdim,k,i,3) + xtbj(k,i,2,2);
                        wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034*sj(jdim,k,i,2) + xtbj(k,i,3,2);
                    } else {
                        uub = xtbj(k,i,1,2); vvb = xtbj(k,i,2,2); wwb = xtbj(k,i,3,2);
                    }
                }
                uub = uub - (double)xmach*cq*sj(jdim,k,i,1)*c2/((double)gamma*pb);
                vvb = vvb - (double)xmach*cq*sj(jdim,k,i,2)*c2/((double)gamma*pb);
                wwb = wwb - (double)xmach*cq*sj(jdim,k,i,3)*c2/((double)gamma*pb);
                qj0(k,i,1,3) = (double)gamma*pb/c2;
                qj0(k,i,2,3) = uub; qj0(k,i,3,3) = vvb; qj0(k,i,4,3) = wwb;
                qj0(k,i,5,3) = pb;
                bcj(k,i,2)   = 1.0;
                double f23 = 0.0;
                int j2 = std::max(1, jdim-2);
                if (j2 == 1) f23 = 0.0;
                double z1 =  -2.0 - 1.5*f23;
                double z2 =       + 0.5*f23;
                double z3 = +(2.0 +     f23);
                qj0(k,i,1,4) = z1*q(jdim1,k,i,1)+z2*q(j2,k,i,1)+z3*qj0(k,i,1,3);
                qj0(k,i,2,4) = z1*q(jdim1,k,i,2)+z2*q(j2,k,i,2)+z3*qj0(k,i,2,3);
                qj0(k,i,3,4) = z1*q(jdim1,k,i,3)+z2*q(j2,k,i,3)+z3*qj0(k,i,3,3);
                qj0(k,i,4,4) = z1*q(jdim1,k,i,4)+z2*q(j2,k,i,4)+z3*qj0(k,i,4,3);
                qj0(k,i,5,4) = z1*q(jdim1,k,i,5)+z2*q(j2,k,i,5)+z3*qj0(k,i,5,3);
            }
        }

        if (ivmx >= 2) {
            if (level >= lglobal && ntime != 0) {
                if (iwf[1] == 0) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        vj0(k,i,1,3) = 0.; vj0(k,i,1,4) = 0.;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        int j = jdim - 1;
                        double c2b  = (double)cbar / (double)tinf;
                        double c2bp = c2b + 1.0;
                        double uu = std::sqrt(
                            (q(j,k,i,2)-qj0(k,i,2,3))*(q(j,k,i,2)-qj0(k,i,2,3)) +
                            (q(j,k,i,3)-qj0(k,i,3,3))*(q(j,k,i,3)-qj0(k,i,3,3)) +
                            (q(j,k,i,4)-qj0(k,i,4,3))*(q(j,k,i,4)-qj0(k,i,4,3)) );
                        double dist = (ivmx == 2) ? snjm(j,k,i) : std::abs(smin(j,k,i));
                        double tt   = (double)gamma*qj0(k,i,5,3)/qj0(k,i,1,3);
                        double fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                        double xnplus;
                        if (rc <= 20.24) {
                            xnplus = std::sqrt(rc);
                        } else if (rc <= 435.) {
                            xnplus = a0[0]+a0[1]*rc+(a0[2]*rc)*rc+(a0[3]*rc)*rc*rc
                                   +(a0[4]*rc)*rc*rc*rc+(a0[5]*rc)*rc*rc*rc*rc
                                   +(a0[6]*rc)*rc*rc*rc*rc*rc;
                        } else if (rc <= 4000.) {
                            xnplus = a1[0]+a1[1]*rc+(a1[2]*rc)*rc+(a1[3]*rc)*rc*rc
                                   +(a1[4]*rc)*rc*rc*rc;
                        } else {
                            xnplus = a2[0]+a2[1]*rc+a2[2]*rc*rc;
                        }
                        double xnplussav = xnplus;
                        if (xnplus > 10.) {
                            bool conv = false;
                            for (int num = 1; num <= 10; num++) {
                                double f_val = rc/xnplus - 2.44*std::log(xnplus) - 5.2;
                                double dfdn  = -rc/(xnplus*xnplus) - 2.44/xnplus;
                                double delta = -f_val/dfdn;
                                xnplus = std::abs(xnplus + delta);
                                if (std::abs(delta) < 1.e-3) { conv = true; break; }
                            }
                            if (!conv) xnplus = xnplussav;
                        }
                        double dudy = uu / dist;
                        double xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/
                                      (q(j,k,i,1)*dudy)*(double)xmach/(double)reue - 1.);
                        vj0(k,i,1,3) = xmut; vj0(k,i,1,4) = 0.;
                    }
                }
            }
        }

        if (level >= lglobal && ntime != 0) {
            if (ivmx == 4 || ivmx == 5) {
                for (int i = ista; i <= iend1; i++)
                for (int k = ksta; k <= kend1; k++) {
                    tj0(k,i,1,3) = -tursav(jdim-1,k,i,1);
                    tj0(k,i,2,3) = -tursav(jdim-1,k,i,2);
                    tj0(k,i,1,4) = 2.*tj0(k,i,1,3) - tursav(jdim-1,k,i,1);
                    tj0(k,i,2,4) = 2.*tj0(k,i,2,3) - tursav(jdim-1,k,i,2);
                }
            }
            if (ivmx >= 6) {
                double c2b  = (double)cbar / (double)tinf;
                double c2bp = c2b + 1.0;
                double re   = (double)reue / (double)xmach;
                double beta1;
                if      (ivmx == 6)                              beta1 = .075;
                else if (ivmx == 7 || ivmx == 30 || ivmx == 40) beta1 = .075;
                else if (ivmx == 8 || ivmx == 12 || ivmx == 14) beta1 = .83;
                else if (ivmx == 72)                             beta1 = .0708;
                else                                             beta1 = .075;
                int j = jdim - 1;
                if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        double tt   = (double)gamma*qj0(k,i,5,3)/qj0(k,i,1,3);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                        tj0(k,i,1,3) = 2.*(2.*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)
                                     - tursav(j,k,i,1);
                        tj0(k,i,2,3) = -tursav(j,k,i,2);
                        tj0(k,i,1,4) = 2.*tj0(k,i,1,3) - tursav(jdim-1,k,i,1);
                        tj0(k,i,2,4) = 2.*tj0(k,i,2,3) - tursav(jdim-1,k,i,2);
                    }
                } else if (ivmx == 15) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = 0.5*(tursav(j,k,i,2)+tursav(j-1,k,i,2))/3.;
                        tj0(k,i,1,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,1);
                        tj0(k,i,2,3) = -tursav(j,k,i,2);
                        tj0(k,i,1,4) = 0.; tj0(k,i,2,4) = 0.;
                    }
                } else if (ivmx == 16) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        tj0(k,i,1,3) = -tursav(j,k,i,1);
                        tj0(k,i,2,3) = -tursav(j,k,i,2);
                        tj0(k,i,1,4) = 0.; tj0(k,i,2,4) = 0.;
                    }
                } else if (ivmx == 70) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        double tke1, tke2;
                        for (int iv = 1; iv <= 6; iv++) {
                            tj0(k,i,iv,3) = -tursav(j,k,i,iv);
                            tj0(k,i,iv,4) = std::numeric_limits<double>::max();
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        tke1 = -0.5*(tursav(j,k,i,1)+tursav(j,k,i,2)+tursav(j,k,i,3));
                        tke2 = -0.5*(tursav(j-1,k,i,1)+tursav(j-1,k,i,2)+tursav(j-1,k,i,3));
                        double dkdy = 0.5*(tke1+tke2)/3.*0.5;
                        tj0(k,i,7,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,7);
                        tj0(k,i,7,4) = 0.0;
                    }
                } else if (ivmx == 72 && issglrrw2012 != 6) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        for (int iv = 1; iv <= 6; iv++) {
                            tj0(k,i,iv,3) = -tursav(j,k,i,iv);
                            tj0(k,i,iv,4) = 2*tj0(k,i,iv,3) - tursav(j,k,i,iv);
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        tj0(k,i,7,3) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,7);
                        tj0(k,i,7,4) = 2*tj0(k,i,7,3) - tursav(j,k,i,7);
                    }
                } else if (ivmx == 72 && issglrrw2012 == 6) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        for (int iv = 1; iv <= 7; iv++) {
                            tj0(k,i,iv,3) = -tursav(j,k,i,iv);
                            tj0(k,i,iv,4) = 2*tj0(k,i,iv,3) - tursav(j,k,i,iv);
                        }
                    }
                } else if (ivmx == 25) {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        tj0(k,i,1,3) = 0.; tj0(k,i,2,3) = 0.;
                        tj0(k,i,1,4) = 0.; tj0(k,i,2,4) = 0.;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++)
                    for (int k = ksta; k <= kend1; k++) {
                        double tt   = (double)gamma*qj0(k,i,5,3)/qj0(k,i,1,3);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        tj0(k,i,1,3) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,1);
                        tj0(k,i,2,3) = -tursav(j,k,i,2);
                        tj0(k,i,1,4) = 2.*tj0(k,i,1,3) - tursav(jdim-1,k,i,1);
                        tj0(k,i,2,4) = 2.*tj0(k,i,2,3) - tursav(jdim-1,k,i,2);
                    }
                }
                if (ivmx == 30) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            int kk = k - ksta + 1;
                            if (iuse3 != 2) {
                                tj0(k,i,3,3) = tursav(j,k,i,3);
                            } else {
                                double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                                double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                                double dist = std::abs(smin(j,k,i));
                                double hee  = (1.+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                                double rhoee= std::pow((double)gamma*q(j,k,i,5), 1./(double)gamma);
                                double uee2 = ((double)gm1*rhoee*hee-(double)gamma*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                                double uee  = (uee2 > 0.) ? std::sqrt(uee2) : (double)xmach;
                                tj0(k,i,3,3) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
                                               uee*bcdata(kk,ii,ip,3)*(double)reue/(double)xmach;
                            }
                            tj0(k,i,3,4) = 2.*tj0(k,i,3,3) - tursav(jdim-1,k,i,3);
                        }
                    }
                } else if (ivmx == 40) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            int kk = k - ksta + 1;
                            tj0(k,i,3,3) = tursav(j,k,i,3);
                            tj0(k,i,4,3) = tursav(j,k,i,4);
                        }
                    }
                }
            }
        }

    } // end nface==4


    // **************************************************************************
    // k=1 boundary   viscous wall with T & cq specified   bctype 2004
    // **************************************************************************
    if (nface == 5) {

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;
                double cq  = bcdata(jj,ii,ip,2);
                double pb  = q(j,1,i,5);
                double dpb = q(j,2,i,5) - q(j,1,i,5);
                pb = pb - dpb/2.0;
                if (pb <= 0.0) pb = q(j,1,i,5);
                double c2 = (double)gamma * q(j,1,i,5) / q(j,1,i,1);
                if (bcdata(jj,ii,ip,1) > 0.) {
                    c2 = bcdata(jj,ii,ip,1);
                } else if (bcdata(jj,ii,ip,1) < 0.) {
                    c2 = 1.e0 + (double)gm1 * 0.5e0 * (double)xmach * (double)xmach;
                } else {
                    double xm2;
                    if (iuns > 0 && irelv > 0) {
                        xm2 = (q(j,1,i,2)-xtbk(j,i,1,1))*(q(j,1,i,2)-xtbk(j,i,1,1))
                            + (q(j,1,i,3)-xtbk(j,i,2,1))*(q(j,1,i,3)-xtbk(j,i,2,1))
                            + (q(j,1,i,4)-xtbk(j,i,3,1))*(q(j,1,i,4)-xtbk(j,i,3,1));
                    } else {
                        xm2 = q(j,1,i,2)*q(j,1,i,2)
                            + q(j,1,i,3)*q(j,1,i,3)
                            + q(j,1,i,4)*q(j,1,i,4);
                    }
                    xm2 = xm2 / c2;
                    c2  = c2 * (1. + 0.5*(double)gm1*xm2);
                }
                double uub, vvb, wwb;
                if (iuse3 == 1) {
                    uub = (double)uub_2034;
                    vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034 * sk(j,1,i,3);
                    wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034 * sk(j,1,i,2);
                } else {
                    uub = 0.; vvb = 0.; wwb = 0.;
                }
                if (iuns > 0 && irelv > 0) {
                    if (iuse3 == 1) {
                        uub = (double)uub_2034 + xtbk(j,i,1,1);
                        vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034*sk(j,1,i,3) + xtbk(j,i,2,1);
                        wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034*sk(j,1,i,2) + xtbk(j,i,3,1);
                    } else {
                        uub = xtbk(j,i,1,1); vvb = xtbk(j,i,2,1); wwb = xtbk(j,i,3,1);
                    }
                }
                uub = uub + (double)xmach*cq*sk(j,1,i,1)*c2/((double)gamma*pb);
                vvb = vvb + (double)xmach*cq*sk(j,1,i,2)*c2/((double)gamma*pb);
                wwb = wwb + (double)xmach*cq*sk(j,1,i,3)*c2/((double)gamma*pb);
                qk0(j,i,1,1) = (double)gamma*pb/c2;
                qk0(j,i,2,1) = uub; qk0(j,i,3,1) = vvb; qk0(j,i,4,1) = wwb;
                qk0(j,i,5,1) = pb;
                bck(j,i,1)   = 1.0;
                double f23 = 0.0;
                int k2 = std::min(2, kdim1);
                if (k2 == 1) f23 = 0.0;
                double z1 =   2.0 + 1.5*f23;
                double z2 =       - 0.5*f23;
                double z3 = -(2.0 +     f23);
                qk0(j,i,1,2) = z1*q(j,1,i,1) + z2*q(j,k2,i,1) + z3*qk0(j,i,1,1);
                qk0(j,i,2,2) = z1*q(j,1,i,2) + z2*q(j,k2,i,2) + z3*qk0(j,i,2,1);
                qk0(j,i,3,2) = z1*q(j,1,i,3) + z2*q(j,k2,i,3) + z3*qk0(j,i,3,1);
                qk0(j,i,4,2) = z1*q(j,1,i,4) + z2*q(j,k2,i,4) + z3*qk0(j,i,4,1);
                qk0(j,i,5,2) = z1*q(j,1,i,5) + z2*q(j,k2,i,5) + z3*qk0(j,i,5,1);
            }
        }

        if (ivmx >= 2) {
            if (level >= lglobal && ntime != 0) {
                if (iwf[2] == 0) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        vk0(j,i,1,1) = 0.; vk0(j,i,1,2) = 0.;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        int k = 1;
                        double c2b  = (double)cbar / (double)tinf;
                        double c2bp = c2b + 1.0;
                        double uu = std::sqrt(
                            (q(j,k,i,2)-qk0(j,i,2,1))*(q(j,k,i,2)-qk0(j,i,2,1)) +
                            (q(j,k,i,3)-qk0(j,i,3,1))*(q(j,k,i,3)-qk0(j,i,3,1)) +
                            (q(j,k,i,4)-qk0(j,i,4,1))*(q(j,k,i,4)-qk0(j,i,4,1)) );
                        double dist = (ivmx == 2) ? snk0(j,k,i) : std::abs(smin(j,k,i));
                        double tt   = (double)gamma*qk0(j,i,5,1)/qk0(j,i,1,1);
                        double fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                        double xnplus;
                        if (rc <= 20.24) {
                            xnplus = std::sqrt(rc);
                        } else if (rc <= 435.) {
                            xnplus = a0[0]+a0[1]*rc+(a0[2]*rc)*rc+(a0[3]*rc)*rc*rc
                                   +(a0[4]*rc)*rc*rc*rc+(a0[5]*rc)*rc*rc*rc*rc
                                   +(a0[6]*rc)*rc*rc*rc*rc*rc;
                        } else if (rc <= 4000.) {
                            xnplus = a1[0]+a1[1]*rc+(a1[2]*rc)*rc+(a1[3]*rc)*rc*rc
                                   +(a1[4]*rc)*rc*rc*rc;
                        } else {
                            xnplus = a2[0]+a2[1]*rc+a2[2]*rc*rc;
                        }
                        double xnplussav = xnplus;
                        if (xnplus > 10.) {
                            bool conv = false;
                            for (int num = 1; num <= 10; num++) {
                                double f_val = rc/xnplus - 2.44*std::log(xnplus) - 5.2;
                                double dfdn  = -rc/(xnplus*xnplus) - 2.44/xnplus;
                                double delta = -f_val/dfdn;
                                xnplus = std::abs(xnplus + delta);
                                if (std::abs(delta) < 1.e-3) { conv = true; break; }
                            }
                            if (!conv) xnplus = xnplussav;
                        }
                        double dudy = uu / dist;
                        double xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/
                                      (q(j,k,i,1)*dudy)*(double)xmach/(double)reue - 1.);
                        vk0(j,i,1,1) = xmut; vk0(j,i,1,2) = 0.;
                    }
                }
            }
        }

        if (level >= lglobal && ntime != 0) {
            if (ivmx == 4 || ivmx == 5) {
                for (int i = ista; i <= iend1; i++)
                for (int j = jsta; j <= jend1; j++) {
                    tk0(j,i,1,1) = -tursav(j,1,i,1);
                    tk0(j,i,2,1) = -tursav(j,1,i,2);
                    tk0(j,i,1,2) = 2.*tk0(j,i,1,1) - tursav(j,1,i,1);
                    tk0(j,i,2,2) = 2.*tk0(j,i,2,1) - tursav(j,1,i,2);
                }
            }
            if (ivmx >= 6) {
                double c2b  = (double)cbar / (double)tinf;
                double c2bp = c2b + 1.0;
                double re   = (double)reue / (double)xmach;
                double beta1;
                if      (ivmx == 6)                              beta1 = .075;
                else if (ivmx == 7 || ivmx == 30 || ivmx == 40) beta1 = .075;
                else if (ivmx == 8 || ivmx == 12 || ivmx == 14) beta1 = .83;
                else if (ivmx == 72)                             beta1 = .0708;
                else                                             beta1 = .075;
                int k = 1;
                if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*qk0(j,i,5,1)/qk0(j,i,1,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                        tk0(j,i,1,1) = 2.*(2.*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)
                                     - tursav(j,k,i,1);
                        tk0(j,i,2,1) = -tursav(j,k,i,2);
                        tk0(j,i,1,2) = 2.*tk0(j,i,1,1) - tursav(j,1,i,1);
                        tk0(j,i,2,2) = 2.*tk0(j,i,2,1) - tursav(j,1,i,2);
                    }
                } else if (ivmx == 15) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = 0.5*(tursav(j,k,i,2)+tursav(j,k+1,i,2))/3.;
                        tk0(j,i,1,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,1);
                        tk0(j,i,2,1) = -tursav(j,k,i,2);
                        tk0(j,i,1,2) = 0.; tk0(j,i,2,2) = 0.;
                    }
                } else if (ivmx == 16) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        tk0(j,i,1,1) = -tursav(j,k,i,1);
                        tk0(j,i,2,1) = -tursav(j,k,i,2);
                        tk0(j,i,1,2) = 0.; tk0(j,i,2,2) = 0.;
                    }
                } else if (ivmx == 70) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tke1, tke2;
                        for (int iv = 1; iv <= 6; iv++) {
                            tk0(j,i,iv,1) = -tursav(j,k,i,iv);
                            tk0(j,i,iv,2) = std::numeric_limits<double>::max();
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        tke1 = -0.5*(tursav(j,k,i,1)+tursav(j,k,i,2)+tursav(j,k,i,3));
                        tke2 = -0.5*(tursav(j,k+1,i,1)+tursav(j,k+1,i,2)+tursav(j,k+1,i,3));
                        double dkdy = 0.5*(tke1+tke2)/3.*0.5;
                        tk0(j,i,7,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,7);
                        tk0(j,i,7,2) = 0.0;
                    }
                } else if (ivmx == 72 && issglrrw2012 != 6) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        for (int iv = 1; iv <= 6; iv++) {
                            tk0(j,i,iv,1) = -tursav(j,k,i,iv);
                            tk0(j,i,iv,2) = 2*tk0(j,i,iv,1) - tursav(j,k,i,iv);
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        tk0(j,i,7,1) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,7);
                        tk0(j,i,7,2) = 2*tk0(j,i,7,1) - tursav(j,k,i,7);
                    }
                } else if (ivmx == 72 && issglrrw2012 == 6) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        for (int iv = 1; iv <= 7; iv++) {
                            tk0(j,i,iv,1) = -tursav(j,k,i,iv);
                            tk0(j,i,iv,2) = 2*tk0(j,i,iv,1) - tursav(j,k,i,iv);
                        }
                    }
                } else if (ivmx == 25) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        tk0(j,i,1,1) = 0.; tk0(j,i,2,1) = 0.;
                        tk0(j,i,1,2) = 0.; tk0(j,i,2,2) = 0.;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*qk0(j,i,5,1)/qk0(j,i,1,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        tk0(j,i,1,1) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,1);
                        tk0(j,i,2,1) = -tursav(j,k,i,2);
                        tk0(j,i,1,2) = 2.*tk0(j,i,1,1) - tursav(j,1,i,1);
                        tk0(j,i,2,2) = 2.*tk0(j,i,2,1) - tursav(j,1,i,2);
                    }
                }
                if (ivmx == 30) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            if (iuse3 != 2) {
                                tk0(j,i,3,1) = tursav(j,k,i,3);
                            } else {
                                double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                                double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                                double dist = std::abs(smin(j,k,i));
                                double hee  = (1.+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                                double rhoee= std::pow((double)gamma*q(j,k,i,5), 1./(double)gamma);
                                double uee2 = ((double)gm1*rhoee*hee-(double)gamma*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                                double uee  = (uee2 > 0.) ? std::sqrt(uee2) : (double)xmach;
                                tk0(j,i,3,1) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
                                               uee*bcdata(jj,ii,ip,3)*(double)reue/(double)xmach;
                            }
                            tk0(j,i,3,2) = 2.*tk0(j,i,3,1) - tursav(j,1,i,3);
                        }
                    }
                } else if (ivmx == 40) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            tk0(j,i,3,1) = tursav(j,k,i,3);
                            tk0(j,i,4,1) = tursav(j,k,i,4);
                        }
                    }
                }
            }
        }

    } // end nface==5


    // **************************************************************************
    // k=kdim boundary   viscous wall with T & cq specified   bctype 2004
    // **************************************************************************
    if (nface == 6) {

        for (int i = ista; i <= iend1; i++) {
            int ii = i - ista + 1;
            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;
                double cq  = bcdata(jj,ii,ip,2);
                double pb  = q(j,kdim1,i,5);
                double dpb = q(j,kdim1,i,5) - q(j,kdim1-1,i,5);
                pb = pb + dpb/2.0;
                if (pb <= 0.0) pb = q(j,kdim1,i,5);
                double c2 = (double)gamma * q(j,kdim1,i,5) / q(j,kdim1,i,1);
                if (bcdata(jj,ii,ip,1) > 0.) {
                    c2 = bcdata(jj,ii,ip,1);
                } else if (bcdata(jj,ii,ip,1) < 0.) {
                    c2 = 1.e0 + (double)gm1 * 0.5e0 * (double)xmach * (double)xmach;
                } else {
                    double xm2;
                    if (iuns > 0 && irelv > 0) {
                        xm2 = (q(j,kdim1,i,2)-xtbk(j,i,1,2))*(q(j,kdim1,i,2)-xtbk(j,i,1,2))
                            + (q(j,kdim1,i,3)-xtbk(j,i,2,2))*(q(j,kdim1,i,3)-xtbk(j,i,2,2))
                            + (q(j,kdim1,i,4)-xtbk(j,i,3,2))*(q(j,kdim1,i,4)-xtbk(j,i,3,2));
                    } else {
                        xm2 = q(j,kdim1,i,2)*q(j,kdim1,i,2)
                            + q(j,kdim1,i,3)*q(j,kdim1,i,3)
                            + q(j,kdim1,i,4)*q(j,kdim1,i,4);
                    }
                    xm2 = xm2 / c2;
                    c2  = c2 * (1. + 0.5*(double)gm1*xm2);
                }
                double uub, vvb, wwb;
                if (iuse3 == 1) {
                    uub = (double)uub_2034;
                    vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034 * sk(j,kdim,i,3);
                    wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034 * sk(j,kdim,i,2);
                } else {
                    uub = 0.; vvb = 0.; wwb = 0.;
                }
                if (iuns > 0 && irelv > 0) {
                    if (iuse3 == 1) {
                        uub = (double)uub_2034 + xtbk(j,i,1,2);
                        vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034*sk(j,kdim,i,3) + xtbk(j,i,2,2);
                        wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034*sk(j,kdim,i,2) + xtbk(j,i,3,2);
                    } else {
                        uub = xtbk(j,i,1,2); vvb = xtbk(j,i,2,2); wwb = xtbk(j,i,3,2);
                    }
                }
                uub = uub - (double)xmach*cq*sk(j,kdim,i,1)*c2/((double)gamma*pb);
                vvb = vvb - (double)xmach*cq*sk(j,kdim,i,2)*c2/((double)gamma*pb);
                wwb = wwb - (double)xmach*cq*sk(j,kdim,i,3)*c2/((double)gamma*pb);
                qk0(j,i,1,3) = (double)gamma*pb/c2;
                qk0(j,i,2,3) = uub; qk0(j,i,3,3) = vvb; qk0(j,i,4,3) = wwb;
                qk0(j,i,5,3) = pb;
                bck(j,i,2)   = 1.0;
                double f23 = 0.0;
                int k2 = std::max(1, kdim-2);
                if (k2 == 1) f23 = 0.0;
                double z1 =  -2.0 - 1.5*f23;
                double z2 =       + 0.5*f23;
                double z3 = +(2.0 +     f23);
                qk0(j,i,1,4) = z1*q(j,kdim1,i,1)+z2*q(j,k2,i,1)+z3*qk0(j,i,1,3);
                qk0(j,i,2,4) = z1*q(j,kdim1,i,2)+z2*q(j,k2,i,2)+z3*qk0(j,i,2,3);
                qk0(j,i,3,4) = z1*q(j,kdim1,i,3)+z2*q(j,k2,i,3)+z3*qk0(j,i,3,3);
                qk0(j,i,4,4) = z1*q(j,kdim1,i,4)+z2*q(j,k2,i,4)+z3*qk0(j,i,4,3);
                qk0(j,i,5,4) = z1*q(j,kdim1,i,5)+z2*q(j,k2,i,5)+z3*qk0(j,i,5,3);
            }
        }

        if (ivmx >= 2) {
            if (level >= lglobal && ntime != 0) {
                if (iwf[2] == 0) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        vk0(j,i,1,3) = 0.; vk0(j,i,1,4) = 0.;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        int k = kdim - 1;
                        double c2b  = (double)cbar / (double)tinf;
                        double c2bp = c2b + 1.0;
                        double uu = std::sqrt(
                            (q(j,k,i,2)-qk0(j,i,2,3))*(q(j,k,i,2)-qk0(j,i,2,3)) +
                            (q(j,k,i,3)-qk0(j,i,3,3))*(q(j,k,i,3)-qk0(j,i,3,3)) +
                            (q(j,k,i,4)-qk0(j,i,4,3))*(q(j,k,i,4)-qk0(j,i,4,3)) );
                        double dist = (ivmx == 2) ? snkm(j,k,i) : std::abs(smin(j,k,i));
                        double tt   = (double)gamma*qk0(j,i,5,3)/qk0(j,i,1,3);
                        double fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                        double xnplus;
                        if (rc <= 20.24) {
                            xnplus = std::sqrt(rc);
                        } else if (rc <= 435.) {
                            xnplus = a0[0]+a0[1]*rc+(a0[2]*rc)*rc+(a0[3]*rc)*rc*rc
                                   +(a0[4]*rc)*rc*rc*rc+(a0[5]*rc)*rc*rc*rc*rc
                                   +(a0[6]*rc)*rc*rc*rc*rc*rc;
                        } else if (rc <= 4000.) {
                            xnplus = a1[0]+a1[1]*rc+(a1[2]*rc)*rc+(a1[3]*rc)*rc*rc
                                   +(a1[4]*rc)*rc*rc*rc;
                        } else {
                            xnplus = a2[0]+a2[1]*rc+a2[2]*rc*rc;
                        }
                        double xnplussav = xnplus;
                        if (xnplus > 10.) {
                            bool conv = false;
                            for (int num = 1; num <= 10; num++) {
                                double f_val = rc/xnplus - 2.44*std::log(xnplus) - 5.2;
                                double dfdn  = -rc/(xnplus*xnplus) - 2.44/xnplus;
                                double delta = -f_val/dfdn;
                                xnplus = std::abs(xnplus + delta);
                                if (std::abs(delta) < 1.e-3) { conv = true; break; }
                            }
                            if (!conv) xnplus = xnplussav;
                        }
                        double dudy = uu / dist;
                        double xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/
                                      (q(j,k,i,1)*dudy)*(double)xmach/(double)reue - 1.);
                        vk0(j,i,1,3) = xmut; vk0(j,i,1,4) = 0.;
                    }
                }
            }
        }

        if (level >= lglobal && ntime != 0) {
            if (ivmx == 4 || ivmx == 5) {
                for (int i = ista; i <= iend1; i++)
                for (int j = jsta; j <= jend1; j++) {
                    tk0(j,i,1,3) = -tursav(j,kdim-1,i,1);
                    tk0(j,i,2,3) = -tursav(j,kdim-1,i,2);
                    tk0(j,i,1,4) = 2.*tk0(j,i,1,3) - tursav(j,kdim-1,i,1);
                    tk0(j,i,2,4) = 2.*tk0(j,i,2,3) - tursav(j,kdim-1,i,2);
                }
            }
            if (ivmx >= 6) {
                double c2b  = (double)cbar / (double)tinf;
                double c2bp = c2b + 1.0;
                double re   = (double)reue / (double)xmach;
                double beta1;
                if      (ivmx == 6)                              beta1 = .075;
                else if (ivmx == 7 || ivmx == 30 || ivmx == 40) beta1 = .075;
                else if (ivmx == 8 || ivmx == 12 || ivmx == 14) beta1 = .83;
                else if (ivmx == 72)                             beta1 = .0708;
                else                                             beta1 = .075;
                int k = kdim - 1;
                if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*qk0(j,i,5,3)/qk0(j,i,1,3);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                        tk0(j,i,1,3) = 2.*(2.*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)
                                     - tursav(j,k,i,1);
                        tk0(j,i,2,3) = -tursav(j,k,i,2);
                        tk0(j,i,1,4) = 2.*tk0(j,i,1,3) - tursav(j,kdim-1,i,1);
                        tk0(j,i,2,4) = 2.*tk0(j,i,2,3) - tursav(j,kdim-1,i,2);
                    }
                } else if (ivmx == 15) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = 0.5*(tursav(j,k,i,2)+tursav(j,k-1,i,2))/3.;
                        tk0(j,i,1,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,1);
                        tk0(j,i,2,3) = -tursav(j,k,i,2);
                        tk0(j,i,1,4) = 0.; tk0(j,i,2,4) = 0.;
                    }
                } else if (ivmx == 16) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        tk0(j,i,1,3) = -tursav(j,k,i,1);
                        tk0(j,i,2,3) = -tursav(j,k,i,2);
                        tk0(j,i,1,4) = 0.; tk0(j,i,2,4) = 0.;
                    }
                } else if (ivmx == 70) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tke1, tke2;
                        for (int iv = 1; iv <= 6; iv++) {
                            tk0(j,i,iv,3) = -tursav(j,k,i,iv);
                            tk0(j,i,iv,4) = std::numeric_limits<double>::max();
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        tke1 = -0.5*(tursav(j,k,i,1)+tursav(j,k,i,2)+tursav(j,k,i,3));
                        tke2 = -0.5*(tursav(j,k-1,i,1)+tursav(j,k-1,i,2)+tursav(j,k-1,i,3));
                        double dkdy = 0.5*(tke1+tke2)/3.*0.5;
                        tk0(j,i,7,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,7);
                        tk0(j,i,7,4) = 0.0;
                    }
                } else if (ivmx == 72 && issglrrw2012 != 6) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        for (int iv = 1; iv <= 6; iv++) {
                            tk0(j,i,iv,3) = -tursav(j,k,i,iv);
                            tk0(j,i,iv,4) = 2*tk0(j,i,iv,3) - tursav(j,k,i,iv);
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        tk0(j,i,7,3) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,7);
                        tk0(j,i,7,4) = 2*tk0(j,i,7,3) - tursav(j,k,i,7);
                    }
                } else if (ivmx == 72 && issglrrw2012 == 6) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        for (int iv = 1; iv <= 7; iv++) {
                            tk0(j,i,iv,3) = -tursav(j,k,i,iv);
                            tk0(j,i,iv,4) = 2*tk0(j,i,iv,3) - tursav(j,k,i,iv);
                        }
                    }
                } else if (ivmx == 25) {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        tk0(j,i,1,3) = 0.; tk0(j,i,2,3) = 0.;
                        tk0(j,i,1,4) = 0.; tk0(j,i,2,4) = 0.;
                    }
                } else {
                    for (int i = ista; i <= iend1; i++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*qk0(j,i,5,3)/qk0(j,i,1,3);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        tk0(j,i,1,3) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,1);
                        tk0(j,i,2,3) = -tursav(j,k,i,2);
                        tk0(j,i,1,4) = 2.*tk0(j,i,1,3) - tursav(j,kdim-1,i,1);
                        tk0(j,i,2,4) = 2.*tk0(j,i,2,3) - tursav(j,kdim-1,i,2);
                    }
                }
                if (ivmx == 30) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            if (iuse3 != 2) {
                                tk0(j,i,3,3) = tursav(j,k,i,3);
                            } else {
                                double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                                double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                                double dist = std::abs(smin(j,k,i));
                                double hee  = (1.+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                                double rhoee= std::pow((double)gamma*q(j,k,i,5), 1./(double)gamma);
                                double uee2 = ((double)gm1*rhoee*hee-(double)gamma*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                                double uee  = (uee2 > 0.) ? std::sqrt(uee2) : (double)xmach;
                                tk0(j,i,3,3) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
                                               uee*bcdata(jj,ii,ip,3)*(double)reue/(double)xmach;
                            }
                            tk0(j,i,3,4) = 2.*tk0(j,i,3,3) - tursav(j,kdim-1,i,3);
                        }
                    }
                } else if (ivmx == 40) {
                    for (int i = ista; i <= iend1; i++) {
                        int ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            tk0(j,i,3,3) = tursav(j,k,i,3);
                            tk0(j,i,4,3) = tursav(j,k,i,4);
                        }
                    }
                }
            }
        }

    } // end nface==6


    // **************************************************************************
    // i=1 boundary   viscous wall with T & cq specified   bctype 2004
    // **************************************************************************
    if (nface == 1) {

        for (int k = ksta; k <= kend1; k++) {
            int kk = k - ksta + 1;
            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;
                double cq  = bcdata(jj,kk,ip,2);
                double pb  = q(j,k,1,5);
                double dpb = q(j,k,2,5) - q(j,k,1,5);
                pb = pb - dpb/2.0;
                if (pb <= 0.0) pb = q(j,k,1,5);
                double c2 = (double)gamma * q(j,k,1,5) / q(j,k,1,1);
                if (bcdata(jj,kk,ip,1) > 0.) {
                    c2 = bcdata(jj,kk,ip,1);
                } else if (bcdata(jj,kk,ip,1) < 0.) {
                    c2 = 1.e0 + (double)gm1 * 0.5e0 * (double)xmach * (double)xmach;
                } else {
                    double xm2;
                    if (iuns > 0 && irelv > 0) {
                        xm2 = (q(j,k,1,2)-xtbi(j,k,1,1))*(q(j,k,1,2)-xtbi(j,k,1,1))
                            + (q(j,k,1,3)-xtbi(j,k,2,1))*(q(j,k,1,3)-xtbi(j,k,2,1))
                            + (q(j,k,1,4)-xtbi(j,k,3,1))*(q(j,k,1,4)-xtbi(j,k,3,1));
                    } else {
                        xm2 = q(j,k,1,2)*q(j,k,1,2)
                            + q(j,k,1,3)*q(j,k,1,3)
                            + q(j,k,1,4)*q(j,k,1,4);
                    }
                    xm2 = xm2 / c2;
                    c2  = c2 * (1. + 0.5*(double)gm1*xm2);
                }
                double uub, vvb, wwb;
                if (iuse3 == 1) {
                    uub = (double)uub_2034;
                    vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034 * si(j,k,1,3);
                    wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034 * si(j,k,1,2);
                } else {
                    uub = 0.; vvb = 0.; wwb = 0.;
                }
                if (iuns > 0 && irelv > 0) {
                    if (iuse3 == 1) {
                        uub = (double)uub_2034 + xtbi(j,k,1,1);
                        vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034*si(j,k,1,3) + xtbi(j,k,2,1);
                        wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034*si(j,k,1,2) + xtbi(j,k,3,1);
                    } else {
                        uub = xtbi(j,k,1,1); vvb = xtbi(j,k,2,1); wwb = xtbi(j,k,3,1);
                    }
                }
                uub = uub + (double)xmach*cq*si(j,k,1,1)*c2/((double)gamma*pb);
                vvb = vvb + (double)xmach*cq*si(j,k,1,2)*c2/((double)gamma*pb);
                wwb = wwb + (double)xmach*cq*si(j,k,1,3)*c2/((double)gamma*pb);
                qi0(j,k,1,1) = (double)gamma*pb/c2;
                qi0(j,k,2,1) = uub; qi0(j,k,3,1) = vvb; qi0(j,k,4,1) = wwb;
                qi0(j,k,5,1) = pb;
                bci(j,k,1)   = 1.0;
                double f23 = 0.0;
                int i2 = std::min(2, idim1);
                if (i2 == 1) f23 = 0.0;
                double z1 =   2.0 + 1.5*f23;
                double z2 =       - 0.5*f23;
                double z3 = -(2.0 +     f23);
                qi0(j,k,1,2) = z1*q(j,k,1,1) + z2*q(j,k,i2,1) + z3*qi0(j,k,1,1);
                qi0(j,k,2,2) = z1*q(j,k,1,2) + z2*q(j,k,i2,2) + z3*qi0(j,k,2,1);
                qi0(j,k,3,2) = z1*q(j,k,1,3) + z2*q(j,k,i2,3) + z3*qi0(j,k,3,1);
                qi0(j,k,4,2) = z1*q(j,k,1,4) + z2*q(j,k,i2,4) + z3*qi0(j,k,4,1);
                qi0(j,k,5,2) = z1*q(j,k,1,5) + z2*q(j,k,i2,5) + z3*qi0(j,k,5,1);
            }
        }

        if (ivmx >= 2) {
            if (level >= lglobal && ntime != 0) {
                if (iwf[0] == 0) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        vi0(j,k,1,1) = 0.; vi0(j,k,1,2) = 0.;
                    }
                } else {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        int i = 1;
                        double c2b  = (double)cbar / (double)tinf;
                        double c2bp = c2b + 1.0;
                        double uu = std::sqrt(
                            (q(j,k,i,2)-qi0(j,k,2,1))*(q(j,k,i,2)-qi0(j,k,2,1)) +
                            (q(j,k,i,3)-qi0(j,k,3,1))*(q(j,k,i,3)-qi0(j,k,3,1)) +
                            (q(j,k,i,4)-qi0(j,k,4,1))*(q(j,k,i,4)-qi0(j,k,4,1)) );
                        double dist = (ivmx == 2) ? sni0(j,k,i) : std::abs(smin(j,k,i));
                        double tt   = (double)gamma*qi0(j,k,5,1)/qi0(j,k,1,1);
                        double fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                        double xnplus;
                        if (rc <= 20.24) {
                            xnplus = std::sqrt(rc);
                        } else if (rc <= 435.) {
                            xnplus = a0[0]+a0[1]*rc+(a0[2]*rc)*rc+(a0[3]*rc)*rc*rc
                                   +(a0[4]*rc)*rc*rc*rc+(a0[5]*rc)*rc*rc*rc*rc
                                   +(a0[6]*rc)*rc*rc*rc*rc*rc;
                        } else if (rc <= 4000.) {
                            xnplus = a1[0]+a1[1]*rc+(a1[2]*rc)*rc+(a1[3]*rc)*rc*rc
                                   +(a1[4]*rc)*rc*rc*rc;
                        } else {
                            xnplus = a2[0]+a2[1]*rc+a2[2]*rc*rc;
                        }
                        double xnplussav = xnplus;
                        if (xnplus > 10.) {
                            bool conv = false;
                            for (int num = 1; num <= 10; num++) {
                                double f_val = rc/xnplus - 2.44*std::log(xnplus) - 5.2;
                                double dfdn  = -rc/(xnplus*xnplus) - 2.44/xnplus;
                                double delta = -f_val/dfdn;
                                xnplus = std::abs(xnplus + delta);
                                if (std::abs(delta) < 1.e-3) { conv = true; break; }
                            }
                            if (!conv) xnplus = xnplussav;
                        }
                        double dudy = uu / dist;
                        double xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/
                                      (q(j,k,i,1)*dudy)*(double)xmach/(double)reue - 1.);
                        vi0(j,k,1,1) = xmut; vi0(j,k,1,2) = 0.;
                    }
                }
            }
        }

        if (level >= lglobal && ntime != 0) {
            if (ivmx == 4 || ivmx == 5) {
                for (int k = ksta; k <= kend1; k++)
                for (int j = jsta; j <= jend1; j++) {
                    ti0(j,k,1,1) = -tursav(j,k,1,1);
                    ti0(j,k,2,1) = -tursav(j,k,1,2);
                    ti0(j,k,1,2) = 2.*ti0(j,k,1,1) - tursav(j,k,1,1);
                    ti0(j,k,2,2) = 2.*ti0(j,k,2,1) - tursav(j,k,1,2);
                }
            }
            if (ivmx >= 6) {
                double c2b  = (double)cbar / (double)tinf;
                double c2bp = c2b + 1.0;
                double re   = (double)reue / (double)xmach;
                double beta1;
                if      (ivmx == 6)                              beta1 = .075;
                else if (ivmx == 7 || ivmx == 30 || ivmx == 40) beta1 = .075;
                else if (ivmx == 8 || ivmx == 12 || ivmx == 14) beta1 = .83;
                else if (ivmx == 72)                             beta1 = .0708;
                else                                             beta1 = .075;
                int i = 1;
                if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*qi0(j,k,5,1)/qi0(j,k,1,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                        ti0(j,k,1,1) = 2.*(2.*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)
                                     - tursav(j,k,i,1);
                        ti0(j,k,2,1) = -tursav(j,k,i,2);
                        ti0(j,k,1,2) = 2.*ti0(j,k,1,1) - tursav(j,k,1,1);
                        ti0(j,k,2,2) = 2.*ti0(j,k,2,1) - tursav(j,k,1,2);
                    }
                } else if (ivmx == 15) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = 0.5*(tursav(j,k,i,2)+tursav(j,k,i+1,2))/3.;
                        ti0(j,k,1,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,1);
                        ti0(j,k,2,1) = -tursav(j,k,i,2);
                        ti0(j,k,1,2) = 0.; ti0(j,k,2,2) = 0.;
                    }
                } else if (ivmx == 16) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        ti0(j,k,1,1) = -tursav(j,k,i,1);
                        ti0(j,k,2,1) = -tursav(j,k,i,2);
                        ti0(j,k,1,2) = 0.; ti0(j,k,2,2) = 0.;
                    }
                } else if (ivmx == 70) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tke1, tke2;
                        for (int iv = 1; iv <= 6; iv++) {
                            ti0(j,k,iv,1) = -tursav(j,k,i,iv);
                            ti0(j,k,iv,2) = std::numeric_limits<double>::max();
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        tke1 = -0.5*(tursav(j,k,i,1)+tursav(j,k,i,2)+tursav(j,k,i,3));
                        tke2 = -0.5*(tursav(j,k,i+1,1)+tursav(j,k,i+1,2)+tursav(j,k,i+1,3));
                        double dkdy = 0.5*(tke1+tke2)/3.*0.5;
                        ti0(j,k,7,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,7);
                        ti0(j,k,7,2) = 0.0;
                    }
                } else if (ivmx == 72 && issglrrw2012 != 6) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        for (int iv = 1; iv <= 6; iv++) {
                            ti0(j,k,iv,1) = -tursav(j,k,i,iv);
                            ti0(j,k,iv,2) = 2.*ti0(j,k,iv,1) - tursav(j,k,i,iv);
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        ti0(j,k,7,1) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,7);
                        ti0(j,k,7,2) = 2.*ti0(j,k,7,1) - tursav(j,k,i,7);
                    }
                } else if (ivmx == 72 && issglrrw2012 == 6) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        for (int iv = 1; iv <= 7; iv++) {
                            ti0(j,k,iv,1) = -tursav(j,k,i,iv);
                            ti0(j,k,iv,2) = 2.*ti0(j,k,iv,1) - tursav(j,k,i,iv);
                        }
                    }
                } else if (ivmx == 25) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        ti0(j,k,1,1) = 0.; ti0(j,k,2,1) = 0.;
                        ti0(j,k,1,2) = 0.; ti0(j,k,2,2) = 0.;
                    }
                } else {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*qi0(j,k,5,1)/qi0(j,k,1,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        ti0(j,k,1,1) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,1);
                        ti0(j,k,2,1) = -tursav(j,k,i,2);
                        ti0(j,k,1,2) = 2.*ti0(j,k,1,1) - tursav(j,k,1,1);
                        ti0(j,k,2,2) = 2.*ti0(j,k,2,1) - tursav(j,k,1,2);
                    }
                }
                if (ivmx == 30) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            if (iuse3 != 2) {
                                ti0(j,k,3,1) = tursav(j,k,i,3);
                            } else {
                                double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                                double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                                double dist = std::abs(smin(j,k,i));
                                double hee  = (1.+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                                double rhoee= std::pow((double)gamma*q(j,k,i,5), 1./(double)gamma);
                                double uee2 = ((double)gm1*rhoee*hee-(double)gamma*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                                double uee  = (uee2 > 0.) ? std::sqrt(uee2) : (double)xmach;
                                ti0(j,k,3,1) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
                                               uee*bcdata(jj,kk,ip,3)*(double)reue/(double)xmach;
                            }
                            ti0(j,k,3,2) = 2.*ti0(j,k,3,1) - tursav(j,k,1,3);
                        }
                    }
                } else if (ivmx == 40) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            ti0(j,k,3,1) = tursav(j,k,i,3);
                            ti0(j,k,4,1) = tursav(j,k,i,4);
                        }
                    }
                }
            }
        }

    } // end nface==1


    // **************************************************************************
    // i=idim boundary   viscous wall with T & cq specified   bctype 2004
    // **************************************************************************
    if (nface == 2) {

        for (int k = ksta; k <= kend1; k++) {
            int kk = k - ksta + 1;
            for (int j = jsta; j <= jend1; j++) {
                int jj = j - jsta + 1;
                double cq  = bcdata(jj,kk,ip,2);
                double pb  = q(j,k,idim1,5);
                double dpb = q(j,k,idim1,5) - q(j,k,idim1-1,5);
                pb = pb + dpb/2.0;
                if (pb <= 0.0) pb = q(j,k,idim1,5);
                double c2 = (double)gamma * q(j,k,idim1,5) / q(j,k,idim1,1);
                if (bcdata(jj,kk,ip,1) > 0.) {
                    c2 = bcdata(jj,kk,ip,1);
                } else if (bcdata(jj,kk,ip,1) < 0.) {
                    c2 = 1.e0 + (double)gm1 * 0.5e0 * (double)xmach * (double)xmach;
                } else {
                    double xm2;
                    if (iuns > 0 && irelv > 0) {
                        xm2 = (q(j,k,idim1,2)-xtbi(j,k,1,2))*(q(j,k,idim1,2)-xtbi(j,k,1,2))
                            + (q(j,k,idim1,3)-xtbi(j,k,2,2))*(q(j,k,idim1,3)-xtbi(j,k,2,2))
                            + (q(j,k,idim1,4)-xtbi(j,k,3,2))*(q(j,k,idim1,4)-xtbi(j,k,3,2));
                    } else {
                        xm2 = q(j,k,idim1,2)*q(j,k,idim1,2)
                            + q(j,k,idim1,3)*q(j,k,idim1,3)
                            + q(j,k,idim1,4)*q(j,k,idim1,4);
                    }
                    xm2 = xm2 / c2;
                    c2  = c2 * (1. + 0.5*(double)gm1*xm2);
                }
                double uub, vvb, wwb;
                if (iuse3 == 1) {
                    uub = (double)uub_2034;
                    vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034 * si(j,k,idim,3);
                    wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034 * si(j,k,idim,2);
                } else {
                    uub = 0.; vvb = 0.; wwb = 0.;
                }
                if (iuns > 0 && irelv > 0) {
                    if (iuse3 == 1) {
                        uub = (double)uub_2034 + xtbi(j,k,1,2);
                        vvb = (double)vvb_2034 + (double)vvb_xaxisrot2034*si(j,k,idim,3) + xtbi(j,k,2,2);
                        wwb = (double)wwb_2034 - (double)vvb_xaxisrot2034*si(j,k,idim,2) + xtbi(j,k,3,2);
                    } else {
                        uub = xtbi(j,k,1,2); vvb = xtbi(j,k,2,2); wwb = xtbi(j,k,3,2);
                    }
                }
                uub = uub - (double)xmach*cq*si(j,k,idim,1)*c2/((double)gamma*pb);
                vvb = vvb - (double)xmach*cq*si(j,k,idim,2)*c2/((double)gamma*pb);
                wwb = wwb - (double)xmach*cq*si(j,k,idim,3)*c2/((double)gamma*pb);
                qi0(j,k,1,3) = (double)gamma*pb/c2;
                qi0(j,k,2,3) = uub; qi0(j,k,3,3) = vvb; qi0(j,k,4,3) = wwb;
                qi0(j,k,5,3) = pb;
                bci(j,k,2)   = 1.0;
                double f23 = 0.0;
                int i2 = std::max(1, idim-2);
                if (i2 == 1) f23 = 0.0;
                double z1 =  -2.0 - 1.5*f23;
                double z2 =       + 0.5*f23;
                double z3 = +(2.0 +     f23);
                qi0(j,k,1,4) = z1*q(j,k,idim1,1)+z2*q(j,k,i2,1)+z3*qi0(j,k,1,3);
                qi0(j,k,2,4) = z1*q(j,k,idim1,2)+z2*q(j,k,i2,2)+z3*qi0(j,k,2,3);
                qi0(j,k,3,4) = z1*q(j,k,idim1,3)+z2*q(j,k,i2,3)+z3*qi0(j,k,3,3);
                qi0(j,k,4,4) = z1*q(j,k,idim1,4)+z2*q(j,k,i2,4)+z3*qi0(j,k,4,3);
                qi0(j,k,5,4) = z1*q(j,k,idim1,5)+z2*q(j,k,i2,5)+z3*qi0(j,k,5,3);
            }
        }

        if (ivmx >= 2) {
            if (level >= lglobal && ntime != 0) {
                if (iwf[0] == 0) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        vi0(j,k,1,3) = 0.; vi0(j,k,1,4) = 0.;
                    }
                } else {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        int i = idim - 1;
                        double c2b  = (double)cbar / (double)tinf;
                        double c2bp = c2b + 1.0;
                        double uu = std::sqrt(
                            (q(j,k,i,2)-qi0(j,k,2,3))*(q(j,k,i,2)-qi0(j,k,2,3)) +
                            (q(j,k,i,3)-qi0(j,k,3,3))*(q(j,k,i,3)-qi0(j,k,3,3)) +
                            (q(j,k,i,4)-qi0(j,k,4,3))*(q(j,k,i,4)-qi0(j,k,4,3)) );
                        double dist = (ivmx == 2) ? snim(j,k,i) : std::abs(smin(j,k,i));
                        double tt   = (double)gamma*qi0(j,k,5,3)/qi0(j,k,1,3);
                        double fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                        double xnplus;
                        if (rc <= 20.24) {
                            xnplus = std::sqrt(rc);
                        } else if (rc <= 435.) {
                            xnplus = a0[0]+a0[1]*rc+(a0[2]*rc)*rc+(a0[3]*rc)*rc*rc
                                   +(a0[4]*rc)*rc*rc*rc+(a0[5]*rc)*rc*rc*rc*rc
                                   +(a0[6]*rc)*rc*rc*rc*rc*rc;
                        } else if (rc <= 4000.) {
                            xnplus = a1[0]+a1[1]*rc+(a1[2]*rc)*rc+(a1[3]*rc)*rc*rc
                                   +(a1[4]*rc)*rc*rc*rc;
                        } else {
                            xnplus = a2[0]+a2[1]*rc+a2[2]*rc*rc;
                        }
                        double xnplussav = xnplus;
                        if (xnplus > 10.) {
                            bool conv = false;
                            for (int num = 1; num <= 10; num++) {
                                double f_val = rc/xnplus - 2.44*std::log(xnplus) - 5.2;
                                double dfdn  = -rc/(xnplus*xnplus) - 2.44/xnplus;
                                double delta = -f_val/dfdn;
                                xnplus = std::abs(xnplus + delta);
                                if (std::abs(delta) < 1.e-3) { conv = true; break; }
                            }
                            if (!conv) xnplus = xnplussav;
                        }
                        double dudy = uu / dist;
                        double xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/
                                      (q(j,k,i,1)*dudy)*(double)xmach/(double)reue - 1.);
                        vi0(j,k,1,3) = xmut; vi0(j,k,1,4) = 0.;
                    }
                }
            }
        }

        if (level >= lglobal && ntime != 0) {
            if (ivmx == 4 || ivmx == 5) {
                for (int k = ksta; k <= kend1; k++)
                for (int j = jsta; j <= jend1; j++) {
                    ti0(j,k,1,3) = -tursav(j,k,idim-1,1);
                    ti0(j,k,2,3) = -tursav(j,k,idim-1,2);
                    ti0(j,k,1,4) = 2.*ti0(j,k,1,3) - tursav(j,k,idim-1,1);
                    ti0(j,k,2,4) = 2.*ti0(j,k,2,3) - tursav(j,k,idim-1,2);
                }
            }
            if (ivmx >= 6) {
                double c2b  = (double)cbar / (double)tinf;
                double c2bp = c2b + 1.0;
                double re   = (double)reue / (double)xmach;
                double beta1;
                if      (ivmx == 6)                              beta1 = .075;
                else if (ivmx == 7 || ivmx == 30 || ivmx == 40) beta1 = .075;
                else if (ivmx == 8 || ivmx == 12 || ivmx == 14) beta1 = .83;
                else if (ivmx == 72)                             beta1 = .0708;
                else                                             beta1 = .075;
                int i = idim - 1;
                if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*qi0(j,k,5,3)/qi0(j,k,1,3);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                        ti0(j,k,1,3) = 2.*(2.*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)
                                     - tursav(j,k,i,1);
                        ti0(j,k,2,3) = -tursav(j,k,i,2);
                        ti0(j,k,1,4) = 2.*ti0(j,k,1,3) - tursav(j,k,idim-1,1);
                        ti0(j,k,2,4) = 2.*ti0(j,k,2,3) - tursav(j,k,idim-1,2);
                    }
                } else if (ivmx == 15) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dkdy = 0.5*(tursav(j,k,i,2)+tursav(j,k,i-1,2))/3.;
                        ti0(j,k,1,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,1);
                        ti0(j,k,2,3) = -tursav(j,k,i,2);
                        ti0(j,k,1,4) = 0.; ti0(j,k,2,4) = 0.;
                    }
                } else if (ivmx == 16) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        ti0(j,k,1,3) = -tursav(j,k,i,1);
                        ti0(j,k,2,3) = -tursav(j,k,i,2);
                        ti0(j,k,1,4) = 0.; ti0(j,k,2,4) = 0.;
                    }
                } else if (ivmx == 70) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tke1, tke2;
                        for (int iv = 1; iv <= 6; iv++) {
                            ti0(j,k,iv,3) = -tursav(j,k,i,iv);
                            ti0(j,k,iv,4) = std::numeric_limits<double>::max();
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        tke1 = -0.5*(tursav(j,k,i,1)+tursav(j,k,i,2)+tursav(j,k,i,3));
                        tke2 = -0.5*(tursav(j,k,i-1,1)+tursav(j,k,i-1,2)+tursav(j,k,i-1,3));
                        double dkdy = 0.5*(tke1+tke2)/3.*0.5;
                        ti0(j,k,7,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i)) - tursav(j,k,i,7);
                        ti0(j,k,7,4) = 0.0;
                    }
                } else if (ivmx == 72 && issglrrw2012 != 6) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        for (int iv = 1; iv <= 6; iv++) {
                            ti0(j,k,iv,3) = -tursav(j,k,i,iv);
                            ti0(j,k,iv,4) = 2.*ti0(j,k,iv,3) - tursav(j,k,i,iv);
                        }
                        double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        ti0(j,k,7,3) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,7);
                        ti0(j,k,7,4) = 2.*ti0(j,k,7,3) - tursav(j,k,i,7);
                    }
                } else if (ivmx == 72 && issglrrw2012 == 6) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        for (int iv = 1; iv <= 7; iv++) {
                            ti0(j,k,iv,3) = -tursav(j,k,i,iv);
                            ti0(j,k,iv,4) = 2.*ti0(j,k,iv,3) - tursav(j,k,i,iv);
                        }
                    }
                } else if (ivmx == 25) {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        ti0(j,k,1,3) = 0.; ti0(j,k,2,3) = 0.;
                        ti0(j,k,1,4) = 0.; ti0(j,k,2,4) = 0.;
                    }
                } else {
                    for (int k = ksta; k <= kend1; k++)
                    for (int j = jsta; j <= jend1; j++) {
                        double tt   = (double)gamma*qi0(j,k,5,3)/qi0(j,k,1,3);
                        double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        double dist = std::abs(smin(j,k,i));
                        ti0(j,k,1,3) = 2.*(60.*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))
                                     - tursav(j,k,i,1);
                        ti0(j,k,2,3) = -tursav(j,k,i,2);
                        ti0(j,k,1,4) = 2.*ti0(j,k,1,3) - tursav(j,k,idim-1,1);
                        ti0(j,k,2,4) = 2.*ti0(j,k,2,3) - tursav(j,k,idim-1,2);
                    }
                }
                if (ivmx == 30) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            if (iuse3 != 2) {
                                ti0(j,k,3,3) = tursav(j,k,i,3);
                            } else {
                                double tt   = (double)gamma*q(j,k,i,5)/q(j,k,i,1);
                                double fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                                double dist = std::abs(smin(j,k,i));
                                double hee  = (1.+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                                double rhoee= std::pow((double)gamma*q(j,k,i,5), 1./(double)gamma);
                                double uee2 = ((double)gm1*rhoee*hee-(double)gamma*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                                double uee  = (uee2 > 0.) ? std::sqrt(uee2) : (double)xmach;
                                ti0(j,k,3,3) = tursav(j,k,i,3) + 2.*dist*q(j,k,i,1)/fnu*
                                               uee*bcdata(jj,kk,ip,3)*(double)reue/(double)xmach;
                            }
                            ti0(j,k,3,4) = 2.*ti0(j,k,3,3) - tursav(j,k,idim-1,3);
                        }
                    }
                } else if (ivmx == 40) {
                    for (int k = ksta; k <= kend1; k++) {
                        int kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            int jj = j - jsta + 1;
                            ti0(j,k,3,3) = tursav(j,k,i,3);
                            ti0(j,k,4,3) = tursav(j,k,i,4);
                        }
                    }
                }
            }
        }

    } // end nface==2

} // end bc2004()

} // namespace bc2004_ns
