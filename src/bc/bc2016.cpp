// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2016.h"
#include "bc.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace bc2016_ns {

// COMMON block aliases
static inline float& gamma_   = cmn_fluid.gamma;
static inline float& gm1      = cmn_fluid.gm1;
static inline float& gp1      = cmn_fluid.gp1;
static inline float& gm1g     = cmn_fluid.gm1g;
static inline float& gp1g     = cmn_fluid.gp1g;
static inline float& ggm1     = cmn_fluid.ggm1;
static inline float& pr       = cmn_fluid2.pr;
static inline float& prt      = cmn_fluid2.prt;
static inline float& cbar     = cmn_fluid2.cbar;
static inline float& xmach    = cmn_info.xmach;
static inline int32_t& ivmx   = cmn_maxiv.ivmx;
static inline int32_t& level  = cmn_mgrd.level;
static inline int32_t& lglobal= cmn_mgrd.lglobal;
static inline float& reue     = cmn_reyue.reue;
static inline float& tinf     = cmn_reyue.tinf;
static inline float& time_    = cmn_unst.time;

// ---------------------------------------------------------------------------
// bc: wrapper calling bc_ns::bc
// ---------------------------------------------------------------------------
void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
        int& nwork, double& cl, FortranArray1DRef<int> nou,
        FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
        int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg,
        FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0,
        FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim,
        FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
        FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek,
        FortranArray3DRef<int> lwdat, int& myid,
        FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg,
        FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, int& nummem)
{
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf,
              ibufdim, maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg,
              nblg, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
              ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek,
              lwdat, myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

// ---------------------------------------------------------------------------
// bc2016: solid wall BC (viscous wall) with T & cq specified (bctype 2016)
// ---------------------------------------------------------------------------
void bc2016(int& jdim, int& kdim, int& idim,
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
            int& nbuf, int& ibufdim, int& myid, int& nummem)
{
    // data statement (for use with wall functions)
    static const double a0[8] = {0.0,
        2.354039, 0.1179840, -4.2899192e-04, 2.0404148e-06,
        -5.1775775e-09, 6.2687308e-12, -2.916958e-15};
    static const double a1[6] = {0.0,
        5.777191, 6.8756983e-02, -7.1582745e-06, 1.5594904e-09,
        -1.4865778e-13};
    static const double a2[4] = {0.0,
        31.08654, 5.0429072e-02, -2.0072314e-8};

    // iuse3 is an implicit integer variable, never declared/assigned, defaults to 0
    int iuse3 = 0;

    // local variables
    double twopi, pb, dpb, c2, xm2, uub, vvb, wwb;
    double cq, cqu, sjetx, sjety, sjetz, rfreq, cqt, tcqt, sjet;
    double f23, z1, z2, z3;
    double c2b, c2bp, uu, dist, tt, fnuw, rc, xnplus, xnplussav;
    double f_val, dfdn, delta, dudy, xmut;
    double fnu, dkdy, re, beta1;
    double hee, rhoee, uee2, uee;
    int jdim1, kdim1, idim1, jend1, kend1, iend1;
    int ip, i, j, k, ii, jj, kk, j2, k2, i2, num;

    twopi = 2.0 * std::acos(-1.0);
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    jend1 = jend - 1;
    kend1 = kend - 1;
    iend1 = iend - 1;

    // this bc makes use of only one plane of data
    ip = 1;



    // **************************************************************************
    // j=1 boundary   viscous wall with T & cq specified   bctype 2016
    // **************************************************************************
    if (nface == 3) {

    for (i = ista; i <= iend1; i++) {
        ii = i - ista + 1;

        for (k = ksta; k <= kend1; k++) {
            kk = k - ksta + 1;
            cq    = bcdata(kk,ii,ip,2);
            cqu   = bcdata(kk,ii,ip,3);
            sjetx = bcdata(kk,ii,ip,4);
            sjety = bcdata(kk,ii,ip,5);
            sjetz = bcdata(kk,ii,ip,6);
            rfreq = bcdata(kk,ii,ip,7);
            pb    = q(1,k,i,5);
            dpb   = q(2,k,i,5) - q(1,k,i,5);
            pb    = pb - dpb / 2.0;
            if ((float)pb <= 0.0f) pb = q(1,k,i,5);
            c2    = (double)gamma_ * q(1,k,i,5) / q(1,k,i,1);
            if ((float)bcdata(kk,ii,ip,1) > 0.0f) {
                c2 = bcdata(kk,ii,ip,1);
            } else if ((float)bcdata(kk,ii,ip,1) < 0.0f) {
                c2 = 1.0 + (double)gm1 * 0.5 * (double)xmach * (double)xmach;
            } else {
                if (iuns > 0 && irelv > 0) {
                    xm2 = (q(1,k,i,2)-xtbj(k,i,1,1))*(q(1,k,i,2)-xtbj(k,i,1,1))
                        + (q(1,k,i,3)-xtbj(k,i,2,1))*(q(1,k,i,3)-xtbj(k,i,2,1))
                        + (q(1,k,i,4)-xtbj(k,i,3,1))*(q(1,k,i,4)-xtbj(k,i,3,1));
                } else {
                    xm2 = q(1,k,i,2)*q(1,k,i,2)+q(1,k,i,3)*q(1,k,i,3)+q(1,k,i,4)*q(1,k,i,4);
                }
                xm2 = xm2 / c2;
                c2  = c2 * (1.0 + 0.5*(double)gm1*xm2);
            }
            uub = 0.0; vvb = 0.0; wwb = 0.0;
            if (iuns > 0 && irelv > 0) {
                uub = xtbj(k,i,1,1); vvb = xtbj(k,i,2,1); wwb = xtbj(k,i,3,1);
            }
            cqt = cq + cqu * std::sin(twopi * rfreq * (double)time_);
            if ((float)rfreq < 0.0f) {
                tcqt = cqu * (double)time_;
                if (std::abs(tcqt) < std::abs(cq)) { cqt = tcqt; } else { cqt = cq; }
            }
            sjet = std::sqrt(sjetx*sjetx + sjety*sjety + sjetz*sjetz);
            if (sjet != 0.0) {
                sjetx = sjetx/sjet; sjety = sjety/sjet; sjetz = sjetz/sjet;
            } else {
                sjetx = sj(1,k,i,1); sjety = sj(1,k,i,2); sjetz = sj(1,k,i,3);
            }
            uub = uub + sjetx*(double)xmach*cqt*c2/((double)gamma_*pb);
            vvb = vvb + sjety*(double)xmach*cqt*c2/((double)gamma_*pb);
            wwb = wwb + sjetz*(double)xmach*cqt*c2/((double)gamma_*pb);
            qj0(k,i,1,1) = (double)gamma_*pb/c2;
            qj0(k,i,2,1) = uub; qj0(k,i,3,1) = vvb; qj0(k,i,4,1) = wwb;
            qj0(k,i,5,1) = pb;
            bcj(k,i,1)   = 1.0;
            f23 = 0.0;
            j2  = std::min(2, jdim1);
            if (j2 == 1) f23 = 0.0;
            z1 = 2.0+1.5*f23; z2 = -0.5*f23; z3 = -(2.0+f23);
            qj0(k,i,1,2) = z1*q(1,k,i,1)+z2*q(j2,k,i,1)+z3*qj0(k,i,1,1);
            qj0(k,i,2,2) = z1*q(1,k,i,2)+z2*q(j2,k,i,2)+z3*qj0(k,i,2,1);
            qj0(k,i,3,2) = z1*q(1,k,i,3)+z2*q(j2,k,i,3)+z3*qj0(k,i,3,1);
            qj0(k,i,4,2) = z1*q(1,k,i,4)+z2*q(j2,k,i,4)+z3*qj0(k,i,4,1);
            qj0(k,i,5,2) = z1*q(1,k,i,5)+z2*q(j2,k,i,5)+z3*qj0(k,i,5,1);
        } // k loop
    } // i loop


    if (ivmx >= 2) {
    if (level >= lglobal && ntime != 0) {
        if (cmn_wallfun.iwf[1] == 0) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                vj0(k,i,1,1) = 0.0; vj0(k,i,1,2) = 0.0;
            }
        } else {
            // ***Wall Function begin
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                j = 1;
                c2b  = (double)cbar/(double)tinf; c2bp = c2b+1.0;
                uu = std::sqrt((q(j,k,i,2)-qj0(k,i,2,1))*(q(j,k,i,2)-qj0(k,i,2,1))
                             + (q(j,k,i,3)-qj0(k,i,3,1))*(q(j,k,i,3)-qj0(k,i,3,1))
                             + (q(j,k,i,4)-qj0(k,i,4,1))*(q(j,k,i,4)-qj0(k,i,4,1)));
                if (ivmx == 2) { dist = snj0(j,k,i); } else { dist = std::abs(smin(j,k,i)); }
                tt   = (double)gamma_*qj0(k,i,5,1)/qj0(k,i,1,1);
                fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                if ((float)rc <= 20.24f) {
                    xnplus = std::sqrt(rc);
                } else if ((float)rc <= 435.0f) {
                    xnplus = a0[1]+a0[2]*rc+(a0[3]*rc)*rc+(a0[4]*rc)*rc*rc
                           +(a0[5]*rc)*rc*rc*rc+(a0[6]*rc)*rc*rc*rc*rc+(a0[7]*rc)*rc*rc*rc*rc*rc;
                } else if ((float)rc <= 4000.0f) {
                    xnplus = a1[1]+a1[2]*rc+(a1[3]*rc)*rc+(a1[4]*rc)*rc*rc+(a1[5]*rc)*rc*rc*rc;
                } else {
                    xnplus = a2[1]+a2[2]*rc+a2[3]*rc*rc;
                }
                xnplussav = xnplus;
                if ((float)xnplus > 10.0f) {
                    bool conv1190 = false;
                    for (num = 1; num <= 10; num++) {
                        f_val = rc/xnplus - 2.44*std::log(xnplus) - 5.2;
                        dfdn  = -rc/(xnplus*xnplus) - 2.44/xnplus;
                        delta = -f_val/dfdn;
                        xnplus = std::abs(xnplus+delta);
                        if (std::abs((float)delta) < 1.e-3f) { conv1190 = true; break; }
                    }
                    if (!conv1190) xnplus = xnplussav;
                }
                dudy = uu/dist;
                xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/(q(j,k,i,1)*dudy)*(double)xmach/(double)reue-1.0);
                vj0(k,i,1,1) = xmut; vj0(k,i,1,2) = 0.0;
            }
            // ***Wall Function end
        }
    }
    } // ivmx >= 2

    // only need to do advanced model turbulence B.C.s on finest grid
    if (level >= lglobal && ntime != 0) {
    if (ivmx == 4 || ivmx == 5) {
        for (i = ista; i <= iend1; i++)
        for (k = ksta; k <= kend1; k++) {
            tj0(k,i,1,1) = -tursav(1,k,i,1); tj0(k,i,2,1) = -tursav(1,k,i,2);
            tj0(k,i,1,2) = 2.0*tj0(k,i,1,1)-tursav(1,k,i,1);
            tj0(k,i,2,2) = 2.0*tj0(k,i,2,1)-tursav(1,k,i,2);
        }
    }
    if (ivmx >= 6) {
        c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
        re  = (double)reue/(double)xmach;
        if (ivmx == 6) { beta1 = 0.075; }
        else if (ivmx == 7 || ivmx == 30 || ivmx == 40) { beta1 = 0.075; }
        else if (ivmx == 8 || ivmx == 12 || ivmx == 14) { beta1 = 0.83; }
        else { beta1 = 0.075; }
        j = 1;
        if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                tt   = (double)gamma_*qj0(k,i,5,1)/qj0(k,i,1,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                tj0(k,i,1,1) = 2.0*(2.0*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)-tursav(j,k,i,1);
                tj0(k,i,2,1) = -tursav(j,k,i,2);
                tj0(k,i,1,2) = 2.0*tj0(k,i,1,1)-tursav(1,k,i,1);
                tj0(k,i,2,2) = 2.0*tj0(k,i,2,1)-tursav(1,k,i,2);
            }
        } else if (ivmx == 15) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = 0.5*(tursav(j,k,i,2)+tursav(j+1,k,i,2))/3.0;
                tj0(k,i,1,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i))-tursav(j,k,i,1);
                tj0(k,i,2,1) = -tursav(j,k,i,2);
                tj0(k,i,1,2) = 0.0; tj0(k,i,2,2) = 0.0;
            }
        } else if (ivmx == 16) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                tj0(k,i,1,1) = -tursav(j,k,i,1); tj0(k,i,2,1) = -tursav(j,k,i,2);
                tj0(k,i,1,2) = 0.0; tj0(k,i,2,2) = 0.0;
            }
        } else {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                tt   = (double)gamma_*qj0(k,i,5,1)/qj0(k,i,1,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dist = std::abs(smin(j,k,i));
                tj0(k,i,1,1) = 2.0*(60.0*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))-tursav(j,k,i,1);
                tj0(k,i,2,1) = -tursav(j,k,i,2);
                tj0(k,i,1,2) = 2.0*tj0(k,i,1,1)-tursav(1,k,i,1);
                tj0(k,i,2,2) = 2.0*tj0(k,i,2,1)-tursav(1,k,i,2);
            }
        }
        // for trans model, intermittency at wall is either zero flux or set by user
        if (ivmx == 30) {
            for (i = ista; i <= iend1; i++) { ii = i-ista+1;
            for (k = ksta; k <= kend1; k++) { kk = k-ksta+1;
                if (iuse3 == 0) {
                    tj0(k,i,3,1) = tursav(j,k,i,3);
                } else {
                    tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                    fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                    dist = std::abs(smin(j,k,i));
                    hee  = (1.0+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                    rhoee = std::pow((double)gamma_*q(j,k,i,5), 1.0/(double)gamma_);
                    uee2 = ((double)gm1*rhoee*hee-(double)gamma_*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                    if ((float)uee2 > 0.0f) { uee = std::sqrt(uee2); } else { uee = (double)xmach; }
                    tj0(k,i,3,1) = tursav(j,k,i,3)+2.0*dist*q(j,k,i,1)/fnu
                                  *uee*bcdata(kk,ii,ip,3)*(double)reue/(double)xmach;
                }
                tj0(k,i,3,2) = 2.0*tj0(k,i,3,1)-tursav(1,k,i,3);
            }}
        } else if (ivmx == 40) {
            for (i = ista; i <= iend1; i++) { ii = i-ista+1;
            for (k = ksta; k <= kend1; k++) { kk = k-ksta+1;
                tj0(k,i,3,1) = tursav(j,k,i,3); tj0(k,i,4,1) = tursav(j,k,i,4);
            }}
        }
    } // ivmx >= 6
    } // level >= lglobal && ntime != 0

    } // nface == 3



    // **************************************************************************
    // j=jdim boundary   viscous wall with T & cq specified   bctype 2016
    // **************************************************************************
    if (nface == 4) {

    for (i = ista; i <= iend1; i++) {
        ii = i - ista + 1;
        for (k = ksta; k <= kend1; k++) {
            kk = k - ksta + 1;
            cq    = bcdata(kk,ii,ip,2); cqu   = bcdata(kk,ii,ip,3);
            sjetx = bcdata(kk,ii,ip,4); sjety = bcdata(kk,ii,ip,5);
            sjetz = bcdata(kk,ii,ip,6); rfreq = bcdata(kk,ii,ip,7);
            pb    = q(jdim1,k,i,5);
            dpb   = q(jdim1,k,i,5) - q(jdim1-1,k,i,5);
            pb    = pb + dpb/2.0;
            if ((float)pb <= 0.0f) pb = q(jdim1,k,i,5);
            c2    = (double)gamma_*q(jdim1,k,i,5)/q(jdim1,k,i,1);
            if ((float)bcdata(kk,ii,ip,1) > 0.0f) {
                c2 = bcdata(kk,ii,ip,1);
            } else if ((float)bcdata(kk,ii,ip,1) < 0.0f) {
                c2 = 1.0+(double)gm1*0.5*(double)xmach*(double)xmach;
            } else {
                if (iuns > 0 && irelv > 0) {
                    xm2 = (q(jdim1,k,i,2)-xtbj(k,i,1,2))*(q(jdim1,k,i,2)-xtbj(k,i,1,2))
                        + (q(jdim1,k,i,3)-xtbj(k,i,2,2))*(q(jdim1,k,i,3)-xtbj(k,i,2,2))
                        + (q(jdim1,k,i,4)-xtbj(k,i,3,2))*(q(jdim1,k,i,4)-xtbj(k,i,3,2));
                } else {
                    xm2 = q(jdim1,k,i,2)*q(jdim1,k,i,2)+q(jdim1,k,i,3)*q(jdim1,k,i,3)
                         +q(jdim1,k,i,4)*q(jdim1,k,i,4);
                }
                xm2 = xm2/c2; c2 = c2*(1.0+0.5*(double)gm1*xm2);
            }
            uub = 0.0; vvb = 0.0; wwb = 0.0;
            if (iuns > 0 && irelv > 0) {
                uub = xtbj(k,i,1,2); vvb = xtbj(k,i,2,2); wwb = xtbj(k,i,3,2);
            }
            cqt = cq + cqu*std::sin(twopi*rfreq*(double)time_);
            if ((float)rfreq < 0.0f) {
                tcqt = cqu*(double)time_;
                if (std::abs(tcqt) < std::abs(cq)) { cqt = tcqt; } else { cqt = cq; }
            }
            sjet = std::sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz);
            if (sjet != 0.0) {
                sjetx = sjetx/sjet; sjety = sjety/sjet; sjetz = sjetz/sjet;
            } else {
                sjetx = -sj(jdim,k,i,1); sjety = -sj(jdim,k,i,2); sjetz = -sj(jdim,k,i,3);
            }
            uub = uub + sjetx*(double)xmach*cqt*c2/((double)gamma_*pb);
            vvb = vvb + sjety*(double)xmach*cqt*c2/((double)gamma_*pb);
            wwb = wwb + sjetz*(double)xmach*cqt*c2/((double)gamma_*pb);
            qj0(k,i,1,3) = (double)gamma_*pb/c2;
            qj0(k,i,2,3) = uub; qj0(k,i,3,3) = vvb; qj0(k,i,4,3) = wwb;
            qj0(k,i,5,3) = pb;
            bcj(k,i,2)   = 1.0;
            f23 = 0.0;
            j2  = std::max(1, jdim-2);
            if (j2 == 1) f23 = 0.0;
            z1 = -2.0-1.5*f23; z2 = 0.5*f23; z3 = 2.0+f23;
            qj0(k,i,1,4) = z1*q(jdim1,k,i,1)+z2*q(j2,k,i,1)+z3*qj0(k,i,1,3);
            qj0(k,i,2,4) = z1*q(jdim1,k,i,2)+z2*q(j2,k,i,2)+z3*qj0(k,i,2,3);
            qj0(k,i,3,4) = z1*q(jdim1,k,i,3)+z2*q(j2,k,i,3)+z3*qj0(k,i,3,3);
            qj0(k,i,4,4) = z1*q(jdim1,k,i,4)+z2*q(j2,k,i,4)+z3*qj0(k,i,4,3);
            qj0(k,i,5,4) = z1*q(jdim1,k,i,5)+z2*q(j2,k,i,5)+z3*qj0(k,i,5,3);
        }
    }

    if (ivmx >= 2) {
    if (level >= lglobal && ntime != 0) {
        if (cmn_wallfun.iwf[1] == 0) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                vj0(k,i,1,3) = 0.0; vj0(k,i,1,4) = 0.0;
            }
        } else {
            // ***Wall Function begin
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                j = jdim-1;
                c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
                uu = std::sqrt((q(j,k,i,2)-qj0(k,i,2,3))*(q(j,k,i,2)-qj0(k,i,2,3))
                             + (q(j,k,i,3)-qj0(k,i,3,3))*(q(j,k,i,3)-qj0(k,i,3,3))
                             + (q(j,k,i,4)-qj0(k,i,4,3))*(q(j,k,i,4)-qj0(k,i,4,3)));
                if (ivmx == 2) { dist = snjm(j,k,i); } else { dist = std::abs(smin(j,k,i)); }
                tt   = (double)gamma_*qj0(k,i,5,3)/qj0(k,i,1,3);
                fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                if ((float)rc <= 20.24f) {
                    xnplus = std::sqrt(rc);
                } else if ((float)rc <= 435.0f) {
                    xnplus = a0[1]+a0[2]*rc+(a0[3]*rc)*rc+(a0[4]*rc)*rc*rc
                           +(a0[5]*rc)*rc*rc*rc+(a0[6]*rc)*rc*rc*rc*rc+(a0[7]*rc)*rc*rc*rc*rc*rc;
                } else if ((float)rc <= 4000.0f) {
                    xnplus = a1[1]+a1[2]*rc+(a1[3]*rc)*rc+(a1[4]*rc)*rc*rc+(a1[5]*rc)*rc*rc*rc;
                } else {
                    xnplus = a2[1]+a2[2]*rc+a2[3]*rc*rc;
                }
                xnplussav = xnplus;
                if ((float)xnplus > 10.0f) {
                    bool conv1290 = false;
                    for (num = 1; num <= 10; num++) {
                        f_val = rc/xnplus-2.44*std::log(xnplus)-5.2;
                        dfdn  = -rc/(xnplus*xnplus)-2.44/xnplus;
                        delta = -f_val/dfdn;
                        xnplus = std::abs(xnplus+delta);
                        if (std::abs((float)delta) < 1.e-3f) { conv1290 = true; break; }
                    }
                    if (!conv1290) xnplus = xnplussav;
                }
                dudy = uu/dist;
                xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/(q(j,k,i,1)*dudy)*(double)xmach/(double)reue-1.0);
                vj0(k,i,1,3) = xmut; vj0(k,i,1,4) = 0.0;
            }
            // ***Wall Function end
        }
    }
    } // ivmx >= 2


    // only need to do advanced model turbulence B.C.s on finest grid (nface==4)
    if (level >= lglobal && ntime != 0) {
    if (ivmx == 4 || ivmx == 5) {
        for (i = ista; i <= iend1; i++)
        for (k = ksta; k <= kend1; k++) {
            tj0(k,i,1,3) = -tursav(jdim-1,k,i,1); tj0(k,i,2,3) = -tursav(jdim-1,k,i,2);
            tj0(k,i,1,4) = 2.0*tj0(k,i,1,3)-tursav(jdim-1,k,i,1);
            tj0(k,i,2,4) = 2.0*tj0(k,i,2,3)-tursav(jdim-1,k,i,2);
        }
    }
    if (ivmx >= 6) {
        c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
        re  = (double)reue/(double)xmach;
        if (ivmx == 6) { beta1 = 0.075; }
        else if (ivmx == 7 || ivmx == 30 || ivmx == 40) { beta1 = 0.075; }
        else if (ivmx == 8 || ivmx == 12 || ivmx == 14) { beta1 = 0.83; }
        else { beta1 = 0.075; }
        j = jdim-1;
        if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                tt   = (double)gamma_*qj0(k,i,5,3)/qj0(k,i,1,3);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                tj0(k,i,1,3) = 2.0*(2.0*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)-tursav(j,k,i,1);
                tj0(k,i,2,3) = -tursav(j,k,i,2);
                tj0(k,i,1,4) = 2.0*tj0(k,i,1,3)-tursav(jdim-1,k,i,1);
                tj0(k,i,2,4) = 2.0*tj0(k,i,2,3)-tursav(jdim-1,k,i,2);
            }
        } else if (ivmx == 15) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = 0.5*(tursav(j,k,i,2)+tursav(j-1,k,i,2))/3.0;
                tj0(k,i,1,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i))-tursav(j,k,i,1);
                tj0(k,i,2,3) = -tursav(j,k,i,2);
                tj0(k,i,1,4) = 0.0; tj0(k,i,2,4) = 0.0;
            }
        } else if (ivmx == 16) {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                tj0(k,i,1,3) = -tursav(j,k,i,1); tj0(k,i,2,3) = -tursav(j,k,i,2);
                tj0(k,i,1,4) = 0.0; tj0(k,i,2,4) = 0.0;
            }
        } else {
            for (i = ista; i <= iend1; i++)
            for (k = ksta; k <= kend1; k++) {
                tt   = (double)gamma_*qj0(k,i,5,3)/qj0(k,i,1,3);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dist = std::abs(smin(j,k,i));
                tj0(k,i,1,3) = 2.0*(60.0*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))-tursav(j,k,i,1);
                tj0(k,i,2,3) = -tursav(j,k,i,2);
                tj0(k,i,1,4) = 2.0*tj0(k,i,1,3)-tursav(jdim-1,k,i,1);
                tj0(k,i,2,4) = 2.0*tj0(k,i,2,3)-tursav(jdim-1,k,i,2);
            }
        }
        if (ivmx == 30) {
            for (i = ista; i <= iend1; i++) { ii = i-ista+1;
            for (k = ksta; k <= kend1; k++) { kk = k-ksta+1;
                if (iuse3 == 0) {
                    tj0(k,i,3,3) = tursav(j,k,i,3);
                } else {
                    tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                    fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                    dist = std::abs(smin(j,k,i));
                    hee  = (1.0+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                    rhoee = std::pow((double)gamma_*q(j,k,i,5), 1.0/(double)gamma_);
                    uee2 = ((double)gm1*rhoee*hee-(double)gamma_*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                    if ((float)uee2 > 0.0f) { uee = std::sqrt(uee2); } else { uee = (double)xmach; }
                    tj0(k,i,3,3) = tursav(j,k,i,3)+2.0*dist*q(j,k,i,1)/fnu
                                  *uee*bcdata(kk,ii,ip,3)*(double)reue/(double)xmach;
                }
                tj0(k,i,3,4) = 2.0*tj0(k,i,3,3)-tursav(jdim-1,k,i,3);
            }}
        } else if (ivmx == 40) {
            for (i = ista; i <= iend1; i++) { ii = i-ista+1;
            for (k = ksta; k <= kend1; k++) { kk = k-ksta+1;
                tj0(k,i,3,3) = tursav(j,k,i,3); tj0(k,i,4,3) = tursav(j,k,i,4);
            }}
        }
    } // ivmx >= 6
    } // level >= lglobal && ntime != 0

    } // nface == 4



    // **************************************************************************
    // k=1 boundary   viscous wall with T & cq specified   bctype 2016
    // **************************************************************************
    if (nface == 5) {

    for (i = ista; i <= iend1; i++) {
        ii = i - ista + 1;
        for (j = jsta; j <= jend1; j++) {
            jj = j - jsta + 1;
            cq    = bcdata(jj,ii,ip,2); cqu   = bcdata(jj,ii,ip,3);
            sjetx = bcdata(jj,ii,ip,4); sjety = bcdata(jj,ii,ip,5);
            sjetz = bcdata(jj,ii,ip,6); rfreq = bcdata(jj,ii,ip,7);
            pb    = q(j,1,i,5);
            dpb   = q(j,2,i,5) - q(j,1,i,5);
            pb    = pb - dpb/2.0;
            if ((float)pb <= 0.0f) pb = q(j,1,i,5);
            c2    = (double)gamma_*q(j,1,i,5)/q(j,1,i,1);
            if ((float)bcdata(jj,ii,ip,1) > 0.0f) {
                c2 = bcdata(jj,ii,ip,1);
            } else if ((float)bcdata(jj,ii,ip,1) < 0.0f) {
                c2 = 1.0+(double)gm1*0.5*(double)xmach*(double)xmach;
            } else {
                if (iuns > 0 && irelv > 0) {
                    xm2 = (q(j,1,i,2)-xtbk(j,i,1,1))*(q(j,1,i,2)-xtbk(j,i,1,1))
                        + (q(j,1,i,3)-xtbk(j,i,2,1))*(q(j,1,i,3)-xtbk(j,i,2,1))
                        + (q(j,1,i,4)-xtbk(j,i,3,1))*(q(j,1,i,4)-xtbk(j,i,3,1));
                } else {
                    xm2 = q(j,1,i,2)*q(j,1,i,2)+q(j,1,i,3)*q(j,1,i,3)+q(j,1,i,4)*q(j,1,i,4);
                }
                xm2 = xm2/c2; c2 = c2*(1.0+0.5*(double)gm1*xm2);
            }
            uub = 0.0; vvb = 0.0; wwb = 0.0;
            if (iuns > 0 && irelv > 0) {
                uub = xtbk(j,i,1,1); vvb = xtbk(j,i,2,1); wwb = xtbk(j,i,3,1);
            }
            cqt = cq + cqu*std::sin(twopi*rfreq*(double)time_);
            if ((float)rfreq < 0.0f) {
                tcqt = cqu*(double)time_;
                if (std::abs(tcqt) < std::abs(cq)) { cqt = tcqt; } else { cqt = cq; }
            }
            sjet = std::sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz);
            if (sjet != 0.0) {
                sjetx = sjetx/sjet; sjety = sjety/sjet; sjetz = sjetz/sjet;
            } else {
                sjetx = sk(j,1,i,1); sjety = sk(j,1,i,2); sjetz = sk(j,1,i,3);
            }
            uub = uub + sjetx*(double)xmach*cqt*c2/((double)gamma_*pb);
            vvb = vvb + sjety*(double)xmach*cqt*c2/((double)gamma_*pb);
            wwb = wwb + sjetz*(double)xmach*cqt*c2/((double)gamma_*pb);
            qk0(j,i,1,1) = (double)gamma_*pb/c2;
            qk0(j,i,2,1) = uub; qk0(j,i,3,1) = vvb; qk0(j,i,4,1) = wwb;
            qk0(j,i,5,1) = pb;
            bck(j,i,1)   = 1.0;
            f23 = 0.0;
            k2  = std::min(2, kdim1);
            if (k2 == 1) f23 = 0.0;
            z1 = 2.0+1.5*f23; z2 = -0.5*f23; z3 = -(2.0+f23);
            qk0(j,i,1,2) = z1*q(j,1,i,1)+z2*q(j,k2,i,1)+z3*qk0(j,i,1,1);
            qk0(j,i,2,2) = z1*q(j,1,i,2)+z2*q(j,k2,i,2)+z3*qk0(j,i,2,1);
            qk0(j,i,3,2) = z1*q(j,1,i,3)+z2*q(j,k2,i,3)+z3*qk0(j,i,3,1);
            qk0(j,i,4,2) = z1*q(j,1,i,4)+z2*q(j,k2,i,4)+z3*qk0(j,i,4,1);
            qk0(j,i,5,2) = z1*q(j,1,i,5)+z2*q(j,k2,i,5)+z3*qk0(j,i,5,1);
        }
    }

    if (ivmx >= 2) {
    if (level >= lglobal && ntime != 0) {
        if (cmn_wallfun.iwf[2] == 0) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                vk0(j,i,1,1) = 0.0; vk0(j,i,1,2) = 0.0;
            }
        } else {
            // ***Wall Function begin
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                k = 1;
                c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
                uu = std::sqrt((q(j,k,i,2)-qk0(j,i,2,1))*(q(j,k,i,2)-qk0(j,i,2,1))
                             + (q(j,k,i,3)-qk0(j,i,3,1))*(q(j,k,i,3)-qk0(j,i,3,1))
                             + (q(j,k,i,4)-qk0(j,i,4,1))*(q(j,k,i,4)-qk0(j,i,4,1)));
                if (ivmx == 2) { dist = snk0(j,k,i); } else { dist = std::abs(smin(j,k,i)); }
                tt   = (double)gamma_*qk0(j,i,5,1)/qk0(j,i,1,1);
                fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                if ((float)rc <= 20.24f) {
                    xnplus = std::sqrt(rc);
                } else if ((float)rc <= 435.0f) {
                    xnplus = a0[1]+a0[2]*rc+(a0[3]*rc)*rc+(a0[4]*rc)*rc*rc
                           +(a0[5]*rc)*rc*rc*rc+(a0[6]*rc)*rc*rc*rc*rc+(a0[7]*rc)*rc*rc*rc*rc*rc;
                } else if ((float)rc <= 4000.0f) {
                    xnplus = a1[1]+a1[2]*rc+(a1[3]*rc)*rc+(a1[4]*rc)*rc*rc+(a1[5]*rc)*rc*rc*rc;
                } else {
                    xnplus = a2[1]+a2[2]*rc+a2[3]*rc*rc;
                }
                xnplussav = xnplus;
                if ((float)xnplus > 10.0f) {
                    bool conv1390 = false;
                    for (num = 1; num <= 10; num++) {
                        f_val = rc/xnplus-2.44*std::log(xnplus)-5.2;
                        dfdn  = -rc/(xnplus*xnplus)-2.44/xnplus;
                        delta = -f_val/dfdn;
                        xnplus = std::abs(xnplus+delta);
                        if (std::abs((float)delta) < 1.e-3f) { conv1390 = true; break; }
                    }
                    if (!conv1390) xnplus = xnplussav;
                }
                dudy = uu/dist;
                xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/(q(j,k,i,1)*dudy)*(double)xmach/(double)reue-1.0);
                vk0(j,i,1,1) = xmut; vk0(j,i,1,2) = 0.0;
            }
            // ***Wall Function end
        }
    }
    } // ivmx >= 2


    // only need to do advanced model turbulence B.C.s on finest grid (nface==5)
    if (level >= lglobal && ntime != 0) {
    if (ivmx == 4 || ivmx == 5) {
        for (i = ista; i <= iend1; i++)
        for (j = jsta; j <= jend1; j++) {
            tk0(j,i,1,1) = -tursav(j,1,i,1); tk0(j,i,2,1) = -tursav(j,1,i,2);
            tk0(j,i,1,2) = 2.0*tk0(j,i,1,1)-tursav(j,1,i,1);
            tk0(j,i,2,2) = 2.0*tk0(j,i,2,1)-tursav(j,1,i,2);
        }
    }
    if (ivmx >= 6) {
        c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
        re  = (double)reue/(double)xmach;
        if (ivmx == 6) { beta1 = 0.075; }
        else if (ivmx == 7 || ivmx == 30 || ivmx == 40) { beta1 = 0.075; }
        else if (ivmx == 8 || ivmx == 12 || ivmx == 14) { beta1 = 0.83; }
        else { beta1 = 0.075; }
        k = 1;
        if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*qk0(j,i,5,1)/qk0(j,i,1,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                tk0(j,i,1,1) = 2.0*(2.0*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)-tursav(j,k,i,1);
                tk0(j,i,2,1) = -tursav(j,k,i,2);
                tk0(j,i,1,2) = 2.0*tk0(j,i,1,1)-tursav(j,1,i,1);
                tk0(j,i,2,2) = 2.0*tk0(j,i,2,1)-tursav(j,1,i,2);
            }
        } else if (ivmx == 15) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = 0.5*(tursav(j,k,i,2)+tursav(j,k+1,i,2))/3.0;
                tk0(j,i,1,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i))-tursav(j,k,i,1);
                tk0(j,i,2,1) = -tursav(j,k,i,2);
                tk0(j,i,1,2) = 0.0; tk0(j,i,2,2) = 0.0;
            }
        } else if (ivmx == 16) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                tk0(j,i,1,1) = -tursav(j,k,i,1); tk0(j,i,2,1) = -tursav(j,k,i,2);
                tk0(j,i,1,2) = 0.0; tk0(j,i,2,2) = 0.0;
            }
        } else {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*qk0(j,i,5,1)/qk0(j,i,1,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dist = std::abs(smin(j,k,i));
                tk0(j,i,1,1) = 2.0*(60.0*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))-tursav(j,k,i,1);
                tk0(j,i,2,1) = -tursav(j,k,i,2);
                tk0(j,i,1,2) = 2.0*tk0(j,i,1,1)-tursav(j,1,i,1);
                tk0(j,i,2,2) = 2.0*tk0(j,i,2,1)-tursav(j,1,i,2);
            }
        }
        if (ivmx == 30) {
            for (i = ista; i <= iend1; i++) { ii = i-ista+1;
            for (j = jsta; j <= jend1; j++) { jj = j-jsta+1;
                if (iuse3 == 0) {
                    tk0(j,i,3,1) = tursav(j,k,i,3);
                } else {
                    tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                    fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                    dist = std::abs(smin(j,k,i));
                    hee  = (1.0+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                    rhoee = std::pow((double)gamma_*q(j,k,i,5), 1.0/(double)gamma_);
                    uee2 = ((double)gm1*rhoee*hee-(double)gamma_*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                    if ((float)uee2 > 0.0f) { uee = std::sqrt(uee2); } else { uee = (double)xmach; }
                    tk0(j,i,3,1) = tursav(j,k,i,3)+2.0*dist*q(j,k,i,1)/fnu
                                  *uee*bcdata(jj,ii,ip,3)*(double)reue/(double)xmach;
                }
                tk0(j,i,3,2) = 2.0*tk0(j,i,3,1)-tursav(j,1,i,3);
            }}
        } else if (ivmx == 40) {
            for (i = ista; i <= iend1; i++) { ii = i-ista+1;
            for (j = jsta; j <= jend1; j++) { jj = j-jsta+1;
                tk0(j,i,3,1) = tursav(j,k,i,3); tk0(j,i,4,1) = tursav(j,k,i,4);
            }}
        }
    } // ivmx >= 6
    } // level >= lglobal && ntime != 0

    } // nface == 5



    // **************************************************************************
    // k=kdim boundary   viscous wall with T & cq specified   bctype 2016
    // **************************************************************************
    if (nface == 6) {

    for (i = ista; i <= iend1; i++) {
        ii = i - ista + 1;
        for (j = jsta; j <= jend1; j++) {
            jj = j - jsta + 1;
            cq    = bcdata(jj,ii,ip,2); cqu   = bcdata(jj,ii,ip,3);
            sjetx = bcdata(jj,ii,ip,4); sjety = bcdata(jj,ii,ip,5);
            sjetz = bcdata(jj,ii,ip,6); rfreq = bcdata(jj,ii,ip,7);
            pb    = q(j,kdim1,i,5);
            dpb   = q(j,kdim1,i,5) - q(j,kdim1-1,i,5);
            pb    = pb + dpb/2.0;
            if ((float)pb <= 0.0f) pb = q(j,kdim1,i,5);
            c2    = (double)gamma_*q(j,kdim1,i,5)/q(j,kdim1,i,1);
            if ((float)bcdata(jj,ii,ip,1) > 0.0f) {
                c2 = bcdata(jj,ii,ip,1);
            } else if ((float)bcdata(jj,ii,ip,1) < 0.0f) {
                c2 = 1.0+(double)gm1*0.5*(double)xmach*(double)xmach;
            } else {
                if (iuns > 0 && irelv > 0) {
                    xm2 = (q(j,kdim1,i,2)-xtbk(j,i,1,2))*(q(j,kdim1,i,2)-xtbk(j,i,1,2))
                        + (q(j,kdim1,i,3)-xtbk(j,i,2,2))*(q(j,kdim1,i,3)-xtbk(j,i,2,2))
                        + (q(j,kdim1,i,4)-xtbk(j,i,3,2))*(q(j,kdim1,i,4)-xtbk(j,i,3,2));
                } else {
                    xm2 = q(j,kdim1,i,2)*q(j,kdim1,i,2)+q(j,kdim1,i,3)*q(j,kdim1,i,3)
                         +q(j,kdim1,i,4)*q(j,kdim1,i,4);
                }
                xm2 = xm2/c2; c2 = c2*(1.0+0.5*(double)gm1*xm2);
            }
            uub = 0.0; vvb = 0.0; wwb = 0.0;
            if (iuns > 0 && irelv > 0) {
                uub = xtbk(j,i,1,2); vvb = xtbk(j,i,2,2); wwb = xtbk(j,i,3,2);
            }
            cqt = cq + cqu*std::sin(twopi*rfreq*(double)time_);
            if ((float)rfreq < 0.0f) {
                tcqt = cqu*(double)time_;
                if (std::abs(tcqt) < std::abs(cq)) { cqt = tcqt; } else { cqt = cq; }
            }
            sjet = std::sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz);
            if (sjet != 0.0) {
                sjetx = sjetx/sjet; sjety = sjety/sjet; sjetz = sjetz/sjet;
            } else {
                sjetx = -sk(j,kdim,i,1); sjety = -sk(j,kdim,i,2); sjetz = -sk(j,kdim,i,3);
            }
            uub = uub + sjetx*(double)xmach*cqt*c2/((double)gamma_*pb);
            vvb = vvb + sjety*(double)xmach*cqt*c2/((double)gamma_*pb);
            wwb = wwb + sjetz*(double)xmach*cqt*c2/((double)gamma_*pb);
            qk0(j,i,1,3) = (double)gamma_*pb/c2;
            qk0(j,i,2,3) = uub; qk0(j,i,3,3) = vvb; qk0(j,i,4,3) = wwb;
            qk0(j,i,5,3) = pb;
            bck(j,i,2)   = 1.0;
            f23 = 0.0;
            k2  = std::max(1, kdim-2);
            if (k2 == 1) f23 = 0.0;
            z1 = -2.0-1.5*f23; z2 = 0.5*f23; z3 = 2.0+f23;
            qk0(j,i,1,4) = z1*q(j,kdim1,i,1)+z2*q(j,k2,i,1)+z3*qk0(j,i,1,3);
            qk0(j,i,2,4) = z1*q(j,kdim1,i,2)+z2*q(j,k2,i,2)+z3*qk0(j,i,2,3);
            qk0(j,i,3,4) = z1*q(j,kdim1,i,3)+z2*q(j,k2,i,3)+z3*qk0(j,i,3,3);
            qk0(j,i,4,4) = z1*q(j,kdim1,i,4)+z2*q(j,k2,i,4)+z3*qk0(j,i,4,3);
            qk0(j,i,5,4) = z1*q(j,kdim1,i,5)+z2*q(j,k2,i,5)+z3*qk0(j,i,5,3);
        }
    }

    if (ivmx >= 2) {
    if (level >= lglobal && ntime != 0) {
        if (cmn_wallfun.iwf[2] == 0) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                vk0(j,i,1,3) = 0.0; vk0(j,i,1,4) = 0.0;
            }
        } else {
            // ***Wall Function begin
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                k = kdim-1;
                c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
                uu = std::sqrt((q(j,k,i,2)-qk0(j,i,2,3))*(q(j,k,i,2)-qk0(j,i,2,3))
                             + (q(j,k,i,3)-qk0(j,i,3,3))*(q(j,k,i,3)-qk0(j,i,3,3))
                             + (q(j,k,i,4)-qk0(j,i,4,3))*(q(j,k,i,4)-qk0(j,i,4,3)));
                if (ivmx == 2) { dist = snkm(j,k,i); } else { dist = std::abs(smin(j,k,i)); }
                tt   = (double)gamma_*qk0(j,i,5,3)/qk0(j,i,1,3);
                fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                if ((float)rc <= 20.24f) {
                    xnplus = std::sqrt(rc);
                } else if ((float)rc <= 435.0f) {
                    xnplus = a0[1]+a0[2]*rc+(a0[3]*rc)*rc+(a0[4]*rc)*rc*rc
                           +(a0[5]*rc)*rc*rc*rc+(a0[6]*rc)*rc*rc*rc*rc+(a0[7]*rc)*rc*rc*rc*rc*rc;
                } else if ((float)rc <= 4000.0f) {
                    xnplus = a1[1]+a1[2]*rc+(a1[3]*rc)*rc+(a1[4]*rc)*rc*rc+(a1[5]*rc)*rc*rc*rc;
                } else {
                    xnplus = a2[1]+a2[2]*rc+a2[3]*rc*rc;
                }
                xnplussav = xnplus;
                if ((float)xnplus > 10.0f) {
                    bool conv1490 = false;
                    for (num = 1; num <= 10; num++) {
                        f_val = rc/xnplus-2.44*std::log(xnplus)-5.2;
                        dfdn  = -rc/(xnplus*xnplus)-2.44/xnplus;
                        delta = -f_val/dfdn;
                        xnplus = std::abs(xnplus+delta);
                        if (std::abs((float)delta) < 1.e-3f) { conv1490 = true; break; }
                    }
                    if (!conv1490) xnplus = xnplussav;
                }
                dudy = uu/dist;
                xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/(q(j,k,i,1)*dudy)*(double)xmach/(double)reue-1.0);
                vk0(j,i,1,3) = xmut; vk0(j,i,1,4) = 0.0;
            }
            // ***Wall Function end
        }
    }
    } // ivmx >= 2


    // only need to do advanced model turbulence B.C.s on finest grid (nface==6)
    if (level >= lglobal && ntime != 0) {
    if (ivmx == 4 || ivmx == 5) {
        for (i = ista; i <= iend1; i++)
        for (j = jsta; j <= jend1; j++) {
            tk0(j,i,1,3) = -tursav(j,kdim-1,i,1); tk0(j,i,2,3) = -tursav(j,kdim-1,i,2);
            tk0(j,i,1,4) = 2.0*tk0(j,i,1,3)-tursav(j,kdim-1,i,1);
            tk0(j,i,2,4) = 2.0*tk0(j,i,2,3)-tursav(j,kdim-1,i,2);
        }
    }
    if (ivmx >= 6) {
        c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
        re  = (double)reue/(double)xmach;
        if (ivmx == 6) { beta1 = 0.075; }
        else if (ivmx == 7 || ivmx == 30 || ivmx == 40) { beta1 = 0.075; }
        else if (ivmx == 8 || ivmx == 12 || ivmx == 14) { beta1 = 0.83; }
        else { beta1 = 0.075; }
        k = kdim-1;
        if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*qk0(j,i,5,3)/qk0(j,i,1,3);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                tk0(j,i,1,3) = 2.0*(2.0*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)-tursav(j,k,i,1);
                tk0(j,i,2,3) = -tursav(j,k,i,2);
                tk0(j,i,1,4) = 2.0*tk0(j,i,1,3)-tursav(j,kdim-1,i,1);
                tk0(j,i,2,4) = 2.0*tk0(j,i,2,3)-tursav(j,kdim-1,i,2);
            }
        } else if (ivmx == 15) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = 0.5*(tursav(j,k,i,2)+tursav(j,k-1,i,2))/3.0;
                tk0(j,i,1,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i))-tursav(j,k,i,1);
                tk0(j,i,2,3) = -tursav(j,k,i,2);
                tk0(j,i,1,4) = 0.0; tk0(j,i,2,4) = 0.0;
            }
        } else if (ivmx == 16) {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                tk0(j,i,1,3) = -tursav(j,k,i,1); tk0(j,i,2,3) = -tursav(j,k,i,2);
                tk0(j,i,1,4) = 0.0; tk0(j,i,2,4) = 0.0;
            }
        } else {
            for (i = ista; i <= iend1; i++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*qk0(j,i,5,3)/qk0(j,i,1,3);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dist = std::abs(smin(j,k,i));
                tk0(j,i,1,3) = 2.0*(60.0*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))-tursav(j,k,i,1);
                tk0(j,i,2,3) = -tursav(j,k,i,2);
                tk0(j,i,1,4) = 2.0*tk0(j,i,1,3)-tursav(j,kdim-1,i,1);
                tk0(j,i,2,4) = 2.0*tk0(j,i,2,3)-tursav(j,kdim-1,i,2);
            }
        }
        if (ivmx == 30) {
            for (i = ista; i <= iend1; i++) { ii = i-ista+1;
            for (j = jsta; j <= jend1; j++) { jj = j-jsta+1;
                if (iuse3 == 0) {
                    tk0(j,i,3,3) = tursav(j,k,i,3);
                } else {
                    tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                    fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                    dist = std::abs(smin(j,k,i));
                    hee  = (1.0+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                    rhoee = std::pow((double)gamma_*q(j,k,i,5), 1.0/(double)gamma_);
                    uee2 = ((double)gm1*rhoee*hee-(double)gamma_*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                    if ((float)uee2 > 0.0f) { uee = std::sqrt(uee2); } else { uee = (double)xmach; }
                    tk0(j,i,3,3) = tursav(j,k,i,3)+2.0*dist*q(j,k,i,1)/fnu
                                  *uee*bcdata(jj,ii,ip,3)*(double)reue/(double)xmach;
                }
                tk0(j,i,3,4) = 2.0*tk0(j,i,3,3)-tursav(j,kdim-1,i,3);
            }}
        } else if (ivmx == 40) {
            for (i = ista; i <= iend1; i++) { ii = i-ista+1;
            for (j = jsta; j <= jend1; j++) { jj = j-jsta+1;
                tk0(j,i,3,3) = tursav(j,k,i,3); tk0(j,i,4,3) = tursav(j,k,i,4);
            }}
        }
    } // ivmx >= 6
    } // level >= lglobal && ntime != 0

    } // nface == 6



    // **************************************************************************
    // i=1 boundary   viscous wall with T & cq specified   bctype 2016
    // **************************************************************************
    if (nface == 1) {

    for (k = ksta; k <= kend1; k++) {
        kk = k - ksta + 1;
        for (j = jsta; j <= jend1; j++) {
            jj = j - jsta + 1;
            cq    = bcdata(jj,kk,ip,2); cqu   = bcdata(jj,kk,ip,3);
            sjetx = bcdata(jj,kk,ip,4); sjety = bcdata(jj,kk,ip,5);
            sjetz = bcdata(jj,kk,ip,6); rfreq = bcdata(jj,kk,ip,7);
            pb    = q(j,k,1,5);
            dpb   = q(j,k,2,5) - q(j,k,1,5);
            pb    = pb - dpb/2.0;
            if ((float)pb <= 0.0f) pb = q(j,k,1,5);
            c2    = (double)gamma_*q(j,k,1,5)/q(j,k,1,1);
            if ((float)bcdata(jj,kk,ip,1) > 0.0f) {
                c2 = bcdata(jj,kk,ip,1);
            } else if ((float)bcdata(jj,kk,ip,1) < 0.0f) {
                c2 = 1.0+(double)gm1*0.5*(double)xmach*(double)xmach;
            } else {
                if (iuns > 0 && irelv > 0) {
                    xm2 = (q(j,k,1,2)-xtbi(j,k,1,1))*(q(j,k,1,2)-xtbi(j,k,1,1))
                        + (q(j,k,1,3)-xtbi(j,k,2,1))*(q(j,k,1,3)-xtbi(j,k,2,1))
                        + (q(j,k,1,4)-xtbi(j,k,3,1))*(q(j,k,1,4)-xtbi(j,k,3,1));
                } else {
                    xm2 = q(j,k,1,2)*q(j,k,1,2)+q(j,k,1,3)*q(j,k,1,3)+q(j,k,1,4)*q(j,k,1,4);
                }
                xm2 = xm2/c2; c2 = c2*(1.0+0.5*(double)gm1*xm2);
            }
            uub = 0.0; vvb = 0.0; wwb = 0.0;
            if (iuns > 0 && irelv > 0) {
                uub = xtbi(j,k,1,1); vvb = xtbi(j,k,2,1); wwb = xtbi(j,k,3,1);
            }
            cqt = cq + cqu*std::sin(twopi*rfreq*(double)time_);
            if ((float)rfreq < 0.0f) {
                tcqt = cqu*(double)time_;
                if (std::abs(tcqt) < std::abs(cq)) { cqt = tcqt; } else { cqt = cq; }
            }
            sjet = std::sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz);
            if (sjet != 0.0) {
                sjetx = sjetx/sjet; sjety = sjety/sjet; sjetz = sjetz/sjet;
            } else {
                sjetx = si(j,k,1,1); sjety = si(j,k,1,2); sjetz = si(j,k,1,3);
            }
            uub = uub + sjetx*(double)xmach*cqt*c2/((double)gamma_*pb);
            vvb = vvb + sjety*(double)xmach*cqt*c2/((double)gamma_*pb);
            wwb = wwb + sjetz*(double)xmach*cqt*c2/((double)gamma_*pb);
            qi0(j,k,1,1) = (double)gamma_*pb/c2;
            qi0(j,k,2,1) = uub; qi0(j,k,3,1) = vvb; qi0(j,k,4,1) = wwb;
            qi0(j,k,5,1) = pb;
            bci(j,k,1)   = 1.0;
            f23 = 0.0;
            i2  = std::min(2, idim1);
            if (i2 == 1) f23 = 0.0;
            z1 = 2.0+1.5*f23; z2 = -0.5*f23; z3 = -(2.0+f23);
            qi0(j,k,1,2) = z1*q(j,k,1,1)+z2*q(j,k,i2,1)+z3*qi0(j,k,1,1);
            qi0(j,k,2,2) = z1*q(j,k,1,2)+z2*q(j,k,i2,2)+z3*qi0(j,k,2,1);
            qi0(j,k,3,2) = z1*q(j,k,1,3)+z2*q(j,k,i2,3)+z3*qi0(j,k,3,1);
            qi0(j,k,4,2) = z1*q(j,k,1,4)+z2*q(j,k,i2,4)+z3*qi0(j,k,4,1);
            qi0(j,k,5,2) = z1*q(j,k,1,5)+z2*q(j,k,i2,5)+z3*qi0(j,k,5,1);
        }
    }

    if (ivmx >= 2) {
    if (level >= lglobal && ntime != 0) {
        if (cmn_wallfun.iwf[0] == 0) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                vi0(j,k,1,1) = 0.0; vi0(j,k,1,2) = 0.0;
            }
        } else {
            // ***Wall Function begin
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                i = 1;
                c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
                uu = std::sqrt((q(j,k,i,2)-qi0(j,k,2,1))*(q(j,k,i,2)-qi0(j,k,2,1))
                             + (q(j,k,i,3)-qi0(j,k,3,1))*(q(j,k,i,3)-qi0(j,k,3,1))
                             + (q(j,k,i,4)-qi0(j,k,4,1))*(q(j,k,i,4)-qi0(j,k,4,1)));
                if (ivmx == 2) { dist = sni0(j,k,i); } else { dist = std::abs(smin(j,k,i)); }
                tt   = (double)gamma_*qi0(j,k,5,1)/qi0(j,k,1,1);
                fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                if ((float)rc <= 20.24f) {
                    xnplus = std::sqrt(rc);
                } else if ((float)rc <= 435.0f) {
                    xnplus = a0[1]+a0[2]*rc+(a0[3]*rc)*rc+(a0[4]*rc)*rc*rc
                           +(a0[5]*rc)*rc*rc*rc+(a0[6]*rc)*rc*rc*rc*rc+(a0[7]*rc)*rc*rc*rc*rc*rc;
                } else if ((float)rc <= 4000.0f) {
                    xnplus = a1[1]+a1[2]*rc+(a1[3]*rc)*rc+(a1[4]*rc)*rc*rc+(a1[5]*rc)*rc*rc*rc;
                } else {
                    xnplus = a2[1]+a2[2]*rc+a2[3]*rc*rc;
                }
                xnplussav = xnplus;
                if ((float)xnplus > 10.0f) {
                    bool conv1590 = false;
                    for (num = 1; num <= 10; num++) {
                        f_val = rc/xnplus-2.44*std::log(xnplus)-5.2;
                        dfdn  = -rc/(xnplus*xnplus)-2.44/xnplus;
                        delta = -f_val/dfdn;
                        xnplus = std::abs(xnplus+delta);
                        if (std::abs((float)delta) < 1.e-3f) { conv1590 = true; break; }
                    }
                    if (!conv1590) xnplus = xnplussav;
                }
                dudy = uu/dist;
                xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/(q(j,k,i,1)*dudy)*(double)xmach/(double)reue-1.0);
                vi0(j,k,1,1) = xmut; vi0(j,k,1,2) = 0.0;
            }
            // ***Wall Function end
        }
    }
    } // ivmx >= 2


    // only need to do advanced model turbulence B.C.s on finest grid (nface==1)
    if (level >= lglobal && ntime != 0) {
    if (ivmx == 4 || ivmx == 5) {
        for (k = ksta; k <= kend1; k++)
        for (j = jsta; j <= jend1; j++) {
            ti0(j,k,1,1) = -tursav(j,k,1,1); ti0(j,k,2,1) = -tursav(j,k,1,2);
            ti0(j,k,1,2) = 2.0*ti0(j,k,1,1)-tursav(j,k,1,1);
            ti0(j,k,2,2) = 2.0*ti0(j,k,2,1)-tursav(j,k,1,2);
        }
    }
    if (ivmx >= 6) {
        c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
        re  = (double)reue/(double)xmach;
        if (ivmx == 6) { beta1 = 0.075; }
        else if (ivmx == 7 || ivmx == 30 || ivmx == 40) { beta1 = 0.075; }
        else if (ivmx == 8 || ivmx == 12 || ivmx == 14) { beta1 = 0.83; }
        else { beta1 = 0.075; }
        i = 1;
        if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*qi0(j,k,5,1)/qi0(j,k,1,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                ti0(j,k,1,1) = 2.0*(2.0*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)-tursav(j,k,i,1);
                ti0(j,k,2,1) = -tursav(j,k,i,2);
                ti0(j,k,1,2) = 2.0*ti0(j,k,1,1)-tursav(j,k,1,1);
                ti0(j,k,2,2) = 2.0*ti0(j,k,2,1)-tursav(j,k,1,2);
            }
        } else if (ivmx == 15) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = 0.5*(tursav(j,k,i,2)+tursav(j,k,i+1,2))/3.0;
                ti0(j,k,1,1) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i))-tursav(j,k,i,1);
                ti0(j,k,2,1) = -tursav(j,k,i,2);
                ti0(j,k,1,2) = 0.0; ti0(j,k,2,2) = 0.0;
            }
        } else if (ivmx == 16) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                ti0(j,k,1,1) = -tursav(j,k,i,1); ti0(j,k,2,1) = -tursav(j,k,i,2);
                ti0(j,k,1,2) = 0.0; ti0(j,k,2,2) = 0.0;
            }
        } else {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*qi0(j,k,5,1)/qi0(j,k,1,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dist = std::abs(smin(j,k,i));
                ti0(j,k,1,1) = 2.0*(60.0*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))-tursav(j,k,i,1);
                ti0(j,k,2,1) = -tursav(j,k,i,2);
                ti0(j,k,1,2) = 2.0*ti0(j,k,1,1)-tursav(j,k,1,1);
                ti0(j,k,2,2) = 2.0*ti0(j,k,2,1)-tursav(j,k,1,2);
            }
        }
        if (ivmx == 30) {
            for (k = ksta; k <= kend1; k++) { kk = k-ksta+1;
            for (j = jsta; j <= jend1; j++) { jj = j-jsta+1;
                if (iuse3 == 0) {
                    ti0(j,k,3,1) = tursav(j,k,i,3);
                } else {
                    tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                    fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                    dist = std::abs(smin(j,k,i));
                    hee  = (1.0+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                    rhoee = std::pow((double)gamma_*q(j,k,i,5), 1.0/(double)gamma_);
                    uee2 = ((double)gm1*rhoee*hee-(double)gamma_*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                    if ((float)uee2 > 0.0f) { uee = std::sqrt(uee2); } else { uee = (double)xmach; }
                    ti0(j,k,3,1) = tursav(j,k,i,3)+2.0*dist*q(j,k,i,1)/fnu
                                  *uee*bcdata(jj,kk,ip,3)*(double)reue/(double)xmach;
                }
                ti0(j,k,3,2) = 2.0*ti0(j,k,3,1)-tursav(j,k,1,3);
            }}
        } else if (ivmx == 40) {
            for (k = ksta; k <= kend1; k++) { kk = k-ksta+1;
            for (j = jsta; j <= jend1; j++) { jj = j-jsta+1;
                ti0(j,k,3,1) = tursav(j,k,i,3); ti0(j,k,4,1) = tursav(j,k,i,4);
            }}
        }
    } // ivmx >= 6
    } // level >= lglobal && ntime != 0

    } // nface == 1



    // **************************************************************************
    // i=idim boundary   viscous wall with T & cq specified   bctype 2016
    // **************************************************************************
    if (nface == 2) {

    for (k = ksta; k <= kend1; k++) {
        kk = k - ksta + 1;
        for (j = jsta; j <= jend1; j++) {
            jj = j - jsta + 1;
            cq    = bcdata(jj,kk,ip,2); cqu   = bcdata(jj,kk,ip,3);
            sjetx = bcdata(jj,kk,ip,4); sjety = bcdata(jj,kk,ip,5);
            sjetz = bcdata(jj,kk,ip,6); rfreq = bcdata(jj,kk,ip,7);
            pb    = q(j,k,idim1,5);
            dpb   = q(j,k,idim1,5) - q(j,k,idim1-1,5);
            pb    = pb + dpb/2.0;
            if ((float)pb <= 0.0f) pb = q(j,k,idim1,5);
            c2    = (double)gamma_*q(j,k,idim1,5)/q(j,k,idim1,1);
            if ((float)bcdata(jj,kk,ip,1) > 0.0f) {
                c2 = bcdata(jj,kk,ip,1);
            } else if ((float)bcdata(jj,kk,ip,1) < 0.0f) {
                c2 = 1.0+(double)gm1*0.5*(double)xmach*(double)xmach;
            } else {
                if (iuns > 0 && irelv > 0) {
                    xm2 = (q(j,k,idim1,2)-xtbi(j,k,1,2))*(q(j,k,idim1,2)-xtbi(j,k,1,2))
                        + (q(j,k,idim1,3)-xtbi(j,k,2,2))*(q(j,k,idim1,3)-xtbi(j,k,2,2))
                        + (q(j,k,idim1,4)-xtbi(j,k,3,2))*(q(j,k,idim1,4)-xtbi(j,k,3,2));
                } else {
                    xm2 = q(j,k,idim1,2)*q(j,k,idim1,2)+q(j,k,idim1,3)*q(j,k,idim1,3)
                         +q(j,k,idim1,4)*q(j,k,idim1,4);
                }
                xm2 = xm2/c2; c2 = c2*(1.0+0.5*(double)gm1*xm2);
            }
            uub = 0.0; vvb = 0.0; wwb = 0.0;
            if (iuns > 0 && irelv > 0) {
                uub = xtbi(j,k,1,2); vvb = xtbi(j,k,2,2); wwb = xtbi(j,k,3,2);
            }
            cqt = cq + cqu*std::sin(twopi*rfreq*(double)time_);
            if ((float)rfreq < 0.0f) {
                tcqt = cqu*(double)time_;
                if (std::abs(tcqt) < std::abs(cq)) { cqt = tcqt; } else { cqt = cq; }
            }
            sjet = std::sqrt(sjetx*sjetx+sjety*sjety+sjetz*sjetz);
            if (sjet != 0.0) {
                sjetx = sjetx/sjet; sjety = sjety/sjet; sjetz = sjetz/sjet;
            } else {
                sjetx = -si(j,k,idim,1); sjety = -si(j,k,idim,2); sjetz = -si(j,k,idim,3);
            }
            uub = uub + sjetx*(double)xmach*cqt*c2/((double)gamma_*pb);
            vvb = vvb + sjety*(double)xmach*cqt*c2/((double)gamma_*pb);
            wwb = wwb + sjetz*(double)xmach*cqt*c2/((double)gamma_*pb);
            qi0(j,k,1,3) = (double)gamma_*pb/c2;
            qi0(j,k,2,3) = uub; qi0(j,k,3,3) = vvb; qi0(j,k,4,3) = wwb;
            qi0(j,k,5,3) = pb;
            bci(j,k,2)   = 1.0;
            f23 = 0.0;
            i2  = std::max(1, idim-2);
            if (i2 == 1) f23 = 0.0;
            z1 = -2.0-1.5*f23; z2 = 0.5*f23; z3 = 2.0+f23;
            qi0(j,k,1,4) = z1*q(j,k,idim1,1)+z2*q(j,k,i2,1)+z3*qi0(j,k,1,3);
            qi0(j,k,2,4) = z1*q(j,k,idim1,2)+z2*q(j,k,i2,2)+z3*qi0(j,k,2,3);
            qi0(j,k,3,4) = z1*q(j,k,idim1,3)+z2*q(j,k,i2,3)+z3*qi0(j,k,3,3);
            qi0(j,k,4,4) = z1*q(j,k,idim1,4)+z2*q(j,k,i2,4)+z3*qi0(j,k,4,3);
            qi0(j,k,5,4) = z1*q(j,k,idim1,5)+z2*q(j,k,i2,5)+z3*qi0(j,k,5,3);
        }
    }

    if (ivmx >= 2) {
    if (level >= lglobal && ntime != 0) {
        if (cmn_wallfun.iwf[0] == 0) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                vi0(j,k,1,3) = 0.0; vi0(j,k,1,4) = 0.0;
            }
        } else {
            // ***Wall Function begin
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                i = idim-1;
                c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
                uu = std::sqrt((q(j,k,i,2)-qi0(j,k,2,3))*(q(j,k,i,2)-qi0(j,k,2,3))
                             + (q(j,k,i,3)-qi0(j,k,3,3))*(q(j,k,i,3)-qi0(j,k,3,3))
                             + (q(j,k,i,4)-qi0(j,k,4,3))*(q(j,k,i,4)-qi0(j,k,4,3)));
                if (ivmx == 2) { dist = snim(j,k,i); } else { dist = std::abs(smin(j,k,i)); }
                tt   = (double)gamma_*qi0(j,k,5,3)/qi0(j,k,1,3);
                fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                rc   = q(j,k,i,1)*uu*dist/fnuw*(double)reue/(double)xmach;
                if ((float)rc <= 20.24f) {
                    xnplus = std::sqrt(rc);
                } else if ((float)rc <= 435.0f) {
                    xnplus = a0[1]+a0[2]*rc+(a0[3]*rc)*rc+(a0[4]*rc)*rc*rc
                           +(a0[5]*rc)*rc*rc*rc+(a0[6]*rc)*rc*rc*rc*rc+(a0[7]*rc)*rc*rc*rc*rc*rc;
                } else if ((float)rc <= 4000.0f) {
                    xnplus = a1[1]+a1[2]*rc+(a1[3]*rc)*rc+(a1[4]*rc)*rc*rc+(a1[5]*rc)*rc*rc*rc;
                } else {
                    xnplus = a2[1]+a2[2]*rc+a2[3]*rc*rc;
                }
                xnplussav = xnplus;
                if ((float)xnplus > 10.0f) {
                    bool conv1690 = false;
                    for (num = 1; num <= 10; num++) {
                        f_val = rc/xnplus-2.44*std::log(xnplus)-5.2;
                        dfdn  = -rc/(xnplus*xnplus)-2.44/xnplus;
                        delta = -f_val/dfdn;
                        xnplus = std::abs(xnplus+delta);
                        if (std::abs((float)delta) < 1.e-3f) { conv1690 = true; break; }
                    }
                    if (!conv1690) xnplus = xnplussav;
                }
                dudy = uu/dist;
                xmut = fnuw*(xnplus*xnplus*fnuw/(dist*dist)/(q(j,k,i,1)*dudy)*(double)xmach/(double)reue-1.0);
                vi0(j,k,1,3) = xmut; vi0(j,k,1,4) = 0.0;
            }
            // ***Wall Function end
        }
    }
    } // ivmx >= 2


    // only need to do advanced model turbulence B.C.s on finest grid (nface==2)
    if (level >= lglobal && ntime != 0) {
    if (ivmx == 4 || ivmx == 5) {
        for (k = ksta; k <= kend1; k++)
        for (j = jsta; j <= jend1; j++) {
            ti0(j,k,1,3) = -tursav(j,k,idim-1,1); ti0(j,k,2,3) = -tursav(j,k,idim-1,2);
            ti0(j,k,1,4) = 2.0*ti0(j,k,1,3)-tursav(j,k,idim-1,1);
            ti0(j,k,2,4) = 2.0*ti0(j,k,2,3)-tursav(j,k,idim-1,2);
        }
    }
    if (ivmx >= 6) {
        c2b = (double)cbar/(double)tinf; c2bp = c2b+1.0;
        re  = (double)reue/(double)xmach;
        if (ivmx == 6) { beta1 = 0.075; }
        else if (ivmx == 7 || ivmx == 30 || ivmx == 40) { beta1 = 0.075; }
        else if (ivmx == 8 || ivmx == 12 || ivmx == 14) { beta1 = 0.83; }
        else { beta1 = 0.075; }
        i = idim-1;
        if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*qi0(j,k,5,3)/qi0(j,k,1,3);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = std::sqrt(tursav(j,k,i,2))/std::abs(smin(j,k,i));
                ti0(j,k,1,3) = 2.0*(2.0*fnu/(q(j,k,i,1)*re*re)*dkdy*dkdy)-tursav(j,k,i,1);
                ti0(j,k,2,3) = -tursav(j,k,i,2);
                ti0(j,k,1,4) = 2.0*ti0(j,k,1,3)-tursav(j,k,idim-1,1);
                ti0(j,k,2,4) = 2.0*ti0(j,k,2,3)-tursav(j,k,idim-1,2);
            }
        } else if (ivmx == 15) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dkdy = 0.5*(tursav(j,k,i,2)+tursav(j,k,i-1,2))/3.0;
                ti0(j,k,1,3) = 2.0*dkdy/(re*smin(j,k,i))/(re*smin(j,k,i))-tursav(j,k,i,1);
                ti0(j,k,2,3) = -tursav(j,k,i,2);
                ti0(j,k,1,4) = 0.0; ti0(j,k,2,4) = 0.0;
            }
        } else if (ivmx == 16) {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                ti0(j,k,1,3) = -tursav(j,k,i,1); ti0(j,k,2,3) = -tursav(j,k,i,2);
                ti0(j,k,1,4) = 0.0; ti0(j,k,2,4) = 0.0;
            }
        } else {
            for (k = ksta; k <= kend1; k++)
            for (j = jsta; j <= jend1; j++) {
                tt   = (double)gamma_*qi0(j,k,5,3)/qi0(j,k,1,3);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                dist = std::abs(smin(j,k,i));
                ti0(j,k,1,3) = 2.0*(60.0*fnu/(re*re*q(j,k,i,1)*beta1*dist*dist))-tursav(j,k,i,1);
                ti0(j,k,2,3) = -tursav(j,k,i,2);
                ti0(j,k,1,4) = 2.0*ti0(j,k,1,3)-tursav(j,k,idim-1,1);
                ti0(j,k,2,4) = 2.0*ti0(j,k,2,3)-tursav(j,k,idim-1,2);
            }
        }
        if (ivmx == 30) {
            for (k = ksta; k <= kend1; k++) { kk = k-ksta+1;
            for (j = jsta; j <= jend1; j++) { jj = j-jsta+1;
                if (iuse3 == 0) {
                    ti0(j,k,3,3) = tursav(j,k,i,3);
                } else {
                    tt   = (double)gamma_*q(j,k,i,5)/q(j,k,i,1);
                    fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                    dist = std::abs(smin(j,k,i));
                    hee  = (1.0+0.5*(double)gm1*(double)xmach*(double)xmach)/(double)gm1;
                    rhoee = std::pow((double)gamma_*q(j,k,i,5), 1.0/(double)gamma_);
                    uee2 = ((double)gm1*rhoee*hee-(double)gamma_*q(j,k,i,5))/(0.5*(double)gm1*rhoee);
                    if ((float)uee2 > 0.0f) { uee = std::sqrt(uee2); } else { uee = (double)xmach; }
                    ti0(j,k,3,3) = tursav(j,k,i,3)+2.0*dist*q(j,k,i,1)/fnu
                                  *uee*bcdata(jj,kk,ip,3)*(double)reue/(double)xmach;
                }
                ti0(j,k,3,4) = 2.0*ti0(j,k,3,3)-tursav(j,k,idim-1,3);
            }}
        } else if (ivmx == 40) {
            for (k = ksta; k <= kend1; k++) { kk = k-ksta+1;
            for (j = jsta; j <= jend1; j++) { jj = j-jsta+1;
                ti0(j,k,3,3) = tursav(j,k,i,3); ti0(j,k,4,3) = tursav(j,k,i,4);
            }}
        }
    } // ivmx >= 6
    } // level >= lglobal && ntime != 0

    } // nface == 2

    return;
} // bc2016

} // namespace bc2016_ns
