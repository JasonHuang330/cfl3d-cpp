// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "csurf.h"
#include <cmath>
#include <cstdint>

namespace csurf_ns {

void csurf(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si, FortranArray4DRef<double> q, int& ics, FortranArray3DRef<double> ub, FortranArray3DRef<double> vb, FortranArray3DRef<double> wb, FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui, FortranArray3DRef<double> vol, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> blank, FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi, int& iuns, int& ncs, FortranArray2DRef<int> icsinfo, double& sx, double& sy, double& sz, double& stot, double& pav, double& ptav, double& tav, double& ttav, double& xmav, double& fmdot, double& cfxp, double& cfyp, double& cfzp, double& cfdp, double& cflp, double& cftp, double& cfxv, double& cfyv, double& cfzv, double& cfdv, double& cflv, double& cftv, double& cfxmom, double& cfymom, double& cfzmom, double& cfdmom, double& cflmom, double& cftmom, double& cfxtot, double& cfytot, double& cfztot, double& cfdtot, double& cfltot, double& cfttot, int& maxcs, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0)
{
    // COMMON block aliases
    float& gamma  = cmn_fluid.gamma;
    float& gm1    = cmn_fluid.gm1;
    float& sref   = cmn_fsum.sref;
    float& xmach  = cmn_info.xmach;
    float& alpha  = cmn_info.alpha;
    float& beta   = cmn_info.beta;
    float& p0     = cmn_ivals.p0;
    float& reue   = cmn_reyue.reue;
    int32_t* ivisc = cmn_reyue.ivisc; // 1-based: ivisc[0]=ivisc(1), [1]=ivisc(2), [2]=ivisc(3)

    // Local variables
    int ist, ifn, jst, jfn, kst, kfn, iwall;
    double fnorm;
    double cpc, cst, cosa, sina, cosb, sinb;
    int k, kc, kcm1, m, mm, kk;
    int j, jc, jcm1, jj;
    int i, ic, icm1, ii;
    double aa, bb, dcp, q1, q2, q3, q4, p1, t1, xm1, tt, term1, pt, cmdot;
    double urel, vrel, wrel, tau, vnorm;


    // Extract icsinfo values
    ist   = icsinfo(ics,2);
    ifn   = icsinfo(ics,3) - 1;
    jst   = icsinfo(ics,4);
    jfn   = icsinfo(ics,5) - 1;
    kst   = icsinfo(ics,6);
    kfn   = icsinfo(ics,7) - 1;
    iwall = icsinfo(ics,8);
    fnorm = (double)icsinfo(ics,9);
    if (icsinfo(ics,9) == 0) fnorm = 1.0;

    cpc  = 2.e0 / ((double)gamma * (double)xmach * (double)xmach);
    cst  = 4.0 / ((double)reue * (double)xmach);

    cosa = std::cos((double)alpha);
    sina = std::sin((double)alpha);
    cosb = std::cos((double)beta);
    sinb = std::sin((double)beta);

    sx     = 0.e0;
    sy     = 0.e0;
    sz     = 0.e0;
    pav    = 0.e0;
    ptav   = 0.e0;
    tav    = 0.e0;
    ttav   = 0.e0;
    xmav   = 0.e0;
    fmdot  = 0.e0;
    cfxp   = 0.e0;
    cfyp   = 0.e0;
    cfzp   = 0.e0;
    cfxv   = 0.e0;
    cfyv   = 0.e0;
    cfzv   = 0.e0;
    cfxmom = 0.e0;
    cfymom = 0.e0;
    cfzmom = 0.e0;

    //******************************************************************
    // forces on k=constant surfaces
    //******************************************************************

    if (kst == kfn+1) {
        k    = kst;
        kc   = k;
        kcm1 = k - 1;
        m    = 0;
        kk   = 0; // initialize to avoid uninitialized use

        if (k == 1) {
            k    = k;
            kc   = 1;
            kcm1 = 1;
            m    = 1;
            mm   = 1;
        }

        if (k == kdim) {
            k    = kdim;
            kc   = kdim - 1;
            kcm1 = kdim - 1;
            m    = 2;
            mm   = 3;
        }

        if (k == 1 && iwall == 1) {
            kk = 1;
        }

        if (k == kdim && iwall == 1) {
            kk = 2;
        }

        for (int ii_i = ist; ii_i <= ifn; ii_i++) {
        for (int jj_j = jst; jj_j <= jfn; jj_j++) {
            i = ii_i;
            j = jj_j;

            // area integration
            sx = sx + sk(j,k,i,1)*sk(j,k,i,4);
            sy = sy + sk(j,k,i,2)*sk(j,k,i,4);
            sz = sz + sk(j,k,i,3)*sk(j,k,i,4);

            // pressure force integration
            if (m != 0) {
                aa  = 1.0 + bck(j,i,m);
                bb  = 1.0 - bck(j,i,m);
                dcp = -(0.5*(aa*qk0(j,i,5,mm) + bb*q(j,kc,i,5))/(double)p0
                       -1.e0)*cpc*sk(j,k,i,4);
            } else {
                dcp = -(0.5*(q(j,kc,i,5)+q(j,kcm1,i,5))/(double)p0
                       -1.e0)*cpc*sk(j,k,i,4);
            }
            cfxp = cfxp - dcp*sk(j,k,i,1);
            cfyp = cfyp - dcp*sk(j,k,i,2);
            cfzp = cfzp - dcp*sk(j,k,i,3);

            // mass and momentum integration
            if (iwall == 0) {
                if (m != 0) {
                    aa = 1.0 + bck(j,i,m);
                    bb = 1.0 - bck(j,i,m);
                    q1 = 0.5*(aa*qk0(j,i,1,mm) + bb*q(j,kc,i,1));
                    q2 = 0.5*(aa*qk0(j,i,2,mm) + bb*q(j,kc,i,2)) / (double)xmach;
                    q3 = 0.5*(aa*qk0(j,i,3,mm) + bb*q(j,kc,i,3)) / (double)xmach;
                    q4 = 0.5*(aa*qk0(j,i,4,mm) + bb*q(j,kc,i,4)) / (double)xmach;
                    p1 = 0.5*(aa*qk0(j,i,5,mm) + bb*q(j,kc,i,5)) * (double)gamma;
                } else {
                    q1 = 0.5*(q(j,kc,i,1) + q(j,kcm1,i,1));
                    q2 = 0.5*(q(j,kc,i,2) + q(j,kcm1,i,2)) / (double)xmach;
                    q3 = 0.5*(q(j,kc,i,3) + q(j,kcm1,i,3)) / (double)xmach;
                    q4 = 0.5*(q(j,kc,i,4) + q(j,kcm1,i,4)) / (double)xmach;
                    p1 = 0.5*(q(j,kc,i,5) + q(j,kcm1,i,5)) * (double)gamma;
                }
                t1    = p1/q1;
                xm1   = std::sqrt((double)xmach*(double)xmach*(q2*q2+q3*q3+q4*q4)/t1);
                tt    = t1*(1.0+0.5*(double)gm1*xm1*xm1);
                term1 = std::pow(1.0+0.5*(double)gm1*xm1*xm1, (double)gamma/(double)gm1);
                pt    = p1*term1;
                cmdot = (q2*sk(j,k,i,1) + q3*sk(j,k,i,2) +
                         q4*sk(j,k,i,3))*q1*sk(j,k,i,4);
                fmdot  = fmdot  + cmdot;
                cfxmom = cfxmom + 2.0*q2*cmdot;
                cfymom = cfymom + 2.0*q3*cmdot;
                cfzmom = cfzmom + 2.0*q4*cmdot;

                pav  = pav  + p1*cmdot;
                ptav = ptav + pt*cmdot;
                tav  = tav  + t1*cmdot;
                ttav = ttav + tt*cmdot;
                xmav = xmav + xm1*cmdot;
            }

            // skin friction integration
            if (iwall == 1 && ivisc[2] > 0) {
                urel = ub(j,kc,i);
                vrel = vb(j,kc,i);
                wrel = wb(j,kc,i);
                if (iuns > 0) {
                    urel = ub(j,kc,i) - xtbk(j,i,1,kk);
                    vrel = vb(j,kc,i) - xtbk(j,i,2,kk);
                    wrel = wb(j,kc,i) - xtbk(j,i,3,kk);
                }
                tau   = vmuk(j,i,kk)*cst/vol(j,kc,i)*sk(j,k,i,4)*sk(j,k,i,4);
                vnorm = (urel*sk(j,k,i,1)+vrel*sk(j,k,i,2)+wrel*sk(j,k,i,3));
                cfxv  = cfxv + tau*(urel-vnorm*sk(j,k,i,1));
                cfyv  = cfyv + tau*(vrel-vnorm*sk(j,k,i,2));
                cfzv  = cfzv + tau*(wrel-vnorm*sk(j,k,i,3));
            }

        } // j loop
        } // i loop

        if (std::abs((float)fmdot) < 0.00001f) {
            pav  = 0.e0;
            ptav = 0.e0;
            tav  = 0.e0;
            ttav = 0.e0;
            xmav = 0.e0;
        } else {
            pav  = pav  / fmdot;
            ptav = ptav / fmdot;
            tav  = tav  / fmdot;
            ttav = ttav / fmdot;
            xmav = xmav / fmdot;
        }

        sx     = sx*fnorm;
        sy     = sy*fnorm;
        sz     = sz*fnorm;
        fmdot  = fmdot*fnorm;
        cfxp   = cfxp*fnorm/(double)sref;
        cfyp   = cfyp*fnorm/(double)sref;
        cfzp   = cfzp*fnorm/(double)sref;
        cfxv   = cfxv/(double)sref;
        cfyv   = cfyv/(double)sref;
        cfzv   = cfzv/(double)sref;
        cfxmom = cfxmom*fnorm/(double)sref;
        cfymom = cfymom*fnorm/(double)sref;
        cfzmom = cfzmom*fnorm/(double)sref;

        stot   = std::sqrt(sx*sx+sy*sy+sz*sz);
        cfxtot = cfxp+cfxv+cfxmom;
        cfytot = cfyp+cfyv+cfymom;
        cfztot = cfzp+cfzv+cfzmom;
    } // end k=constant surface block


    //******************************************************************
    // forces on j=constant surfaces
    //******************************************************************

    if (jst == jfn+1) {
        j    = jst;
        jc   = j;
        jcm1 = j - 1;
        m    = 0;
        jj   = 0; // initialize to avoid uninitialized use

        if (j == 1) {
            j    = j;
            jc   = 1;
            jcm1 = 1;
            m    = 1;
            mm   = 1;
        }

        if (j == jdim) {
            j    = jdim;
            jc   = jdim - 1;
            jcm1 = jdim - 1;
            m    = 2;
            mm   = 3;
        }

        if (j == 1 && iwall == 1) {
            jj = 1;
        }

        if (j == jdim && iwall == 1) {
            jj = 2;
        }

        for (int ii_i = ist; ii_i <= ifn; ii_i++) {
        for (int kk_k = kst; kk_k <= kfn; kk_k++) {
            i = ii_i;
            k = kk_k;

            // area integration
            sx = sx + sj(j,k,i,1)*sj(j,k,i,4);
            sy = sy + sj(j,k,i,2)*sj(j,k,i,4);
            sz = sz + sj(j,k,i,3)*sj(j,k,i,4);

            // pressure force integration
            if (m != 0) {
                aa  = 1.0 + bcj(k,i,m);
                bb  = 1.0 - bcj(k,i,m);
                dcp = -(0.5*(aa*qj0(k,i,5,mm) + bb*q(jc,k,i,5))/(double)p0
                       -1.e0)*cpc*sj(j,k,i,4);
            } else {
                dcp = -(0.5*(q(jc,k,i,5)+q(jcm1,k,i,5))/(double)p0
                       -1.e0)*cpc*sj(j,k,i,4);
            }
            cfxp = cfxp - dcp*sj(j,k,i,1);
            cfyp = cfyp - dcp*sj(j,k,i,2);
            cfzp = cfzp - dcp*sj(j,k,i,3);

            // mass and momentum integration
            if (iwall == 0) {
                if (m != 0) {
                    aa = 1.0 + bcj(k,i,m);
                    bb = 1.0 - bcj(k,i,m);
                    q1 = 0.5*(aa*qj0(k,i,1,mm) + bb*q(jc,k,i,1));
                    q2 = 0.5*(aa*qj0(k,i,2,mm) + bb*q(jc,k,i,2)) / (double)xmach;
                    q3 = 0.5*(aa*qj0(k,i,3,mm) + bb*q(jc,k,i,3)) / (double)xmach;
                    q4 = 0.5*(aa*qj0(k,i,4,mm) + bb*q(jc,k,i,4)) / (double)xmach;
                    p1 = 0.5*(aa*qj0(k,i,5,mm) + bb*q(jc,k,i,5)) * (double)gamma;
                } else {
                    q1 = 0.5*(q(jc,k,i,1) + q(jcm1,k,i,1));
                    q2 = 0.5*(q(jc,k,i,2) + q(jcm1,k,i,2)) / (double)xmach;
                    q3 = 0.5*(q(jc,k,i,3) + q(jcm1,k,i,3)) / (double)xmach;
                    q4 = 0.5*(q(jc,k,i,4) + q(jcm1,k,i,4)) / (double)xmach;
                    p1 = 0.5*(q(jc,k,i,5) + q(jcm1,k,i,5)) * (double)gamma;
                }
                t1    = p1/q1;
                xm1   = std::sqrt((double)xmach*(double)xmach*(q2*q2+q3*q3+q4*q4)/t1);
                tt    = t1*(1.0+0.5*(double)gm1*xm1*xm1);
                term1 = std::pow(1.0+0.5*(double)gm1*xm1*xm1, (double)gamma/(double)gm1);
                pt    = p1*term1;
                cmdot = (q2*sj(j,k,i,1) + q3*sj(j,k,i,2) +
                         q4*sj(j,k,i,3))*q1*sj(j,k,i,4);
                fmdot  = fmdot  + cmdot;
                cfxmom = cfxmom + 2.0*q2*cmdot;
                cfymom = cfymom + 2.0*q3*cmdot;
                cfzmom = cfzmom + 2.0*q4*cmdot;

                pav  = pav  + p1*cmdot;
                ptav = ptav + pt*cmdot;
                tav  = tav  + t1*cmdot;
                ttav = ttav + tt*cmdot;
                xmav = xmav + xm1*cmdot;
            }

            // skin friction integration
            if (iwall == 1 && ivisc[1] > 0) {
                urel = ub(jc,k,i);
                vrel = vb(jc,k,i);
                wrel = wb(jc,k,i);
                if (iuns > 0) {
                    urel = ub(jc,k,i) - xtbj(k,i,1,jj);
                    vrel = vb(jc,k,i) - xtbj(k,i,2,jj);
                    wrel = wb(jc,k,i) - xtbj(k,i,3,jj);
                }
                tau   = vmuj(k,i,jj)*cst/vol(jc,k,i)*sj(j,k,i,4)*sj(j,k,i,4);
                vnorm = (urel*sj(j,k,i,1)+vrel*sj(j,k,i,2)+wrel*sj(j,k,i,3));
                cfxv  = cfxv + tau*(urel-vnorm*sj(j,k,i,1));
                cfyv  = cfyv + tau*(vrel-vnorm*sj(j,k,i,2));
                cfzv  = cfzv + tau*(wrel-vnorm*sj(j,k,i,3));
            }

        } // k loop
        } // i loop

        if (std::abs((float)fmdot) < 0.00001f) {
            pav  = 0.e0;
            ptav = 0.e0;
            tav  = 0.e0;
            ttav = 0.e0;
            xmav = 0.e0;
        } else {
            pav  = pav  / fmdot;
            ptav = ptav / fmdot;
            tav  = tav  / fmdot;
            ttav = ttav / fmdot;
            xmav = xmav / fmdot;
        }

        sx     = sx*fnorm;
        sy     = sy*fnorm;
        sz     = sz*fnorm;
        fmdot  = fmdot*fnorm;
        cfxp   = cfxp*fnorm/(double)sref;
        cfyp   = cfyp*fnorm/(double)sref;
        cfzp   = cfzp*fnorm/(double)sref;
        cfxv   = cfxv/(double)sref;
        cfyv   = cfyv/(double)sref;
        cfzv   = cfzv/(double)sref;
        cfxmom = cfxmom*fnorm/(double)sref;
        cfymom = cfymom*fnorm/(double)sref;
        cfzmom = cfzmom*fnorm/(double)sref;

        stot   = std::sqrt(sx*sx+sy*sy+sz*sz);
        cfxtot = cfxp+cfxv+cfxmom;
        cfytot = cfyp+cfyv+cfymom;
        cfztot = cfzp+cfzv+cfzmom;

    } // end j=constant surface block


    //******************************************************************
    // forces on i=constant surfaces
    //******************************************************************

    if (ist == ifn+1) {
        i    = ist;
        ic   = i;
        icm1 = i - 1;
        ii   = 1;
        m    = 0;

        if (i == 1) {
            ic   = 1;
            icm1 = 1;
            if (iwall == 1) ii = 1;
            m    = 1;
            mm   = 1;
        }

        if (i == idim) {
            ic   = idim - 1;
            icm1 = idim - 1;
            if (iwall == 1) ii = 2;
            m    = 2;
            mm   = 3;
        }

        for (int jj_j = jst; jj_j <= jfn; jj_j++) {
        for (int kk_k = kst; kk_k <= kfn; kk_k++) {
            j = jj_j;
            k = kk_k;

            // area integration
            sx = sx + si(j,k,i,1)*si(j,k,i,4);
            sy = sy + si(j,k,i,2)*si(j,k,i,4);
            sz = sz + si(j,k,i,3)*si(j,k,i,4);

            // pressure force integration
            if (m != 0) {
                aa  = 1.0 + bci(j,k,m);
                bb  = 1.0 - bci(j,k,m);
                dcp = -(0.5*(aa*qi0(j,k,5,mm) + bb*q(j,k,ic,5))/(double)p0
                       -1.e0)*cpc*si(j,k,i,4);
            } else {
                dcp = -(0.5*(q(j,k,ic,5)+q(j,k,icm1,5))/(double)p0
                       -1.e0)*cpc*si(j,k,i,4);
            }
            cfxp = cfxp - dcp*si(j,k,i,1);
            cfyp = cfyp - dcp*si(j,k,i,2);
            cfzp = cfzp - dcp*si(j,k,i,3);

            // mass and momentum integration
            if (iwall == 0) {
                if (m != 0) {
                    aa = 1.0 + bci(j,k,m);
                    bb = 1.0 - bci(j,k,m);
                    q1 = 0.5*(aa*qi0(j,k,1,mm) + bb*q(j,k,ic,1));
                    q2 = 0.5*(aa*qi0(j,k,2,mm) + bb*q(j,k,ic,2)) / (double)xmach;
                    q3 = 0.5*(aa*qi0(j,k,3,mm) + bb*q(j,k,ic,3)) / (double)xmach;
                    q4 = 0.5*(aa*qi0(j,k,4,mm) + bb*q(j,k,ic,4)) / (double)xmach;
                    p1 = 0.5*(aa*qi0(j,k,5,mm) + bb*q(j,k,ic,5)) * (double)gamma;
                } else {
                    q1 = 0.5*(q(j,k,ic,1) + q(j,k,icm1,1));
                    q2 = 0.5*(q(j,k,ic,2) + q(j,k,icm1,2)) / (double)xmach;
                    q3 = 0.5*(q(j,k,ic,3) + q(j,k,icm1,3)) / (double)xmach;
                    q4 = 0.5*(q(j,k,ic,4) + q(j,k,icm1,4)) / (double)xmach;
                    p1 = 0.5*(q(j,k,ic,5) + q(j,k,icm1,5)) * (double)gamma;
                }
                t1    = p1/q1;
                xm1   = std::sqrt((double)xmach*(double)xmach*(q2*q2+q3*q3+q4*q4)/t1);
                tt    = t1*(1.0+0.5*(double)gm1*xm1*xm1);
                term1 = std::pow(1.0+0.5*(double)gm1*xm1*xm1, (double)gamma/(double)gm1);
                pt    = p1*term1;
                cmdot = (q2*si(j,k,i,1) + q3*si(j,k,i,2) +
                         q4*si(j,k,i,3))*q1*si(j,k,i,4);
                fmdot  = fmdot  + cmdot;
                cfxmom = cfxmom + 2.0*q2*cmdot;
                cfymom = cfymom + 2.0*q3*cmdot;
                cfzmom = cfzmom + 2.0*q4*cmdot;

                pav  = pav  + p1*cmdot;
                ptav = ptav + pt*cmdot;
                tav  = tav  + t1*cmdot;
                ttav = ttav + tt*cmdot;
                xmav = xmav + xm1*cmdot;
            }

            // skin friction integration
            if (iwall == 1 && ivisc[0] > 0) {
                urel = ub(j,k,ic);
                vrel = vb(j,k,ic);
                wrel = wb(j,k,ic);
                if (iuns > 0) {
                    urel = ub(j,k,ic) - xtbi(j,k,1,ii);
                    vrel = vb(j,k,ic) - xtbi(j,k,2,ii);
                    wrel = wb(j,k,ic) - xtbi(j,k,3,ii);
                }
                tau   = vmui(j,k,ii)*cst/vol(j,k,ic)*si(j,k,i,4)*si(j,k,i,4);
                vnorm = (urel*si(j,k,i,1)+vrel*si(j,k,i,2)+wrel*si(j,k,i,3));
                cfxv  = cfxv + tau*(urel-vnorm*si(j,k,i,1));
                cfyv  = cfyv + tau*(vrel-vnorm*si(j,k,i,2));
                cfzv  = cfzv + tau*(wrel-vnorm*si(j,k,i,3));
            }

        } // k loop
        } // j loop

        if (std::abs((float)fmdot) < 0.00001f) {
            pav  = 0.e0;
            ptav = 0.e0;
            tav  = 0.e0;
            ttav = 0.e0;
            xmav = 0.e0;
        } else {
            pav  = pav  / fmdot;
            ptav = ptav / fmdot;
            tav  = tav  / fmdot;
            ttav = ttav / fmdot;
            xmav = xmav / fmdot;
        }

        sx     = sx*fnorm;
        sy     = sy*fnorm;
        sz     = sz*fnorm;
        fmdot  = fmdot*fnorm;
        cfxp   = cfxp*fnorm/(double)sref;
        cfyp   = cfyp*fnorm/(double)sref;
        cfzp   = cfzp*fnorm/(double)sref;
        cfxv   = cfxv/(double)sref;
        cfyv   = cfyv/(double)sref;
        cfzv   = cfzv/(double)sref;
        cfxmom = cfxmom*fnorm/(double)sref;
        cfymom = cfymom*fnorm/(double)sref;
        cfzmom = cfzmom*fnorm/(double)sref;

        stot   = std::sqrt(sx*sx+sy*sy+sz*sz);
        cfxtot = cfxp+cfxv+cfxmom;
        cfytot = cfyp+cfyv+cfymom;
        cfztot = cfzp+cfzv+cfzmom;

    } // end i=constant surface block

    // calculate total magnitude and lift and drag components
    cftp   = std::sqrt(cfxp*cfxp+cfyp*cfyp+cfzp*cfzp);
    cftv   = std::sqrt(cfxv*cfxv+cfyv*cfyv+cfzv*cfzv);
    cftmom = std::sqrt(cfxmom*cfxmom+cfymom*cfymom+cfzmom*cfzmom);
    cfttot = std::sqrt(cfxtot*cfxtot+cfytot*cfytot+cfztot*cfztot);

    cfdp   =  cfxp*cosa*cosb+cfyp*cosa*sinb+cfzp*sina;
    cflp   = -cfxp*sina*cosb-cfyp*sina*sinb+cfzp*cosa;
    cfdv   =  cfxv*cosa*cosb+cfyv*cosa*sinb+cfzv*sina;
    cflv   = -cfxv*sina*cosb-cfyv*sina*sinb+cfzv*cosa;
    cfdmom =  cfxmom*cosa*cosb+cfymom*cosa*sinb+cfzmom*sina;
    cflmom = -cfxmom*sina*cosb-cfymom*sina*sinb+cfzmom*cosa;
    cfdtot =  cfxtot*cosa*cosb+cfytot*cosa*sinb+cfztot*sina;
    cfltot = -cfxtot*sina*cosb-cfytot*sina*sinb+cfztot*cosa;

} // end csurf

} // namespace csurf_ns
