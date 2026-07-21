// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2026.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace bc2026_ns {

void bc2026(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0, FortranArray4DRef<double> sj,
            FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci, FortranArray4DRef<double> xtbj,
            FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi,
            FortranArray4DRef<double> atbj, FortranArray4DRef<double> atbk,
            FortranArray4DRef<double> atbi, int& ista, int& iend,
            int& jsta, int& jend, int& ksta, int& kend, int& nface,
            FortranArray4DRef<double> tursav, FortranArray4DRef<double> tj0,
            FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> smin, FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0, int& mdim, int& ndim,
            FortranArray4DRef<double> bcdata, char (&filname)[80],
            int& iuns, int& irelv, FortranArray3DRef<double> snj0,
            FortranArray3DRef<double> snk0, FortranArray3DRef<double> sni0,
            int& ntime, FortranArray3DRef<double> snjm,
            FortranArray3DRef<double> snkm, FortranArray3DRef<double> snim,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid, int& nummem)
{
    // COMMON block aliases
    float& gamma    = cmn_fluid.gamma;
    float& gm1      = cmn_fluid.gm1;
    float& gp1      = cmn_fluid.gp1;
    float& gm1g     = cmn_fluid.gm1g;
    float& gp1g     = cmn_fluid.gp1g;
    float& ggm1     = cmn_fluid.ggm1;
    float& pr       = cmn_fluid2.pr;
    float& prt      = cmn_fluid2.prt;
    float& cbar     = cmn_fluid2.cbar;
    int32_t& ivmx   = cmn_maxiv.ivmx;
    int32_t& levt   = cmn_mgrd.levt;
    int32_t& kode   = cmn_mgrd.kode;
    int32_t& mode   = cmn_mgrd.mode;
    int32_t& ncyc   = cmn_mgrd.ncyc;
    int32_t& mtt    = cmn_mgrd.mtt;
    int32_t& icyc   = cmn_mgrd.icyc;
    int32_t& level  = cmn_mgrd.level;
    int32_t& lglobal= cmn_mgrd.lglobal;
    float& reue     = cmn_reyue.reue;
    float& tinf     = cmn_reyue.tinf;
    int32_t* ivisc  = cmn_reyue.ivisc;  // 1-based in Fortran, use ivisc[0..2]
    int32_t& isklton= cmn_sklton.isklton;
    int32_t* iwf    = cmn_wallfun.iwf;
    float& time     = cmn_unst.time;
    float& cfltau   = cmn_unst.cfltau;
    int32_t& ntstep = cmn_unst.ntstep;
    int32_t& ita    = cmn_unst.ita;
    int32_t& iunst  = cmn_unst.iunst;
    float& cfltau0  = cmn_unst.cfltau0;
    float& cfltauMax= cmn_unst.cfltaumax;
    float& radtodeg = cmn_conversion.radtodeg;
    float* tur10    = cmn_ivals.tur10;  // 1-based in Fortran, use tur10[l-1]

    // Local variables
    double twopi, thetajet;
    double vmag, rfreq, sideangj;
    double sxa, sya, sza, sxb, syb, szb;
    double uub, vvb, wwb;
    double sjetx1, sjety1, sjetz1;
    double snorma, snormb, dotp;
    double sjetx, sjety, sjetz;
    double t11, t12, t13, t14;
    int jdim1, kdim1, idim1;
    int jend1, kend1, iend1;
    int ipp, itrflg1, ip, ii, kk, jj, l;

    twopi = 2.0 * std::acos(-1.0);
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    jend1 = jend - 1;
    kend1 = kend - 1;
    iend1 = iend - 1;



    // ******************************************************************************
    //      j=1 boundary        sweeping jet                             bctype 2026
    // ******************************************************************************
    if (nface == 3) {
        // check to see if turbulence data is input (itrflg1 = 1) or
        // if freestream values are to be used (itrflg1 = 0)
        ipp     = 1;
        itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,10) > -1.e10f) itrflg1 = 1;

        for (ip = 1; ip <= 2; ip++) {
        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;

            for (int k = ksta; k <= kend1; k++) {
                kk = k - ksta + 1;
                vmag     = bcdata(kk,ii,ip,1);
                rfreq    = bcdata(kk,ii,ip,2);
                sideangj = bcdata(kk,ii,ip,3) / (double)radtodeg;
                sxa      = bcdata(kk,ii,ip,4);
                sya      = bcdata(kk,ii,ip,5);
                sza      = bcdata(kk,ii,ip,6);
                sxb      = bcdata(kk,ii,ip,7);
                syb      = bcdata(kk,ii,ip,8);
                szb      = bcdata(kk,ii,ip,9);

                // surface velocities
                uub = 0.;
                vvb = 0.;
                wwb = 0.;

                // for dynamic mesh, set velocity at wall to grid velocity at wall
                // if irelv > 0; otherwise, set to zero
                if (iuns > 0 && irelv > 0) {
                    uub = xtbj(k,i,1,1);
                    vvb = xtbj(k,i,2,1);
                    wwb = xtbj(k,i,3,1);
                }

                // find angle of sweeping jet relative to frame pointing up in z
                thetajet = sideangj * std::sin(twopi * rfreq * (double)time);
                sjetx1 = 0.;
                sjety1 = std::sin(thetajet);
                sjetz1 = std::cos(thetajet);

                // normalize vectors
                snorma = std::sqrt(sxa*sxa + sya*sya + sza*sza);
                sxa = sxa / snorma;
                sya = sya / snorma;
                sza = sza / snorma;
                snormb = std::sqrt(sxb*sxb + syb*syb + szb*szb);
                sxb = sxb / snormb;
                syb = syb / snormb;
                szb = szb / snormb;

                // check to make sure vectors are normal to each other
                dotp = sxa*sxb + sya*syb + sza*szb;
                if (std::abs(dotp) > 1.e-8) {
                    // error!
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " stopping in bc2026 - dot product not zero: %g", dotp);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }

                // rotate vector to required position and find u,v,w
                sjetx = sxb*sjetx1 + (sya*szb - sza*syb)*sjety1 + sxa*sjetz1;
                sjety = syb*sjetx1 + (sza*sxb - sxa*szb)*sjety1 + sya*sjetz1;
                sjetz = szb*sjetx1 + (sxa*syb - sya*sxb)*sjety1 + sza*sjetz1;
                uub = uub + sjetx*vmag;
                vvb = vvb + sjety*vmag;
                wwb = wwb + sjetz*vmag;

                qj0(k,i,1,ip) = q(1,k,i,1);
                qj0(k,i,2,ip) = uub;
                qj0(k,i,3,ip) = vvb;
                qj0(k,i,4,ip) = wwb;
                qj0(k,i,5,ip) = q(1,k,i,5);

                bcj(k,i,1) = 0.0;
            } // k loop
        } // i loop
        } // ip loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int k = ksta; k <= kend1; k++) {
                vj0(k,i,1,1) = vist3d(1,k,i);
                vj0(k,i,1,2) = vist3d(1,k,i);
            }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                    ii = i - ista + 1;
                    for (int k = ksta; k <= kend1; k++) {
                        kk = k - ksta + 1;
                        ip  = 1;
                        t11 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(kk,ii,ip,9+l);
                        ip  = 2;
                        t12 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(kk,ii,ip,9+l);
                        tj0(k,i,l,1) = t11;
                        tj0(k,i,l,2) = t12;
                    }
                }
                } // l loop
            }
        }

    } // end nface==3



    // ******************************************************************************
    //      j=jdim boundary     sweeping jet                             bctype 2026
    // ******************************************************************************
    if (nface == 4) {
        // check to see if turbulence data is input (itrflg1 = 1) or
        // if freestream values are to be used (itrflg1 = 0)
        ipp     = 1;
        itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,10) > -1.e10f) itrflg1 = 1;

        for (ip = 1; ip <= 2; ip++) {
        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;

            for (int k = ksta; k <= kend1; k++) {
                kk = k - ksta + 1;
                vmag     = bcdata(kk,ii,ip,1);
                rfreq    = bcdata(kk,ii,ip,2);
                sideangj = bcdata(kk,ii,ip,3) / (double)radtodeg;
                sxa      = bcdata(kk,ii,ip,4);
                sya      = bcdata(kk,ii,ip,5);
                sza      = bcdata(kk,ii,ip,6);
                sxb      = bcdata(kk,ii,ip,7);
                syb      = bcdata(kk,ii,ip,8);
                szb      = bcdata(kk,ii,ip,9);

                // surface velocities
                uub = 0.;
                vvb = 0.;
                wwb = 0.;

                // for dynamic mesh, set velocity at wall to grid velocity at wall
                // if irelv > 0; otherwise, set to zero
                if (iuns > 0 && irelv > 0) {
                    uub = xtbj(k,i,1,2);
                    vvb = xtbj(k,i,2,2);
                    wwb = xtbj(k,i,3,2);
                }

                // find angle of sweeping jet relative to frame pointing up in z
                thetajet = sideangj * std::sin(twopi * rfreq * (double)time);
                sjetx1 = 0.;
                sjety1 = std::sin(thetajet);
                sjetz1 = std::cos(thetajet);

                // normalize vectors
                snorma = std::sqrt(sxa*sxa + sya*sya + sza*sza);
                sxa = sxa / snorma;
                sya = sya / snorma;
                sza = sza / snorma;
                snormb = std::sqrt(sxb*sxb + syb*syb + szb*szb);
                sxb = sxb / snormb;
                syb = syb / snormb;
                szb = szb / snormb;

                // check to make sure vectors are normal to each other
                dotp = sxa*sxb + sya*syb + sza*szb;
                if (std::abs(dotp) > 1.e-8) {
                    // error!
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " stopping in bc2026 - dot product not zero: %g", dotp);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }

                // rotate vector to required position and find u,v,w
                sjetx = sxb*sjetx1 + (sya*szb - sza*syb)*sjety1 + sxa*sjetz1;
                sjety = syb*sjetx1 + (sza*sxb - sxa*szb)*sjety1 + sya*sjetz1;
                sjetz = szb*sjetx1 + (sxa*syb - sya*sxb)*sjety1 + sza*sjetz1;
                uub = uub + sjetx*vmag;
                vvb = vvb + sjety*vmag;
                wwb = wwb + sjetz*vmag;

                qj0(k,i,1,ip+2) = q(jdim1,k,i,1);
                qj0(k,i,2,ip+2) = uub;
                qj0(k,i,3,ip+2) = vvb;
                qj0(k,i,4,ip+2) = wwb;
                qj0(k,i,5,ip+2) = q(jdim1,k,i,5);

                bcj(k,i,2) = 0.0;
            } // k loop
        } // i loop
        } // ip loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int k = ksta; k <= kend1; k++) {
                vj0(k,i,1,3) = vist3d(jdim1,k,i);
                vj0(k,i,1,4) = vist3d(jdim1,k,i);
            }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                    ii = i - ista + 1;
                    for (int k = ksta; k <= kend1; k++) {
                        kk = k - ksta + 1;
                        ip  = 1;
                        t13 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(kk,ii,ip,9+l);
                        ip  = 2;
                        t14 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(kk,ii,ip,9+l);
                        tj0(k,i,l,3) = t13;
                        tj0(k,i,l,4) = t14;
                    }
                }
                } // l loop
            }
        }

    } // end nface==4



    // ******************************************************************************
    //      k=1 boundary        sweeping jet                             bctype 2026
    // ******************************************************************************
    if (nface == 5) {
        // check to see if turbulence data is input (itrflg1 = 1) or
        // if freestream values are to be used (itrflg1 = 0)
        ipp     = 1;
        itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,10) > -1.e10f) itrflg1 = 1;

        for (ip = 1; ip <= 2; ip++) {
        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;
                vmag     = bcdata(jj,ii,ip,1);
                rfreq    = bcdata(jj,ii,ip,2);
                sideangj = bcdata(jj,ii,ip,3) / (double)radtodeg;
                sxa      = bcdata(jj,ii,ip,4);
                sya      = bcdata(jj,ii,ip,5);
                sza      = bcdata(jj,ii,ip,6);
                sxb      = bcdata(jj,ii,ip,7);
                syb      = bcdata(jj,ii,ip,8);
                szb      = bcdata(jj,ii,ip,9);

                // surface velocities
                uub = 0.;
                vvb = 0.;
                wwb = 0.;

                // for dynamic mesh, set velocity at wall to grid velocity at wall
                // if irelv > 0; otherwise, set to zero
                if (iuns > 0 && irelv > 0) {
                    uub = xtbk(j,i,1,1);
                    vvb = xtbk(j,i,2,1);
                    wwb = xtbk(j,i,3,1);
                }

                // find angle of sweeping jet relative to frame pointing up in z
                thetajet = sideangj * std::sin(twopi * rfreq * (double)time);
                sjetx1 = 0.;
                sjety1 = std::sin(thetajet);
                sjetz1 = std::cos(thetajet);

                // normalize vectors
                snorma = std::sqrt(sxa*sxa + sya*sya + sza*sza);
                sxa = sxa / snorma;
                sya = sya / snorma;
                sza = sza / snorma;
                snormb = std::sqrt(sxb*sxb + syb*syb + szb*szb);
                sxb = sxb / snormb;
                syb = syb / snormb;
                szb = szb / snormb;

                // check to make sure vectors are normal to each other
                dotp = sxa*sxb + sya*syb + sza*szb;
                if (std::abs(dotp) > 1.e-8) {
                    // error!
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " stopping in bc2026 - dot product not zero: %g", dotp);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }

                // rotate vector to required position and find u,v,w
                sjetx = sxb*sjetx1 + (sya*szb - sza*syb)*sjety1 + sxa*sjetz1;
                sjety = syb*sjetx1 + (sza*sxb - sxa*szb)*sjety1 + sya*sjetz1;
                sjetz = szb*sjetx1 + (sxa*syb - sya*sxb)*sjety1 + sza*sjetz1;
                uub = uub + sjetx*vmag;
                vvb = vvb + sjety*vmag;
                wwb = wwb + sjetz*vmag;

                qk0(j,i,1,ip) = q(j,1,i,1);
                qk0(j,i,2,ip) = uub;
                qk0(j,i,3,ip) = vvb;
                qk0(j,i,4,ip) = wwb;
                qk0(j,i,5,ip) = q(j,1,i,5);

                bck(j,i,1) = 0.0;
            } // j loop
        } // i loop
        } // ip loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int j = jsta; j <= jend1; j++) {
                vk0(j,i,1,1) = vist3d(j,1,i);
                vk0(j,i,1,2) = vist3d(j,1,i);
            }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                    ii = i - ista + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        jj = j - jsta + 1;
                        ip  = 1;
                        t11 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,ii,ip,9+l);
                        ip  = 2;
                        t12 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,ii,ip,9+l);
                        tk0(j,i,l,1) = t11;
                        tk0(j,i,l,2) = t12;
                    }
                }
                } // l loop
            }
        }

    } // end nface==5



    // ******************************************************************************
    //      k=kdim boundary     sweeping jet                             bctype 2026
    // ******************************************************************************
    if (nface == 6) {
        // check to see if turbulence data is input (itrflg1 = 1) or
        // if freestream values are to be used (itrflg1 = 0)
        ipp     = 1;
        itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,10) > -1.e10f) itrflg1 = 1;

        for (ip = 1; ip <= 2; ip++) {
        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;
                vmag     = bcdata(jj,ii,ip,1);
                rfreq    = bcdata(jj,ii,ip,2);
                sideangj = bcdata(jj,ii,ip,3) / (double)radtodeg;
                sxa      = bcdata(jj,ii,ip,4);
                sya      = bcdata(jj,ii,ip,5);
                sza      = bcdata(jj,ii,ip,6);
                sxb      = bcdata(jj,ii,ip,7);
                syb      = bcdata(jj,ii,ip,8);
                szb      = bcdata(jj,ii,ip,9);

                // surface velocities
                uub = 0.;
                vvb = 0.;
                wwb = 0.;

                // for dynamic mesh, set velocity at wall to grid velocity at wall
                // if irelv > 0; otherwise, set to zero
                if (iuns > 0 && irelv > 0) {
                    uub = xtbk(j,i,1,2);
                    vvb = xtbk(j,i,2,2);
                    wwb = xtbk(j,i,3,2);
                }

                // find angle of sweeping jet relative to frame pointing up in z
                thetajet = sideangj * std::sin(twopi * rfreq * (double)time);
                sjetx1 = 0.;
                sjety1 = std::sin(thetajet);
                sjetz1 = std::cos(thetajet);

                // normalize vectors
                snorma = std::sqrt(sxa*sxa + sya*sya + sza*sza);
                sxa = sxa / snorma;
                sya = sya / snorma;
                sza = sza / snorma;
                snormb = std::sqrt(sxb*sxb + syb*syb + szb*szb);
                sxb = sxb / snormb;
                syb = syb / snormb;
                szb = szb / snormb;

                // check to make sure vectors are normal to each other
                dotp = sxa*sxb + sya*syb + sza*szb;
                if (std::abs(dotp) > 1.e-8) {
                    // error!
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " stopping in bc2026 - dot product not zero: %g", dotp);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }

                // rotate vector to required position and find u,v,w
                sjetx = sxb*sjetx1 + (sya*szb - sza*syb)*sjety1 + sxa*sjetz1;
                sjety = syb*sjetx1 + (sza*sxb - sxa*szb)*sjety1 + sya*sjetz1;
                sjetz = szb*sjetx1 + (sxa*syb - sya*sxb)*sjety1 + sza*sjetz1;
                uub = uub + sjetx*vmag;
                vvb = vvb + sjety*vmag;
                wwb = wwb + sjetz*vmag;

                qk0(j,i,1,ip+2) = q(j,kdim1,i,1);
                qk0(j,i,2,ip+2) = uub;
                qk0(j,i,3,ip+2) = vvb;
                qk0(j,i,4,ip+2) = wwb;
                qk0(j,i,5,ip+2) = q(j,kdim1,i,5);

                bck(j,i,2) = 0.0;
            } // j loop
        } // i loop
        } // ip loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
            for (int j = jsta; j <= jend1; j++) {
                vk0(j,i,1,3) = vist3d(j,kdim1,i);
                vk0(j,i,1,4) = vist3d(j,kdim1,i);
            }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                for (int i = ista; i <= iend1; i++) {
                    ii = i - ista + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        jj = j - jsta + 1;
                        ip  = 1;
                        t13 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,ii,ip,9+l);
                        ip  = 2;
                        t14 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,ii,ip,9+l);
                        tk0(j,i,l,3) = t13;
                        tk0(j,i,l,4) = t14;
                    }
                }
                } // l loop
            }
        }

    } // end nface==6



    // ******************************************************************************
    //      i=1 boundary        sweeping jet                             bctype 2026
    // ******************************************************************************
    if (nface == 1) {
        // check to see if turbulence data is input (itrflg1 = 1) or
        // if freestream values are to be used (itrflg1 = 0)
        ipp     = 1;
        itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,10) > -1.e10f) itrflg1 = 1;

        for (ip = 1; ip <= 2; ip++) {
        for (int k = ksta; k <= kend1; k++) {
            kk = k - ksta + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;
                vmag     = bcdata(jj,kk,ip,1);
                rfreq    = bcdata(jj,kk,ip,2);
                sideangj = bcdata(jj,kk,ip,3) / (double)radtodeg;
                sxa      = bcdata(jj,kk,ip,4);
                sya      = bcdata(jj,kk,ip,5);
                sza      = bcdata(jj,kk,ip,6);
                sxb      = bcdata(jj,kk,ip,7);
                syb      = bcdata(jj,kk,ip,8);
                szb      = bcdata(jj,kk,ip,9);

                // surface velocities
                uub = 0.;
                vvb = 0.;
                wwb = 0.;

                // for dynamic mesh, set velocity at wall to grid velocity at wall
                // if irelv > 0; otherwise, set to zero
                if (iuns > 0 && irelv > 0) {
                    uub = xtbi(j,k,1,1);
                    vvb = xtbi(j,k,2,1);
                    wwb = xtbi(j,k,3,1);
                }

                // find angle of sweeping jet relative to frame pointing up in z
                thetajet = sideangj * std::sin(twopi * rfreq * (double)time);
                sjetx1 = 0.;
                sjety1 = std::sin(thetajet);
                sjetz1 = std::cos(thetajet);

                // normalize vectors
                snorma = std::sqrt(sxa*sxa + sya*sya + sza*sza);
                sxa = sxa / snorma;
                sya = sya / snorma;
                sza = sza / snorma;
                snormb = std::sqrt(sxb*sxb + syb*syb + szb*szb);
                sxb = sxb / snormb;
                syb = syb / snormb;
                szb = szb / snormb;

                // check to make sure vectors are normal to each other
                dotp = sxa*sxb + sya*syb + sza*szb;
                if (std::abs(dotp) > 1.e-8) {
                    // error!
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " stopping in bc2026 - dot product not zero: %g", dotp);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }

                // rotate vector to required position and find u,v,w
                sjetx = sxb*sjetx1 + (sya*szb - sza*syb)*sjety1 + sxa*sjetz1;
                sjety = syb*sjetx1 + (sza*sxb - sxa*szb)*sjety1 + sya*sjetz1;
                sjetz = szb*sjetx1 + (sxa*syb - sya*sxb)*sjety1 + sza*sjetz1;
                uub = uub + sjetx*vmag;
                vvb = vvb + sjety*vmag;
                wwb = wwb + sjetz*vmag;

                qi0(j,k,1,ip) = q(j,k,1,1);
                qi0(j,k,2,ip) = uub;
                qi0(j,k,3,ip) = vvb;
                qi0(j,k,4,ip) = wwb;
                qi0(j,k,5,ip) = q(j,k,1,5);

                bci(j,k,1) = 0.0;
            } // j loop
        } // k loop
        } // ip loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
            for (int j = jsta; j <= jend1; j++) {
                vi0(j,k,1,1) = vist3d(j,k,1);
                vi0(j,k,1,2) = vist3d(j,k,1);
            }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    kk = k - ksta + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        jj = j - jsta + 1;
                        ip  = 1;
                        t11 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,kk,ip,9+l);
                        ip  = 2;
                        t12 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,kk,ip,9+l);
                        ti0(j,k,l,1) = t11;
                        ti0(j,k,l,2) = t12;
                    }
                }
                } // l loop
            }
        }

    } // end nface==1



    // ******************************************************************************
    //      i=idim boundary     sweeping jet                             bctype 2026
    // ******************************************************************************
    if (nface == 2) {
        // check to see if turbulence data is input (itrflg1 = 1) or
        // if freestream values are to be used (itrflg1 = 0)
        ipp     = 1;
        itrflg1 = 0;
        if ((float)bcdata(1,1,ipp,10) > -1.e10f) itrflg1 = 1;

        for (ip = 1; ip <= 2; ip++) {
        for (int k = ksta; k <= kend1; k++) {
            kk = k - ksta + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;
                vmag     = bcdata(jj,kk,ip,1);
                rfreq    = bcdata(jj,kk,ip,2);
                sideangj = bcdata(jj,kk,ip,3) / (double)radtodeg;
                sxa      = bcdata(jj,kk,ip,4);
                sya      = bcdata(jj,kk,ip,5);
                sza      = bcdata(jj,kk,ip,6);
                sxb      = bcdata(jj,kk,ip,7);
                syb      = bcdata(jj,kk,ip,8);
                szb      = bcdata(jj,kk,ip,9);

                // surface velocities
                uub = 0.;
                vvb = 0.;
                wwb = 0.;

                // for dynamic mesh, set velocity at wall to grid velocity at wall
                // if irelv > 0; otherwise, set to zero
                if (iuns > 0 && irelv > 0) {
                    uub = xtbi(j,k,1,2);
                    vvb = xtbi(j,k,2,2);
                    wwb = xtbi(j,k,3,2);
                }

                // find angle of sweeping jet relative to frame pointing up in z
                thetajet = sideangj * std::sin(twopi * rfreq * (double)time);
                sjetx1 = 0.;
                sjety1 = std::sin(thetajet);
                sjetz1 = std::cos(thetajet);

                // normalize vectors
                snorma = std::sqrt(sxa*sxa + sya*sya + sza*sza);
                sxa = sxa / snorma;
                sya = sya / snorma;
                sza = sza / snorma;
                snormb = std::sqrt(sxb*sxb + syb*syb + szb*szb);
                sxb = sxb / snormb;
                syb = syb / snormb;
                szb = szb / snormb;

                // check to make sure vectors are normal to each other
                dotp = sxa*sxb + sya*syb + sza*szb;
                if (std::abs(dotp) > 1.e-8) {
                    // error!
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, " stopping in bc2026 - dot product not zero: %g", dotp);
                    int m1 = -1;
                    termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                }

                // rotate vector to required position and find u,v,w
                sjetx = sxb*sjetx1 + (sya*szb - sza*syb)*sjety1 + sxa*sjetz1;
                sjety = syb*sjetx1 + (sza*sxb - sxa*szb)*sjety1 + sya*sjetz1;
                sjetz = szb*sjetx1 + (sxa*syb - sya*sxb)*sjety1 + sza*sjetz1;
                uub = uub + sjetx*vmag;
                vvb = vvb + sjety*vmag;
                wwb = wwb + sjetz*vmag;

                qi0(j,k,1,ip+2) = q(j,k,idim1,1);
                qi0(j,k,2,ip+2) = uub;
                qi0(j,k,3,ip+2) = vvb;
                qi0(j,k,4,ip+2) = wwb;
                qi0(j,k,5,ip+2) = q(j,k,idim1,5);

                bci(j,k,2) = 0.0;
            } // j loop
        } // k loop
        } // ip loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
            for (int j = jsta; j <= jend1; j++) {
                vi0(j,k,1,3) = vist3d(j,k,idim1);
                vi0(j,k,1,4) = vist3d(j,k,idim1);
            }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    kk = k - ksta + 1;
                    for (int j = jsta; j <= jend1; j++) {
                        jj = j - jsta + 1;
                        ip  = 1;
                        t13 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,kk,ip,9+l);
                        ip  = 2;
                        t14 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,kk,ip,9+l);
                        ti0(j,k,l,3) = t13;
                        ti0(j,k,l,4) = t14;
                    }
                }
                } // l loop
            }
        }

    } // end nface==2

} // end bc2026()



void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
        int& nwork, double& cl, FortranArray1DRef<int> nou,
        FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
        int& maxbl, int& maxgr, int& maxseg,
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
    int iuns  = cmn_unst.iunst;

    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        // Build array views from w using lw offsets
        FortranArray4DRef<double> q    (&w(lw(1,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray4DRef<double> qj0  (&w(lw(2,ibl)),  kdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qk0  (&w(lw(3,ibl)),  jdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qi0  (&w(lw(4,ibl)),  jdim,   kdim,   5,      4);
        FortranArray4DRef<double> sj   (&w(lw(5,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> sk   (&w(lw(6,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> si   (&w(lw(7,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray3DRef<double> bcj  (&w(lw(8,ibl)),  kdim,   idim-1, 2);
        FortranArray3DRef<double> bck  (&w(lw(9,ibl)),  jdim,   idim-1, 2);
        FortranArray3DRef<double> bci  (&w(lw(10,ibl)), jdim,   kdim,   2);
        FortranArray4DRef<double> xtbj (&w(lw(11,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbk (&w(lw(12,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> xtbi (&w(lw(13,ibl)), jdim,   kdim,   3,      2);
        FortranArray4DRef<double> atbj (&w(lw(14,ibl)), kdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbk (&w(lw(15,ibl)), jdim,   idim-1, 3,      2);
        FortranArray4DRef<double> atbi (&w(lw(16,ibl)), jdim,   kdim,   3,      2);

        // Build array views from wk using lw2 offsets
        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&wk(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&wk(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&wk(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&wk(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&wk(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&wk(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // smin, snj0, snk0, sni0, snjm, snkm, snim from lw2(9..15)
        FortranArray3DRef<double> smin  (&wk(lw2(9,ibl)),  jdim-1, kdim-1, idim-1);
        FortranArray3DRef<double> snj0  (&wk(lw2(10,ibl)), jdim-1, kdim-1, idim-1);
        FortranArray3DRef<double> snk0  (&wk(lw2(11,ibl)), jdim-1, kdim-1, idim-1);
        FortranArray3DRef<double> sni0  (&wk(lw2(12,ibl)), jdim-1, kdim-1, idim-1);
        FortranArray3DRef<double> snjm  (&wk(lw2(13,ibl)), jdim-1, kdim-1, idim-1);
        FortranArray3DRef<double> snkm  (&wk(lw2(14,ibl)), jdim-1, kdim-1, idim-1);
        FortranArray3DRef<double> snim  (&wk(lw2(15,ibl)), jdim-1, kdim-1, idim-1);

        // irelv: set based on itrans/irotat/idefrm flags
        int irelv = 0;
        if (itrans(ibl) > 0 || irotat(ibl) > 0 || idefrm(ibl) > 0) irelv = 1;

        // Loop over j-face BC segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) == 2026) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                // bcdata from wk at offset given by bcfilej
                // bcdata from wk at offset given by bcfilej
                int ibcf  = bcfilej(ibl,iseg,1);
                // filname from bcfiles
                char filname[80];
                if (ibcf >= 1 && ibcf <= mxbcfil) {
                    std::memcpy(filname, bcfiles(ibcf), 80);
                } else {
                    std::memset(filname, ' ', 80);
                }
                FortranArray4DRef<double> bcdata(&wk(lw2(9,ibl)), mdim, ndim, 2, 12);
                bc2026(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, smin, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns, irelv,
                       snj0, snk0, sni0, ntime, snjm, snkm, snim,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // Loop over k-face BC segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 2026) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                int ibcf  = bcfilek(ibl,iseg,1);
                char filname[80];
                if (ibcf >= 1 && ibcf <= mxbcfil) {
                    std::memcpy(filname, bcfiles(ibcf), 80);
                } else {
                    std::memset(filname, ' ', 80);
                }
                FortranArray4DRef<double> bcdata(&wk(lw2(9,ibl)), mdim, ndim, 2, 12);
                bc2026(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, smin, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns, irelv,
                       snj0, snk0, sni0, ntime, snjm, snkm, snim,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // Loop over i-face BC segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 2026) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                int mdim  = lwdat(ibl,iseg,3);
                int ndim  = lwdat(ibl,iseg,4);
                int ibcf  = bcfilei(ibl,iseg,1);
                char filname[80];
                if (ibcf >= 1 && ibcf <= mxbcfil) {
                    std::memcpy(filname, bcfiles(ibcf), 80);
                } else {
                    std::memset(filname, ' ', 80);
                }
                FortranArray4DRef<double> bcdata(&wk(lw2(9,ibl)), mdim, ndim, 2, 12);
                bc2026(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci, xtbj, xtbk, xtbi, atbj, atbk, atbi,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, smin, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns, irelv,
                       snj0, snk0, sni0, ntime, snjm, snkm, snim,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

    } // end ibl loop

} // end bc()

} // namespace bc2026_ns
