// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "moddefl.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace moddefl_ns {

// Inline helpers for CFL3D non-standard intrinsics
static inline double ccmin(double a, double b) { return std::min(a, b); }
static inline double ccabs(double x) { return std::fabs(x); }

void moddefl(FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn,
             FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> perturb,
             double& cmyt, double& cnwt, double& xorig, double& yorig, double& zorig,
             int& maxaes, int& nmds, int& irbtrim, int& maxbl, int& myid)
{
    // COMMON block aliases
    int32_t& ndefrm   = cmn_elastic.ndefrm;
    int32_t& naesrf   = cmn_elastic.naesrf;
    int32_t& iexp     = cmn_zero.iexp;
    float&   time     = cmn_unst.time;
    float&   dt       = cmn_info.dt;
    int32_t& ntt      = cmn_info.ntt;
    float&   xmach    = cmn_info.xmach;
    float&   tmass    = cmn_rbstmt2.tmass;
    float&   qinfrb   = cmn_rbstmt2.qinfrb;
    float&   greflrb  = cmn_rbstmt2.greflrb;
    float&   gaccel   = cmn_rbstmt2.gaccel;
    float&   crefrb   = cmn_rbstmt2.crefrb;
    float&   areat    = cmn_rbstmt2.areat;
    float&   dmtrmn   = cmn_trim.dmtrmn;
    float&   dmtrmnm  = cmn_trim.dmtrmnm;
    float&   dlcln    = cmn_trim.dlcln;
    float&   dlclnm   = cmn_trim.dlclnm;
    float&   trtol    = cmn_trim.trtol;
    float&   cmy      = cmn_trim.cmy;
    float&   cnw      = cmn_trim.cnw;
    float&   alf0     = cmn_trim.alf0;
    float&   alf1     = cmn_trim.alf1;
    float&   dtrmsmx  = cmn_trim.dtrmsmx;
    float&   dtrmsmn  = cmn_trim.dtrmsmn;
    float&   dalfmx   = cmn_trim.dalfmx;
    float&   ddtmx    = cmn_trim.ddtmx;
    float&   ddtrm0   = cmn_trim.ddtrm0;
    float&   ddtrm1   = cmn_trim.ddtrm1;
    int32_t& itrmt    = cmn_trim.itrmt;
    int32_t& itrminc  = cmn_trim.itrminc;
    float*   fp_flat  = cmn_trim.fp;   // fp(4,4) column-major float[16]
    float*   tp_flat  = cmn_trim.tp;   // tp(4,4) column-major float[16]
    float&   relax    = cmn_trim.relax;
    float&   r33      = cmn_trim1.r33;
    float&   r44      = cmn_trim1.r44;
    int32_t& ittr     = cmn_trim1.ittr;
    float&   a11      = cmn_trim1.a11;
    float&   a12      = cmn_trim1.a12;
    float&   a22      = cmn_trim1.a22;
    float&   r11      = cmn_trim1.r11;
    float&   r22      = cmn_trim1.r22;
    float&   xmc0     = cmn_motionmc.xmc0;

    // fp(i,j) column-major: fp(i,j) = fp_flat[(j-1)*4 + (i-1)]
    // fp(1,1)=fp_flat[0], fp(1,2)=fp_flat[4], fp(2,1)=fp_flat[1], fp(2,2)=fp_flat[5]

    double dtimemx = -std::log(std::pow(10.0, (double)(-iexp)));

    for (int iaes = 1; iaes <= naesrf; iaes++) {

        int    nmodes = (int)aesrfdat(5, iaes);
        double grefl  = aesrfdat(2, iaes);
        double uinf   = aesrfdat(3, iaes);
        double qinf   = aesrfdat(4, iaes);
        double ainf   = uinf / (double)xmach;
        double times  = (double)time * grefl / ainf;
        double dts    = (double)dt * grefl / ainf;

        for (int nm = 1; nm <= nmodes; nm++) {

            int moddfl = (int)perturb(nm, iaes, 1);

            // fixed mode
            if (moddfl < 0) {
                xs(2*nm-1, iaes) = 0.;
                xs(2*nm,   iaes) = 0.;
            }

            // harmonic modal oscillation
            if (moddfl == 1) {
                double amp   = perturb(nm, iaes, 2);
                double freqp = perturb(nm, iaes, 3);
                double t0    = perturb(nm, iaes, 4);
                xs(2*nm-1, iaes) = amp * std::sin(freqp * (times - t0));
                xs(2*nm,   iaes) = amp * freqp * std::cos(freqp * (times - t0));
            }

            // Gaussian pulse modal deflection (half is the pulse half-life)
            if (moddfl == 2) {
                double amp     = perturb(nm, iaes, 2);
                double half    = perturb(nm, iaes, 3);
                double t0      = perturb(nm, iaes, 4);
                double constv  = std::log(2.0) / (half * half);
                double dtime   = constv * (times - t0) * (times - t0);
                dtime          = ccmin(dtime, dtimemx);
                double expterm = std::exp(-dtime);
                xs(2*nm-1, iaes) = amp * expterm;
                xs(2*nm,   iaes) = -2. * constv * (times - t0) *
                                   xs(2*nm-1, iaes);
            }

            // step pulse modal deflection
            if (moddfl == 3) {
                double amp = perturb(nm, iaes, 2);
                double t0  = perturb(nm, iaes, 4);
                if ((float)times < (float)(t0 - dts/2.)) {
                    xs(2*nm-1, iaes) = 0.;
                    xs(2*nm,   iaes) = 0.;
                } else if ((float)times > (float)(t0 - dts/2.) &&
                           (float)times < (float)(t0 + dts/2.)) {
                    xs(2*nm-1, iaes) = amp;
                    xs(2*nm,   iaes) = amp / dts;
                } else {
                    xs(2*nm-1, iaes) = amp;
                    xs(2*nm,   iaes) = 0.;
                }
            }

            // rigid body mode
            if (moddfl == 4) {

                if (ntt / itrminc * itrminc == ntt) {
                    if ((std::abs((float)dlcln) < (float)trtol) &&
                        (std::abs((float)dmtrmn) < (float)trtol) &&
                        (itrmt > 0)) return;
                    itrmt  = itrmt + 1;
                    alf0   = alf1;
                    ddtrm0 = (float)xs(2*nm-1, iaes);
//                 call sqrtcumm(itrmt,dcl,dcm,ddclda,ddcmda,da,dd
//    .                         ,dlcln,dmtrmn,tp,fp,a11,a12,a22
//    .                         ,r11,r22,r33,r44,alf1,alf0,ddtrm1
//    .                         ,ddtrm0)

                    // Level 1 g flight
                    dmtrmnm = dmtrmn;
                    qinfrb  = (float)qinf;

                    // assumes ialph = 0
                    // crefrb was set equal to cref in init_rb
                    cmy    = (float)(cmyt + cnwt * (xorig - (double)xmc0) / (double)crefrb);
                    dmtrmn = cmy;
                    cnw    = (float)cnwt;
                    dlclnm = dlcln;
                    dlcln  = (float)(cnwt - ((double)tmass * (double)gaccel) /
                                     ((double)qinfrb * (double)areat * grefl * grefl));
                    if ((std::abs((float)dlcln) < (float)trtol) &&
                        (std::abs((float)dmtrmn) < (float)trtol)) return;

                    // fp(1,1)=fp_flat[0], fp(1,2)=fp_flat[4]
                    // fp(2,1)=fp_flat[1], fp(2,2)=fp_flat[5]
                    alf1   = (float)((double)alf0  + (double)relax *
                                     (-(double)fp_flat[0] * (double)dlcln
                                      -(double)fp_flat[4] * (double)dmtrmn));
                    ddtrm1 = (float)((double)ddtrm0 + (double)relax *
                                     (-(double)fp_flat[1] * (double)dlcln
                                      -(double)fp_flat[5] * (double)dmtrmn));
                    if ((float)ddtrm1 > (float)dtrmsmx) ddtrm1 = dtrmsmx;
                    if ((float)ddtrm1 < (float)dtrmsmn) ddtrm1 = dtrmsmn;
                    if (std::abs((float)(ddtrm1 - ddtrm0)) > (float)ddtmx) {
                        ddtrm1 = (float)((double)ddtrm0 + (double)ddtmx *
                                         ((double)ddtrm1 - (double)ddtrm0) /
                                         ccabs((double)ddtrm1 - (double)ddtrm0));
                    }
                    xs(2*nm-1, iaes) = (double)ddtrm1;
                    if (std::abs((float)(alf1 - alf0)) > (float)dalfmx) {
                        alf1 = (float)((double)alf0 + (double)dalfmx *
                                       ((double)alf1 - (double)alf0) /
                                       ccabs((double)alf1 - (double)alf0));
                    }
                    if (myid == 0) {
                        fortran_write_unit(79 + myid,
                            "%8d%8d %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f\n",
                            (int)itrmt, (int)ntt,
                            (double)cmyt, (double)cmy, (double)cnwt, (double)dlcln,
                            (double)ddtrm1 * 180. / 3.14159,
                            (double)alf1   * 180. / 3.14159,
                            (double)xorig, (double)xmc0, (double)crefrb);
                    }
                } else {
                    xs(2*nm-1, iaes) = xxn(2*nm-1, iaes);
                }
            }

        } // end do nm
    } // end do iaes

    return;
}

void sqrtcumm(int& ittr,
              FortranArray1DRef<double> dcl, FortranArray1DRef<double> dcm,
              FortranArray1DRef<double> ddclda, FortranArray1DRef<double> ddcmda,
              FortranArray1DRef<double> da, FortranArray1DRef<double> dd,
              double& dlcln, double& dmtrmn,
              FortranArray2DRef<double> tp, FortranArray2DRef<double> fp,
              double& a11, double& a12, double& a22,
              double& r11, double& r22, double& r33, double& r44,
              double& alf1, double& alf0, double& ddtrm1, double& ddtrm0)
{
    dcl(ittr) = dlcln;
    dcm(ittr) = dmtrmn;
    if (ittr == 1) {
        da(ittr) = 0.;
        dd(ittr) = 0.;
    } else {
        da(ittr)     = alf1 - alf0;
        dd(ittr)     = ddtrm1 - ddtrm0;
        ddclda(ittr) = (dcl(ittr) - dcl(ittr-1));
        ddcmda(ittr) = (dcm(ittr) - dcm(ittr-1));
        double a1 = da(ittr-1) * da(ittr-1);
        double a2 = da(ittr-1) * dd(ittr-1);
        double a3 = dd(ittr-1) * dd(ittr-1);
        double r1 = ddclda(ittr) * da(ittr-1);
        double r2 = ddclda(ittr) * dd(ittr-1);
        double r3 = ddcmda(ittr) * da(ittr-1);
        double r4 = ddcmda(ittr) * dd(ittr-1);
        a11 = a11 + a1;
        a12 = a12 + a2;
        a22 = a22 + a3;
        r11 = r11 + r1;
        r22 = r22 + r2;
        r33 = r33 + r3;
        r44 = r44 + r4;
        double dtrm = a11 * a22 - a12 * a12;
        if (ittr > 2) {
            if (dtrm != 0.) {
                tp(1,1) = ( a22*r11 - a12*r22) / dtrm;
                tp(1,2) = (-a12*r11 + a11*r22) / dtrm;
                tp(2,1) = ( a22*r33 - a12*r44) / dtrm;
                tp(2,2) = (-a12*r33 + a11*r44) / dtrm;
            }
            double dtr = tp(1,1)*tp(2,2) - tp(1,2)*tp(2,1);
            if (dtr != 0.) {
                fp(1,1) =  tp(2,2) / dtr;
                fp(1,2) = -tp(1,2) / dtr;
                fp(2,1) = -tp(2,1) / dtr;
                fp(2,2) =  tp(1,1) / dtr;
            }
        }
    }
    return;
}

} // namespace moddefl_ns
