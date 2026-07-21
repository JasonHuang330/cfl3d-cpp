// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// avgint.cpp — faithful C++ translation of avgint.F (CFL3D ronnie).
// Create interpolation data on coarser levels by averaging finer-level data.
// All arrays are 1-based flat (nptc / nptf), so element (ll) -> (ll-1).
#include "ron_common.h"

void avgint(double* xiec,double* etac,int* nblkc,int nptc,
            double* xief,double* etaf,int* nblkf,int nptf,
            int j1c,int j2c,int k1c,int k2c,int j1f,int j2f,int k1f,int k2f)
{
    (void)nptc; (void)nptf;
    for (int kc = k1c; kc <= k2c-1; ++kc) {
        int kf  = (kc-1)*2 + 1;
        int kfp = kf + 1;
        kfp = i_min(kfp, k2f-1);
        for (int jc = j1c; jc <= j2c-1; ++jc) {
            int jf  = (jc-1)*2 + 1;
            int jfp = jf + 1;
            jfp = i_min(jfp, j2f-1);
            int llc     = (j2c-j1c)*(kc-k1c)  + (jc -j1c+1);
            int llf     = (j2f-j1f)*(kf-k1f)  + (jf -j1f+1);
            int llfjp   = (j2f-j1f)*(kf-k1f)  + (jfp-j1f+1);
            int llfkp   = (j2f-j1f)*(kfp-k1f) + (jf -j1f+1);
            int llfjpkp = (j2f-j1f)*(kfp-k1f) + (jfp-j1f+1);

            // all finer points must lie in same block for averaging;
            // otherwise, use only data from lower corner point
            if (nblkf[llfjp-1]   != nblkf[llf-1]) llfjp   = llf;
            if (nblkf[llfkp-1]   != nblkf[llf-1]) llfkp   = llf;
            if (nblkf[llfjpkp-1] != nblkf[llf-1]) llfjpkp = llf;

            double xiefavg = 0.25*(xief[llf-1]   + xief[llfjp-1]
                                 + xief[llfkp-1] + xief[llfjpkp-1]);
            double etafavg = 0.25*(etaf[llf-1]   + etaf[llfjp-1]
                                 + etaf[llfkp-1] + etaf[llfjpkp-1]);
            int    jfhat   = f_int(xiefavg);
            double xiefhat = xiefavg - jfhat;
            int    kfhat   = f_int(etafavg);
            double etafhat = etafavg - kfhat;
            int    jchat   = (jfhat+1)/2;
            int    kchat   = (kfhat+1)/2;
            double xiechat, etachat;
            if (j1f == 1 && j2f == 2) {
                // 2d case
                xiechat = xiefhat;
            } else {
                if (jfhat/2*2 == jfhat) {
                    xiechat = 0.5*(1.0+xiefhat);
                } else {
                    xiechat = 0.5*xiefhat;
                }
            }
            if (k1f == 1 && k2f == 2) {
                // 2d case
                etachat = etafhat;
            } else {
                if (kfhat/2*2 == kfhat) {
                    etachat = 0.5*(1.0+etafhat);
                } else {
                    etachat = 0.5*etafhat;
                }
            }
            xiec[llc-1]  = jchat + xiechat;
            etac[llc-1]  = kchat + etachat;
            nblkc[llc-1] = nblkf[llf-1];
        }
    }
}
