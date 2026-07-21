// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// arc.cpp — faithful C++ translation of arc.F (CFL3D ronnie).
// Arc-length correction to generalized coordinates near a boundary when the
// shearing correction has failed.
//
//   x1/y1/z1(jdim1,kdim1,msub1)   at (j,k,l) -> (j-1)+(k-1)*jdim1+(l-1)*jdim1*kdim1
//   x2/y2/z2(jdim1,kdim1,msub2)   at (j,k,l) -> (j-1)+(k-1)*jdim1+(l-1)*jdim1*kdim1
//   sxie/seta(jdim1,kdim1,msub1)  same 3-index layout
//   sxie2/seta2(jdim1,kdim1,msub2) same 3-index layout
//   xie2s/eta2s(jdim1,kdim1)      at (j,k)   -> (j-1)+(k-1)*jdim1
//   xie2/eta2/mblkpt(npt)         1-based flat, (ll) -> (ll-1)
//   jjmax1/kkmax1(msub1), jjmax2/kkmax2(msub2)
#include "ron_common.h"

void arc(int jdim1,int kdim1,int msub1,int msub2,
         int* jjmax1,int* kkmax1,int lmax1,double* x1,double* y1,double* z1,
         int limit0,int* jjmax2,int* kkmax2,
         double* x2,double* y2,double* z2,double* xie2,double* eta2,int* mblkpt,
         int ifit,int itmax,int jcorr,int kcorr,
         double* sxie,double* seta,double* sxie2,double* seta2,
         double* xie2s,double* eta2s,
         int intmx,int icheck,int* nblkj,int* nblkk,int* jmm,int* kmm,
         int j21,int j22,int k21,int k22,
         int npt,int* xif1,int* xif2,int* etf1,int* etf2,
         int* nou,FStr* bou,int nbuf,int ibufdim,int* mblk2nd,int maxbl)
{
    (void)msub2; (void)ifit; (void)intmx; (void)icheck; (void)nblkj; (void)nblkk;
    (void)npt; (void)nbuf;
    int myid = 0;   // rechk/topol2 sequential build

    auto I3  = [&](int j,int k,int l){
        return (long)(j-1) + (long)(k-1)*jdim1 + (long)(l-1)*jdim1*kdim1; };
    auto I2  = [&](int j,int k){ return (long)(j-1) + (long)(k-1)*jdim1; };

    char line[160];

    int l2    = 1;
    int jmax2 = jjmax2[l2-1];
    int kmax2 = kkmax2[l2-1];

    // arc lengths for (unexpanded) "from" grid(s)
    for (int l = 1; l <= lmax1; ++l) {
        int jmax1 = jjmax1[l-1];
        int kmax1 = kkmax1[l-1];
        for (int j = 2; j <= jmax1-1; ++j) {
            seta[I3(j-1,1,l)] = 0.0;
            for (int k = 3; k <= kmax1-1; ++k) {
                double dx = x1[I3(j,k,l)]-x1[I3(j,k-1,l)];
                double dy = y1[I3(j,k,l)]-y1[I3(j,k-1,l)];
                double dz = z1[I3(j,k,l)]-z1[I3(j,k-1,l)];
                seta[I3(j-1,k-1,l)] = seta[I3(j-1,k-2,l)]
                                    + std::sqrt(dx*dx+dy*dy+dz*dz);
            }
        }
        for (int k = 2; k <= kmax1-1; ++k) {
            sxie[I3(1,k-1,l)] = 0.0;
            for (int j = 3; j <= jmax1-1; ++j) {
                double dx = x1[I3(j,k,l)]-x1[I3(j-1,k,l)];
                double dy = y1[I3(j,k,l)]-y1[I3(j-1,k,l)];
                double dz = z1[I3(j,k,l)]-z1[I3(j-1,k,l)];
                sxie[I3(j-1,k-1,l)] = sxie[I3(j-2,k-1,l)]
                                    + std::sqrt(dx*dx+dy*dy+dz*dz);
            }
        }
    }

    // arc lengths for "to" grid
    for (int j = 1; j <= jmax2; ++j) {
        seta2[I3(j,1,1)] = 0.0;
        for (int k = 2; k <= kmax2; ++k) {
            double dx = x2[I3(j,k,1)]-x2[I3(j,k-1,1)];
            double dy = y2[I3(j,k,1)]-y2[I3(j,k-1,1)];
            double dz = z2[I3(j,k,1)]-z2[I3(j,k-1,1)];
            seta2[I3(j,k,1)] = seta2[I3(j,k-1,1)] + std::sqrt(dx*dx+dy*dy+dz*dz);
        }
    }
    for (int k = 1; k <= kmax2; ++k) {
        sxie2[I3(1,k,1)] = 0.0;
        for (int j = 2; j <= jmax2; ++j) {
            double dx = x2[I3(j,k,1)]-x2[I3(j-1,k,1)];
            double dy = y2[I3(j,k,1)]-y2[I3(j-1,k,1)];
            double dz = z2[I3(j,k,1)]-z2[I3(j-1,k,1)];
            sxie2[I3(j,k,1)] = sxie2[I3(j-1,k,1)] + std::sqrt(dx*dx+dy*dy+dz*dz);
        }
    }

    if (kcorr == 1) {
        // try arc length correction near eta=0
        for (int j = j21; j <= j22-1; ++j) {
            int km = kmm[j-1];
            if (km < k21) continue;   // go to 6901

            double test1 = 0.5;
            double test2 = 0.75;
            int kup  = (int)(test1*kmm[j-1]);
            int kup2 = (int)(test2*kmm[j-1]);
            if (kup2 <= k21) continue;   // go to 6901

            for (int k = k21; k <= kup2; ++k) {
                // center of "to" grid cell j,k
                double sc = 0.25*(seta2[I3(j,k,1)]+seta2[I3(j+1,k,1)]
                                 +seta2[I3(j+1,k+1,1)]+seta2[I3(j,k+1,1)]);
                int ll = (j22-j21)*(k-k21) + (j-j21+1);
                int jp   = (int)xie2[ll-1];
                double xie = xie2[ll-1]-jp;
                int lp   = mblkpt[ll-1];
                int jmax = jjmax1[lp-1]; (void)jmax;
                int kmax = kkmax1[lp-1]; (void)kmax;
                int kp;
                if (k == k21) {
                    kp = 1;
                } else {
                    kp = (int)eta2s[I2(j,k-1)]+1;
                }
                double eta = kp;
                int lsrch = 2;
                double xiet, etat;
                topol2(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,xie,eta,seta,
                       limit0,sc,jp,kp,lp,lsrch,itmax,xiet,etat,xif1,xif2,
                       etf1,etf2,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl);
                eta2s[I2(j,k)] = etat;
            }

            // blend corrected/uncorrected coordinates
            for (int k = k21; k <= kup; ++k) {
                int ll = (j22-j21)*(k-k21) + (j-j21+1);
                eta2[ll-1] = eta2s[I2(j,k)];
            }
            if (kup2 > kup) {
                for (int k = kup; k <= kup2; ++k) {
                    double phi1 = 1.0-f_float(k-kup)/f_float(kup2-kup);
                    double phi2 = 1.0-phi1;
                    int ll = (j22-j21)*(k-k21) + (j-j21+1);
                    eta2[ll-1] = eta2[ll-1]*phi2+eta2s[I2(j,k)]*phi1;
                }
            }

            // check for monotonicity of corrected coordinates
            for (int k = k21+1; k <= km; ++k) {
                int ll  = (j22-j21)*(k-k21)   + (j-j21+1);
                int ll1 = (j22-j21)*(k-1-k21) + (j-j21+1);
                if (eta2[ll-1] <= eta2[ll1-1]) {
                    std::snprintf(line,sizeof(line),
                        "  WARNING: corrected eta values are"
                        " not monotonic at j,k = %5d,%5d",j,k);
                    bou_put(nou,bou,ibufdim,4,line);
                }
            }
        }
    }

    if (jcorr == 1) {
        // try arc length correction near xie=0
        for (int k = k21; k <= k22-1; ++k) {
            int jm = jmm[k-1];
            if (jm < j21) continue;   // go to 7901

            double test1 = 0.5;
            double test2 = 0.75;
            int jup  = (int)(test1*jmm[k-1]);
            int jup2 = (int)(test2*jmm[k-1]);
            if (jup2 <= j21) continue;   // go to 7901

            for (int j = j21; j <= jup2; ++j) {
                double sc = 0.25*(sxie2[I3(j,k,1)]+sxie2[I3(j+1,k,1)]
                                 +sxie2[I3(j+1,k+1,1)]+sxie2[I3(j,k+1,1)]);
                int ll = (j22-j21)*(k-k21) + (j-j21+1);
                int kp   = (int)eta2[ll-1];
                double eta = eta2[ll-1]-kp;
                int lp   = mblkpt[ll-1];
                int jmax = jjmax1[lp-1]; (void)jmax;
                int kmax = kkmax1[lp-1]; (void)kmax;
                int jp;
                if (j == j21) {
                    jp = 1;
                } else {
                    jp = (int)xie2s[I2(j-1,k)]+1;
                }
                double xie = jp;
                int lsrch = 1;
                double xiet, etat;
                topol2(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,xie,eta,sxie,
                       limit0,sc,jp,kp,lp,lsrch,itmax,xiet,etat,xif1,xif2,
                       etf1,etf2,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl);
                xie2s[I2(j,k)] = xiet;
            }

            for (int j = j21; j <= jup; ++j) {
                int ll = (j22-j21)*(k-k21) + (j-j21+1);
                xie2[ll-1] = xie2s[I2(j,k)];
            }
            for (int j = jup; j <= jup2; ++j) {
                double phi1 = 1.0-f_float(j-jup)/f_float(jup2-jup);
                double phi2 = 1.0-phi1;
                int ll = (j22-j21)*(k-k21) + (j-j21+1);
                xie2[ll-1] = xie2[ll-1]*phi2+xie2s[I2(j,k)]*phi1;
            }

            // check for monotonicity of corrected coordinates
            for (int j = j21+1; j <= jm; ++j) {
                int ll  = (j22-j21)*(k-k21) + (j-j21+1);
                int ll1 = (j22-j21)*(k-k21) + (j-1-j21+1);
                if (xie2[ll-1] <= xie2[ll1-1]) {
                    std::snprintf(line,sizeof(line),
                        "  WARNING: corrected xie values are"
                        " not monotonic at j,k = %5d,%5d",j,k);
                    bou_put(nou,bou,ibufdim,4,line);
                }
            }
        }
    }
}
