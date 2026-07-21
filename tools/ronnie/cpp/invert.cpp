// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// invert.cpp — faithful C++ translation of invert.F (CFL3D ronnie).
// Determine generalized coordinates of the cell centers of the "to" grid in
// terms of the generalized coordinate system(s) defined on the "from" grid(s).
//
// COMMON mapping:
//   /sklt1/ isklt1   -> sklt1_.isklt1
//   /areas/ ap(3),imaxa -> areas_.ap[0..2], areas_.imaxa   (BY POSITION)
//   /tacos/ iretry   -> tacos_.iretry
//   /tracer/ itrace  -> tracer_.itrace
//   /filenam/ grdmov -> filenam_.grdmov
//   (invert.F declares NO /tol/.)
//
// Array dims (column-major, 1-based):
//   x1/y1/z1, x1mid.., x1mide.. : (jdim1,kdim1,msub1)
//   x2/y2/z2, sxie2/seta2       : (jdim1,kdim1,msub2)
//   sxie/seta                   : (jdim1,kdim1,msub1)
//   jimage/kimage               : (msub1,jdim1,kdim1)
//   xie2/eta2/mblkpt            : (npt)
//   xie2f/eta2f/mblkptf         : (nptf)
//   temp                        : (jdim1*kdim1)
//   xie2s/eta2s                 : (jdim1,kdim1)
//   nblkj(jdim1),nblkk(kdim1),jmm(kdim1),kmm(jdim1)
//   jjmax1/kkmax1/lout/xif1/xif2/etf1/etf2 : (msub1)
//   jjmax2/kkmax2 : (msub2)
//   bou(ibufdim,nbuf), nou(nbuf)
#include "ron_common.h"

// forward decls (add to ron_common.h if not present elsewhere)
void extra(int jdim,int kdim,int nsub,int l,double* x,double* y,double* z,
           int j,int k,int jl,int jr,double& x5,double& y5,double& z5,
           int& icase,int ifit);
void extrae(int jdim,int kdim,int nsub,int l,double* x,double* y,double* z,
            int j,int k,int kl,int kr,double& x7,double& y7,double& z7,
            int& icase,int ifit);
void direct(double x5,double x6,double x7,double x8,
            double y5,double y6,double y7,double y8,
            double z5,double z6,double z7,double z8,
            double& a1,double& a2,double& a3,int& imaxa,
            int* nou,FStr* bou,int nbuf,int ibufdim);
void newfit(int ifits,int& ifit,int icount);
void topol(int jdim,int kdim,int nsub,int* jjmax,int* kkmax,int lmax,int& l,
           double* x,double* y,double* z,
           double* xmid,double* ymid,double* zmid,
           double* xmide,double* ymide,double* zmide,int limit,
           double xc,double yc,double zc,double& xie,double& eta,
           int* jimage,int* kimage,int ifit,int itmax,int& igap,int& iok,
           int* lout,int ic0,int itoss0,int jto,int kto,int iself,
           int* xif1,int* xif2,int* etf1,int* etf2,
           int* nou,FStr* bou,int nbuf,int ibufdim,int myid);
void shear(int ishear,int& istop,int& iout,int igap,int jdim1,int kdim1,
           int msub1,int msub2,int* jjmax1,int* kkmax1,int lmax1,
           double* x1,double* y1,double* z1,double* x1mid,double* y1mid,
           double* z1mid,double* x1mide,double* y1mide,double* z1mide,
           int limit0,int* jjmax2,int* kkmax2,double* x2,double* y2,double* z2,
           double* xie2,double* eta2,int* mblkpt,double* temp,int* jimage,
           int* kimage,int ifit,int itmax,double& xc,double& yc,double& zc,
           double* sxie2,double* seta2,int& jcorr,int& kcorr,int intmx,
           int icheck,int* nblkj,int* nblkk,int* jmm,int* kmm,int mcxie,
           int mceta,int* lout,int j21,int j22,int k21,int k22,int npt,
           int ic0,int iorph,int itoss0,int* xif1,int* xif2,int* etf1,
           int* etf2,int iself,int ifiner,int* nou,FStr* bou,int nbuf,
           int ibufdim,int myid,int* mblk2nd,int maxbl);
void arc(int jdim1,int kdim1,int msub1,int msub2,int* jjmax1,int* kkmax1,
         int lmax1,double* x1,double* y1,double* z1,int limit0,int* jjmax2,
         int* kkmax2,double* x2,double* y2,double* z2,double* xie2,
         double* eta2,int* mblkpt,int ifit,int itmax,int jcorr,int kcorr,
         double* sxie,double* seta,double* sxie2,double* seta2,double* xie2s,
         double* eta2s,int intmx,int icheck,int* nblkj,int* nblkk,int* jmm,
         int* kmm,int j21,int j22,int k21,int k22,int npt,int* xif1,int* xif2,
         int* etf1,int* etf2,int* nou,FStr* bou,int nbuf,int ibufdim,
         int* mblk2nd,int maxbl);
void avgint(double* xie2,double* eta2,int* mblkpt,int npt,double* xie2f,
            double* eta2f,int* mblkptf,int nptf,int j21,int j22,int k21,
            int k22,int xi1f,int xi2f,int et1f,int et2f);
void diagnos(int istop,int iout,int igap,int jdim1,int kdim1,int msub1,
             int msub2,int* jjmax1,int* kkmax1,int lmax1,double* x1,double* y1,
             double* z1,double* x1mid,double* y1mid,double* z1mid,
             double* x1mide,double* y1mide,double* z1mide,double* sxie,
             double* seta,double* sxie2,double* seta2,double* xie2s,
             double* eta2s,int* jjmax2,int* kkmax2,double* x2,double* y2,
             double* z2,double* xie2,double* eta2,int* mblkpt,int icheck,
             int intmx,double xc,double yc,double zc,int ifit,int j21,int j22,
             int k21,int k22,int npt,int ic0,int iorph,int* xif1,int* xif2,
             int* etf1,int* etf2,int itoss0,int iself,int* nou,FStr* bou,
             int nbuf,int ibufdim,int myid);

void invert(int jdim1,int kdim1,int msub1,int msub2,int* jjmax1,int* kkmax1,
            int lmax1,double* x1,double* y1,double* z1,double* x1mid,
            double* y1mid,double* z1mid,double* x1mide,double* y1mide,
            double* z1mide,int limit0,int* jjmax2,int* kkmax2,double* x2,
            double* y2,double* z2,double* xie2,double* eta2,int* mblkpt,
            double* temp,int* jimage,int* kimage,int ifit,int itmax,
            double* sxie,double* seta,double* sxie2,double* seta2,
            double* xie2s,double* eta2s,int intmx,int icheck,int* nblkj,
            int* nblkk,int* jmm,int* kmm,int mcxie,int mceta,int* lout,
            int j21,int j22,int k21,int k22,int npt,int ic0,int iorph,
            int itoss0,int ncall,int ioutpt,int* xif1,int* xif2,int* etf1,
            int* etf2,int iself,int ifiner,double* xie2f,double* eta2f,
            int* mblkptf,int nptf,int xi1f,int xi2f,int et1f,int et2f,
            int iavg,int* nou,FStr* bou,int nbuf,int ibufdim,int myid,
            int* mblk2nd,int maxbl)
{
    (void)temp; (void)mcxie; (void)mceta;
    // storage-order accessors ---------------------------------------------
    const int d1 = jdim1, d2 = kdim1;
    #define X1(i,j,k)    x1   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y1(i,j,k)    y1   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z1(i,j,k)    z1   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define X2(i,j,k)    x2   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y2(i,j,k)    y2   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z2(i,j,k)    z2   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define X1MID(i,j,k)  x1mid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y1MID(i,j,k)  y1mid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z1MID(i,j,k)  z1mid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define X1MIDE(i,j,k) x1mide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y1MIDE(i,j,k) y1mide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z1MIDE(i,j,k) z1mide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    // xie2/eta2/mblkpt(npt): 1-based
    #define XIE2(i)   xie2  [(i)-1]
    #define ETA2(i)   eta2  [(i)-1]
    #define MBLKPT(i) mblkpt[(i)-1]

    int l2    = 1;
    int jmax2 = jjmax2[l2-1];
    int kmax2 = kkmax2[l2-1];
    int istop = 0;
    int iout  = 0;
    int igap  = 0;
    // xc,yc,zc: Fortran function-scope locals holding the last cell center;
    // read by shear/diagnos. Declared here (before any goto) so the value
    // survives 'goto L552'.
    double xc=0,yc=0,zc=0;
    if (ifit < 0) iout = 1;
    ifit = f_iabs(ifit);

    // itrace history control (kept for fidelity; here forced to -1)
    tracer_.itrace = -1;

    if (sklt1_.isklt1 > 0) {
        nou[4-1] = i_min(nou[4-1]+1,ibufdim);
        bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
            "     beginning inversion for xie, eta of cell centers");
        if (iavg == 0 || ifiner == 0) {
            if (ifit == 1) {
                nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                bou[(nou[4-1]-1)+(4-1)*ibufdim].assign("       using bi-linear fit");
            }
            if (ifit == 2) {
                nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                bou[(nou[4-1]-1)+(4-1)*ibufdim].assign("       using bi-quadratic fit");
            }
            if (ifit == 3) {
                nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                    "       using quadratic fit in xie, linear fit in eta");
            }
            if (ifit == 4) {
                nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                    "       using linear fit in xie, quadratic fit in eta");
            }
        } else if (ifiner > 0) {
            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
            char buf[128];
            std::snprintf(buf,sizeof(buf),
                "     using averages of finer level: interpolation number%4d",ifiner);
            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(buf);
        }
    }

    // compute cell edge midpoints of "from" grid(s) via quadratic least squares
    for (int l = 1; l <= lmax1; ++l) {
        int jmax1 = jjmax1[l-1];
        int kmax1 = kkmax1[l-1];
        for (int j = 1; j <= jmax1-1; ++j) {
            for (int k = 1; k <= kmax1; ++k) {
                int jl = 2;
                int jr = jmax1-2;
                double x5,y5,z5; int icase;
                extra(jdim1,kdim1,msub1,l,x1,y1,z1,j,k,jl,jr,x5,y5,z5,icase,ifit);
                X1MID(j,k,l) = x5;
                Y1MID(j,k,l) = y5;
                Z1MID(j,k,l) = z5;
            }
        }
        for (int k = 1; k <= kmax1-1; ++k) {
            for (int j = 1; j <= jmax1; ++j) {
                int kl = 2;
                int kr = kmax1-2;
                double x7,y7,z7; int icase;
                extrae(jdim1,kdim1,msub1,l,x1,y1,z1,j,k,kl,kr,x7,y7,z7,icase,ifit);
                X1MIDE(j,k,l) = x7;
                Y1MIDE(j,k,l) = y7;
                Z1MIDE(j,k,l) = z7;
            }
        }
    }

    // don't run inversion for a coarser level: coarse interpolation data are
    // always obtained by averaging fine level data.
    if (ifiner > 0 && iavg > 0) goto L552;

    {
    int nblsav = 1;
    int jl     = 1;
    int jr     = jmax2-1;
    int kl     = 1;
    int kr     = kmax2-1;

    // reused across cells: x5,y5,z5 (this cell) and x6,y6,z6 (prev in j-strip)
    double x5=0,y5=0,z5=0,x6=0,y6=0,z6=0,x7=0,y7=0,z7=0,x8=0,y8=0,z8=0;
    // (xc,yc,zc declared at function scope above.)
    int jcorr=0, kcorr=0;

    // loop over all "to" cells
    for (int j = j21; j <= j22-1; ++j) {
        for (int k = k22-1; k >= k21; --k) {

            // compute cell edge midpoints of "to" grid via quadratic LS
            if (k == k22-1) {
                int kcall = k+1;
                int icase;
                extra(jdim1,kdim1,msub2,l2,x2,y2,z2,j,kcall,jl,jr,x6,y6,z6,icase,ifit);
            }
            {
                int icase;
                extra(jdim1,kdim1,msub2,l2,x2,y2,z2,j,k,jl,jr,x5,y5,z5,icase,ifit);
                extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,j,k,kl,kr,x7,y7,z7,icase,ifit);
                int jcall = j+1;
                extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,k,kl,kr,x8,y8,z8,icase,ifit);
            }

            // compute normalized directed areas/unit normals of "to" cell
            if (itoss0 == 0) {
                double a1,a2,a3; int imaxa;
                direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,
                       a1,a2,a3,imaxa,nou,bou,nbuf,ibufdim);
                areas_.ap[1-1] = a1;
                areas_.ap[2-1] = a2;
                areas_.ap[3-1] = a3;
                areas_.imaxa   = imaxa;
            }

            int ifits  = ifit;
            tacos_.iretry = 0;
            int icount = 1;
            int iorfn  = 0;

            double xiet=0,etat=0;
            int l1=1;
            int iok=0;

        L17085: ;   // retry with alternate fit

            // compute center of "to" cell, consistent with ifit
            if (ifit == 1) {   // bi-linear
                xc = 0.25*( X2(j,k,l2) + X2(j+1,k,l2)
                          + X2(j+1,k+1,l2) + X2(j,k+1,l2) );
                yc = 0.25*( Y2(j,k,l2) + Y2(j+1,k,l2)
                          + Y2(j+1,k+1,l2) + Y2(j,k+1,l2) );
                zc = 0.25*( Z2(j,k,l2) + Z2(j+1,k,l2)
                          + Z2(j+1,k+1,l2) + Z2(j,k+1,l2) );
            }
            if (ifit == 2) {   // bi-quadratic
                xc = 0.5*( x5 + x6 + x7 + x8 )
                    -0.25*( X2(j,k,l2)     + X2(j+1,k,l2)
                          + X2(j+1,k+1,l2) + X2(j,k+1,l2) );
                yc = 0.5*( y5 + y6 + y7 + y8 )
                    -0.25*( Y2(j,k,l2)     + Y2(j+1,k,l2)
                          + Y2(j+1,k+1,l2) + Y2(j,k+1,l2) );
                zc = 0.5*( z5 + z6 + z7 + z8 )
                    -0.25*( Z2(j,k,l2)     + Z2(j+1,k,l2)
                          + Z2(j+1,k+1,l2) + Z2(j,k+1,l2) );
            }
            if (ifit == 3) {   // quadratic in xie, linear in eta
                xc = 0.5*(x5 + x6);
                yc = 0.5*(y5 + y6);
                zc = 0.5*(z5 + z6);
            }
            if (ifit == 4) {   // linear in xie, quadratic in eta
                xc = 0.5*(x7 + x8);
                yc = 0.5*(y7 + y8);
                zc = 0.5*(z7 + z8);
            }

            // set starting point for the search routine.
            //   first time..... use solution at previous point, or search for
            //                    minimum distance if xie=eta=0.
            //   subsequent..... use solution at same point from previous time.
            // xiet,etat: gen. coords in expanded "from" grid
            // xie2,eta2: gen. coords in original "from" grid
            if (tacos_.iretry == 0) {
                if (ncall == 1) {
                    if (k < k22-1) {
                        int ll = (j22-j21)*(k+1-k21) + (j-j21+1);
                        if (MBLKPT(ll) != 0) {
                            xiet = XIE2(ll) + 1.0;
                            etat = ETA2(ll) + 1.0;
                            l1   = MBLKPT(ll);
                        } else {
                            xiet = 0.0;
                            etat = 0.0;
                            l1   = nblsav;
                            l1   = i_max(l1,1);
                        }
                    } else if (j == j21) {
                        xiet = 0.0;
                        etat = 0.0;
                        l1   = 1;
                    } else if (j > j21) {
                        int ll = (j22-j21)*(k-k21) + (j-1-j21+1);
                        if (MBLKPT(ll) != 0) {
                            xiet = XIE2(ll) + 1.0;
                            etat = ETA2(ll) + 1.0;
                            l1   = MBLKPT(ll);
                        } else {
                            xiet = 0.0;
                            etat = 0.0;
                            l1   = nblsav;
                            l1   = i_max(l1,1);
                        }
                    }
                } else {
                    int ll = (j22-j21)*(k-k21) + (j-j21+1);
                    xiet = XIE2(ll) + 1.0;
                    etat = ETA2(ll) + 1.0;
                    l1   = MBLKPT(ll);
                    l1   = i_max(l1,1);
                }
            } else {
                int ll = (j22-j21)*(k-k21) + (j-j21+1);
                if (MBLKPT(ll) != 0) {
                    xiet = XIE2(ll) + 1.0;
                    etat = ETA2(ll) + 1.0;
                    l1   = MBLKPT(ll);
                } else {
                    xiet = 0.0;
                    etat = 0.0;
                    l1   = nblsav;
                    l1   = i_max(l1,1);
                }
            }

        L17086: ;

            // search over "from" blocks to find xie,eta for "to" center xc,yc,zc
            topol(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,l1,x1,y1,z1,
                  x1mid,y1mid,z1mid,x1mide,y1mide,z1mide,limit0,xc,yc,zc,
                  xiet,etat,jimage,kimage,ifit,itmax,igap,iok,lout,ic0,
                  itoss0,j,k,iself,xif1,xif2,etf1,etf2,nou,bou,nbuf,
                  ibufdim,myid);

            // search unsuccessful... try an alternative polynomial fit
            if (tacos_.iretry == 0 && iok != 1) {
                if (icount < 4) {
                    newfit(ifits,ifit,icount);
                    icount = icount + 1;
                    goto L17085;
                } else {
                    if (iorph <= 0) {
                        if (ifiner > 0) {
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                                "        all attempts to find generalized coordinates");
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            {
                                char b[128];
                                std::snprintf(b,sizeof(b),
                                    "        of cell center j,k =%4d,%4d have been unsuccessful...",
                                    j,k);
                                bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(b);
                            }
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            {
                                char b[128];
                                std::snprintf(b,sizeof(b),
                                    "        will use  averages of finer level: interpolation number%4d",
                                    ifiner);
                                bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(b);
                            }
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                                "        for ALL points on this interface");
                            istop = 1;
                            goto L552;
                        } else {
                            // build titlptchgrd (character*14)
                            char titl[15];
                            int len1;
                            if (icheck > 99) {
                                len1 = 13;
                                std::snprintf(titl,sizeof(titl),"patch_p3d.%3d",icheck);
                            } else if (icheck > 9) {
                                len1 = 12;
                                std::snprintf(titl,sizeof(titl),"patch_p3d.%2d",icheck);
                            } else {
                                len1 = 11;
                                std::snprintf(titl,sizeof(titl),"patch_p3d.%1d",icheck);
                            }
                            for (int i = len1+1; i <= 14; ++i) titl[i-1] = ' ';
                            titl[14] = '\0';

                            nou[1-1] = i_min(nou[1-1]+1,ibufdim);
                            {
                                std::string s = "  program terminated in dynamic patching routines - see file "
                                                + std::string(filenam_.grdmov.c(),60);
                                bou[(nou[1-1]-1)+(1-1)*ibufdim].assign(s);
                            }
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                                "        stopping...all attempts to find generalized coordinates of");
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            {
                                char b[160];
                                std::snprintf(b,sizeof(b),
                                    "        cell center j,k =%4d,%4d have been unsuccessful...check%.*s",
                                    j,k,len1,titl);
                                bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(b);
                            }
                            istop = 1;
                            iout  = 1;
                            goto L552;
                        }
                    } else {
                        // mark orphan point as interpolated from block "0"
                        iok = 1;
                        int ll = (j22-j21)*(k-k21) + (j-j21+1);
                        nblsav      = MBLKPT(ll);
                        XIE2(ll)    = 1.0;
                        ETA2(ll)    = 1.0;
                        MBLKPT(ll)  = 0;
                        x6 = x5;
                        y6 = y5;
                        z6 = z5;
                        goto L1001;
                    }
                }
            }

            // search successful with current polynomial fit
            if (iok == 1) {
                int ll = (j22-j21)*(k-k21) + (j-j21+1);
                XIE2(ll)   = xiet-1.0;
                ETA2(ll)   = etat-1.0;
                MBLKPT(ll) = l1;
                // ifit1 = ifit;  (dead: value unused downstream)
            }

            if (iorph > 0) {
                // search "successful", but only because it found what might
                // otherwise be an orphan point in the expanded "from" grid.
                int ll = (j22-j21)*(k-k21) + (j-j21+1);
                if (XIE2(ll) < 1.0 || XIE2(ll) > jjmax1[l1-1]-2 ||
                    ETA2(ll) < 1.0 || ETA2(ll) > kkmax1[l1-1]-2) {
                    iorfn = iorfn + 1;
                    if (iorfn <= 1) {
                        // redo the search starting with a minimum distance
                        // search to confirm this is really an orphan
                        iok    = 0;
                        icount = 1;
                        tacos_.iretry = 0;
                        xiet   = 0.0;
                        etat   = 0.0;
                        goto L17086;
                    } else {
                        // mark the point as an orphan
                        iok = 1;
                        nblsav     = MBLKPT(ll);
                        XIE2(ll)   = 1.0;
                        ETA2(ll)   = 1.0;
                        MBLKPT(ll) = 0;
                        x6 = x5;
                        y6 = y5;
                        z6 = z5;
                        goto L1001;
                    }
                }
            }

            // if current polynomial fit is not the input one, retry with input
            if (ifit != ifits) {
                tacos_.iretry = 1;
                ifit = ifits;
                goto L17085;
            }

            // if second try with input ifit not successful, keep the ifit
            // for which the search succeeded.
            if (tacos_.iretry == 1) {
                if (iok != 1) {
                    if (sklt1_.isklt1 > 0) {
                        nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                        {
                            char b[160];
                            std::snprintf(b,sizeof(b),
                                "         iterations using original fit not successful at j,k=  %11d %11d",
                                j,k);
                            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(b);
                        }
                        if (ifit == 1) {
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                                "          used bilinear fit instead");
                        }
                        if (ifit == 2) {
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                                "          used biquadratic fit instead");
                        }
                        if (ifit == 3) {
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                                "          used quadratic fit in xie, linear fit in eta instead");
                        }
                        if (ifit == 4) {
                            nou[4-1] = i_min(nou[4-1]+1,ibufdim);
                            bou[(nou[4-1]-1)+(4-1)*ibufdim].assign(
                                "          used linear fit in xie, quadratic fit in eta instead");
                        }
                    }
                }
            }

            x6 = x5;
            y6 = y5;
            z6 = z5;

        L1001: ;
        } // k loop (1001)
    } // j loop (1000)

    // check and, if necessary, correct values of xie and/or eta near
    // boundaries.
    //   ishear > 0: try a shearing correction to make "to"/"from" boundaries
    //               coincident; if shearing fails, try an arc length correction.
    //   ishear < 0: use only the arc length correction.
    {
        int ishear = 1;

        shear(ishear,istop,iout,igap,jdim1,kdim1,msub1,msub2,jjmax1,kkmax1,
              lmax1,x1,y1,z1,x1mid,y1mid,z1mid,x1mide,y1mide,z1mide,limit0,
              jjmax2,kkmax2,x2,y2,z2,xie2,eta2,mblkpt,temp,jimage,kimage,ifit,
              itmax,xc,yc,zc,sxie2,seta2,jcorr,kcorr,
              intmx,icheck,nblkj,nblkk,jmm,kmm,mcxie,mceta,lout,j21,j22,k21,
              k22,npt,ic0,iorph,itoss0,xif1,xif2,etf1,etf2,iself,ifiner,
              nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl);

        if (jcorr != 0 || kcorr != 0)
            arc(jdim1,kdim1,msub1,msub2,jjmax1,kkmax1,lmax1,x1,y1,z1,limit0,
                jjmax2,kkmax2,x2,y2,z2,xie2,eta2,mblkpt,ifit,itmax,jcorr,kcorr,
                sxie,seta,sxie2,seta2,xie2s,eta2s,intmx,icheck,nblkj,nblkk,
                jmm,kmm,j21,j22,k21,k22,npt,xif1,xif2,etf1,etf2,
                nou,bou,nbuf,ibufdim,mblk2nd,maxbl);
    }
    } // end block guarded by (ifiner>0 && iavg>0) skip

L552:

    // for coarse level interfaces where the standard search algorithm has
    // failed, use average of finer level data.
    if (ifiner > 0 && istop == 1) {
        avgint(xie2,eta2,mblkpt,npt,xie2f,eta2f,mblkptf,nptf,
               j21,j22,k21,k22,xi1f,xi2f,et1f,et2f);
        istop = 0;
    }

    // perform diagnostic checks on xie,eta found for "to" cell centers
    if (ioutpt > 0 || istop > 0) {
        // xc,yc,zc: last-computed cell center (Fortran function-scope locals).
        diagnos(istop,iout,igap,jdim1,kdim1,msub1,msub2,jjmax1,kkmax1,lmax1,
                x1,y1,z1,x1mid,y1mid,z1mid,x1mide,y1mide,z1mide,sxie,seta,
                sxie2,seta2,xie2s,eta2s,jjmax2,kkmax2,x2,y2,z2,xie2,eta2,
                mblkpt,icheck,intmx,xc,yc,zc,ifit,j21,j22,k21,k22,npt,
                ic0,iorph,xif1,xif2,etf1,etf2,itoss0,iself,
                nou,bou,nbuf,ibufdim,myid);
    }
    return;

    #undef X1
    #undef Y1
    #undef Z1
    #undef X2
    #undef Y2
    #undef Z2
    #undef XIE2
    #undef ETA2
    #undef MBLKPT
}
