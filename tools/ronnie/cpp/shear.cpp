// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// shear.cpp — faithful C++ translation of shear.F (CFL3D ronnie).
// Determine generalized coordinates of cell-edge midpoints on the xie=0 and
// eta=0 boundaries and apply the shearing (or arc-length) correction needed to
// render viscous "to"/"from" boundaries coincident.
//
// Sequential real build: real(x)->x, float(i)->double(i).
// Commons: /sklt1/ isklt1 -> sklt1_.isklt1
//          /areas/ ap(3),imaxa -> areas_.ap[0..2], areas_.imaxa
//          /tacos/ iretry -> tacos_.iretry
//          /filenam/ grdmov -> filenam_.grdmov
// Array layouts (column-major, 1-based):
//   x1/y1/z1, x1mid/..., x1mide/... : (jdim1,kdim1,msub1)
//   x2/y2/z2 : (jdim1,kdim1,msub2)
//   sxie2/seta2 : (jdim1,kdim1,msub2)
//   jimage/kimage : (msub1,jdim1,kdim1)
//   xie2/eta2/mblkpt : (npt);  temp : (jdim1*kdim1)
//   jjmax1/kkmax1/lout/xif1/xif2/etf1/etf2 : (msub1); jjmax2/kkmax2 : (msub2)
//   nblkj(jdim1),nblkk(kdim1),jmm(kdim1),kmm(jdim1); mblk2nd(maxbl)
//   bou(ibufdim,nbuf), nou(nbuf)
//
// When mmcxie=mmceta=0 (point-matched test) both correction blocks are skipped
// immediately (mceta==0 -> goto 670; mcxie==0 -> goto 770), but the full
// logic is reproduced faithfully.
#include "ron_common.h"

// forward decls (mirror the Fortran argument lists)
void extra(int jdim,int kdim,int nsub,int l,
           double* x,double* y,double* z,int jcell,int kcell,int jcl,int jcr,
           double& x5,double& y5,double& z5,int& icase,int ifit);
void extrae(int jdim,int kdim,int nsub,int l,
            double* x,double* y,double* z,int jcell,int kcell,int kcl,int kcr,
            double& x7,double& y7,double& z7,int& icase,int ifit);
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
           double* xie2,double* eta2,int* mblkpt,double* temp,
           int* jimage,int* kimage,int ifit,int itmax,double& xc,double& yc,
           double& zc,double* sxie2,double* seta2,int& jcorr,int& kcorr,
           int intmx,int icheck,int* nblkj,int* nblkk,int* jmm,int* kmm,
           int mcxie,int mceta,int* lout,int j21,int j22,int k21,int k22,
           int npt,int ic0,int iorph,int itoss0,int* xif1,int* xif2,
           int* etf1,int* etf2,int iself,int ifiner,
           int* nou,FStr* bou,int nbuf,int ibufdim,int myid,
           int* mblk2nd,int maxbl)
{
    (void)intmx; (void)npt; (void)mblk2nd; (void)maxbl;
    const int d1 = jdim1, d2 = kdim1;
    #define X1(i,j,k)    x1   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y1(i,j,k)    y1   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z1(i,j,k)    z1   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define X2(i,j,k)    x2   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y2(i,j,k)    y2   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z2(i,j,k)    z2   [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define SXIE2(i,j,k) sxie2[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define SETA2(i,j,k) seta2[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    // 1-d arrays (1-based)
    #define XIE2(i)   xie2[(i)-1]
    #define ETA2(i)   eta2[(i)-1]
    #define MBLKPT(i) mblkpt[(i)-1]
    #define TEMP(i)   temp[(i)-1]
    #define JJMAX1(i) jjmax1[(i)-1]
    #define KKMAX1(i) kkmax1[(i)-1]
    #define JJMAX2(i) jjmax2[(i)-1]
    #define KKMAX2(i) kkmax2[(i)-1]
    #define NBLKJ(i)  nblkj[(i)-1]
    #define NBLKK(i)  nblkk[(i)-1]
    #define JMM(i)    jmm[(i)-1]
    #define KMM(i)    kmm[(i)-1]

    int  l2    = 1;
    int  jmax2 = JJMAX2(l2);
    int  kmax2 = KKMAX2(l2);
    int  jl = 1;
    int  jr = jmax2-1;
    int  kl = 1;
    int  kr = kmax2-1;

    int j,k,ll,l1,jcall,kcall,icase;
    double x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;
    double a1,a2,a3;
    int imaxa;
    int ifits,iretry_local,icount,iok;
    double xiet,etat;
    double etat1=0.,xiet1=0.,xiec,etac,fact,fact1,fact2,dis0,dis1;
    int l11=0,ifit1=0;
    int ncoin,mcoin;
    int js,ks,km,jm,kmaxck,jmaxck;
    double beta,red;

    if (sklt1_.isklt1 > 0) {
        bou_put(nou,bou,ibufdim,4,"     beginning check of boundary values");
    }

    //*** correct boundary values near eta=0 ***
    kcorr = 0;
    ncoin = 0;

    if (mceta == 0) {
        if (sklt1_.isklt1 > 0)
            bou_put(nou,bou,ibufdim,4,"       eta=0 boundaries not rendered coincident");
        goto L670;
    }
    if (ishear >= 0) {
        if (sklt1_.isklt1 > 0)
            bou_put(nou,bou,ibufdim,4,"       eta=0 boundaries being rendered coincident via shearing method");
    } else {
        if (sklt1_.isklt1 > 0)
            bou_put(nou,bou,ibufdim,4,"       eta=0 boundaries being rendered coincident via arc length method");
    }

    // loop over all "to" cells on eta=0 boundary
    for (j = j21; j <= j22-1; ++j) {
        // compute edge midpoints of first layer of "to" grid cells
        jcall = j;
        kcall = k21;
        extra(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,kcall,jl,jr,x5,y5,z5,icase,ifit);
        if (j == j21) {
            extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,kcall,kl,kr,x7,y7,z7,icase,ifit);
        }
        kcall = k21+1;
        extra(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,kcall,jl,jr,x6,y6,z6,icase,ifit);
        jcall = j+1;
        kcall = k21;
        extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,kcall,kl,kr,x8,y8,z8,icase,ifit);

        // compute normalized directed areas/unit normals of "to" cell
        if (itoss0 == 0) {
            direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,a1,a2,a3,imaxa,nou,bou,nbuf,ibufdim);
            areas_.ap[0] = a1;
            areas_.ap[1] = a2;
            areas_.ap[2] = a3;
        }

        ifits  = ifit;
        tacos_.iretry = 0;
        icount = 1;
    L18085:
        // compute midpoint of "to" cell on eta=0 boundary, consistent with ifit
        if (ifit == 1 || ifit == 4) {
            xc = 0.5*(X2(j,k21,l2)+X2(j+1,k21,l2));
            yc = 0.5*(Y2(j,k21,l2)+Y2(j+1,k21,l2));
            zc = 0.5*(Z2(j,k21,l2)+Z2(j+1,k21,l2));
        } else {
            xc = x5;
            yc = y5;
            zc = z5;
        }

        ll   = (j-j21+1);
        xiet = XIE2(ll) + 1.;
        etat = ETA2(ll) + 1.;
        l1   = MBLKPT(ll);

        // search over "from" blocks for xie,eta of "to" cell midpoint on eta=0
        topol(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,l1,x1,y1,z1,
              x1mid,y1mid,z1mid,x1mide,y1mide,z1mide,limit0,xc,yc,zc,
              xiet,etat,jimage,kimage,ifit,itmax,igap,iok,lout,ic0,
              itoss0,j,1,iself,xif1,xif2,etf1,etf2,nou,bou,nbuf,ibufdim,myid);

        // search routine unsuccessful...try alternative polynomial fit
        if (tacos_.iretry == 0 && iok != 1) {
            if (icount < 4) {
                newfit(ifits,ifit,icount);
                icount = icount + 1;
                goto L18085;
            } else {
                // don't shear orphan points
                if (iorph > 0) {
                    if (sklt1_.isklt1 > 0)
                        bou_put(nou,bou,ibufdim,4," quitting boundary correction because points have been flagged as orphans");
                    return;
                }
                if (ifiner > 0) {
                    bou_put(nou,bou,ibufdim,4,"        all attempts to find generalized coordinates");
                    {
                        char b[160];
                        std::snprintf(b,sizeof(b),"        of cell center j,k =%s,%s have been unsuccessful...",
                                      f_iw(j,4).c_str(), f_iw(1,4).c_str());
                        bou_put(nou,bou,ibufdim,4,b);
                        std::snprintf(b,sizeof(b),"        will use  averages of finer level: interpolation number%s",
                                      f_iw(ifiner,4).c_str());
                        bou_put(nou,bou,ibufdim,4,b);
                    }
                    bou_put(nou,bou,ibufdim,4,"        for ALL points on this interface");
                } else {
                    int len1;
                    char titl[15];
                    if (icheck > 99)      { len1 = 13; std::snprintf(titl,sizeof(titl),"patch_p3d.%3d",icheck); }
                    else if (icheck > 9)  { len1 = 12; std::snprintf(titl,sizeof(titl),"patch_p3d.%2d",icheck); }
                    else                  { len1 = 11; std::snprintf(titl,sizeof(titl),"patch_p3d.%1d",icheck); }
                    {
                        std::string s = " program terminated in dynamic patching routines - see file "
                                        + std::string(filenam_.grdmov.c(),60);
                        bou_put(nou,bou,ibufdim,1,s);
                    }
                    bou_put(nou,bou,ibufdim,4,"         stopping...all attempts to find generalized coordinates");
                    {
                        char b[200];
                        std::snprintf(b,sizeof(b),"         of cell midpoint j,k = %s%s have been unsuccessful...check%.*s",
                                      f_iw(j,0).c_str(), f_iw(1,0).c_str(), len1, titl);
                        bou_put(nou,bou,ibufdim,4,b);
                    }
                }
                istop = 1;
                iout  = 1;
                return;
            }
        }

        // search routine successful with current polynomial fit
        if (iok == 1) {
            etat1 = etat;
            xiet1 = xiet;
            l11   = l1;
            ifit1 = ifit;
        }

        // try again with input polynomial fit
        if (ifit != ifits) {
            tacos_.iretry = 1;
            ifit = ifits;
            goto L18085;
        }

        // if second try with input value of ifit not successful, use the ifit
        // value for which the search routine was successful
        if (tacos_.iretry == 1) {
            if (iok != 1) {
                etat = etat1;
                xiet = xiet1;
                l1   = l11;
                if (sklt1_.isklt1 > 0) {
                    char b[160];
                    std::snprintf(b,sizeof(b),"         iterations using original fit  not successful at j,k= %s%s",
                                  f_iw(j,0).c_str(), f_iw(k21,0).c_str());
                    bou_put(nou,bou,ibufdim,4,b);
                    if (ifit == 1) bou_put(nou,bou,ibufdim,4,"          used bilinear fit instead");
                    if (ifit == 2) bou_put(nou,bou,ibufdim,4,"          used biquadratic fit instead");
                    if (ifit == 3) bou_put(nou,bou,ibufdim,4,"          used quadratic fit in xie, linear fit in eta instead");
                    if (ifit == 4) bou_put(nou,bou,ibufdim,4,"          used linear fit in xie, quadratic fit in eta instead");
                }
            }
        }

        // block locations for k=k21 boundary points stored in nblkj(j)
        NBLKJ(j) = l1;
        ifit     = ifits;
        etac     = etat - 2.0;
        fact     = 1.;

        mcoin = mceta*(j22-j21)/100;
        if (etac > 1.) {
            ncoin = ncoin+1;
            if (ncoin >= mcoin) {
                if (sklt1_.isklt1 > 0) {
                    bou_put(nou,bou,ibufdim,4,"      quitting check of boundary...");
                    bou_put(nou,bou,ibufdim,4,"      \"to\" and \"from\" eta=0 lines apparently not supposed to be coincident at this boundary");
                }
                goto L670;
            }
        }

        if (etat < 2.0) {
            // correction for widely varying cell sizes near the boundary
            js   = f_int(xiet);
            dis1 = (X1(js,3,l1)-X1(js,2,l1))*(X1(js,3,l1)-X1(js,2,l1))
                 + (Y1(js,3,l1)-Y1(js,2,l1))*(Y1(js,3,l1)-Y1(js,2,l1))
                 + (Z1(js,3,l1)-Z1(js,2,l1))*(Z1(js,3,l1)-Z1(js,2,l1));
            dis0 = (X1(js,2,l1)-X1(js,1,l1))*(X1(js,2,l1)-X1(js,1,l1))
                 + (Y1(js,2,l1)-Y1(js,1,l1))*(Y1(js,2,l1)-Y1(js,1,l1))
                 + (Z1(js,2,l1)-Z1(js,1,l1))*(Z1(js,2,l1)-Z1(js,1,l1));
            fact1 = std::sqrt(dis0/dis1);
            js   = js + 1;
            dis1 = (X1(js,3,l1)-X1(js,2,l1))*(X1(js,3,l1)-X1(js,2,l1))
                 + (Y1(js,3,l1)-Y1(js,2,l1))*(Y1(js,3,l1)-Y1(js,2,l1))
                 + (Z1(js,3,l1)-Z1(js,2,l1))*(Z1(js,3,l1)-Z1(js,2,l1));
            dis0 = (X1(js,2,l1)-X1(js,1,l1))*(X1(js,2,l1)-X1(js,1,l1))
                 + (Y1(js,2,l1)-Y1(js,1,l1))*(Y1(js,2,l1)-Y1(js,1,l1))
                 + (Z1(js,2,l1)-Z1(js,1,l1))*(Z1(js,2,l1)-Z1(js,1,l1));
            fact2 = std::sqrt(dis0/dis1);
            fact  = 0.5*( fact1 + fact2 );
            etac  = etac*fact;
        }

        TEMP(j) = etac;

        for (k = k21; k <= k22-1; ++k) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (MBLKPT(ll) == l1) {
                SETA2(j,k,1) = ETA2(ll);
                if (ETA2(ll) < 1.) SETA2(j,k,1) = (ETA2(ll)-1.)*fact+1.;
            }
        }

        // find km, last k-point on current j-line in same "from" block as k=1 pt
        KMM(j) = 0;
        for (k = k21; k <= k22-1; ++k) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (MBLKPT(ll) != NBLKJ(j)) break;
            KMM(j) = k;
        }
        x7 = x8;
        y7 = y8;
        z7 = z8;
    } // 2000

    // shearing correction (eta)
    if (ishear < 0) {
        kcorr = 1;
        goto L670;
    }

    beta = 25.;
    for (j = j21; j <= j22-1; ++j) {
        kmaxck = KKMAX1(NBLKJ(j));
        km     = KMM(j);
        if (km <= 1) continue;
        bool quit = false;
        for (k = k21; k <= km; ++k) {
            fact = 1. - ((double)(km-k)-0.5)/(double)(km-k21);
            red  = std::exp(-beta*fact*fact);
            SETA2(j,k,1) = SETA2(j,k,1) - TEMP(j)*red;

            if (SETA2(j,k,1) < 1. || SETA2(j,k,1) > (double)(kmaxck-2)) {
                if (sklt1_.isklt1 > 0) {
                    char b[160];
                    std::snprintf(b,sizeof(b),"        quitting eta shearing correction at j,k=%s%s",
                                  f_iw(j,0).c_str(), f_iw(k,0).c_str());
                    bou_put(nou,bou,ibufdim,4,b);
                    std::snprintf(b,sizeof(b),"        eta2(j,k),temp(j)= %s %s",
                                  f_glist(SETA2(j,k,1)).c_str(), f_glist(TEMP(j)).c_str());
                    bou_put(nou,bou,ibufdim,4,b);
                    bou_put(nou,bou,ibufdim,4,"        will try arc length correction near eta boundary ");
                }
                kcorr = 1;
                quit = true;
                break;
            }
        }
        if (quit) goto L670;

        for (k = k21; k <= k22-1; ++k) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (MBLKPT(ll) == NBLKJ(j)) ETA2(ll) = SETA2(j,k,1);
        }
    }
L670:

    //*** correct boundary values near xie=0 ***
    jcorr = 0;
    ncoin = 0;

    if (mcxie == 0) {
        if (sklt1_.isklt1 > 0)
            bou_put(nou,bou,ibufdim,4,"       xie=0 boundaries not rendered coincident");
        goto L770;
    }
    if (ishear >= 0) {
        if (sklt1_.isklt1 > 0)
            bou_put(nou,bou,ibufdim,4,"       xie=0 boundaries being rendered coincident via shearing method");
    } else {
        if (sklt1_.isklt1 > 0)
            bou_put(nou,bou,ibufdim,4,"       xie=0 boundaries being rendered coincident via arc length method");
    }

    for (k = k21; k <= k22-1; ++k) {
        kcall = k;
        jcall = j21;
        extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,kcall,kl,kr,x7,y7,z7,icase,ifit);
        if (k == k21) {
            extra(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,kcall,jl,jr,x5,y5,z5,icase,ifit);
        }
        jcall = j21+1;
        extrae(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,kcall,kl,kr,x8,y8,z8,icase,ifit);
        jcall = j21;
        kcall = k+1;
        extra(jdim1,kdim1,msub2,l2,x2,y2,z2,jcall,kcall,jl,jr,x6,y6,z6,icase,ifit);

        if (itoss0 == 0) {
            direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,a1,a2,a3,imaxa,nou,bou,nbuf,ibufdim);
            areas_.ap[0] = a1;
            areas_.ap[1] = a2;
            areas_.ap[2] = a3;
        }

        ifits  = ifit;
        tacos_.iretry = 0;
        icount = 1;
    L19085:
        if (ifit == 1 || ifit == 3) {
            xc = 0.5*(X2(j21,k+1,l2)+X2(j21,k,l2));
            yc = 0.5*(Y2(j21,k+1,l2)+Y2(j21,k,l2));
            zc = 0.5*(Z2(j21,k+1,l2)+Z2(j21,k,l2));
        } else {
            xc = x7;
            yc = y7;
            zc = z7;
        }

        ll   = (j22-j21)*(k-k21) + 1;
        xiet = XIE2(ll) + 1.;
        etat = ETA2(ll) + 1.;
        l1   = MBLKPT(ll);

        topol(jdim1,kdim1,msub1,jjmax1,kkmax1,lmax1,l1,x1,y1,z1,
              x1mid,y1mid,z1mid,x1mide,y1mide,z1mide,limit0,xc,yc,zc,
              xiet,etat,jimage,kimage,ifit,itmax,igap,iok,lout,ic0,
              itoss0,1,k,iself,xif1,xif2,etf1,etf2,nou,bou,nbuf,ibufdim,myid);

        if (tacos_.iretry == 0 && iok != 1) {
            if (icount < 4) {
                newfit(ifits,ifit,icount);
                icount = icount + 1;
                goto L19085;
            } else {
                if (iorph > 0) {
                    if (sklt1_.isklt1 > 0)
                        bou_put(nou,bou,ibufdim,4," quitting boundary correction because points have been flagged as orphans");
                    return;
                }
                if (ifiner > 0) {
                    bou_put(nou,bou,ibufdim,4,"        all attempts to find generalized coordinates");
                    {
                        char b[160];
                        std::snprintf(b,sizeof(b),"        of cell center j,k =%s,%s have been unsuccessful...",
                                      f_iw(1,4).c_str(), f_iw(k,4).c_str());
                        bou_put(nou,bou,ibufdim,4,b);
                        std::snprintf(b,sizeof(b),"        will use  averages of finer level: interpolation number%s",
                                      f_iw(ifiner,4).c_str());
                        bou_put(nou,bou,ibufdim,4,b);
                    }
                    bou_put(nou,bou,ibufdim,4,"        for ALL points on this interface");
                } else {
                    int len1;
                    char titl[15];
                    if (icheck > 99)      { len1 = 13; std::snprintf(titl,sizeof(titl),"patch_p3d.%3d",icheck); }
                    else if (icheck > 9)  { len1 = 12; std::snprintf(titl,sizeof(titl),"patch_p3d.%2d",icheck); }
                    else                  { len1 = 11; std::snprintf(titl,sizeof(titl),"patch_p3d.%1d",icheck); }
                    {
                        std::string s = " program terminated in dynamic patching routines - see file "
                                        + std::string(filenam_.grdmov.c(),60);
                        bou_put(nou,bou,ibufdim,1,s);
                    }
                    bou_put(nou,bou,ibufdim,4,"         stopping...all attempts to find generalized coordinates");
                    {
                        char b[200];
                        std::snprintf(b,sizeof(b),"         of cell center j,k =  1,k,have been unsuccessful...check%.*s",
                                      len1, titl);
                        bou_put(nou,bou,ibufdim,4,b);
                    }
                }
                istop = 1;
                iout  = 1;
                return;
            }
        }

        if (iok == 1) {
            xiet1 = xiet;
            etat1 = etat;
            l11   = l1;
            ifit1 = ifit;
        }

        if (ifit != ifits) {
            tacos_.iretry = 1;
            ifit = ifits;
            goto L19085;
        }

        if (tacos_.iretry == 1) {
            if (iok != 1) {
                xiet = xiet1;
                etat = etat1;
                l1   = l11;
                if (sklt1_.isklt1 > 0) {
                    char b[160];
                    std::snprintf(b,sizeof(b),"         iterations using original fit not successful at j,k= %s%s",
                                  f_iw(j21,0).c_str(), f_iw(k,0).c_str());
                    bou_put(nou,bou,ibufdim,4,b);
                    if (ifit == 1) bou_put(nou,bou,ibufdim,4,"          used bilinear fit instead");
                    if (ifit == 2) bou_put(nou,bou,ibufdim,4,"          used biquadratic fit instead");
                    if (ifit == 3) bou_put(nou,bou,ibufdim,4,"          used quadratic fit in xie, linear fit in eta instead");
                    if (ifit == 4) bou_put(nou,bou,ibufdim,4,"          used linear fit in xie, quadratic fit in eta instead");
                }
            }
        }

        // block locations for j=j21 boundary points stored in nblkk(k)
        NBLKK(k) = l1;
        ifit     = ifits;
        xiec     = xiet - 2.0;
        fact     = 1.;

        mcoin = mcxie*(k22-k21)/100;
        if (xiec > 1.) {
            ncoin = ncoin+1;
            if (ncoin >= mcoin) {
                if (sklt1_.isklt1 > 0) {
                    bou_put(nou,bou,ibufdim,4,"      quitting check of boundary...");
                    bou_put(nou,bou,ibufdim,4,"      \"to\" and \"from\" xie=0 lines apparently not supposed to be coincident at this boundary");
                }
                goto L770;
            }
        }

        if (xiet < 2.0) {
            ks   = f_int(etat);
            dis1 = (X1(3,ks,l1)-X1(2,ks,l1))*(X1(3,ks,l1)-X1(2,ks,l1))
                 + (Y1(3,ks,l1)-Y1(2,ks,l1))*(Y1(3,ks,l1)-Y1(2,ks,l1))
                 + (Z1(3,ks,l1)-Z1(2,ks,l1))*(Z1(3,ks,l1)-Z1(2,ks,l1));
            dis0 = (X1(2,ks,l1)-X1(1,ks,l1))*(X1(2,ks,l1)-X1(1,ks,l1))
                 + (Y1(2,ks,l1)-Y1(1,ks,l1))*(Y1(2,ks,l1)-Y1(1,ks,l1))
                 + (Z1(2,ks,l1)-Z1(1,ks,l1))*(Z1(2,ks,l1)-Z1(1,ks,l1));
            fact1 = std::sqrt(dis0/dis1);
            ks   = ks + 1;
            dis1 = (X1(3,ks,l1)-X1(2,ks,l1))*(X1(3,ks,l1)-X1(2,ks,l1))
                 + (Y1(3,ks,l1)-Y1(2,ks,l1))*(Y1(3,ks,l1)-Y1(2,ks,l1))
                 + (Z1(3,ks,l1)-Z1(2,ks,l1))*(Z1(3,ks,l1)-Z1(2,ks,l1));
            dis0 = (X1(2,ks,l1)-X1(1,ks,l1))*(X1(2,ks,l1)-X1(1,ks,l1))
                 + (Y1(2,ks,l1)-Y1(1,ks,l1))*(Y1(2,ks,l1)-Y1(1,ks,l1))
                 + (Z1(2,ks,l1)-Z1(1,ks,l1))*(Z1(2,ks,l1)-Z1(1,ks,l1));
            fact2 = std::sqrt(dis0/dis1);
            fact  = 0.5*( fact1 + fact2 );
            xiec  = xiec*fact;
        }

        TEMP(k) = xiec;

        for (j = j21; j <= j22-1; ++j) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (MBLKPT(ll) == l1) {
                SXIE2(j,k,1) = XIE2(ll);
                if (XIE2(ll) < (double)j21) SXIE2(j,k,1) = (XIE2(ll)-1.)*fact+1.;
            }
        }

        // find jm, last j-point on current k-line in same "from" block as j=1 pt
        JMM(k) = 0;
        for (j = j21; j <= j22-1; ++j) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (MBLKPT(ll) != NBLKK(k)) break;
            JMM(k) = j;
        }
        x5 = x6;
        y5 = y6;
        z5 = z6;
    } // 4000

    // shearing correction (xie)
    if (ishear < 0) {
        jcorr = 1;
        goto L770;
    }

    beta = 25.;
    for (k = k21; k <= k22-1; ++k) {
        jmaxck = JJMAX1(NBLKK(k));
        jm     = JMM(k);
        if (jm <= j21) continue;
        bool quit = false;
        for (j = j21; j <= jm; ++j) {
            fact = 1. - ((double)(jm-j)-0.5)/(double)(jm-j21);
            red  = std::exp(-beta*fact*fact);
            SXIE2(j,k,1) = SXIE2(j,k,1) - TEMP(k)*red;

            if (SXIE2(j,k,1) < 1. || SXIE2(j,k,1) > (double)(jmaxck-2)) {
                if (sklt1_.isklt1 > 0) {
                    char b[160];
                    std::snprintf(b,sizeof(b),"        quitting xie shearing correction at j,k=%s%s",
                                  f_iw(j,0).c_str(), f_iw(k,0).c_str());
                    bou_put(nou,bou,ibufdim,4,b);
                    std::snprintf(b,sizeof(b),"        xie2(j,k),temp(k)= %s %s",
                                  f_glist(SXIE2(j,k,1)).c_str(), f_glist(TEMP(k)).c_str());
                    bou_put(nou,bou,ibufdim,4,b);
                    bou_put(nou,bou,ibufdim,4,"        will try arc length correction near xie boundary");
                }
                jcorr = 1;
                quit = true;
                break;
            }
        }
        if (quit) goto L770;

        for (j = j21; j <= j22-1; ++j) {
            ll = (j22-j21)*(k-k21) + (j-j21+1);
            if (MBLKPT(ll) == NBLKK(k)) XIE2(ll) = SXIE2(j,k,1);
        }
    }
L770:
    return;

    #undef X1
    #undef Y1
    #undef Z1
    #undef X2
    #undef Y2
    #undef Z2
    #undef SXIE2
    #undef SETA2
    #undef XIE2
    #undef ETA2
    #undef MBLKPT
    #undef TEMP
    #undef JJMAX1
    #undef KKMAX1
    #undef JJMAX2
    #undef KKMAX2
    #undef NBLKJ
    #undef NBLKK
    #undef JMM
    #undef KMM
}
