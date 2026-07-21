// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// topol.cpp — faithful C++ translation of topol.F (CFL3D ronnie).
// Search appropriate "from" blocks for current "to" cell center (xc,yc,zc);
// determine xie,eta of that cell center via Newton inversion (xe/xe2).
//
// /tol/ epsc,epsc0,epsreen,epscoll -> tol_.epsc (slot0), tol_.epsc0 (slot1).
// Array dims (column-major, 1-based):
//   x/y/z, xmid/ymid/zmid, xmide/ymide/zmide : (jdim,kdim,nsub)
//       (i,j,k) -> [(i-1)+(j-1)*jdim+(k-1)*jdim*kdim]
//   jimage/kimage : (nsub,jdim,kdim)
//       (i,j,k) -> [(i-1)+(j-1)*nsub+(k-1)*nsub*jdim]
//   jjmax/kkmax/lout/xif1/xif2/etf1/etf2 : (nsub)
//   jfroz/kfroz : local (itmax)
//   bou(ibufdim,nbuf), nou(nbuf)
//
// NOTE on numeric fidelity:
//   jp=xie, kp=eta assign REAL->INT (Fortran truncates toward zero) -> f_int.
//   jinc=abs(xie) etc. compute real abs then truncate to int -> f_int(fabs()).
#include "ron_common.h"

// Group D / B forward decls (declared here if not yet in ron_common.h):
void dsmin(int jdim,int kdim,int nsub,int* jjmax,int* kkmax,int lmax,
           double* x,double* y,double* z,double xc,double yc,double zc,
           int& jp,int& kp,int& l,int* lout,int iflag,
           int* xif1,int* xif2,int* etf1,int* etf2);
void xe(int jdim,int kdim,int nsub,int l,double* x,double* y,double* z,
        double* xmid,double* ymid,double* zmid,double* xmide,double* ymide,
        double* zmide,int jp,int kp,double xc,double yc,double zc,
        double& xie,double& eta,int& imiss,int ifit,int ic0,
        int* nou,FStr* bou,int nbuf,int ibufdim,int myid);
void xe2(double y1,double y2,double y3,double y4,double y5,double y6,
         double y7,double y8,double yc,double z1,double z2,double z3,
         double z4,double z5,double z6,double z7,double z8,double zc,
         double& xie,double& eta,int& imiss,int ifit,
         int* nou,FStr* bou,int nbuf,int ibufdim,int myid);

void topol(int jdim,int kdim,int nsub,int* jjmax,int* kkmax,int lmax,int& l,
           double* x,double* y,double* z,
           double* xmid,double* ymid,double* zmid,
           double* xmide,double* ymide,double* zmide,int limit,
           double xc,double yc,double zc,double& xie,double& eta,
           int* jimage,int* kimage,int ifit,int itmax,int& igap,int& iok,
           int* lout,int ic0,int itoss0,int jto,int kto,int iself,
           int* xif1,int* xif2,int* etf1,int* etf2,
           int* nou,FStr* bou,int nbuf,int ibufdim,int myid)
{
    (void)nou; (void)bou; (void)nbuf; (void)ibufdim;
    const int dg1 = jdim, dg2 = kdim;      // x/y/z leading dims
    #define X(i,j,k)    x   [((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    #define Y(i,j,k)    y   [((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    #define Z(i,j,k)    z   [((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    #define XMID(i,j,k) xmid[((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    #define YMID(i,j,k) ymid[((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    #define ZMID(i,j,k) zmid[((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    #define XMIDE(i,j,k) xmide[((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    #define YMIDE(i,j,k) ymide[((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    #define ZMIDE(i,j,k) zmide[((i)-1)+((j)-1)*dg1+((k)-1)*dg1*dg2]
    // jimage/kimage (nsub,jdim,kdim)
    #define JIMAGE(i,j,k) jimage[((i)-1)+((j)-1)*nsub+((k)-1)*nsub*dg1]
    #define KIMAGE(i,j,k) kimage[((i)-1)+((j)-1)*nsub+((k)-1)*nsub*dg1]

    std::vector<int> jfroz((size_t)itmax+1), kfroz((size_t)itmax+1); // 1-based use

    int iatmpt,ichk,ifroze,ihuge;
    int jp,kp,lsav,jpsav,kpsav,jpc,kpc;
    double xiesav,etasav;
    int jmax,kmax,js,je,ks,ke;
    int ifroz;      // frozen flag (declared, mirrors Fortran ifroz)
    int ibrdir = 0;

    for (int ll = 1; ll <= lmax; ++ll) lout[ll-1] = 0;
    iatmpt = 0;
    ichk   = 0;
    ifroze = 0;
    ifroz  = 0;
    ihuge  = 0;
    jp     = f_int(xie);       // jp = xie  (real->int truncation)
    kp     = f_int(eta);       // kp = eta
    lsav   = 1;
    xiesav = xie;
    etasav = eta;
    jpsav  = jp;
    kpsav  = kp;
    jpc    = jp;
    kpc    = kp;

    // to start search, use solution from last cell if one exists; otherwise
    // start by searching for minimum distance point over all "from" blocks
    if (xie < 1.0 || eta < 1.0) {
        dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
              jp,kp,l,lout,999,xif1,xif2,etf1,etf2);
    }

L999:
    iatmpt = iatmpt + 1;

    // all "from" blocks have been searched if iatmpt > lmax
    if (iatmpt > lmax) goto L1000;

    jmax = jjmax[l-1];
    kmax = kkmax[l-1];
    js   = xif1[l-1];
    je   = xif2[l-1];
    ks   = etf1[l-1];
    ke   = etf2[l-1];

    // avoid starting outside specified search range
    jp = i_min(jp,je-1);
    kp = i_min(kp,ke-1);
    jp = i_max(jp,js);
    kp = i_max(kp,ks);

    for (int intern = 1; intern <= itmax; ++intern) {
        jfroz[intern] = jp;
        kfroz[intern] = kp;

        int imiss = 0;

        // find local xie, eta via Newton iteration in current target cell jp,kp
        if (itoss0 == 0) {
            xe(jdim,kdim,nsub,l,x,y,z,xmid,ymid,zmid,xmide,ymide,zmide,
               jp,kp,xc,yc,zc,xie,eta,imiss,ifit,ic0,
               nou,bou,nbuf,ibufdim,myid);
        } else {
            double x1 = X(jp,kp,l),     y1 = Y(jp,kp,l),     z1 = Z(jp,kp,l);
            double x2 = X(jp+1,kp,l),   y2 = Y(jp+1,kp,l),   z2 = Z(jp+1,kp,l);
            double x4 = X(jp,kp+1,l),   y4 = Y(jp,kp+1,l),   z4 = Z(jp,kp+1,l);
            double x3 = X(jp+1,kp+1,l), y3 = Y(jp+1,kp+1,l), z3 = Z(jp+1,kp+1,l);
            double x5 = XMID(jp,kp,l),  y5 = YMID(jp,kp,l),  z5 = ZMID(jp,kp,l);
            double x6 = XMID(jp,kp+1,l),y6 = YMID(jp,kp+1,l),z6 = ZMID(jp,kp+1,l);
            double x7 = XMIDE(jp,kp,l), y7 = YMIDE(jp,kp,l), z7 = ZMIDE(jp,kp,l);
            double x8 = XMIDE(jp+1,kp,l),y8 = YMIDE(jp+1,kp,l),z8 = ZMIDE(jp+1,kp,l);

            if (itoss0 == 1) {
                // use only y and z equations
                xe2(y1,y2,y3,y4,y5,y6,y7,y8,yc,z1,z2,z3,z4,z5,z6,z7,z8,zc,
                    xie,eta,imiss,ifit,nou,bou,nbuf,ibufdim,myid);
            }
            if (itoss0 == 2) {
                // use only x and z equations
                xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,z1,z2,z3,z4,z5,z6,z7,z8,zc,
                    xie,eta,imiss,ifit,nou,bou,nbuf,ibufdim,myid);
            }
            if (itoss0 == 3) {
                // use only x and y equations
                xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,y1,y2,y3,y4,y5,y6,y7,y8,yc,
                    xie,eta,imiss,ifit,nou,bou,nbuf,ibufdim,myid);
            }
        }

        // check the search did not find the current "to" cell as target
        // (possible for a grid that communicates with itself along a branch
        // cut). If so, move to a cell on the other side of the branch cut.
        if (iself > 0) {
            if (imiss == 0 && (jto == jp-1 && kto == kp-1)) {
                // determine in which direction the branch cut lies
                jpc = jpc - 1;
                jpc = i_max(1,jpc);
                if (jpc != JIMAGE(l,jpc,kp)) ibrdir = 1;
                jpc = jpc + 1;
                jpc = i_min(jpc,jmax-1);
                if (jpc != JIMAGE(l,jpc,kp)) ibrdir = 2;
                kpc = kpc - 1;
                kpc = i_max(1,kpc);
                if (kpc != KIMAGE(l,jp,kpc)) ibrdir = 3;
                kpc = kpc + 1;
                kpc = i_min(kpc,kmax-1);
                if (kpc != KIMAGE(l,jp,kpc)) ibrdir = 4;

                // set cell index to trigger a branch cut jump
                if (ibrdir == 1) jp = jp - 1;
                if (ibrdir == 2) jp = jp + 1;
                if (ibrdir == 3) kp = kp - 1;
                if (ibrdir == 4) kp = kp + 1;

                jp = i_min(jp,jmax-1);
                kp = i_min(kp,kmax-1);
                jp = i_max(1,jp);
                kp = i_max(1,kp);
                jp = JIMAGE(l,jp,kp);
                kp = KIMAGE(l,jp,kp);
                imiss = 1;
                continue;   // go to 5555
            }
        }

        // current target cell correct if imiss = 0
        if (imiss == 0) goto L5556;

        // update current guess for target cell based on Newton result,
        // with max allowable change set by limit
        int jinc,kinc;
        if (xie >= 0.0) jinc = f_int(std::fabs(xie));
        else            jinc = f_int(std::fabs(xie-1.0));
        if (eta >= 0.0) kinc = f_int(std::fabs(eta));
        else            kinc = f_int(std::fabs(eta-1.0));

        jinc = i_min(jinc,limit);
        kinc = i_min(kinc,limit);

        if (xie > 1.0)      jp = jp + jinc;
        else if (xie < 0.0) jp = jp - jinc;
        if (eta > 1.0)      kp = kp + kinc;
        else if (eta < 0.0) kp = kp - kinc;

        // keep within bounds of (expanded) "from" block
        jp = i_min(jp,jmax-1);
        kp = i_min(kp,kmax-1);
        jp = i_max(1,jp);
        kp = i_max(1,kp);

        jmax = jjmax[l-1];
        kmax = kkmax[l-1];
        js   = xif1[l-1];
        je   = xif2[l-1];
        ks   = etf1[l-1];
        ke   = etf2[l-1];

        // avoid cells outside specified search range
        jp = i_min(jp,je-1);
        kp = i_min(kp,ke-1);
        jp = i_max(jp,js);
        kp = i_max(kp,ks);

        // account for any branch cuts
        jpc = JIMAGE(l,jp,kp);
        kpc = KIMAGE(l,jp,kp);
        jp = jpc;
        kp = kpc;

        // search routine off track if local xie or eta become huge...
        // try to get back on track via minimum distance search in current block
        {
            double hugev = 1.e+5;
            if (std::fabs(xie) > hugev || std::fabs(eta) > hugev) {
                ihuge = ihuge + 1;
                if (ihuge > 1) goto L1000;
                dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
                      jp,kp,l,lout,-999,xif1,xif2,etf1,etf2);
                continue;   // go to 5555
            }
        }

        // check for frozen convergence: search keeps returning to the same
        // point without 0 < xie,eta < 1. If frozen, break out via min-dist.
        ifroz = 0;
        for (int ii = 1; ii <= intern; ++ii) {
            int intv = intern - ii + 1;
            if (jp == jfroz[intv] && kp == kfroz[intv]) ifroz = 1;
        }
        if (ifroz == 1) {
            ifroze = ifroze + 1;
            if (ifroze > 1) goto L1000;
            dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
                  jp,kp,l,lout,-999,xif1,xif2,etf1,etf2);
        }
    } // 5555 continue (loop increment)

L1000:
    // search unsuccessful in the current "from" block; if all "from" blocks
    // have not been searched, try another; else use any previously obtained
    // generalized coordinates found in an expanded cell (ichk>0); otherwise
    // exit to try a different basis function and start over.
    if (iatmpt < lmax) {
        // find new starting cell by min-distance over all not-yet-searched blocks
        lout[l-1] = 1;
        dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
              jp,kp,l,lout,999,xif1,xif2,etf1,etf2);
        ifroze = 0;
        ihuge  = 0;
        goto L999;
    } else {
        if (ichk == 0) {
            iok = 0;
            return;
        } else {
            l   = lsav;
            jp  = jpsav;
            kp  = kpsav;
            xie = xiesav;
            eta = etasav;
        }
        if (ichk == 1) {
            // point really does lie in the expanded cell of the "from" block
            // in which it was first found
            goto L5557;
        }
        if (ichk > 1) {
            // point lies in the expanded zone of two or more grids, possibly
            // due to a small gap between blocks on the "from" side.
            igap = 1;
            goto L5557;
        }
    }

L5556:
    // search successful; however, with multiple "from" blocks, check whether:
    // 1) for C-0 interfaces, the "from" cell center matches the "to" cell
    //    center to insure the right "from" block was used;
    // 2) the "to" cell center lies in the overlap region formed when the
    //    "from" blocks are expanded (may have been found in the wrong block;
    //    alternatively, a gap may exist between blocks).
    if (lmax > 1) {
        if (ic0 > 0) {
            int j = jp;
            int k = kp;
            double x1c = 0.25*( X(j,k,l) + X(j+1,k,l)
                              + X(j+1,k+1,l) + X(j,k+1,l) );
            double y1c = 0.25*( Y(j,k,l) + Y(j+1,k,l)
                              + Y(j+1,k+1,l) + Y(j,k+1,l) );
            double z1c = 0.25*( Z(j,k,l) + Z(j+1,k,l)
                              + Z(j+1,k+1,l) + Z(j,k+1,l) );
            if (std::fabs(x1c-xc) > tol_.epsc0 ||
                std::fabs(y1c-yc) > tol_.epsc0 ||
                std::fabs(z1c-zc) > tol_.epsc0) {
                ichk = ichk+1;
                if (ichk == 1) {
                    lsav   = l;
                    jpsav  = jp;
                    kpsav  = kp;
                    xiesav = xie;
                    etasav = eta;
                }
                // find alternate block via min-distance over not-yet-searched
                lout[l-1] = 1;
                dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
                      jp,kp,l,lout,999,xif1,xif2,etf1,etf2);
                ifroze = 0;
                ihuge  = 0;
                goto L999;
            }
        }

        if (jp == 1 || jp == jmax-1 || kp == 1 || kp == kmax-1) {
            ichk = ichk+1;
            if (ichk == 1) {
                lsav   = l;
                jpsav  = jp;
                kpsav  = kp;
                xiesav = xie;
                etasav = eta;
            }
            // find alternate block via min-distance over not-yet-searched
            lout[l-1] = 1;
            dsmin(jdim,kdim,nsub,jjmax,kkmax,lmax,x,y,z,xc,yc,zc,
                  jp,kp,l,lout,999,xif1,xif2,etf1,etf2);
            ifroze = 0;
            ihuge  = 0;
            goto L999;
        }
    }

L5557:
    // completed search for "to" cell center; convert to global xie, eta
    iok = 1;
    xie = xie + jp;
    eta = eta + kp;
    return;

    #undef X
    #undef Y
    #undef Z
    #undef XMID
    #undef YMID
    #undef ZMID
    #undef XMIDE
    #undef YMIDE
    #undef ZMIDE
    #undef JIMAGE
    #undef KIMAGE
}
