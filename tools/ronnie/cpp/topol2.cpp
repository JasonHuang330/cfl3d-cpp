// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// topol2.cpp — faithful C++ translation of topol2.F (CFL3D ronnie).
// One-dimensional equivalent of topol: search in one direction only.
//
// /tol/ epsc,epsc0,epsreen,epscoll  -> tol_.epsc (slot0 = epsc; only epsc used).
// /filenam/ grdmov -> filenam_.grdmov.
// Array dims (column-major, 1-based):
//   sarc(jdim1,kdim1,msub1)  (i,j,k) -> [(i-1)+(j-1)*jdim1+(k-1)*jdim1*kdim1]
//   jjmax1(msub1),kkmax1(msub1),xif1/xif2/etf1/etf2(msub1),mblk2nd(maxbl)
//   bou(ibufdim,nbuf), nou(nbuf)
#include "ron_common.h"

// NOTE: signature matches the canonical declaration in ron_common.h:
//   xie,eta passed BY VALUE (callers do not read them back; the real outputs
//   are xiet,etat), and the sarc array is named `s` there. Fortran passed
//   xie/eta as arguments, but only the entry values are needed here.
void topol2(int jdim1,int kdim1,int msub1,int* jjmax1,int* kkmax1,int lmax1,
            double xie,double eta,double* s,int limit0,double sc,
            int jp,int kp,int lp,int lsrch,int itmax,double& xiet,double& etat,
            int* xif1,int* xif2,int* etf1,int* etf2,
            int* nou,FStr* bou,int nbuf,int ibufdim,int myid,
            int* mblk2nd,int maxbl)
{
    (void)jdim1; (void)kdim1; (void)msub1; (void)jjmax1; (void)kkmax1;
    (void)lmax1; (void)mblk2nd; (void)maxbl;
    double* sarc = s;
    const int d1 = jdim1, d2 = kdim1;   // sarc leading dims
    #define SARC(i,j,k) sarc[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]

    int jmax = jjmax1[lp-1];
    int kmax = kkmax1[lp-1];
    (void)jmax; (void)kmax;
    int js = xif1[lp-1];
    int je = xif2[lp-1];
    int ks = etf1[lp-1];
    int ke = etf2[lp-1];
    int intern = 0;
    int limit  = limit0;

L6901:
    intern = intern + 1;
    if (intern > itmax) {
        nou[1-1] = i_min(nou[1-1]+1,ibufdim);
        {
            std::string s = "  program terminated in dynamic patching routines - see file "
                            + std::string(filenam_.grdmov.c(),60);
            bou[(nou[1-1]-1)+(1-1)*ibufdim].assign(s);
        }
        nou[1-1] = i_min(nou[1-1]+1,ibufdim);
        bou[(nou[1-1]-1)+(1-1)*ibufdim].assign(" stopping...intern > itmax ");
        termn8(myid,-1,ibufdim,nbuf,bou,nou);
    }

    // bi-linear fit in target cell jp,kp in "from" block lp
    double dy2 = SARC(jp+1,kp,lp)   - SARC(jp,kp,lp);
    double dy3 = SARC(jp+1,kp+1,lp) - SARC(jp,kp,lp);
    double dy4 = SARC(jp,kp+1,lp)   - SARC(jp,kp,lp);
    double a1  = SARC(jp,kp,lp);
    double a2  = dy2;
    double a3  = dy4;
    double a4  = dy3-a2-a3;

    if (lsrch == 1) {
        // for fixed eta
        xie = (sc-a1-a3*eta)/(a2+a4*eta);

        // check to ensure point is inside cell
        int imiss = 0;
        if (xie < -tol_.epsc || xie > 1.0+tol_.epsc) {
            imiss = 1;
        }
        if (imiss == 0) goto L6902;

        // try new target cell
        int jpc = jp;
        int jinc = 0;
        if (xie >= 0.0) jinc = f_int(std::fabs(xie));
        if (xie <  0.0) jinc = f_int(std::fabs(xie-1.0));
        if (limit > 3 && intern > 5)  limit = 3;
        if (limit > 1 && intern > 10) limit = 1;
        jinc = i_min(jinc,limit);
        if (xie > 1.0) {
            jpc = jp + jinc;
        } else if (xie < 0.0) {
            jpc = jp - jinc;
        }
        jp = i_min(jpc,je-1);
        jp = i_max(js,jp);
        goto L6901;
    } else if (lsrch == 2) {
        // for fixed xie
        eta = (sc-a1-a2*xie)/(a3+a4*xie);

        int imiss = 0;
        if (eta < -tol_.epsc || eta > 1.0+tol_.epsc) {
            imiss = 1;
        }
        if (imiss == 0) goto L6902;

        // try new target cell
        int kpc = kp;
        int kinc = 0;
        if (eta >= 0.0) kinc = f_int(std::fabs(eta));
        if (eta <  0.0) kinc = f_int(std::fabs(eta-1.0));
        if (limit > 3 && intern > 5)  limit = 3;
        if (limit > 1 && intern > 10) limit = 1;
        kinc = i_min(kinc,limit);
        if (eta > 1.0) {
            kpc = kp + kinc;
        } else if (eta < 0.0) {
            kpc = kp - kinc;
        }
        kp = i_min(kpc,ke-1);
        kp = i_max(ks,kp);
        goto L6901;
    }

L6902:
    xiet = jp + xie;
    etat = kp + eta;
    return;
    #undef SARC
}
