// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// patcher.cpp — faithful C++ translation of patcher.F (CFL3D ronnie).
// Orchestrator: calculate patched-grid interpolation coefficients.
//
// SEQUENTIAL build: all `#if defined DIST_MPI` code (the whole IRecv/Send
// preamble and the MPI branch of the "from"-grid load) is REMOVED; only the
// non-MPI loadgr calls are kept. CMPLX is OFF (`real(x)` -> x).
//
// Array-index mapping (1-based column-major, CONTRACT rule 2):
//   lw(65,maxbl)            (i,nbl)  -> (i-1)+(nbl-1)*65
//   iindx(intmx,6*msub1+9)  (n,j)    -> (n-1)+(j-1)*intmx
//   isav_dpat(intmx,17)     (i,j)    -> (i-1)+(j-1)*intmx
//   isav_dpat_b(intmx,msub1,6)(i,j,k)-> (i-1)+(j-1)*intmx+(k-1)*intmx*msub1
//   windx(mxxe,2)           (i,j)    -> (i-1)+(j-1)*mxxe
//   factjlo/../dx/../dthetx(intmx,msub1)(i,l)-> (i-1)+(l-1)*intmx
//   jdimg/kdimg/idimg(maxbl)(nb)     -> (nb-1)
//   local 2-D work x1..(mptch+2,mptch+2)(j,k)-> (j-1)+(k-1)*(mptch+2)
//   xte/../(mptch+2,mptch+2,msub1)(j,k,l)->(j-1)+(k-1)*(mptch+2)+(l-1)*(mptch+2)^2
//   jimage/kimage(msub1,mptch+2,mptch+2)(l,j,k)->
//        (l-1)+(j-1)*msub1+(k-1)*msub1*(mptch+2)
#include "ron_common.h"

void patcher(int nbl,int* lw,double* w,int mgwk,double* wk,int nwork,int ncall,
             int ioutpt,int it_thro,int maxbl,int msub1,int intmx,int mxxe,int mptch,
             int* jdimg,int* kdimg,int* idimg,double* windx,int nintr,int* iindx,
             int* llimit,int* iitmax,int* mmcxie,int* mmceta,int* ncheck,int* iifit,
             int* mblkpt,int* iic0,int* iiorph,int* iitoss,int* ifiner,
             double* factjlo,double* factjhi,double* factklo,double* factkhi,
             double* dx,double* dy,double* dz,double* dthetx,double* dthety,double* dthetz,
             int* isav_dpat,int* isav_dpat_b,
             double* xte,double* yte,double* zte,double* xmi,double* ymi,double* zmi,
             double* xmie,double* ymie,double* zmie,int* jjmax1,int* kkmax1,
             int* jimage,int* kimage,double* xorig,double* yorig,double* zorig,
             int* jte,int* kte,double* sxie,double* seta,double* sxie2,double* seta2,
             double* xie2s,double* eta2s,double* temp,double* x2,double* y2,double* z2,
             int* nblk1,int* nblk2,int* jmm,int* kmm,double* x1,double* y1,double* z1,
             int* lout,int* xif1,int* xif2,int* etf1,int* etf2,int* ireq_ar,
             int myid,int myhost,int mycomm,int* mblk2nd,int* nou,FStr* bou,
             int nbuf,int ibufdim)
{
    (void)wk; (void)nwork; (void)nintr; (void)ncheck; (void)ireq_ar;
    (void)myhost; (void)mycomm; (void)nbuf;

    const int m2 = mptch+2;
    // 1-based accessors
    auto LW      = [&](int i,int nb){ return lw[(long)(i-1)+(long)(nb-1)*65]; };
    auto IINDX   = [&](int n,int j){ return iindx[(long)(n-1)+(long)(j-1)*intmx]; };
    auto ISDP    = [&](int i,int j){ return isav_dpat[(long)(i-1)+(long)(j-1)*intmx]; };
    auto ISDPB   = [&](int i,int j,int k){
        return isav_dpat_b[(long)(i-1)+(long)(j-1)*intmx+(long)(k-1)*intmx*msub1]; };
    auto FJLO    = [&](int i,int l)->double&{ return factjlo[(long)(i-1)+(long)(l-1)*intmx]; };
    auto FJHI    = [&](int i,int l)->double&{ return factjhi[(long)(i-1)+(long)(l-1)*intmx]; };
    auto FKLO    = [&](int i,int l)->double&{ return factklo[(long)(i-1)+(long)(l-1)*intmx]; };
    auto FKHI    = [&](int i,int l)->double&{ return factkhi[(long)(i-1)+(long)(l-1)*intmx]; };
    auto DXv     = [&](int i,int l){ return dx[(long)(i-1)+(long)(l-1)*intmx]; };
    auto DYv     = [&](int i,int l){ return dy[(long)(i-1)+(long)(l-1)*intmx]; };
    auto DZv     = [&](int i,int l){ return dz[(long)(i-1)+(long)(l-1)*intmx]; };
    auto DTX     = [&](int i,int l){ return dthetx[(long)(i-1)+(long)(l-1)*intmx]; };
    auto DTY     = [&](int i,int l){ return dthety[(long)(i-1)+(long)(l-1)*intmx]; };
    auto DTZ     = [&](int i,int l){ return dthetz[(long)(i-1)+(long)(l-1)*intmx]; };
    // windx(lst,col): pointer to base element windx(lst,col)
    auto WINDX   = [&](int lst,int col){ return &windx[(long)(lst-1)+(long)(col-1)*mxxe]; };

    char line[200];

    // scalar arrays jjmax2(1),kkmax2(1)
    int jjmax2[1], kkmax2[1];

    const int iavg = 0;

    //********** set tolerances, etc. (10.**(-iexp) is machine zero) **********
    const int iexp = zero_.iexp;
    (void)ncall;

    // expansion factors for "from" grids
    double factj = 0.01;
    double factk = 0.01;

    // convergence tolerance for generalized coordinates
    tol_.epsc = f_max(1.0e-07, std::pow(10.0,(double)(-iexp+1)));
    // threshold for collapsed boundaries  (patcher names slot3 `epscoll` = tol_.eps)
    tol_.eps  = f_max(1.0e-10, std::pow(10.0,(double)(-iexp+1)));
    // threshold for reentrant boundaries
    tol_.epsreen = f_max(1.0e-09, std::pow(10.0,(double)(-iexp+1)));
    // threshold for C-0 continuous interfaces
    tol_.epsc0 = f_max(1.0e-07, std::pow(10.0,(double)(-iexp+1)));

    sklt1_.isklt1 = 0;
    if (ioutpt > 0) sklt1_.isklt1 = 1;
    const int isklt1 = sklt1_.isklt1;

    // (it_thro.gt.1 -> go to 909; the DIST_MPI preamble between here and 909 is
    //  removed entirely for the sequential build, so nothing to skip.)
    (void)it_thro;

// 909 continue
    const int maxdcnt = is_dpatch_.maxdcnt;

    for (int lcnt = 1; lcnt <= maxdcnt; ++lcnt) {
        // nbll is the current (to) block; nbl is block coming into this routine
        int nbll = ISDP(lcnt,1);
        if (nbll != nbl) continue;

        int nd_dest = mblk2nd[nbl-1]; (void)nd_dest;
        int lmax1   = ISDP(lcnt,2);
        int icheck  = ISDP(lcnt,17);
        int lst    = IINDX(icheck,2*lmax1+5);
        int npt    = IINDX(icheck,2*lmax1+4);
        int ifit   = iifit[icheck-1];
        int limit0 = llimit[icheck-1];
        int itmax  = iitmax[icheck-1];
        int mcxie  = mmcxie[icheck-1];
        int mceta  = mmceta[icheck-1];
        int ic0    = iic0[icheck-1];
        int iorph  = iiorph[icheck-1];
        int itoss0 = iitoss[icheck-1];
        int xi1    = IINDX(icheck,2*lmax1+6);
        int xi2    = IINDX(icheck,2*lmax1+7);
        int et1    = IINDX(icheck,2*lmax1+8);
        int et2    = IINDX(icheck,2*lmax1+9);
        if (mcxie > 100) mcxie = 200;
        if (mceta > 100) mceta = 200;

        if (isklt1 > 0) {
            bou_put(nou,bou,ibufdim,4,"");
            bou_put(nou,bou,ibufdim,4,"");
            std::snprintf(line,sizeof(line),
                " generalized coordinate interpolation number %3d",icheck);
            bou_put(nou,bou,ibufdim,4,line);   // 721
        }

        // ---- describe each "from" sub-face ------------------------------
        for (int l = 1; l <= lmax1; ++l) {   // do 1605
            int mbl   = IINDX(icheck,l+1);
            int xif1l = IINDX(icheck,2*lmax1+9+l);
            int xif2l = IINDX(icheck,3*lmax1+9+l);
            int etf1l = IINDX(icheck,4*lmax1+9+l);
            int etf2l = IINDX(icheck,5*lmax1+9+l);
            int jindexl = IINDX(icheck,l+lmax1+2);

            // from an i=constant surface
            if (jindexl/10 == 1) {
                jjmax1[l-1] = jdimg[mbl-1];
                kkmax1[l-1] = kdimg[mbl-1];
                if (isklt1 > 0) {
                    int isurf = (jindexl == 11) ? 1 : idimg[mbl-1];
                    std::snprintf(line,sizeof(line),
                        " interpolation from surface i =%3d of block%3d"
                        "  (j: %3d to %3d, k: %3d to %3d)",
                        isurf,mbl,xif1l,xif2l,etf1l,etf2l);   // 1957
                    bou_put(nou,bou,ibufdim,4,line);
                }
            }
            // from a j=constant surface
            if (jindexl/10 == 2) {
                jjmax1[l-1] = kdimg[mbl-1];
                kkmax1[l-1] = idimg[mbl-1];
                if (isklt1 > 0) {
                    int isurf = (jindexl == 21) ? 1 : jdimg[mbl-1];
                    std::snprintf(line,sizeof(line),
                        " interpolation from surface j =%3d of block%3d"
                        "  (k: %3d to %3d, i: %3d to %3d)",
                        isurf,mbl,xif1l,xif2l,etf1l,etf2l);   // 1958
                    bou_put(nou,bou,ibufdim,4,line);
                }
            }
            // from a k=constant surface
            if (jindexl/10 == 3) {
                jjmax1[l-1] = jdimg[mbl-1];
                kkmax1[l-1] = idimg[mbl-1];
                if (isklt1 > 0) {
                    int isurf = (jindexl == 31) ? 1 : kdimg[mbl-1];
                    std::snprintf(line,sizeof(line),
                        " interpolation from surface k =%3d of block%3d"
                        "  (j: %3d to %3d, i: %3d to %3d)",
                        isurf,mbl,xif1l,xif2l,etf1l,etf2l);   // 1959
                    bou_put(nou,bou,ibufdim,4,line);
                }
            }

            if (jjmax1[l-1] > mptch || kkmax1[l-1] > mptch) {
                std::snprintf(line,sizeof(line),
                    " program terminated in dynamic patching routines"
                    " - see file %-60s", filenam_.grdmov.c());
                bou_put(nou,bou,ibufdim,1,line);
                std::snprintf(line,sizeof(line),
                    " stopping ... mptch =%11d too small ",mptch);
                bou_put(nou,bou,ibufdim,1,line);
                termn8(myid,-1,ibufdim,nbuf,bou,nou);
            }
        } // 1605

        // ---- describe the "to" surface ----------------------------------
        int jindex2 = IINDX(icheck,2*lmax1+3);
        // to an i=constant surface
        if (jindex2/10 == 1) {
            jjmax2[0] = jdimg[nbl-1];
            kkmax2[0] = kdimg[nbl-1];
            if (isklt1 > 0) {
                int isurf = (jindex2 == 11) ? 1 : idimg[nbl-1];
                std::snprintf(line,sizeof(line),
                    "                 to surface i =%3d of block%3d"
                    "  (j: %3d to %3d, k: %3d to %3d)",
                    isurf,nbl,xi1,xi2,et1,et2);   // 1857
                bou_put(nou,bou,ibufdim,4,line);
            }
        }
        // to a j=constant surface
        if (jindex2/10 == 2) {
            jjmax2[0] = kdimg[nbl-1];
            kkmax2[0] = idimg[nbl-1];
            if (isklt1 > 0) {
                int isurf = (jindex2 == 21) ? 1 : jdimg[nbl-1];
                std::snprintf(line,sizeof(line),
                    "                 to surface j =%3d of block%3d"
                    "  (k: %3d to %3d, i: %3d to %3d)",
                    isurf,nbl,xi1,xi2,et1,et2);   // 1858
                bou_put(nou,bou,ibufdim,4,line);
            }
        }
        // to a k=constant surface
        if (jindex2/10 == 3) {
            jjmax2[0] = jdimg[nbl-1];
            kkmax2[0] = idimg[nbl-1];
            if (isklt1 > 0) {
                int isurf = (jindex2 == 31) ? 1 : kdimg[nbl-1];
                std::snprintf(line,sizeof(line),
                    "                 to surface k =%3d of block%3d"
                    "  (j: %3d to %3d, i: %3d to %3d)",
                    isurf,nbl,xi1,xi2,et1,et2);   // 1859
                bou_put(nou,bou,ibufdim,4,line);
            }
        }

        if (jjmax2[0] > mptch || kkmax2[0] > mptch) {
            std::snprintf(line,sizeof(line),
                " program terminated in dynamic patching routines"
                " - see file %-60s", filenam_.grdmov.c());
            bou_put(nou,bou,ibufdim,1,line);
            std::snprintf(line,sizeof(line),
                " stopping ... mptch =%11d too small ",mptch);
            bou_put(nou,bou,ibufdim,1,line);
            termn8(myid,-1,ibufdim,nbuf,bou,nou);
        }

        if (isklt1 > 0) {
            bou_put(nou,bou,ibufdim,4,"");
            bou_put(nou,bou,ibufdim,4,
                " note: j and k referred to below are indicies"
                " local to the patch surface");   // 99
        }

        // set-up for "to" side of patch interface
        int jmax2  = jjmax2[0];
        int kmax2  = kkmax2[0];
        int lx2    = LW(10,nbl);
        int ly2    = LW(11,nbl);
        int lz2    = LW(12,nbl);
        int idim2g = idimg[nbl-1];
        int jdim2g = jdimg[nbl-1];
        int kdim2g = kdimg[nbl-1];
        int jindex = IINDX(icheck,2*lmax1+3);

        // load proper grid from 1-d array into 2-d work array.
        loadgr(w,mgwk,lx2,ly2,lz2,jindex,x2,y2,z2,m2,m2,idim2g,jdim2g,kdim2g);

        // check for collapsed grid lines ("to" side)
        if (isklt1 > 0) {
            std::snprintf(line,sizeof(line),
                "    checking for collapsed grid lines on  \"to\" "
                " side in block %3d",nbl);   // 3958
            bou_put(nou,bou,ibufdim,4,line);
        }
        collapse(m2,m2,jmax2,kmax2,x2,y2,z2,nou,bou,nbuf,ibufdim);

        // ---- set-up for "from" side of patch interface ------------------
        int iself = 0;
        for (int l = 1; l <= lmax1; ++l) {   // do 125
            int mbl    = IINDX(icheck,l+1);
            int lx1    = LW(10,mbl);
            int ly1    = LW(11,mbl);
            int lz1    = LW(12,mbl);
            int jmax1  = jjmax1[l-1];
            int kmax1  = kkmax1[l-1];
            int jindexf = IINDX(icheck,l+lmax1+2);
            int idim1  = idimg[mbl-1];
            int jdim1  = jdimg[mbl-1];
            int kdim1  = kdimg[mbl-1];
            // iself = 1 if a block face communicates with itself
            int itest1 = 100*nbl+IINDX(icheck,2*lmax1+3);
            int itest2 = 100*mbl+IINDX(icheck,l+lmax1+2);
            if (itest1 == itest2) iself = 1;

            // load proper grid from 1-d array into 2-d work array
            // (DIST_MPI branch removed; sequential loadgr only)
            loadgr(w,mgwk,lx1,ly1,lz1,jindexf,x1,y1,z1,m2,m2,idim1,jdim1,kdim1);

            // check for collapsed grid lines ("from" side)
            if (isklt1 > 0) {
                std::snprintf(line,sizeof(line),
                    "    checking for collapsed grid lines on "
                    "\"from\" side in block %3d",mbl);   // 3957
                bou_put(nou,bou,ibufdim,4,line);
            }
            collapse(m2,m2,jmax1,kmax1,x1,y1,z1,nou,bou,nbuf,ibufdim);

            // check for branch cuts
            if (isklt1 > 0) {
                std::snprintf(line,sizeof(line),
                    "    checking for branch cuts on \"from\" side"
                    " in block %3d",mbl);   // 98
                bou_put(nou,bou,ibufdim,4,line);
            }

            int xif1l = IINDX(icheck,2*lmax1+9+l);
            int xif2l = IINDX(icheck,3*lmax1+9+l);
            int etf1l = IINDX(icheck,4*lmax1+9+l);
            int etf2l = IINDX(icheck,5*lmax1+9+l);
            rechk(m2,m2,jimage,kimage,msub1,jmax1,kmax1,l,x1,y1,z1,
                  xif1l,xif2l,etf1l,etf2l,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl);

            // expand "from" grid(s) at boundaries
            if (isklt1 > 0) {
                std::snprintf(line,sizeof(line),
                    "    expanding grid boundaries on \"from\" side"
                    " in block %3d",mbl);   // 1755
                bou_put(nou,bou,ibufdim,4,line);
            }

            FJLO(icheck,l) = factj;
            FJHI(icheck,l) = factj;
            FKLO(icheck,l) = factk;
            FKHI(icheck,l) = factk;

            // should not really need any significant expansion for c-0 grids
            if (ic0 > 0) {
                FJLO(icheck,l) = 1.0e-5;
                FJHI(icheck,l) = 1.0e-5;
                FKLO(icheck,l) = 1.0e-5;
                FKHI(icheck,l) = 1.0e-5;
            }

            expand(m2,m2,msub1,jmax1,kmax1,l,x1,y1,z1,xte,yte,zte,
                   FJLO(icheck,l),FJHI(icheck,l),FKLO(icheck,l),FKHI(icheck,l),
                   jmax2,kmax2,x2,y2,z2);
            jte[l-1] = jjmax1[l-1]+2;
            kte[l-1] = kkmax1[l-1]+2;

            // translate/rotate "from" blocks as needed
            if (f_abs(DXv(icheck,l)) > 0.0 ||
                f_abs(DYv(icheck,l)) > 0.0 ||
                f_abs(DZv(icheck,l)) > 0.0) {
                int jjte = jte[l-1];
                int kkte = kte[l-1];
                transp(m2,m2,jjte,kkte,msub1,l,xte,yte,zte,dx,dy,dz,intmx,icheck);
            }
            if (f_abs(DTX(icheck,l)) > 0.0 ||
                f_abs(DTY(icheck,l)) > 0.0 ||
                f_abs(DTZ(icheck,l)) > 0.0) {
                int jjte = jte[l-1];
                int kkte = kte[l-1];
                rotatp(m2,m2,jjte,kkte,msub1,l,xte,yte,zte,dthetx,dthety,dthetz,
                       xorig,yorig,zorig,mbl,maxbl,intmx,icheck);
            }

            // search range on "from" side
            xif1[l-1] = IINDX(icheck,2*lmax1+9+l);
            if (xif1[l-1] > 1) xif1[l-1] = xif1[l-1] + 1;
            xif2[l-1] = IINDX(icheck,3*lmax1+9+l);
            if (xif2[l-1] == jjmax1[l-1]) {
                xif2[l-1] = xif2[l-1] + 2;
            } else {
                xif2[l-1] = xif2[l-1] + 1;
            }
            etf1[l-1] = IINDX(icheck,4*lmax1+9+l);
            if (etf1[l-1] > 1) etf1[l-1] = etf1[l-1] + 1;
            etf2[l-1] = IINDX(icheck,5*lmax1+9+l) + 1;
            if (etf2[l-1] == kkmax1[l-1]) {
                etf2[l-1] = etf2[l-1] + 2;
            } else {
                etf2[l-1] = etf2[l-1] + 1;
            }
        } // 125

        int lstf, nptf, xi1f, xi2f, et1f, et2f;
        if (ifiner[icheck-1] != 0) {
            int icheckf = ifiner[icheck-1];
            int lmax1f  = IINDX(icheckf,1);
            lstf = IINDX(icheckf,2*lmax1f+5);
            nptf = IINDX(icheckf,2*lmax1f+4);
            xi1f = IINDX(icheckf,2*lmax1f+6);
            xi2f = IINDX(icheckf,2*lmax1f+7);
            et1f = IINDX(icheckf,2*lmax1f+8);
            et2f = IINDX(icheckf,2*lmax1f+9);
        } else {
            lstf = IINDX(icheck,2*lmax1+5);
            nptf = IINDX(icheck,2*lmax1+4);
            xi1f = IINDX(icheck,2*lmax1+6);
            xi2f = IINDX(icheck,2*lmax1+7);
            et1f = IINDX(icheck,2*lmax1+8);
            et2f = IINDX(icheck,2*lmax1+9);
        }

        // obtain interpolation coefficients for coarser levels by averaging
        if (ifiner[icheck-1] != 0 && iavg > 0) {
            avgint(WINDX(lst,1),WINDX(lst,2),&mblkpt[lst-1],npt,
                   WINDX(lstf,1),WINDX(lstf,2),&mblkpt[lstf-1],nptf,
                   xi1,xi2,et1,et2,xi1f,xi2f,et1f,et2f);
        }

        // determine projection of x2,y2,z2 onto "from" generalized coords
        invert(m2,m2,msub1,1,jte,kte,lmax1,xte,yte,zte,xmi,ymi,zmi,
               xmie,ymie,zmie,limit0,jjmax2,kkmax2,x2,y2,z2,
               WINDX(lst,1),WINDX(lst,2),&mblkpt[lst-1],temp,jimage,kimage,
               ifit,itmax,sxie,seta,sxie2,seta2,xie2s,eta2s,
               intmx,icheck,nblk1,nblk2,jmm,kmm,mcxie,mceta,
               lout,xi1,xi2,et1,et2,npt,ic0,iorph,itoss0,ncall,ioutpt,
               xif1,xif2,etf1,etf2,iself,ifiner[icheck-1],
               WINDX(lstf,1),WINDX(lstf,2),&mblkpt[lstf-1],nptf,
               xi1f,xi2f,et1f,et2f,iavg,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl);
    }
}
