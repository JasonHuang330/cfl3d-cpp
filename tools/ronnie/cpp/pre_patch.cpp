// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// pre_patch.cpp — set up data arrays for patch bc message passing.
// lw(65,maxbl); jdimg/kdimg/idimg(maxbl); isav_pat(intmax,17);
// isav_pat_b(intmax,nsub1,6); jjmax1/kkmax1/iiint1/iiint2(nsub1);
// iindex(intmax,6*nsub1+9).  All column-major, 1-based; icount is inout.
#include "ron_common.h"

void pre_patch(int nbl,int* lw,int& icount,int ninter,int* iindex,int intmax,
               int nsub1,int* isav_pat,int* isav_pat_b,int* jjmax1,int* kkmax1,
               int* iiint1,int* iiint2,int /*maxbl*/,int* jdimg,int* kdimg,
               int* idimg,int /*ierrflg*/){
    auto LW =[&](int i)->int{ return lw[(i-1)+(nbl-1)*65]; };
    auto II =[&](int r,int c)->int&{ return iindex[(r-1)+(c-1)*intmax]; };
    auto SP =[&](int r,int c)->int&{ return isav_pat[(r-1)+(c-1)*intmax]; };
    // isav_pat_b(intmax,nsub1,6): (r-1)+(l-1)*intmax+(c-1)*intmax*nsub1
    auto SPB=[&](int r,int l,int c)->int&{
        return isav_pat_b[(r-1)+(l-1)*intmax+(c-1)*intmax*nsub1]; };

    if (f_iabs(ninter)>0){
        for(int icheck=1;icheck<=f_iabs(ninter);++icheck){
            int lmax1 = II(icheck,1);
            int nblcc = II(icheck,lmax1+2);
            if (nblcc!=nbl) continue;
            int lst = II(icheck,2*lmax1+5);
            int npt = II(icheck,2*lmax1+4);

            int j21 = II(icheck,2*lmax1+6);
            int j22 = II(icheck,2*lmax1+7);
            int k21 = II(icheck,2*lmax1+8);
            int k22 = II(icheck,2*lmax1+9);

            int lqedge=0,lqedgt=0,lqedgv=0,lqedgb=0,jmax2=0,kmax2=0,mint1=0,mint2=0;
            int toty = II(icheck,2*lmax1+3);
            if (toty/10==1){
                lqedge=LW(4); lqedgt=LW(25); lqedgv=LW(30); lqedgb=LW(33);
                jmax2=jdimg[nbl-1]; kmax2=kdimg[nbl-1];
                if (toty==11){ mint1=1; mint2=2; } else { mint1=3; mint2=4; }
            }
            if (toty/10==2){
                lqedge=LW(2); lqedgt=LW(23); lqedgv=LW(28); lqedgb=LW(31);
                jmax2=kdimg[nbl-1]; kmax2=idimg[nbl-1]-1;
                if (toty==21){ mint1=1; mint2=2; } else { mint1=3; mint2=4; }
            }
            if (toty/10==3){
                lqedge=LW(3); lqedgt=LW(24); lqedgv=LW(29); lqedgb=LW(32);
                jmax2=jdimg[nbl-1]; kmax2=idimg[nbl-1]-1;
                if (toty==31){ mint1=1; mint2=2; } else { mint1=3; mint2=4; }
            }

            icount = icount + 1;
            SP(icount,1) =nbl;
            SP(icount,2) =lmax1;
            SP(icount,3) =j21;
            SP(icount,4) =j22;
            SP(icount,5) =k21;
            SP(icount,6) =k22;
            SP(icount,7) =lqedge;
            SP(icount,8) =lqedgt;
            SP(icount,9) =lqedgv;
            SP(icount,10)=lqedgb;
            SP(icount,11)=jmax2;
            SP(icount,12)=kmax2;
            SP(icount,13)=mint1;
            SP(icount,14)=mint2;
            SP(icount,15)=lst;
            SP(icount,16)=npt;
            SP(icount,17)=icheck;

            for(int l=1;l<=lmax1;++l){
                int mbl   = II(icheck,l+1);
                int mtype = II(icheck,l+lmax1+2);
                if (mtype/10==1){
                    jjmax1[l-1]=jdimg[mbl-1]; kkmax1[l-1]=kdimg[mbl-1];
                    if (mtype==11){ iiint1[l-1]=1; iiint2[l-1]=i_min(2,idimg[mbl-1]-1); }
                    else { iiint1[l-1]=idimg[mbl-1]-1; iiint2[l-1]=i_max(1,idimg[mbl-1]-2); }
                }
                if (mtype/10==2){
                    jjmax1[l-1]=kdimg[mbl-1]; kkmax1[l-1]=idimg[mbl-1];
                    if (mtype==21){ iiint1[l-1]=1; iiint2[l-1]=i_min(2,jdimg[mbl-1]-1); }
                    else { iiint1[l-1]=jdimg[mbl-1]-1; iiint2[l-1]=i_max(1,jdimg[mbl-1]-2); }
                }
                if (mtype/10==3){
                    jjmax1[l-1]=jdimg[mbl-1]; kkmax1[l-1]=idimg[mbl-1];
                    if (mtype==31){ iiint1[l-1]=1; iiint2[l-1]=i_min(2,kdimg[mbl-1]-1); }
                    else { iiint1[l-1]=kdimg[mbl-1]-1; iiint2[l-1]=i_max(1,kdimg[mbl-1]-2); }
                }
                SPB(icount,l,1)=mbl;
                SPB(icount,l,2)=mtype/10;
                SPB(icount,l,3)=jjmax1[l-1];
                SPB(icount,l,4)=kkmax1[l-1];
                SPB(icount,l,5)=iiint1[l-1];
                SPB(icount,l,6)=iiint2[l-1];
            }
        }
    }
}
