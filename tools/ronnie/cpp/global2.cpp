// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// global2.cpp — read dynamic patch input parameters (unit iunit25 output,
// iunit5 input). Faithful translation; ioflag=0 path is the tested one.
//
// Array shapes (column-major, 1-based):
//   iindx(intmx,6*msub1+9), llimit/iitmax/mmcxie/mmceta/iifit/iic0/iiorph/
//   iitoss/ifiner(intmx), ncheck(maxbl), dx.. dthetz(intmx,msub1),
//   idimg/jdimg/kdimg/levelg(maxbl), ncgg/nblg/iemg(maxgr), igridg(maxbl),
//   ifrom/xif1/xif2/etf1/etf2(msub1), nou(nbuf), bou(ibufdim,nbuf).
#include "ron_common.h"

void global2(int maxbl,int maxgr,int msub1,int& nintr,int intmx,int ngrid,
             int* idimg,int* jdimg,int* kdimg,int* levelg,int* ncgg,int* nblg,
             int* iindx,int* llimit,int* iitmax,int* mmcxie,int* mmceta,
             int* ncheck,int* iifit,int* iic0,int* iiorph,int* iitoss,int* ifiner,
             double* dx,double* dy,double* dz,double* dthetx,double* dthety,
             double* dthetz,int myid,int mptch,int mxxe,int icall,int iunit25,
             int* nou,FStr* bou,int ibufdim,int nbuf,int* ifrom,int* xif1,int* etf1,
             int* xif2,int* etf2,int* igridg,int* iemg,int nblock,int ioflag,int imode){
    (void)maxgr;(void)ncheck;
    FILE* f5=FUNIT(unit5_.iunit5);
    FILE* f25=FUNIT(iunit25);
    tc::TextReader tr(f5);

    // accessors (1-based)
    auto IX =[&](int r,int c)->int&{ return iindx[(r-1)+(c-1)*intmx]; };
    auto DX =[&](int r,int c)->double&{ return dx[(r-1)+(c-1)*intmx]; };
    auto DY =[&](int r,int c)->double&{ return dy[(r-1)+(c-1)*intmx]; };
    auto DZ =[&](int r,int c)->double&{ return dz[(r-1)+(c-1)*intmx]; };
    auto DTX=[&](int r,int c)->double&{ return dthetx[(r-1)+(c-1)*intmx]; };
    auto DTY=[&](int r,int c)->double&{ return dthety[(r-1)+(c-1)*intmx]; };
    auto DTZ=[&](int r,int c)->double&{ return dthetz[(r-1)+(c-1)*intmx]; };

    int ierrflg = (icall==0)? -99 : -1;

    nintr = 0;

    if (imode==1){
        tr.line();   // read(iunit5,*,end=999)  -- header (imode==1 only)
    }
    tr.line();       // read(iunit5,*)   skip
    { int tmp[1]; tr.ints(tmp,1); nintr=tmp[0]; }  // read(iunit5,*)nintr

    if (nintr==0) return;

    params_.lnbuf = params_.lnbuf + 2;

    if (ioflag==2){
        std::fprintf(f25,"\n dynamic patch input data\n");
        std::fprintf(f25," nintr\n");
    } else {
        std::fprintf(f25,"\npatch input data\n");
        std::fprintf(f25," ninter\n");
    }

    // ncgmax over grids
    int ncgmax=0;
    for(int igrid=1;igrid<=ngrid;++igrid){
        int ncg=ncgg[igrid-1];
        if (igrid==1) ncgmax=ncg; else ncgmax=i_max(ncgmax,ncg);
    }
    int ntest=(ncgmax+1)*nintr;
    if (ntest>intmx){
        std::fprintf(f25," too many blocks to interpolate; intmx should be at"
                     " least %s\n",f_iw(ntest,4).c_str());
        termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
    }
    // format 36: (13i6)
    std::fprintf(f25,"%s\n",f_iw(nintr,6).c_str());

    tr.line();   // read(iunit5,*)  (blank before iifit table)
    if (ioflag==2){
        std::fprintf(f25,"   int iifit    limit    itmax    mcxie    mceta"
                         "      c-0    iorph    itoss\n");
    } else {
        std::fprintf(f25,"   int iifit    limit    itmax    mcxie    mceta"
                         "      c-0    iorph\n");
    }
    if (nintr>0){
        for(int n=1;n<=nintr;++n) iifit[n-1]=-99;
        for(int n=1;n<=nintr;++n){
            int intv,iif;
            if (ioflag==2){
                int v[9]; tr.ints(v,9);
                intv=v[0]; iif=v[1];
                llimit[intv-1]=v[2]; iitmax[intv-1]=v[3]; mmcxie[intv-1]=v[4];
                mmceta[intv-1]=v[5]; iic0[intv-1]=v[6]; iiorph[intv-1]=v[7];
                iitoss[intv-1]=v[8];
            } else {
                int v[8]; tr.ints(v,8);
                intv=v[0]; iif=v[1];
                llimit[intv-1]=v[2]; iitmax[intv-1]=v[3]; mmcxie[intv-1]=v[4];
                mmceta[intv-1]=v[5]; iic0[intv-1]=v[6]; iiorph[intv-1]=v[7];
                iitoss[intv-1]=0;
            }
            if (iifit[intv-1]!=-99){
                FILE* f11=FUNIT(11);
                if(f11){ std::string g=filenam_.grdmov.s; g.resize(60,' ');
                    std::fprintf(f11," program terminated in dynamic patching "
                        "routines - see file %s\n",g.c_str()); }
                std::fprintf(f25," stopping...attempting to set data for int %d"
                             " more than once\n",intv);
                termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
            }
            iifit[intv-1]=iif;
            if (iifit[intv-1]==0) iifit[intv-1]=1;
            if (iic0[intv-1]>0){
                if (iifit[intv-1]>=0) iifit[intv-1]=1; else iifit[intv-1]=-1;
                mmcxie[intv-1]=0; mmceta[intv-1]=0;
            }
            // format 38: (i6,i6,10i9)
            if (ioflag==2){
                std::fprintf(f25,"%s%s%s%s%s%s%s%s%s\n",
                    f_iw(intv,6).c_str(),f_iw(iifit[intv-1],6).c_str(),
                    f_iw(llimit[intv-1],9).c_str(),f_iw(iitmax[intv-1],9).c_str(),
                    f_iw(mmcxie[intv-1],9).c_str(),f_iw(mmceta[intv-1],9).c_str(),
                    f_iw(iic0[intv-1],9).c_str(),f_iw(iiorph[intv-1],9).c_str(),
                    f_iw(iitoss[intv-1],9).c_str());
            } else {
                std::fprintf(f25,"%s%s%s%s%s%s%s%s\n",
                    f_iw(intv,6).c_str(),f_iw(iifit[intv-1],6).c_str(),
                    f_iw(llimit[intv-1],9).c_str(),f_iw(iitmax[intv-1],9).c_str(),
                    f_iw(mmcxie[intv-1],9).c_str(),f_iw(mmceta[intv-1],9).c_str(),
                    f_iw(iic0[intv-1],9).c_str(),f_iw(iiorph[intv-1],9).c_str());
            }
        }
    }

    int nfbmax=0;
    for(int n=1;n<=nintr;++n) IX(n,1)=-99;

    // format 8146/8147/8148/8149 headers
    auto hdr8146=[&](){
        std::fprintf(f25,"   int    to  xie1  xie2  eta1  eta2   nfb"
                         "  from  from  from  from   ...\n"); };
    auto hdr8147=[&](){
        std::fprintf(f25,"   int    to     xie1     xie2     eta1     eta2"
                         "      nfb\n"); };
    auto hdr8148=[&](){
        std::fprintf(f25,"          dx       dy       dz   dthetx   dthety"
                         "   dthetz\n"); };
    auto hdr8149=[&](){
        std::fprintf(f25,"        from     xie1     xie2     eta1     eta2\n"); };

    if (ioflag==0){
        tr.line();        // read(iunit5,*)
        hdr8146();
    } else if (ioflag==1){
        tr.line(); tr.line();
        hdr8147();
        hdr8149();
    }

    int intv=0,ito=0,xi1=0,xi2=0,et1=0,et2=0,nfb=0;
    for(int n=1;n<=nintr;++n){
        if (ioflag==2){
            tr.line();
            int v[7]; tr.ints(v,7);
            intv=v[0]; ito=v[1]; xi1=v[2]; xi2=v[3]; et1=v[4]; et2=v[5]; nfb=v[6];
        } else if (ioflag==1){
            int v[7]; tr.ints(v,7);
            intv=v[0]; ito=v[1]; xi1=v[2]; xi2=v[3]; et1=v[4]; et2=v[5]; nfb=v[6];
        } else { // ioflag==0
            // read(iunit5,*) int,ito,xi1,xi2,et1,et2,nfb,(ifrom(l),l=1,nfb)
            std::string l=tr.line(); std::stringstream ss(l);
            auto nexti=[&](int& out){ std::string t; if(ss>>t){
                for(char&c:t) if(c==',')c=' '; out=std::atoi(t.c_str()); } else out=0; };
            nexti(intv); nexti(ito); nexti(xi1); nexti(xi2); nexti(et1); nexti(et2);
            nexti(nfb);
            for(int lp=1;lp<=nfb;++lp) nexti(ifrom[lp-1]);
            for(int lp=1;lp<=nfb;++lp){
                xif1[lp-1]=0; xif2[lp-1]=0; etf1[lp-1]=0; etf2[lp-1]=0;
                DX(intv,lp)=0; DY(intv,lp)=0; DZ(intv,lp)=0;
                DTX(intv,lp)=0; DTY(intv,lp)=0; DTZ(intv,lp)=0;
            }
        }
        if (nfb>msub1){
            std::fprintf(f25," msub1 is too small; should be least: %s\n",
                         f_iw(nfb,2).c_str());
            termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
        }
        nfbmax=i_max(nfbmax,nfb);

        if (ioflag==1){
            for(int lp=1;lp<=nfb;++lp){
                int v[5]; tr.ints(v,5);
                ifrom[lp-1]=v[0]; xif1[lp-1]=v[1]; xif2[lp-1]=v[2];
                etf1[lp-1]=v[3]; etf2[lp-1]=v[4];
                DX(intv,lp)=0; DY(intv,lp)=0; DZ(intv,lp)=0;
                DTX(intv,lp)=0; DTY(intv,lp)=0; DTZ(intv,lp)=0;
            }
        } else if (ioflag==2){
            for(int lp=1;lp<=nfb;++lp){
                tr.line();
                int v[5]; tr.ints(v,5);
                ifrom[lp-1]=v[0]; xif1[lp-1]=v[1]; xif2[lp-1]=v[2];
                etf1[lp-1]=v[3]; etf2[lp-1]=v[4];
                tr.line();
                double rv[6]; tr.dbls(rv,6);
                DX(intv,lp)=rv[0]; DTX(intv,lp)=rv[3];
                if (igrdtyp_.ialph>0){
                    DY(intv,lp)=-rv[2]; DTY(intv,lp)=-rv[5];
                    DZ(intv,lp)= rv[1]; DTZ(intv,lp)= rv[4];
                } else {
                    DY(intv,lp)= rv[1]; DTY(intv,lp)= rv[4];
                    DZ(intv,lp)= rv[2]; DTZ(intv,lp)= rv[5];
                }
            }
        }

        if (IX(intv,1)!=-99){
            FILE* f11=FUNIT(11);
            if(f11){ std::string g=filenam_.grdmov.s; g.resize(60,' ');
                std::fprintf(f11," program terminated in dynamic patching "
                    "routines - see file %s\n",g.c_str()); }
            std::fprintf(f25," stopping...attempting to set data for int %d"
                         " more than once\n",intv);
            termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
        }
        IX(intv,1)=nfb;
        for(int lp=1;lp<=nfb;++lp){
            IX(intv,lp+1)     = ifrom[lp-1]/100;
            IX(intv,lp+nfb+2) = ifrom[lp-1]-IX(intv,lp+1)*100;
        }

        IX(intv,nfb+2) = ito/100;
        int nbl = nblg[IX(intv,nfb+2)-1];
        IX(intv,2*nfb+3) = ito - IX(intv,nfb+2)*100;

        // set full ranges if zeros are input ("to" side)
        int toty = IX(intv,2*nfb+3)/10;
        if (toty==1){
            if (xi1==0 && xi2==0){ xi1=1; xi2=jdimg[nbl-1]; }
            if (et1==0 && et2==0){ et1=1; et2=kdimg[nbl-1]; }
        }
        if (toty==2){
            if (xi1==0 && xi2==0){ xi1=1; xi2=kdimg[nbl-1]; }
            if (et1==0 && et2==0){ et1=1; et2=idimg[nbl-1]; }
        }
        if (toty==3){
            if (xi1==0 && xi2==0){ xi1=1; xi2=jdimg[nbl-1]; }
            if (et1==0 && et2==0){ et1=1; et2=idimg[nbl-1]; }
        }

        // set full ranges if zeros input for search range in "from" block(s)
        for(int lp=1;lp<=nfb;++lp){
            int mbl=nblg[IX(intv,lp+1)-1];
            int fty=IX(intv,lp+nfb+2)/10;
            if (fty==1){
                if (xif1[lp-1]==0 && xif2[lp-1]==0){ xif1[lp-1]=1; xif2[lp-1]=jdimg[mbl-1]; }
                if (etf1[lp-1]==0 && etf2[lp-1]==0){ etf1[lp-1]=1; etf2[lp-1]=kdimg[mbl-1]; }
            }
            if (fty==2){
                if (xif1[lp-1]==0 && xif2[lp-1]==0){ xif1[lp-1]=1; xif2[lp-1]=kdimg[mbl-1]; }
                if (etf1[lp-1]==0 && etf2[lp-1]==0){ etf1[lp-1]=1; etf2[lp-1]=idimg[mbl-1]; }
            }
            if (fty==3){
                if (xif1[lp-1]==0 && xif2[lp-1]==0){ xif1[lp-1]=1; xif2[lp-1]=jdimg[mbl-1]; }
                if (etf1[lp-1]==0 && etf2[lp-1]==0){ etf1[lp-1]=1; etf2[lp-1]=idimg[mbl-1]; }
            }
        }

        // format 38: (i6,i6,10i9); format 40: (6x,i6,10i9)
        auto w38=[&](int a,int b,int c,int d,int e,int ff,int g){
            std::fprintf(f25,"%s%s%s%s%s%s%s\n",
                f_iw(a,6).c_str(),f_iw(b,6).c_str(),f_iw(c,9).c_str(),
                f_iw(d,9).c_str(),f_iw(e,9).c_str(),f_iw(ff,9).c_str(),
                f_iw(g,9).c_str()); };
        auto w40=[&](int a,int b,int c,int d,int e){
            std::fprintf(f25,"      %s%s%s%s%s\n",
                f_iw(a,6).c_str(),f_iw(b,9).c_str(),f_iw(c,9).c_str(),
                f_iw(d,9).c_str(),f_iw(e,9).c_str()); };
        auto w39=[&](double a,double b,double c,double d,double e,double ff){
            // format 39: (3x,6f9.4)
            std::fprintf(f25,"   %s%s%s%s%s%s\n",
                f_fw(a,9,4).c_str(),f_fw(b,9,4).c_str(),f_fw(c,9,4).c_str(),
                f_fw(d,9,4).c_str(),f_fw(e,9,4).c_str(),f_fw(ff,9,4).c_str()); };

        if (ioflag==0 || ioflag==1){
            w38(intv,ito,xi1,xi2,et1,et2,nfb);
            for(int lp=1;lp<=nfb;++lp)
                w40(ifrom[lp-1],xif1[lp-1],xif2[lp-1],etf1[lp-1],etf2[lp-1]);
        } else {
            hdr8147();
            w38(intv,ito,xi1,xi2,et1,et2,nfb);
            for(int lp=1;lp<=nfb;++lp){
                hdr8149();
                w40(ifrom[lp-1],xif1[lp-1],xif2[lp-1],etf1[lp-1],etf2[lp-1]);
                hdr8148();
                w39(DX(intv,lp),DY(intv,lp),DZ(intv,lp),
                    DTX(intv,lp),DTY(intv,lp),DTZ(intv,lp));
                DTX(intv,lp)=DTX(intv,lp)/conversion_.radtodeg;
                DTY(intv,lp)=DTY(intv,lp)/conversion_.radtodeg;
                DTZ(intv,lp)=DTZ(intv,lp)/conversion_.radtodeg;
            }
        }

        // check input ranges of xie and eta ("to" side)
        auto errxie=[&](const char* face){
            std::fprintf(f25," incorrect xie range input on this %s "
                         "\"to\" face\n",face);
            termn8(myid,ierrflg,ibufdim,nbuf,bou,nou); };
        auto erreta=[&](const char* face){
            std::fprintf(f25," incorrect eta range input on this %s "
                         "\"to\" face\n",face);
            termn8(myid,ierrflg,ibufdim,nbuf,bou,nou); };
        if (toty==1){
            if (xi1<1 || xi1>jdimg[nbl-1]) errxie("i=constant");
            if (xi2<1 || xi2>jdimg[nbl-1]) errxie("i=constant");
            if (xi2<=xi1) errxie("i=constant");
            if (et1<1 || et1>kdimg[nbl-1]) erreta("i=constant");
            if (et2<1 || et2>kdimg[nbl-1]) erreta("i=constant");
            if (et2<=et1) erreta("i=constant");
        }
        if (toty==2){
            if (xi1<1 || xi1>kdimg[nbl-1]) errxie("j=constant");
            if (xi2<1 || xi2>kdimg[nbl-1]) errxie("j=constant");
            if (xi2<=xi1) errxie("j=constant");
            if (et1<1 || et1>idimg[nbl-1]) erreta("j=constant");
            if (et2<1 || et2>idimg[nbl-1]) erreta("j=constant");
            if (et2<=et1) erreta("j=constant");
        }
        if (toty==3){
            if (xi1<1 || xi1>jdimg[nbl-1]) errxie("k=constant");
            if (xi2<1 || xi2>jdimg[nbl-1]) errxie("k=constant");
            if (xi2<=xi1) errxie("k=constant");
            if (et1<1 || et1>idimg[nbl-1]) erreta("k=constant");
            if (et2<1 || et2>idimg[nbl-1]) erreta("k=constant");
            if (et2<=et1) erreta("k=constant");
        }

        // check input search ranges - "from" block(s)
        for(int lp=1;lp<=nfb;++lp){
            int mbl=nblg[IX(intv,lp+1)-1];
            int fty=IX(intv,lp+nfb+2)/10;
            auto exs=[&](const char* face){
                std::fprintf(f25," incorrect xie search range input on the %s "
                    "\"from\" face%s\n",face,f_iw(lp,12).c_str());
                termn8(myid,ierrflg,ibufdim,nbuf,bou,nou); };
            auto ees=[&](const char* face){
                std::fprintf(f25," incorrect eta search range input on the %s "
                    "\"from\" face%s\n",face,f_iw(lp,12).c_str());
                termn8(myid,ierrflg,ibufdim,nbuf,bou,nou); };
            if (fty==1){
                if (xif1[lp-1]<1 || xif1[lp-1]>jdimg[mbl-1]) exs("i=constant");
                if (xif2[lp-1]<1 || xif2[lp-1]>jdimg[mbl-1]) exs("i=constant");
                if (xif2[lp-1]<=xif1[lp-1]) exs("i=constant");
                if (etf1[lp-1]<1 || etf1[lp-1]>kdimg[mbl-1]) ees("i=constant");
                if (etf2[lp-1]<1 || etf2[lp-1]>kdimg[mbl-1]) ees("i=constant");
                if (etf2[lp-1]<=etf1[lp-1]) ees("i=constant");
            }
            if (fty==2){
                if (xif1[lp-1]<1 || xif1[lp-1]>kdimg[mbl-1]) exs("j=constant");
                if (xif2[lp-1]<1 || xif2[lp-1]>kdimg[mbl-1]) exs("j=constant");
                if (xif2[lp-1]<=xif1[lp-1]) exs("j=constant");
                if (etf1[lp-1]<1 || etf1[lp-1]>idimg[mbl-1]) ees("j=constant");
                if (etf2[lp-1]<1 || etf2[lp-1]>idimg[mbl-1]) ees("j=constant");
                if (etf2[lp-1]<=etf1[lp-1]) ees("j=constant");
            }
            if (fty==3){
                if (xif1[lp-1]<1 || xif1[lp-1]>jdimg[mbl-1]) exs("k=constant");
                if (xif2[lp-1]<1 || xif2[lp-1]>jdimg[mbl-1]) exs("k=constant");
                if (xif2[lp-1]<=xif1[lp-1]) exs("k=constant");
                if (etf1[lp-1]<1 || etf1[lp-1]>idimg[mbl-1]) ees("k=constant");
                if (etf2[lp-1]<1 || etf2[lp-1]>idimg[mbl-1]) ees("k=constant");
                if (etf2[lp-1]<=etf1[lp-1]) ees("k=constant");
            }
        }

        IX(intv,2*nfb+6)=xi1;
        IX(intv,2*nfb+7)=xi2;
        IX(intv,2*nfb+8)=et1;
        IX(intv,2*nfb+9)=et2;
        for(int lp=1;lp<=nfb;++lp){
            IX(intv,2*nfb+9+lp)=xif1[lp-1];
            IX(intv,3*nfb+9+lp)=xif2[lp-1];
            IX(intv,4*nfb+9+lp)=etf1[lp-1];
            IX(intv,5*nfb+9+lp)=etf2[lp-1];
        }
        IX(intv,2*nfb+4)=(xi2-xi1)*(et2-et1);
        if (intv==1){
            IX(intv,2*nfb+5)=1;
        } else {
            int nfb1=IX(intv-1,1);
            IX(intv,2*nfb+5)=IX(intv-1,2*nfb1+5)+IX(intv-1,2*nfb1+4);
        }
    } // end do 1703

    if (imode==0){
        std::fprintf(f25,"\n SUMMARY BY GRIDS\n");
        std::fprintf(f25,"   grid  level  block   jdim   kdim   idim"
                         "   grid pts.\n");
        int lfem=0, lfgm=ncgmax+1, lcgm=1, igptot=0;
        for(int igrid=1;igrid<=ngrid;++igrid){
            int iem=iemg[igrid-1];
            int nbl=nblg[igrid-1];
            int ncg=ncgg[igrid-1];
            igridg[nbl-1]=igrid;
            levelg[nbl-1]=iem+lfgm;
            if (iem==0) lfgm=i_max(lfgm,levelg[nbl-1]);
            if (iem>0)  lfem=i_max(lfem,levelg[nbl-1]);
            int igpts=jdimg[nbl-1]*kdimg[nbl-1]*idimg[nbl-1];
            igptot+=igpts;
            // format 7: (6i7,i12)
            auto w7=[&](int a,int b,int c,int d,int e,int ff,int g){
                std::fprintf(f25,"%s%s%s%s%s%s%s\n",
                    f_iw(a,7).c_str(),f_iw(b,7).c_str(),f_iw(c,7).c_str(),
                    f_iw(d,7).c_str(),f_iw(e,7).c_str(),f_iw(ff,7).c_str(),
                    f_iw(g,12).c_str()); };
            w7(igrid,levelg[nbl-1],nbl,jdimg[nbl-1],kdimg[nbl-1],idimg[nbl-1],igpts);
            if (ncg>0 && iem==0){
                for(int nn=1;nn<=ncg;++nn){
                    nbl=nbl+1;
                    igridg[nbl-1]=igrid;
                    levelg[nbl-1]=levelg[nbl-2]-1;
                    igpts=jdimg[nbl-1]*kdimg[nbl-1]*idimg[nbl-1];
                    igptot+=igpts;
                    w7(igrid,levelg[nbl-1],nbl,jdimg[nbl-1],kdimg[nbl-1],
                       idimg[nbl-1],igpts);
                }
            }
        }
        // format 2021: (/,37x,5hTOTAL,i12)
        std::fprintf(f25,"\n%sTOTAL%s\n",std::string(37,' ').c_str(),
                     f_iw(igptot,12).c_str());

        std::fprintf(f25,"\n SUMMARY BY LEVELS\n");
        std::fprintf(f25,"  level   grid  block\n");
        int lf=lfem; if (lfem==0) lf=lfgm;
        for(int levelc=lf;levelc>=lcgm;--levelc){
            for(int nbl=1;nbl<=nblock;++nbl){
                if (levelc!=levelg[nbl-1]) continue;
                int igrid=igridg[nbl-1];
                std::fprintf(f25,"%s%s%s\n",f_iw(levelg[nbl-1],7).c_str(),
                             f_iw(igrid,7).c_str(),f_iw(nbl,7).c_str());
            }
        }

        if (info_.mseq>(ncgmax+1)){
            std::fprintf(f25," mseq,ncgmax=%s%s\n",f_iw(info_.mseq,4).c_str(),
                         f_iw(ncgmax,4).c_str());
            std::fprintf(f25," error in input, mseq, ncgmax%s%s\n",
                         f_iw(info_.mseq,5).c_str(),f_iw(ncgmax,5).c_str());
            termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
        }

        // format 4008: (/)  -> blank line
        std::fprintf(f25,"\n\n");
        if (lfem!=0)
            std::fprintf(f25," level of finest embedded mesh=%s\n",f_iw(lfem,3).c_str());
        std::fprintf(f25," level of coarsest global mesh=%s\n",f_iw(lcgm,3).c_str());
        std::fprintf(f25," level of finest global mesh  =%s\n",f_iw(lfgm,3).c_str());

        std::fprintf(f25,"\n SUMMARY OF GRID SEQUENCES\n");
        std::fprintf(f25," sequence  starting level    ending level\n");
        for(int m=1;m<=info_.mseq;++m){
            // format 29: (2x,i7,2i16)
            std::fprintf(f25,"  %s%s%s\n",f_iw(m,7).c_str(),
                         f_iw(info_.levelt[m-1],16).c_str(),
                         f_iw(info_.levelb[m-1],16).c_str());
        }
    }

    // generate additional patch data for coarser blocks
    int ntemp=nintr;
    int ifirst=1;
    for(int n=1;n<=nintr;++n){
        nfb=IX(n,1);
        ifiner[n-1]=0;
        if (ncgg[IX(n,nfb+2)-1]>0){
            if (ifirst==1){
                if (ioflag==2){
                    std::fprintf(f25,"\n REARRANGED BLOCK NUMBERS AND DYNAMIC "
                        "PATCH DATA GENERATED FOR COARSER MESHES\n"
                        "               block references to CFL3D internal ordering\n"
                        "                    * ---> coarser level interpolation\n\n");
                } else {
                    std::fprintf(f25,"\n REARRANGED BLOCK NUMBERS AND PATCH DATA "
                        "GENERATED FOR COARSER MESHES\n"
                        "               block references to CFL3D internal ordering\n"
                        "                    * ---> coarser level interpolation\n\n");
                }
                ifirst=0;
            }
            ito=(IX(n,nfb+2)*100)+IX(n,2*nfb+3);
            int itonew=nblg[IX(n,nfb+2)-1]*100+IX(n,2*nfb+3);
            xi1=IX(n,2*nfb+6); xi2=IX(n,2*nfb+7);
            et1=IX(n,2*nfb+8); et2=IX(n,2*nfb+9);
            for(int i=1;i<=nfb;++i)
                ifrom[i-1]=nblg[IX(n,i+1)-1]*100+IX(n,i+nfb+2);
            hdr8147();
            // format 38
            std::fprintf(f25,"%s%s%s%s%s%s%s\n",
                f_iw(n,6).c_str(),f_iw(itonew,6).c_str(),f_iw(xi1,9).c_str(),
                f_iw(xi2,9).c_str(),f_iw(et1,9).c_str(),f_iw(et2,9).c_str(),
                f_iw(nfb,9).c_str());
            for(int lp=1;lp<=nfb;++lp){
                hdr8149();
                std::fprintf(f25,"      %s%s%s%s%s\n",
                    f_iw(ifrom[lp-1],6).c_str(),f_iw(xif1[lp-1],9).c_str(),
                    f_iw(xif2[lp-1],9).c_str(),f_iw(etf1[lp-1],9).c_str(),
                    f_iw(etf2[lp-1],9).c_str());
                if (ioflag==2){
                    hdr8148();
                    std::fprintf(f25,"   %s%s%s%s%s%s\n",
                        f_fw(DX(n,lp),9,4).c_str(),f_fw(DY(n,lp),9,4).c_str(),
                        f_fw(DZ(n,lp),9,4).c_str(),
                        f_fw(DTX(n,lp)*conversion_.radtodeg,9,4).c_str(),
                        f_fw(DTY(n,lp)*conversion_.radtodeg,9,4).c_str(),
                        f_fw(DTZ(n,lp)*conversion_.radtodeg,9,4).c_str());
                }
            }

            for(int m=1;m<=ncgg[IX(n,nfb+2)-1];++m){
                ntemp=ntemp+1;
                int lprev=ntemp-1;
                if (m==1) ifiner[ntemp-1]=n; else ifiner[ntemp-1]=ntemp-1;
                iifit[ntemp-1] =iifit[n-1];
                iic0[ntemp-1]  =iic0[n-1];
                iiorph[ntemp-1]=iiorph[n-1];
                iitoss[ntemp-1]=iitoss[n-1];
                llimit[ntemp-1]=llimit[n-1];
                iitmax[ntemp-1]=iitmax[n-1];
                mmcxie[ntemp-1]=mmcxie[n-1];
                mmceta[ntemp-1]=mmceta[n-1];
                for(int lp=1;lp<=nfb;++lp){
                    DX(ntemp,lp)=DX(n,lp); DY(ntemp,lp)=DY(n,lp); DZ(ntemp,lp)=DZ(n,lp);
                    DTX(ntemp,lp)=DTX(n,lp); DTY(ntemp,lp)=DTY(n,lp); DTZ(ntemp,lp)=DTZ(n,lp);
                }
                IX(ntemp,1)=nfb;
                IX(ntemp,nfb+2)=nblg[IX(n,nfb+2)-1]+m;
                IX(ntemp,nfb+3)=IX(n,nfb+3);
                IX(ntemp,2*nfb+3)=IX(n,2*nfb+3);
                int nfb1=IX(lprev,1);
                int nbl=IX(ntemp,nfb+2); (void)nbl;

                xi1=IX(n,2*nfb+6); xi2=IX(n,2*nfb+7);
                et1=IX(n,2*nfb+8); et2=IX(n,2*nfb+9);
                IX(ntemp,2*nfb+6)=xi1/(2*m)+1;
                IX(ntemp,2*nfb+7)=xi2/(2*m)+1;
                IX(ntemp,2*nfb+8)=et1/(2*m)+1;
                IX(ntemp,2*nfb+9)=et2/(2*m)+1;
                if (IX(n,2*nfb+6)<=2) IX(ntemp,2*nfb+6)=IX(n,2*nfb+6);
                if (IX(n,2*nfb+7)<=2) IX(ntemp,2*nfb+7)=IX(n,2*nfb+7);
                if (IX(n,2*nfb+8)<=2) IX(ntemp,2*nfb+8)=IX(n,2*nfb+8);
                if (IX(n,2*nfb+9)<=2) IX(ntemp,2*nfb+9)=IX(n,2*nfb+9);
                xi1=IX(ntemp,2*nfb+6); xi2=IX(ntemp,2*nfb+7);
                et1=IX(ntemp,2*nfb+8); et2=IX(ntemp,2*nfb+9);

                for(int lp=1;lp<=nfb;++lp){
                    int xif11=IX(n,2*nfb+9+lp);
                    int xif22=IX(n,3*nfb+9+lp);
                    int etf11=IX(n,4*nfb+9+lp);
                    int etf22=IX(n,5*nfb+9+lp);
                    IX(ntemp,2*nfb+9+lp)=xif11/(2*m)+1;
                    IX(ntemp,3*nfb+9+lp)=xif22/(2*m)+1;
                    IX(ntemp,4*nfb+9+lp)=etf11/(2*m)+1;
                    IX(ntemp,5*nfb+9+lp)=etf22/(2*m)+1;
                    if (IX(n,2*nfb+9+lp)<=2) IX(ntemp,2*nfb+9+lp)=IX(n,2*nfb+9+lp);
                    if (IX(n,3*nfb+9+lp)<=2) IX(ntemp,3*nfb+9+lp)=IX(n,3*nfb+9+lp);
                    if (IX(n,4*nfb+9+lp)<=2) IX(ntemp,4*nfb+9+lp)=IX(n,4*nfb+9+lp);
                    if (IX(n,5*nfb+9+lp)<=2) IX(ntemp,5*nfb+9+lp)=IX(n,5*nfb+9+lp);
                }

                IX(ntemp,2*nfb+4)=(xi2-xi1)*(et2-et1);
                IX(ntemp,2*nfb+5)=IX(lprev,2*nfb1+4)+IX(lprev,2*nfb1+5);
                int mxxchk=IX(ntemp,2*nfb+4)+IX(ntemp,2*nfb+5)-1;
                if (icall>0){
                    if (mxxchk>mxxe){
                        FILE* f11=FUNIT(11);
                        if(f11){ std::string g=filenam_.grdmov.s; g.resize(60,' ');
                            std::fprintf(f11," program terminated in dynamic patching "
                                "routines - see file %s\n",g.c_str()); }
                        std::fprintf(f25,"  stopping...parameter mxxe too small\n");
                        std::fprintf(f25,"  must make  mxxe at least %s\n",
                                     f_iw(mxxchk,12).c_str());
                        termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
                    }
                }

                for(int i=1;i<=nfb;++i){
                    IX(ntemp,i+1)=nblg[IX(n,i+1)-1]+m;
                    IX(ntemp,nfb+i+2)=IX(n,nfb+i+2);
                    if (ncgg[IX(n,i+1)-1]<m)
                        IX(ntemp,i+1)=nblg[IX(n,i+1)-1]+ncgg[IX(n,i+1)-1];
                    ifrom[i-1]=(IX(ntemp,i+1)*100)+IX(ntemp,i+nfb+2);
                }
                ito=(IX(ntemp,nfb+2)*100)+IX(ntemp,2*nfb+3);
                hdr8147();
                // format 37: (2h *,i4,i6,10i9)
                std::fprintf(f25," *%s%s%s%s%s%s%s\n",
                    f_iw(ntemp,4).c_str(),f_iw(ito,6).c_str(),f_iw(xi1,9).c_str(),
                    f_iw(xi2,9).c_str(),f_iw(et1,9).c_str(),f_iw(et2,9).c_str(),
                    f_iw(nfb,9).c_str());
                for(int lp=1;lp<=nfb;++lp){
                    hdr8149();
                    std::fprintf(f25,"      %s%s%s%s%s\n",
                        f_iw(ifrom[lp-1],6).c_str(),
                        f_iw(IX(ntemp,2*nfb+9+lp),9).c_str(),
                        f_iw(IX(ntemp,3*nfb+9+lp),9).c_str(),
                        f_iw(IX(ntemp,4*nfb+9+lp),9).c_str(),
                        f_iw(IX(ntemp,5*nfb+9+lp),9).c_str());
                    if (ioflag==2){
                        hdr8148();
                        std::fprintf(f25,"   %s%s%s%s%s%s\n",
                            f_fw(DX(n,lp),9,4).c_str(),f_fw(DY(n,lp),9,4).c_str(),
                            f_fw(DZ(n,lp),9,4).c_str(),
                            f_fw(DTX(n,lp)*conversion_.radtodeg,9,4).c_str(),
                            f_fw(DTY(n,lp)*conversion_.radtodeg,9,4).c_str(),
                            f_fw(DTZ(n,lp)*conversion_.radtodeg,9,4).c_str());
                    }
                }
            }
        }

        // renumber "to"/"from" block numbers to CFL3D ordering
        IX(n,nfb+2)=nblg[IX(n,nfb+2)-1];
        for(int i=1;i<=nfb;++i)
            IX(n,i+1)=nblg[IX(n,i+1)-1];
    } // end do 8400
    nintr=ntemp;

    // check parameter sizes, multigridability
    int mpa=0;
    int mxxchk=0;
    int nfb1=0;
    for(int n=1;n<=nintr;++n){
        if (n>1) nfb1=IX(n-1,1);
        nfb=IX(n,1);
        int nbl=IX(n,nfb+2);
        xi1=IX(n,2*nfb+6); xi2=IX(n,2*nfb+7);
        et1=IX(n,2*nfb+8); et2=IX(n,2*nfb+9);
        int id2=idimg[nbl-1], jd2=jdimg[nbl-1], kd2=kdimg[nbl-1];
        int itop=IX(n,2*nfb+3)/10;
        if (itop==1) mpa=i_max(i_max(jd2,kd2),mpa);
        if (itop==2) mpa=i_max(i_max(id2,kd2),mpa);
        if (itop==3) mpa=i_max(i_max(id2,jd2),mpa);
        IX(n,2*nfb+4)=(xi2-xi1)*(et2-et1);
        for(int lp=1;lp<=nfb;++lp){
            int mbl=IX(n,lp+1);
            int id1=idimg[mbl-1], jd1=jdimg[mbl-1], kd1=kdimg[mbl-1];
            int itopf=(IX(n,lp+nfb+2)-IX(n,lp+nfb+2)/100*100)/10;
            if (itopf==1) mpa=i_max(i_max(jd1,kd1),mpa);
            if (itopf==2) mpa=i_max(i_max(id1,kd1),mpa);
            if (itopf==3) mpa=i_max(i_max(id1,jd1),mpa);
            params_.lmptch=i_max(params_.lmptch,mpa);
        }
        if (icall>0){
            if (mpa>mptch){
                std::fprintf(f25," mptch is too small; should be: %s\n",
                             f_iw(mpa,4).c_str());
                termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
            }
        }
        if (n==1) IX(n,2*nfb+5)=1;
        else IX(n,2*nfb+5)=IX(n-1,2*nfb1+4)+IX(n-1,2*nfb1+5);
        mxxchk=IX(nintr,2*nfb+4)+IX(nintr,2*nfb+5)-1;
    }

    if (icall>0){
        if (mxxchk>mxxe){
            FILE* f11=FUNIT(11);
            if(f11){ std::string g=filenam_.grdmov.s; g.resize(60,' ');
                std::fprintf(f11," program terminated in dynamic patching "
                    "routines - see file %s\n",g.c_str()); }
            std::fprintf(f25,"  stopping...parameter mxxe too small\n");
            std::fprintf(f25,"  must make mxxe at least %s\n",f_iw(mxxchk,12).c_str());
            termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
        }
    }

    // check multigridability - "to" side only
    int istop=0;
    for(int n=1;n<=nintr;++n){
        nfb=IX(n,1);
        int nbl=IX(n,nfb+2);
        xi1=IX(n,2*nfb+6); xi2=IX(n,2*nfb+7);
        et1=IX(n,2*nfb+8); et2=IX(n,2*nfb+9);
        if (levelg[nbl-1]>1){
            if ((double)(xi1/2)==(double)xi1/2.0){
                std::fprintf(f25," begining \"to\" xie index not multigridable"
                    " for interpolation %s\n",f_iw(n,4).c_str());
                istop=1;
            }
            if ((double)(xi2/2)==(double)xi2/2.0){
                std::fprintf(f25," ending \"to\" xie index not multigridable"
                    " for interpolation %s\n",f_iw(n,4).c_str());
                istop=1;
            }
            if (idimg[nbl-1]==2){
                if (et1!=1 || et2!=2){
                    std::fprintf(f25," for 2D cases must have eta1 = 1 and eta2 = 2\n");
                    istop=1;
                }
            } else {
                if ((double)(et1/2)==(double)et1/2.0){
                    std::fprintf(f25," begining \"to\" eta index not multigridable"
                        " for interpolation %s\n",f_iw(n,4).c_str());
                    istop=1;
                }
                if ((double)(et2/2)==(double)et2/2.0){
                    std::fprintf(f25," ending \"to\" eta index not multigridable"
                        " for interpolation %s\n",f_iw(n,4).c_str());
                    istop=1;
                }
            }
        }
    }
    if (istop>0){
        std::fprintf(f25," stopping\n");
        termn8(myid,ierrflg,ibufdim,nbuf,bou,nou);
    }

    // ichk = 0 (diagnostic dump block skipped)

    // adjust patch parameters so dynamic data can be appended (icall==0)
    if (icall==0){
        params_.lmsub1=nfbmax;
        params_.lnsub1=i_max(params_.lnsub1,params_.lmsub1);
        params_.lintmx=nintr;
        if (params_.lintmax>1) params_.lintmax=params_.lintmax+params_.lintmx;
        else params_.lintmax=params_.lintmx;
        params_.lmxxe=0;
        int ntot=0;
        for(int intv2=1;intv2<=nintr;++intv2){
            int nfbl=IX(intv2,1);
            int npts=IX(intv2,2*nfbl+4);
            int nst =IX(intv2,2*nfbl+5);
            ntot=nst+npts-1;
        }
        params_.lmxxe=ntot;
        if (params_.lmaxxe>1) params_.lmaxxe=params_.lmaxxe+params_.lmxxe;
        else params_.lmaxxe=params_.lmxxe;
    }
}
