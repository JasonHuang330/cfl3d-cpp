// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// sizer.cpp — establish parameter sizes for ronnie (preronnie stage, unit 66).
// Mirrors sizer.F: allocates sizing arrays (tracking memuse via umalloc), reads
// the deck header + grid/patch data, calls global2 (icall=0), prints sizes,
// then rewinds iunit5. Sequential, single-precision build.
#include "ron_common.h"

void sizer(int& maxgr,int& maxbl,int& maxxe,int& nsub1,int& intmax,int& mpatch,
           int& nwork,int nbuf0,int ibufdim0,int maxgr0,int maxbl0,int maxxe0,
           int nsub10,int intmax0,int mpatch0,int imode){
    (void)maxxe0;(void)mpatch0;
    FILE* f5=FUNIT(unit5_.iunit5);

    int ierrflg=-99;
    long memuse=0;

    // ---- allocations (sizes tracked; contents used below) ----
    auto AL_D=[&](long words,const char* nm)->std::vector<double>{
        umalloc(words,0,nm,memuse,0); return std::vector<double>((size_t)words,0.0); };
    auto AL_I=[&](long words,const char* nm)->std::vector<int>{
        umalloc(words,1,nm,memuse,0); return std::vector<int>((size_t)words,0); };

    auto dthetx=AL_D((long)intmax0*nsub10,"dthetx");
    auto dthety=AL_D((long)intmax0*nsub10,"dthety");
    auto dthetz=AL_D((long)intmax0*nsub10,"dthetz");
    auto dx    =AL_D((long)intmax0*nsub10,"dx");
    auto dy    =AL_D((long)intmax0*nsub10,"dy");
    auto dz    =AL_D((long)intmax0*nsub10,"dz");
    auto etf1  =AL_I(nsub10,"etf1");
    auto etf2  =AL_I(nsub10,"etf2");
    auto idimg =AL_I(maxbl0,"idimg");
    auto iemg  =AL_I(maxgr0,"iemg");
    auto ifiner=AL_I(intmax0,"ifiner");
    auto ifrom =AL_D(nsub10,"ifrom");      // note: Fortran allocs ifrom as REAL(intflag 0)
    std::vector<int> ifrom_i((size_t)nsub10,0); // integer working copy for global2
    auto igridg=AL_I(maxbl0,"igridg");
    auto iic0  =AL_I(intmax0,"iic0");
    auto iifit =AL_I(intmax0,"iifit");
    auto iiint1=AL_I(nsub10,"iiint1");
    auto iiint2=AL_I(nsub10,"iiint2");
    auto iindex=AL_I((long)intmax0*(6*nsub10+9),"iindex");
    auto iiorph=AL_I(intmax0,"iiorph");
    auto iitmax=AL_I(intmax0,"iitmax");
    auto iitoss=AL_I(intmax0,"iitoss");
    auto isav_pat  =AL_I((long)intmax0*17,"isav_pat");
    auto isav_pat_b=AL_I((long)intmax0*nsub10*6,"isav_pat_b");
    auto jdimg =AL_I(maxbl0,"jdimg");
    auto jjmax1=AL_I(nsub10,"jjmax1");
    auto kdimg =AL_I(maxbl0,"kdimg");
    auto kkmax1=AL_I(nsub10,"kkmax1");
    auto levelg=AL_I(maxbl0,"levelg");
    auto llimit=AL_I(intmax0,"llimit");
    auto mblk2nd=AL_I(maxbl0,"mblk2nd");
    auto mglevg=AL_I(maxbl0,"mglevg");
    auto mmceta=AL_I(intmax0,"mmceta");
    auto mmcxie=AL_I(intmax0,"mmcxie");
    auto nblcg =AL_I(maxbl0,"nblcg");
    auto nblg  =AL_I(maxgr0,"nblg");
    auto ncgg  =AL_I(maxgr0,"ncgg");
    auto ncheck=AL_I(maxbl0,"ncheck");
    auto nemgl =AL_I(maxbl0,"nemgl");
    auto xif1  =AL_I(nsub10,"xif1");
    auto xif2  =AL_I(nsub10,"xif2");

    // open unit 66 (preronnie.out)
    FUNIT(66)=std::fopen("preronnie.out","w");
    FILE* f66=FUNIT(66);

    tc::TextReader tr(f5);
    tr.line();                              // read(iunit5,*)
    auto rd_a60=[&](FStr& dst){ std::string s=tr.line();
        if((int)s.size()>60) s=s.substr(0,60); dst=FStr(80); dst.assign(s); };
    rd_a60(filenam_.grid);                  // read(iunit5,'(a60)')grid
    rd_a60(filenam_.output);                // read(iunit5,'(a60)')output
    rd_a60(filenam_.patch);                 // read(iunit5,'(a60)')patch

    double pi=4.0*std::atan(1.0);
    conversion_.radtodeg=180.0/pi;

    mydist2_.nnodes=1; mydist2_.myhost=0; mydist2_.myid=0; mydist2_.mycomm=0;
    for(int nn=1;nn<=maxbl0;++nn) mblk2nd[nn-1]=mydist2_.myhost;

    // initialize output buffers
    std::vector<int> nou((size_t)nbuf0,0);
    std::vector<FStr> bou((size_t)ibufdim0*nbuf0,FStr(120));

    // default dimensions
    params_.lmaxgr=1; params_.lmaxbl=1; params_.lmxseg=1; params_.lmaxcs=1;
    params_.lnplts=1; params_.lmxbli=1; params_.lmaxxe=1; params_.lnsub1=1;
    params_.lintmx=1; params_.lmxxe=1; params_.liitot=1; params_.isum=1;
    params_.lncycm=1; params_.isum_n=1; params_.lminnode=1; params_.isumi=1;
    params_.isumi_n=1; params_.lmptch=1; params_.lmsub1=1; params_.lintmax=1;
    params_.libufdim=1; params_.lnbuf=1; params_.llbcprd=1; params_.llbcemb=1;
    params_.llbcrad=1;

    // banner + memory line (single precision) to unit 66
    ron_banner(66,0,(double)memuse/1.0e6);
    ron_iofiles(66,filenam_.grid,filenam_.output,filenam_.patch);

    int ioflag=0,itrace=0;
    tr.line();                              // read(iunit5,*)
    { int v[2]; tr.ints(v,2); ioflag=v[0]; itrace=v[1]; }
    // format 398: (/,4x,6hioflag,4x,6hitrace/,2i10)
    std::fprintf(f66,"\n    ioflag    itrace\n%s%s\n",
                 f_iw(ioflag,10).c_str(),f_iw(itrace,10).c_str());

    // read title (20a4) -> 80 chars
    {
        std::string t=tr.line();
        // title stored as 20 real*8 words in Fortran; we keep it as text for
        // output only. Print via FORMAT 11: (2h  ,20a4) = "  " + 80-char field.
        std::string t80=t; t80.resize(80,' ');
        std::fprintf(f66,"\ntitle\n  %s\n",t80.c_str());
    }

    int ngrid=0;
    tr.line();                              // read(iunit5,*)
    { int v[1]; tr.ints(v,1); ngrid=v[0]; }
    std::fprintf(f66,"\ngrid/level data\n");
    std::fprintf(f66," ngrid\n");
    std::fprintf(f66,"%s\n",f_iw(ngrid,6).c_str());  // format 36 (13i6)

    igrdtyp_.ip3dgrd=0;
    if (ngrid<0){ igrdtyp_.ip3dgrd=1; ngrid=f_iabs(ngrid); }
    int nchk=maxgr0-ngrid;
    if (nchk<0){
        std::fprintf(f66," stopping - insufficient maximum number of grids(maxgr)\n");
        termn8(mydist2_.myid,ierrflg,ibufdim0,nbuf0,bou.data(),nou.data());
    }

    tr.line();                              // read(iunit5,*)
    int nbl=0, iemtot=0;
    std::fprintf(f66,"   ncg   iem  idim  jdim  kdim\n");  // format 1631
    int ncgmax=0, iemmax=0;
    for(int igrid=1;igrid<=ngrid;++igrid){
        nbl=nbl+1;
        int ncg,iem,idim,jdim,kdim;
        { int v[5]; tr.ints(v,5); ncg=v[0]; iem=v[1]; idim=v[2]; jdim=v[3]; kdim=v[4]; }
        // format 36 (13i6)
        std::fprintf(f66,"%s%s%s%s%s\n",f_iw(ncg,6).c_str(),f_iw(iem,6).c_str(),
                     f_iw(idim,6).c_str(),f_iw(jdim,6).c_str(),f_iw(kdim,6).c_str());
        iemtot+=iem;
        ncgg[igrid-1]=ncg;
        if (igrid==1){ ncgmax=ncg; iemmax=iem; }
        else { ncgmax=i_max(ncgmax,ncg); iemmax=i_max(iemmax,iem); }
        iemg[igrid-1]=iem;
        nblg[igrid-1]=nbl;
        idimg[nbl-1]=idim; jdimg[nbl-1]=jdim; kdimg[nbl-1]=kdim;
        if (ncg>0){
            if (iem>0){
                std::fprintf(f66," embedded grids must have ncg = 0\n");
                termn8(mydist2_.myid,ierrflg,ibufdim0,nbuf0,bou.data(),nou.data());
            }
            for(int nn=1;nn<=ncg;++nn){
                nbl=nbl+1;
                idimg[nbl-1]=idimg[nbl-2]/2+1;
                jdimg[nbl-1]=jdimg[nbl-2]/2+1;
                kdimg[nbl-1]=kdimg[nbl-2]/2+1;
                if (idimg[nbl-2]<=2) idimg[nbl-1]=idimg[nbl-2];
                int istop=0;
                if ((double)(idimg[nbl-2]/2)==(double)idimg[nbl-2]/2.0 && idim>2){
                    std::fprintf(f66," Cannot create coarser level for idim past%s\n",
                                 f_iw(idimg[nbl-2],6).c_str()); istop=1;
                }
                if ((double)(jdimg[nbl-2]/2)==(double)jdimg[nbl-2]/2.0){
                    std::fprintf(f66," Cannot create coarser level for jdim past%s\n",
                                 f_iw(jdimg[nbl-2],6).c_str()); istop=1;
                }
                if ((double)(kdimg[nbl-2]/2)==(double)kdimg[nbl-2]/2.0){
                    std::fprintf(f66," Cannot create coarser level for kdim past%s\n",
                                 f_iw(kdimg[nbl-2],6).c_str()); istop=1;
                }
                if (istop==1)
                    termn8(mydist2_.myid,ierrflg,ibufdim0,nbuf0,bou.data(),nou.data());
            }
        }
    }
    (void)iemtot;

    info_.mseq=1;
    if (ncgmax>0) info_.mgflag=1;
    if (iemmax>0) info_.mgflag=2;
    if (info_.mseq>1){
        nbl=0;
        for(int igrid=1;igrid<=ngrid;++igrid){
            nbl=nbl+1; int ncg=ncgg[igrid-1]; if(ncg>0) nbl=nbl+ncg;
        }
    }
    for(int m=1;m<=info_.mseq;++m){ mglevg[m-1]=ncgmax+1; nemgl[m-1]=iemmax; }

    int nblock=nbl;
    nchk=maxbl0-nblock;
    if (nchk<0){
        std::fprintf(f66," maxbl,nchk=%s%s\n",f_iw(maxbl0,5).c_str(),f_iw(nchk,5).c_str());
        termn8(mydist2_.myid,ierrflg,ibufdim0,nbuf0,bou.data(),nou.data());
    }

    for(int m=1;m<=info_.mseq;++m){
        info_.levelt[m-1]=ncgmax-(info_.mseq-m)+nemgl[m-1]+1;
        info_.levelb[m-1]=info_.levelt[m-1]-(mglevg[m-1]-1)-nemgl[m-1];
        if (info_.levelb[m-1]<1){
            std::fprintf(f66," error in input, m, levelt, levelb, ncgmax=%s%s%s%s\n",
                f_iw(m,5).c_str(),f_iw(info_.levelt[m-1],5).c_str(),
                f_iw(info_.levelb[m-1],5).c_str(),f_iw(ncgmax,5).c_str());
            termn8(mydist2_.myid,ierrflg,ibufdim0,nbuf0,bou.data(),nou.data());
        }
    }

    params_.lmaxgr=ngrid;
    params_.lmaxbl=nblock;

    int icall1=0, iunit1=66, imode1=0, ninter0=0;
    // copy ifrom(real) placeholder not needed; global2 uses integer ifrom
    global2(maxbl0,maxgr0,nsub10,ninter0,intmax0,ngrid,idimg.data(),
            jdimg.data(),kdimg.data(),levelg.data(),ncgg.data(),nblg.data(),
            iindex.data(),llimit.data(),iitmax.data(),mmcxie.data(),mmceta.data(),
            ncheck.data(),iifit.data(),iic0.data(),iiorph.data(),iitoss.data(),
            ifiner.data(),dx.data(),dy.data(),dz.data(),dthetx.data(),
            dthety.data(),dthetz.data(),mydist2_.myid,mpatch0,maxxe0,icall1,iunit1,
            nou.data(),bou.data(),ibufdim0,nbuf0,ifrom_i.data(),xif1.data(),
            etf1.data(),xif2.data(),etf2.data(),igridg.data(),iemg.data(),
            nblock,ioflag,imode1);

    // work array sizing
    int nstart=1;
    for(int nb=1;nb<=nblock;++nb){
        int j=jdimg[nb-1], k=kdimg[nb-1], i=idimg[nb-1];
        int ns=j*k*i;
        int lw=nstart; lw+=ns; lw+=ns; nstart=lw+ns;
    }
    nwork=nstart;
    nwork=nwork+1;
    maxgr=params_.lmaxgr;
    maxbl=params_.lmaxbl;
    maxxe=params_.lmaxxe;
    nsub1=params_.lnsub1;
    mpatch=params_.lmptch;
    intmax=params_.lintmax;
    int nbuf_out=nbuf0, ibufdim_out=ibufdim0; (void)nbuf_out;(void)ibufdim_out;

    // list-directed writes: leading blank
    std::fprintf(f66," \n");
    std::fprintf(f66," ***********************************************\n");
    std::fprintf(f66," \n");
    std::fprintf(f66,"     PARAMETER SIZES REQUIRED FOR THIS CASE:\n");
    std::fprintf(f66," \n");
    std::fprintf(f66," ***********************************************\n");
    std::fprintf(f66," \n");
    std::fprintf(f66," nwork   = %s\n",f_iw(nwork,10).c_str());
    std::fprintf(f66," intmax  = %s\n",f_iw(intmax,10).c_str());
    std::fprintf(f66," maxxe   = %s\n",f_iw(maxxe,10).c_str());
    std::fprintf(f66," nsub1   = %s\n",f_iw(nsub1,10).c_str());
    std::fprintf(f66," mpatch  = %s\n",f_iw(mpatch,10).c_str());
    std::fprintf(f66," maxbl   = %s\n",f_iw(maxbl,10).c_str());
    std::fprintf(f66," maxgr   = %s\n",f_iw(maxgr,10).c_str());
    std::fprintf(f66," nbuf    = %s\n",f_iw(nbuf0,10).c_str());
    std::fprintf(f66," ibufdim = %s\n",f_iw(ibufdim0,10).c_str());

    nstart=1;
    for(int nb=1;nb<=nblock;++nb){
        int j=jdimg[nb-1], k=kdimg[nb-1], i=idimg[nb-1];
        int ns=j*k*i;
        int lw=nstart; lw+=ns; lw+=ns; nstart=lw+ns;
    }

    std::fprintf(f66,"  \n");
    std::fprintf(f66," *******************************************************\n");
    std::fprintf(f66,"  \n");
    std::fprintf(f66,"             SUMMARY OF STORAGE REQUIREMENTS \n");
    std::fprintf(f66,"  \n");
    // format 790: ('      permanent array w  requires ',i9,' (words)')
    std::fprintf(f66,"      permanent array w  requires %s (words)\n",
                 f_iw((int)nstart,9).c_str());
    std::fprintf(f66,"  \n");
    std::fprintf(f66,"  \n");
    // format 800: ('          >>> Estimate for nwork = ',i9,' <<<')
    std::fprintf(f66,"          >>> Estimate for nwork = %s <<<\n",
                 f_iw(nstart,9).c_str());
    std::fprintf(f66,"  \n");
    std::fprintf(f66," *******************************************************\n");
    std::fprintf(f66,"  \n");

    std::rewind(f5);   // rewind (iunit5)

    // free memory: (no-op — vectors auto-free)

    if (imode==1){
        std::fprintf(f66,"\n memory for preronnie has been deallocated\n");
    } else {
        std::fprintf(stdout,"\npreronnie has completed successfully\n");
        std::fprintf(stdout,"preronnie information has been put in"
                            " file preronnie.out\n\n");
    }
    if(f66){ std::fclose(f66); FUNIT(66)=nullptr; }
}
