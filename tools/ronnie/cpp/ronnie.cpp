// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// ronnie.cpp — C++ translation of ronnie/ronnie.F (subroutine ronnie).
// Top-level orchestrator for the CFL3D generalized patched-grid preprocessor.
//
// Reads the input deck (unit iunit5=55), reads grids, sets up 1-D work indices,
// loops over patch interfaces calling patcher (which drives topol/invert/shear/
// project/global2/diagnos/...), and writes the connectivity file patch.bin
// (unit 2) plus the primary information file ronnie.out (unit 25).
//
// COMMON blocks -> shared structs in ron_common.h.
// Arrays are column-major, 1-based; allocated here as std::vector and passed as
// raw pointers to the workers. Every allocate/umalloc pair is reproduced so the
// reported "memory allocation" (memuse) matches the reference byte-for-byte.
#include "ron_common.h"

void ronnie(int maxbl,int maxgr,int nwork,int intmax,int nsub1,int maxxe,
            int mpatch,int nbuf,int ibufdim){
    int ierrflg = -1;

    FILE* f25 = FUNIT(25);
    FILE* f5  = FUNIT(unit5_.iunit5);
    tc::TextReader tr(f5);

    // ---- memory allocation (track memuse; allocate as vectors) -------------
    long memuse = 0;
    auto AL_D=[&](long words,const char* nm)->std::vector<double>{
        umalloc(words,0,nm,memuse,0); return std::vector<double>((size_t)words,0.0); };
    auto AL_I=[&](long words,const char* nm)->std::vector<int>{
        umalloc(words,1,nm,memuse,0); return std::vector<int>((size_t)words,0); };

    auto dthetx = AL_D((long)intmax*nsub1,"dthetx");
    auto dthety = AL_D((long)intmax*nsub1,"dthety");
    auto dthetz = AL_D((long)intmax*nsub1,"dthetz");
    auto dx     = AL_D((long)intmax*nsub1,"dx");
    auto dy     = AL_D((long)intmax*nsub1,"dy");
    auto dz     = AL_D((long)intmax*nsub1,"dx");
    auto eta2s  = AL_D((long)(mpatch+2)*(mpatch+2),"eta2s");
    auto etf1   = AL_I(nsub1,"etf1");
    auto etf2   = AL_I(nsub1,"etf2");
    auto factjhi= AL_D((long)intmax*nsub1,"factjhi");
    auto factjlo= AL_D((long)intmax*nsub1,"factjlo");
    auto factkhi= AL_D((long)intmax*nsub1,"factkhi");
    auto factklo= AL_D((long)intmax*nsub1,"factklo");
    auto idimg  = AL_I(maxbl,"idimg");
    auto iemg   = AL_I(maxgr,"iemg");
    auto ifiner = AL_I(intmax,"ifiner");
    auto ifrom  = AL_D(nsub1,"ifrom");
    std::vector<int> ifrom_i((size_t)nsub1,0);
    auto igridg = AL_I(maxbl,"igridg");
    auto iic0   = AL_I(intmax,"iic0");
    auto iifit  = AL_I(intmax,"iifit");
    auto iiint1 = AL_I(nsub1,"iiint1");
    auto iiint2 = AL_I(nsub1,"iiint2");
    auto iindex = AL_I((long)intmax*(6*nsub1+9),"iindex");
    auto iiorph = AL_I(intmax,"iiorph");
    auto iitmax = AL_I(intmax,"iitmax");
    auto iitoss = AL_I(intmax,"iitoss");
    auto ireq_ar= AL_I((long)intmax*3,"ireq_ar");
    auto isav_pat  = AL_I((long)intmax*17,"isav_pat");
    auto isav_pat_b= AL_I((long)intmax*nsub1*6,"isav_pat_b");
    auto itest  = AL_I(maxgr,"itest");
    auto jdimg  = AL_I(maxbl,"jdimg");
    auto jimage = AL_I((long)nsub1*(mpatch+2)*(mpatch+2),"jimage");
    auto jjmax1 = AL_I(nsub1,"jjmax1");
    auto jmm    = AL_I(mpatch+2,"jmm");
    auto jte    = AL_I(nsub1,"jte");
    auto jtest  = AL_I(maxgr,"jtest");
    auto kdimg  = AL_I(maxbl,"kdimg");
    auto kimage = AL_I((long)nsub1*(mpatch+2)*(mpatch+2),"kimage");
    auto kkmax1 = AL_I(nsub1,"kkmax1");
    auto kmm    = AL_I(mpatch+2,"kmm");
    auto kte    = AL_I(nsub1,"kte");
    auto ktest  = AL_I(maxgr,"ktest");
    auto levelg = AL_I(maxbl,"levelg");
    auto llimit = AL_I(intmax,"llimit");
    auto lout   = AL_I(nsub1,"lout");
    auto lw     = AL_I((long)65*maxbl,"lw");
    auto lw2    = AL_I((long)43*maxbl,"lw2");
    auto mblk2nd= AL_I(maxbl,"mblk2nd");
    auto mglevg = AL_I(maxbl,"mglevg");
    auto mmceta = AL_I(intmax,"mmceta");
    auto mmcxie = AL_I(intmax,"mmcxie");
    auto nblcg  = AL_I(maxbl,"nblcg");
    auto nblg   = AL_I(maxgr,"nblg");
    auto nblk1  = AL_I(mpatch+2,"nblk1");
    auto nblk2  = AL_I(mpatch+2,"nblk2");
    auto nblkpt = AL_I(maxxe,"nblkpt");
    auto ncgg   = AL_I(maxgr,"ncgg");
    auto ncheck = AL_I(maxbl,"ncheck");
    auto nemgl  = AL_I(maxbl,"nemgl");
    auto seta   = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"seta");
    auto seta2  = AL_D((long)(mpatch+2)*(mpatch+2),"seta2");
    auto sxie   = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"sxie");
    auto sxie2  = AL_D((long)(mpatch+2)*(mpatch+2),"sxie2");
    auto temp   = AL_D((long)(mpatch+2)*(mpatch+2),"temp");
    auto windex = AL_D((long)maxxe*2,"windex");
    auto work   = AL_D(nwork,"work");
    auto x1     = AL_D((long)(mpatch+2)*(mpatch+2),"x1");
    auto x2     = AL_D((long)(mpatch+2)*(mpatch+2),"x2");
    auto xie2s  = AL_D((long)(mpatch+2)*(mpatch+2),"xie2s");
    auto xif1   = AL_I(nsub1,"xif1");
    auto xif2   = AL_I(nsub1,"xif2");
    auto xmi    = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"xmi");
    auto xmie   = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"xmie");
    auto xorig  = AL_D(maxbl," xorig");
    auto xte    = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"xte");
    auto y1     = AL_D((long)(mpatch+2)*(mpatch+2),"y1");
    auto y2     = AL_D((long)(mpatch+2)*(mpatch+2),"y2");
    auto ymi    = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"ymi");
    auto ymie   = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"ymie");
    auto yorig  = AL_D(maxbl,"yorig");
    auto yte    = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"yte");
    auto z1     = AL_D((long)(mpatch+2)*(mpatch+2),"z1");
    auto z2     = AL_D((long)(mpatch+2)*(mpatch+2),"z2");
    auto zmi    = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"zmi");
    auto zmie   = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"zmie");
    auto zorig  = AL_D(maxbl,"zorig");
    auto zte    = AL_D((long)(mpatch+2)*(mpatch+2)*nsub1,"zte");

    FStr string(50); string.assign(" ");
    cputim(0,1,string,mydist2_.myhost,mydist2_.myid,mydist2_.mycomm,25);

    // ---- determine machine zero (10**(-iexp)) ------------------------------
    {
        double compare = 1.0;
        int iexp = 0;
        for (int i=1;i<=20;++i){
            double add=1.0;
            for (int n=1;n<=i;++n) add*=0.1;
            double x11 = compare + add;
            if (x11==compare){ iexp=i-1; break; }
        }
        zero_.iexp = iexp;
    }

    double pi = 4.0*std::atan(1.0);
    conversion_.radtodeg = 180.0/pi;

    // dummy parallel-related values
    mydist2_.nnodes=1; mydist2_.myhost=0; mydist2_.myid=0; mydist2_.mycomm=0;
    for (int nn=1;nn<=maxbl;++nn) mblk2nd[nn-1]=mydist2_.myhost;

    // initialize output buffers
    std::vector<int>  nou((size_t)nbuf,0);
    std::vector<FStr> bou((size_t)ibufdim*nbuf,FStr(120));

    // ---- banner + memory line ---------------------------------------------
    ron_banner(25,1,(double)memuse/1.0e6);

    // input/output files
    tr.line();                                   // read(iunit5,*)
    auto rd_a60=[&](FStr& dst){ std::string s=tr.line();
        if((int)s.size()>60) s=s.substr(0,60); dst=FStr(80); dst.assign(s); };
    rd_a60(filenam_.grid);
    rd_a60(filenam_.output);
    rd_a60(filenam_.patch);
    ron_iofiles(25,filenam_.grid,filenam_.output,filenam_.patch);

    int ioflag=0,itrace=0;
    tr.line();                                   // read(iunit5,*)
    { int v[2]; tr.ints(v,2); ioflag=v[0]; itrace=v[1]; }
    tracer_.itrace = itrace;
    // format 398: (/,6hioflag,4x,6hitrace/,i6,4x,i6)
    std::fprintf(f25,"\nioflag    itrace\n%s    %s\n",
                 f_iw(ioflag,6).c_str(),f_iw(itrace,6).c_str());

    // title (20a4) -> keep as text; store packed into info_.title too
    std::string title_txt = tr.line();
    { std::string t80=title_txt; t80.resize(80,' ');
      // pack into info_.title as 20 4-char words (little-endian order of chars
      // within each REAL word is unobservable except via output, which uses the
      // text; store bytes for completeness).
      for (int i=0;i<20;++i){
          char w[8]={' ',' ',' ',' ',0,0,0,0};
          for (int c=0;c<4;++c) w[c]=t80[i*4+c];
          double d; std::memcpy(&d,w,8); info_.title[i]=d;
      }
    }
    // format 111 then 11: (/5htitle) then (2h  ,20a4)
    { std::string t80=title_txt; t80.resize(80,' ');
      std::fprintf(f25,"\ntitle\n  %s\n",t80.c_str()); }

    tr.line();                                   // read(iunit5,*)
    int ngrid=0;
    { int v[1]; tr.ints(v,1); ngrid=v[0]; }
    std::fprintf(f25,"\ngrid/level data\n");
    std::fprintf(f25," ngrid\n");
    std::fprintf(f25,"%s\n",f_iw(ngrid,6).c_str());   // format 36 (13i6)

    igrdtyp_.ip3dgrd=0;
    if (ngrid<0){ igrdtyp_.ip3dgrd=1; ngrid=f_iabs(ngrid); }
    int nchk=maxgr-ngrid;
    if (nchk<0){
        std::fprintf(f25," stopping - insufficient maximum number of grids(maxgr)\n");
        termn8(mydist2_.myid,ierrflg,ibufdim,nbuf,bou.data(),nou.data());
    }

    tr.line();                                   // read(iunit5,*)
    int nbl=0, iemtot=0;
    std::fprintf(f25,"   ncg   iem  idim  jdim  kdim\n");  // format 1631
    int ncgmax=0, iemmax=0;
    int idim=0,jdim=0,kdim=0;
    for (int igrid=1;igrid<=ngrid;++igrid){
        nbl=nbl+1;
        int ncg,iem;
        { int v[5]; tr.ints(v,5); ncg=v[0]; iem=v[1]; idim=v[2]; jdim=v[3]; kdim=v[4]; }
        std::fprintf(f25,"%s%s%s%s%s\n",f_iw(ncg,6).c_str(),f_iw(iem,6).c_str(),
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
                std::fprintf(f25," embedded grids must have ncg = 0\n");
                termn8(mydist2_.myid,ierrflg,ibufdim,nbuf,bou.data(),nou.data());
            }
            for (int n=1;n<=ncg;++n){
                nbl=nbl+1;
                idimg[nbl-1]=idimg[nbl-2]/2+1;
                jdimg[nbl-1]=jdimg[nbl-2]/2+1;
                kdimg[nbl-1]=kdimg[nbl-2]/2+1;
                if (idimg[nbl-2]<=2) idimg[nbl-1]=idimg[nbl-2];
                int istop=0;
                if ((double)(idimg[nbl-2]/2)==(double)idimg[nbl-2]/2.0 && idim>2){
                    std::fprintf(f25," Cannot create coarser level for idim past%s\n",
                                 f_iw(idimg[nbl-2],6).c_str()); istop=1;
                }
                if ((double)(jdimg[nbl-2]/2)==(double)jdimg[nbl-2]/2.0){
                    std::fprintf(f25," Cannot create coarser level for jdim past%s\n",
                                 f_iw(jdimg[nbl-2],6).c_str()); istop=1;
                }
                if ((double)(kdimg[nbl-2]/2)==(double)kdimg[nbl-2]/2.0){
                    std::fprintf(f25," Cannot create coarser level for kdim past%s\n",
                                 f_iw(kdimg[nbl-2],6).c_str()); istop=1;
                }
                if (istop==1)
                    termn8(mydist2_.myid,ierrflg,ibufdim,nbuf,bou.data(),nou.data());
            }
        }
    }
    (void)iemtot;

    info_.mseq=1;
    if (ncgmax>0) info_.mgflag=1;
    if (iemmax>0) info_.mgflag=2;
    if (info_.mseq>1){
        nbl=0;
        for (int igrid=1;igrid<=ngrid;++igrid){
            nbl=nbl+1; int ncg=ncgg[igrid-1]; if(ncg>0) nbl=nbl+ncg;
        }
    }
    for (int m=1;m<=info_.mseq;++m){ mglevg[m-1]=ncgmax+1; nemgl[m-1]=iemmax; }

    int nblock=nbl;
    nchk=maxbl-nblock;
    if (nchk<0){
        std::fprintf(f25," maxbl,nchk=%s%s\n",f_iw(maxbl,5).c_str(),f_iw(nchk,5).c_str());
        termn8(mydist2_.myid,ierrflg,ibufdim,nbuf,bou.data(),nou.data());
    }

    for (int m=1;m<=info_.mseq;++m){
        info_.levelt[m-1]=ncgmax-(info_.mseq-m)+nemgl[m-1]+1;
        info_.levelb[m-1]=info_.levelt[m-1]-(mglevg[m-1]-1)-nemgl[m-1];
        if (info_.levelb[m-1]<1){
            std::fprintf(f25," error in input, m, levelt, levelb, ncgmax=%s%s%s%s\n",
                f_iw(m,5).c_str(),f_iw(info_.levelt[m-1],5).c_str(),
                f_iw(info_.levelb[m-1],5).c_str(),f_iw(ncgmax,5).c_str());
            termn8(mydist2_.myid,ierrflg,ibufdim,nbuf,bou.data(),nou.data());
        }
    }

    // ---- read patch input + summaries (global2, icall=1, unit 25) ----------
    int ninter=0;
    {
        int icall=1, iunit=25, imode=0;
        global2(maxbl,maxgr,nsub1,ninter,intmax,ngrid,idimg.data(),
                jdimg.data(),kdimg.data(),levelg.data(),ncgg.data(),nblg.data(),
                iindex.data(),llimit.data(),iitmax.data(),mmcxie.data(),mmceta.data(),
                ncheck.data(),iifit.data(),iic0.data(),iiorph.data(),iitoss.data(),
                ifiner.data(),dx.data(),dy.data(),dz.data(),dthetx.data(),
                dthety.data(),dthetz.data(),mydist2_.myid,mpatch,maxxe,icall,iunit,
                nou.data(),bou.data(),ibufdim,nbuf,ifrom_i.data(),xif1.data(),
                etf1.data(),xif2.data(),etf2.data(),igridg.data(),iemg.data(),
                nblock,ioflag,imode);
    }

    // ---- set up 1-D storage indices in work array (lw/lw2) -----------------
    // lw(i,nbl) 1-based: element (i,nbl) -> lw[(i-1)+(nbl-1)*65]
    #define LW(i,nb)  lw [((i)-1)+((nb)-1)*65]
    #define LW2(i,nb) lw2[((i)-1)+((nb)-1)*43]

    int nstart=1;
    for (int nb=1;nb<=nblock;++nb){
        int j=jdimg[nb-1], k=kdimg[nb-1], i=idimg[nb-1];
        int igrid=igridg[nb-1]; (void)igrid;
        for (int ll=1;ll<=9;++ll) LW(ll,nb)=nstart;
        LW(10,nb)=nstart;
        int ns=j*k*i;
        LW(11,nb)=LW(10,nb)+ns;
        ns=j*k*i;
        LW(12,nb)=LW(11,nb)+ns;
        ns=j*k*i;
        nstart=LW(12,nb)+ns;
        for (int ll=13;ll<=65;++ll) LW(ll,nb)=nstart;

        ginfo_.jdim=jdimg[nb-1]; ginfo_.kdim=kdimg[nb-1]; ginfo_.idim=idimg[nb-1];
        LW2(1,nb)=jdimg[nb-1];
        LW2(2,nb)=kdimg[nb-1];
        LW2(3,nb)=idimg[nb-1];
        LW2(4,nb)=nb;
        LW2(5,nb)=ginfo_.jdim/2+1;
        LW2(6,nb)=ginfo_.kdim/2+1;
        LW2(7,nb)=ginfo_.idim/2+1;
        for (int ll=8;ll<=40;++ll) LW2(ll,nb)=0;
    }

    if (nstart>nwork){
        std::fprintf(f25," must increase parameter nwork to %d\n",nstart);
        termn8(mydist2_.myid,ierrflg,ibufdim,nbuf,bou.data(),nou.data());
    }

    // read grids and preliminary setup
    setup(lw.data(),lw2.data(),work.data(),nstart,itest.data(),jtest.data(),
          ktest.data(),maxbl,maxgr,ngrid,ncgg.data(),iemg.data(),nblg.data());

    // ---- auxiliary arrays for patching subroutines -------------------------
    ierrflg=-1;
    int icount_pat=0;
    int levt=info_.levelt[info_.mseq-1];
    std::vector<int> is_pat(6,0), ie_pat(6,0);   // 1..levt (Fortran ie_pat(5))
    for (int levl=1;levl<=levt;++levl){
        is_pat[levl]=icount_pat+1;
        for (int nb=1;nb<=nblock;++nb){
            if (levl!=levelg[nb-1]) continue;
            int icount_pat1=icount_pat;
            pre_patch(nb,lw.data(),icount_pat,ninter,iindex.data(),intmax,nsub1,
                      isav_pat.data(),isav_pat_b.data(),jjmax1.data(),kkmax1.data(),
                      iiint1.data(),iiint2.data(),maxbl,jdimg.data(),kdimg.data(),
                      idimg.data(),ierrflg);
            if (icount_pat>icount_pat1) ie_pat[levl]=icount_pat;
        }
    }
    is_dpatch_.maxdcnt=icount_pat;

    // ---- generalized-coordinate interpolation ------------------------------
    std::fprintf(f25,"\n\n BEGINNING GENERALIZED-COORDINATE INTERPOLATION\n");  // format 91

    int it_thro=0;
    int ncall=1;
    int mgwk=nstart;
    int ioutpt=1;
    int mwork=nwork-nstart;

    for (int nb=1;nb<=nblock;++nb){
        lead(nb,lw.data(),lw2.data(),maxbl);
        it_thro=it_thro+1;
        patcher(nb,lw.data(),work.data(),mgwk,&work[mgwk],mwork,ncall,ioutpt,
                it_thro,maxbl,nsub1,intmax,maxxe,mpatch,jdimg.data(),kdimg.data(),
                idimg.data(),windex.data(),ninter,iindex.data(),llimit.data(),
                iitmax.data(),mmcxie.data(),mmceta.data(),ncheck.data(),iifit.data(),
                nblkpt.data(),iic0.data(),iiorph.data(),iitoss.data(),ifiner.data(),
                factjlo.data(),factjhi.data(),factklo.data(),factkhi.data(),
                dx.data(),dy.data(),dz.data(),dthetx.data(),dthety.data(),
                dthetz.data(),isav_pat.data(),isav_pat_b.data(),
                xte.data(),yte.data(),zte.data(),xmi.data(),ymi.data(),zmi.data(),
                xmie.data(),ymie.data(),zmie.data(),jjmax1.data(),kkmax1.data(),
                jimage.data(),kimage.data(),xorig.data(),yorig.data(),zorig.data(),
                jte.data(),kte.data(),sxie.data(),seta.data(),sxie2.data(),
                seta2.data(),xie2s.data(),eta2s.data(),temp.data(),x2.data(),
                y2.data(),z2.data(),nblk1.data(),nblk2.data(),jmm.data(),kmm.data(),
                x1.data(),y1.data(),z1.data(),lout.data(),xif1.data(),xif2.data(),
                etf1.data(),etf2.data(),ireq_ar.data(),
                mydist2_.myid,mydist2_.myhost,mydist2_.mycomm,mblk2nd.data(),
                nou.data(),bou.data(),nbuf,ibufdim);
        writ_buf(nb,25,nou.data(),bou.data(),nbuf,ibufdim,mydist2_.myhost,
                 mydist2_.myid,mydist2_.mycomm,mblk2nd.data(),maxbl);
        writ_buf(nb,9,nou.data(),bou.data(),nbuf,ibufdim,mydist2_.myhost,
                 mydist2_.myid,mydist2_.mycomm,mblk2nd.data(),maxbl);
    }

    // ---- write out connectivity file (patch.bin, unit 2) -------------------
    #define IINDEX(n,c) iindex[((n)-1)+((c)-1)*intmax]
    FILE* f2 = FUNIT(2);
    {
        int v = ninter;
        tc::urec_write(f2,&v,sizeof(int));               // write(2) ninter
    }
    for (int n=1;n<=ninter;++n){
        int nfb = IINDEX(n,1);
        { int v=IINDEX(n,1); tc::urec_write(f2,&v,sizeof(int)); }        // write(2) iindex(n,1)
        { tc::RecordWriter rw; for(int ll=1;ll<=nfb;++ll) rw.put_i(IINDEX(n,1+ll)); rw.flush(f2); }
        { int v=IINDEX(n,nfb+2); tc::urec_write(f2,&v,sizeof(int)); }    // write(2) iindex(n,nfb+2)
        { tc::RecordWriter rw; for(int ll=1;ll<=nfb;++ll) rw.put_i(IINDEX(n,nfb+2+ll)); rw.flush(f2); }
        { int v=IINDEX(n,2*nfb+3); tc::urec_write(f2,&v,sizeof(int)); }
        { int v=IINDEX(n,2*nfb+4); tc::urec_write(f2,&v,sizeof(int)); }
        { int v=IINDEX(n,2*nfb+5); tc::urec_write(f2,&v,sizeof(int)); }
        { int v=IINDEX(n,2*nfb+6); tc::urec_write(f2,&v,sizeof(int)); }
        { int v=IINDEX(n,2*nfb+7); tc::urec_write(f2,&v,sizeof(int)); }
        { int v=IINDEX(n,2*nfb+8); tc::urec_write(f2,&v,sizeof(int)); }
        { int v=IINDEX(n,2*nfb+9); tc::urec_write(f2,&v,sizeof(int)); }
        int lst = IINDEX(n,2*nfb+5);
        int len = lst + IINDEX(n,2*nfb+4) - 1;
        { tc::RecordWriter rw; for(int nnn=lst;nnn<=len;++nnn) rw.put_i(nblkpt[nnn-1]); rw.flush(f2); }
        // write(2) ((windex(nnn,ll),nnn=lst,len),ll=1,2)  — column-major order
        { tc::RecordWriter rw;
          for (int ll=1;ll<=2;++ll)
              for (int nnn=lst;nnn<=len;++nnn)
                  rw.put_d(windex[(nnn-1)+(ll-1)*maxxe]);
          rw.flush(f2); }
    }

    // angular displacements for ioflag=2
    if (ioflag>1){
        for (int n=1;n<=ninter;++n){
            int nfb=IINDEX(n,1);
            for (int ll=1;ll<=nfb;++ll){
                tc::RecordWriter rw;
                rw.put_d(dthetx[(n-1)+(ll-1)*intmax]);
                rw.put_d(dthety[(n-1)+(ll-1)*intmax]);
                rw.put_d(dthetz[(n-1)+(ll-1)*intmax]);
                rw.flush(f2);
            }
        }
    }

    std::fprintf(f25,"\n");                       // write(25,*)
    // format 9998: (2x,26h..,/,2x,a60,/)
    { std::string pn=filenam_.patch.s; pn.resize(60,' ');
      std::fprintf(f25,"  writing connectivity file:\n  %s\n\n",pn.c_str()); }
    {
        int nfb=IINDEX(ninter,1);
        int lst=IINDEX(ninter,2*nfb+5);
        int len=lst + IINDEX(ninter,2*nfb+4) - 1;
        // format 9997: ('  minimum dimension for parameter MAXXE in CFL3D:',i6)
        std::fprintf(f25,"  minimum dimension for parameter MAXXE in CFL3D:%s\n",
                     f_iw(len,6).c_str());
    }

    string=FStr(50);
    string.assign("    timing for complete run - time in seconds     ");
    cputim(-1,mydist2_.nnodes,string,mydist2_.myhost,mydist2_.myid,mydist2_.mycomm,25);

    if(f2){ std::fflush(f2); std::fclose(f2); FUNIT(2)=nullptr; }

    #undef LW
    #undef LW2
    #undef IINDEX
    return;
}
