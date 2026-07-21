// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// main.cpp — C++ driver for the CFL3D ronnie patched-grid preprocessor.
// Faithful translation of ronnie/main.F (program main).
//
// The Fortran driver:
//   * copies stdin into a temporary input file (unit iunit5=55)
//   * opens grid.bin (unit 1), output=ronnie.out (unit 25), patch.bin (unit 2),
//     and ronnie.error (unit 99)
//   * pre-reads the deck to size maxgr0/maxbl0/intmax0/nsub10
//   * calls sizer(...) then ronnie(...) then termn8
//
// signal() handlers are no-ops here (skipped, as instructed).
#include "ron_common.h"

int main(){
    const int ibufdim = 2000, nbuf = 4;

    unit5_.iunit5 = 55;
    int iunit5 = unit5_.iunit5;

    // ---- copy stdin -> temporary input file (tempz3y2x1.inp) ---------------
    // open(iunit5,file='tempz3y2x1.inp',...); loop read(5,'(a80)') / write(iunit5)
    {
        FILE* tmp = std::fopen("tempz3y2x1.inp","w+");
        FUNIT(iunit5) = tmp;
        char line[8192];
        for (int n=1; n<=9999; ++n){
            if (!std::fgets(line,sizeof(line),stdin)) break;
            // strip newline, then emulate read(5,'(a80)')/write(iunit5,'(a80)'):
            // Fortran reads exactly 80 chars (a80) and writes them back a80.
            size_t L=std::strlen(line);
            while (L>0 && (line[L-1]=='\n'||line[L-1]=='\r')) line[--L]='\0';
            // pad/truncate to 80 like a80 field
            char f80[81];
            for (int i=0;i<80;++i) f80[i] = (i<(int)L)? line[i] : ' ';
            f80[80]='\0';
            // trailing blanks are significant in the file record; write 80 then \n
            std::fwrite(f80,1,80,tmp); std::fputc('\n',tmp);
        }
        std::rewind(tmp);
    }

    FILE* f5 = FUNIT(iunit5);
    tc::TextReader tr(f5);

    tr.line();                                   // read(iunit5,*)

    auto rd_a60 = [&](FStr& dst){ std::string s=tr.line();
        if((int)s.size()>60) s=s.substr(0,60); dst=FStr(80); dst.assign(s); };

    rd_a60(filenam_.grid);                       // read(iunit5,'(a60)')grid
    {
        std::string gn=filenam_.grid.s; size_t p=gn.find_last_not_of(' ');
        gn = (p==std::string::npos)? std::string() : gn.substr(0,p+1);
        FUNIT(1) = std::fopen(gn.c_str(),"rb");  // form='unformatted'
    }

    rd_a60(filenam_.output);                     // read(iunit5,'(a60)')output
    // open, delete, reopen (truncate) the formatted output file:
    {
        std::string on=filenam_.output.s; size_t p=on.find_last_not_of(' ');
        on = (p==std::string::npos)? std::string() : on.substr(0,p+1);
        FUNIT(25) = std::fopen(on.c_str(),"w");
    }
    filenam_.grdmov = filenam_.output;

    rd_a60(filenam_.patch);                      // read(iunit5,'(a60)')patch
    {
        std::string pn=filenam_.patch.s; size_t p=pn.find_last_not_of(' ');
        pn = (p==std::string::npos)? std::string() : pn.substr(0,p+1);
        FUNIT(2) = std::fopen(pn.c_str(),"wb");  // form='unformatted'
    }

    filenam_.errfile = FStr(80); filenam_.errfile.assign("ronnie.error");
    FUNIT(99) = std::fopen("ronnie.error","w");

    // ---- pre-read deck to determine array-size requirements ----------------
    tr.line();                                   // read(iunit5,*)
    int ioflag=0,itrace=0;
    { int v[2]; tr.ints(v,2); ioflag=v[0]; itrace=v[1]; } // read ioflag,itrace
    tr.line();                                   // read(iunit5,*)
    tr.line();                                   // read(iunit5,*)
    int ngrid=0;
    { int v[1]; tr.ints(v,1); ngrid=v[0]; }      // read ngrid
    ngrid = std::abs(ngrid);

    int maxgr0 = ngrid;

    tr.line();                                   // read(iunit5,*)
    int ncgmax=0;
    for (int ll=1; ll<=ngrid; ++ll){
        int v[5]; tr.ints(v,5);                  // ncg,iem,idim,jdim,kdim
        int ncg=v[0];
        ncgmax = i_max(ncgmax,ncg);
    }

    int maxbl0 = ngrid*(ncgmax+1);

    tr.line();                                   // read(iunit5,*)
    int nint0=0;
    { int v[1]; tr.ints(v,1); nint0=v[0]; }      // read nint0

    int intmax0 = nint0*(ncgmax+1);

    tr.line();                                   // read(iunit5,*)
    for (int n=1; n<=nint0; ++n){
        int v[8]; tr.ints(v,8);                  // int,ifit,limit,itmax,mcxie,mceta,ic0,iorph
    }

    int nsub10 = 1;
    if (nint0>0){
        if (ioflag==0){
            tr.line();                           // read(iunit5,*)
            for (int n=1;n<=nint0;++n){
                int v[7]; tr.ints(v,7);          // int,ito,i1,i2,j1,j2,nfb
                nsub10 = i_max(nsub10,v[6]);
            }
        }
        if (ioflag==1){
            tr.line(); tr.line();
            for (int n=1;n<=nint0;++n){
                int v[7]; tr.ints(v,7);
                int nfb=v[6];
                for (int nn=1;nn<=nfb;++nn) tr.line();
                nsub10 = i_max(nsub10,nfb);
            }
        }
        if (ioflag==2){
            for (int n=1;n<=nint0;++n){
                tr.line();
                int v[7]; tr.ints(v,7);
                int nfb=v[6];
                for (int nn=1;nn<=nfb;++nn){ tr.line(); tr.line(); tr.line(); tr.line(); }
                nsub10 = i_max(nsub10,nfb);
            }
        }
    }

    std::rewind(f5);

    int nwork0=1, mpatch0=1, maxxe0=1, ibufdim0=ibufdim, nbuf0=nbuf;
    (void)nwork0;

    int maxgr=0,maxbl=0,maxxe=0,nsub1=0,intmax=0,mpatch=0,nwork=0;
    int imode=1;
    sizer(maxgr,maxbl,maxxe,nsub1,intmax,mpatch,nwork,
          nbuf0,ibufdim0,maxgr0,maxbl0,maxxe0,nsub10,
          intmax0,mpatch0,imode);

    // ---- call ronnie main routine -----------------------------------------
    ronnie(maxbl,maxgr,nwork,intmax,nsub1,maxxe,mpatch,nbuf,ibufdim);

    // ---- normal termination -----------------------------------------------
    if (FUNIT(iunit5)){ std::fclose(FUNIT(iunit5)); FUNIT(iunit5)=nullptr; }
    std::remove("tempz3y2x1.inp");

    {
        std::vector<int> nou((size_t)nbuf,0);
        std::vector<FStr> bou((size_t)ibufdim*nbuf,FStr(120));
        termn8(mydist2_.myid,0,ibufdim,nbuf,bou.data(),nou.data());
    }
    return 0;
}
