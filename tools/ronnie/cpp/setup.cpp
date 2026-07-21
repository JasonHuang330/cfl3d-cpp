// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// setup.cpp — read grids and compute preliminary 1-d array indices (unit 25).
// lw(65,maxbl), lw2(43,maxbl), w(mgwk). itest/jtest/ktest/ncgg/iemg/nblg(maxgr).
#include "ron_common.h"

void setup(int* lw,int* lw2,double* w,int /*mgwk*/,int* itest,int* jtest,int* ktest,
           int maxbl,int /*maxgr*/,int ngrid,int* ncgg,int* iemg,int* nblg){
    FILE* f25=FUNIT(25);
    // format 92: (/,1x,16hREADING IN GRIDS)
    std::fprintf(f25,"\n READING IN GRIDS\n");

    for(int igrid=1;igrid<=ngrid;++igrid){
        int iem=iemg[igrid-1];  (void)iem;
        int nbl=nblg[igrid-1];
        lead(nbl,lw,lw2,maxbl);

        if (igrdtyp_.ip3dgrd==0){
            // cfl3d type: one record of jdum,kdum,idum
            FILE* f1=FUNIT(1);
            int hdr[3]={0,0,0};
            tc::urec_read(f1,hdr,3*sizeof(int));
            int jdum=hdr[0],kdum=hdr[1],idum=hdr[2];
            // format 93: (/,13h reading grid,i3,24h of dimensions (I/J/K) :,3i5)
            std::fprintf(f25,"\n reading grid%s of dimensions (I/J/K) :%s%s%s\n",
                         f_iw(igrid,3).c_str(),f_iw(idum,5).c_str(),
                         f_iw(jdum,5).c_str(),f_iw(kdum,5).c_str());
            if (jdum!=ginfo_.jdim || kdum!=ginfo_.kdim || idum!=ginfo_.idim){
                std::fprintf(f25," stopping....inconsistency in grid data file\n");
                std::exit(1);
            }
            // read(1) (w(i),i=lx,lvis-1) : one record, (lvis-1)-lx+1 doubles
            int lx=ginfo_.lx, lvis=ginfo_.lvis;
            int cnt=lvis-1 - lx + 1;
            // w passed as base pointer to element 1; w(i) -> w[i-1]
            tc::urec_read(f1,&w[lx-1],(size_t)cnt*sizeof(double));
        } else {
            // plot3d type
            FILE* f1=FUNIT(1);
            if (igrid==1){
                int ndum=0;
                tc::urec_read(f1,&ndum,sizeof(int));
                if (ndum!=ngrid){
                    std::fprintf(f25," stopping....ngrid = %d but grid file"
                                 " contains %d grids\n",ngrid,ndum);
                    std::exit(1);
                }
                // read(1) (itest(ll),jtest(ll),ktest(ll),ll=1,ngrid)
                std::vector<int> tmp(3*ngrid);
                tc::urec_read(f1,tmp.data(),(size_t)3*ngrid*sizeof(int));
                for(int ll=1;ll<=ngrid;++ll){
                    itest[ll-1]=tmp[(ll-1)*3+0];
                    jtest[ll-1]=tmp[(ll-1)*3+1];
                    ktest[ll-1]=tmp[(ll-1)*3+2];
                }
            }
            std::fprintf(f25,"\n reading grid%s of dimensions (I/J/K) :%s%s%s\n",
                         f_iw(igrid,3).c_str(),f_iw(itest[igrid-1],5).c_str(),
                         f_iw(jtest[igrid-1],5).c_str(),f_iw(ktest[igrid-1],5).c_str());
            if (jtest[igrid-1]!=ginfo_.jdim || ktest[igrid-1]!=ginfo_.kdim ||
                itest[igrid-1]!=ginfo_.idim){
                std::fprintf(f25," stopping....inconsistency in grid data file\n");
                std::exit(1);
            }
            igrdtyp_.ialph=1;   // same as old ronnie
            int irr=0;
            rp3d(&w[ginfo_.lx-1],&w[ginfo_.ly-1],&w[ginfo_.lz-1],
                 ginfo_.jdim,ginfo_.kdim,ginfo_.idim,igrid,igrdtyp_.ialph,0,irr);
            if (irr!=0){
                // write(11,*) messages then stop
                FILE* f11=FUNIT(11);
                if(f11){
                    std::fprintf(f11," Stopping... error reading grid...\n");
                    std::fprintf(f11," (Common error:  grid not written in same"
                                     " precision that RONNIE was compiled)\n");
                }
                std::exit(1);
            }
        }

        // coarser global grids
        int ncg=ncgg[igrid-1];
        if (ncg>0){
            for(int m=1;m<=ncg;++m){
                nbl=nbl+1;
                int lxc=lw[(10-1)+(nbl-1)*65];
                int lyc=lw[(11-1)+(nbl-1)*65];
                int lzc=lw[(12-1)+(nbl-1)*65];
                collx(&w[ginfo_.lx-1],&w[ginfo_.ly-1],&w[ginfo_.lz-1],
                      &w[lxc-1],&w[lyc-1],&w[lzc-1],
                      ginfo_.jdim,ginfo_.kdim,ginfo_.idim,
                      ginfo_.jj2,ginfo_.kk2,ginfo_.ii2);
                lead(nbl,lw,lw2,maxbl);
            }
        }
    }
}
