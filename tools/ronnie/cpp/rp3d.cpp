// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// rp3d.cpp — read grids in plot3d format (real build; unit 1, unformatted).
// x,y,z are (jdim,kdim,idim). One Fortran READ = one unformatted record holding
// the whole x block, then the whole y (or z) block, then z (or y).
//   ialph == 0 : record order x, y, z   (cfl3d standard)
//   ialph  > 0 : record order x, z, y   then negate y (tlns3d standard)
#include "ron_common.h"

void rp3d(double* x,double* y,double* z,int jdim,int kdim,int idim,
          int /*igrid*/,int ialph,int /*igeom_img*/,int& irr){
    auto X=[&](int j,int k,int i)->double&{
        return x[(j-1)+(k-1)*jdim+(i-1)*jdim*kdim]; };
    auto Y=[&](int j,int k,int i)->double&{
        return y[(j-1)+(k-1)*jdim+(i-1)*jdim*kdim]; };
    auto Z=[&](int j,int k,int i)->double&{
        return z[(j-1)+(k-1)*jdim+(i-1)*jdim*kdim]; };

    size_t n=(size_t)jdim*kdim*idim;
    std::vector<double> buf(3*n);
    FILE* f=FUNIT(1);
    if(!f || !tc::urec_read_try(f,buf.data(),3*n*sizeof(double))){
        irr=1; return;
    }
    // fill in the loop order (((.(j,k,i),i=1,idim),j=1,jdim),k=1,kdim):
    // the record streams a1 = first block, a2 = second block, a3 = third.
    // element order within a block is i fastest, then j, then k.
    size_t p=0;
    // block A -> x
    for(int k=1;k<=kdim;++k) for(int j=1;j<=jdim;++j) for(int i=1;i<=idim;++i)
        X(j,k,i)=buf[p++];
    if(ialph==0){
        for(int k=1;k<=kdim;++k) for(int j=1;j<=jdim;++j) for(int i=1;i<=idim;++i)
            Y(j,k,i)=buf[p++];
        for(int k=1;k<=kdim;++k) for(int j=1;j<=jdim;++j) for(int i=1;i<=idim;++i)
            Z(j,k,i)=buf[p++];
    } else {
        // order x, z, y
        for(int k=1;k<=kdim;++k) for(int j=1;j<=jdim;++j) for(int i=1;i<=idim;++i)
            Z(j,k,i)=buf[p++];
        for(int k=1;k<=kdim;++k) for(int j=1;j<=jdim;++j) for(int i=1;i<=idim;++i)
            Y(j,k,i)=buf[p++];
        for(int i=1;i<=idim;++i) for(int j=1;j<=jdim;++j) for(int k=1;k<=kdim;++k)
            Y(j,k,i)=-Y(j,k,i);
    }
    irr=0;
}
