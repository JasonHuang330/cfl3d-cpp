// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// collx.cpp — restrict x,y,z to a coarser mesh (every other point).
// x,y,z are (jdim,kdim,idim); xx,yy,zz are (jj2,kk2,ii2). Column-major, 1-based.
#include "ron_common.h"

void collx(double* x,double* y,double* z,double* xx,double* yy,double* zz,
           int jdim,int kdim,int idim,int jj2,int kk2,int ii2){
    auto X=[&](int j,int k,int i)->double&{
        return x[(j-1)+(k-1)*jdim+(i-1)*jdim*kdim]; };
    auto Y=[&](int j,int k,int i)->double&{
        return y[(j-1)+(k-1)*jdim+(i-1)*jdim*kdim]; };
    auto Z=[&](int j,int k,int i)->double&{
        return z[(j-1)+(k-1)*jdim+(i-1)*jdim*kdim]; };
    auto XX=[&](int j,int k,int i)->double&{
        return xx[(j-1)+(k-1)*jj2+(i-1)*jj2*kk2]; };
    auto YY=[&](int j,int k,int i)->double&{
        return yy[(j-1)+(k-1)*jj2+(i-1)*jj2*kk2]; };
    auto ZZ=[&](int j,int k,int i)->double&{
        return zz[(j-1)+(k-1)*jj2+(i-1)*jj2*kk2]; };

    int iinc=2;
    if (idim==2) iinc=1;
    int ii=0;
    for(int i=1;i<=idim;i+=iinc){
        ii=ii+1;
        int kk=0;
        for(int k=1;k<=kdim;k+=2){
            kk=kk+1;
            int jj=0;
            for(int j=1;j<=jdim;j+=2){
                jj=jj+1;
                XX(jj,kk,ii)=X(j,k,i);
                YY(jj,kk,ii)=Y(j,k,i);
                ZZ(jj,kk,ii)=Z(j,k,i);
            }
        }
    }
}
