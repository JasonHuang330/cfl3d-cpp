// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// loadgr.cpp — load a patch-surface grid from 1-d storage w into 2-d x,y,z.
// w is the flat grid store; lx,ly,lz are 1-based start indices into w. x,y,z
// are (mdim,ndim). Storage formula matches CONTRACT: l = (i-1)*jdimg*kdimg +
// (k-1)*jdimg + (j-1), value = w(lx+l) (lx is 1-based, w accessed 1-based).
#include "ron_common.h"

void loadgr(double* w,int /*mgwk*/,int lx,int ly,int lz,int jindex,
            double* x,double* y,double* z,int mdim,int /*ndim*/,
            int idimg,int jdimg,int kdimg){
    // w passed as base pointer to element 1; w(lx+l) -> w[(lx+l)-1]
    auto WV=[&](int idx)->double{ return w[idx-1]; };
    auto X=[&](int a,int b)->double&{ return x[(a-1)+(b-1)*mdim]; };
    auto Y=[&](int a,int b)->double&{ return y[(a-1)+(b-1)*mdim]; };
    auto Z=[&](int a,int b)->double&{ return z[(a-1)+(b-1)*mdim]; };

    // patch surface is an i=constant surface
    if (jindex/10==1){
        int i = (jindex==11)? 1 : idimg;
        for(int j=1;j<=jdimg;++j) for(int k=1;k<=kdimg;++k){
            int l=(i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1);
            X(j,k)=WV(lx+l); Y(j,k)=WV(ly+l); Z(j,k)=WV(lz+l);
        }
    }
    // patch surface is a j=constant surface
    if (jindex/10==2){
        int j = (jindex==21)? 1 : jdimg;
        for(int i=1;i<=idimg;++i) for(int k=1;k<=kdimg;++k){
            int l=(i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1);
            X(k,i)=WV(lx+l); Y(k,i)=WV(ly+l); Z(k,i)=WV(lz+l);
        }
    }
    // patch surface is a k=constant surface
    if (jindex/10==3){
        int k = (jindex==31)? 1 : kdimg;
        for(int i=1;i<=idimg;++i) for(int j=1;j<=jdimg;++j){
            int l=(i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1);
            X(j,i)=WV(lx+l); Y(j,i)=WV(ly+l); Z(j,i)=WV(lz+l);
        }
    }
}
