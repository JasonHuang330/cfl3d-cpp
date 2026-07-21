// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// xe.cpp — faithful C++ translation of xe.F (CFL3D ronnie).
// Select the proper pair of coordinates (2 of x,y,z) to use for the inversion
// in the "to" cell, project "from" cell nodes onto the "to" cell plane, then
// dispatch to xe2 with the appropriate coordinate pair.
//
// Sequential real build: real(x)->x.
// /areas/ ap(3),imaxa -> areas_.ap[0..2], areas_.imaxa.
// Arrays x/y/z, xmid/ymid/zmid, xmide/ymide/zmide are (jdim,kdim,nsub):
//   (i,j,k) -> [(i-1)+(j-1)*jdim+(k-1)*jdim*kdim].
// jcell=jp, kcell=kp (1-based cell indices).  bou(ibufdim,nbuf), nou(nbuf).
#include "ron_common.h"

void direct(double x5,double x6,double x7,double x8,
            double y5,double y6,double y7,double y8,
            double z5,double z6,double z7,double z8,
            double& a1,double& a2,double& a3,int& imaxa,
            int* nou,FStr* bou,int nbuf,int ibufdim);
void project(double x0,double y0,double z0,double xc,double yc,double zc,
             double a1,double a2,double a3,
             double& xcp,double& ycp,double& zcp);
void xe2(double x1,double x2,double x3,double x4,double x5,double x6,
         double x7,double x8,double xc,double y1,double y2,double y3,
         double y4,double y5,double y6,double y7,double y8,double yc,
         double& xie,double& eta,int& imiss,int ifit,
         int* nou,FStr* bou,int nbuf,int ibufdim,int myid);

void xe(int jdim,int kdim,int nsub,int l,double* x,double* y,double* z,
        double* xmid,double* ymid,double* zmid,double* xmide,double* ymide,
        double* zmide,int jcell,int kcell,double xc,double yc,double zc,
        double& xie,double& eta,int& imiss,int ifit,int ic0,
        int* nou,FStr* bou,int nbuf,int ibufdim,int myid)
{
    (void)nsub;
    const int d1 = jdim, d2 = kdim;
    #define X(i,j,k)     x    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y(i,j,k)     y    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z(i,j,k)     z    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define XMID(i,j,k)  xmid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define YMID(i,j,k)  ymid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define ZMID(i,j,k)  zmid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define XMIDE(i,j,k) xmide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define YMIDE(i,j,k) ymide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define ZMIDE(i,j,k) zmide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]

    double x1 = X(jcell,kcell,l);
    double y1 = Y(jcell,kcell,l);
    double z1 = Z(jcell,kcell,l);
    double x2 = X(jcell+1,kcell,l);
    double y2 = Y(jcell+1,kcell,l);
    double z2 = Z(jcell+1,kcell,l);
    double x4 = X(jcell,kcell+1,l);
    double y4 = Y(jcell,kcell+1,l);
    double z4 = Z(jcell,kcell+1,l);
    double x3 = X(jcell+1,kcell+1,l);
    double y3 = Y(jcell+1,kcell+1,l);
    double z3 = Z(jcell+1,kcell+1,l);
    double x5 = XMID(jcell,kcell,l);
    double y5 = YMID(jcell,kcell,l);
    double z5 = ZMID(jcell,kcell,l);
    double x6 = XMID(jcell,kcell+1,l);
    double y6 = YMID(jcell,kcell+1,l);
    double z6 = ZMID(jcell,kcell+1,l);
    double x7 = XMIDE(jcell,kcell,l);
    double y7 = YMIDE(jcell,kcell,l);
    double z7 = ZMIDE(jcell,kcell,l);
    double x8 = XMIDE(jcell+1,kcell,l);
    double y8 = YMIDE(jcell+1,kcell,l);
    double z8 = ZMIDE(jcell+1,kcell,l);

    // compute normalized directed areas/unit normals of "from" cell
    double a1,a2,a3;
    int itoss;
    direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,
           a1,a2,a3,itoss,nou,bou,nbuf,ibufdim);
    double a[3];
    a[1-1] = a1;
    a[2-1] = a2;
    a[3-1] = a3;

    // check inner product of unit normals for compatible orientation before
    // projecting "from" cell nodes onto plane defined by "to" cell
    double tol  = 0.1;
    double prod = areas_.ap[1-1]*a[1-1] + areas_.ap[2-1]*a[2-1]
                + areas_.ap[3-1]*a[3-1];

    if (std::fabs(prod) > tol) {
        if (ic0 == 0) {
            project(xc,yc,zc,x1,y1,z1,areas_.ap[0],areas_.ap[1],areas_.ap[2],x1,y1,z1);
            project(xc,yc,zc,x2,y2,z2,areas_.ap[0],areas_.ap[1],areas_.ap[2],x2,y2,z2);
            project(xc,yc,zc,x3,y3,z3,areas_.ap[0],areas_.ap[1],areas_.ap[2],x3,y3,z3);
            project(xc,yc,zc,x4,y4,z4,areas_.ap[0],areas_.ap[1],areas_.ap[2],x4,y4,z4);
            project(xc,yc,zc,x5,y5,z5,areas_.ap[0],areas_.ap[1],areas_.ap[2],x5,y5,z5);
            project(xc,yc,zc,x6,y6,z6,areas_.ap[0],areas_.ap[1],areas_.ap[2],x6,y6,z6);
            project(xc,yc,zc,x7,y7,z7,areas_.ap[0],areas_.ap[1],areas_.ap[2],x7,y7,z7);
            project(xc,yc,zc,x8,y8,z8,areas_.ap[0],areas_.ap[1],areas_.ap[2],x8,y8,z8);
        }
    } else {
        imiss = 1;
        xie = xie + 1.;
        eta = eta + 1.;
        return;
    }

    // use best coordinates for inversion... toss out equation for which
    // directed area ap(i) is maximum.
    itoss = areas_.imaxa;

    if (itoss == 1) {
        // use only y and z equations
        xe2(y1,y2,y3,y4,y5,y6,y7,y8,yc,z1,z2,z3,z4,z5,
            z6,z7,z8,zc,xie,eta,imiss,ifit,nou,bou,nbuf,ibufdim,myid);
    }

    if (itoss == 2) {
        // use only x and z equations
        xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,z1,z2,z3,z4,z5,
            z6,z7,z8,zc,xie,eta,imiss,ifit,nou,bou,nbuf,ibufdim,myid);
    }

    if (itoss == 3) {
        // use only x and y equations
        xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,y1,y2,y3,y4,y5,
            y6,y7,y8,yc,xie,eta,imiss,ifit,nou,bou,nbuf,ibufdim,myid);
    }
    return;
    #undef X
    #undef Y
    #undef Z
    #undef XMID
    #undef YMID
    #undef ZMID
    #undef XMIDE
    #undef YMIDE
    #undef ZMIDE
}
