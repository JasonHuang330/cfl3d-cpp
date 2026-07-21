// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "mms.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace mms_ns {

// Helper: ccabs for double (absolute value)
static inline double ccabs(double x) { return std::abs(x); }

void exact_norm(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
    FortranArray4DRef<double> q, FortranArray3DRef<double> tursav,
    FortranArray3DRef<double> vol, FortranArray4DRef<double> res,
    FortranArray3DRef<double> smin,
    int& iexact_trunc, int& iexact_disc)
{
    double xmach = (double)cmn_info.xmach;
    double gamma = (double)cmn_fluid.gamma;
    double reue  = (double)cmn_reyue.reue;

    FortranArray1D<double> fvec(6);
    double u_max=0., w_max=0., p_max=0., t_max=0.;
    int i_u_max=0, j_u_max=0, k_u_max=0;
    int i_w_max=0, j_w_max=0, k_w_max=0;
    int i_p_max=0, j_p_max=0, k_p_max=0;
    int i_t_max=0, j_t_max=0, k_t_max=0;
    double sum_u=0., sum_w=0., sum_p=0., sum_t=0.;

    int i_convert_q = 0;
    int i_forcing = 0;
    int i_gradient = 0;
    int nvar=6;
    for (int i=1; i<=idim-1; i++) {
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double distf=smin(j,k,i);
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);

                double delta_u=ccabs(q(j,k,i,2)-fvec(2))/xmach;
                double delta_w=ccabs(q(j,k,i,4)-fvec(4))/xmach;
                double delta_p=ccabs(q(j,k,i,5)-fvec(5))*gamma/(xmach*xmach);
                double delta_t=ccabs(tursav(j,k,i)-fvec(6))/reue;
                if (delta_u > u_max) {
                    u_max=delta_u; i_u_max=i; j_u_max=j; k_u_max=k;
                }
                if (delta_w > w_max) {
                    w_max=delta_w; i_w_max=i; j_w_max=j; k_w_max=k;
                }
                if (delta_p > p_max) {
                    p_max=delta_p; i_p_max=i; j_p_max=j; k_p_max=k;
                }
                if (delta_t > t_max) {
                    t_max=delta_t; i_t_max=i; j_t_max=j; k_t_max=k;
                }
                sum_u=sum_u+delta_u;
                sum_w=sum_w+delta_w;
                sum_p=sum_p+delta_p;
                sum_t=sum_t+delta_t;
            }
        }
    }
    sum_u=sum_u/((jdim-1)*(kdim-1)*(idim-1));
    sum_w=sum_w/((jdim-1)*(kdim-1)*(idim-1));
    sum_p=sum_p/((jdim-1)*(kdim-1)*(idim-1));
    sum_t=sum_t/((jdim-1)*(kdim-1)*(idim-1));

    FILE* f201 = fortran_get_unit(201);
    fprintf(f201,"# discretization error (compared to exact soln)\n");
    fprintf(f201,"variables=\"N\",\"1/N\",\"sqrt(1/N)\",");
    fprintf(f201,"\"u_max\",\"i_u_max\",\"j_u_max\",\"k_u_max\",\"sum_u\",\n");
    fprintf(f201,"\"w_max\",\"i_w_max\",\"j_w_max\",\"k_w_max\",\"sum_w\",");
    fprintf(f201,"\"p_max\",\"i_p_max\",\"j_p_max\",\"k_p_max\",\"sum_p\",\n");
    fprintf(f201,"\"t_max\",\"i_t_max\",\"j_t_max\",\"k_t_max\",\"sum_t\"\n");
    int inumtot=(jdim-1)*(kdim-1);
    double hhhh=1./(double)((jdim-1)*(kdim-1));
    double sqhhhh=sqrt(hhhh);
    fprintf(f201,"%9d%18.8e%18.8e%18.8e%8d%8d%8d%18.8e%18.8e%8d%8d%8d%18.8e%18.8e%8d%8d%8d%18.8e%18.8e%5d%5d%5d%18.8e\n",
        inumtot,hhhh,sqhhhh,
        u_max,i_u_max,j_u_max,k_u_max,sum_u,
        w_max,i_w_max,j_w_max,k_w_max,sum_w,
        p_max,i_p_max,j_p_max,k_p_max,sum_p,
        t_max,i_t_max,j_t_max,k_t_max,sum_t);

    // residuals
    double res_1_max=0., res_2_max=0., res_3_max=0., res_4_max=0., res_5_max=0.;
    double sum_1_res=0., sum_2_res=0., sum_3_res=0., sum_4_res=0., sum_5_res=0.;
    int i_res_1_max=0, j_res_1_max=0, k_res_1_max=0;
    int i_res_2_max=0, j_res_2_max=0, k_res_2_max=0;
    int i_res_3_max=0, j_res_3_max=0, k_res_3_max=0;
    int i_res_4_max=0, j_res_4_max=0, k_res_4_max=0;
    int i_res_5_max=0, j_res_5_max=0, k_res_5_max=0;
    double res_1_max_alone=0., vol_1_max_alone=0.;
    double res_2_max_alone=0., vol_2_max_alone=0.;
    double res_3_max_alone=0., vol_3_max_alone=0.;
    double res_4_max_alone=0., vol_4_max_alone=0.;
    double res_5_max_alone=0., vol_5_max_alone=0.;
    {
        int i=1;
        for (int j=3; j<=jdim-3; j++) {
            for (int k=3; k<=kdim-3; k++) {
                if (ccabs(res(j,k,i,1)/vol(j,k,i)) > res_1_max) {
                    res_1_max = ccabs(res(j,k,i,1)/vol(j,k,i));
                    i_res_1_max=i; j_res_1_max=j; k_res_1_max=k;
                    res_1_max_alone=res(j,k,i,1); vol_1_max_alone=vol(j,k,i);
                }
                if (ccabs(res(j,k,i,2)/vol(j,k,i)) > res_2_max) {
                    res_2_max = ccabs(res(j,k,i,2)/vol(j,k,i));
                    i_res_2_max=i; j_res_2_max=j; k_res_2_max=k;
                    res_2_max_alone=res(j,k,i,2); vol_2_max_alone=vol(j,k,i);
                }
                if (ccabs(res(j,k,i,3)/vol(j,k,i)) > res_3_max) {
                    res_3_max = ccabs(res(j,k,i,3)/vol(j,k,i));
                    i_res_3_max=i; j_res_3_max=j; k_res_3_max=k;
                    res_3_max_alone=res(j,k,i,3); vol_3_max_alone=vol(j,k,i);
                }
                if (ccabs(res(j,k,i,4)/vol(j,k,i)) > res_4_max) {
                    res_4_max = ccabs(res(j,k,i,4)/vol(j,k,i));
                    i_res_4_max=i; j_res_4_max=j; k_res_4_max=k;
                    res_4_max_alone=res(j,k,i,4); vol_4_max_alone=vol(j,k,i);
                }
                if (ccabs(res(j,k,i,5)/vol(j,k,i)) > res_5_max) {
                    res_5_max = ccabs(res(j,k,i,5)/vol(j,k,i));
                    i_res_5_max=i; j_res_5_max=j; k_res_5_max=k;
                    res_5_max_alone=res(j,k,i,5); vol_5_max_alone=vol(j,k,i);
                }
                sum_1_res=sum_1_res+ccabs(res(j,k,i,1)/vol(j,k,i));
                sum_2_res=sum_2_res+ccabs(res(j,k,i,2)/vol(j,k,i));
                sum_3_res=sum_3_res+ccabs(res(j,k,i,3)/vol(j,k,i));
                sum_4_res=sum_4_res+ccabs(res(j,k,i,4)/vol(j,k,i));
                sum_5_res=sum_5_res+ccabs(res(j,k,i,5)/vol(j,k,i));
            }
        }
    }
    sum_1_res=sum_1_res/((jdim-5)*(kdim-5)*(idim-1));
    sum_2_res=sum_2_res/((jdim-5)*(kdim-5)*(idim-1));
    sum_3_res=sum_3_res/((jdim-5)*(kdim-5)*(idim-1));
    sum_4_res=sum_4_res/((jdim-5)*(kdim-5)*(idim-1));
    sum_5_res=sum_5_res/((jdim-5)*(kdim-5)*(idim-1));

    FILE* f301 = fortran_get_unit(301);
    fprintf(f301,"# mean flow residuals - truncation error if iexact_trunc=1\n");
    fprintf(f301,"variables=\"N\",\"1/N\",\"sqrt(1/N)\",\n");
    fprintf(f301,"\"res1_max\",\"i_res1_max\",\"j_res1_max\",");
    fprintf(f301,"\"k_res1_max\",\"res1_L1norm\",\"res_1_max_alone\",");
    fprintf(f301,"\"vol_1_max_alone\",\n");
    fprintf(f301,"\"res2_max\",\"i_res2_max\",\"j_res2_max\",");
    fprintf(f301,"\"k_res2_max\",\"res2_L1norm\",\"res_2_max_alone\",");
    fprintf(f301,"\"vol_2_max_alone\",\n");
    fprintf(f301,"\"res3_max\",\"i_res3_max\",\"j_res3_max\",");
    fprintf(f301,"\"k_res3_max\",\"res3_L1norm\",\"res_3_max_alone\",");
    fprintf(f301,"\"vol_3_max_alone\",\n");
    fprintf(f301,"\"res4_max\",\"i_res4_max\",\"j_res4_max\",");
    fprintf(f301,"\"k_res4_max\",\"res4_L1norm\",\"res_4_max_alone\",");
    fprintf(f301,"\"vol_4_max_alone\",\n");
    fprintf(f301,"\"res5_max\",\"i_res5_max\",\"j_res5_max\",");
    fprintf(f301,"\"k_res5_max\",\"res5_L1norm\",\"res_5_max_alone\",");
    fprintf(f301,"\"vol_5_max_alone\"\n");
    inumtot=(jdim-1)*(kdim-1);
    hhhh=1./(double)((jdim-1)*(kdim-1));
    sqhhhh=sqrt(hhhh);
    fprintf(f301,"%9d%18.8e%18.8e"
        "%18.8e%5d%5d%5d%18.8e%18.8e%18.8e"
        "%18.8e%5d%5d%5d%18.8e%18.8e%18.8e"
        "%18.8e%5d%5d%5d%18.8e%18.8e%18.8e"
        "%18.8e%5d%5d%5d%18.8e%18.8e%18.8e"
        "%18.8e%5d%5d%5d%18.8e%18.8e%18.8e\n",
        inumtot,hhhh,sqhhhh,
        res_1_max,i_res_1_max,j_res_1_max,k_res_1_max,sum_1_res,res_1_max_alone,vol_1_max_alone,
        res_2_max,i_res_2_max,j_res_2_max,k_res_2_max,sum_2_res,res_2_max_alone,vol_2_max_alone,
        res_3_max,i_res_3_max,j_res_3_max,k_res_3_max,sum_3_res,res_3_max_alone,vol_3_max_alone,
        res_4_max,i_res_4_max,j_res_4_max,k_res_4_max,sum_4_res,res_4_max_alone,vol_4_max_alone,
        res_5_max,i_res_5_max,j_res_5_max,k_res_5_max,sum_5_res,res_5_max_alone,vol_5_max_alone);
}


void exact_turb_res(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> vol, FortranArray3DRef<double> res,
    int& iexact_trunc, int& iexact_disc)
{
    double res_max=0.;
    double sum_res=0.;
    double res_max_alone=0., vol_max_alone=0.;
    int i_res_max=0, j_res_max=0, k_res_max=0;
    {
        int i=1;
        for (int j=3; j<=jdim-3; j++) {
            for (int k=3; k<=kdim-3; k++) {
                if (ccabs(res(j,k,i)) > res_max) {
                    res_max = ccabs(res(j,k,i));
                    i_res_max=i; j_res_max=j; k_res_max=k;
                    res_max_alone=res(j,k,i);
                    vol_max_alone=vol(j,k,i);
                }
                sum_res=sum_res+ccabs(res(j,k,i));
            }
        }
    }
    sum_res=sum_res/((jdim-5)*(kdim-5)*(idim-1));

    FILE* f302 = fortran_get_unit(302);
    fprintf(f302,"# turbulence residuals - truncation error if iexact_trunc=1\n");
    fprintf(f302,"variables=\"N\",\"1/N\",\"sqrt(1/N)\",\n");
    fprintf(f302,"\"res_max\",\"i_res_max\",\"j_res_max\",");
    fprintf(f302,"\"k_res_max\",\"res_L1norm\",\"res_max_alone\",");
    fprintf(f302,"\"vol_max_alone\"\n");
    int inumtot=(jdim-1)*(kdim-1);
    double hhhh=1./(double)((jdim-1)*(kdim-1));
    double sqhhhh=sqrt(hhhh);
    fprintf(f302,"%9d%18.8e%18.8e%18.8e%5d%5d%5d%18.8e%18.8e%18.8e\n",
        inumtot,hhhh,sqhhhh,
        res_max,i_res_max,j_res_max,k_res_max,sum_res,
        res_max_alone,vol_max_alone);
}

void exact_flow_force(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
    FortranArray3DRef<double> vol, FortranArray4DRef<double> res,
    FortranArray3DRef<double> smin,
    int& iexact_trunc, int& iexact_disc)
{
    FortranArray1D<double> fvec(6);
    int i_convert_q = 0;
    int i_forcing = 1;
    int i_gradient = 0;
    int nvar=6;
    for (int i=1; i<=idim-1; i++) {
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double distf=smin(j,k,i);
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) {
                    res(j,k,i,n)=res(j,k,i,n)-vol(j,k,i)*fvec(n);
                }
            }
        }
    }
}

void exact_turb_force(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
    FortranArray3DRef<double> vol, FortranArray3DRef<double> res,
    FortranArray3DRef<double> smin,
    int& iexact_trunc, int& iexact_disc)
{
    FortranArray1D<double> fvec(6);
    int i_convert_q = 0;
    int i_forcing = 1;
    int i_gradient = 0;
    int nvar=6;
    for (int i=1; i<=idim-1; i++) {
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double distf=smin(j,k,i);
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                res(j,k,i)=res(j,k,i)+fvec(6);
            }
        }
    }
}

void exact_flow_q(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
    FortranArray4DRef<double> q,
    int& iexact_trunc, int& iexact_disc)
{
    FortranArray1D<double> fvec(6);
    int i_convert_q = 0;
    int i_forcing = 0;
    int i_gradient = 0;
    int nvar=6;
    double distf=0.;
    for (int i=1; i<=idim-1; i++) {
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) {
                    q(j,k,i,n)=fvec(n);
                }
            }
        }
    }
}

void zero_resid_ring(int& jdim, int& kdim, int& idim,
    FortranArray4DRef<double> res,
    int& jj, int& kk, int& ii, int& numeq, int& numrows,
    int& iexact_trunc, int& iexact_disc)
{
    for (int n=1; n<=numeq; n++) {
        for (int i=1; i<=idim-1; i++) {
            for (int j=1; j<=numrows; j++) {
                for (int k=1; k<=kdim-1; k++) {
                    res(j,k,i,n)=0.0;
                }
            }
            for (int j=jdim-numrows; j<=jdim-1; j++) {
                for (int k=1; k<=kdim-1; k++) {
                    res(j,k,i,n)=0.0;
                }
            }
            for (int j=1; j<=jdim-1; j++) {
                for (int k=1; k<=numrows; k++) {
                    res(j,k,i,n)=0.0;
                }
            }
            for (int j=1; j<=jdim-1; j++) {
                for (int k=kdim-numrows; k<=kdim-1; k++) {
                    res(j,k,i,n)=0.0;
                }
            }
        }
    }
}



void exact_flow_q_ring(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
    FortranArray4DRef<double> q,
    int& iexact_trunc, int& iexact_disc)
{
    FortranArray1D<double> fvec(6);
    int i_convert_q = 0;
    int i_forcing = 0;
    int i_gradient = 0;
    int nvar=6;
    double distf=0.0;
    for (int i=1; i<=idim-1; i++) {
        // j=1,2
        for (int j=1; j<=2; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) q(j,k,i,n)=fvec(n);
            }
        }
        // j=jdim-2,jdim-1
        for (int j=jdim-2; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) q(j,k,i,n)=fvec(n);
            }
        }
        // k=1,2
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=2; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) q(j,k,i,n)=fvec(n);
            }
        }
        // k=kdim-2,kdim-1
        for (int j=1; j<=jdim-1; j++) {
            for (int k=kdim-2; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) q(j,k,i,n)=fvec(n);
            }
        }
    }
}

void exact_turb_q(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
    FortranArray3DRef<double> q, FortranArray3DRef<double> smin,
    FortranArray3DRef<double> vist3d,
    int& iexact_trunc, int& iexact_disc)
{
    FortranArray1D<double> fvec(6);
    int i_convert_q = 0;
    int i_forcing = 0;
    int i_gradient = 0;
    int nvar=6;
    for (int i=1; i<=idim-1; i++) {
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double distf = smin(j,k,i);
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                q(j,k,i)=fvec(6);
                vist3d(j,k,i)=xmut;
            }
        }
    }
}



void exact_turb_q_ring(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
    FortranArray3DRef<double> q, FortranArray3DRef<double> smin,
    FortranArray3DRef<double> vist3d,
    int& iexact_trunc, int& iexact_disc)
{
    FortranArray1D<double> fvec(6);
    int i_convert_q = 0;
    int i_forcing = 0;
    int i_gradient = 0;
    int nvar=6;
    for (int i=1; i<=idim-1; i++) {
        // j=1,2
        for (int j=1; j<=2; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double distf = smin(j,k,i);
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                q(j,k,i)=fvec(6); vist3d(j,k,i)=xmut;
            }
        }
        // j=jdim-2,jdim-1
        for (int j=jdim-2; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double distf = smin(j,k,i);
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                q(j,k,i)=fvec(6); vist3d(j,k,i)=xmut;
            }
        }
        // k=1,2
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=2; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double distf = smin(j,k,i);
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                q(j,k,i)=fvec(6); vist3d(j,k,i)=xmut;
            }
        }
        // k=kdim-2,kdim-1
        for (int j=1; j<=jdim-1; j++) {
            for (int k=kdim-2; k<=kdim-1; k++) {
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double distf = smin(j,k,i);
                double xx=0.125*(x(j  ,k  ,i  )+x(j+1,k  ,i  )+
                                 x(j  ,k+1,i  )+x(j  ,k  ,i+1)+
                                 x(j+1,k+1,i  )+x(j+1,k  ,i+1)+
                                 x(j  ,k+1,i+1)+x(j+1,k+1,i+1));
                double yy=0.125*(y(j  ,k  ,i  )+y(j+1,k  ,i  )+
                                 y(j  ,k+1,i  )+y(j  ,k  ,i+1)+
                                 y(j+1,k+1,i  )+y(j+1,k  ,i+1)+
                                 y(j  ,k+1,i+1)+y(j+1,k+1,i+1));
                double zz=0.125*(z(j  ,k  ,i  )+z(j+1,k  ,i  )+
                                 z(j  ,k+1,i  )+z(j  ,k  ,i+1)+
                                 z(j+1,k+1,i  )+z(j+1,k  ,i+1)+
                                 z(j  ,k+1,i+1)+z(j+1,k+1,i+1));
                double xmut;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                q(j,k,i)=fvec(6); vist3d(j,k,i)=xmut;
            }
        }
    }
}



// Helper macro for cell-center coordinate computation
#define CELL_CENTER_XYZ(xx,yy,zz,xarr,yarr,zarr,j,k,i) \
    xx=0.125*(xarr(j  ,k  ,i  )+xarr(j+1,k  ,i  )+ \
              xarr(j  ,k+1,i  )+xarr(j  ,k  ,i+1)+ \
              xarr(j+1,k+1,i  )+xarr(j+1,k  ,i+1)+ \
              xarr(j  ,k+1,i+1)+xarr(j+1,k+1,i+1)); \
    yy=0.125*(yarr(j  ,k  ,i  )+yarr(j+1,k  ,i  )+ \
              yarr(j  ,k+1,i  )+yarr(j  ,k  ,i+1)+ \
              yarr(j+1,k+1,i  )+yarr(j+1,k  ,i+1)+ \
              yarr(j  ,k+1,i+1)+yarr(j+1,k+1,i+1)); \
    zz=0.125*(zarr(j  ,k  ,i  )+zarr(j+1,k  ,i  )+ \
              zarr(j  ,k+1,i  )+zarr(j  ,k  ,i+1)+ \
              zarr(j+1,k+1,i  )+zarr(j+1,k  ,i+1)+ \
              zarr(j  ,k+1,i+1)+zarr(j+1,k+1,i+1));

void exact_flow_bc(int& jdim, int& kdim, int& idim,
    FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
    FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
    int& nface,
    FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
    FortranArray3DRef<double> smin,
    FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0,
    FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0,
    int& nummem,
    int& iexact_trunc, int& iexact_disc)
{
    int i2d    = cmn_twod.i2d;
    int ivmx   = cmn_maxiv.ivmx;
    int level  = cmn_mgrd.level;
    int lglobal= cmn_mgrd.lglobal;

    FortranArray1D<double> fvec(6);
    int i_convert_q = 0;
    int i_forcing = 0;
    int i_gradient = 0;
    int nvar=6;

    // QI0 lo
    if (i2d != 1) {
    if (nface==1) {
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                bci(j,k,1) = 0.0;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx1,yy1,zz1,xx2,yy2,zz2,xx,yy,zz,distf,xmut;
                {int i=1; double distf1=smin(j,k,i);
                CELL_CENTER_XYZ(xx1,yy1,zz1,x,y,z,j,k,i)
                i=2; double distf2=smin(j,k,i);
                CELL_CENTER_XYZ(xx2,yy2,zz2,x,y,z,j,k,i)
                xx=2.*xx1-xx2; yy=2.*yy1-yy2; zz=2.*zz1-zz2;
                distf=2.*distf1-distf2;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qi0(j,k,n,1)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) ti0(j,k,n,1)=fvec(6);
                if (ivmx >= 2) vi0(j,k,1,1)=xmut;
                xx=3.*xx1-2.*xx2; yy=3.*yy1-2.*yy2; zz=3.*zz1-2.*zz2;
                distf=3.*distf1-2.*distf2;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qi0(j,k,n,2)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) ti0(j,k,n,2)=fvec(6);
                if (ivmx >= 2) vi0(j,k,1,2)=xmut;}
            }
        }
    }
    // QI0 hi
    if (nface==2) {
        for (int j=1; j<=jdim-1; j++) {
            for (int k=1; k<=kdim-1; k++) {
                bci(j,k,2) = 0.0;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx1,yy1,zz1,xx2,yy2,zz2,xx,yy,zz,distf,xmut;
                {int i=idim-1; double distf1=smin(j,k,i);
                CELL_CENTER_XYZ(xx1,yy1,zz1,x,y,z,j,k,i)
                i=idim-2; double distf2=smin(j,k,i);
                CELL_CENTER_XYZ(xx2,yy2,zz2,x,y,z,j,k,i)
                xx=2.*xx1-xx2; yy=2.*yy1-yy2; zz=2.*zz1-zz2;
                distf=2.*distf1-distf2;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qi0(j,k,n,3)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) ti0(j,k,n,3)=fvec(6);
                if (ivmx >= 2) vi0(j,k,1,3)=xmut;
                xx=3.*xx1-2.*xx2; yy=3.*yy1-2.*yy2; zz=3.*zz1-2.*zz2;
                distf=3.*distf1-2.*distf2;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qi0(j,k,n,4)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) ti0(j,k,n,4)=fvec(6);
                if (ivmx >= 2) vi0(j,k,1,4)=xmut;}
            }
        }
    }
    } // end if (i2d != 1)


    // QJ0 lo
    if (nface == 3) {
        for (int i=1; i<=idim-1; i++) {
            for (int k=1; k<=kdim-1; k++) {
                bcj(k,i,1) = 0.0;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx1,yy1,zz1,xx2,yy2,zz2,xx,yy,zz,distf,xmut;
                {int j=1; double distf1=smin(j,k,i);
                CELL_CENTER_XYZ(xx1,yy1,zz1,x,y,z,j,k,i)
                j=2; double distf2=smin(j,k,i);
                CELL_CENTER_XYZ(xx2,yy2,zz2,x,y,z,j,k,i)
                xx=2.*xx1-xx2; yy=2.*yy1-yy2; zz=2.*zz1-zz2;
                distf=2.*distf1-distf2;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qj0(k,i,n,1)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) tj0(k,i,n,1)=fvec(6);
                if (ivmx >= 2) vj0(k,i,1,1)=xmut;
                xx=3.*xx1-2.*xx2; yy=3.*yy1-2.*yy2; zz=3.*zz1-2.*zz2;
                distf=3.*distf1-2.*distf2;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qj0(k,i,n,2)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) tj0(k,i,n,2)=fvec(6);
                if (ivmx >= 2) vj0(k,i,1,2)=xmut;}
            }
        }
    }
    // QJ0 hi
    if (nface == 4) {
        for (int i=1; i<=idim-1; i++) {
            for (int k=1; k<=kdim-1; k++) {
                bcj(k,i,2) = 0.0;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx1,yy1,zz1,xx2,yy2,zz2,xx,yy,zz,distf,xmut;
                {int j=jdim-1; double distf1=smin(j,k,i);
                CELL_CENTER_XYZ(xx1,yy1,zz1,x,y,z,j,k,i)
                j=jdim-2; double distf2=smin(j,k,i);
                CELL_CENTER_XYZ(xx2,yy2,zz2,x,y,z,j,k,i)
                xx=2.*xx1-xx2; yy=2.*yy1-yy2; zz=2.*zz1-zz2;
                distf=2.*distf1-distf2;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qj0(k,i,n,3)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) tj0(k,i,n,3)=fvec(6);
                if (ivmx >= 2) vj0(k,i,1,3)=xmut;
                xx=3.*xx1-2.*xx2; yy=3.*yy1-2.*yy2; zz=3.*zz1-2.*zz2;
                distf=3.*distf1-2.*distf2;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qj0(k,i,n,4)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) tj0(k,i,n,4)=fvec(6);
                if (ivmx >= 2) vj0(k,i,1,4)=xmut;}
            }
        }
    }
    // QK0 lo
    if (nface==5) {
        for (int i=1; i<=idim-1; i++) {
            for (int j=1; j<=jdim-1; j++) {
                bck(j,i,1) = 0.0;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx1,yy1,zz1,xx2,yy2,zz2,xx,yy,zz,distf,xmut;
                {int k=1; double distf1=smin(j,k,i);
                CELL_CENTER_XYZ(xx1,yy1,zz1,x,y,z,j,k,i)
                k=2; double distf2=smin(j,k,i);
                CELL_CENTER_XYZ(xx2,yy2,zz2,x,y,z,j,k,i)
                xx=2.*xx1-xx2; yy=2.*yy1-yy2; zz=2.*zz1-zz2;
                distf=2.*distf1-distf2;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qk0(j,i,n,1)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) tk0(j,i,n,1)=fvec(6);
                if (ivmx >= 2) vk0(j,i,1,1)=xmut;
                xx=3.*xx1-2.*xx2; yy=3.*yy1-2.*yy2; zz=3.*zz1-2.*zz2;
                distf=3.*distf1-2.*distf2;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qk0(j,i,n,2)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) tk0(j,i,n,2)=fvec(6);
                if (ivmx >= 2) vk0(j,i,1,2)=xmut;}
            }
        }
    }
    // QK0 hi
    if (nface==6) {
        for (int i=1; i<=idim-1; i++) {
            for (int j=1; j<=jdim-1; j++) {
                bck(j,i,2) = 0.0;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                double xx1,yy1,zz1,xx2,yy2,zz2,xx,yy,zz,distf,xmut;
                {int k=kdim-1; double distf1=smin(j,k,i);
                CELL_CENTER_XYZ(xx1,yy1,zz1,x,y,z,j,k,i)
                k=kdim-2; double distf2=smin(j,k,i);
                CELL_CENTER_XYZ(xx2,yy2,zz2,x,y,z,j,k,i)
                xx=2.*xx1-xx2; yy=2.*yy1-yy2; zz=2.*zz1-zz2;
                distf=2.*distf1-distf2;
                FortranArray1DRef<double> fref = fvec.ref();
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qk0(j,i,n,3)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) tk0(j,i,n,3)=fvec(6);
                if (ivmx >= 2) vk0(j,i,1,3)=xmut;
                xx=3.*xx1-2.*xx2; yy=3.*yy1-2.*yy2; zz=3.*zz1-2.*zz2;
                distf=3.*distf1-2.*distf2;
                for (int n=1; n<=nvar; n++) fvec(n)=0.;
                analytic_compressible(xx,yy,zz,nvar,fref,
                    i_convert_q,i_forcing,i_gradient,distf,xmut,
                    iexact_trunc,iexact_disc);
                for (int n=1; n<=5; n++) qk0(j,i,n,4)=fvec(n);
                if (level >= lglobal && ivmx >= 4)
                    for (int n=1; n<=nummem; n++) tk0(j,i,n,4)=fvec(6);
                if (ivmx >= 2) vk0(j,i,1,4)=xmut;}
            }
        }
    }
}

#undef CELL_CENTER_XYZ



void analytic_compressible(double& x, double& y, double& z, int& neq,
    FortranArray1DRef<double> q, int& i_convert_q, int& i_forcing,
    int& i_gradient, double& distf, double& xmut,
    int& iexact_trunc, int& iexact_disc)
{
    double xmach  = (double)cmn_info.xmach;
    double gamma  = (double)cmn_fluid.gamma;
    double gm1    = (double)cmn_fluid.gm1;
    double pr     = (double)cmn_fluid2.pr;
    double prt    = (double)cmn_fluid2.prt;
    double reue   = (double)cmn_reyue.reue;
    int    ivmx   = cmn_maxiv.ivmx;
    int    iexp   = cmn_zero.iexp;

    double xminn = std::pow(10.0, (double)(-iexp+1));

    // Initialize q
    for (int n=1; n<=neq; n++) q(n)=0.;

    // Spalart Constants
    double sig  = 0.66667;
    double cv1  = 7.1;
    double ct3  = 1.2;
    double ct4  = 0.5;
    double vkar = 0.41;
    double cb1  = 0.1355;
    double cb2  = 0.622;
    double cw1  = cb1/(vkar*vkar)+(1.+cb2)/sig;
    double cw2  = 0.3;
    double cw3  = 2.0;

    // Other variables
    double xmre   = xmach / reue;
    double xmre_s = xmre / sig;
    double cgp    = xmre /(gm1*pr);
    double cgpt   = xmre /(gm1*prt);

    // Initialize
    double u,ux,uy,uz,uxx,uyy,uzz,uxy,uxz,uyz;
    double v,vx,vy,vz,vxx,vyy,vzz,vxy,vxz,vyz;
    double w,wx,wy,wz,wxx,wyy,wzz,wxy,wxz,wyz;
    double p,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz;
    double t,tx,ty,tz,txx,tyy,tzz,txy,txz,tyz;
    // vt and derivatives declared at function scope for use in final section
    double vt=0.,vtx=0.,vty=0.,vtz=0.,vtxx=0.,vtyy=0.,vtzz=0.,vtxy=0.,vtxz=0.,vtyz=0.;
    exact_zero(u,ux,uy,uz,uxx,uyy,uzz,uxy,uxz,uyz);
    exact_zero(v,vx,vy,vz,vxx,vyy,vzz,vxy,vxz,vyz);
    exact_zero(w,wx,wy,wz,wxx,wyy,wzz,wxy,wxz,wyz);
    exact_zero(p,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz);
    exact_zero(t,tx,ty,tz,txx,tyy,tzz,txy,txz,tyz);

    // Choose manufactured solution
    int ims_sol = std::max(iexact_trunc, iexact_disc);
    if (ims_sol == 999) {
        double u_scale=xmach;
        double w_scale=xmach;
        double p_scale=1.0/gamma;
        FortranArray1D<double> u_polyf(13), w_polyf(13), p_polyf(13);
        u_polyf(1)=1.; u_polyf(2)=1.;
        for (int n=3; n<=13; n++) u_polyf(n)=0.;
        w_polyf(1)=1.; w_polyf(2)=1.;
        for (int n=3; n<=13; n++) w_polyf(n)=0.;
        p_polyf(1)=1.; p_polyf(2)=1.;
        for (int n=3; n<=13; n++) p_polyf(n)=0.;
        FortranArray1DRef<double> upref=u_polyf.ref();
        FortranArray1DRef<double> wpref=w_polyf.ref();
        FortranArray1DRef<double> ppref=p_polyf.ref();
        exact_polyfg(u_scale,upref,x,y,z,u,ux,uy,uz,uxx,uyy,uzz,uxy,uxz,uyz);
        exact_polyfg(w_scale,wpref,x,y,z,w,wx,wy,wz,wxx,wyy,wzz,wxy,wxz,wyz);
        exact_polyfg(p_scale,ppref,x,y,z,p,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz);
    } else if (ims_sol==1 || ims_sol==2 || ims_sol==4) {
        double xlisbon_p_scale = 1.0/gamma;
        exact_lisbon_ms2_u(xmach, x, z, u, ux, uz, uxx, uzz, uxz);
        exact_lisbon_ms2_v(xmach, x, z, w, wx, wz, wxx, wzz, wxz);
        exact_lisbon_ms2_p(xlisbon_p_scale, x, z, p, px, pz);
        p = p + xlisbon_p_scale;
    } else {
        fprintf(fortran_get_unit(6)," MFG soln not chosen\n");
    }

    // Enforce constant adiabatic wall temperature (2D).
    total_temperature_2d(t,tx,tz,txx,tzz,txz,
                         u,ux,uz,uxx,uzz,uxz,
                         w,wx,wz,wxx,wzz,wxz);

    // Density
    double rho  = gamma*( p/t );
    double rhox = gamma*( px/t - p*tx/(t*t) );
    double rhoy = gamma*( py/t - p*ty/(t*t) );
    double rhoz = gamma*( pz/t - p*tz/(t*t) );

    // Molecular viscosity
    double xmu    = viscosity_law(t);
    double dxmudt = dviscosity_law(t);
    double xmux  = dxmudt*tx;
    double xmuy  = dxmudt*ty;
    double xmuz  = dxmudt*tz;

    // Kinematic viscosity
    double xnu   = xmu/rho;
    double xnux  = (xmux - xnu*rhox)/rho;
    double xnuy  = (xmuy - xnu*rhoy)/rho;
    double xnuz  = (xmuz - xnu*rhoz)/rho;

    // Continuity equation
    double d    = ux + vy + wz;
    FortranArray1D<double> f(neq);
    for (int n=1; n<=neq; n++) f(n)=0.;
    f(1) = rho*d + u*rhox + v*rhoy + w*rhoz;

    // Inviscid momentum equation
    f(2) = rho*( u*ux + v*uy + w*uz ) + px + u*f(1);
    f(3) = rho*( u*vx + v*vy + w*vz ) + py + v*f(1);
    f(4) = rho*( u*wx + v*wy + w*wz ) + pz + w*f(1);

    // Inviscid energy equation
    double h  =  t/gm1 + 0.5*(u*u + v*v + w*w);
    double hx = tx/gm1 + u*ux + v*vx + w*wx;
    double hy = ty/gm1 + u*uy + v*vy + w*wy;
    double hz = tz/gm1 + u*uz + v*vz + w*wz;
    double rhoh  = rho*h;
    double rhohx = rhox*h + rho*hx;
    double rhohy = rhoy*h + rho*hy;
    double rhohz = rhoz*h + rho*hz;
    f(5) = u*rhohx + v*rhohy + w*rhohz + rhoh*d;

    // s2 = twice mean-strain-rate tensor
    double s2_xx = ux + ux;
    double s2_xy = vx + uy;
    double s2_xz = wx + uz;
    double s2_yy = vy + vy;
    double s2_yz = wy + vz;
    double s2_zz = wz + wz;

    // Laminar momentum equations
    double xla_u = uxx + uyy + uzz;
    double xla_v = vxx + vyy + vzz;
    double xla_w = wxx + wyy + wzz;
    double d3   = 1./3.*d;
    double d3x  = 1./3.*( uxx + vxy + wxz );
    double d3y  = 1./3.*( uxy + vyy + wyz );
    double d3z  = 1./3.*( uxz + vyz + wzz );

    // Momentum equations - first term
    double vf2 = - xmre*xmu*( xla_u + d3x );
    double vf3 = - xmre*xmu*( xla_v + d3y );
    double vf4 = - xmre*xmu*( xla_w + d3z );

    // Momentum equations - second term
    vf2 = vf2 - xmre*( xmux*s2_xx + xmuy*s2_xy + xmuz*s2_xz );
    vf3 = vf3 - xmre*( xmux*s2_xy + xmuy*s2_yy + xmuz*s2_yz );
    vf4 = vf4 - xmre*( xmux*s2_xz + xmuy*s2_yz + xmuz*s2_zz );
    vf2 = vf2 + xmre*2.0*d3*xmux;
    vf3 = vf3 + xmre*2.0*d3*xmuy;
    vf4 = vf4 + xmre*2.0*d3*xmuz;
    f(2) = f(2) + vf2;
    f(3) = f(3) + vf3;
    f(4) = f(4) + vf4;

    // Laminar energy equation
    double xla_t = txx + tyy + tzz;
    double vf5  = u*vf2 + v*vf3 + w*vf4
       - cgp*( xmu*xla_t + xmux*tx + xmuy*ty + xmuz*tz )
       - xmre*xmu*( ux*s2_xx + uy*s2_xy + uz*s2_xz
                  + vx*s2_xy + vy*s2_yy + vz*s2_yz
                  + wx*s2_xz + wy*s2_yz + wz*s2_zz - 2.0*d3*d );
    f(5) = f(5) + vf5;


    // Spalart turbulence model equation
    if (ivmx == 5) {
        exact_zero(vt,vtx,vty,vtz,vtxx,vtyy,vtzz,vtxy,vtxz,vtyz);
        double term = 1.0;
        if (ims_sol == 1) {
            exact_lisbon_ms1_vt(term, x, z, vt, vtx, vtz, vtxx, vtzz);
        } else if (ims_sol == 2) {
            exact_lisbon_ms2_vt(term, x, z, vt, vtx, vtz, vtxx, vtzz);
        } else if (ims_sol == 4) {
            exact_lisbon_ms4_vt(term, x, z, vt, vtx, vtz, vtxx, vtzz);
        } else {
            fprintf(fortran_get_unit(6)," Error... must have ims_sol=1,2, or 4 with SA\n");
        }
        double chi  =   vt/xnu;
        double chix = ( vtx - chi*xnux )/xnu;
        double chiy = ( vty - chi*xnuy )/xnu;
        double chiz = ( vtz - chi*xnuz )/xnu;
        double fv1  = chi*chi*chi / (chi*chi*chi+cv1*cv1*cv1);
        double fv1x = 3.0*chix*( 1.0 - fv1 )*chi*chi/(chi*chi*chi+cv1*cv1*cv1);
        double fv1y = 3.0*chiy*( 1.0 - fv1 )*chi*chi/(chi*chi*chi+cv1*cv1*cv1);
        double fv1z = 3.0*chiz*( 1.0 - fv1 )*chi*chi/(chi*chi*chi+cv1*cv1*cv1);
        xmut  =  rho*fv1*vt;
        double xmutx = rhox*fv1*vt + rho*( fv1x*vt + fv1*vtx );
        double xmuty = rhoy*fv1*vt + rho*( fv1y*vt + fv1*vty );
        double xmutz = rhoz*fv1*vt + rho*( fv1z*vt + fv1*vtz );
        // Momentum equations - first term
        double tf2 = - xmre*xmut*( xla_u + d3x );
        double tf3 = - xmre*xmut*( xla_v + d3y );
        double tf4 = - xmre*xmut*( xla_w + d3z );
        // Momentum equations - second term
        tf2 = tf2 - xmre*( xmutx*s2_xx + xmuty*s2_xy + xmutz*s2_xz );
        tf3 = tf3 - xmre*( xmutx*s2_xy + xmuty*s2_yy + xmutz*s2_yz );
        tf4 = tf4 - xmre*( xmutx*s2_xz + xmuty*s2_yz + xmutz*s2_zz );
        tf2 = tf2 + xmre*2.0*d3*xmutx;
        tf3 = tf3 + xmre*2.0*d3*xmuty;
        tf4 = tf4 + xmre*2.0*d3*xmutz;
        f(2) = f(2) + tf2;
        f(3) = f(3) + tf3;
        f(4) = f(4) + tf4;
        // Energy equation
        f(5) = f(5) + u*tf2 + v*tf3 + w*tf4
           - cgpt*( xmut*xla_t + xmutx*tx + xmuty*ty + xmutz*tz )
           - xmre*xmut*( ux*s2_xx + uy*s2_xy + uz*s2_xz
                       + vx*s2_xy + vy*s2_yy + vz*s2_yz
                       + wx*s2_xz + wy*s2_yz + wz*s2_zz - 2.0*d3*d );

        if (distf != 0.) {
            double t_conv=1.0, t_diff1=1.0, t_diff2=1.0, t_diff3=1.0;
            double t_prod=1.0, t_dest=1.0;
            double distance = distf;
            double fv2 = 1.0 - chi/(1.0+chi*fv1);
            double ft2 = ct3 * std::exp(-ct4*chi*chi);
            double vkar2 = vkar * vkar;
            double bot = (vkar*distance)*(vkar*distance);
            double s   = std::sqrt((wy-vz)*(wy-vz)+(uz-wx)*(uz-wx)+(vx-uy)*(vx-uy));
            double sw  = s + xmre*vt*fv2/bot;
            sw  = std::max(sw,xminn);
            double rr  = xmre*vt/(bot*sw);
            if(rr > 10.0) rr = 10.0;
            double gg  = rr + cw2*(std::pow(rr,6.0)-rr);
            gg = std::max(gg,xminn);
            double fw  = gg * std::pow((1.0+std::pow(cw3,6.0))/(std::pow(gg,6.0)+std::pow(cw3,6.0)),1./6.);
            double xla_vt = vtxx + vtyy + vtzz;
            double f_conv = u*vtx + v*vty + w*vtz;
            double f_prod = cb1*(1.0-ft2)*s*vt;
            double f_dest = xmre*( cb1*((1.-ft2)*fv2+ft2)/vkar2 -
                            cw1*fw )*(vt/distance)*(vt/distance);
            double cb21 = t_diff2 + t_diff3*cb2;
            double f_diff = xmre_s*( ( t_diff1*xnu + t_diff2*vt )*xla_vt
                                   + ( t_diff1*xnux + cb21*vtx )*vtx
                                   + ( t_diff1*xnuy + cb21*vty )*vty
                                   + ( t_diff1*xnuz + cb21*vtz )*vtz );
            q(6) = vt;
            f(6) = t_conv*f_conv
                 - ( t_prod*f_prod + t_dest*f_dest )
                 - f_diff;
        } else {
            q(6) = vt;
            f(6) = 0.0;
        }
    }

    // Return requested result in q
    q(1) = rho;
    q(2) = u;
    q(3) = v;
    q(4) = w;
    q(5) = p;
    if ( ivmx == 5 ) q(6) = vt;
    if (i_convert_q == 1) {
        q(1) = rho;
        q(2) = rho*u;
        q(3) = rho*v;
        q(4) = rho*w;
        q(5) = p/gm1 + 0.5*rho*(u*u + v*v + w*w);
    } else if (i_forcing == 1) {
        q(1) = f(1);
        q(2) = f(2);
        q(3) = f(3);
        q(4) = f(4);
        q(5) = f(5);
        if ( ivmx == 5 ) q(6) = f(6);
    } else if (i_gradient == 1) {
        q(1) = rhox;
        q(2) = ux;
        q(3) = vx;
        q(4) = wx;
        q(5) = px;
        if ( ivmx == 5 ) q(6) = vtx;
    } else if (i_gradient == 2) {
        q(1) = rhoy;
        q(2) = uy;
        q(3) = vy;
        q(4) = wy;
        q(5) = py;
        if ( ivmx == 5 ) q(6) = vty;
    } else if (i_gradient == 3) {
        q(1) = rhoz;
        q(2) = uz;
        q(3) = vz;
        q(4) = wz;
        q(5) = pz;
        if ( ivmx == 5 ) q(6) = vtz;
    }
}



void exact_zero(double& q, double& qx, double& qy, double& qz,
    double& qxx, double& qyy, double& qzz, double& qxy, double& qxz, double& qyz)
{
    q = 0.0; qx = 0.0; qy = 0.0; qz = 0.0;
    qxx = 0.0; qyy = 0.0; qzz = 0.0;
    qxy = 0.0; qxz = 0.0; qyz = 0.0;
}

void total_temperature_2d(double& t, double& tx, double& ty,
    double& txx, double& tyy, double& txy,
    double& u, double& ux, double& uy, double& uxx, double& uyy, double& uxy,
    double& v, double& vx, double& vy, double& vxx, double& vyy, double& vxy)
{
    double xmach = (double)cmn_info.xmach;
    double gm1   = (double)cmn_fluid.gm1;
    // correction for consistency with CFL3D's TWTYPE=-1 for BC2004:
    double d = gm1;
    double t0 = 1.0 + 0.5*d*xmach*xmach;
    t = t0 - 0.5*d*( u*u + v*v );
    tx = -d*( u*ux + v*vx );
    ty = -d*( u*uy + v*vy );
    txx = -d*( u*uxx + ux*ux + v*vxx + vx*vx );
    tyy = -d*( u*uyy + uy*uy + v*vyy + vy*vy );
    txy = -d*( u*uxy + uy*ux + v*vxy + vy*vx );
}

double viscosity_law(double& t)
{
    double tinf = (double)cmn_reyue.tinf;
    double suth_const=198.6;
    double cstar=suth_const/tinf;
    return (1.0+cstar)/(t+cstar)*std::pow(t,1.5);
}

double dviscosity_law(double& t)
{
    double tinf = (double)cmn_reyue.tinf;
    double suth_const=198.6;
    double cstar=suth_const/tinf;
    double tstar  = t+cstar;
    double tstari = 1.0/tstar;
    double visc = (1.0+cstar)*tstari*std::pow(t,1.5);
    return visc*( 0.5 + 1.5*cstar/t )*tstari;
}

void exact_lisbon_ms2_forcing(double& beta, double& x, double& y,
    FortranArray1DRef<double> f)
{
    double emu=1.e-6;
    f(1) = beta*( dudxms(x,y) + dvdyms(x,y) );
    f(2) = smxsams(x,y);
    f(3) = 0.0;
    f(4) = smysams(x,y);
    f(5) = ssams(x,y)/emu;
}

void exact_lisbon_ms2_u(double& scale, double& x, double& y,
    double& q, double& qx, double& qy, double& qxx, double& qyy, double& qxy)
{
    q   = q   + scale*ums(x,y);
    qx  = qx  + scale*dudxms(x,y);
    qy  = qy  + scale*dudyms(x,y);
    qxx = qxx + scale*dudx2ms(x,y);
    qyy = qyy + scale*dudy2ms(x,y);
    qxy = qxy + scale*dudxyms(x,y);
}

void exact_lisbon_ms2_v(double& scale, double& x, double& y,
    double& q, double& qx, double& qy, double& qxx, double& qyy, double& qxy)
{
    q   = q   + scale*vms(x,y);
    qx  = qx  + scale*dvdxms(x,y);
    qy  = qy  + scale*dvdyms(x,y);
    qxx = qxx + scale*dvdx2ms(x,y);
    qyy = qyy + scale*dvdy2ms(x,y);
    qxy = qxy + scale*dvdxyms(x,y);
}

void exact_lisbon_ms2_p(double& scale, double& x, double& y,
    double& q, double& qx, double& qy)
{
    q  = q  + scale*pms(x,y);
    qx = qx + scale*dpdxms(x,y);
    qy = qy + scale*dpdyms(x,y);
}

void exact_lisbon_ms1_vt(double& scale, double& x, double& y,
    double& q, double& qx, double& qy, double& qxx, double& qyy)
{
    double emu=1.e-6;
    q   = q   + scale*em1ms(x,y)/emu;
    qx  = qx  + scale*dem1dxms(x,y)/emu;
    qy  = qy  + scale*dem1dyms(x,y)/emu;
    qxx = qxx + scale*dem1dx2ms(x,y)/emu;
    qyy = qyy + scale*dem1dy2ms(x,y)/emu;
}

void exact_lisbon_ms2_vt(double& scale, double& x, double& y,
    double& q, double& qx, double& qy, double& qxx, double& qyy)
{
    double emu=1.e-6;
    q   = q   + scale*em2ms(x,y)/emu;
    qx  = qx  + scale*dem2dxms(x,y)/emu;
    qy  = qy  + scale*dem2dyms(x,y)/emu;
    qxx = qxx + scale*dem2dx2ms(x,y)/emu;
    qyy = qyy + scale*dem2dy2ms(x,y)/emu;
}

void exact_lisbon_ms4_vt(double& scale, double& x, double& y,
    double& q, double& qx, double& qy, double& qxx, double& qyy)
{
    double emu=1.e-6;
    q   = q   + scale*em4ms(x,y)/emu;
    qx  = qx  + scale*dem4dxms(x,y)/emu;
    qy  = qy  + scale*dem4dyms(x,y)/emu;
    qxx = qxx + scale*dem4dx2ms(x,y)/emu;
    qyy = qyy + scale*dem4dy2ms(x,y)/emu;
}

double ums(double& x, double& y)
{
    double sigma=4.0;
    double eta=sigma*y/x;
    return user_derf(eta);
}

double dudxms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    return -u1*sigma*y*2.0/(x*x*std::sqrt(pi))*std::exp(-eta*eta);
}

double dudyms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    return u1*sigma*2.0/(x*std::sqrt(pi))*std::exp(-eta*eta);
}

double dudx2ms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double t1=4.0*u1/std::sqrt(pi);
    double x2=x*x;
    return t1*(eta/x2)*std::exp(-eta*eta)*(1.0-eta*eta);
}

double dudy2ms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double t1=4.0*u1/std::sqrt(pi);
    double x2=sigma/x;
    return -t1*x2*x2*eta*std::exp(-eta*eta);
}

double dudxyms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double t1=2.0*u1/std::sqrt(pi)*sigma/x/x;
    return t1*std::exp(-eta*eta)*(2.0*eta*eta-1.0);
}

double dudx2yms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double eta2=eta*eta;
    double t1=4.0*u1/std::sqrt(pi)*sigma/x/x/x;
    return t1*std::exp(-eta2)*(1.0+2.0*eta2*(eta2-2.5));
}

double dudy3ms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double x2=sigma/x;
    double t1=4.0*u1/std::sqrt(pi)*x2*x2*x2;
    return t1*std::exp(-eta*eta)*(2.0*eta*eta-1.0);
}



double vms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double t1=u1/(sigma*std::sqrt(pi));
    return t1*(1.0-std::exp(-eta*eta));
}

double dvdxms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    return -u1*sigma*y*y*2.0/(x*x*x*std::sqrt(pi))*std::exp(-eta*eta);
}

double dvdyms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    return u1*sigma*y*2.0/(x*x*std::sqrt(pi))*std::exp(-eta*eta);
}

double dvdx2ms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double t1=2.0*u1/std::sqrt(pi);
    double x2=sigma*y*y/(x*x*x*x);
    return t1*x2*std::exp(-eta*eta)*(3.0-2.0*eta*eta);
}

double dvdy2ms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double t1=2.0*u1/std::sqrt(pi);
    double x2=sigma/(x*x);
    return t1*x2*std::exp(-eta*eta)*(1.0-2.0*eta*eta);
}

double dvdxyms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double t1=4.0*u1/std::sqrt(pi)*sigma*y/x/x/x;
    return t1*std::exp(-eta*eta)*(eta*eta-1.0);
}

double dvdxy2ms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double eta2=eta*eta;
    double t1=4.0*u1/std::sqrt(pi)*sigma/x/x/x;
    return t1*std::exp(-eta2)*(2.0*eta2*(2.5-eta2)-1.0);
}

double dvdx3ms(double& x, double& y)
{
    double sigma=4.0;
    double pi=std::acos(-1.);
    double u1=1.0;
    double eta=sigma*y/x;
    double eta2=eta*eta;
    double t1=4.0*u1/std::sqrt(pi);
    double x2=eta*y/(x*x*x*x);
    return t1*x2*std::exp(-eta2)*(2.0*eta2*(4.5-eta2)-6.0);
}

double pms(double& x, double& y)
{
    return 0.5*std::log(2.0*x-x*x+0.25)*std::log(4.0*y*y*y-3.0*y*y+1.25);
}

double dpdxms(double& x, double& y)
{
    return 0.5*(2.0-2.0*x)/(2.0*x-x*x+0.25)
           *std::log(4.0*y*y*y-3.0*y*y+1.25);
}

double dpdyms(double& x, double& y)
{
    return 0.5*(12.0*y*y-6.0*y)/(4.0*y*y*y-3.0*y*y+1.25)
           *std::log(2.0*x-x*x+0.25);
}

double strainms(double& x, double& y)
{
    double u1x=dudxms(x,y);
    double u1y=dudyms(x,y);
    double u2x=dvdxms(x,y);
    double u2y=dvdyms(x,y);
    double val=2.0*(u1x*u1x+u2y*u2y)+(u1y+u2x)*(u1y+u2x);
    return std::sqrt(val);
}

double vortms(double& x, double& y)
{
    double u1y=dudyms(x,y);
    double u2x=dvdxms(x,y);
    double val=(u1y-u2x)*(u1y-u2x);
    return std::sqrt(val);
}

double dvodxms(double& x, double& y)
{
    double u1y=dudyms(x,y);
    double u2x=dvdxms(x,y);
    double vort=std::copysign(1.0, u1y-u2x);
    return vort*(dudxyms(x,y)-dvdx2ms(x,y));
}

double dvodyms(double& x, double& y)
{
    double u1y=dudyms(x,y);
    double u2x=dvdxms(x,y);
    double vort=std::copysign(1.0, u1y-u2x);
    return vort*(dudy2ms(x,y)-dvdxyms(x,y));
}



double dvodx2ms(double& x, double& y)
{
    double u1y=dudyms(x,y);
    double u2x=dvdxms(x,y);
    double vort=std::copysign(1.0, u1y-u2x);
    return vort*(dudx2yms(x,y)-dvdx3ms(x,y));
}

double dvody2ms(double& x, double& y)
{
    double u1y=dudyms(x,y);
    double u2x=dvdxms(x,y);
    double vort=std::copysign(1.0, u1y-u2x);
    return vort*(dudy3ms(x,y)-dvdxy2ms(x,y));
}

double smxsams(double& x, double& y)
{
    double emu=1.e-6;
    double u1=ums(x,y);
    double u2=vms(x,y);
    double em=eddysams(x,y);
    double cv1=u1*dudxms(x,y);
    double cv2=u2*dudyms(x,y);
    double cvx=cv1+cv2;
    double df1=-(emu+em)*(dudx2ms(x,y)+dudy2ms(x,y));
    double df2=-2.0*desadxms(x,y)*dudxms(x,y);
    double df3=-desadyms(x,y)*(dudyms(x,y)+dvdxms(x,y));
    double dfx=df1+df2+df3;
    double dpx=dpdxms(x,y);
    return cvx+dfx+dpx;
}

double smysams(double& x, double& y)
{
    double emu=1.e-6;
    double u1=ums(x,y);
    double u2=vms(x,y);
    double em=eddysams(x,y);
    double cv1=u1*dvdxms(x,y);
    double cv2=u2*dvdyms(x,y);
    double cvy=cv1+cv2;
    double df1=-(emu+em)*(dvdx2ms(x,y)+dvdy2ms(x,y));
    double df2=-2.0*desadyms(x,y)*dvdyms(x,y);
    double df3=-desadxms(x,y)*(dudyms(x,y)+dvdxms(x,y));
    double dfy=df1+df2+df3;
    double dpy=dpdyms(x,y);
    return cvy+dfy+dpy;
}

double smxmtms(double& x, double& y)
{
    double emu=1.e-6;
    double u1=ums(x,y);
    double u2=vms(x,y);
    double em=eddymtms(x,y);
    double cv1=u1*dudxms(x,y);
    double cv2=u2*dudyms(x,y);
    double cvx=cv1+cv2;
    double df1=-(emu+em)*(dudx2ms(x,y)+dudy2ms(x,y));
    double df2=-2.0*demtdxms(x,y)*dudxms(x,y);
    double df3=-demtdyms(x,y)*(dudyms(x,y)+dvdxms(x,y));
    double dfx=df1+df2+df3;
    double dpx=dpdxms(x,y);
    return cvx+dfx+dpx;
}

double smymtms(double& x, double& y)
{
    double emu=1.e-6;
    double u1=ums(x,y);
    double u2=vms(x,y);
    double em=eddymtms(x,y);
    double cv1=u1*dvdxms(x,y);
    double cv2=u2*dvdyms(x,y);
    double cvy=cv1+cv2;
    double df1=-(emu+em)*(dvdx2ms(x,y)+dvdy2ms(x,y));
    double df2=-2.0*demtdyms(x,y)*dvdyms(x,y);
    double df3=-demtdxms(x,y)*(dudyms(x,y)+dvdxms(x,y));
    double dfy=df1+df2+df3;
    double dpy=dpdyms(x,y);
    return cvy+dfy+dpy;
}

double em1ms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double t0=std::sqrt(2.0)*std::exp(0.5)*emmax;
    return t0*eta*std::exp(-eta*eta);
}

double em2ms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double t1=std::exp(1.0)*emmax;
    return t1*eta*eta*std::exp(-eta*eta);
}

double em4ms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double t1=0.25*std::exp(2.)*emmax;
    return t1*eta*eta*eta*eta*std::exp(-eta*eta);
}

double dem1dxms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double t0=std::sqrt(2.0)*std::exp(0.5)*emmax;
    double eta=sigeta*y/x;
    double eta_x = -eta/x;
    return t0*eta_x*std::exp(-eta*eta)*( 1.0 - 2.0*eta*eta );
}

double dem2dxms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double t1=2.0*std::exp(1.0)*emmax/x;
    return t1*eta*eta*std::exp(-eta*eta)*(eta*eta-1.0);
}

double dem4dxms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double t1=0.5*std::exp(2.)*emmax/x;
    return t1*eta*eta*eta*eta*std::exp(-eta*eta)*(eta*eta-2.0);
}

double dem1dyms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double t0=std::sqrt(2.0)*std::exp(0.5)*emmax;
    double eta   = sigeta*y/x;
    double eta_y = sigeta/x;
    return t0*eta_y*std::exp(-eta*eta)*( 1.0 - 2.0*eta*eta );
}

double dem2dyms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double t1=2.0*std::exp(1.0)*emmax*sigeta/x;
    return t1*eta*std::exp(-eta*eta)*(1.0-eta*eta);
}

double dem4dyms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double t1=0.5*std::exp(2.)*emmax*sigeta/x;
    return t1*eta*eta*eta*std::exp(-eta*eta)*(2.0-eta*eta);
}



double dem1dx2ms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double t0=std::sqrt(2.0)*std::exp(0.5)*emmax;
    double eta    = sigeta*y/x;
    double eta_x  = -eta/x;
    double eta_xx = -2.0*eta_x/x;
    return t0*std::exp(-eta*eta)*(        eta_xx*(  1.0 - 2.0*eta*eta )
                           + eta*eta_x*eta_x*( -6.0 + 4.0*eta*eta ) );
}

double dem2dx2ms(double& x, double& y)
{
    double sigeta=1.e1;
    double eta=sigeta*y/x;
    double et2=eta*eta;
    double t1=2.0*em2ms(x,y)/x/x;
    return t1*(et2*(2.0*et2-7.0)+3.0);
}

double dem4dx2ms(double& x, double& y)
{
    double sigeta=1.e1;
    double eta=sigeta*y/x;
    double et2=eta*eta;
    double t1=2.0*em4ms(x,y)/x/x;
    return t1*(et2*(2.0*et2-11.0)+10.0);
}

double dem1dy2ms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double t0=std::sqrt(2.0)*std::exp(0.5)*emmax;
    double eta    = sigeta*y/x;
    double eta_y  = sigeta/x;
    double eta_yy = 0.0;
    return t0*std::exp(-eta*eta)*(        eta_yy*(  1.0 - 2.0*eta*eta )
                           + eta*eta_y*eta_y*( -6.0 + 4.0*eta*eta ) );
}

double dem2dy2ms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double et2=eta*eta;
    double t1=2.0*std::exp(1.0)*emmax*(sigeta/x)*(sigeta/x)*std::exp(-et2);
    return t1*(et2*(2.0*et2-5.0)+1.0);
}

double dem4dy2ms(double& x, double& y)
{
    double sigeta=1.e1;
    double emmax=1.e-3;
    double eta=sigeta*y/x;
    double et2=eta*eta;
    double t1=0.5*std::exp(2.)*emmax*et2*(sigeta/x)*(sigeta/x)*std::exp(-et2);
    return t1*(et2*(2.0*et2-9.0)+6.0);
}

double ssams(double& x, double& y)
{
    double emu=1.e-6;
    double yloc=std::max(1.e-40,y);
    yloc=yloc*yloc;
    double u1  =ums(x,y);
    double u2  =vms(x,y);
    double emsa=em2ms(x,y);
    double desadx  =dem2dxms(x,y);
    double desady  =dem2dyms(x,y);
    double desadx2 =dem2dx2ms(x,y);
    double desady2 =dem2dy2ms(x,y);
    double cv1=u1*desadx+u2*desady;
    double sig=2./3.;
    double df1=(emu+emsa)*(desadx2+desady2);
    df1=df1+1.622*(desadx*desadx+desady*desady);
    df1=-df1/sig;
    double emt=emsa/emu;
    double fv1=emsa*emsa*emsa/(emsa*emsa*emsa+(7.1*emu)*(7.1*emu)*(7.1*emu));
    double fv2=1.0-emt/(emt*fv1+1.0);
    double dyp=emsa/0.41/0.41/yloc;
    double str=vortms(x,y)+dyp*fv2;
    double cp1=-0.1355*str*emsa;
    double cw1=0.1355/0.41/0.41+1.622/sig;
    double r=dyp/str;
    double g=r*(1.0+0.3*(std::pow(r,5.0)-1.0));
    double fw=65.0/(std::pow(g,6.0)+64.0);
    fw=g*(std::pow(fw,1.0/6.0));
    double cd1=cw1*fw*emsa*emsa/yloc;
    return cv1+df1+cp1+cd1;
}

double eddysams(double& x, double& y)
{
    double emsa=em2ms(x,y);
    return fv1sams(emsa)*emsa;
}

double desadxms(double& x, double& y)
{
    double emsa=em2ms(x,y);
    return (fv1sams(emsa)+dfv1sams(emsa)*emsa)*dem2dxms(x,y);
}

double desadyms(double& x, double& y)
{
    double emsa=em2ms(x,y);
    return (fv1sams(emsa)+dfv1sams(emsa)*emsa)*dem2dyms(x,y);
}

double fv1sams(double& em)
{
    double emu=1.e-6;
    double cv1=(7.1*emu)*(7.1*emu)*(7.1*emu);
    return em*em*em/(em*em*em+cv1);
}

double dfv1sams(double& em)
{
    double emu=1.e-6;
    double cv1=(7.1*emu)*(7.1*emu)*(7.1*emu);
    return 3.0*em*em*cv1/(em*em*em+cv1)/(em*em*em+cv1);
}



double smtms(double& x, double& y)
{
    double emu=1.e-6;
    double u1  =ums(x,y);
    double u2  =vms(x,y);
    double emmt=em2ms(x,y);
    double emtu=eddymtms(x,y);
    double demtdx  =dem2dxms(x,y);
    double demtdy  =dem2dyms(x,y);
    double demtdx2 =dem2dx2ms(x,y);
    double demtdy2 =dem2dy2ms(x,y);
    double cv1=u1*demtdx+u2*demtdy;
    double sig=1.0;
    double df1=-(emu+emmt/sig)*(demtdx2+demtdy2);
    df1=df1-(demtdx*demtdx+demtdy*demtdy)/sig;
    double fd1=(emtu+emu)/(emmt+emu);
    double str=strainms(x,y);
    double cp1=-0.144*fd1*str*emmt;
    double u1dx  =dudxms(x,y);
    double u1dy  =dudyms(x,y);
    double u1dx2 =dudx2ms(x,y);
    double u1dy2 =dudy2ms(x,y);
    double u1dxy2=dudxyms(x,y);
    double u2dx  =dvdxms(x,y);
    double u2dy  =dvdyms(x,y);
    double u2dx2 =dvdx2ms(x,y);
    double u2dy2 =dvdy2ms(x,y);
    double u2dxy2=dvdxyms(x,y);
    double dsdx=4.0*(u1dx*u1dx2+u2dy*u2dxy2)+2.0*(u1dy+u2dx)*(u1dxy2+u2dx2);
    double dsdy=4.0*(u1dx*u1dxy2+u2dy*u2dy2)+2.0*(u1dy+u2dx)*(u1dy2+u2dxy2);
    double dsdt=(dsdx*dsdx+dsdy*dsdy)/(4.0*str*str);
    double ekep=emmt*emmt*dsdt/str/str;
    double epbb=demtdx*demtdx+demtdy*demtdy;
    double ep1e=7.0*epbb*std::tanh(ekep/(7.0*epbb));
    double cd1 =1.862*ep1e;
    return cv1+df1+cp1+cd1;
}

double eddymtms(double& x, double& y)
{
    double emmt=em2ms(x,y);
    return d2mtms(emmt)*emmt;
}

double demtdxms(double& x, double& y)
{
    double emmt=em2ms(x,y);
    return (d2mtms(emmt)+dd2mtms(emmt)*emmt)*dem2dxms(x,y);
}

double demtdyms(double& x, double& y)
{
    double emmt=em2ms(x,y);
    return (d2mtms(emmt)+dd2mtms(emmt)*emmt)*dem2dyms(x,y);
}

double d2mtms(double& em)
{
    double emu=1.e-6;
    double cv1 =13.0*0.41*emu;
    double cdum=(em/cv1)*(em/cv1);
    return 1.0-std::exp(-cdum);
}

double dd2mtms(double& em)
{
    double emu=1.e-6;
    double cv1    =13.0*0.41*emu;
    double cdum   =em/cv1;
    return 2.0*cdum*std::exp(-cdum*cdum)/cv1;
}

double user_derfc(double& x)
{
    const double pv= 1.26974899965115684e+01;
    const double ph= 6.10399733098688199e+00;
    const double p0= 2.96316885199227378e-01;
    const double p1= 1.81581125134637070e-01;
    const double p2= 6.81866451424939493e-02;
    const double p3= 1.56907543161966709e-02;
    const double p4= 2.21290116681517573e-03;
    const double p5= 1.91395813098742864e-04;
    const double p6= 9.71013284010551623e-06;
    const double p7= 1.66642447174307753e-07;
    const double q0= 6.12158644495538758e-02;
    const double q1= 5.50942780056002085e-01;
    const double q2= 1.53039662058770397e+00;
    const double q3= 2.99957952311300634e+00;
    const double q4= 4.95867777128246701e+00;
    const double q5= 7.41471251099335407e+00;
    const double q6= 1.04765104356545238e+01;
    const double q7= 1.48455557345597957e+01;
    double y=x*x;
    y=std::exp(-y)*x*(p7/(y+q7)+p6/(y+q6)
      +p5/(y+q5)+p4/(y+q4)+p3/(y+q3)
      +p2/(y+q2)+p1/(y+q1)+p0/(y+q0));
    if(x < ph) y=y+2.0/(std::exp(pv*x)+1.0);
    return y;
}

double user_derf(double& x)
{
    const double p0= 1.12837916709551257e+00;
    const double p1=-3.76126389031833602e-01;
    const double p2= 1.12837916706621301e-01;
    const double p3=-2.68661698447642378e-02;
    const double p4= 5.22387877685618101e-03;
    const double p5=-8.49202435186918470e-04;
    double y=std::abs(x);
    if(y > 0.125) {
        return std::copysign(1.0-user_derfc(y),x);
    } else {
        y=x*x;
        return (((((p5*y+p4)*y+p3)*y+p2)*y+p1)*y+p0)*x;
    }
}

void exact_polyfg(double& scale, FortranArray1DRef<double> polyfc,
    double& x, double& y, double& z,
    double& q, double& qx, double& qy, double& qz,
    double& qxx, double& qyy, double& qzz,
    double& qxy, double& qxz, double& qyz)
{
    double c     = polyfc( 1);
    double cx    = polyfc( 2);
    double cy    = polyfc( 3);
    double cz    = polyfc( 4);
    double cxx   = polyfc( 5);
    double cyy   = polyfc( 6);
    double czz   = polyfc( 7);
    double cxxx  = polyfc( 8);
    double cyyy  = polyfc( 9);
    double czzz  = polyfc(10);
    double cxxxx = polyfc(11);
    double cyyyy = polyfc(12);
    double czzzz = polyfc(13);

    q = q + scale*( c
         + cx*x + 0.5*cxx*(x*x) + 1./6.*cxxx*(x*x*x)
         + cy*y + 0.5*cyy*(y*y) + 1./6.*cyyy*(y*y*y)
         + cz*z + 0.5*czz*(z*z) + 1./6.*czzz*(z*z*z)
         + 1./12.*cxxxx*(x*x*x*x)
         + 1./12.*cyyyy*(y*y*y*y)
         + 1./12.*czzzz*(z*z*z*z) );

    qx =  qx + scale*( cx + cxx*x + 0.5*cxxx*(x*x) +
       1./3.*cxxxx*(x*x*x) );
    qy =  qy + scale*( cy + cyy*y + 0.5*cyyy*(y*y) +
       1./3.*cyyyy*(y*y*y) );
    qz =  qz + scale*( cz + czz*z + 0.5*czzz*(z*z) +
       1./3.*czzzz*(z*z*z) );

    qxx = qxx + scale*( cxx + cxxx*x + cxxxx*(x*x) );
    qyy = qyy + scale*( cyy + cyyy*y + cyyyy*(y*y) );
    qzz = qzz + scale*( czz + czzz*z + czzzz*(z*z) );

    qxy = qxy + scale*( 0.0 );
    qxz = qxz + scale*( 0.0 );
    qyz = qyz + scale*( 0.0 );
}


} // namespace mms_ns
