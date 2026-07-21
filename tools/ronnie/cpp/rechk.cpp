// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// rechk.cpp — faithful C++ translation of rechk.F (CFL3D ronnie).
// Check for branch cuts; build jimage/kimage image-index arrays.
//
// /tol/ declared in rechk.F as (epsc,epsc0,eps,epscoll): its `eps` is the 3rd
// slot = tol_.epsreen (map BY POSITION, CONTRACT rule 3). /sklt1/ -> isklt1.
//
//   x1/y1/z1(mdim1,ndim1)          at (j,k)   -> (j-1)+(k-1)*mdim1
//   jimage/kimage(msub1,mdim1,ndim1) at (l,j,k)->
//        (l-1)+(j-1)*msub1+(k-1)*msub1*mdim1
#include "ron_common.h"

void rechk(int mdim1,int ndim1,int* jimage,int* kimage,int msub1,
           int jmax1,int kmax1,int l,double* x1,double* y1,double* z1,
           int xie1,int xie2,int eta1,int eta2,
           int* nou,FStr* bou,int nbuf,int ibufdim,int myid,int* mblk2nd,int maxbl)
{
    (void)ndim1; (void)nbuf; (void)myid; (void)mblk2nd; (void)maxbl;
    const double eps    = tol_.epsreen;   // 3rd /tol/ slot by position
    const int    isklt1 = sklt1_.isklt1;

    auto IX = [&](int j,int k){ return (long)(j-1) + (long)(k-1)*mdim1; };
    auto JI = [&](int ll,int j,int k)->int& {
        return jimage[(long)(ll-1) + (long)(j-1)*msub1 + (long)(k-1)*msub1*mdim1]; };
    auto KI = [&](int ll,int j,int k)->int& {
        return kimage[(long)(ll-1) + (long)(j-1)*msub1 + (long)(k-1)*msub1*mdim1]; };

    char line[200];

    // default to no branch cuts
    for (int j = 1; j <= jmax1+1; ++j) {
        for (int k = 1; k <= kmax1+1; ++k) {
            JI(l,j,k) = j;
            KI(l,j,k) = k;
        }
    }

    // check for C-grid type branch cut on k=1 and k=kmax1 boundaries
    if (!(eta1 > 1 && eta2 < kmax1)) {   // if not skipping (go to 11)
        for (int k = 1; k <= kmax1; k += (kmax1-1)) {
            int kcheck = 0;
            for (int js = xie1; js <= xie2-1; ++js) {
                for (int jq = js+1; jq <= xie2; ++jq) {
                    double delxyz = ccabs(x1[IX(js,k)]-x1[IX(jq,k)])
                                  + ccabs(y1[IX(js,k)]-y1[IX(jq,k)])
                                  + ccabs(z1[IX(js,k)]-z1[IX(jq,k)]);
                    if (delxyz < eps) {
                        kcheck = kcheck+1;
                        if (k == 1) {
                            JI(l,js+1,1) = jq;   KI(l,js+1,1) = 2;
                            JI(l,jq,1)   = js+1; KI(l,jq,1)   = 2;
                        }
                        if (k == kmax1) {
                            JI(l,js+1,kmax1+1) = jq;   KI(l,js+1,kmax1+1) = kmax1;
                            JI(l,jq,kmax1+1)   = js+1; KI(l,jq,kmax1+1)   = kmax1;
                        }
                    }
                }
            }
            if (kcheck == 0) continue;   // go to 10
            if (k == 1) {
                // correction for border points between start/end of cuts
                for (int j = xie1+1; j <= xie2-1; ++j) {
                    if (JI(l,j+1,1) == j+1) {
                        JI(l, JI(l,j,1), 1) = JI(l,j,1);
                        KI(l, JI(l,j,1), 1) = 1;
                        JI(l,j,1) = j;
                        KI(l,j,1) = 1;
                        break;   // go to 1809
                    }
                }
                // image point at ends of boundary of extended grid
                if (JI(l,xie1+1,1) == xie2) {
                    JI(l,xie1,1)   = xie2+1; KI(l,xie1,1)   = 2;
                    JI(l,xie2+1,1) = xie1;   KI(l,xie2+1,1) = 2;
                }
            }
            if (k == kmax1) {
                for (int j = xie1+1; j <= xie2-1; ++j) {
                    if (JI(l,j+1,kmax1+1) == j+1) {
                        JI(l, JI(l,j,kmax1+1), kmax1+1) = JI(l,j,kmax1+1);
                        KI(l, JI(l,j,kmax1+1), kmax1+1) = kmax1+1;
                        JI(l,j,kmax1+1) = j;
                        KI(l,j,kmax1+1) = kmax1+1;
                        break;   // go to 1808
                    }
                }
                // image point at left boundary of extended grid
                if (JI(l,xie1+1,kmax1+1) == xie2) {
                    JI(l,xie1,kmax1+1)   = xie2+1; KI(l,xie1,kmax1+1)   = kmax1;
                    JI(l,xie2+1,kmax1+1) = xie1;   KI(l,xie2+1,kmax1+1) = kmax1;
                }
            }
        }
    }
    // 11 continue

    // check for C-grid type branch cut on j=1 and j=jmax1 boundaries
    if (!(xie1 > 1 && xie2 < jmax1)) {   // if not skipping (go to 21)
        for (int j = 1; j <= jmax1; j += (jmax1-1)) {
            int jcheck = 0;
            for (int ks = eta1; ks <= eta2-1; ++ks) {
                for (int kq = ks+1; kq <= eta2; ++kq) {
                    double delxyz = ccabs(x1[IX(j,ks)]-x1[IX(j,kq)])
                                  + ccabs(y1[IX(j,ks)]-y1[IX(j,kq)])
                                  + ccabs(z1[IX(j,ks)]-z1[IX(j,kq)]);
                    if (delxyz < eps) {
                        jcheck = jcheck+1;
                        if (j == 1) {
                            JI(l,1,ks+1) = 2;    KI(l,1,ks+1) = kq;
                            JI(l,1,kq)   = 2;    KI(l,1,kq)   = ks+1;
                        }
                        if (j == jmax1) {
                            JI(l,jmax1+1,ks+1) = jmax1; KI(l,jmax1+1,ks+1) = kq;
                            JI(l,jmax1+1,kq)   = jmax1; KI(l,jmax1+1,kq)   = ks+1;
                        }
                        break;   // go to 2712 (exit inner kq loop)
                    }
                }
            }
            if (jcheck == 0) continue;   // go to 20
            if (j == 1) {
                for (int k = eta1+1; k <= eta2-1; ++k) {
                    if (KI(l,1,k+1) == k+1) {
                        JI(l,1, KI(l,1,k)) = 1;
                        KI(l,1, KI(l,1,k)) = KI(l,1,k);
                        JI(l,1,k) = 1;
                        KI(l,1,k) = k;
                        break;   // go to 2809
                    }
                }
                if (KI(l,1,eta1+1) == eta2) {
                    JI(l,1,eta1)   = 2; KI(l,1,eta1)   = eta2+1;
                    JI(l,1,eta2+1) = 2; KI(l,1,eta2+1) = eta1;
                }
            }
            if (j == jmax1) {
                for (int k = eta1+1; k <= eta2-1; ++k) {
                    if (KI(l,jmax1+1,k+1) == k+1) {
                        JI(l,jmax1+1, KI(l,jmax1+1,k)) = jmax1+1;
                        KI(l,jmax1+1, KI(l,jmax1+1,k)) = KI(l,jmax1+1,k);
                        JI(l,jmax1+1,k) = jmax1+1;
                        KI(l,jmax1+1,k) = k;
                        break;   // go to 2808
                    }
                }
                if (KI(l,jmax1+1,eta1+1) == eta2) {
                    JI(l,jmax1+1,eta1)   = jmax1; KI(l,jmax1+1,eta1)   = eta2+1;
                    JI(l,jmax1+1,eta2+1) = jmax1; KI(l,jmax1+1,eta2+1) = eta1;
                }
            }
        }
    }
    // 21 continue

    // check for O-grid type branch cut on k=1 and k=kmax1 boundaries
    {
        int jcount = 0;
        if (!(eta1 > 1 && eta2 < kmax1)) {   // if not skipping (go to 31)
            for (int j = xie1; j <= xie2; ++j) {
                double delxyz = ccabs(x1[IX(j,1)]-x1[IX(j,kmax1)])
                              + ccabs(y1[IX(j,1)]-y1[IX(j,kmax1)])
                              + ccabs(z1[IX(j,1)]-z1[IX(j,kmax1)]);
                if (delxyz < eps) jcount = jcount+1;
            }
            if (jcount > 1 && jcount != xie2-xie1+1) {
                if (isklt1 > 1) {
                    std::snprintf(line,sizeof(line),
                        " WARNING...O-type branch cut does not extend over"
                        " entire k=constant  boundary, as assumed");
                    bou_put(nou,bou,ibufdim,4,line);
                }
            }
            if (jcount > 1) {
                for (int j = xie1; j <= xie2+1; ++j) {
                    JI(l,j,1)       = j; KI(l,j,1)       = kmax1;
                    JI(l,j,kmax1+1) = j; KI(l,j,kmax1+1) = 2;
                }
            }
        }
    }
    // 31 continue

    // check for O-grid type branch cut on j=1 and j=jmax1 boundaries
    {
        int kcount = 0;
        if (!(xie1 > 1 && xie2 < jmax1)) {   // if not skipping (go to 41)
            for (int k = eta1; k <= eta2; ++k) {
                double delxyz = ccabs(x1[IX(1,k)]-x1[IX(jmax1,k)])
                              + ccabs(y1[IX(1,k)]-y1[IX(jmax1,k)])
                              + ccabs(z1[IX(1,k)]-z1[IX(jmax1,k)]);
                if (delxyz < eps) kcount = kcount+1;
            }
            if (kcount > 1 && kcount != eta2-eta1+1) {
                if (isklt1 > 0) {
                    std::snprintf(line,sizeof(line),
                        " WARNING...O-type branch cut does not extend over"
                        " entire k=constant  boundary, as assumed");
                    bou_put(nou,bou,ibufdim,4,line);
                }
            }
            if (kcount > 1) {
                for (int k = eta1; k <= eta2+1; ++k) {
                    JI(l,1,k)       = jmax1; KI(l,1,k)       = k;
                    JI(l,jmax1+1,k) = 2;     KI(l,jmax1+1,k) = k;
                }
            }
        }
    }
    // 41 continue

    // print out branch cut topology (j direction)
    for (int j = xie1; j <= xie2+1; j += (xie2+1-xie1)) {
        int iflag = 0;
        int jpr   = j;
        if (j == jmax1+1) jpr = jmax1;
        for (int k = eta1+1; k <= eta2; ++k) {
            if (JI(l,j,k) != j || KI(l,j,k) != k) iflag = 1;
        }
        int iflago = 0;
        for (int k = eta1+1; k <= eta2; ++k) {
            if (JI(l,j,k) == j || KI(l,j,k) != k) iflago = 1;
        }
        if (iflag != 0 && iflago == 0 && isklt1 > 0) {
            std::snprintf(line,sizeof(line),
                " %s%3d%s","     j=",jpr,
                " is an O-type (periodic) branch cut boundary");
            bou_put(nou,bou,ibufdim,4,line);
        }
        if (iflag != 0 && iflago != 0 && isklt1 > 0) {
            std::snprintf(line,sizeof(line),
                " %s%3d%s","     j=",jpr,
                " is a C-type branch cut boundary");
            bou_put(nou,bou,ibufdim,4,line);
        }
    }
    // print out branch cut topology (k direction)
    for (int k = eta1; k <= eta2+1; k += (eta2+1-eta1)) {
        int kpr = k;
        if (k == kmax1+1) kpr = kmax1;
        int iflag = 0;
        for (int j = xie1+1; j <= xie2; ++j) {
            if (JI(l,j,k) != j || KI(l,j,k) != k) iflag = 1;
        }
        int iflago = 0;
        for (int j = xie1+1; j <= xie2; ++j) {
            if (JI(l,j,k) != j || KI(l,j,k) == k) iflago = 1;
        }
        if (iflag != 0 && iflago == 0 && isklt1 > 0) {
            std::snprintf(line,sizeof(line),
                " %s%3d%s","     k=",kpr,
                " is an O-type (periodic) branch cut boundary");
            bou_put(nou,bou,ibufdim,4,line);
        }
        if (iflag != 0 && iflago != 0 && isklt1 > 0) {
            std::snprintf(line,sizeof(line),
                " %s%3d%s","     k=",kpr,
                " is a C-type branch cut boundary");
            bou_put(nou,bou,ibufdim,4,line);
        }
    }
}
