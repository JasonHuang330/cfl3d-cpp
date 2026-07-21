// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// collapse.cpp — faithful C++ translation of collapse.F (CFL3D ronnie).
// Check for collapsed points in grid, and expand any collapsed lines detected.
//
// /tol/ is declared here as (epsc,epsc0,epsreen,eps); by POSITION `eps` = the
// 4th slot = tol_.eps (CONTRACT rule 3). /sklt1/ -> sklt1_.isklt1.
//   xo/yo/zo(jdim,kdim) at (j,k) -> (j-1)+(k-1)*jdim
#include "ron_common.h"

void collapse(int jdim,int kdim,int jmaxo,int kmaxo,
              double* xo,double* yo,double* zo,
              int* nou,FStr* bou,int nbuf,int ibufdim)
{
    (void)kdim; (void)nbuf;
    const double eps    = tol_.eps;      // 4th /tol/ slot by position
    const int    isklt1 = sklt1_.isklt1;
    auto IX = [&](int j,int k){ return (long)(j-1) + (long)(k-1)*jdim; };

    const double fact  = 1.0e5;
    const double epss1 = fact*eps;
    const double epss2 = 1.0e-01;

    char line[160];

    //******************************
    //     check along lines k=const
    //******************************
    for (int k = 1; k <= kmaxo; ++k) {
        int kc     = k;
        int j1     = 1;
        int ncount = 0;
        bool polar = false;
    L110:
        {
            int j2 = j1 + 1 + ncount;
            if (f_abs(xo[IX(j2,k)]-xo[IX(j1,k)]) < eps &&
                f_abs(yo[IX(j2,k)]-yo[IX(j1,k)]) < eps &&
                f_abs(zo[IX(j2,k)]-zo[IX(j1,k)]) < eps) {
                ncount = ncount+1;
                if (j1 == 1 && j2 == jmaxo) { polar = true; goto L109; }  // go to 109
                goto L110;
            }
            if (ncount > 0) {
                j2 = j2-1;
                if (isklt1 > 0) {
                    std::snprintf(line,sizeof(line),
                        " %s%s%3d%s%3d%s%3d",
                        "     ","collapsed boundary on k =",kc,
                        "  between j =",j1,"  and j =",j2);
                    bou_put(nou,bou,ibufdim,4,line);
                }
                // 1) collapsed segments only over part of the k=const line
                double sx0,sy0,sz0;
                if (j1 > 1) {
                    sx0 = xo[IX(j1,k)]-xo[IX(j1-1,k)];
                    sy0 = yo[IX(j1,k)]-yo[IX(j1-1,k)];
                    sz0 = zo[IX(j1,k)]-zo[IX(j1-1,k)];
                } else {
                    sx0 = xo[IX(j2+1,k)]-xo[IX(j2,k)];
                    sy0 = yo[IX(j2+1,k)]-yo[IX(j2,k)];
                    sz0 = zo[IX(j2+1,k)]-zo[IX(j2,k)];
                }
                double xoo = xo[IX(j1,k)];
                double yoo = yo[IX(j1,k)];
                double zoo = zo[IX(j1,k)];
                int jc;
                if (j1 > 1 && j2 < jmaxo) {
                    jc = j1 + ncount/2;
                } else if (j1 == 1) {
                    jc = j1;
                } else {
                    jc = j2;
                }
                for (int jj = j1; jj <= j2; ++jj) {
                    double fact1 = f_float(jj-jc)*epss1;
                    xo[IX(jj,k)] = xoo+fact1*sx0;
                    yo[IX(jj,k)] = yoo+fact1*sy0;
                    zo[IX(jj,k)] = zoo+fact1*sz0;
                }
            }
            j1     = j1 + ncount + 1;
            ncount = 0;
            if (j1 < jmaxo-1) goto L110;
            continue;   // go to 111 (next k)
        }
    L109:
        if (polar) {
            // 2) collapsed over entire k=const line - polar singularity
            int j1 = 1, j2 = jmaxo;   // preserve values at branch (j1==1,j2==jmaxo)
            if (isklt1 > 0) {
                std::snprintf(line,sizeof(line),
                    " %s%s%3d%s%3d%s%3d",
                    "     ","collapsed boundary on k =",kc,
                    "  between j =",j1,"  and j =",j2);
                bou_put(nou,bou,ibufdim,4,line);
            }
            int m = 1;
            if (k == kmaxo) m = -1;
            for (int j = 1; j <= jmaxo; ++j) {
                double sx = xo[IX(j,k+m)]-xo[IX(j,k)];
                double sy = yo[IX(j,k+m)]-yo[IX(j,k)];
                double sz = zo[IX(j,k+m)]-zo[IX(j,k)];
                xo[IX(j,k)] = xo[IX(j,k)]+sx*epss2;
                yo[IX(j,k)] = yo[IX(j,k)]+sy*epss2;
                zo[IX(j,k)] = zo[IX(j,k)]+sz*epss2;
            }
        }
        // falls through to 111 continue (next k)
    }

    //******************************
    //     check along lines j=const
    //******************************
    for (int j = 1; j <= jmaxo; ++j) {
        int jc     = j;
        int k1     = 1;
        int ncount = 0;
        bool polar = false;
    L1101:
        {
            int k2 = k1 + 1 + ncount;
            if (f_abs(xo[IX(j,k2)]-xo[IX(j,k1)]) < eps &&
                f_abs(yo[IX(j,k2)]-yo[IX(j,k1)]) < eps &&
                f_abs(zo[IX(j,k2)]-zo[IX(j,k1)]) < eps) {
                ncount = ncount+1;
                if (k1 == 1 && k2 == kmaxo) { polar = true; goto L1009; }
                goto L1101;
            }
            if (ncount > 0) {
                k2 = k2-1;
                if (isklt1 > 0) {
                    std::snprintf(line,sizeof(line),
                        " %s%s%3d%s%3d%s%3d",
                        "     ","collapsed boundary on j =",jc,
                        "  between k =",k1,"  and k =",k2);
                    bou_put(nou,bou,ibufdim,4,line);
                }
                double sx0,sy0,sz0;
                if (k1 > 1) {
                    sx0 = xo[IX(j,k1)]-xo[IX(j,k1-1)];
                    sy0 = yo[IX(j,k1)]-yo[IX(j,k1-1)];
                    sz0 = zo[IX(j,k1)]-zo[IX(j,k1-1)];
                } else {
                    sx0 = xo[IX(j,k2+1)]-xo[IX(j,k2)];
                    sy0 = yo[IX(j,k2+1)]-yo[IX(j,k2)];
                    sz0 = zo[IX(j,k2+1)]-zo[IX(j,k2)];
                }
                double xoo = xo[IX(j,k1)];
                double yoo = yo[IX(j,k1)];
                double zoo = zo[IX(j,k1)];
                int kc;
                if (k1 > 1 && k2 < kmaxo) {
                    kc = k1 + ncount/2;
                } else if (k1 == 1) {
                    kc = k1;
                } else {
                    kc = k2;
                }
                for (int kk = k1; kk <= k2; ++kk) {
                    double fact1 = f_float(kk-kc)*epss1;
                    xo[IX(j,kk)] = xoo+fact1*sx0;
                    yo[IX(j,kk)] = yoo+fact1*sy0;
                    zo[IX(j,kk)] = zoo+fact1*sz0;
                }
            }
            k1     = k1 + ncount + 1;
            ncount = 0;
            if (k1 < kmaxo-1) goto L1101;
            continue;   // go to 1111 (next j)
        }
    L1009:
        if (polar) {
            int k1 = 1, k2 = kmaxo;
            if (isklt1 > 0) {
                std::snprintf(line,sizeof(line),
                    " %s%s%3d%s%3d%s%3d",
                    "     ","collapsed boundary on j =",jc,
                    "  between k =",k1,"  and k =",k2);
                bou_put(nou,bou,ibufdim,4,line);
            }
            int m = 1;
            if (j == jmaxo) m = -1;
            for (int k = 1; k <= kmaxo; ++k) {
                double sx = xo[IX(j+m,k)]-xo[IX(j,k)];
                double sy = yo[IX(j+m,k)]-yo[IX(j,k)];
                double sz = zo[IX(j+m,k)]-zo[IX(j,k)];
                xo[IX(j,k)] = xo[IX(j,k)]+sx*epss2;
                yo[IX(j,k)] = yo[IX(j,k)]+sy*epss2;
                zo[IX(j,k)] = zo[IX(j,k)]+sz*epss2;
            }
        }
    }
}
