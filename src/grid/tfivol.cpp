// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "tfivol.h"
#include <cmath>
#include <algorithm>

namespace tfivol_ns {

void tfivol(int& idim, int& jdim, int& kdim, int& iskp, int& jskp, int& kskp,
            FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax,
            FortranArray1DRef<int> kskmax, FortranArray2DRef<int> iskip,
            FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip,
            int& isktyp, FortranArray3DRef<double> x, FortranArray3DRef<double> y,
            FortranArray3DRef<double> z, FortranArray3DRef<double> x1,
            FortranArray3DRef<double> y1, FortranArray3DRef<double> z1,
            FortranArray3DRef<double> x2, FortranArray3DRef<double> y2,
            FortranArray3DRef<double> z2, FortranArray3DRef<double> x3,
            FortranArray3DRef<double> y3, FortranArray3DRef<double> z3,
            FortranArray3DRef<double> x4, FortranArray3DRef<double> y4,
            FortranArray3DRef<double> z4, FortranArray3DRef<double> x5,
            FortranArray3DRef<double> y5, FortranArray3DRef<double> z5,
            FortranArray3DRef<double> x6, FortranArray3DRef<double> y6,
            FortranArray3DRef<double> z6, FortranArray3DRef<double> x7,
            FortranArray3DRef<double> y7, FortranArray3DRef<double> z7,
            FortranArray3DRef<double> arci, FortranArray3DRef<double> arcj,
            FortranArray3DRef<double> arck, FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
            int& myid, int& maxbl, int& nbl)
{
    // local arrays (1-indexed, size 2)
    double phi[3], psi[3], omg[3];

    // tolerance for switch to linear blending function
    // (10.**(-iexp) is machine zero)
    double tol = std::max(1.e-07, std::pow(10.0, (double)(-cmn_zero.iexp + 1)));

    int i1, j1, k1, is, ie, js, je, ks, ke, i, j, k;
    double eta, zeta, xi;
    double denomj, denomk, denomi;

    if (std::abs(isktyp) == 1) {
        for (i1 = 1; i1 <= idim - iskp; i1 += iskp) {
            is = i1;
            ie = i1 + iskp;
            for (j1 = 1; j1 <= jdim - jskp; j1 += jskp) {
                js = j1;
                je = j1 + jskp;
                for (k1 = 1; k1 <= kdim - kskp; k1 += kskp) {
                    ks = k1;
                    ke = k1 + kskp;

                    for (j = js; j <= je; j++) {
                        for (k = ks; k <= ke; k++) {
                            for (i = is; i <= ie; i++) {
                                denomj = (arcj(je,k,i) - arcj(js,k,i));
                                if ((float)denomj < (float)tol) {
                                    eta = 0.;
                                } else {
                                    eta = (arcj(j,k,i) - arcj(js,k,i)) / denomj;
                                }
                                psi[1] = eta;
                                psi[2] = 1. - eta;

                                denomk = (arck(j,ke,i) - arck(j,ks,i));
                                if ((float)denomk < (float)tol) {
                                    zeta = 0.;
                                } else {
                                    zeta = (arck(j,k,i) - arck(j,ks,i)) / denomk;
                                }
                                omg[1] = zeta;
                                omg[2] = 1. - zeta;

                                denomi = (arci(j,k,ie) - arci(j,k,is));
                                if ((float)denomi < (float)tol) {
                                    xi = 0.;
                                } else {
                                    xi = (arci(j,k,i) - arci(j,k,is)) / denomi;
                                }
                                phi[1] = xi;
                                phi[2] = 1. - xi;

                                x1(j,k,i) = psi[2]*x(js,k,i) + psi[1]*x(je,k,i);
                                y1(j,k,i) = psi[2]*y(js,k,i) + psi[1]*y(je,k,i);
                                z1(j,k,i) = psi[2]*z(js,k,i) + psi[1]*z(je,k,i);

                                x2(j,k,i) = omg[2]*x(j,ks,i) + omg[1]*x(j,ke,i);
                                y2(j,k,i) = omg[2]*y(j,ks,i) + omg[1]*y(j,ke,i);
                                z2(j,k,i) = omg[2]*z(j,ks,i) + omg[1]*z(j,ke,i);

                                x3(j,k,i) = phi[2]*x(j,k,is) + phi[1]*x(j,k,ie);
                                y3(j,k,i) = phi[2]*y(j,k,is) + phi[1]*y(j,k,ie);
                                z3(j,k,i) = phi[2]*z(j,k,is) + phi[1]*z(j,k,ie);

                                x4(j,k,i) = psi[2]*omg[2]*x(js,ks,i)
                                           +psi[2]*omg[1]*x(js,ke,i)
                                           +psi[1]*omg[2]*x(je,ks,i)
                                           +psi[1]*omg[1]*x(je,ke,i);
                                y4(j,k,i) = psi[2]*omg[2]*y(js,ks,i)
                                           +psi[2]*omg[1]*y(js,ke,i)
                                           +psi[1]*omg[2]*y(je,ks,i)
                                           +psi[1]*omg[1]*y(je,ke,i);
                                z4(j,k,i) = psi[2]*omg[2]*z(js,ks,i)
                                           +psi[2]*omg[1]*z(js,ke,i)
                                           +psi[1]*omg[2]*z(je,ks,i)
                                           +psi[1]*omg[1]*z(je,ke,i);

                                x5(j,k,i) = psi[2]*phi[2]*x(js,k,is)
                                           +psi[2]*phi[1]*x(js,k,ie)
                                           +psi[1]*phi[2]*x(je,k,is)
                                           +psi[1]*phi[1]*x(je,k,ie);
                                y5(j,k,i) = psi[2]*phi[2]*y(js,k,is)
                                           +psi[2]*phi[1]*y(js,k,ie)
                                           +psi[1]*phi[2]*y(je,k,is)
                                           +psi[1]*phi[1]*y(je,k,ie);
                                z5(j,k,i) = psi[2]*phi[2]*z(js,k,is)
                                           +psi[2]*phi[1]*z(js,k,ie)
                                           +psi[1]*phi[2]*z(je,k,is)
                                           +psi[1]*phi[1]*z(je,k,ie);

                                x6(j,k,i) = omg[2]*phi[2]*x(j,ks,is)
                                           +omg[2]*phi[1]*x(j,ks,ie)
                                           +omg[1]*phi[2]*x(j,ke,is)
                                           +omg[1]*phi[1]*x(j,ke,ie);
                                y6(j,k,i) = omg[2]*phi[2]*y(j,ks,is)
                                           +omg[2]*phi[1]*y(j,ks,ie)
                                           +omg[1]*phi[2]*y(j,ke,is)
                                           +omg[1]*phi[1]*y(j,ke,ie);
                                z6(j,k,i) = omg[2]*phi[2]*z(j,ks,is)
                                           +omg[2]*phi[1]*z(j,ks,ie)
                                           +omg[1]*phi[2]*z(j,ke,is)
                                           +omg[1]*phi[1]*z(j,ke,ie);

                                x7(j,k,i) = psi[2]*omg[2]*phi[2]*x(js,ks,is)
                                           +psi[2]*omg[2]*phi[1]*x(js,ks,ie)
                                           +psi[2]*omg[1]*phi[2]*x(js,ke,is)
                                           +psi[2]*omg[1]*phi[1]*x(js,ke,ie)
                                           +psi[1]*omg[2]*phi[2]*x(je,ks,is)
                                           +psi[1]*omg[2]*phi[1]*x(je,ks,ie)
                                           +psi[1]*omg[1]*phi[2]*x(je,ke,is)
                                           +psi[1]*omg[1]*phi[1]*x(je,ke,ie);
                                y7(j,k,i) = psi[2]*omg[2]*phi[2]*y(js,ks,is)
                                           +psi[2]*omg[2]*phi[1]*y(js,ks,ie)
                                           +psi[2]*omg[1]*phi[2]*y(js,ke,is)
                                           +psi[2]*omg[1]*phi[1]*y(js,ke,ie)
                                           +psi[1]*omg[2]*phi[2]*y(je,ks,is)
                                           +psi[1]*omg[2]*phi[1]*y(je,ks,ie)
                                           +psi[1]*omg[1]*phi[2]*y(je,ke,is)
                                           +psi[1]*omg[1]*phi[1]*y(je,ke,ie);
                                z7(j,k,i) = psi[2]*omg[2]*phi[2]*z(js,ks,is)
                                           +psi[2]*omg[2]*phi[1]*z(js,ks,ie)
                                           +psi[2]*omg[1]*phi[2]*z(js,ke,is)
                                           +psi[2]*omg[1]*phi[1]*z(js,ke,ie)
                                           +psi[1]*omg[2]*phi[2]*z(je,ks,is)
                                           +psi[1]*omg[2]*phi[1]*z(je,ks,ie)
                                           +psi[1]*omg[1]*phi[2]*z(je,ke,is)
                                           +psi[1]*omg[1]*phi[1]*z(je,ke,ie);
                            }
                        }
                    }

                    for (i = is+1; i <= ie-1; i++) {
                        for (j = js+1; j <= je-1; j++) {
                            for (k = ks+1; k <= ke-1; k++) {
                                x(j,k,i) = x1(j,k,i) + x2(j,k,i) + x3(j,k,i)
                                          - x4(j,k,i) - x5(j,k,i) - x6(j,k,i)
                                          + x7(j,k,i);
                                y(j,k,i) = y1(j,k,i) + y2(j,k,i) + y3(j,k,i)
                                          - y4(j,k,i) - y5(j,k,i) - y6(j,k,i)
                                          + y7(j,k,i);
                                z(j,k,i) = z1(j,k,i) + z2(j,k,i) + z3(j,k,i)
                                          - z4(j,k,i) - z5(j,k,i) - z6(j,k,i)
                                          + z7(j,k,i);
                            }
                        }
                    }

                } // k1
            } // j1
        } // i1
    } else {
        for (i1 = 1; i1 <= iskmax(nbl) - 1; i1++) {
            is = iskip(nbl, i1);
            ie = iskip(nbl, i1+1);
            for (j1 = 1; j1 <= jskmax(nbl) - 1; j1++) {
                js = jskip(nbl, j1);
                je = jskip(nbl, j1+1);
                for (k1 = 1; k1 <= kskmax(nbl) - 1; k1++) {
                    ks = kskip(nbl, k1);
                    ke = kskip(nbl, k1+1);

                    for (j = js; j <= je; j++) {
                        for (k = ks; k <= ke; k++) {
                            for (i = is; i <= ie; i++) {
                                denomj = (arcj(je,k,i) - arcj(js,k,i));
                                if ((float)denomj < (float)tol) {
                                    eta = 0.;
                                } else {
                                    eta = (arcj(j,k,i) - arcj(js,k,i)) / denomj;
                                }
                                psi[1] = eta;
                                psi[2] = 1. - eta;

                                denomk = (arck(j,ke,i) - arck(j,ks,i));
                                if ((float)denomk < (float)tol) {
                                    zeta = 0.;
                                } else {
                                    zeta = (arck(j,k,i) - arck(j,ks,i)) / denomk;
                                }
                                omg[1] = zeta;
                                omg[2] = 1. - zeta;

                                denomi = (arci(j,k,ie) - arci(j,k,is));
                                if ((float)denomi < (float)tol) {
                                    xi = 0.;
                                } else {
                                    xi = (arci(j,k,i) - arci(j,k,is)) / denomi;
                                }
                                phi[1] = xi;
                                phi[2] = 1. - xi;

                                x1(j,k,i) = psi[2]*x(js,k,i) + psi[1]*x(je,k,i);
                                y1(j,k,i) = psi[2]*y(js,k,i) + psi[1]*y(je,k,i);
                                z1(j,k,i) = psi[2]*z(js,k,i) + psi[1]*z(je,k,i);

                                x2(j,k,i) = omg[2]*x(j,ks,i) + omg[1]*x(j,ke,i);
                                y2(j,k,i) = omg[2]*y(j,ks,i) + omg[1]*y(j,ke,i);
                                z2(j,k,i) = omg[2]*z(j,ks,i) + omg[1]*z(j,ke,i);

                                x3(j,k,i) = phi[2]*x(j,k,is) + phi[1]*x(j,k,ie);
                                y3(j,k,i) = phi[2]*y(j,k,is) + phi[1]*y(j,k,ie);
                                z3(j,k,i) = phi[2]*z(j,k,is) + phi[1]*z(j,k,ie);

                                x4(j,k,i) = psi[2]*omg[2]*x(js,ks,i)
                                           +psi[2]*omg[1]*x(js,ke,i)
                                           +psi[1]*omg[2]*x(je,ks,i)
                                           +psi[1]*omg[1]*x(je,ke,i);
                                y4(j,k,i) = psi[2]*omg[2]*y(js,ks,i)
                                           +psi[2]*omg[1]*y(js,ke,i)
                                           +psi[1]*omg[2]*y(je,ks,i)
                                           +psi[1]*omg[1]*y(je,ke,i);
                                z4(j,k,i) = psi[2]*omg[2]*z(js,ks,i)
                                           +psi[2]*omg[1]*z(js,ke,i)
                                           +psi[1]*omg[2]*z(je,ks,i)
                                           +psi[1]*omg[1]*z(je,ke,i);

                                x5(j,k,i) = psi[2]*phi[2]*x(js,k,is)
                                           +psi[2]*phi[1]*x(js,k,ie)
                                           +psi[1]*phi[2]*x(je,k,is)
                                           +psi[1]*phi[1]*x(je,k,ie);
                                y5(j,k,i) = psi[2]*phi[2]*y(js,k,is)
                                           +psi[2]*phi[1]*y(js,k,ie)
                                           +psi[1]*phi[2]*y(je,k,is)
                                           +psi[1]*phi[1]*y(je,k,ie);
                                z5(j,k,i) = psi[2]*phi[2]*z(js,k,is)
                                           +psi[2]*phi[1]*z(js,k,ie)
                                           +psi[1]*phi[2]*z(je,k,is)
                                           +psi[1]*phi[1]*z(je,k,ie);

                                x6(j,k,i) = omg[2]*phi[2]*x(j,ks,is)
                                           +omg[2]*phi[1]*x(j,ks,ie)
                                           +omg[1]*phi[2]*x(j,ke,is)
                                           +omg[1]*phi[1]*x(j,ke,ie);
                                y6(j,k,i) = omg[2]*phi[2]*y(j,ks,is)
                                           +omg[2]*phi[1]*y(j,ks,ie)
                                           +omg[1]*phi[2]*y(j,ke,is)
                                           +omg[1]*phi[1]*y(j,ke,ie);
                                z6(j,k,i) = omg[2]*phi[2]*z(j,ks,is)
                                           +omg[2]*phi[1]*z(j,ks,ie)
                                           +omg[1]*phi[2]*z(j,ke,is)
                                           +omg[1]*phi[1]*z(j,ke,ie);

                                x7(j,k,i) = psi[2]*omg[2]*phi[2]*x(js,ks,is)
                                           +psi[2]*omg[2]*phi[1]*x(js,ks,ie)
                                           +psi[2]*omg[1]*phi[2]*x(js,ke,is)
                                           +psi[2]*omg[1]*phi[1]*x(js,ke,ie)
                                           +psi[1]*omg[2]*phi[2]*x(je,ks,is)
                                           +psi[1]*omg[2]*phi[1]*x(je,ks,ie)
                                           +psi[1]*omg[1]*phi[2]*x(je,ke,is)
                                           +psi[1]*omg[1]*phi[1]*x(je,ke,ie);
                                y7(j,k,i) = psi[2]*omg[2]*phi[2]*y(js,ks,is)
                                           +psi[2]*omg[2]*phi[1]*y(js,ks,ie)
                                           +psi[2]*omg[1]*phi[2]*y(js,ke,is)
                                           +psi[2]*omg[1]*phi[1]*y(js,ke,ie)
                                           +psi[1]*omg[2]*phi[2]*y(je,ks,is)
                                           +psi[1]*omg[2]*phi[1]*y(je,ks,ie)
                                           +psi[1]*omg[1]*phi[2]*y(je,ke,is)
                                           +psi[1]*omg[1]*phi[1]*y(je,ke,ie);
                                z7(j,k,i) = psi[2]*omg[2]*phi[2]*z(js,ks,is)
                                           +psi[2]*omg[2]*phi[1]*z(js,ks,ie)
                                           +psi[2]*omg[1]*phi[2]*z(js,ke,is)
                                           +psi[2]*omg[1]*phi[1]*z(js,ke,ie)
                                           +psi[1]*omg[2]*phi[2]*z(je,ks,is)
                                           +psi[1]*omg[2]*phi[1]*z(je,ks,ie)
                                           +psi[1]*omg[1]*phi[2]*z(je,ke,is)
                                           +psi[1]*omg[1]*phi[1]*z(je,ke,ie);
                            }
                        }
                    }

                    for (i = is+1; i <= ie-1; i++) {
                        for (j = js+1; j <= je-1; j++) {
                            for (k = ks+1; k <= ke-1; k++) {
                                x(j,k,i) = x1(j,k,i) + x2(j,k,i) + x3(j,k,i)
                                          - x4(j,k,i) - x5(j,k,i) - x6(j,k,i)
                                          + x7(j,k,i);
                                y(j,k,i) = y1(j,k,i) + y2(j,k,i) + y3(j,k,i)
                                          - y4(j,k,i) - y5(j,k,i) - y6(j,k,i)
                                          + y7(j,k,i);
                                z(j,k,i) = z1(j,k,i) + z2(j,k,i) + z3(j,k,i)
                                          - z4(j,k,i) - z5(j,k,i) - z6(j,k,i)
                                          + z7(j,k,i);
                            }
                        }
                    }

                } // k1
            } // j1
        } // i1
    } // if abs(isktyp)==1

    return;
}

} // namespace tfivol_ns
