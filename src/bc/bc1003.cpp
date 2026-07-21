// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc1003.h"
#include <cmath>
#include <algorithm>

namespace bc1003_ns {


// rie1d: 1D Riemann characteristic inflow/outflow boundary condition solver
// t array layout:
//   t(jv, 1..5) = interior q state: rho, u, v, w, p
//   t(jv, 6..8) = outward normal vector components (snx, sny, snz)
//   t(jv, 20)   = face area (magnitude of normal)
//   t(jv, 18..19) = x,z coords (for iipv=1 case, not used in rie1d)
//   t(jv, 21..23) = scratch (reserved for rie1d)
void rie1d(int& jvdim, FortranArray2DRef<double> t, int& jv, double& cl)
{
    double p0   = (double)cmn_ivals.p0;
    double rho0 = (double)cmn_ivals.rho0;
    double c0   = (double)cmn_ivals.c0;
    double u0   = (double)cmn_ivals.u0;
    double v0   = (double)cmn_ivals.v0;
    double w0   = (double)cmn_ivals.w0;

    const double gam = 1.4;
    const double gm1 = 0.4;

    for (int jv_ = 1; jv_ <= jv; jv_++) {
        double rho = t(jv_, 1);
        double u   = t(jv_, 2);
        double v   = t(jv_, 3);
        double w   = t(jv_, 4);
        double p   = t(jv_, 5);
        double snx = t(jv_, 6);
        double sny = t(jv_, 7);
        double snz = t(jv_, 8);

        double ds = std::sqrt(snx*snx + sny*sny + snz*snz);
        double nx = snx / ds;
        double ny = sny / ds;
        double nz = snz / ds;

        double ci  = std::sqrt(gam * p / rho);
        double qni = u*nx + v*ny + w*nz;
        double qne = u0*nx + v0*ny + w0*nz;

        double rplus  = qni + 2.0*ci / gm1;
        double rminus = qne - 2.0*c0 / gm1;

        double qnb = 0.5*(rplus + rminus);
        double cb  = 0.25*gm1*(rplus - rminus);
        if (cb < 0.0) cb = 0.0;

        double rho_b, u_b, v_b, w_b, p_b;
        if (qnb >= 0.0) {
            // outflow: entropy from interior
            double s_int = p / std::pow(rho, gam);
            rho_b = std::pow(cb*cb / (gam * s_int), 1.0/gm1);
            p_b   = rho_b * cb*cb / gam;
            u_b   = (u  - qni*nx) + qnb*nx;
            v_b   = (v  - qni*ny) + qnb*ny;
            w_b   = (w  - qni*nz) + qnb*nz;
        } else {
            // inflow: entropy from freestream
            double s_ext = p0 / std::pow(rho0, gam);
            rho_b = std::pow(cb*cb / (gam * s_ext), 1.0/gm1);
            p_b   = rho_b * cb*cb / gam;
            u_b   = (u0 - qne*nx) + qnb*nx;
            v_b   = (v0 - qne*ny) + qnb*ny;
            w_b   = (w0 - qne*nz) + qnb*nz;
        }

        t(jv_, 1) = rho_b;
        t(jv_, 2) = u_b;
        t(jv_, 3) = v_b;
        t(jv_, 4) = w_b;
        t(jv_, 5) = p_b;
    }
}


// bc1003: characteristic inflow/outflow boundary conditions (bctype 1003)
void bc1003(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0,
            FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj,
            FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            FortranArray2DRef<double> t,
            int& jvdim, int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0,
            FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0,
            FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& iuns,
            FortranArray3DRef<double> x,
            FortranArray3DRef<double> z,
            double& cl,
            FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim,
            int& myid, int& nummem)
{
    // COMMON block aliases
    int& level   = cmn_mgrd.level;
    int& lglobal = cmn_mgrd.lglobal;
    int* ivisc   = cmn_reyue.ivisc;   // 0-based
    float* tur10 = cmn_ivals.tur10;   // 0-based
    int& iipv    = cmn_info.iipv;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    double ubar;
    double xi_, zi_, xo_, zo_;


    // **************************************************************************
    //      j=1 boundary                inflow/outflow                   bctype 1003
    // **************************************************************************
    if (nface == 3) {
        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista)*(kend - ksta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    t(js + k - ksta, l) = q(1, k, i, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    t(js + k - ksta, 5 + l) = -sj(1, k, i, l);
                }
            }
            for (int k = ksta; k <= kend1; k++) {
                t(js + k - ksta, 20) = -sj(1, k, i, 5);
            }
            if (iipv == 1) {
                for (int k = ksta; k <= kend1; k++) {
                    xi_ = 0.25*(x(1,k,i)+x(1,k+1,i) + x(2,k,i)+x(2,k+1,i));
                    zi_ = 0.25*(z(1,k,i)+z(1,k+1,i) + z(2,k,i)+z(2,k+1,i));
                    xo_ = 0.5*(x(1,k,i)+x(1,k+1,i));
                    zo_ = 0.5*(z(1,k,i)+z(1,k+1,i));
                    t(js + k - ksta, 18) = xo_ + (xo_ - xi_);
                    t(js + k - ksta, 19) = zo_ + (zo_ - zi_);
                }
            }
        }
        int jv = (kend - ksta)*(iend - ista);
        rie1d(jvdim, t, jv, cl);
        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista)*(kend - ksta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    qj0(k, i, l, 1) = t(k - ksta + js, l);
                    qj0(k, i, l, 2) = qj0(k, i, l, 1);
                    bcj(k, i, 1)    = 0.0;
                }
            }
        }
        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k, i, 1, 1) = vist3d(1, k, i);
                    vj0(k, i, 1, 2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int k = ksta; k <= kend1; k++) {
                            ubar = -(qj0(k,i,2,1)*sj(1,k,i,1) +
                                     qj0(k,i,3,1)*sj(1,k,i,2) +
                                     qj0(k,i,4,1)*sj(1,k,i,3));
                            if ((float)ubar < 0.f) {
                                tj0(k, i, l, 1) = (double)tur10[l-1];
                                tj0(k, i, l, 2) = (double)tur10[l-1];
                            } else {
                                tj0(k, i, l, 1) = tursav(1, k, i, l);
                                tj0(k, i, l, 2) = tj0(k, i, l, 1);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==3


    // **************************************************************************
    //      j=jdim boundary             inflow/outflow                   bctype 1003
    // **************************************************************************
    if (nface == 4) {
        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista)*(kend - ksta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    t(js + k - ksta, l) = q(jdim1, k, i, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    t(js + k - ksta, 5 + l) = sj(jdim, k, i, l);
                }
            }
            for (int k = ksta; k <= kend1; k++) {
                t(js + k - ksta, 20) = sj(jdim, k, i, 5);
            }
            if (iipv == 1) {
                for (int k = ksta; k <= kend1; k++) {
                    xi_ = 0.25*(x(jdim ,k,i)+x(jdim,k+1,i)
                              + x(jdim1,k,i)+x(jdim1,k+1,i));
                    zi_ = 0.25*(z(jdim ,k,i)+z(jdim ,k+1,i)
                              + z(jdim1,k,i)+z(jdim1,k+1,i));
                    xo_ = 0.5*(x(jdim,k,i)+x(jdim,k+1,i));
                    zo_ = 0.5*(z(jdim,k,i)+z(jdim,k+1,i));
                    t(js + k - ksta, 18) = xo_ + (xo_ - xi_);
                    t(js + k - ksta, 19) = zo_ + (zo_ - zi_);
                }
            }
        }
        int jv = (kend - ksta)*(iend - ista);
        rie1d(jvdim, t, jv, cl);
        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista)*(kend - ksta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int k = ksta; k <= kend1; k++) {
                    qj0(k, i, l, 3) = t(k - ksta + js, l);
                    qj0(k, i, l, 4) = qj0(k, i, l, 3);
                    bcj(k, i, 2)    = 0.0;
                }
            }
        }
        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k, i, 1, 3) = vist3d(jdim-1, k, i);
                    vj0(k, i, 1, 4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int k = ksta; k <= kend1; k++) {
                            ubar = qj0(k,i,2,3)*sj(jdim,k,i,1) +
                                   qj0(k,i,3,3)*sj(jdim,k,i,2) +
                                   qj0(k,i,4,3)*sj(jdim,k,i,3);
                            if ((float)ubar < 0.f) {
                                tj0(k, i, l, 3) = (double)tur10[l-1];
                                tj0(k, i, l, 4) = (double)tur10[l-1];
                            } else {
                                tj0(k, i, l, 3) = tursav(jdim-1, k, i, l);
                                tj0(k, i, l, 4) = tj0(k, i, l, 3);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==4


    // **************************************************************************
    //      k=1 boundary                inflow/outflow                   bctype 1003
    // **************************************************************************
    if (nface == 5) {
        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista)*(jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    t(js + jj - jsta, l) = q(jj, 1, i, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    t(js + jj - jsta, 5 + l) = -sk(jj, 1, i, l);
                }
            }
            for (int jj = jsta; jj <= jend1; jj++) {
                t(js + jj - jsta, 20) = -sk(jj, 1, i, 5);
            }
            if (iipv == 1) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    xi_ = 0.25*(x(jj,1,i)+x(jj+1,1,i) + x(jj,2,i)+x(jj+1,2,i));
                    zi_ = 0.25*(z(jj,1,i)+z(jj+1,1,i) + z(jj,2,i)+z(jj+1,2,i));
                    xo_ = 0.5*(x(jj,1,i)+x(jj+1,1,i));
                    zo_ = 0.5*(z(jj,1,i)+z(jj+1,1,i));
                    t(js + jj - jsta, 18) = xo_ + (xo_ - xi_);
                    t(js + jj - jsta, 19) = zo_ + (zo_ - zi_);
                }
            }
        }
        int jv = (jend - jsta)*(iend - ista);
        rie1d(jvdim, t, jv, cl);
        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista)*(jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    qk0(j, i, l, 1) = t(j - jsta + js, l);
                    qk0(j, i, l, 2) = qk0(j, i, l, 1);
                    bck(j, i, 1)    = 0.0;
                }
            }
        }
        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j, i, 1, 1) = vist3d(j, 1, i);
                    vk0(j, i, 1, 2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int j = jsta; j <= jend1; j++) {
                            ubar = -(qk0(j,i,2,1)*sk(j,1,i,1) +
                                     qk0(j,i,3,1)*sk(j,1,i,2) +
                                     qk0(j,i,4,1)*sk(j,1,i,3));
                            if ((float)ubar < 0.f) {
                                tk0(j, i, l, 1) = (double)tur10[l-1];
                                tk0(j, i, l, 2) = (double)tur10[l-1];
                            } else {
                                tk0(j, i, l, 1) = tursav(j, 1, i, l);
                                tk0(j, i, l, 2) = tk0(j, i, l, 1);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==5


    // **************************************************************************
    //      k=kdim boundary             inflow/outflow                   bctype 1003
    // **************************************************************************
    if (nface == 6) {
        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista)*(jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    t(js + jj - jsta, l) = q(jj, kdim1, i, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    t(js + jj - jsta, 5 + l) = sk(jj, kdim, i, l);
                }
            }
            for (int jj = jsta; jj <= jend1; jj++) {
                t(js + jj - jsta, 20) = sk(jj, kdim, i, 5);
            }
            if (iipv == 1) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    xi_ = 0.25*(x(jj,kdim ,i)+x(jj+1,kdim ,i)
                              + x(jj,kdim1,i)+x(jj+1,kdim1,i));
                    zi_ = 0.25*(z(jj,kdim ,i)+z(jj+1,kdim ,i)
                              + z(jj,kdim1,i)+z(jj+1,kdim1,i));
                    xo_ = 0.5*(x(jj,kdim,i)+x(jj+1,kdim,i));
                    zo_ = 0.5*(z(jj,kdim,i)+z(jj+1,kdim,i));
                    t(js + jj - jsta, 18) = xo_ + (xo_ - xi_);
                    t(js + jj - jsta, 19) = zo_ + (zo_ - zi_);
                }
            }
        }
        int jv = (jend - jsta)*(iend - ista);
        rie1d(jvdim, t, jv, cl);
        for (int i = ista; i <= iend1; i++) {
            int js = (i - ista)*(jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    qk0(j, i, l, 3) = t(j - jsta + js, l);
                    qk0(j, i, l, 4) = qk0(j, i, l, 3);
                    bck(j, i, 2)    = 0.0;
                }
            }
        }
        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j, i, 1, 3) = vist3d(j, kdim-1, i);
                    vk0(j, i, 1, 4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        for (int j = jsta; j <= jend1; j++) {
                            ubar = qk0(j,i,2,3)*sk(j,kdim,i,1) +
                                   qk0(j,i,3,3)*sk(j,kdim,i,2) +
                                   qk0(j,i,4,3)*sk(j,kdim,i,3);
                            if ((float)ubar < 0.f) {
                                tk0(j, i, l, 3) = (double)tur10[l-1];
                                tk0(j, i, l, 4) = (double)tur10[l-1];
                            } else {
                                tk0(j, i, l, 3) = tursav(j, kdim-1, i, l);
                                tk0(j, i, l, 4) = tk0(j, i, l, 3);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==6


    // **************************************************************************
    //      i=1 boundary                inflow/outflow                   bctype 1003
    // **************************************************************************
    if (nface == 1) {
        for (int k = ksta; k <= kend1; k++) {
            int js = (k - ksta)*(jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    t(js + jj - jsta, l) = q(jj, k, 1, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    t(js + jj - jsta, 5 + l) = -si(jj, k, 1, l);
                }
            }
            for (int jj = jsta; jj <= jend1; jj++) {
                t(js + jj - jsta, 20) = -si(jj, k, 1, 5);
            }
        }
        int jv = (jend - jsta)*(kend - ksta);
        rie1d(jvdim, t, jv, cl);
        for (int k = ksta; k <= kend1; k++) {
            int js = (k - ksta)*(jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    qi0(j, k, l, 1) = t(j - jsta + js, l);
                    qi0(j, k, l, 2) = qi0(j, k, l, 1);
                    bci(j, k, 1)    = 0.0;
                }
            }
        }
        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j, k, 1, 1) = vist3d(j, k, 1);
                    vi0(j, k, 1, 2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        for (int j = jsta; j <= jend1; j++) {
                            ubar = -(qi0(j,k,2,1)*si(j,k,1,1) +
                                     qi0(j,k,3,1)*si(j,k,1,2) +
                                     qi0(j,k,4,1)*si(j,k,1,3));
                            if ((float)ubar < 0.f) {
                                ti0(j, k, l, 1) = (double)tur10[l-1];
                                ti0(j, k, l, 2) = (double)tur10[l-1];
                            } else {
                                ti0(j, k, l, 1) = tursav(j, k, 1, l);
                                ti0(j, k, l, 2) = ti0(j, k, l, 1);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==1


    // **************************************************************************
    //      i=idim boundary             inflow/outflow                   bctype 1003
    // **************************************************************************
    if (nface == 2) {
        for (int k = ksta; k <= kend1; k++) {
            int js = (k - ksta)*(jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    t(js + jj - jsta, l) = q(jj, k, idim1, l);
                }
            }
            for (int l = 1; l <= 3; l++) {
                for (int jj = jsta; jj <= jend1; jj++) {
                    t(js + jj - jsta, 5 + l) = si(jj, k, idim, l);
                }
            }
            for (int jj = jsta; jj <= jend1; jj++) {
                t(js + jj - jsta, 20) = si(jj, k, idim, 5);
            }
        }
        int jv = (jend - jsta)*(kend - ksta);
        rie1d(jvdim, t, jv, cl);
        for (int k = ksta; k <= kend1; k++) {
            int js = (k - ksta)*(jend - jsta) + 1;
            for (int l = 1; l <= 5; l++) {
                for (int j = jsta; j <= jend1; j++) {
                    qi0(j, k, l, 3) = t(j - jsta + js, l);
                    qi0(j, k, l, 4) = qi0(j, k, l, 3);
                    bci(j, k, 2)    = 0.0;
                }
            }
        }
        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j, k, 1, 3) = vist3d(j, k, idim-1);
                    vi0(j, k, 1, 4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (int l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        for (int j = jsta; j <= jend1; j++) {
                            ubar = qi0(j,k,2,3)*si(j,k,idim,1) +
                                   qi0(j,k,3,3)*si(j,k,idim,2) +
                                   qi0(j,k,4,3)*si(j,k,idim,3);
                            if ((float)ubar < 0.f) {
                                ti0(j, k, l, 3) = (double)tur10[l-1];
                                ti0(j, k, l, 4) = (double)tur10[l-1];
                            } else {
                                ti0(j, k, l, 3) = tursav(j, k, idim-1, l);
                                ti0(j, k, l, 4) = ti0(j, k, l, 3);
                            }
                        }
                    }
                }
            }
        }
    } // end nface==2
} // end bc1003


// bc: dispatcher for bctype 1003 boundary conditions
void bc(int& ntime, int& nbl,
        FortranArray2DRef<int> lw,
        FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w,
        int& mgwk,
        FortranArray1DRef<double> wk,
        int& nwork,
        double& cl,
        FortranArray1DRef<int> nou,
        FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim,
        int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans,
        FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm,
        FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg,
        FortranArray1DRef<int> nbci0,
        FortranArray1DRef<int> nbcj0,
        FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim,
        FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim,
        FortranArray4DRef<int> ibcinfo,
        FortranArray4DRef<int> jbcinfo,
        FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei,
        FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek,
        FortranArray3DRef<int> lwdat,
        int& myid,
        FortranArray1DRef<int> idimg,
        FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg,
        FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil,
        int& nummem)
{
    int iuns = 0;

    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        // Build array views from w (lw offsets)
        FortranArray4DRef<double> q    (&w(lw(1,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray4DRef<double> qj0  (&w(lw(2,ibl)),  kdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qk0  (&w(lw(3,ibl)),  jdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qi0  (&w(lw(4,ibl)),  jdim,   kdim,   5,      4);
        FortranArray4DRef<double> sj   (&w(lw(5,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> sk   (&w(lw(6,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> si   (&w(lw(7,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray3DRef<double> bcj  (&w(lw(8,ibl)),  kdim,   idim-1, 2);
        FortranArray3DRef<double> bck  (&w(lw(9,ibl)),  jdim,   idim-1, 2);
        FortranArray3DRef<double> bci  (&w(lw(10,ibl)), jdim,   kdim,   2);
        // lw(11..16) = xtbj, xtbk, xtbi, atbj, atbk, atbi (not needed for bc1003)
        // lw(17) = x, lw(18) = y, lw(19) = z
        FortranArray3DRef<double> x    (&w(lw(17,ibl)), jdim,   kdim,   idim);
        FortranArray3DRef<double> z    (&w(lw(19,ibl)), jdim,   kdim,   idim);

        // Build array views from wk (lw2 offsets)
        FortranArray4DRef<double> tursav(&wk(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0   (&wk(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0   (&wk(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0   (&wk(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&wk(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0   (&wk(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0   (&wk(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0   (&wk(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // Compute jvdim: maximum face size for t work array
        int jvdim = std::max(jdim * kdim,
                    std::max(jdim * (idim-1),
                             kdim * (idim-1)));
        FortranArray2D<double> t_arr(jvdim, 23);
        FortranArray2DRef<double> t = t_arr.ref();

        // Loop over j-face boundary segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl, iseg, 1, 1) == 1003) {
                int ista  = jbcinfo(ibl, iseg, 1, 2);
                int iend  = jbcinfo(ibl, iseg, 1, 3);
                int jsta  = jbcinfo(ibl, iseg, 1, 4);
                int jend  = jbcinfo(ibl, iseg, 1, 5);
                int ksta  = jbcinfo(ibl, iseg, 1, 6);
                int kend  = jbcinfo(ibl, iseg, 1, 7);
                int nface = jbcinfo(ibl, iseg, 1, 8);
                bc1003(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci,
                       ista, iend, jsta, jend, ksta, kend,
                       t, jvdim, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, x, z, cl, nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // Loop over k-face boundary segments
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl, iseg, 1, 1) == 1003) {
                int ista  = kbcinfo(ibl, iseg, 1, 2);
                int iend  = kbcinfo(ibl, iseg, 1, 3);
                int jsta  = kbcinfo(ibl, iseg, 1, 4);
                int jend  = kbcinfo(ibl, iseg, 1, 5);
                int ksta  = kbcinfo(ibl, iseg, 1, 6);
                int kend  = kbcinfo(ibl, iseg, 1, 7);
                int nface = kbcinfo(ibl, iseg, 1, 8);
                bc1003(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci,
                       ista, iend, jsta, jend, ksta, kend,
                       t, jvdim, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, x, z, cl, nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // Loop over i-face boundary segments
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl, iseg, 1, 1) == 1003) {
                int ista  = ibcinfo(ibl, iseg, 1, 2);
                int iend  = ibcinfo(ibl, iseg, 1, 3);
                int jsta  = ibcinfo(ibl, iseg, 1, 4);
                int jend  = ibcinfo(ibl, iseg, 1, 5);
                int ksta  = ibcinfo(ibl, iseg, 1, 6);
                int kend  = ibcinfo(ibl, iseg, 1, 7);
                int nface = ibcinfo(ibl, iseg, 1, 8);
                bc1003(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci,
                       ista, iend, jsta, jend, ksta, kend,
                       t, jvdim, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       iuns, x, z, cl, nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }
    } // end ibl loop
}

} // namespace bc1003_ns
