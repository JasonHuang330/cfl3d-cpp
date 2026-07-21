// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rotate.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace rotate_ns {

void rotate(int& jdim, int& kdim, int& idim,
            FortranArray2DRef<double> t,
            FortranArray3DRef<double> tti,
            FortranArray3DRef<double> ttj,
            FortranArray3DRef<double> ttk,
            FortranArray2DRef<double> t1,
            FortranArray3DRef<double> x,
            FortranArray3DRef<double> y,
            FortranArray3DRef<double> z,
            int& nbl, int& irot,
            double& rfreqr,
            double& omegx, double& omegy, double& omegz,
            double& xorg, double& yorg, double& zorg,
            double& thetax, double& thetay, double& thetaz,
            double& thxold, double& thyold, double& thzold,
            int& iupdat, double& time2,
            FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim)
{
    // COMMON /sklton/
    int& isklton = cmn_sklton.isklton;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    //
    // ft modulates the rotation
    // dfdt is the time derivative of ft
    // d2fdt2 is the second time derivative of ft
    //
    double ft = 0.0, dfdt = 0.0, d2fdt2 = 0.0;
    double expt;

    if (irot == 0) {
        return;
    } else if (irot == 1) {
        ft     = time2;
        dfdt   = 1.0;
        d2fdt2 = 0.0;
    } else if (irot == 2) {
        ft     = std::sin(rfreqr * time2);
        dfdt   = rfreqr * std::cos(rfreqr * time2);
        d2fdt2 = -(rfreqr) * (rfreqr) * std::sin(rfreqr * time2);
    } else if (irot == 3) {
        expt   = std::exp(-rfreqr * time2);
        ft     = 1.0 - expt;
        dfdt   = rfreqr * expt;
        d2fdt2 = -(rfreqr) * (rfreqr) * expt;
    } else if (irot == 99) {
        ft     = 0.0;
        dfdt   = 1.0;
        d2fdt2 = 0.0;
    }

    // Raw pointers for flat access to x, y, z (column-major: x(j,k,i))
    // Flat 0-based index for x(izz,1,i): (izz-1) + (i-1)*jdim*kdim
    double* xptr = &x(1, 1, 1);
    double* yptr = &y(1, 1, 1);
    double* zptr = &z(1, 1, 1);


    if (std::abs((float)omegx) > 0.0f) {
        //
        // rotate about an axis parallel to the x-axis
        //
        double theold;
        if (time2 != 0.0) {
            theold = thxold;
        } else {
            theold = 0.0e0;
        }
        //
        // calculate rotated y and z at grid points
        // t1(1)=x (unaltered) t1(2)=y t1(3)=z
        //
        double theta, dthedt, d2thedt2;
        if (irot != 99) {
            theta = omegx * ft;
        } else {
            theta = thetax;
        }
        dthedt   = omegx * dfdt;
        d2thedt2 = omegx * d2fdt2;
        double dtheta = theta - theold;
        double ca = std::cos(dtheta);
        double sa = std::sin(dtheta);
        int n = jdim * kdim;
        for (int i = 1; i <= idim; i++) {
            int js = jdim * kdim * (i - 1) + 1;
            for (int izz = 1; izz <= n; izz++) {
                int flat = izz - 1 + (i - 1) * jdim * kdim;
                t1(izz + js - 1, 1) =   xptr[flat];
                t1(izz + js - 1, 2) =  (yptr[flat] - yorg) * ca
                                       - (zptr[flat] - zorg) * sa + yorg;
                t1(izz + js - 1, 3) =  (yptr[flat] - yorg) * sa
                                       + (zptr[flat] - zorg) * ca + zorg;
            }
        }
        //
        // update grid to new position
        //
        if (iupdat > 0) {
            thetax = theta;
            for (int i = 1; i <= idim; i++) {
                int js = jdim * kdim * (i - 1) + 1;
                for (int izz = 1; izz <= n; izz++) {
                    int flat = izz - 1 + (i - 1) * jdim * kdim;
                    xptr[flat] = t1(izz + js - 1, 1);
                    yptr[flat] = t1(izz + js - 1, 2);
                    zptr[flat] = t1(izz + js - 1, 3);
                }
            }
        }
        //
        // calculate increment to speed of grid points
        // t(1)=dx/dt t(2)=dy/dt t(3)=dz/dt
        //
        for (int i = 1; i <= idim; i++) {
            int js = jdim * kdim * (i - 1) + 1;
            for (int izz = 1; izz <= n; izz++) {
                int flat = izz - 1 + (i - 1) * jdim * kdim;
                t(izz + js - 1, 1) = t(izz + js - 1, 1) + 0.0e0;
                t(izz + js - 1, 2) = t(izz + js - 1, 2) - (zptr[flat] - zorg) * dthedt;
                t(izz + js - 1, 3) = t(izz + js - 1, 3) + (yptr[flat] - yorg) * dthedt;
            }
        }
        // i0/idim boundaries
        n = jdim * kdim;
        { int i_var = 1;
          for (int ii = 1; ii <= 2; ii++) {
            for (int izz = 1; izz <= n; izz++) {
                tti(izz, 1, ii) = tti(izz, 1, ii) + 0.0e0;
                tti(izz, 2, ii) = tti(izz, 2, ii) + 0.0e0;
                tti(izz, 3, ii) = tti(izz, 3, ii) + 0.0e0;
            }
            i_var = i_var + idim - 1;
          }
        }
        // j0/jdim boundaries
        n = idim * kdim;
        { int j_var = 1;
          for (int jj = 1; jj <= 2; jj++) {
            for (int izz = 1; izz <= n; izz++) {
                ttj(izz, 1, jj) = ttj(izz, 1, jj) + 0.0e0;
                ttj(izz, 2, jj) = ttj(izz, 2, jj) + 0.0e0;
                ttj(izz, 3, jj) = ttj(izz, 3, jj) + 0.0e0;
            }
            j_var = j_var + jdim - 1;
          }
        }
        // k0/kdim boundaries
        n = jdim * idim;
        { int k_var = 1;
          for (int kk = 1; kk <= 2; kk++) {
            for (int izz = 1; izz <= n; izz++) {
                ttk(izz, 1, kk) = ttk(izz, 1, kk) + 0.0e0;
                ttk(izz, 2, kk) = ttk(izz, 2, kk) + 0.0e0;
                ttk(izz, 3, kk) = ttk(izz, 3, kk) + 0.0e0;
            }
            k_var = k_var + kdim - 1;
          }
        }

    } else if (std::abs((float)omegy) > 0.0f) {
        //
        // rotate about an axis parallel to the y-axis
        //
        double theold;
        if (time2 != 0.0) {
            theold = thyold;
        } else {
            theold = 0.0e0;
        }
        //
        // calculate rotated x and z at grid points
        // t1(1)=x t1(2)=y (unaltered) t1(3)=z
        //
        double theta, dthedt, d2thedt2;
        if (irot != 99) {
            theta = omegy * ft;
        } else {
            theta = thetay;
        }
        dthedt   = omegy * dfdt;
        d2thedt2 = omegy * d2fdt2;
        double dtheta = theta - theold;
        double ca = std::cos(dtheta);
        double sa = std::sin(dtheta);
        int n = jdim * kdim;
        for (int i = 1; i <= idim; i++) {
            int js = jdim * kdim * (i - 1) + 1;
            for (int izz = 1; izz <= n; izz++) {
                int flat = izz - 1 + (i - 1) * jdim * kdim;
                t1(izz + js - 1, 1) =  (xptr[flat] - xorg) * ca
                                       + (zptr[flat] - zorg) * sa + xorg;
                t1(izz + js - 1, 2) =   yptr[flat];
                t1(izz + js - 1, 3) = -(xptr[flat] - xorg) * sa
                                       + (zptr[flat] - zorg) * ca + zorg;
            }
        }
        //
        // update grid to new position
        //
        if (iupdat > 0) {
            thetay = theta;
            for (int i = 1; i <= idim; i++) {
                int js = jdim * kdim * (i - 1) + 1;
                for (int izz = 1; izz <= n; izz++) {
                    int flat = izz - 1 + (i - 1) * jdim * kdim;
                    xptr[flat] = t1(izz + js - 1, 1);
                    yptr[flat] = t1(izz + js - 1, 2);
                    zptr[flat] = t1(izz + js - 1, 3);
                }
            }
        }
        //
        // calculate increment to speed of grid points
        // t(1)=dx/dt t(2)=dy/dt t(3)=dz/dt
        //
        for (int i = 1; i <= idim; i++) {
            int js = jdim * kdim * (i - 1) + 1;
            for (int izz = 1; izz <= n; izz++) {
                int flat = izz - 1 + (i - 1) * jdim * kdim;
                t(izz + js - 1, 1) = t(izz + js - 1, 1) + (zptr[flat] - zorg) * dthedt;
                t(izz + js - 1, 2) = t(izz + js - 1, 2) + 0.0e0;
                t(izz + js - 1, 3) = t(izz + js - 1, 3) - (xptr[flat] - xorg) * dthedt;
            }
        }
        // i0/idim boundaries
        n = jdim * kdim;
        { int i_var = 1;
          for (int ii = 1; ii <= 2; ii++) {
            for (int izz = 1; izz <= n; izz++) {
                tti(izz, 1, ii) = tti(izz, 1, ii) + 0.0e0;
                tti(izz, 2, ii) = tti(izz, 2, ii) + 0.0e0;
                tti(izz, 3, ii) = tti(izz, 3, ii) + 0.0e0;
            }
            i_var = i_var + idim - 1;
          }
        }
        // j0/jdim boundaries
        n = idim * kdim;
        { int j_var = 1;
          for (int jj = 1; jj <= 2; jj++) {
            for (int izz = 1; izz <= n; izz++) {
                ttj(izz, 1, jj) = ttj(izz, 1, jj) + 0.0e0;
                ttj(izz, 2, jj) = ttj(izz, 2, jj) + 0.0e0;
                ttj(izz, 3, jj) = ttj(izz, 3, jj) + 0.0e0;
            }
            j_var = j_var + jdim - 1;
          }
        }
        // k0/kdim boundaries
        n = jdim * idim;
        { int k_var = 1;
          for (int kk = 1; kk <= 2; kk++) {
            for (int izz = 1; izz <= n; izz++) {
                ttk(izz, 1, kk) = ttk(izz, 1, kk) + 0.0e0;
                ttk(izz, 2, kk) = ttk(izz, 2, kk) + 0.0e0;
                ttk(izz, 3, kk) = ttk(izz, 3, kk) + 0.0e0;
            }
            k_var = k_var + kdim - 1;
          }
        }

    } else if (std::abs((float)omegz) > 0.0f) {
        //
        // rotate about an axis parallel to the z-axis
        //
        double theold;
        if (time2 != 0.0) {
            theold = thzold;
        } else {
            theold = 0.0e0;
        }
        //
        // calculate rotated x and y at grid points
        // t1(1)=x t1(2)=y t1(3)=z (unaltered)
        //
        double theta, dthedt, d2thedt2;
        if (irot != 99) {
            theta = omegz * ft;
        } else {
            theta = thetaz;
        }
        dthedt   = omegz * dfdt;
        d2thedt2 = omegz * d2fdt2;
        double dtheta = theta - theold;
        double ca = std::cos(dtheta);
        double sa = std::sin(dtheta);
        int n = jdim * kdim;
        for (int i = 1; i <= idim; i++) {
            int js = jdim * kdim * (i - 1) + 1;
            for (int izz = 1; izz <= n; izz++) {
                int flat = izz - 1 + (i - 1) * jdim * kdim;
                t1(izz + js - 1, 1) =  (xptr[flat] - xorg) * ca
                                       - (yptr[flat] - yorg) * sa + xorg;
                t1(izz + js - 1, 2) =  (xptr[flat] - xorg) * sa
                                       + (yptr[flat] - yorg) * ca + yorg;
                t1(izz + js - 1, 3) =   zptr[flat];
            }
        }
        //
        // update grid to new position
        //
        if (iupdat > 0) {
            thetaz = theta;
            for (int i = 1; i <= idim; i++) {
                int js = jdim * kdim * (i - 1) + 1;
                for (int izz = 1; izz <= n; izz++) {
                    int flat = izz - 1 + (i - 1) * jdim * kdim;
                    xptr[flat] = t1(izz + js - 1, 1);
                    yptr[flat] = t1(izz + js - 1, 2);
                    zptr[flat] = t1(izz + js - 1, 3);
                }
            }
        }
        //
        // calculate increment to speed of grid points
        // t(1)=dx/dt t(2)=dy/dt t(3)=dz/dt
        //
        for (int i = 1; i <= idim; i++) {
            int js = jdim * kdim * (i - 1) + 1;
            for (int izz = 1; izz <= n; izz++) {
                int flat = izz - 1 + (i - 1) * jdim * kdim;
                t(izz + js - 1, 1) = t(izz + js - 1, 1) - (yptr[flat] - yorg) * dthedt;
                t(izz + js - 1, 2) = t(izz + js - 1, 2) + (xptr[flat] - xorg) * dthedt;
                t(izz + js - 1, 3) = t(izz + js - 1, 3) + 0.0e0;
            }
        }
        // i0/idim boundaries
        n = jdim * kdim;
        { int i_var = 1;
          for (int ii = 1; ii <= 2; ii++) {
            for (int izz = 1; izz <= n; izz++) {
                tti(izz, 1, ii) = tti(izz, 1, ii) + 0.0e0;
                tti(izz, 2, ii) = tti(izz, 2, ii) + 0.0e0;
                tti(izz, 3, ii) = tti(izz, 3, ii) + 0.0e0;
            }
            i_var = i_var + idim - 1;
          }
        }
        // j0/jdim boundaries
        n = idim * kdim;
        { int j_var = 1;
          for (int jj = 1; jj <= 2; jj++) {
            for (int izz = 1; izz <= n; izz++) {
                ttj(izz, 1, jj) = ttj(izz, 1, jj) + 0.0e0;
                ttj(izz, 2, jj) = ttj(izz, 2, jj) + 0.0e0;
                ttj(izz, 3, jj) = ttj(izz, 3, jj) + 0.0e0;
            }
            j_var = j_var + jdim - 1;
          }
        }
        // k0/kdim boundaries
        n = jdim * idim;
        { int k_var = 1;
          for (int kk = 1; kk <= 2; kk++) {
            for (int izz = 1; izz <= n; izz++) {
                ttk(izz, 1, kk) = ttk(izz, 1, kk) + 0.0e0;
                ttk(izz, 2, kk) = ttk(izz, 2, kk) + 0.0e0;
                ttk(izz, 3, kk) = ttk(izz, 3, kk) + 0.0e0;
            }
            k_var = k_var + kdim - 1;
          }
        }

    } else {
        if (isklton > 0) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            // format 101: "WARNING: this block has zero rotational displacement"
            std::snprintf(bou(nou(1), 1), 120,
                          " WARNING: this block has zero rotational"
                          " displacement");
        }
    }

    return;
} // end rotate

} // namespace rotate_ns
