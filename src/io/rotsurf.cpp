// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rotsurf.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace rotsurf_ns {

void rotsurf(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti, int& nbl, int& idef, double& xorg, double& yorg, double& zorg, double& omegx, double& omegy, double& omegz, double& thetax, double& thetay, double& thetaz, double& rfreqr, int& ici, int& icf, int& jci, int& jcf, int& kci, int& kcf, double& time, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray3DRef<double> wkj, FortranArray3DRef<double> wkk, FortranArray3DRef<double> wki)
{
    int& isklton = cmn_sklton.isklton;

    // determine which block face is to be rotated
    int isurf = 0;
    int ii = 0, jj = 0, kk = 0, ll = 0;

    if (ici == icf) {
        isurf = 1;
        ii    = ici;
        ll    = 1;
        if (ici == idim) ll = 2;
    } else if (jci == jcf) {
        isurf = 2;
        jj    = jci;
        ll    = 1;
        if (jci == jdim) ll = 2;
    } else if (kci == kcf) {
        isurf = 3;
        kk    = kci;
        ll    = 1;
        if (kci == kdim) ll = 2;
    }

    // ft modulates the rotation
    // dfdt is the time derivative of ft
    // d2fdt2 is the second time derivative of ft

    double ft = 0.0, dfdt = 0.0, d2fdt2 = 0.0;

    if (idef == 0) {
        return;
    } else if (idef == 2) {
        ft     = std::sin(rfreqr * time);
        dfdt   = rfreqr * std::cos(rfreqr * time);
        d2fdt2 = -(rfreqr) * (rfreqr) * std::sin(rfreqr * time);
    } else if (idef == 999) {
        return;
    }

    if (std::abs((float)omegx) > 0.0f) {

        // rotate about an axis parallel to the x-axis
        // calculate rotated y and z surface points
        // delt(1)=0 (unaltered) delt(2)=dy delt(3)=dz

        double theold;
        if (time != 0.) {
            theold = thetax;
        } else {
            theold = 0.0;
        }

        double theta    = omegx * ft;
        double dthedt   = omegx * dfdt;
        double d2thedt2 = omegx * d2fdt2;
        double dtheta   = theta - theold;
        double ca = std::cos(dtheta);
        double sa = std::sin(dtheta);

        if (isurf == 1) {
            for (int j = jci; j <= jcf; j++) {
                for (int k = kci; k <= kcf; k++) {
                    double tempy           = y(j, k, ii);
                    double tempz           = z(j, k, ii);
                    delti(j, k, 2, ll) = wki(j, k, ll) * ((tempy - yorg) * ca
                                       - (tempz - zorg) * sa + yorg - tempy)
                                       + delti(j, k, 2, ll);
                    delti(j, k, 3, ll) = wki(j, k, ll) * ((tempy - yorg) * sa
                                       + (tempz - zorg) * ca + zorg - tempz)
                                       + delti(j, k, 3, ll);
                    wki(j, k, ll)      = 0.;
                }
            }
        } else if (isurf == 2) {
            for (int k = kci; k <= kcf; k++) {
                for (int i = ici; i <= icf; i++) {
                    double tempy           = y(jj, k, i);
                    double tempz           = z(jj, k, i);
                    deltj(k, i, 2, ll) = wkj(k, i, ll) * ((tempy - yorg) * ca
                                       - (tempz - zorg) * sa + yorg - tempy)
                                       + deltj(k, i, 2, ll);
                    deltj(k, i, 3, ll) = wkj(k, i, ll) * ((tempy - yorg) * sa
                                       + (tempz - zorg) * ca + zorg - tempz)
                                       + deltj(k, i, 3, ll);
                    wkj(k, i, ll)      = 0.;
                }
            }
        } else if (isurf == 3) {
            for (int j = jci; j <= jcf; j++) {
                for (int i = ici; i <= icf; i++) {
                    double tempy           = y(j, kk, i);
                    double tempz           = z(j, kk, i);
                    deltk(j, i, 2, ll) = wkk(j, i, ll) * ((tempy - yorg) * ca
                                       - (tempz - zorg) * sa + yorg - tempy)
                                       + deltk(j, i, 2, ll);
                    deltk(j, i, 3, ll) = wkk(j, i, ll) * ((tempy - yorg) * sa
                                       + (tempz - zorg) * ca + zorg - tempz)
                                       + deltk(j, i, 3, ll);
                    wkk(j, i, ll)      = 0.;
                }
            }
        }

        thetax = theta;

    } else if (std::abs((float)omegy) > 0.0f) {

        // rotate about an axis parallel to the y-axis
        // calculate rotated x and z surface points
        // delt(1)=dx delt(2)=0 (unaltered) delt(3)=dz

        double theold;
        if (time != 0.) {
            theold = thetay;
        } else {
            theold = 0.0;
        }

        double theta    = omegy * ft;
        double dthedt   = omegy * dfdt;
        double d2thedt2 = omegy * d2fdt2;
        double dtheta   = theta - theold;
        double ca = std::cos(dtheta);
        double sa = std::sin(dtheta);

        if (isurf == 1) {
            for (int j = jci; j <= jcf; j++) {
                for (int k = kci; k <= kcf; k++) {
                    double tempx           = x(j, k, ii);
                    double tempz           = z(j, k, ii);
                    delti(j, k, 1, ll) = wki(j, k, ll) * ((tempx - xorg) * ca
                                       + (tempz - zorg) * sa + xorg - tempx)
                                       + delti(j, k, 1, ll);
                    delti(j, k, 3, ll) = wki(j, k, ll) * (-(tempx - xorg) * sa
                                       + (tempz - zorg) * ca + zorg - tempz)
                                       + delti(j, k, 3, ll);
                    wki(j, k, ll)      = 0.;
                }
            }
        } else if (isurf == 2) {
            for (int k = kci; k <= kcf; k++) {
                for (int i = ici; i <= icf; i++) {
                    double tempx           = x(jj, k, i);
                    double tempz           = z(jj, k, i);
                    deltj(k, i, 1, ll) = wkj(k, i, ll) * ((tempx - xorg) * ca
                                       + (tempz - zorg) * sa + xorg - tempx)
                                       + deltj(k, i, 1, ll);
                    deltj(k, i, 3, ll) = wkj(k, i, ll) * (-(tempx - xorg) * sa
                                       + (tempz - zorg) * ca + zorg - tempz)
                                       + deltj(k, i, 3, ll);
                    wkj(k, i, ll)      = 0.;
                }
            }
        } else if (isurf == 3) {
            for (int j = jci; j <= jcf; j++) {
                for (int i = ici; i <= icf; i++) {
                    double tempx           = x(j, kk, i);
                    double tempz           = z(j, kk, i);
                    deltk(j, i, 1, ll) = wkk(j, i, ll) * ((tempx - xorg) * ca
                                       + (tempz - zorg) * sa + xorg - tempx)
                                       + deltk(j, i, 1, ll);
                    deltk(j, i, 3, ll) = wkk(j, i, ll) * (-(tempx - xorg) * sa
                                       + (tempz - zorg) * ca + zorg - tempz)
                                       + deltk(j, i, 3, ll);
                    wkk(j, i, ll)      = 0.;
                }
            }
        }

        thetay = theta;

    } else if (std::abs((float)omegz) > 0.0f) {

        // rotate about an axis parallel to the z-axis
        // calculate rotated x and y surface points
        // delt(1)=dx delt(2)=dy delt(3)=0 (unaltered)

        double theold;
        if (time != 0.) {
            theold = thetaz;
        } else {
            theold = 0.0;
        }

        double theta    = omegz * ft;
        double dthedt   = omegz * dfdt;
        double d2thedt2 = omegz * d2fdt2;
        double dtheta   = theta - theold;
        double ca = std::cos(dtheta);
        double sa = std::sin(dtheta);

        if (isurf == 1) {
            for (int j = jci; j <= jcf; j++) {
                for (int k = kci; k <= kcf; k++) {
                    double tempx           = x(j, k, ii);
                    double tempy           = y(j, k, ii);
                    delti(j, k, 1, ll) = wki(j, k, ll) * ((tempx - xorg) * ca
                                       - (tempy - yorg) * sa + xorg - tempx)
                                       + delti(j, k, 1, ll);
                    delti(j, k, 2, ll) = wki(j, k, ll) * ((tempx - xorg) * sa
                                       + (tempy - yorg) * ca + yorg - tempy)
                                       + delti(j, k, 2, ll);
                    wki(j, k, ll)      = 0.;
                }
            }
        } else if (isurf == 2) {
            for (int k = kci; k <= kcf; k++) {
                for (int i = ici; i <= icf; i++) {
                    double tempx           = x(jj, k, i);
                    double tempy           = y(jj, k, i);
                    deltj(k, i, 1, ll) = wkj(k, i, ll) * ((tempx - xorg) * ca
                                       - (tempy - yorg) * sa + xorg - tempx)
                                       + deltj(k, i, 1, ll);
                    deltj(k, i, 2, ll) = wkj(k, i, ll) * ((tempx - xorg) * sa
                                       + (tempy - yorg) * ca + yorg - tempy)
                                       + deltj(k, i, 2, ll);
                    wkj(k, i, ll)      = 0.;
                }
            }
        } else if (isurf == 3) {
            for (int j = jci; j <= jcf; j++) {
                for (int i = ici; i <= icf; i++) {
                    double tempx           = x(j, kk, i);
                    double tempy           = y(j, kk, i);
                    deltk(j, i, 1, ll) = wkk(j, i, ll) * ((tempx - xorg) * ca
                                       - (tempy - yorg) * sa + xorg - tempx)
                                       + deltk(j, i, 1, ll);
                    deltk(j, i, 2, ll) = wkk(j, i, ll) * ((tempx - xorg) * sa
                                       + (tempy - yorg) * ca + yorg - tempy)
                                       + deltk(j, i, 2, ll);
                    wkk(j, i, ll)      = 0.;
                }
            }
        }

        thetaz = theta;

    } else {

        if (isklton > 0) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "%-40s%-21s",
                "WARNING: this block has zero rotational",
                " surface displacement");
        }

    }

    return;
}

} // namespace rotsurf_ns
