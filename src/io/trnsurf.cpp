// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "trnsurf.h"
#include <cmath>

namespace trnsurf_ns {

void trnsurf(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti, int& nbl, int& idef, double& xorg, double& yorg, double& zorg, double& utran, double& vtran, double& wtran, double& rfreqt, int& ici, int& icf, int& jci, int& jcf, int& kci, int& kcf, double& time, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray3DRef<double> wkj, FortranArray3DRef<double> wkk, FortranArray3DRef<double> wki, double& xorg0, double& yorg0, double& zorg0)
{
    // common /sklton/ isklton
    int& isklton = cmn_sklton.isklton;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    // determine which block face is to be translated
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

    // ft modulates the translation
    // dfdt is the time derivative of ft
    // d2fdt2 is the second time derivative of ft
    double ft     = 0.0;
    double dfdt   = 0.0;
    double d2fdt2 = 0.0;

    if (idef == 0) {
        return;
    } else if (idef == 1) {
        ft     = std::sin(rfreqt * time);
        dfdt   = rfreqt * std::cos(rfreqt * time);
        d2fdt2 = -(rfreqt) * (rfreqt) * std::sin(rfreqt * time);
    } else if (idef == 999) {
        return;
    }

    double xold = xorg;
    double yold = yorg;
    double zold = zorg;

    double xnew = utran * ft + xorg0;
    double ynew = vtran * ft + yorg0;
    double znew = wtran * ft + zorg0;

    double dx = xnew - xold;
    double dy = ynew - yold;
    double dz = znew - zold;

    if (isurf == 1) {
        for (int j = jci; j <= jcf; j++) {
            for (int k = kci; k <= kcf; k++) {
                delti(j, k, 1, ll) = wki(j, k, ll) * dx + delti(j, k, 1, ll);
                delti(j, k, 2, ll) = wki(j, k, ll) * dy + delti(j, k, 2, ll);
                delti(j, k, 3, ll) = wki(j, k, ll) * dz + delti(j, k, 3, ll);
                wki(j, k, ll)      = 0.;
            }
        }
    } else if (isurf == 2) {
        for (int k = kci; k <= kcf; k++) {
            for (int i = ici; i <= icf; i++) {
                deltj(k, i, 1, ll) = wkj(k, i, ll) * dx + deltj(k, i, 1, ll);
                deltj(k, i, 2, ll) = wkj(k, i, ll) * dy + deltj(k, i, 2, ll);
                deltj(k, i, 3, ll) = wkj(k, i, ll) * dz + deltj(k, i, 3, ll);
                wkj(k, i, ll)      = 0.;
            }
        }
    } else if (isurf == 3) {
        for (int j = jci; j <= jcf; j++) {
            for (int i = ici; i <= icf; i++) {
                deltk(j, i, 1, ll) = wkk(j, i, ll) * dx + deltk(j, i, 1, ll);
                deltk(j, i, 2, ll) = wkk(j, i, ll) * dy + deltk(j, i, 2, ll);
                deltk(j, i, 3, ll) = wkk(j, i, ll) * dz + deltk(j, i, 3, ll);
                wkk(j, i, ll)      = 0.;
            }
        }
    }

    xorg = xorg + dx;
    yorg = yorg + dy;
    zorg = zorg + dz;

    return;
}

} // namespace trnsurf_ns
