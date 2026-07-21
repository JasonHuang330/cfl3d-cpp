// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rotateqb.h"
#include "rotate.h"
#include "rotateq.h"
#include <cmath>

namespace rotateqb_ns {

void rotate(int& jdim, int& kdim, int& idim, FortranArray2DRef<double> t, FortranArray3DRef<double> tti, FortranArray3DRef<double> ttj, FortranArray3DRef<double> ttk, FortranArray2DRef<double> t1, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, int& nbl, int& irot, double& rfreqr, double& omegx, double& omegy, double& omegz, double& xorg, double& yorg, double& zorg, double& thetax, double& thetay, double& thetaz, double& thxold, double& thyold, double& thzold, int& iupdat, double& time2, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim) {
    rotate_ns::rotate(jdim, kdim, idim, t, tti, ttj, ttk, t1, x, y, z, nbl, irot, rfreqr, omegx, omegy, omegz, xorg, yorg, zorg, thetax, thetay, thetaz, thxold, thyold, thzold, iupdat, time2, nou, bou, nbuf, ibufdim);
}

void rotateq(int& jd, int& kd, int& id, FortranArray4DRef<double> q, FortranArray4DRef<double> qrot, int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend, double& dthtx, double& dthty, double& dthtz) {
    rotateq_ns::rotateq(jd, kd, id, q, qrot, ista, iend, jsta, jend, ksta, kend, dthtx, dthty, dthtz);
}

void rotateqb(int& nbl, double& dthtx, double& dthty, double& dthtz, int& maxbl, int& iitot, FortranArray1DRef<int> ibcg, FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg, FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg, FortranArray3DRef<double> qb) {
//
//***********************************************************************
//     Purpose:  Rotate solution in qb array through angle dthtx/y/z for
//     chimera scheme with rotating grids
//***********************************************************************
//
    int iset = 1;
//
//     lsta = lbg(nbl)
//     lend = lsta-1
    int lsta = lig(nbl);
    int lend = lsta + iipntsg(nbl) - 1;
//
    double ca, sa;

    if (std::abs((float)dthtx) > 0.) {
//
//        rotate qb about an axis parallel to the x-axis
//
        ca = std::cos(dthtx);
        sa = std::sin(dthtx);
//
//        if (ibpntsg(nbl,1).gt.0) then
//           lend = lsta+ibpntsg(nbl,1)-1
        for (int l = lsta; l <= lend; l++) {
//           qb3 = qb(ibcg(l),3,iset)
//           qb(ibcg(l),3,iset) = qb3*ca - qb(ibcg(l),4,iset)*sa
//           qb(ibcg(l),4,iset) = qb3*sa + qb(ibcg(l),4,iset)*ca
            double qb3 = qb(l, 3, iset);
            qb(l, 3, iset) = qb3 * ca - qb(l, 4, iset) * sa;
            qb(l, 4, iset) = qb3 * sa + qb(l, 4, iset) * ca;
        }
//        end if
//
    }
//
    if (std::abs((float)dthty) > 0.) {
//
//        rotate qb about an axis parallel to the y-axis
//
        ca = std::cos(dthty);
        sa = std::sin(dthty);
//
//        if (ibpntsg(nbl,1).gt.0) then
//           lend = lsta+ibpntsg(nbl,1)-1
        for (int l = lsta; l <= lend; l++) {
//           qb2 = qb(ibcg(l),2,iset)
//           qb(ibcg(l),2,iset) =  qb2*ca + qb(ibcg(l),4,iset)*sa
//           qb(ibcg(l),4,iset) = -qb2*sa + qb(ibcg(l),4,iset)*ca
            double qb2 = qb(l, 2, iset);
            qb(l, 2, iset) =  qb2 * ca + qb(l, 4, iset) * sa;
            qb(l, 4, iset) = -qb2 * sa + qb(l, 4, iset) * ca;
        }
//        end if
//
    }
//
    if (std::abs((float)dthtz) > 0.) {
//
//        rotate qb about an axis parallel to the z-axis
//
        ca = std::cos(dthtz);
        sa = std::sin(dthtz);
//
//        if (ibpntsg(nbl,1).gt.0) then
//           lend = lsta+ibpntsg(nbl,1)-1
        for (int l = lsta; l <= lend; l++) {
//           qb2 = qb(ibcg(l),2,iset)
//           qb(ibcg(l),2,iset) = qb2*ca - qb(ibcg(l),3,iset)*sa
//           qb(ibcg(l),3,iset) = qb2*sa + qb(ibcg(l),3,iset)*ca
            double qb2 = qb(l, 2, iset);
            qb(l, 2, iset) = qb2 * ca - qb(l, 3, iset) * sa;
            qb(l, 3, iset) = qb2 * sa + qb(l, 3, iset) * ca;
        }
//        end if
//
    }
//
    return;
}

} // namespace rotateqb_ns
