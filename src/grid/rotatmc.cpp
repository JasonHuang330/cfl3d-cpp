// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rotatmc.h"
#include <cmath>

namespace rotatmc_ns {

void rotatmc(int& irot, double& rfreqr, double& omegx, double& omegy, double& omegz, double& xorg, double& yorg, double& zorg, double& thetax, double& thetay, double& thetaz, double& xmc, double& ymc, double& zmc, int& iupdat, double& time2)
{
//
//     ft modulates the rotation
//
      double ft, expt, theold, theta, dtheta, ca, sa;
      double xml, yml, zml;

      if (irot == 0) {
         return;
      } else if (irot == 1) {
         ft     = time2;
      } else if (irot == 2) {
         ft     = std::sin(rfreqr*time2);
      } else if (irot == 3) {
         expt   = std::exp(-rfreqr*time2);
         ft     = 1.-expt;
      }

      if (std::fabs((double)(omegx)) > 0.0) {

//***********************************************************************
//        rotate about an axis parallel to the x-axis
//***********************************************************************
//
         if (time2 != 0.) {
            theold = thetax;
         } else {
            theold = 0.e0;
         }
//
         theta    = omegx*ft;
         dtheta   = theta - theold;
         ca = std::cos(dtheta);
         sa = std::sin(dtheta);
//
         if (iupdat > 0) {
            thetax = theta;
            yml = ymc;
            zml = zmc;
            xmc = xmc;
            ymc = (yml-yorg)*ca-(zml-zorg)*sa+yorg;
            zmc = (yml-yorg)*sa+(zml-zorg)*ca+zorg;
         }
//
      } else if (std::fabs((double)(omegy)) > 0.0) {
//
//***********************************************************************
//        rotate about an axis parallel to the y-axis
//***********************************************************************
//
         if (time2 != 0.) {
            theold = thetay;
         } else {
            theold = 0.e0;
         }
//
         theta    = omegy*ft;
         dtheta   = theta - theold;
         ca = std::cos(dtheta);
         sa = std::sin(dtheta);
//
         if (iupdat > 0) {
            thetay = theta;
            xml = xmc;
            zml = zmc;
            xmc =  (xml-xorg)*ca+(zml-zorg)*sa+xorg;
            ymc =   ymc;
            zmc = -(xml-xorg)*sa+(zml-zorg)*ca+zorg;
         }
//
      } else if (std::fabs((double)(omegz)) > 0.0) {
//
//***********************************************************************
//        rotate about an axis parallel to the z-axis
//***********************************************************************
//
         if (time2 != 0.) {
            theold = thetaz;
         } else {
            theold = 0.e0;
         }
//
         theta    = omegz*ft;
         dtheta   = theta - theold;
         ca = std::cos(dtheta);
         sa = std::sin(dtheta);
//
         if (iupdat > 0) {
            thetaz = theta;
            xml = xmc;
            yml = ymc;
            xmc = (xml-xorg)*ca-(yml-yorg)*sa+xorg;
            ymc = (xml-xorg)*sa+(yml-yorg)*ca+yorg;
            zmc = zmc;
         }
//
      }
//
      return;
}

} // namespace rotatmc_ns
