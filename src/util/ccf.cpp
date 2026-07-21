// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ccf.h"
#include <cmath>

namespace ccf_ns {

void ccf(double& x, double& y, double& ca, double& sa, double& cl, double& xm, double& uf, double& wf, double& cf, double& pi)
{
//
//   far field boundary condition (point vortex)
//
    double xe = ca*x + sa*y;
    double ye = -sa*x + ca*y;
    double beta = std::sqrt(1.e0 - xm*xm);
    double re = xe*xe + (ye*beta)*(ye*beta);
    double fact = cl*beta*xm / (re*4.e0*pi);
    double ue = xm + fact*ye;
    double ve = -fact*xe;
    uf = ca*ue - sa*ve;
    wf = sa*ue + ca*ve;
    cf = std::sqrt(1.e0 + 0.2e0*(xm*xm - uf*uf - wf*wf));
    return;
}

} // namespace ccf_ns
