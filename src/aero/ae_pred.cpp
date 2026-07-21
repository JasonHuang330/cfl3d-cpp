// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ae_pred.h"
#include "moddefl.h"

namespace ae_pred_ns {

void ae_pred(FortranArray2DRef<double> aesrfdat, FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi, FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcnm, FortranArray2DRef<double> xs, FortranArray2DRef<double> xxn, FortranArray2DRef<double> x0, FortranArray3DRef<double> perturb, double& cmyt, double& cnwt, double& xorig, double& yorig, double& zorig, int& nmds, int& maxaes, int& irbtrim, int& maxbl, int& myid)
{
    // common /elastic/ ndefrm, naesrf
    int& ndefrm = cmn_elastic.ndefrm;
    int& naesrf = cmn_elastic.naesrf;

    // modal displacement and velocity prediction via aeroelastic
    // equations of motion
    for (int iaes = 1; iaes <= naesrf; iaes++) {
        int nmodes = (int)aesrfdat(5, iaes);
        int iskyhk = (int)aesrfdat(1, iaes);
        for (int n = 1; n <= 2 * nmodes; n++) {
            // don't update if the modal time variation is specified
            int moddfl = (int)perturb((n + 1) / 2, iaes, 1);
            if (moddfl == 0) {
                xs(n, iaes) = 0.;
                for (int j = 1; j <= 2 * nmodes; j++) {
                    xs(n, iaes) = xs(n, iaes) + stm(n, j, iaes) * xxn(j, iaes)
                                + .5 * stmi(n, j, iaes) * (3. * gforcn(j, iaes)
                                - gforcnm(j, iaes));
                }
            }
        }
    }

    // modal displacement and velocity as specified functions of time
    moddefl_ns::moddefl(xs, xxn, aesrfdat, perturb, cmyt, cnwt, xorig, yorig,
                        zorig, maxaes, nmds, irbtrim, maxbl, myid);

    return;
}

} // namespace ae_pred_ns
