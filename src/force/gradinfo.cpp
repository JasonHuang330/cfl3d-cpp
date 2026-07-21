// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "gradinfo.h"
#include "runtime/fortran_io.h"

namespace gradinfo_ns {

void gradinfo(double& delh, int& iunit)
{
    // COMMON /complx/
    float& xmach_img    = cmn_complx.xmach_img;
    float& alpha_img    = cmn_complx.alpha_img;
    float& beta_img     = cmn_complx.beta_img;
    float& reue_img     = cmn_complx.reue_img;
    float& tinf_img     = cmn_complx.tinf_img;
    float& geom_img     = cmn_complx.geom_img;
    float& surf_img     = cmn_complx.surf_img;
    float& xrotrate_img = cmn_complx.xrotrate_img;
    float& yrotrate_img = cmn_complx.yrotrate_img;
    float& zrotrate_img = cmn_complx.zrotrate_img;

    // COMMON /igrdtyp/
    int32_t& ialph = cmn_igrdtyp.ialph;

    if (xmach_img > 0.f) {
        delh = (double)xmach_img;
        if (iunit > 0) {
            fortran_write_unit(iunit,
                "                    derivatives with respect to Mach number\n");
        }
    } else if (alpha_img > 0.f) {
        delh = (double)alpha_img;
        if (iunit > 0) {
            fortran_write_unit(iunit,
                "                       derivatives with respect to alpha\n");
        }
    } else if (beta_img > 0.f) {
        delh = (double)beta_img;
        if (iunit > 0) {
            fortran_write_unit(iunit,
                "                       derivatives with respect to beta\n");
        }
    } else if (reue_img > 0.f) {
        delh = (double)reue_img;
        if (iunit > 0) {
            fortran_write_unit(iunit,
                "                       derivatives with respect to reue\n");
        }
    } else if (tinf_img > 0.f) {
        delh = (double)tinf_img;
        if (iunit > 0) {
            fortran_write_unit(iunit,
                "                       derivatives with respect to Tinf\n");
        }
    } else if (geom_img > 0.f) {
        delh = (double)geom_img;
        if (iunit > 0) {
            fortran_write_unit(iunit,
                "                      derivatives with respect to geometry\n");
        }
    } else if (surf_img > 0.f) {
        delh = (double)surf_img;
        if (iunit > 0) {
            fortran_write_unit(iunit,
                "                      derivatives with respect to (surface) geometry\n");
        }
    } else if (xrotrate_img > 0.f) {
        delh = (double)xrotrate_img;
        if (iunit > 0) {
            fortran_write_unit(iunit,
                "                      derivatives with respect to omega_x (roll rate)\n");
        }
    } else if (yrotrate_img > 0.f) {
        delh = (double)yrotrate_img;
        if (iunit > 0) {
            if (ialph == 0) {
                fortran_write_unit(iunit,
                    "                      derivatives with respect to omega_y (pitch rate)\n");
            } else {
                fortran_write_unit(iunit,
                    "                      derivatives with respect to omega_y (yaw rate)\n");
            }
        }
    } else if (zrotrate_img > 0.f) {
        delh = (double)zrotrate_img;
        if (iunit > 0) {
            if (ialph == 0) {
                fortran_write_unit(iunit,
                    "                      derivatives with respect to omega_z (yaw rate)\n");
            } else {
                fortran_write_unit(iunit,
                    "                      derivatives with respect to omega_z (pitch rate)\n");
            }
        }
    }

    return;
}

} // namespace gradinfo_ns
