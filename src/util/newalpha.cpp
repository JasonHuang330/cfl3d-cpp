// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "newalpha.h"
#include "runtime/fortran_io.h"
#include "ccomplex.h"
#include "init_mast.h"
#include <cmath>
#include <cstdio>

namespace newalpha_ns {

void newalpha(int& ncycmax, FortranArray1DRef<double> rms, FortranArray1DRef<double> clw, int& myid, int& myhost, int& mycomm)
{
    // COMMON block aliases
    int32_t& ialphit  = cmn_alphait.ialphit;
    float&   cltarg   = cmn_alphait.cltarg;
    float&   rlxalph  = cmn_alphait.rlxalph;
    float&   dalim    = cmn_alphait.dalim;
    float&   dalpha   = cmn_alphait.dalpha;
    int32_t& icycupdt = cmn_alphait.icycupdt;

    float&   alpha    = cmn_info.alpha;
    int32_t& ntt      = cmn_info.ntt;

    float&   radtodeg = cmn_conversion.radtodeg;

    // cltarg...target cl
    // rlxalph..relaxation parameter for updating alpha (default = 1)
    // dalim....limit on magnitude of alpha change (default = 0.2 deg)
    // dalpha...alpha change
    // ioalph...if > 0, output alpha history to cfl3d.alpha
    int ioalph = 1;

    // update alpha based on error in cl;
    if (myid == myhost) {
        double rlx = (double)rlxalph / (double)radtodeg;
        dalpha = (float)(-rlx * ((double)clw(ntt) - (double)cltarg));
        if (std::fabs((double)dalpha) > (double)dalim / (double)radtodeg) {
            double a1 = 1.0;
            double b1 = (double)dalpha;
            double sgn = ccomplex_ns::ccsignrc(a1, b1);
            dalpha = (float)(sgn * ((double)dalim / (double)radtodeg));
        }
        alpha = alpha + dalpha;

        // output to cfl3d.alpha
        if (ioalph > 0) {
            FILE* f27 = fortran_get_unit(27);
            fprintf(f27, "   %6d%16.8E%16.8E%16.8E%16.8E\n",
                    (int)ntt,
                    (double)std::log10((double)rms(ntt)),
                    (double)clw(ntt),
                    (double)alpha * (double)radtodeg,
                    (double)dalpha * (double)radtodeg);
        }
    }

    // reset freestream for new alpha
    init_mast_ns::init_mast();

    return;
}

} // namespace newalpha_ns
