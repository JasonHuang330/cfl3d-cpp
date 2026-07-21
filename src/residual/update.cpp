// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "update.h"
#include "af3f.h"
#include "mms.h"
#include "termn8.h"
#include "ccomplex.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <algorithm>

namespace update_ns {

void update(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray2DRef<double> qj0, FortranArray2DRef<double> qk0, FortranArray2DRef<double> qi0, FortranArray2DRef<double> sj, FortranArray2DRef<double> sk, FortranArray2DRef<double> si, FortranArray1DRef<double> vol, FortranArray1DRef<double> dtj, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> blank, FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z, FortranArray4DRef<double> res, FortranArray1DRef<double> wk0, FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui, FortranArray1DRef<double> wk, int& nwork, int& nbl, int& iover, FortranArray4DRef<double> vk0, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, FortranArray3DRef<double> volk0, FortranArray4DRef<double> tursav, FortranArray4DRef<double> tk0, FortranArray3DRef<double> cmuv, FortranArray1DRef<int> iadvance, int& nummem, FortranArray4DRef<double> ux)
{
    // COMMON block references
    int& ichk         = cmn_chk.ichk;
    int& isklton      = cmn_sklton.isklton;
    int& i2d          = cmn_twod.i2d;
    int& iexact_trunc = cmn_mms.iexact_trunc;
    int& iexact_disc  = cmn_mms.iexact_disc;
    int& iexact_ring  = cmn_mms.iexact_ring;
    int& iforcev0     = cmn_axisym.iforcev0;
    int& iupdatemean  = cmn_iupdate.iupdatemean;

    if (iupdatemean == 0) return;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    double term = 1.0;
    if (i2d == 1 || iforcev0 == 1) term = 0.0;

    double alpq = -0.2;
    double phiq = 1.0 / 0.5;
    double betq = 1.0 + alpq * phiq;

    if (isklton > 0) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::memset(bou(nou(1), 1), ' ', 120);
        bou(nou(1), 1)[0] = '\0';
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
            " updating with 3-D 3-factor AF scheme for block%6d", nbl);
    }

    // assemble and solve matrix equations
    double resd = 0.0;
    int nwk0_val = 1;
    {
        // Create flat/reshaped views for af3f which uses different array dimensions
        int idim1 = idim - 1;
        FortranArray2DRef<double> q_2d(&q(1,1,1,1), jdim*kdim*idim, 5);
        FortranArray3DRef<double> qj0_3d(&qj0(1,1), kdim*idim1, 5, 1);
        FortranArray3DRef<double> qk0_3d(&qk0(1,1), jdim*idim1, 5, 1);
        FortranArray3DRef<double> qi0_3d(&qi0(1,1), jdim*kdim, 5, 1);
        FortranArray2DRef<double> res_2d(&res(1,1,1,1), jdim*kdim*idim1, 5);
        FortranArray3DRef<double> x_3d(&x(1), jdim, kdim, idim);
        FortranArray3DRef<double> y_3d(&y(1), jdim, kdim, idim);
        FortranArray3DRef<double> z_3d(&z(1), jdim, kdim, idim);
        af3f_ns::af3f(nbl, jdim, kdim, idim, q_2d, vol, qj0_3d, qk0_3d, qi0_3d, dtj,
                      sj, sk, si, res_2d, vist3d, x_3d, y_3d, z_3d, blank, vmuk, resd,
                      wk, nwork, wk0, nwk0_val, iover,
                      vk0, bcj, bck, bci, nou, bou, nbuf, ibufdim, myid,
                      mblk2nd, maxbl, volk0, tursav, tk0, cmuv, iadvance,
                      nummem, ux);
    }

    // update 3-d solution
    int nt   = jdim * kdim;
    int nplq = std::min(idim1, 999000 / nt);
    int npl  = nplq;

    for (int i = 1; i <= idim1; i += nplq) {
        if (i + npl - 1 > idim1) npl = idim1 - i + 1;
        int n = nt * npl - jdim - 1;

        for (int izz = 1; izz <= n; izz++) {
            // ensure positivity of density
            double t1 = res(izz, 1, i, 1) / q(izz, 1, i, 1);
            double t2 = res(izz, 1, i, 1) / (betq + ccomplex_ns::ccabs(t1) * phiq);
            bool cond1 = ((float)t1 < (float)alpq);
            res(izz, 1, i, 1) = ccomplex_ns::ccvmgt(t2, res(izz, 1, i, 1), cond1);

            // ensure positivity of pressure
            t1 = res(izz, 1, i, 5) / q(izz, 1, i, 5);
            t2 = res(izz, 1, i, 5) / (betq + ccomplex_ns::ccabs(t1) * phiq);
            bool cond2 = ((float)t1 < (float)alpq);
            res(izz, 1, i, 5) = ccomplex_ns::ccvmgt(t2, res(izz, 1, i, 5), cond2);

            q(izz, 1, i, 1) =  q(izz, 1, i, 1) + res(izz, 1, i, 1);
            q(izz, 1, i, 2) =  q(izz, 1, i, 2) + res(izz, 1, i, 2);
            q(izz, 1, i, 3) = (q(izz, 1, i, 3) + res(izz, 1, i, 3)) * term;
            q(izz, 1, i, 4) =  q(izz, 1, i, 4) + res(izz, 1, i, 4);
            q(izz, 1, i, 5) =  q(izz, 1, i, 5) + res(izz, 1, i, 5);
        }
    }

    if (iexact_ring == 1) {
        // overwrite ring of exact MMS values in outer 2 rows of grid
        FortranArray3DRef<double> x_3d(&x(1), jdim, kdim, idim);
        FortranArray3DRef<double> y_3d(&y(1), jdim, kdim, idim);
        FortranArray3DRef<double> z_3d(&z(1), jdim, kdim, idim);
        mms_ns::exact_flow_q_ring(jdim, kdim, idim, x_3d, y_3d, z_3d, q,
                                  iexact_trunc, iexact_disc);
    }

    if (ichk == 1) {
        double epsz = 1.0e-05;
        double epss = 1.0e+03;
        for (int i = 1; i <= idim1; i++) {
            for (int k = 1; k <= kdim1; k++) {
                for (int j = 1; j <= jdim1; j++) {
                    if ((float)q(j, k, i, 5) < (float)epsz ||
                        (float)q(j, k, i, 1) < (float)epsz ||
                        (float)q(j, k, i, 5) > (float)epss ||
                        (float)q(j, k, i, 1) > (float)epss) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " stopping in update on block %d", nbl);
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " *neg. (or large) d/p*(j,k,i,q)=%5d%5d%5d%12.5e%12.5e%12.5e%12.5e%12.5e",
                            j, k, i,
                            (float)q(j, k, i, 1),
                            (float)q(j, k, i, 2),
                            (float)q(j, k, i, 3),
                            (float)q(j, k, i, 4),
                            (float)q(j, k, i, 5));
                        int ierrflg = -1;
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
    }
}

} // namespace update_ns
