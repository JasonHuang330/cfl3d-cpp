// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "initvist.h"
#include <cstdio>
#include "mms.h"
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <cmath>
#include <cstdint>

namespace initvist_ns {

// Inline helpers for ccabs/ccmin/ccmax on double rvalues
static inline double ccabs(double x) { return std::fabs(x); }
static inline double ccmin(double a, double b) { return (a < b) ? a : b; }
static inline double ccmax(double a, double b) { return (a > b) ? a : b; }

void init(int& nbl, int& jdim, int& kdim, int& idim,
          FortranArray4DRef<double> q,
          FortranArray4DRef<double> qj0,
          FortranArray4DRef<double> qk0,
          FortranArray4DRef<double> qi0,
          FortranArray4DRef<double> tj0,
          FortranArray4DRef<double> tk0,
          FortranArray4DRef<double> ti0,
          FortranArray3DRef<double> vol,
          FortranArray3DRef<double> volj0,
          FortranArray3DRef<double> volk0,
          FortranArray3DRef<double> voli0,
          int& nummem,
          FortranArray3DRef<double> x,
          FortranArray3DRef<double> y,
          FortranArray3DRef<double> z,
          FortranArray1DRef<int> nou,
          FortranArray2DRef<char[120]> bou,
          int& nbuf,
          int& ibufdim,
          int& iflagprnt)
{
    (void)nbl; (void)jdim; (void)kdim; (void)idim;
    (void)q; (void)qj0; (void)qk0; (void)qi0;
    (void)tj0; (void)tk0; (void)ti0;
    (void)vol; (void)volj0; (void)volk0; (void)voli0;
    (void)nummem;
    (void)x; (void)y; (void)z;
    (void)nou; (void)bou;
    (void)nbuf; (void)ibufdim; (void)iflagprnt;
}

void initvist(int& nbl, int& jdim, int& kdim, int& idim,
              FortranArray3DRef<double> vist3d,
              FortranArray4DRef<double> tursav,
              FortranArray3DRef<double> smin,
              FortranArray3DRef<double> cmuv,
              int& nummem,
              FortranArray3DRef<double> x,
              FortranArray3DRef<double> y,
              FortranArray3DRef<double> z)
{
    (void)nbl;

    // COMMON block aliases
    float  scal_ic      = cmn_initfac.scal_ic;
    float  rho0         = cmn_ivals.rho0;
    float* tur10        = cmn_ivals.tur10;   // 0-based: tur10[0]..tur10[6]
    float  reue         = cmn_reyue.reue;
    int*   ivisc        = cmn_reyue.ivisc;   // 0-based: ivisc[0]..ivisc[2]
    int    iexact_trunc = cmn_mms.iexact_trunc;
    int    iexact_disc  = cmn_mms.iexact_disc;
    int    issglrrw2012 = cmn_reystressmodel.issglrrw2012;

    // Local variables
    int    iscal;
    double refac = 0.0;
    int    i, j, k;

    // Scaling factor for distance-based I.C.s
    if (scal_ic == 0.f) {
        // Use freestream ICs everywhere
        iscal = 0;
    } else {
        // Scale turb-like IC in BL; lower scal_ic yields thinner layer
        refac = (double)reue / (double)scal_ic;
        iscal = 1;
    }

    // ivisc == 8,9,13,14: initialize cmuv
    if (ivisc[0] == 8  || ivisc[1] == 8  || ivisc[2] == 8  ||
        ivisc[0] == 9  || ivisc[1] == 9  || ivisc[2] == 9  ||
        ivisc[0] == 13 || ivisc[1] == 13 || ivisc[2] == 13 ||
        ivisc[0] == 14 || ivisc[1] == 14 || ivisc[2] == 14) {
        for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++)
                    cmuv(j,k,i) = -.09;
    }

    // For field eqn turbulence models:
    if (ivisc[0] >= 4 || ivisc[1] >= 4 || ivisc[2] >= 4) {
        if (ivisc[0] == 4 || ivisc[1] == 4 || ivisc[2] == 4) {
            for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                    for (j = 1; j <= jdim-1; j++) {
                        tursav(j,k,i,1) = (double)tur10[0];
                        vist3d(j,k,i)   = 0.09 * (double)tur10[0];
                    }
        } else if (ivisc[0] == 5 || ivisc[1] == 5 || ivisc[2] == 5) {
            for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                    for (j = 1; j <= jdim-1; j++) {
                        tursav(j,k,i,1) = (double)tur10[0];
                        vist3d(j,k,i)   = (double)tur10[0] *
                            (std::pow((double)tur10[0], 3.0) /
                             (std::pow((double)tur10[0], 3.0) + 357.911));
                    }
        } else if (ivisc[0] == 16 || ivisc[1] == 16 || ivisc[2] == 16) {
            for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                    for (j = 1; j <= jdim-1; j++) {
                        tursav(j,k,i,1) = (double)tur10[0];
                        tursav(j,k,i,2) = (double)tur10[1];
                        vist3d(j,k,i)   = 0.54772 * (double)tur10[0] /
                                          std::sqrt((double)tur10[1]);
                    }
        } else if (ivisc[0] == 11 || ivisc[1] == 11 || ivisc[2] == 11 ||
                   ivisc[0] == 10 || ivisc[1] == 10 || ivisc[2] == 10 ||
                   ivisc[0] ==  9 || ivisc[1] ==  9 || ivisc[2] ==  9 ||
                   ivisc[0] == 13 || ivisc[1] == 13 || ivisc[2] == 13 ||
                   ivisc[0] == 15 || ivisc[1] == 15 || ivisc[2] == 15) {
            for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                    for (j = 1; j <= jdim-1; j++) {
                        double zk1, zk2, zk3, zk4;
                        double ep1, ep2, ep3, ep4;
                        double exponent;

                        zk1 = (double)tur10[1];
                        exponent = -471.0 * (ccabs(smin(j,k,i)) * refac) + .47;
                        if (exponent < -40.0) {
                            zk2 = 0.0;
                        } else {
                            zk2 = std::pow(10.0, exponent);
                        }
                        exponent = -37.5 * (ccabs(smin(j,k,i)) * refac) - 3.7;
                        if (exponent < -40.0) {
                            zk3 = 0.0;
                        } else {
                            zk3 = std::pow(10.0, exponent);
                        }
                        zk4 = 6.7 * (ccabs(smin(j,k,i)) * refac);
                        tursav(j,k,i,2) = ccmin(zk2, zk3);
                        tursav(j,k,i,2) = ccmax(tursav(j,k,i,2), zk1);
                        tursav(j,k,i,2) = ccmin(tursav(j,k,i,2), zk4);
                        tursav(j,k,i,2) = tursav(j,k,i,2) * (double)iscal +
                                          (double)tur10[1] * (double)(1 - iscal);

                        ep1 = (double)tur10[0];
                        exponent = -555.0 * (ccabs(smin(j,k,i)) * refac) - 6.0;
                        if (exponent < -20.0) {
                            ep2 = 0.0;
                        } else {
                            ep2 = std::pow(10.0, exponent);
                        }
                        exponent = -280.0 * (ccabs(smin(j,k,i)) * refac) - 9.2;
                        if (exponent < -20.0) {
                            ep3 = 0.0;
                        } else {
                            ep3 = std::pow(10.0, exponent);
                        }
                        exponent = 13333.0 * (ccabs(smin(j,k,i)) * refac) - 9.8;
                        if (exponent > 20.0) {
                            ep4 = 1.e20;
                        } else {
                            ep4 = std::pow(10.0, exponent);
                        }
                        tursav(j,k,i,1) = ccmin(ep2, ep3);
                        tursav(j,k,i,1) = ccmax(tursav(j,k,i,1), ep1);
                        tursav(j,k,i,1) = ccmin(tursav(j,k,i,1), ep4);
                        tursav(j,k,i,1) = tursav(j,k,i,1) * (double)iscal +
                                          (double)tur10[0] * (double)(1 - iscal);
                        vist3d(j,k,i) = .09 * (double)rho0 *
                                        tursav(j,k,i,2) * tursav(j,k,i,2) /
                                        tursav(j,k,i,1);
                    }
        } else if (ivisc[0] == 25 || ivisc[1] == 25 || ivisc[2] == 25) {
            // do nothing
            for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                    for (j = 1; j <= jdim-1; j++) {
                        tursav(j,k,i,1) = 0.0;
                        vist3d(j,k,i)   = 0.0;
                    }
        } else if (ivisc[0] == 72 || ivisc[1] == 72 || ivisc[2] == 72) {
            if (issglrrw2012 == 6) {
                for (i = 1; i <= idim-1; i++)
                    for (k = 1; k <= kdim-1; k++)
                        for (j = 1; j <= jdim-1; j++) {
                            double zkinf;
                            tursav(j,k,i,1) = (double)tur10[0];
                            tursav(j,k,i,2) = (double)tur10[1];
                            tursav(j,k,i,3) = (double)tur10[2];
                            tursav(j,k,i,4) = (double)tur10[3];
                            tursav(j,k,i,5) = (double)tur10[4];
                            tursav(j,k,i,6) = (double)tur10[5];
                            tursav(j,k,i,7) = (double)tur10[6];
                            zkinf = -((double)tur10[0] + (double)tur10[1] +
                                      (double)tur10[2]) / 2.0;
                            vist3d(j,k,i) = (double)rho0 * zkinf *
                                            ((double)tur10[6] * (double)tur10[6]);
                        }
            } else {
                for (i = 1; i <= idim-1; i++)
                    for (k = 1; k <= kdim-1; k++)
                        for (j = 1; j <= jdim-1; j++) {
                            double zkinf;
                            tursav(j,k,i,1) = (double)tur10[0];
                            tursav(j,k,i,2) = (double)tur10[1];
                            tursav(j,k,i,3) = (double)tur10[2];
                            tursav(j,k,i,4) = (double)tur10[3];
                            tursav(j,k,i,5) = (double)tur10[4];
                            tursav(j,k,i,6) = (double)tur10[5];
                            tursav(j,k,i,7) = (double)tur10[6];
                            zkinf = -((double)tur10[0] + (double)tur10[1] +
                                      (double)tur10[2]) / 2.0;
                            vist3d(j,k,i) = (double)rho0 * zkinf / (double)tur10[6];
                        }
            }
        } else {
            double const1 = 45.8;
            double const2 = 1.68;
            double smax   = const2 / (2.0 * const1);
            double tmax   = -const1 * smax * smax + const2 * smax;
            double cmu    = 1.0;
            if (ivisc[0] ==  8 || ivisc[1] ==  8 || ivisc[2] ==  8 ||
                ivisc[0] == 12 || ivisc[1] == 12 || ivisc[2] == 12 ||
                ivisc[0] == 14 || ivisc[1] == 14 || ivisc[2] == 14)
                cmu = .09;
            double v3dset = cmu * (double)rho0 * (double)tur10[1] / (double)tur10[0];
            for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                    for (j = 1; j <= jdim-1; j++) {
                        double zk1, zk2;
                        double om1, om2, v3d;
                        zk1 = (double)tur10[1];
                        zk2 = -const1 * (ccabs(smin(j,k,i)) * refac) *
                                        (ccabs(smin(j,k,i)) * refac) +
                              const2 * ccabs(smin(j,k,i)) * refac;
                        tursav(j,k,i,2) = ccmax(zk1, zk2);
                        tursav(j,k,i,2) = tursav(j,k,i,2) * (double)iscal +
                                          (double)tur10[1] * (double)(1 - iscal);
                        v3d = tursav(j,k,i,2) * 100.0 / tmax;
                        v3d = ccmax(v3d, v3dset);
                        om1 = -12444.0 * (ccabs(smin(j,k,i)) * refac) + .54;
                        om2 = cmu * tursav(j,k,i,2) / v3d;
                        tursav(j,k,i,1) = ccmax(om1, om2);
                        tursav(j,k,i,1) = tursav(j,k,i,1) * (double)iscal +
                                          (double)tur10[0] * (double)(1 - iscal);
                        vist3d(j,k,i) = cmu * (double)rho0 *
                                        tursav(j,k,i,2) / tursav(j,k,i,1);
                    }
        }

        // Special for 3-eqn transition model ivisc=30
        if (ivisc[0] == 30 || ivisc[1] == 30 || ivisc[2] == 30) {
            for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                    for (j = 1; j <= jdim-1; j++)
                        tursav(j,k,i,3) = (double)tur10[2];
        }

        // Special for 4-eqn transition model ivisc=40
        if (ivisc[0] == 40 || ivisc[1] == 40 || ivisc[2] == 40) {
            for (i = 1; i <= idim-1; i++)
                for (k = 1; k <= kdim-1; k++)
                    for (j = 1; j <= jdim-1; j++) {
                        tursav(j,k,i,3) = (double)tur10[2];
                        tursav(j,k,i,4) = (double)tur10[3];
                    }
        }

    } else {
        for (i = 1; i <= idim-1; i++)
            for (k = 1; k <= kdim-1; k++)
                for (j = 1; j <= jdim-1; j++)
                    vist3d(j,k,i) = 0.0;
    }

    // Overwrite with exact soln if doing MMS
    if (iexact_trunc != 0 || iexact_disc != 0) {
        // tursav is 4D but exact_turb_q expects 3D q — pass as 3D view of same data
        FortranArray3DRef<double> tursav3d(
            reinterpret_cast<double*>(&tursav(1,1,1,1)),
            jdim, kdim, idim);
        mms_ns::exact_turb_q(jdim, kdim, idim, x, y, z,
                             tursav3d, smin, vist3d,
                             iexact_trunc, iexact_disc);
    }
}

} // namespace initvist_ns
