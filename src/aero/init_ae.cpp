// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "init_ae.h"
#include "init.h"
#include <cmath>

namespace init_ae_ns {

void init(int& nbl, int& jdim, int& kdim, int& idim,
          FortranArray4DRef<double> q,
          FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0,
          FortranArray4DRef<double> tj0, FortranArray4DRef<double> tk0, FortranArray4DRef<double> ti0,
          FortranArray3DRef<double> vol,
          FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0,
          int& nummem,
          FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z,
          FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
          int& nbuf, int& ibufdim, int& iflagprnt)
{
    init_ns::init(nbl, jdim, kdim, idim, q, qj0, qk0, qi0, tj0, tk0, ti0,
                  vol, volj0, volk0, voli0, nummem, x, y, z,
                  nou, bou, nbuf, ibufdim, iflagprnt);
}

void init_ae(int& nmds, int& maxaes,
             FortranArray2DRef<double> x0,
             FortranArray2DRef<double> xxn,
             FortranArray2DRef<double> temp,
             FortranArray3DRef<double> bmat,
             FortranArray3DRef<double> stm,
             FortranArray3DRef<double> stmi,
             FortranArray2DRef<double> gforcn,
             FortranArray2DRef<double> gforcnm,
             FortranArray2DRef<double> freq,
             FortranArray2DRef<double> damp,
             FortranArray2DRef<double> gmass,
             FortranArray2DRef<double> aesrfdat,
             FortranArray2DRef<double> xs,
             FortranArray2DRef<double> gforcs)
{
    // COMMON block references
    int32_t& naesrf = cmn_elastic.naesrf;
    float& xmach    = cmn_info.xmach;
    float& dt       = cmn_info.dt;

    for (int iaes = 1; iaes <= naesrf; iaes++) {

        int    iskyhk  = (int)aesrfdat(1, iaes);
        double grefl   = aesrfdat(2, iaes);
        double uinf    = aesrfdat(3, iaes);
        double qinf    = aesrfdat(4, iaes);
        int    nmodes  = (int)aesrfdat(5, iaes);
        double ainf    = uinf / (double)xmach;
        double rgrefl  = 1.0 / grefl;

        for (int m = 1; m <= 2 * nmodes; m++) {
            for (int n = 1; n <= 2 * nmodes; n++) {
                stm(n, m, iaes)    = 0.0;
                stmi(n, m, iaes)   = 0.0;
                bmat(n, m, iaes)   = 0.0;
                gforcn(m, iaes)    = 0.0;
                gforcnm(m, iaes)   = 0.0;
                xxn(m, iaes)       = x0(m, iaes);
            }
        }

        for (int i = 1; i <= 2 * nmodes - 1; i += 2) {

            int j = (i + 1) / 2;

            double aa = -freq(j, iaes) * damp(j, iaes);
            double bb =  freq(j, iaes) * std::sqrt(1.0 - damp(j, iaes) * damp(j, iaes));
            double ec = std::exp(aa * (double)dt * grefl / ainf) * std::cos(bb * (double)dt * grefl / ainf);
            double es = std::exp(aa * (double)dt * grefl / ainf) * std::sin(bb * (double)dt * grefl / ainf);

            stm(i,   i,   iaes) =   ec - aa * es / bb;
            stm(i,   i+1, iaes) =   es / bb;
            stm(i+1, i,   iaes) = -(aa * aa + bb * bb) * es / bb;
            stm(i+1, i+1, iaes) =   ec + aa * es / bb;

            stmi(i,   i,   iaes) = (2.0 * aa * (ec - 1.0)
                                   + (bb - aa * aa / bb) * es)
                                   / (aa * aa + bb * bb);
            stmi(i,   i+1, iaes) = (aa * es / bb - ec + 1.0)
                                   / (aa * aa + bb * bb);
            stmi(i+1, i,   iaes) = -aa * es / bb + ec - 1.0;
            stmi(i+1, i+1, iaes) =  es / bb;

            bmat(i+1, i+1, iaes) = 1.0 / gmass(j, iaes);

        }

        // overwrite stmi with matrix product stmi*bmat
        // (theta*B in the reference's notation)
        for (int j = 1; j <= 2 * nmodes; j++) {
            for (int i = 1; i <= 2 * nmodes; i++) {
                temp(i, j) = 0.0;
                for (int k = 1; k <= 2 * nmodes; k++) {
                    temp(i, j) = temp(i, j)
                                 + stmi(i, k, iaes) * bmat(k, j, iaes);
                }
            }
        }
        for (int j = 1; j <= 2 * nmodes; j++) {
            for (int i = 1; i <= 2 * nmodes; i++) {
                stmi(i, j, iaes) = temp(i, j);
            }
        }

    } // end do iaes

    return;
}

} // namespace init_ae_ns
