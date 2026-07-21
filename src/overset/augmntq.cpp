// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "augmntq.h"
#include "runtime/common_blocks.h"

namespace augmntq_ns {

void augmntq(FortranArray4DRef<double> q, int& jdim, int& kdim, int& idim, int& nbl, int& ldim, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qq, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& ighost)
{
    // COMMON /twod/
    int& i2d = cmn_twod.i2d;

    // load interior q data into qq array
    for (int l = 1; l <= ldim; l++) {
        for (int i = 1; i <= idim-1; i++) {
            int ii = i + 1;
            for (int k = 1; k <= kdim-1; k++) {
                int kk = k + 1;
                for (int j = 1; j <= jdim-1; j++) {
                    int jj = j + 1;
                    qq(jj,kk,ii,l) = q(j,k,i,l);
                }
            }
        }
    }

    // fill in block faces (excluding for now the edges and corners)
    // if ghost-cell data is not available, use duplicate interior
    // data to fill in the block face info.
    //
    // note: eddy viscosity and turbulence quantities are always
    // stored at cell centers, so must override the bcj/k/i arrays

    double visflg = 1.;
    if (ldim != 5) visflg = 0.;

    if (ighost != 0) {

        // calculate and store cell-face center data into qq

        for (int l = 1; l <= ldim; l++) {
            for (int i = 1; i <= idim-1; i++) {
                int ii = i + 1;
                for (int k = 1; k <= kdim-1; k++) {
                    int kk = k + 1;
                    double aa = 1. + bcj(k,i,1)*visflg;
                    double bb = 1. - bcj(k,i,1)*visflg;
                    double cc = 1. + bcj(k,i,2)*visflg;
                    double dd = 1. - bcj(k,i,2)*visflg;
                    qq(1,kk,ii,l)      = (aa*qj0(k,i,l,1) +
                                          bb*q(1,k,i,l))*0.5;
                    qq(jdim+1,kk,ii,l) = (cc*qj0(k,i,l,3) +
                                          dd*q(jdim-1,k,i,l))*0.5;
                }
            }
            for (int i = 1; i <= idim-1; i++) {
                int ii = i + 1;
                for (int j = 1; j <= jdim-1; j++) {
                    int jj = j + 1;
                    double aa = 1. + bck(j,i,1)*visflg;
                    double bb = 1. - bck(j,i,1)*visflg;
                    double cc = 1. + bck(j,i,2)*visflg;
                    double dd = 1. - bck(j,i,2)*visflg;
                    qq(jj,1,ii,l)      = (aa*qk0(j,i,l,1) +
                                          bb*q(j,1,i,l))*0.5;
                    qq(jj,kdim+1,ii,l) = (cc*qk0(j,i,l,3) +
                                          dd*q(j,kdim-1,i,l))*0.5;
                }
            }
            for (int k = 1; k <= kdim-1; k++) {
                int kk = k + 1;
                for (int j = 1; j <= jdim-1; j++) {
                    int jj = j + 1;
                    double aa = 1. + bci(j,k,1)*visflg;
                    double bb = 1. - bci(j,k,1)*visflg;
                    double cc = 1. + bci(j,k,2)*visflg;
                    double dd = 1. - bci(j,k,2)*visflg;
                    qq(jj,kk,1,l)      = (aa*qi0(j,k,l,1) +
                                          bb*q(j,k,1,l))*0.5;
                    qq(jj,kk,idim+1,l) = (cc*qi0(j,k,l,3) +
                                          dd*q(j,k,idim-1,l))*0.5;
                }
            }
        }

    } else {

        // load duplicate interior data in block faces (1st order approx)

        for (int l = 1; l <= ldim; l++) {
            for (int i = 1; i <= idim-1; i++) {
                int ii = i + 1;
                for (int k = 1; k <= kdim-1; k++) {
                    int kk = k + 1;
                    qq(1,kk,ii,l)      = qq(2,kk,ii,l);
                    qq(jdim+1,kk,ii,l) = qq(jdim,kk,ii,l);
                }
            }
            for (int i = 1; i <= idim-1; i++) {
                int ii = i + 1;
                for (int j = 1; j <= jdim-1; j++) {
                    int jj = j + 1;
                    qq(jj,1,ii,l)      = qq(jj,2,ii,l);
                    qq(jj,kdim+1,ii,l) = qq(jj,kdim,ii,l);
                }
            }
            for (int k = 1; k <= kdim-1; k++) {
                int kk = k + 1;
                for (int j = 1; j <= jdim-1; j++) {
                    int jj = j + 1;
                    qq(jj,kk,1,l)      = qq(jj,kk,2,l);
                    qq(jj,kk,idim+1,l) = qq(jj,kk,idim,l);
                }
            }
        }
    }

    // fill in edge values and corner values via extrapolation/averaging

    int itwo = 2;
    if (i2d != 0) itwo = 1;

    for (int l = 1; l <= ldim; l++) {

        // fill in j-k edge values via extrapolation in both j and k
        // directions and averaging

        for (int i = 2; i <= idim; i++) {
            // k loop: k=1 and k=kdim+1 (step kdim)
            for (int k = 1; k <= kdim+1; k += kdim) {
                int k1 = 1;
                int k2 = 2;
                if (k == kdim+1) {
                    k1 = -1;
                    k2 = -2;
                }
                // j loop: j=1 and j=jdim+1 (step jdim)
                for (int j = 1; j <= jdim+1; j += jdim) {
                    int j1 = 1;
                    int j2 = 2;
                    if (j == jdim+1) {
                        j1 = -1;
                        j2 = -2;
                    }
                    double dqj = qq(j+j1,k,i,l) - qq(j+j2,k,i,l);
                    double qj  = qq(j+j1,k,i,l) + dqj;
                    double dqk = qq(j,k+k1,i,l) - qq(j,k+k2,i,l);
                    double qk  = qq(j,k+k1,i,l) + dqk;
                    qq(j,k,i,l) = (qj + qk)*0.5;
                }
            }
        }

        // fill in i-j edge values via extrapolation in both i and j
        // directions and averaging

        // i loop: i=1 and i=idim+1 (step idim)
        for (int i = 1; i <= idim+1; i += idim) {
            int i1 = 1;
            int i2 = itwo;
            if (i == idim+1) {
                i1 = -1;
                i2 = -itwo;
            }
            for (int k = 2; k <= kdim; k++) {
                // j loop: j=1 and j=jdim+1 (step jdim)
                for (int j = 1; j <= jdim+1; j += jdim) {
                    int j1 = 1;
                    int j2 = 2;
                    if (j == jdim+1) {
                        j1 = -1;
                        j2 = -2;
                    }
                    double dqj = qq(j+j1,k,i,l) - qq(j+j2,k,i,l);
                    double qj  = qq(j+j1,k,i,l) + dqj;
                    double dqi = qq(j,k,i+i1,l) - qq(j,k,i+i2,l);
                    double qi  = qq(j,k,i+i1,l) + dqi;
                    qq(j,k,i,l) = (qj + qi)*0.5;
                }
            }
        }

        // fill in i-k edge values via extrapolation in both i and k
        // directions and averaging

        // i loop: i=1 and i=idim+1 (step idim)
        for (int i = 1; i <= idim+1; i += idim) {
            int i1 = 1;
            int i2 = itwo;
            if (i == idim+1) {
                i1 = -1;
                i2 = -itwo;
            }
            // k loop: k=1 and k=kdim+1 (step kdim)
            for (int k = 1; k <= kdim+1; k += kdim) {
                int k1 = 1;
                int k2 = 2;
                if (k == kdim+1) {
                    k1 = -1;
                    k2 = -2;
                }
                for (int j = 2; j <= jdim; j++) {
                    double dqi = qq(j,k,i+i1,l) - qq(j,k,i+i2,l);
                    double qi  = qq(j,k,i+i1,l) + dqi;
                    double dqk = qq(j,k+k1,i,l) - qq(j,k+k2,i,l);
                    double qk  = qq(j,k+k1,i,l) + dqk;
                    qq(j,k,i,l) = (qi + qk)*0.5;
                }
            }
        }

        // fill in corner values via extrapolation in all 3 directions
        // and averaging

        // k loop: k=1 and k=kdim+1 (step kdim)
        for (int k = 1; k <= kdim+1; k += kdim) {
            int k1 = 1;
            int k2 = 2;
            if (k == kdim+1) {
                k1 = -1;
                k2 = -2;
            }
            // j loop: j=1 and j=jdim+1 (step jdim)
            for (int j = 1; j <= jdim+1; j += jdim) {
                int j1 = 1;
                int j2 = 2;
                if (j == jdim+1) {
                    j1 = -1;
                    j2 = -2;
                }
                // i loop: i=1 and i=idim+1 (step idim)
                for (int i = 1; i <= idim+1; i += idim) {
                    int i1 = 1;
                    int i2 = itwo;
                    if (i == idim+1) {
                        i1 = -1;
                        i2 = -itwo;
                    }
                    double dqi = qq(j,k,i+i1,l) - qq(j,k,i+i2,l);
                    double qi  = qq(j,k,i+i1,l) + dqi;
                    double dqj = qq(j+j1,k,i,l) - qq(j+j2,k,i,l);
                    double qj  = qq(j+j1,k,i,l) + dqj;
                    double dqk = qq(j,k+k1,i,l) - qq(j,k+k2,i,l);
                    double qk  = qq(j,k+k1,i,l) + dqk;
                    qq(j,k,i,l) = (qi + qj + qk)/3.;
                }
            }
        }
    }

    return;
}

} // namespace augmntq_ns
