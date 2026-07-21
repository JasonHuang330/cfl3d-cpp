// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "add2x.h"
#include "fill.h"
#include "termn8.h"
#include "ccomplex.h"
#include <cstring>
#include <algorithm>
#include <cmath>

namespace add2x_ns {

void add2x(FortranArray4DRef<double> q, FortranArray4DRef<double> qc,
           int& jdim, int& kdim, int& idim, int& jj2, int& kk2, int& ii2,
           FortranArray4DRef<double> q1, FortranArray4DRef<double> dq,
           FortranArray4DRef<double> wq, FortranArray3DRef<double> wqj,
           FortranArray3DRef<double> wqjk, int& js, int& ks, int& is,
           int& je, int& ke, int& ie, int& ipass, int& nbl, int& nblc,
           FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim, int& ll, int& myid)
{
    // COMMON /mgrd/
    int& mode = cmn_mgrd.mode;

    // Note:  mode=1 must ONLY be used for ll=5 (primitive variables):
    if (mode != 0 && ll != 5) {
        nou(1) = std::min(nou(1) + 1, ibufdim);
        std::snprintf(bou(nou(1), 1), 120,
                      "mode must = 0 when ll .ne. 5 in addx");
        int ierrflg = -1;
        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
    }

    //      jdim,kdim,idim  finer mesh
    //      jj2,kk2,ii2     coarser mesh
    //      js,ks,is        coarser mesh starting indices
    //      je,ke,ie        coarser mesh ending indices

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;
    int jdim2 = jdim - 2;
    int kdim2 = kdim - 2;
    int idim2 = idim - 2;
    int jjl   = jj2 - 1;
    int kkl   = kk2 - 1;
    int iil   = ii2 - 1;

    //      semi-coarsening / directional refinement
    int nsi  = (idim - 1) / (ie - is);
    int ista = is - 1;
    int iend = ie;
    if (nsi == 1) {
        ista = is;
        iend = ie - 1;
    }

    // Initialize dq to zero
    int nn = jdim * kdim;
    for (int i = 1; i <= idim1; i++) {
        for (int n = 1; n <= ll; n++) {
            for (int izz = 1; izz <= nn; izz++) {
                dq(izz, 1, i, n) = 0.0e0;
            }
        }
    }

    if (ipass == 1) {
        //      wq=qc     mode=0    coarser grid
        //      wq=qc-q1  mode=1    coarser grid
        if (mode == 0) {
            nn = jj2 * kk2;
            for (int i = 1; i <= ii2; i++) {
                for (int n = 1; n <= ll; n++) {
                    for (int izz = 1; izz <= nn; izz++) {
                        wq(izz, 1, i, n) = qc(izz, 1, i, n);
                    }
                }
            }
        } else {
            nn = jj2 * kk2;
            for (int i = 1; i <= ii2; i++) {
                for (int n = 1; n <= ll; n++) {
                    for (int izz = 1; izz <= nn; izz++) {
                        wq(izz, 1, i, n) = qc(izz, 1, i, n) - q1(izz, 1, i, n);
                    }
                }
            }
        }
    }

    for (int n = 1; n <= ll; n++) {

        //      interpolate in j onto finer mesh
        int ksta = std::max(ks - 1, 1);
        for (int i = ista; i <= iend; i++) {
            for (int k = ksta; k <= ke; k++) {
                int jj = 0;
                for (int j = js; j <= je - 1; j++) {
                    int jm = std::max(1, j - 1);
                    int jp = std::min(jjl, j + 1);
                    jj = jj + 1;
                    wqj(jj, k, i) = 0.75e0 * wq(j, k, i, n) + 0.25e0 * wq(jm, k, i, n);
                    jj = jj + 1;
                    wqj(jj, k, i) = 0.75e0 * wq(j, k, i, n) + 0.25e0 * wq(jp, k, i, n);
                }
            }
        }

        //      interpolate in k onto finer mesh
        for (int i = ista; i <= iend; i++) {
            int kk = 0;
            for (int k = ks; k <= ke - 1; k++) {
                int km = std::max(1, k - 1);
                int kp = std::min(kkl, k + 1);
                kk = kk + 1;
                for (int j = 1; j <= jdim1; j++) {
                    wqjk(j, kk, i) = 0.75e0 * wqj(j, k, i) + 0.25e0 * wqj(j, km, i);
                }
                kk = kk + 1;
                for (int j = 1; j <= jdim1; j++) {
                    wqjk(j, kk, i) = 0.75e0 * wqj(j, k, i) + 0.25e0 * wqj(j, kp, i);
                }
            }
        }

        //      interpolate in i
        int np = jdim * kdim1 - 1;
        if (nsi == 2) {
            int ii = 0;
            for (int i = is; i <= ie - 1; i++) {
                int im = std::max(1, i - 1);
                int ip = std::min(iil, i + 1);
                ii = ii + 1;
                for (int izz = 1; izz <= np; izz++) {
                    dq(izz, 1, ii, n) = 0.75e0 * wqjk(izz, 1, i) + 0.25e0 * wqjk(izz, 1, im);
                }
                ii = ii + 1;
                for (int izz = 1; izz <= np; izz++) {
                    dq(izz, 1, ii, n) = 0.75e0 * wqjk(izz, 1, i) + 0.25e0 * wqjk(izz, 1, ip);
                }
            }
        } else {
            int ii = 0;
            for (int i = is; i <= ie - 1; i++) {
                ii = ii + 1;
                for (int izz = 1; izz <= np; izz++) {
                    dq(izz, 1, ii, n) = wqjk(izz, 1, i);
                }
            }
        }

    } // end do 80 n=1,ll

    //      q=dq    mode=0      fine grid    interpolation
    //      q=q+dq  mode=1      fine grid    correction
    if (mode == 0) {
        for (int i = 1; i <= idim1; i++) {
            nn = jdim * kdim - jdim - 1;
            for (int n = 1; n <= ll; n++) {
                for (int izz = 1; izz <= nn; izz++) {
                    q(izz, 1, i, n) = dq(izz, 1, i, n);
                }
            }
        }
    } else {
        for (int i = 1; i <= idim1; i++) {
            nn = jdim * kdim - jdim - 1;

            //         update density and pressure to ensure positivity
            double alpq = -.2;
            double phiq = 1. / 0.5;
            double betq = 1. + alpq * phiq;
            for (int izz = 1; izz <= nn; izz++) {
                double t1 = dq(izz, 1, i, 1) / q(izz, 1, i, 1);
                double t2 = dq(izz, 1, i, 1) / (betq + std::abs(t1) * phiq);
                bool cond1 = (t1 < alpq);
                dq(izz, 1, i, 1) = ccomplex_ns::ccvmgt(t2, dq(izz, 1, i, 1), cond1);

                t1 = dq(izz, 1, i, 5) / q(izz, 1, i, 5);
                t2 = dq(izz, 1, i, 5) / (betq + std::abs(t1) * phiq);
                bool cond5 = (t1 < alpq);
                dq(izz, 1, i, 5) = ccomplex_ns::ccvmgt(t2, dq(izz, 1, i, 5), cond5);
            }

            //         update primitive variables
            for (int izz = 1; izz <= nn; izz++) {
                q(izz, 1, i, 1) = q(izz, 1, i, 1) + dq(izz, 1, i, 1);
                q(izz, 1, i, 2) = q(izz, 1, i, 2) + dq(izz, 1, i, 2);
                q(izz, 1, i, 3) = q(izz, 1, i, 3) + dq(izz, 1, i, 3);
                q(izz, 1, i, 4) = q(izz, 1, i, 4) + dq(izz, 1, i, 4);
                q(izz, 1, i, 5) = q(izz, 1, i, 5) + dq(izz, 1, i, 5);
            }
        }
    }

    fill_ns::fill(jdim, kdim, idim, q, ll);
}

} // namespace add2x_ns
