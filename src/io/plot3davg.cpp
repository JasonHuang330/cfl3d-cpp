// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "plot3davg.h"
#include "plot3d.h"
#include "cctogp.h"
#include "runtime/common_blocks.h"
#include "runtime/fortran_array.h"
#include <cmath>
#include <cstdint>

namespace plot3davg_ns {

void plot3d(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> xw, FortranArray4DRef<double> blank2, FortranArray3DRef<double> blank, FortranArray4DRef<double> xg, int& iflag, FortranArray3DRef<double> vist3d, int& iover, int& nblk, FortranArray1DRef<int> nmap, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& ifunc, int& iplot, int& jdw, int& kdw, int& idw, int& nplots, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay, int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d, int& nblock, FortranArray1DRef<int> nblkpt, FortranArray4DRef<double> xv, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, int& nset)
{
    plot3d_ns::plot3d(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                      q, qi0, qj0, qk0, x, y, z, xw, blank2, blank, xg,
                      iflag, vist3d, iover, nblk, nmap, bcj, bck, bci,
                      vj0, vk0, vi0, ifunc, iplot, jdw, kdw, idw,
                      nplots, jdimg, kdimg, idimg, nblcg, jsg, ksg, isg,
                      jeg, keg, ieg, ninter, iindex, intmax, nsub1,
                      maxxe, nblkk, nbli, limblk, isva, nblon, mxbli,
                      thetay, maxbl, maxgr, myid, myhost, mycomm,
                      mblk2nd, inpl3d, nblock, nblkpt, xv, sj, sk, si,
                      vol, nset);
}

void plot3davg(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> xw, FortranArray4DRef<double> blank2, FortranArray3DRef<double> blank, FortranArray4DRef<double> xg, int& iflag, FortranArray3DRef<double> vist3d, int& iover, int& nblk, FortranArray1DRef<int> nmap, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& ifunc, int& iplot, int& jdw, int& kdw, int& idw, int& nplots, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay, int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d, int& nblock, FortranArray1DRef<int> nblkpt, FortranArray4DRef<double> xv, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, int& nset, FortranArray4DRef<double> qavg, FortranArray4DRef<double> q2avg, int& nt, int& movabs)
{
    // COMMON block aliases
    float& radtodeg   = cmn_conversion.radtodeg;
    float& xmach      = cmn_info.xmach;
    float& alpha      = cmn_info.alpha;
    float& reue       = cmn_reyue.reue;
    float& time       = cmn_unst.time;
    float& xnumavg    = cmn_avgdata.xnumavg;
    float& xnumavg2   = cmn_avgdata.xnumavg2;

    // Local variables
    int jdim1, kdim1, idim1;
    int jw, kw, iw;
    int ldw;
    double alphaw, xmachw, alphww, reuew, timew;
    int maxnum, minnum, maxpos, minpos;

    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    // initialize xw, xv and xg arrays
    jw = (j2 - j1) / j3 + 1;
    kw = (k2 - k1) / k3 + 1;
    iw = (i2 - i1) / i3 + 1;

    for (int j = 1; j <= jw; j++) {
        for (int k = 1; k <= kw; k++) {
            for (int i = 1; i <= iw; i++) {
                for (int l = 1; l <= 5; l++) {
                    xw(j, k, i, l) = 0.;
                    xv(j, k, i, l) = 0.;
                }
                for (int l = 1; l <= 4; l++) {
                    xg(j, k, i, l) = 0.;
                }
            }
        }
    }

    // assign single precision scalars
    alphaw = (double)radtodeg * ((double)alpha + (double)thetay(nblk));
    xmachw = (double)xmach;
    alphww = alphaw;
    reuew  = (double)reue;
    timew  = (double)time;

    // determine q values at grid points and load into xv array
    jdw = (j2 - j1) / j3 + 1;
    kdw = (k2 - k1) / k3 + 1;
    idw = (i2 - i1) / i3 + 1;
    ldw = 5;

    cctogp_ns::cctogp(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                      q, qi0, qj0, qk0, jdw, kdw, idw, xv, ldw);

    // get iteration-averaged Q values
    // note: the qavg values are kept as primitive variables

    // maxnum contains the max vals of each of the variables
    maxnum = 0;
    minnum = 99999999;

    // maxpos contains the indices of the maximum vals (4 is block #)
    maxpos = 0;
    minpos = 0;

    for (int i = 1; i <= idim; i++) {
        for (int j = 1; j <= jdim; j++) {
            for (int k = 1; k <= kdim; k++) {
                qavg(j, k, i, 1) = (qavg(j, k, i, 1) * ((double)xnumavg - 1.) +
                                    xv(j, k, i, 1)) / (double)xnumavg;
                qavg(j, k, i, 2) = (qavg(j, k, i, 2) * ((double)xnumavg - 1.) +
                                    xv(j, k, i, 2)) / (double)xnumavg;
                qavg(j, k, i, 3) = (qavg(j, k, i, 3) * ((double)xnumavg - 1.) +
                                    xv(j, k, i, 3)) / (double)xnumavg;
                qavg(j, k, i, 4) = (qavg(j, k, i, 4) * ((double)xnumavg - 1.) +
                                    xv(j, k, i, 4)) / (double)xnumavg;
                qavg(j, k, i, 5) = (qavg(j, k, i, 5) * ((double)xnumavg - 1.) +
                                    xv(j, k, i, 5)) / (double)xnumavg;

                q2avg(j, k, i, 1) = (q2avg(j, k, i, 1) * ((double)xnumavg2 - 1.) +
                                     xv(j, k, i, 1) * xv(j, k, i, 1)) / (double)xnumavg2;
                q2avg(j, k, i, 2) = (q2avg(j, k, i, 2) * ((double)xnumavg2 - 1.) +
                                     xv(j, k, i, 2) * xv(j, k, i, 2)) / (double)xnumavg2;
                q2avg(j, k, i, 3) = (q2avg(j, k, i, 3) * ((double)xnumavg2 - 1.) +
                                     xv(j, k, i, 3) * xv(j, k, i, 3)) / (double)xnumavg2;
                q2avg(j, k, i, 4) = (q2avg(j, k, i, 4) * ((double)xnumavg2 - 1.) +
                                     xv(j, k, i, 4) * xv(j, k, i, 4)) / (double)xnumavg2;
                q2avg(j, k, i, 5) = (q2avg(j, k, i, 5) * ((double)xnumavg2 - 1.) +
                                     xv(j, k, i, 5) * xv(j, k, i, 5)) / (double)xnumavg2;
            }
        }
    }

    return;
}

} // namespace plot3davg_ns
