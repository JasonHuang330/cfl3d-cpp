// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "af3f.h"
#include "termn8.h"
#include "gfluxl.h"
#include "amafj.h"
#include "abcjz.h"
#include "vlutr.h"
#include "vlutrp.h"
#include "swafj.h"
#include "diagj.h"
#include "l2norm.h"
#include "hfluxl.h"
#include "amafk.h"
#include "hfluxv.h"
#include "abckz.h"
#include "swafk.h"
#include "diagk.h"
#include "ffluxl.h"
#include "amafi.h"
#include "abciz.h"
#include "swafi.h"
#include "diagi.h"
#include "runtime/fortran_io.h"
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cmath>

namespace af3f_ns {

void af3f(int& nbl, int& jdim, int& kdim, int& idim,
          FortranArray2DRef<double> q, FortranArray1DRef<double> vol,
          FortranArray3DRef<double> qj0, FortranArray3DRef<double> qk0,
          FortranArray3DRef<double> qi0, FortranArray1DRef<double> dtj,
          FortranArray2DRef<double> sj, FortranArray2DRef<double> sk,
          FortranArray2DRef<double> si, FortranArray2DRef<double> res,
          FortranArray3DRef<double> vist3d, FortranArray3DRef<double> x,
          FortranArray3DRef<double> y, FortranArray3DRef<double> z,
          FortranArray3DRef<double> blank, FortranArray3DRef<double> vmuk,
          double& resid, FortranArray1DRef<double> wk, int& nwork,
          FortranArray1DRef<double> wk0, int& nwk0, int& iover,
          FortranArray4DRef<double> vk0, FortranArray3DRef<double> bcj,
          FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
          FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
          int& nbuf, int& ibufdim, int& myid,
          FortranArray1DRef<int> mblk2nd, int& maxbl,
          FortranArray3DRef<double> volk0, FortranArray4DRef<double> tursav,
          FortranArray4DRef<double> tk0, FortranArray3DRef<double> cmuv,
          FortranArray1DRef<int> iadvance, int& nummem,
          FortranArray4DRef<double> ux)
{
    // COMMON block references
    int32_t& isklton  = cmn_sklton.isklton;
    int32_t& ita      = cmn_unst.ita;
    int32_t* ifds     = cmn_fvfds.ifds;   // 0-based C array, ifds(1..3) -> ifds[0..2]
    float*   rkap     = cmn_info.rkap;    // 0-based C array, rkap(1..3) -> rkap[0..2]
    int32_t* idiag    = cmn_info.idiag;   // 0-based C array, idiag(1..3) -> idiag[0..2]
    int32_t* ivisc    = cmn_reyue.ivisc;  // 0-based C array, ivisc(1..3) -> ivisc[0..2]
    int32_t& icyc     = cmn_mgrd.icyc;
    int32_t& ncyc     = cmn_mgrd.ncyc;
    int32_t& i2d      = cmn_twod.i2d;

    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    if (isklton > 0) {
        if (std::abs(ita) == 2) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                          "   using second order time differencing");
        } else {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                          "   using first order time differencing");
        }
    }

    //
    // ****************** spatially-split af ***********************
    //
    //      idiag(l)=1  diagonal inversion for l direction
    //      idiag(l)=0  block inversion
    //
    // cc     implicit   J direction
    //
    int imult = 1;
    if (ifds[1] == 0 || idiag[1] == 0) {
        if (isklton > 0) nou(1) = std::min(nou(1) + 1, ibufdim);
        if (isklton > 0) std::snprintf(bou(nou(1), 1), 120,
                          "   5x5 block inversion in J-direction");
        int nvtq = std::min(999000, nwork / 145);
        int n    = jdim * kdim;
        int nplq = std::min(idim1, nvtq / n);
        if (imult == 0) nplq = 1;
        int npl = nplq;
        if (npl < 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                          " insufficient memory - stopping in af3f");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        int iperd = 0;
        for (int i = 1; i <= idim1; i += nplq) {
            if (i + npl - 1 > idim1) npl = idim1 - i + 1;

            nvtq = npl * jdim * kdim;
            int iwka = nvtq * 20 + 1;
            int iwkb = nvtq * 45 + 1;
            int iwkc = nvtq * 70 + 1;
            int iwkd = nvtq * 95 + 1;
            int iwke = nvtq * 120 + 1;

            // gfluxl declares q(jdim,kdim,idim,5), res(jdim,kdim,idim-1,5),
            // qj0(kdim,idim-1,5,4), sj(jdim*kdim,idim-1,5)
            FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
            FortranArray4DRef<double> q4d(&q(1, 1), jdim, kdim, idim, 5);
            FortranArray4DRef<double> qj04d(&qj0(1, 1, 1), kdim, idim - 1, 5, 4);
            FortranArray3DRef<double> sj3d(&sj(1, 1), jdim * kdim, idim - 1, 5);
            FortranArray3DRef<double> dgp_j(&wk(iwkd), nvtq, 5, 1);
            FortranArray3DRef<double> dgm_j(&wk(iwke), nvtq, 5, 1);
            FortranArray2DRef<double> wk2d(&wk(1), nvtq, 145);
            double rkap2 = (double)rkap[1];
            int idf1 = 1;

            gfluxl_ns::gfluxl(i, npl, rkap2, idf1, jdim, kdim, idim,
                               res4d, q4d, qj04d, sj3d, dgp_j, dgm_j, wk2d, nvtq);

            FortranArray4DRef<double> aj(&wk(iwka), nvtq, 5, 5, 1);
            FortranArray4DRef<double> bj(&wk(iwkb), nvtq, 5, 5, 1);
            FortranArray4DRef<double> cj(&wk(iwkc), nvtq, 5, 5, 1);
            // amafj declares q(jdim,kdim,idim,5), dtj(jdim,kdim,idim-1)
            FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim - 1);
            FortranArray4DRef<double> dgp_j4d(&wk(iwkd), nvtq, 5, 1, 1);
            FortranArray4DRef<double> dgm_j4d(&wk(iwke), nvtq, 5, 1, 1);

            amafj_ns::amafj(i, npl, jdim, kdim, idim, q4d, aj, bj, cj,
                            dtj3d, wk2d, nvtq, dgp_j4d, dgm_j4d);

            if (iover == 1) {
                FortranArray5DRef<double> aj5d(&wk(iwka), kdim-1, npl, jdim, 5, 5);
                FortranArray5DRef<double> bj5d(&wk(iwkb), kdim-1, npl, jdim, 5, 5);
                FortranArray5DRef<double> cj5d(&wk(iwkc), kdim-1, npl, jdim, 5, 5);
                abcjz_ns::abcjz(i, npl, jdim, kdim, idim, aj5d, bj5d, cj5d, blank);
            }

            int nvmax = npl * kdim1;
            if (iperd == 0) {
                int il1 = 1;
                vlutr_ns::vlutr(nvmax, nvmax, jdim, il1, jdim1, aj, bj, cj,
                                nou, bou, nbuf, ibufdim);
            } else {
                if (isklton == 1) nou(1) = std::min(nou(1) + 1, ibufdim);
                if (isklton == 1) std::snprintf(bou(nou(1), 1), 120,
                                   " periodic matrix equation in J-direction");
                FortranArray4DRef<double> gj(&wk(iwkd), nvtq, 5, 5, 1);
                FortranArray4DRef<double> hj(&wk(iwke), nvtq, 5, 5, 1);
                int il1 = 1;
                vlutrp_ns::vlutrp(nvmax, nvmax, jdim, il1, jdim1, aj, bj, cj, gj, hj);
            }

            FortranArray3DRef<double> wk3d_f(&wk(1), nvtq, 145, 1);
            FortranArray4DRef<double> gj_sw(&wk(iwkd), nvtq, 5, 5, 1);
            FortranArray4DRef<double> hj_sw(&wk(iwke), nvtq, 5, 5, 1);
            swafj_ns::swafj(i, npl, jdim, kdim, idim, aj, bj, cj,
                            wk3d_f, nvtq, res4d, iperd, gj_sw, hj_sw);
        }
    } else {
        // diagonal inversion in J-direction
        if (isklton > 0) nou(1) = std::min(nou(1) + 1, ibufdim);
        if (isklton > 0) std::snprintf(bou(nou(1), 1), 120,
                          "   diagonal inversion in J-direction");
        int nvtq = std::min(999000, nwork / 35);
        int n    = jdim * kdim;
        int nplq = std::min(idim1, nvtq / n);
        if (imult == 0) nplq = 1;
        int npl  = nplq;
        if (npl < 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                          " insufficient memory - stopping in af3f");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        int iperd = 0;
        if (isklton > 0 && iperd == 1) nou(1) = std::min(nou(1) + 1, ibufdim);
        if (isklton > 0 && iperd == 1)
            std::snprintf(bou(nou(1), 1), 120,
                          "  periodic matrix equation in J-direction");
        for (int i = 1; i <= idim1; i += nplq) {
            if (i + npl - 1 > idim1) npl = idim1 - i + 1;
            // diagj declares q(jdim,kdim,idim,5), res(jdim,kdim,idim-1,5),
            // sj(jdim,kdim,idim-1,5), dtj(jdim,kdim,idim-1), vol(jdim,kdim,idim-1)
            FortranArray4DRef<double> q4d(&q(1, 1), jdim, kdim, idim, 5);
            FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
            FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim - 1);
            FortranArray4DRef<double> sj4d(&sj(1, 1), jdim, kdim, idim - 1, 5);
            FortranArray2DRef<double> wk2d(&wk(1), npl * kdim * jdim, 35);
            FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim - 1);
            diagj_ns::diagj(i, npl, jdim, kdim, idim, q4d, res4d, dtj3d, sj4d,
                            wk2d, iperd, vol3d, vist3d, blank, iover);
        }
    }

    //
    // l2-norm of delta q
    //
    if (icyc <= 2 || icyc == ncyc) {
        // l2norm declares res(jdim,kdim,idim-1,5), vol(jdim,kdim,idim-1)
        FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
        FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim - 1);
        int irdq1 = 1;
        l2norm_ns::l2norm(nbl, icyc, resid, irdq1, jdim, kdim, idim, res4d, vol3d);
    }


    //
    // cc     implicit   K direction
    //
    if (ifds[2] == 0 || idiag[2] == 0) {
        if (isklton > 0) nou(1) = std::min(nou(1) + 1, ibufdim);
        if (isklton > 0) std::snprintf(bou(nou(1), 1), 120,
                          "   5x5 block inversion in K-direction");
        int nvtq = std::min(999000, nwork / 145);
        int n    = jdim * kdim;
        int nplq = std::min(idim1, nvtq / n);
        if (imult == 0) nplq = 1;
        int npl  = nplq;
        if (npl < 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                          " insufficient memory - stopping in af3f");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        int jz = npl * jdim;
        (void)jz; // jz computed but not used further in this scope
        for (int il = 1; il <= idim1; il += nplq) {
            if (il + npl - 1 > idim1) npl = idim1 - il + 1;
            int i = idim1 - il - npl + 2;

            nvtq = npl * jdim * kdim;
            int iwka = nvtq * 20 + 1;
            int iwkb = nvtq * 45 + 1;
            int iwkc = nvtq * 70 + 1;
            int iwkd = nvtq * 95 + 1;
            int iwke = nvtq * 120 + 1;

            // hfluxl declares q(jdim,kdim,idim,5), res(jdim,kdim,idim-1,5),
            // qk0(jdim,idim-1,5,4), sk(jdim*kdim,idim-1,5)
            FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
            FortranArray4DRef<double> q4d(&q(1, 1), jdim, kdim, idim, 5);
            FortranArray4DRef<double> qk04d(&qk0(1, 1, 1), jdim, idim - 1, 5, 4);
            FortranArray3DRef<double> sk3d(&sk(1, 1), jdim * kdim, idim - 1, 5);
            FortranArray3DRef<double> dhp_k(&wk(iwkd), nvtq, 5, 1);
            FortranArray3DRef<double> dhm_k(&wk(iwke), nvtq, 5, 1);
            FortranArray2DRef<double> wk2d(&wk(1), nvtq, 145);
            double rkap3 = (double)rkap[2];
            int idf1 = 1;

            hfluxl_ns::hfluxl(i, npl, rkap3, idf1, jdim, kdim, idim,
                               res4d, q4d, qk04d, sk3d, dhp_k, dhm_k, wk2d, nvtq);

            FortranArray4DRef<double> ak(&wk(iwka), nvtq, 5, 5, 1);
            FortranArray4DRef<double> bk(&wk(iwkb), nvtq, 5, 5, 1);
            FortranArray4DRef<double> ck(&wk(iwkc), nvtq, 5, 5, 1);
            // amafk declares q(jdim,kdim,idim,5), dtj(jdim,kdim,idim-1)
            FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim - 1);
            FortranArray4DRef<double> dhp_k4d(&wk(iwkd), nvtq, 5, 1, 1);
            FortranArray4DRef<double> dhm_k4d(&wk(iwke), nvtq, 5, 1, 1);

            amafk_ns::amafk(i, npl, jdim, kdim, idim, q4d, ak, bk, ck,
                            dtj3d, wk2d, nvtq, dhp_k4d, dhm_k4d);

            if (ivisc[2] > 0) {
                // Call viscous LHS.
                FortranArray2DRef<double> wk02d(&wk0(1), nwk0, 1);
                // hfluxv declares vol(jdim,kdim,idim-1)
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim - 1);
                FortranArray2DRef<double> t_hv(&wk(1), nvtq, 145);
                int iadv_nbl = iadvance(nbl);
                hfluxv_ns::hfluxv(i, npl, jdim, kdim, idim, idf1,
                                  ak, bk, ck, res4d, q4d, qk04d, sk3d, vol3d,
                                  t_hv, nvtq, wk02d, vist3d, vmuk, vk0, bck,
                                  tursav, tk0, cmuv, volk0,
                                  nou, bou, nbuf, ibufdim, iadv_nbl, nummem, ux);
            }

            if (iover == 1) {
                FortranArray5DRef<double> ak5d(&wk(iwka), jdim-1, npl, kdim-1, 5, 5);
                FortranArray5DRef<double> bk5d(&wk(iwkb), jdim-1, npl, kdim-1, 5, 5);
                FortranArray5DRef<double> ck5d(&wk(iwkc), jdim-1, npl, kdim-1, 5, 5);
                abckz_ns::abckz(i, npl, jdim, kdim, idim, ak5d, bk5d, ck5d, blank);
            }

            int imw = 0;
            int nvmax = npl * (jdim - 1) / (imw + 1);
            int nrec  = kdim1 * (imw + 1);
            int il1 = 1;
            vlutr_ns::vlutr(nvmax, nvmax, nrec, il1, nrec, ak, bk, ck,
                            nou, bou, nbuf, ibufdim);
            // swafk declares q(jdim,kdim,idim,5), res(jdim,kdim,idim-1,5),
            // dtj(jdim,kdim,idim-1)
            FortranArray3DRef<double> dtj3d_sw(&dtj(1), jdim, kdim, idim - 1);
            FortranArray3DRef<double> f3d_sw(&wk(1), nvtq, 145, 1);
            swafk_ns::swafk(i, npl, jdim, kdim, idim, q4d, ak, bk, ck,
                            dtj3d_sw, f3d_sw, nvtq, res4d, imw);
        }
    } else {
        // diagonal inversion in K-direction
        if (isklton > 0) nou(1) = std::min(nou(1) + 1, ibufdim);
        if (isklton > 0) std::snprintf(bou(nou(1), 1), 120,
                          "   diagonal inversion in K-direction");
        int nvtq = std::min(999000, nwork / 35);
        int n    = jdim * kdim;
        int nplq = std::min(idim1, nvtq / n);
        if (imult == 0) nplq = 1;
        int npl  = nplq;
        if (npl < 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                          " insufficient memory - stopping in af3f");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        for (int i = 1; i <= idim1; i += nplq) {
            if (i + npl - 1 > idim1) npl = idim1 - i + 1;
            // diagk declares q(jdim,kdim,idim,5), res(jdim,kdim,idim-1,5),
            // sk(jdim,kdim,idim-1,5), dtj(jdim,kdim,idim-1), vol(jdim,kdim,idim-1)
            FortranArray4DRef<double> q4d(&q(1, 1), jdim, kdim, idim, 5);
            FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
            FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim - 1);
            FortranArray4DRef<double> sk4d(&sk(1, 1), jdim, kdim, idim - 1, 5);
            FortranArray2DRef<double> wk2d(&wk(1), npl * (jdim - 1) * kdim, 35);
            FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim - 1);
            diagk_ns::diagk(i, npl, jdim, kdim, idim, q4d, res4d, dtj3d, sk4d,
                            wk2d, vol3d, vist3d, blank, iover);
        }
    }

    //
    // l2-norm of delta q
    //
    if (icyc <= 2 || icyc == ncyc) {
        // l2norm declares res(jdim,kdim,idim-1,5), vol(jdim,kdim,idim-1)
        FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
        FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim - 1);
        int irdq1 = 1;
        l2norm_ns::l2norm(nbl, icyc, resid, irdq1, jdim, kdim, idim, res4d, vol3d);
    }


    //
    // cc     implicit in I direction
    //
    if (i2d != 1) {
        if (ifds[0] == 0 || idiag[0] == 0) {
            if (isklton > 0) nou(1) = std::min(nou(1) + 1, ibufdim);
            if (isklton > 0) std::snprintf(bou(nou(1), 1), 120,
                              "   5x5 block inversion in I-direction");
            int n    = jdim * idim;
            int nvtq = std::min(999000, nwork / 145);
            int nplq = std::min(kdim1, nvtq / n);
            if (imult == 0) nplq = 1;
            int npl  = nplq;
            if (npl < 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                              " insufficient memory - stopping in af3f");
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
            for (int k = 1; k <= kdim1; k += nplq) {
                if (k + npl - 1 > kdim1) npl = kdim1 - k + 1;

                nvtq = npl * jdim * idim;
                int iwka = nvtq * 20 + 1;
                int iwkb = nvtq * 45 + 1;
                int iwkc = nvtq * 70 + 1;
                int iwkd = nvtq * 95 + 1;
                int iwke = nvtq * 120 + 1;

                // ffluxl declares q(jdim,kdim,idim,5), res(jdim,kdim,idim-1,5),
                // qi0(jdim,kdim,5,4), si(jdim*kdim,idim,5)
                FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
                FortranArray4DRef<double> q4d(&q(1, 1), jdim, kdim, idim, 5);
                FortranArray4DRef<double> qi04d(&qi0(1, 1, 1), jdim, kdim, 5, 4);
                FortranArray3DRef<double> si3d(&si(1, 1), jdim * kdim, idim, 5);
                FortranArray3DRef<double> dfp_i(&wk(iwkd), nvtq, 5, 1);
                FortranArray3DRef<double> dfm_i(&wk(iwke), nvtq, 5, 1);
                FortranArray2DRef<double> wk2d(&wk(1), nvtq, 145);
                double rkap1 = (double)rkap[0];
                int idf1 = 1;

                ffluxl_ns::ffluxl(k, npl, rkap1, idf1, jdim, kdim, idim,
                                   res4d, q4d, qi04d, si3d, dfp_i, dfm_i, wk2d, nvtq);

                int imw = 0;
                FortranArray4DRef<double> ai(&wk(iwka), nvtq, 5, 5, 1);
                FortranArray4DRef<double> bi(&wk(iwkb), nvtq, 5, 5, 1);
                FortranArray4DRef<double> ci(&wk(iwkc), nvtq, 5, 5, 1);
                // amafi declares q(jdim,kdim,idim,5), dtj(jdim,kdim,idim-1)
                FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim - 1);
                FortranArray4DRef<double> dfp_i4d(&wk(iwkd), nvtq, 5, 1, 1);
                FortranArray4DRef<double> dfm_i4d(&wk(iwke), nvtq, 5, 1, 1);

                amafi_ns::amafi(k, npl, jdim, kdim, idim, q4d, ai, bi, ci,
                                dtj3d, wk2d, nvtq, dfp_i4d, dfm_i4d, imw);

                if (iover == 1) {
                    FortranArray5DRef<double> ai5d(&wk(iwka), jdim-1, npl, idim-1, 5, 5);
                    FortranArray5DRef<double> bi5d(&wk(iwkb), jdim-1, npl, idim-1, 5, 5);
                    FortranArray5DRef<double> ci5d(&wk(iwkc), jdim-1, npl, idim-1, 5, 5);
                    abciz_ns::abciz(k, npl, jdim, kdim, idim, ai5d, bi5d, ci5d, blank);
                }

                int nvmax = npl * (jdim - 1) / (imw + 1);
                int nrec  = idim1 * (imw + 1);
                int il1 = 1;
                vlutr_ns::vlutr(nvmax, nvmax, nrec, il1, nrec, ai, bi, ci,
                                nou, bou, nbuf, ibufdim);

                // swafi declares q(jdim,kdim,idim,5), res(jdim,kdim,idim-1,5),
                // dtj(jdim,kdim,idim-1)
                FortranArray3DRef<double> dtj3d_sw(&dtj(1), jdim, kdim, idim - 1);
                FortranArray3DRef<double> f3d_sw(&wk(1), nvtq, 145, 1);
                swafi_ns::swafi(k, npl, jdim, kdim, idim, q4d, ai, bi, ci,
                                dtj3d_sw, f3d_sw, nvtq, res4d, imw);
            }
        } else {
            // diagonal inversion in I-direction
            if (isklton > 0) nou(1) = std::min(nou(1) + 1, ibufdim);
            if (isklton > 0) std::snprintf(bou(nou(1), 1), 120,
                              "   diagonal inversion in I-direction");
            int n    = jdim * idim;
            int nvtq = std::min(999000, nwork / 35);
            int nplq = std::min(kdim1, nvtq / n);
            if (imult == 0) nplq = 1;
            int npl  = nplq;
            if (npl < 1) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                              " insufficient memory - stopping in af3f");
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
            for (int k = 1; k <= kdim1; k += nplq) {
                if (k + npl - 1 > kdim1) npl = kdim1 - k + 1;
                // diagi declares q(jdim,kdim,idim,5), res(jdim,kdim,idim-1,5),
                // si(jdim,kdim,idim,5), dtj(jdim,kdim,idim-1), vol(jdim,kdim,idim-1)
                FortranArray4DRef<double> q4d(&q(1, 1), jdim, kdim, idim, 5);
                FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
                FortranArray3DRef<double> dtj3d(&dtj(1), jdim, kdim, idim - 1);
                FortranArray4DRef<double> si4d(&si(1, 1), jdim, kdim, idim, 5);
                FortranArray2DRef<double> wk2d(&wk(1), npl * (jdim - 1) * idim, 35);
                FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim - 1);
                diagi_ns::diagi(k, npl, jdim, kdim, idim, q4d, res4d, dtj3d, si4d,
                                wk2d, vol3d, vist3d, blank, iover);
            }
        }

        //
        // l2-norm of delta q
        //
        if (icyc <= 2 || icyc == ncyc) {
            // l2norm declares res(jdim,kdim,idim-1,5), vol(jdim,kdim,idim-1)
            FortranArray4DRef<double> res4d(&res(1, 1), jdim, kdim, idim - 1, 5);
            FortranArray3DRef<double> vol3d(&vol(1), jdim, kdim, idim - 1);
            int irdq1 = 1;
            l2norm_ns::l2norm(nbl, icyc, resid, irdq1, jdim, kdim, idim, res4d, vol3d);
        }
    }

    //
    // ****************** spatially-split af ***********************
    //
    return;
} // end af3f

} // namespace af3f_ns
