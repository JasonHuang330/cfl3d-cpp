// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "hfluxr.h"
#include "xlim.h"
#include "fluxp.h"
#include "fluxm.h"
#include "fhat.h"
#include "fmaps.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace hfluxr_ns {

void hfluxr(int& i, int& npl, double& xkap, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qk0, FortranArray3DRef<double> sk,
            FortranArray2DRef<double> t, int& nvtq, int& nv, int& nfa,
            FortranArray1DRef<double> wfa, FortranArray1DRef<int> iwfa,
            FortranArray1DRef<int> kbctyp, int& isf, int& nbl,
            FortranArray3DRef<double> bck, FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
            int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& idef)
{
    // COMMON block aliases
    float& gamma     = cmn_fluid.gamma;
    int32_t& ichk    = cmn_chk.ichk;
    int32_t* ncell   = cmn_cpurate.ncell;   // 1-based: ncell[level-1]
    int32_t& level   = cmn_mgrd.level;
    int32_t& i2d     = cmn_twod.i2d;
    int32_t* ifds    = cmn_fvfds.ifds;      // 1-based: ifds[3-1]
    int32_t* iflim   = cmn_info.iflim;      // 1-based: iflim[3-1]
    int32_t& isklton = cmn_sklton.isklton;

    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    int n   = npl * jdim * kdim;
    int nr  = n - jdim;
    int nr2 = n - 2 * jdim;

    if ((float)xkap < -2.0f) {
        //
        // first order
        //
        // fill temp arrays with interior values
        //
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= nr; izz++) {
                t(izz + jdim, 20 + l) = q(izz, 1, i, l);
                t(izz,        25 + l) = q(izz, 1, i, l);
            }
        }
        //
        // fill temp arrays with boundary values
        //
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii  = i + ipl - 1;
            int jk  = 1 + (ipl - 1) * jdim * kdim;
            int jk2 = 1 + (ipl - 1) * jdim * kdim + kdim1 * jdim;
            for (int l = 1; l <= 5; l++) {
                for (int izz = 1; izz <= jdim; izz++) {
                    t(izz + jk  - 1, 20 + l) = qk0(izz, ii, l, 1);
                    t(izz + jk2 - 1, 25 + l) = qk0(izz, ii, l, 3);
                }
            }
        }
    } else {
        //
        // higher order
        //
        for (int l = 1; l <= 5; l++) {
            //
            // calculate interior gradients
            //
            for (int izz = 1; izz <= nr2; izz++) {
                t(izz + jdim, 1) = q(izz + jdim, 1, i, l) - q(izz, 1, i, l);
            }
            //
            // edge gradients - left boundary
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                int jc = (ipl - 1) * kdim * jdim;
                for (int jj = 1; jj <= jdim1; jj++) {
                    jc = jc + 1;
                    t(jc, 1) = (1.0 - bck(jj, ii, 1)) * (q(jj, 1, ii, l) - qk0(jj, ii, l, 1))
                                   + bck(jj, ii, 1) * qk0(jj, ii, l, 2);
                }
            }
            //
            // edge gradients - right boundary
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                int jc = (ipl - 1) * kdim * jdim + kdim1 * jdim;
                for (int jj = 1; jj <= jdim1; jj++) {
                    jc = jc + 1;
                    t(jc, 1) = (1.0 - bck(jj, ii, 2)) * (qk0(jj, ii, l, 3) - q(jj, kdim1, ii, l))
                                   + bck(jj, ii, 2) * qk0(jj, ii, l, 4);
                }
            }
            //
            // zero out gradients for dummy interfaces at j=jdim
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int jvl = (ipl - 1) * jdim * kdim;
                for (int k = 1; k <= kdim; k++) {
                    t(k * jdim + jvl, 1) = 0.0;
                }
            }
            //
            // gradient limiting - cell interface interpolations
            //
            for (int izz = 1; izz <= nr; izz++) {
                t(izz, 2) = t(izz + jdim, 1);
            }
            {
                int ifl = iflim[3 - 1];
                int ncells;
                if (ifl == 4) {
                    if (i2d == 1) {
                        ncells = (int)std::sqrt((float)ncell[level - 1]);
                    } else {
                        ncells = (int)std::pow((float)ncell[level - 1], 1.0f / 3.0f);
                    }
                } else {
                    ncells = kdim1;
                }
                FortranArray1DRef<double> t_col1(&t(1, 1), t.size(1), t.lbound(1));
                FortranArray1DRef<double> t_col2(&t(1, 2), t.size(1), t.lbound(1));
                FortranArray1DRef<double> q_col(&q(1, 1, i, l), n, 1);
                xlim_ns::xlim(xkap, nr, t_col1, t_col2, q_col, ifl, ncells, l);
            }
            for (int izz = 1; izz <= nr; izz++) {
                t(izz + jdim, 20 + l) = q(izz, 1, i, l) + t(izz, 2);
                t(izz,        25 + l) = q(izz, 1, i, l) - t(izz, 1);
            }
        } // end l loop (do 1000)
    } // end if xkap


    // do 2000 l=1,5
    for (int l = 1; l <= 5; l++) {
        //
        // edge values - left boundary
        //
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii = i + ipl - 1;
            int jl = 0;
            for (int jj = 1; jj <= jdim1; jj++) {
                jl = jl + 1;
                t(jl, 1) = qk0(jj, ii, l, 1) - qk0(jj, ii, l, 2);
                t(jl, 2) = q(jj, 1, ii, l)   - qk0(jj, ii, l, 1);
                t(jl, 3) = qk0(jj, ii, l, 1);
            }
            {
                int ifl = iflim[3 - 1];
                int ncells;
                if (ifl == 4) {
                    if (i2d == 1) {
                        ncells = (int)std::sqrt((float)ncell[level - 1]);
                    } else {
                        ncells = (int)std::pow((float)ncell[level - 1], 1.0f / 3.0f);
                    }
                } else {
                    ncells = kdim1;
                }
                FortranArray1DRef<double> t_col1(&t(1, 1), t.size(1), t.lbound(1));
                FortranArray1DRef<double> t_col2(&t(1, 2), t.size(1), t.lbound(1));
                FortranArray1DRef<double> t_col3(&t(1, 3), t.size(1), t.lbound(1));
                xlim_ns::xlim(xkap, jdim1, t_col1, t_col2, t_col3, ifl, ncells, l);
            }
            int jc = (ipl - 1) * (kdim * jdim);
            jl = 0;
            for (int jj = 1; jj <= jdim1; jj++) {
                jl = jl + 1;
                jc = jc + 1;
                t(jc, 20 + l) = (1.0 - bck(jj, ii, 1)) * (qk0(jj, ii, l, 1) + t(jl, 2))
                                    + bck(jj, ii, 1) * qk0(jj, ii, l, 1);
                t(jc, 25 + l) = (1.0 - bck(jj, ii, 1)) * t(jc, 25 + l)
                                    + bck(jj, ii, 1) * qk0(jj, ii, l, 1);
            }
        } // end ipl (do 1300)
        //
        // edge values - right boundary
        //
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii = i + ipl - 1;
            int jl = 0;
            for (int jj = 1; jj <= jdim1; jj++) {
                jl = jl + 1;
                t(jl, 2) = qk0(jj, ii, l, 4) - qk0(jj, ii, l, 3);
                t(jl, 1) = qk0(jj, ii, l, 3) - q(jj, kdim1, ii, l);
                t(jl, 3) = qk0(jj, ii, l, 3);
            }
            {
                int ifl = iflim[3 - 1];
                int ncells;
                if (ifl == 4) {
                    if (i2d == 1) {
                        ncells = (int)std::sqrt((float)ncell[level - 1]);
                    } else {
                        ncells = (int)std::pow((float)ncell[level - 1], 1.0f / 3.0f);
                    }
                } else {
                    ncells = kdim1;
                }
                FortranArray1DRef<double> t_col1(&t(1, 1), t.size(1), t.lbound(1));
                FortranArray1DRef<double> t_col2(&t(1, 2), t.size(1), t.lbound(1));
                FortranArray1DRef<double> t_col3(&t(1, 3), t.size(1), t.lbound(1));
                xlim_ns::xlim(xkap, jdim1, t_col1, t_col2, t_col3, ifl, ncells, l);
            }
            int jc = (ipl - 1) * (kdim * jdim) + kdim1 * jdim;
            jl = 0;
            for (int jj = 1; jj <= jdim1; jj++) {
                jl = jl + 1;
                jc = jc + 1;
                t(jc, 20 + l) = (1.0 - bck(jj, ii, 2)) * t(jc, 20 + l)
                                    + bck(jj, ii, 2) * qk0(jj, ii, l, 3);
                t(jc, 25 + l) = (1.0 - bck(jj, ii, 2)) * (qk0(jj, ii, l, 3) - t(jl, 1))
                                    + bck(jj, ii, 2) * qk0(jj, ii, l, 3);
            }
        } // end ipl (do 1800)
        //
        // fill end points for safety
        t(jdim,              20 + l) = t(jdim - 1,              20 + l);
        t(jdim * kdim * npl, 25 + l) = t(jdim * kdim * npl - 1, 25 + l);
    } // end l=1,5 (do 2000)




    // ichk check
    if (ichk == 1) {
        double epsz = 1.0e-03;
        double epss = 1.0e+03;
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii = i + ipl - 1;
            for (int k = 1; k <= kdim; k++) {
                for (int j = 1; j <= jdim; j++) {
                    int ic = (ipl - 1) * jdim * kdim + (k - 1) * jdim + j;
                    if ((float)t(ic, 21) < (float)epsz ||
                        (float)t(ic, 25) < (float)epsz ||
                        (float)t(ic, 21) > (float)epss ||
                        (float)t(ic, 25) > (float)epss) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " on block %d", nbl);
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " stopping in hflux left - small (or large) density and/or pressure at ");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " j,k,i,t(21),t(25) = %d %d %d %g %g",
                            j, k, ii, (float)t(ic, 21), (float)t(ic, 25));
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }
                    if ((float)t(ic, 26) < (float)epsz ||
                        (float)t(ic, 30) < (float)epsz ||
                        (float)t(ic, 26) > (float)epss ||
                        (float)t(ic, 30) > (float)epss) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " on block %d", nbl);
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " stopping in hflux right - small  (or large) density and/or pressure at ");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " j,k,i,t(26),t(30) = %d %d %d %g %g",
                            j, k, ii, (float)t(ic, 26), (float)t(ic, 30));
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
    }

    // ifds(3) branch: flux splitting method
    if (ifds[3 - 1] == 0) {
        //
        if (isklton > 0 && i == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   computing inviscid fluxes, K-direction - flux-vector splitting - block%4d", nbl);
        }
        //
        {
            FortranArray1DRef<double> sk_i1(&sk(1, i, 1), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i2(&sk(1, i, 2), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i3(&sk(1, i, 3), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i4(&sk(1, i, 4), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i5(&sk(1, i, 5), sk.size(1), sk.lbound(1));
            FortranArray2DRef<double> t_21(&t(1, 21), nvtq, 5);
            FortranArray2DRef<double> t_31(&t(1, 31), nvtq, 5);
            fluxp_ns::fluxp(sk_i1, sk_i2, sk_i3, sk_i4, sk_i5,
                            t_21, t_31, n, t, n, nvtq, nou, bou, nbuf, ibufdim);
        }
        //
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 20 + l) = t(izz, 30 + l);
            }
        }
        //
        {
            FortranArray1DRef<double> sk_i1(&sk(1, i, 1), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i2(&sk(1, i, 2), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i3(&sk(1, i, 3), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i4(&sk(1, i, 4), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i5(&sk(1, i, 5), sk.size(1), sk.lbound(1));
            FortranArray2DRef<double> t_26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_31(&t(1, 31), nvtq, 5);
            fluxm_ns::fluxm(sk_i1, sk_i2, sk_i3, sk_i4, sk_i5,
                            t_26, t_31, n, t, n, nvtq, nou, bou, nbuf, ibufdim);
        }
        //
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= n; izz++) {
                t(izz, 30 + l) = t(izz, 30 + l) + t(izz, 20 + l);
            }
        }

    } else if (ifds[3 - 1] == 1) {
        //
        if (isklton > 0 && i == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   computing inviscid fluxes, K-direction - flux-difference splitting");
        }
        //
        {
            FortranArray1DRef<double> sk_i1(&sk(1, i, 1), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i2(&sk(1, i, 2), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i3(&sk(1, i, 3), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i4(&sk(1, i, 4), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i5(&sk(1, i, 5), sk.size(1), sk.lbound(1));
            FortranArray2DRef<double> t_31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_21(&t(1, 21), nvtq, 5);
            fhat_ns::fhat(sk_i1, sk_i2, sk_i3, sk_i4, sk_i5,
                          t_31, t_26, t_21, n, nvtq);
        }

    } else {
        //
        if (isklton > 0 && i == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            std::snprintf(bou(nou(1), 1), 120,
                "   computing inviscid fluxes, K-direction - MAPS+ flux splitting");
        }
        //
        {
            FortranArray1DRef<double> sk_i1(&sk(1, i, 1), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i2(&sk(1, i, 2), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i3(&sk(1, i, 3), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i4(&sk(1, i, 4), sk.size(1), sk.lbound(1));
            FortranArray1DRef<double> sk_i5(&sk(1, i, 5), sk.size(1), sk.lbound(1));
            FortranArray2DRef<double> t_31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_21(&t(1, 21), nvtq, 5);
            fmaps_ns::fmaps(sk_i1, sk_i2, sk_i3, sk_i4, sk_i5,
                            t_31, t_26, t_21, n, nvtq);
        }
    }

    //
    // conservative at embedded boundaries - wfa array contains fluxes
    // from a finer grid
    //
    if (nfa > 0) {
        int iis  = i;
        int iie  = i + npl - 1;
        int lfcc = iwfa(7);
        for (int ifa = 1; ifa <= nfa; ifa++) {
            int ic   = (ifa - 1) * 7;
            int jsb  = iwfa(ic + 1);
            int ksb  = iwfa(ic + 2);
            int isb  = iwfa(ic + 3);
            int jeb  = iwfa(ic + 4);
            int keb  = iwfa(ic + 5);
            int ieb  = iwfa(ic + 6);
            int ifts = iwfa(ic + 7);

            if (ifts != lfcc) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " inconsistent summations - stopping");
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }

            if (isklton > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "  installing accumulated fluxes on k=%3d at js,je,is,ie=%4d%4d%4d%4d",
                    ksb, jsb, jeb, isb, ieb);
            }

            //
            // loop over all planes in I-direction
            //
            lfcc = ifts;
            for (int ii = 1; ii <= idim1; ii++) {
                //
                // skip planes not in embedded region
                //
                if (ii >= isb && ii < ieb) {
                    //
                    // check for ii in region to be updated on this pass
                    //
                    if (iis <= ii && ii <= iie) {
                        int ipl = ii - i + 1;
                        for (int l = 1; l <= 5; l++) {
                            int lfc = ifts + (l - 1) * (jeb - jsb) * (ieb - isb)
                                           + (ii - isb) * (jeb - jsb);
                            int loc = (ipl - 1) * jdim * kdim + (ksb - 1) * jdim + jsb;
                            for (int j = jsb; j <= jeb - 1; j++) {
                                t(loc, 30 + l) = wfa(lfc);
                                loc  = loc  + 1;
                                lfc  = lfc  + 1;
                                lfcc = lfcc + 1;
                            }
                        }
                    } else {
                        //
                        // increment counter for planes not updated on this pass
                        //
                        lfcc = lfcc + 5 * (jeb - jsb);
                    }
                }
            } // end ii loop (do 754)
        } // end ifa loop (do 7055)
    }

    //
    for (int l = 1; l <= 5; l++) {
        for (int izz = 1; izz <= nr; izz++) {
            res(izz, 1, i, l) = res(izz, 1, i, l) + t(izz + jdim, 30 + l) - t(izz, 30 + l);
        }
    }

    //
    // geometric conservation law terms for deforming grids
    //
    if (idef > 0) {
        double oogmo = 1.0 / ((double)gamma - 1.0);
        for (int izz = 1; izz <= nr; izz++) {
            t(izz, 41) = -(sk(izz + jdim, i, 5) * sk(izz + jdim, i, 4)
                         - sk(izz,         i, 5) * sk(izz,         i, 4));
        }
        for (int ii = i; ii <= i + npl - 1; ii++) {
            int ijk2 = (ii - i + 1) * jdim * kdim - jdim;
            for (int izz = 1; izz <= jdim; izz++) {
                t(izz + ijk2, 41) = -(sk(izz + ijk2,       i, 5) * sk(izz + ijk2,       i, 4)
                                    - sk(izz + ijk2 - jdim, i, 5) * sk(izz + ijk2 - jdim, i, 4));
            }
        }
        for (int izz = 1; izz <= nr; izz++) {
            t(izz, 36) = q(izz, 1, i, 1);
            t(izz, 37) = q(izz, 1, i, 1) * q(izz, 1, i, 2);
            t(izz, 38) = q(izz, 1, i, 1) * q(izz, 1, i, 3);
            t(izz, 39) = q(izz, 1, i, 1) * q(izz, 1, i, 4);
            t(izz, 40) = q(izz, 1, i, 5) * oogmo
                       + 0.5 * q(izz, 1, i, 1) * (q(izz, 1, i, 2) * q(izz, 1, i, 2)
                                                  + q(izz, 1, i, 3) * q(izz, 1, i, 3)
                                                  + q(izz, 1, i, 4) * q(izz, 1, i, 4));
        }
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= nr; izz++) {
                res(izz, 1, i, l) = res(izz, 1, i, l) + t(izz, 35 + l) * t(izz, 41);
            }
        }
    }

    //
    // store finer-grid fluxes for enforcing conservation on coarser meshes
    //
    if (isf == 1) {
        if (kbctyp(1) == 21) {
            //
            if (isklton > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "  installing k fluxes into qk0 for edge k=0 for block%3d", nbl);
            }
            //
            // left boundary
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                for (int l = 1; l <= 5; l++) {
                    int jk = (ipl - 1) * jdim * kdim;
                    for (int j = 1; j <= jdim1; j++) {
                        jk = jk + 1;
                        qk0(j, ii, l, 2) = t(jk, 30 + l);
                    }
                }
            }
        }
        if (kbctyp(2) == 21) {
            //
            if (isklton > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120,
                    "  installing k fluxes into qk0 for edge k=kdim for block%3d", nbl);
            }
            //
            // right boundary
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                for (int l = 1; l <= 5; l++) {
                    int jk = (ipl - 1) * jdim * kdim + jdim * kdim1;
                    for (int j = 1; j <= jdim1; j++) {
                        jk = jk + 1;
                        qk0(j, ii, l, 4) = t(jk, 30 + l);
                    }
                }
            }
        }
    }

    return;
}

} // namespace hfluxr_ns
