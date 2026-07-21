// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "ffluxr.h"
#include "runtime/fortran_io.h"
#include "xlim.h"
#include "termn8.h"
#include "fluxp.h"
#include "fluxm.h"
#include "fhat.h"
#include "fmaps.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace ffluxr_ns {

void ffluxr(int& k, int& npl, double& xkap, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qi0, FortranArray3DRef<double> si,
            FortranArray2DRef<double> t, int& nvtq, int& nv, int& nfa,
            FortranArray1DRef<double> wfa, FortranArray1DRef<int> iwfa,
            FortranArray1DRef<int> ibctyp, int& isf, int& nbl,
            FortranArray3DRef<double> bci, FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
            int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& idef)
{
    // COMMON block aliases
    float& gamma   = cmn_fluid.gamma;
    int32_t& level = cmn_mgrd.level;
    int32_t& i2d   = cmn_twod.i2d;
    int32_t& ichk  = cmn_chk.ichk;
    // ifds(1) -> cmn_fvfds.ifds[0]  (0-based)
    // iflim(1) -> cmn_info.iflim[0] (0-based)
    // ncell(level) -> cmn_cpurate.ncell[level-1] (0-based)
    int32_t& isklton = cmn_sklton.isklton;

    // Local variables
    int jdim1, idim1, kdim1;
    int jv, n, nr;
    int i, l, izz, kpl, kk, jc, jl, jk, jk1, kv;
    int ifl, ncells;
    double epsz, epss, oogmo;
    int kc;
    // embedded boundary vars
    int kks, kke, lfcc, ifa, ic, jsb, ksb, isb, jeb, keb, ieb, ifts;
    int lfc, loc;
    // termn8 needs int& for literal -1
    int minus1 = -1;



    jdim1 = jdim - 1;
    idim1 = idim - 1;
    kdim1 = kdim - 1;

    // Mirror Fortran dummy-arg contract: res(jdim,kdim,idim-1,5). res is indexed
    // res(izz,k,i,l) with SEPARATE izz and k (k = chunk-start plane), so the k
    // stride must be jdim, not jdim*kdim. A caller that folded (j,k) into a single
    // dim (jdim*kdim,1,...) only works when k==1 (npl==kdim); for the chunked path
    // (npl!=kdim, k>1) it lands accumulations on the wrong k-plane. Re-view here so
    // the callee owns the correct shape regardless of how the caller dimensioned it.
    res = FortranArray4DRef<double>(res.data(), jdim, kdim, idim1, 5);

    if (npl == kdim1 && nvtq >= jdim * kdim * idim) npl = kdim;

    jv = npl * jdim;
    n  = jv * idim;
    nr = n - jv;

    if (npl != kdim) {
        // fill temp arrays with metrics if all k planes are not done at once
        for (i = 1; i <= idim; i++) {
            int jk_  = (i - 1) * jv + 1;
            int jk1_ = (k - 1) * jdim + 1;
            for (l = 1; l <= 5; l++) {
                for (izz = 1; izz <= jv; izz++) {
                    t(izz + jk_ - 1, 14 + l) = si(izz + jk1_ - 1, i, l);
                }
            }
        }
    }

    // all orders
    // fill temp arrays with interior values
    for (l = 1; l <= 5; l++) {
        if (npl == kdim) {
            for (izz = 1; izz <= nr; izz++) {
                t(izz + jv, 20 + l) = q(izz, 1, 1, l);
                t(izz,      25 + l) = q(izz, 1, 1, l);
            }
        } else {
            for (i = 1; i <= idim1; i++) {
                jk = (i - 1) * jv + 1;
                for (izz = 1; izz <= jv; izz++) {
                    t(izz + jv + jk - 1, 20 + l) = q(izz, k, i, l);
                    t(izz + jk - 1,      25 + l) = q(izz, k, i, l);
                }
            }
        }

        // fill temp array with boundary values
        for (izz = 1; izz <= jv; izz++) {
            t(izz,      20 + l) = qi0(izz, k, l, 1);
            t(izz + nr, 25 + l) = qi0(izz, k, l, 3);
        }
    }

    if (xkap > -2.0) {
        // higher order
        for (l = 1; l <= 5; l++) {
            // gradients on interior planes
            for (izz = 1; izz <= nr - jv; izz++) {
                t(izz + jv, 1) = t(izz + jv, 25 + l) - t(izz, 25 + l);
            }

            // initialize gradients to zero on planes i=1 and i=idim
            kv = jv * idim1;
            for (izz = 1; izz <= jv; izz++) {
                t(izz,      1) = 0.0;
                t(kv + izz, 1) = 0.0;
            }

            // edge gradients - left boundary
            for (kpl = 1; kpl <= npl; kpl++) {
                kk = k + kpl - 1;
                jc = (kpl - 1) * jdim;
                for (int jj = 1; jj <= jdim1; jj++) {
                    jc = jc + 1;
                    t(jc, 1) = (1.0 - bci(jj, kk, 1)) * (q(jj, kk, 1, l) - qi0(jj, kk, l, 1))
                                   + bci(jj, kk, 1) * qi0(jj, kk, l, 2);
                }
            }

            // edge gradients - right boundary
            for (kpl = 1; kpl <= npl; kpl++) {
                kk = k + kpl - 1;
                jc = idim1 * jv + (kpl - 1) * jdim;
                for (int jj = 1; jj <= jdim1; jj++) {
                    jc = jc + 1;
                    t(jc, 1) = (1.0 - bci(jj, kk, 2)) * (qi0(jj, kk, l, 3) - q(jj, kk, idim1, l))
                                   + bci(jj, kk, 2) * qi0(jj, kk, l, 4);
                }
            }

            // gradient limiting - cell interface interpolations
            for (izz = 1; izz <= nr; izz++) {
                t(izz, 2) = t(izz + jv, 1);
            }

            ifl = cmn_info.iflim[0];
            if (ifl == 4) {
                if (i2d == 1) {
                    ncells = (int)std::sqrt((float)cmn_cpurate.ncell[level - 1]);
                } else {
                    ncells = (int)std::pow((float)cmn_cpurate.ncell[level - 1], 1.0f / 3.0f);
                }
            } else {
                ncells = idim1;
            }
            {
                FortranArray1DRef<double> t_col1(&t(1, 1), nr);
                FortranArray1DRef<double> t_col2(&t(1, 2), nr);
                FortranArray1DRef<double> t_col25l(&t(1, 25 + l), nr);
                xlim_ns::xlim(xkap, nr, t_col1, t_col2, t_col25l, ifl, ncells, l);
            }

            for (izz = 1; izz <= nr; izz++) {
                t(izz + jv, 20 + l) = t(izz + jv, 20 + l) + t(izz, 2);
                t(izz,      25 + l) = t(izz,      25 + l) - t(izz, 1);
            }
        }
    }



    // loop 2000: edge values
    for (l = 1; l <= 5; l++) {
        // edge values - left boundary
        for (kpl = 1; kpl <= npl; kpl++) {
            kk = k + kpl - 1;
            jl = 0;
            for (int jj = 1; jj <= jdim1; jj++) {
                jl = jl + 1;
                t(jl, 1) = qi0(jj, kk, l, 1) - qi0(jj, kk, l, 2);
                t(jl, 2) = q(jj, kk, 1, l)   - qi0(jj, kk, l, 1);
                t(jl, 3) = qi0(jj, kk, l, 1);
            }
            ifl = cmn_info.iflim[0];
            if (ifl == 4) {
                if (i2d == 1) {
                    ncells = (int)std::sqrt((float)cmn_cpurate.ncell[level - 1]);
                } else {
                    ncells = (int)std::pow((float)cmn_cpurate.ncell[level - 1], 1.0f / 3.0f);
                }
            } else {
                ncells = idim1;
            }
            {
                FortranArray1DRef<double> t_col1(&t(1, 1), jdim1);
                FortranArray1DRef<double> t_col2(&t(1, 2), jdim1);
                FortranArray1DRef<double> t_col3(&t(1, 3), jdim1);
                xlim_ns::xlim(xkap, jdim1, t_col1, t_col2, t_col3, ifl, ncells, l);
            }
            jc = (kpl - 1) * jdim;
            jl = 0;
            for (int jj = 1; jj <= jdim1; jj++) {
                jl = jl + 1;
                jc = jc + 1;
                t(jc, 20 + l) = (1.0 - bci(jj, kk, 1)) * (qi0(jj, kk, l, 1) + t(jl, 2))
                                    + bci(jj, kk, 1) * qi0(jj, kk, l, 1);
                t(jc, 25 + l) = (1.0 - bci(jj, kk, 1)) * t(jc, 25 + l)
                                    + bci(jj, kk, 1) * qi0(jj, kk, l, 1);
            }
        }

        // edge values - right boundary
        for (kpl = 1; kpl <= npl; kpl++) {
            kk = k + kpl - 1;
            jl = 0;
            for (int jj = 1; jj <= jdim1; jj++) {
                jl = jl + 1;
                t(jl, 2) = qi0(jj, kk, l, 4) - qi0(jj, kk, l, 3);
                t(jl, 1) = qi0(jj, kk, l, 3) - q(jj, kk, idim1, l);
                t(jl, 3) = qi0(jj, kk, l, 3);
            }
            ifl = cmn_info.iflim[0];
            if (ifl == 4) {
                if (i2d == 1) {
                    ncells = (int)std::sqrt((float)cmn_cpurate.ncell[level - 1]);
                } else {
                    ncells = (int)std::pow((float)cmn_cpurate.ncell[level - 1], 1.0f / 3.0f);
                }
            } else {
                ncells = idim1;
            }
            {
                FortranArray1DRef<double> t_col1(&t(1, 1), jdim1);
                FortranArray1DRef<double> t_col2(&t(1, 2), jdim1);
                FortranArray1DRef<double> t_col3(&t(1, 3), jdim1);
                xlim_ns::xlim(xkap, jdim1, t_col1, t_col2, t_col3, ifl, ncells, l);
            }
            jc = idim1 * jv + (kpl - 1) * jdim;
            jl = 0;
            for (int jj = 1; jj <= jdim1; jj++) {
                jl = jl + 1;
                jc = jc + 1;
                t(jc, 20 + l) = (1.0 - bci(jj, kk, 2)) * t(jc, 20 + l)
                                    + bci(jj, kk, 2) * qi0(jj, kk, l, 3);
                t(jc, 25 + l) = (1.0 - bci(jj, kk, 2)) * (qi0(jj, kk, l, 3) - t(jl, 1))
                                    + bci(jj, kk, 2) * qi0(jj, kk, l, 3);
            }
        }

        // fill end points for safety
        for (kpl = 1; kpl <= npl; kpl++) {
            jc = (kpl - 1) * jdim + jdim;
            int jl2 = idim1 * jv + (kpl - 1) * jdim + jdim;
            t(jc,  20 + l) = t(jc  - 1, 20 + l);
            t(jl2, 25 + l) = t(jl2 - 1, 25 + l);
        }
    }



    // ichk check
    if (ichk == 1) {
        epsz = 1.0e-03;
        epss = 1.0e+03;
        for (int ipl = 1; ipl <= npl; ipl++) {
            kk = k + ipl - 1;
            for (i = 1; i <= idim; i++) {
                for (int j = 1; j <= jdim; j++) {
                    kc = jv * (i - 1) + (ipl - 1) * jdim + j;
                    if (t(kc, 21) < epsz || t(kc, 25) < epsz ||
                        t(kc, 21) > epss || t(kc, 25) > epss) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " on block %d", nbl);
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " stopping in fflux left - small (or large) density and/or pressure at ");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " j,k,i,t(21),t(25) = %d %d %d %g %g",
                            j, kk, i, (float)t(kc, 21), (float)t(kc, 25));
                        termn8_ns::termn8(myid, minus1, ibufdim, nbuf, bou, nou);
                    }

                    if (t(kc, 26) < epsz || t(kc, 30) < epsz ||
                        t(kc, 26) > epss || t(kc, 30) > epss) {
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120, " on block %d", nbl);
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " stopping in fflux right - small (or large) density and/or pressure at ");
                        nou(1) = std::min(nou(1) + 1, ibufdim);
                        std::snprintf(bou(nou(1), 1), 120,
                            " j,k,i,t(26),t(30) = %d %d %d %g %g",
                            j, kk, i, (float)t(kc, 26), (float)t(kc, 30));
                        termn8_ns::termn8(myid, minus1, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
    }



    if (cmn_fvfds.ifds[0] == 0) {
        // flux-vector splitting
        if (isklton > 0 && k == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            // format 185: "   computing inviscid fluxes, I-direction - flux-vector splitting - block  N"
            std::snprintf(bou(nou(1), 1), 120,
                "   computing inviscid fluxes, I-direction - flux-vector splitting - block%4d", nbl);
        }

        if (npl == kdim) {
            FortranArray1DRef<double> si_col1(&si(1, 1, 1), n);
            FortranArray1DRef<double> si_col2(&si(1, 1, 2), n);
            FortranArray1DRef<double> si_col3(&si(1, 1, 3), n);
            FortranArray1DRef<double> si_col4(&si(1, 1, 4), n);
            FortranArray1DRef<double> si_col5(&si(1, 1, 5), n);
            FortranArray2DRef<double> t_col21(&t(1, 21), nvtq, 5);
            FortranArray2DRef<double> t_col31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_full(&t(1, 1), nvtq, nv);
            fluxp_ns::fluxp(si_col1, si_col2, si_col3, si_col4, si_col5,
                            t_col21, t_col31, n, t_full, n, nvtq,
                            nou, bou, nbuf, ibufdim);
        } else {
            FortranArray1DRef<double> t_col15(&t(1, 15), n);
            FortranArray1DRef<double> t_col16(&t(1, 16), n);
            FortranArray1DRef<double> t_col17(&t(1, 17), n);
            FortranArray1DRef<double> t_col18(&t(1, 18), n);
            FortranArray1DRef<double> t_col19(&t(1, 19), n);
            FortranArray2DRef<double> t_col21(&t(1, 21), nvtq, 5);
            FortranArray2DRef<double> t_col31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_full(&t(1, 1), nvtq, nv);
            fluxp_ns::fluxp(t_col15, t_col16, t_col17, t_col18, t_col19,
                            t_col21, t_col31, n, t_full, n, nvtq,
                            nou, bou, nbuf, ibufdim);
        }

        for (l = 1; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                t(izz, 20 + l) = t(izz, 30 + l);
            }
        }

        if (npl == kdim) {
            FortranArray1DRef<double> si_col1(&si(1, 1, 1), n);
            FortranArray1DRef<double> si_col2(&si(1, 1, 2), n);
            FortranArray1DRef<double> si_col3(&si(1, 1, 3), n);
            FortranArray1DRef<double> si_col4(&si(1, 1, 4), n);
            FortranArray1DRef<double> si_col5(&si(1, 1, 5), n);
            FortranArray2DRef<double> t_col26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_col31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_full(&t(1, 1), nvtq, nv);
            fluxm_ns::fluxm(si_col1, si_col2, si_col3, si_col4, si_col5,
                            t_col26, t_col31, n, t_full, n, nvtq,
                            nou, bou, nbuf, ibufdim);
        } else {
            FortranArray1DRef<double> t_col15(&t(1, 15), n);
            FortranArray1DRef<double> t_col16(&t(1, 16), n);
            FortranArray1DRef<double> t_col17(&t(1, 17), n);
            FortranArray1DRef<double> t_col18(&t(1, 18), n);
            FortranArray1DRef<double> t_col19(&t(1, 19), n);
            FortranArray2DRef<double> t_col26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_col31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_full(&t(1, 1), nvtq, nv);
            fluxm_ns::fluxm(t_col15, t_col16, t_col17, t_col18, t_col19,
                            t_col26, t_col31, n, t_full, n, nvtq,
                            nou, bou, nbuf, ibufdim);
        }

        for (l = 1; l <= 5; l++) {
            for (izz = 1; izz <= n; izz++) {
                t(izz, 30 + l) = t(izz, 30 + l) + t(izz, 20 + l);
            }
        }

    } else if (cmn_fvfds.ifds[0] == 1) {
        // flux-difference splitting
        if (isklton > 0 && k == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            // format 184: "   computing inviscid fluxes, I-direction - flux-difference splitting"
            std::snprintf(bou(nou(1), 1), 120,
                "   computing inviscid fluxes, I-direction - flux-difference splitting");
        }
        if (npl == kdim) {
            FortranArray1DRef<double> si_col1(&si(1, 1, 1), n);
            FortranArray1DRef<double> si_col2(&si(1, 1, 2), n);
            FortranArray1DRef<double> si_col3(&si(1, 1, 3), n);
            FortranArray1DRef<double> si_col4(&si(1, 1, 4), n);
            FortranArray1DRef<double> si_col5(&si(1, 1, 5), n);
            FortranArray2DRef<double> t_col31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_col26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_col21(&t(1, 21), nvtq, 5);
            fhat_ns::fhat(si_col1, si_col2, si_col3, si_col4, si_col5,
                          t_col31, t_col26, t_col21, n, nvtq);
        } else {
            FortranArray1DRef<double> t_col15(&t(1, 15), n);
            FortranArray1DRef<double> t_col16(&t(1, 16), n);
            FortranArray1DRef<double> t_col17(&t(1, 17), n);
            FortranArray1DRef<double> t_col18(&t(1, 18), n);
            FortranArray1DRef<double> t_col19(&t(1, 19), n);
            FortranArray2DRef<double> t_col31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_col26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_col21(&t(1, 21), nvtq, 5);
            fhat_ns::fhat(t_col15, t_col16, t_col17, t_col18, t_col19,
                          t_col31, t_col26, t_col21, n, nvtq);
        }

    } else {
        // MAPS+ flux splitting
        if (isklton > 0 && k == 1) {
            nou(1) = std::min(nou(1) + 1, ibufdim);
            // format 183: "   computing inviscid fluxes, I-direction - MAPS+ flux splitting"
            std::snprintf(bou(nou(1), 1), 120,
                "   computing inviscid fluxes, I-direction - MAPS+ flux splitting");
        }
        if (npl == kdim) {
            FortranArray1DRef<double> si_col1(&si(1, 1, 1), n);
            FortranArray1DRef<double> si_col2(&si(1, 1, 2), n);
            FortranArray1DRef<double> si_col3(&si(1, 1, 3), n);
            FortranArray1DRef<double> si_col4(&si(1, 1, 4), n);
            FortranArray1DRef<double> si_col5(&si(1, 1, 5), n);
            FortranArray2DRef<double> t_col31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_col26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_col21(&t(1, 21), nvtq, 5);
            fmaps_ns::fmaps(si_col1, si_col2, si_col3, si_col4, si_col5,
                            t_col31, t_col26, t_col21, n, nvtq);
        } else {
            FortranArray1DRef<double> t_col15(&t(1, 15), n);
            FortranArray1DRef<double> t_col16(&t(1, 16), n);
            FortranArray1DRef<double> t_col17(&t(1, 17), n);
            FortranArray1DRef<double> t_col18(&t(1, 18), n);
            FortranArray1DRef<double> t_col19(&t(1, 19), n);
            FortranArray2DRef<double> t_col31(&t(1, 31), nvtq, 5);
            FortranArray2DRef<double> t_col26(&t(1, 26), nvtq, 5);
            FortranArray2DRef<double> t_col21(&t(1, 21), nvtq, 5);
            fmaps_ns::fmaps(t_col15, t_col16, t_col17, t_col18, t_col19,
                            t_col31, t_col26, t_col21, n, nvtq);
        }
    }



    // conservative at embedded boundaries - wfa array contains fluxes from a finer grid
    if (nfa > 0) {
        kks  = k;
        kke  = k + npl - 1;
        lfcc = iwfa(7);
        for (ifa = 1; ifa <= nfa; ifa++) {
            ic   = (ifa - 1) * 7;
            jsb  = iwfa(ic + 1);
            ksb  = iwfa(ic + 2);
            isb  = iwfa(ic + 3);
            jeb  = iwfa(ic + 4);
            keb  = iwfa(ic + 5);
            ieb  = iwfa(ic + 6);
            ifts = iwfa(ic + 7);

            if (ifts != lfcc) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                std::snprintf(bou(nou(1), 1), 120, " inconsistent summations - stopping");
                termn8_ns::termn8(myid, minus1, ibufdim, nbuf, bou, nou);
            }

            if (isklton > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 269: "  installing accumulated fluxes on i=  N at js,je,ks,ke=  N  N  N  N"
                std::snprintf(bou(nou(1), 1), 120,
                    "  installing accumulated fluxes on i=%3d at js,je,ks,ke=%4d%4d%4d%4d",
                    isb, jsb, jeb, ksb, keb);
            }

            // loop over all planes in K-direction
            lfcc = ifts;
            for (kk = 1; kk <= kdim1; kk++) {
                // skip planes not in embedded region
                if (kk >= ksb && kk < keb) {
                    // check for kk in region to be updated on this pass (npl planes of data)
                    if (kks <= kk && kk <= kke) {
                        int kpl2 = kk - k + 1;
                        for (l = 1; l <= 5; l++) {
                            lfc = ifts + (l - 1) * (keb - ksb) * (jeb - jsb)
                                       + (kk - ksb) * (jeb - jsb);
                            loc = (kpl2 - 1) * jdim + (isb - 1) * jdim * npl + jsb;
                            for (int j = jsb; j <= jeb - 1; j++) {
                                t(loc, 30 + l) = wfa(lfc);
                                loc  = loc  + 1;
                                lfc  = lfc  + 1;
                                lfcc = lfcc + 1;
                            }
                        }
                    } else {
                        // increment counter for planes not updated on this pass
                        lfcc = lfcc + 5 * (jeb - jsb);
                    }
                }
            }
        }
    }



    if (nbl == 3) {
        for (int col=20; col<=35; col++){
            long badi=-1;
            for (long a=1;a<=(long)n;a++){ double v=t(a,col); if(std::isnan(v)||std::isinf(v)){badi=a;break;} }
        }
        std::fflush(stderr);
    }

    // accumulate residuals
    if (npl == kdim) {
        for (l = 1; l <= 5; l++) {
            for (izz = 1; izz <= nr; izz++) {
                res(izz, 1, 1, l) = res(izz, 1, 1, l) + t(izz + jv, 30 + l) - t(izz, 30 + l);
            }
        }

        // geometric conservation law terms for deforming grids
        if (idef > 0) {
            oogmo = 1.0 / (gamma - 1.0);
            for (izz = 1; izz <= nr; izz++) {
                t(izz, 36) = q(izz, 1, 1, 1);
                t(izz, 37) = q(izz, 1, 1, 1) * q(izz, 1, 1, 2);
                t(izz, 38) = q(izz, 1, 1, 1) * q(izz, 1, 1, 3);
                t(izz, 39) = q(izz, 1, 1, 1) * q(izz, 1, 1, 4);
                t(izz, 40) = q(izz, 1, 1, 5) * oogmo
                           + 0.5 * q(izz, 1, 1, 1) * (q(izz, 1, 1, 2) * q(izz, 1, 1, 2)
                           + q(izz, 1, 1, 3) * q(izz, 1, 1, 3)
                           + q(izz, 1, 1, 4) * q(izz, 1, 1, 4));
                t(izz, 41) = -(si(izz + jv, 1, 5) * si(izz + jv, 1, 4)
                             - si(izz,      1, 5) * si(izz,      1, 4));
            }
            for (l = 1; l <= 5; l++) {
                for (izz = 1; izz <= nr; izz++) {
                    res(izz, 1, 1, l) = res(izz, 1, 1, l) + t(izz, 35 + l) * t(izz, 41);
                }
            }
        }
    } else {
        for (i = 1; i <= idim1; i++) {
            jk = (i - 1) * jv;
            for (l = 1; l <= 5; l++) {
                for (izz = 1; izz <= jv; izz++) {
                    res(izz, k, i, l) = res(izz, k, i, l)
                                      + t(izz + jk + jv, 30 + l)
                                      - t(izz + jk,      30 + l);
                }
            }
        }

        // geometric conservation law terms for deforming grids
        if (idef > 0) {
            oogmo = 1.0 / (gamma - 1.0);
            int jk1_ = (k - 1) * jdim;
            for (i = 1; i <= idim1; i++) {
                for (izz = 1; izz <= jv; izz++) {
                    int izz1 = izz + jk1_;
                    t(izz, 36) = q(izz, k, i, 1);
                    t(izz, 37) = q(izz, k, i, 1) * q(izz, k, i, 2);
                    t(izz, 38) = q(izz, k, i, 1) * q(izz, k, i, 3);
                    t(izz, 39) = q(izz, k, i, 1) * q(izz, k, i, 4);
                    t(izz, 40) = q(izz, k, i, 5) * oogmo
                               + 0.5 * q(izz, k, i, 1) * (q(izz, k, i, 2) * q(izz, k, i, 2)
                               + q(izz, k, i, 3) * q(izz, k, i, 3)
                               + q(izz, k, i, 4) * q(izz, k, i, 4));
                    t(izz, 41) = -(si(izz1, i + 1, 5) * si(izz1, i + 1, 4)
                                 - si(izz1, i,     5) * si(izz1, i,     4));
                }

                for (l = 1; l <= 5; l++) {
                    for (izz = 1; izz <= jv; izz++) {
                        res(izz, k, i, l) = res(izz, k, i, l)
                                          + t(izz, 35 + l) * t(izz, 41);
                    }
                }
            }
        }
    }



    // store finer-grid fluxes for enforcing conservation on coarser meshes
    if (isf == 1) {
        if (ibctyp(1) == 21) {
            if (isklton > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 379: "  installing i fluxes into qi0 for edge i=0 for block  N"
                std::snprintf(bou(nou(1), 1), 120,
                    "  installing i fluxes into qi0 for edge i=0 for block%3d", nbl);
            }

            // left boundary
            for (kpl = 1; kpl <= npl; kpl++) {
                kk = k + kpl - 1;
                for (l = 1; l <= 5; l++) {
                    jk = (kpl - 1) * jdim;
                    for (int j = 1; j <= jdim1; j++) {
                        jk = jk + 1;
                        qi0(j, kk, l, 2) = t(jk, 30 + l);
                    }
                }
            }
        }
        if (ibctyp(2) == 21) {
            if (isklton > 0) {
                nou(1) = std::min(nou(1) + 1, ibufdim);
                // format 389: "  installing i fluxes into qi0 for edge i=idim for block  N"
                std::snprintf(bou(nou(1), 1), 120,
                    "  installing i fluxes into qi0 for edge i=idim for block%3d", nbl);
            }

            // right boundary
            for (kpl = 1; kpl <= npl; kpl++) {
                kk = k + kpl - 1;
                for (l = 1; l <= 5; l++) {
                    jk = idim1 * jv + (kpl - 1) * jdim;
                    for (int j = 1; j <= jdim1; j++) {
                        jk = jk + 1;
                        qi0(j, kk, l, 4) = t(jk, 30 + l);
                    }
                }
            }
        }
    }

    return;
}

} // namespace ffluxr_ns
