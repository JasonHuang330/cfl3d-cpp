// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "gfluxr.h"
#include "runtime/fortran_io.h"
#include "xlim.h"
#include "fluxp.h"
#include "fluxm.h"
#include "fhat.h"
#include "fmaps.h"
#include "termn8.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace gfluxr_ns {

void gfluxr(int& i, int& npl, double& xkap, int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> res, FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0, FortranArray3DRef<double> sj,
            FortranArray2DRef<double> t, int& nvtq, int& nv, int& nfa,
            FortranArray1DRef<double> wfa, FortranArray1DRef<int> iwfa,
            FortranArray1DRef<int> jbctyp, int& isf, int& nbl,
            FortranArray3DRef<double> bcj, FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim,
            int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& idef)
{
    // COMMON block references
    float& gamma     = cmn_fluid.gamma;
    int32_t& ichk    = cmn_chk.ichk;
    int32_t& isklton = cmn_sklton.isklton;
    int32_t& i2d     = cmn_twod.i2d;
    int32_t& level   = cmn_mgrd.level;

    // Derived dimensions
    int idim1 = idim - 1;
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;

    int kg = kdim * npl;
    int n  = jdim * kdim * npl;
    int nr = n - 1;


    if ((float)xkap < -2.0f) {
        //
        // first order
        //
        // fill temp arrays with interior values
        //
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= nr; izz++) {
                t(izz+1, 20+l) = q(izz, 1, i, l);
                t(izz,   25+l) = q(izz, 1, i, l);
            }
            //
            // fill temp arrays with boundary values
            //
            int jj = 1 - jdim;
            for (int kk = 1; kk <= kg; kk++) {
                jj  = jj + jdim;
                int jj2 = jj + jdim1;
                t(jj,  20+l) = qj0(kk, i, l, 1);
                t(jj2, 25+l) = qj0(kk, i, l, 3);
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
            for (int izz = 1; izz <= nr; izz++) {
                t(izz+1, 1) = q(izz+1, 1, i, l) - q(izz, 1, i, l);
            }
            //
            // edge gradients - left boundary
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                int jc = (ipl-1)*kdim*jdim + 1 - jdim;
                for (int kk = 1; kk <= kdim1; kk++) {
                    jc = jc + jdim;
                    t(jc, 1) = (1.0 - bcj(kk,ii,1)) * (q(1,kk,ii,l) - qj0(kk,ii,l,1))
                                   + bcj(kk,ii,1) * qj0(kk,ii,l,2);
                }
            }
            //
            // edge gradients - right boundary
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                int jc = (ipl-1)*kdim*jdim;
                for (int kk = 1; kk <= kdim1; kk++) {
                    jc = jc + jdim;
                    t(jc, 1) = (1.0 - bcj(kk,ii,2)) * (qj0(kk,ii,l,3) - q(jdim1,kk,ii,l))
                                   + bcj(kk,ii,2) * qj0(kk,ii,l,4);
                }
            }
            //
            // zero out gradients of dummy interfaces at k=kdim
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int kvl = (ipl-1)*jdim*kdim + kdim1*jdim;
                for (int j = 1; j <= jdim; j++) {
                    t(j+kvl, 1) = 0.0;
                }
            }
            //
            // gradient limiting - cell interface interpolations
            //
            for (int izz = 1; izz <= nr; izz++) {
                t(izz, 2) = t(izz+1, 1);
            }

            int ifl = cmn_info.iflim[1];
            int ncells;
            if (ifl == 4) {
                if (i2d == 1) {
                    ncells = (int)std::sqrt((float)cmn_cpurate.ncell[level-1]);
                } else {
                    ncells = (int)std::pow((float)cmn_cpurate.ncell[level-1], 1.0f/3.0f);
                }
            } else {
                ncells = jdim1;
            }
            {
                FortranArray1DRef<double> t_col1(&t(1,1), nvtq, 1);
                FortranArray1DRef<double> t_col2(&t(1,2), nvtq, 1);
                FortranArray1DRef<double> q_il(&q(1,1,i,l), nr, 1);
                xlim_ns::xlim(xkap, nr, t_col1, t_col2, q_il, ifl, ncells, l);
            }

            for (int izz = 1; izz <= nr; izz++) {
                t(izz+1, 20+l) = q(izz, 1, i, l) + t(izz, 2);
                t(izz,   25+l) = q(izz, 1, i, l) - t(izz, 1);
            }
        } // end l=1,5
    } // end if xkap


    // do 2000 l=1,5
    for (int l = 1; l <= 5; l++) {
        //
        // edge values - left boundary
        //
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii = i + ipl - 1;
            int kl = 0;
            for (int kk = 1; kk <= kdim1; kk++) {
                kl = kl + 1;
                t(kl, 1) = qj0(kk,ii,l,1) - qj0(kk,ii,l,2);
                t(kl, 2) = q(1,kk,ii,l)   - qj0(kk,ii,l,1);
                t(kl, 3) = qj0(kk,ii,l,1);
            }
            int ifl = cmn_info.iflim[1];
            int ncells;
            if (ifl == 4) {
                if (i2d == 1) {
                    ncells = (int)std::sqrt((float)cmn_cpurate.ncell[level-1]);
                } else {
                    ncells = (int)std::pow((float)cmn_cpurate.ncell[level-1], 1.0f/3.0f);
                }
            } else {
                ncells = jdim1;
            }
            {
                FortranArray1DRef<double> t_col1(&t(1,1), nvtq, 1);
                FortranArray1DRef<double> t_col2(&t(1,2), nvtq, 1);
                FortranArray1DRef<double> t_col3(&t(1,3), nvtq, 1);
                xlim_ns::xlim(xkap, kdim1, t_col1, t_col2, t_col3, ifl, ncells, l);
            }
            int jc = (ipl-1)*(kdim*jdim) + 1 - jdim;
            kl = 0;
            for (int kk = 1; kk <= kdim1; kk++) {
                kl = kl + 1;
                jc = jc + jdim;
                t(jc, 20+l) = (1.0 - bcj(kk,ii,1)) * (qj0(kk,ii,l,1) + t(kl,2))
                                  + bcj(kk,ii,1) * qj0(kk,ii,l,1);
                t(jc, 25+l) = (1.0 - bcj(kk,ii,1)) * t(jc, 25+l)
                                  + bcj(kk,ii,1) * qj0(kk,ii,l,1);
            }
        }
        //
        // edge values - right boundary
        //
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii = i + ipl - 1;
            int kl = 0;
            for (int kk = 1; kk <= kdim1; kk++) {
                kl = kl + 1;
                t(kl, 2) = qj0(kk,ii,l,4) - qj0(kk,ii,l,3);
                t(kl, 1) = qj0(kk,ii,l,3) - q(jdim1,kk,ii,l);
                t(kl, 3) = qj0(kk,ii,l,3);
            }
            int ifl = cmn_info.iflim[1];
            int ncells;
            if (ifl == 4) {
                if (i2d == 1) {
                    ncells = (int)std::sqrt((float)cmn_cpurate.ncell[level-1]);
                } else {
                    ncells = (int)std::pow((float)cmn_cpurate.ncell[level-1], 1.0f/3.0f);
                }
            } else {
                ncells = jdim1;
            }
            {
                FortranArray1DRef<double> t_col1(&t(1,1), nvtq, 1);
                FortranArray1DRef<double> t_col2(&t(1,2), nvtq, 1);
                FortranArray1DRef<double> t_col3(&t(1,3), nvtq, 1);
                xlim_ns::xlim(xkap, kdim1, t_col1, t_col2, t_col3, ifl, ncells, l);
            }
            int jc = (ipl-1)*(kdim*jdim);
            kl = 0;
            for (int kk = 1; kk <= kdim1; kk++) {
                kl = kl + 1;
                jc = jc + jdim;
                t(jc, 20+l) = (1.0 - bcj(kk,ii,2)) * t(jc, 20+l)
                                  + bcj(kk,ii,2) * qj0(kk,ii,l,3);
                t(jc, 25+l) = (1.0 - bcj(kk,ii,2)) * (qj0(kk,ii,l,3) - t(kl,1))
                                  + bcj(kk,ii,2) * qj0(kk,ii,l,3);
            }
        }
        //
        // fill end point for safety
        t(n, 25+l) = t(nr, 25+l);

    } // end do 2000 l=1,5


    // ichk check
    if (ichk == 1) {
        double epsz = 1.0e-03;
        double epss = 1.0e+03;
        for (int ipl = 1; ipl <= npl; ipl++) {
            int ii = i + ipl - 1;
            for (int k = 1; k <= kdim; k++) {
                for (int j = 1; j <= jdim; j++) {
                    int ic = jdim*kdim*(ipl-1) + (k-1)*jdim + j;
                    if ((float)t(ic,21) < (float)epsz ||
                        (float)t(ic,25) < (float)epsz ||
                        (float)t(ic,21) > (float)epss ||
                        (float)t(ic,25) > (float)epss) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " on block %d", nbl);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " stopping in gflux left - small (or large) density and/or pressure at ");
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " j,k,i,t(21),t(25) = %d %d %d %g %g",
                            j, k, ii, (float)t(ic,21), (float)t(ic,25));
                        int ierrflg = -1;
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }

                    if ((float)t(ic,26) < (float)epsz ||
                        (float)t(ic,30) < (float)epsz ||
                        (float)t(ic,26) > (float)epss ||
                        (float)t(ic,30) > (float)epss) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " on block %d", nbl);
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " stopping in gflux right - small  (or large) density and/or pressure at ");
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120,
                            " j,k,i,t(26),t(30) = %d %d %d %g %g",
                            j, k, ii, (float)t(ic,26), (float)t(ic,30));
                        int ierrflg = -1;
                        termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
                    }
                }
            }
        }
    }

    int jkpro = jdim*kdim*npl - jdim;

    if (cmn_fvfds.ifds[1] == 0) {
        //
        if (isklton > 0 && i == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "   computing inviscid fluxes, J-direction - flux-vector splitting - block%4d", nbl);
        }
        //
        {
            FortranArray1DRef<double> sj_i1(&sj(1,i,1), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i2(&sj(1,i,2), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i3(&sj(1,i,3), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i4(&sj(1,i,4), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i5(&sj(1,i,5), jdim*kdim*npl, 1);
            FortranArray2DRef<double> t_col21(&t(1,21), nvtq, 5);
            FortranArray2DRef<double> t_col31(&t(1,31), nvtq, 5);
            fluxp_ns::fluxp(sj_i1, sj_i2, sj_i3, sj_i4, sj_i5,
                            t_col21, t_col31, jkpro, t, jkpro, nvtq,
                            nou, bou, nbuf, ibufdim);
        }
        //
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= jkpro; izz++) {
                t(izz, 20+l) = t(izz, 30+l);
            }
        }
        //
        {
            FortranArray1DRef<double> sj_i1(&sj(1,i,1), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i2(&sj(1,i,2), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i3(&sj(1,i,3), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i4(&sj(1,i,4), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i5(&sj(1,i,5), jdim*kdim*npl, 1);
            FortranArray2DRef<double> t_col26(&t(1,26), nvtq, 5);
            FortranArray2DRef<double> t_col31(&t(1,31), nvtq, 5);
            fluxm_ns::fluxm(sj_i1, sj_i2, sj_i3, sj_i4, sj_i5,
                            t_col26, t_col31, jkpro, t, jkpro, nvtq,
                            nou, bou, nbuf, ibufdim);
        }
        //
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= jkpro; izz++) {
                t(izz, 30+l) = t(izz, 30+l) + t(izz, 20+l);
            }
        }

    } else if (cmn_fvfds.ifds[1] == 1) {
        //
        if (isklton > 0 && i == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "   computing inviscid fluxes, J-direction - flux-difference splitting");
        }
        //
        {
            FortranArray1DRef<double> sj_i1(&sj(1,i,1), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i2(&sj(1,i,2), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i3(&sj(1,i,3), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i4(&sj(1,i,4), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i5(&sj(1,i,5), jdim*kdim*npl, 1);
            FortranArray2DRef<double> t_col31(&t(1,31), nvtq, 5);
            FortranArray2DRef<double> t_col26(&t(1,26), nvtq, 5);
            FortranArray2DRef<double> t_col21(&t(1,21), nvtq, 5);
            fhat_ns::fhat(sj_i1, sj_i2, sj_i3, sj_i4, sj_i5,
                          t_col31, t_col26, t_col21, jkpro, nvtq);
        }

    } else {
        //
        if (isklton > 0 && i == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "   computing inviscid fluxes, J-direction - MAPS+ flux splitting");
        }
        //
        {
            FortranArray1DRef<double> sj_i1(&sj(1,i,1), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i2(&sj(1,i,2), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i3(&sj(1,i,3), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i4(&sj(1,i,4), jdim*kdim*npl, 1);
            FortranArray1DRef<double> sj_i5(&sj(1,i,5), jdim*kdim*npl, 1);
            FortranArray2DRef<double> t_col31(&t(1,31), nvtq, 5);
            FortranArray2DRef<double> t_col26(&t(1,26), nvtq, 5);
            FortranArray2DRef<double> t_col21(&t(1,21), nvtq, 5);
            fmaps_ns::fmaps(sj_i1, sj_i2, sj_i3, sj_i4, sj_i5,
                            t_col31, t_col26, t_col21, jkpro, nvtq);
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
            int ic   = (ifa-1)*7;
            int jsb  = iwfa(ic+1);
            int ksb  = iwfa(ic+2);
            int isb  = iwfa(ic+3);
            int jeb  = iwfa(ic+4);
            int keb  = iwfa(ic+5);
            int ieb  = iwfa(ic+6);
            int ifts = iwfa(ic+7);

            if (ifts != lfcc) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " inconsistent summations - stopping");
                int ierrflg = -1;
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }

            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   installing accumulated fluxes on j=%3d at ks,ke,is,ie=%4d%4d%4d%4d",
                    jsb, ksb, keb, isb, ieb);
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
                            int lfc = ifts + (l-1)*(keb-ksb)*(ieb-isb) + (ii-isb)*(keb-ksb);
                            int loc = (ipl-1)*jdim*kdim + (ksb-1)*jdim + jsb;
                            for (int k = ksb; k <= keb-1; k++) {
                                t(loc, 30+l) = wfa(lfc);
                                loc  = loc  + jdim;
                                lfc  = lfc  + 1;
                                lfcc = lfcc + 1;
                            }
                        }
                    } else {
                        //
                        // increment counter for planes not updated on this pass
                        //
                        lfcc = lfcc + 5*(keb-ksb);
                    }
                }
            } // end do 754
        } // end do 7055
    }

    for (int l = 1; l <= 5; l++) {
        for (int izz = 1; izz <= nr; izz++) {
            res(izz, 1, i, l) = t(izz+1, 30+l) - t(izz, 30+l);
        }
    }

    //
    // geometric conservation law terms for deforming grids
    //
    if (idef > 0) {
        double oogmo = 1.0 / ((double)gamma - 1.0);
        for (int izz = 2; izz <= nr; izz++) {
            t(izz, 41) = -(sj(izz+1,i,5)*sj(izz+1,i,4)
                         - sj(izz,  i,5)*sj(izz,  i,4));
        }
        for (int izz = 1; izz <= nr; izz += jdim) {
            t(izz+jdim1, 41) = -(sj(izz+jdim1,  i,5)*sj(izz+jdim1,  i,4)
                                - sj(izz+jdim1-1,i,5)*sj(izz+jdim1-1,i,4));
            t(izz,       41) = -(sj(izz+1,i,5)*sj(izz+1,i,4)
                                - sj(izz,  i,5)*sj(izz,  i,4));
        }
        for (int izz = 1; izz <= nr; izz++) {
            t(izz, 36) = q(izz, 1, i, 1);
            t(izz, 37) = q(izz, 1, i, 1) * q(izz, 1, i, 2);
            t(izz, 38) = q(izz, 1, i, 1) * q(izz, 1, i, 3);
            t(izz, 39) = q(izz, 1, i, 1) * q(izz, 1, i, 4);
            t(izz, 40) = q(izz, 1, i, 5) * oogmo
                       + 0.5 * q(izz, 1, i, 1) * (q(izz, 1, i, 2)*q(izz, 1, i, 2)
                                                  + q(izz, 1, i, 3)*q(izz, 1, i, 3)
                                                  + q(izz, 1, i, 4)*q(izz, 1, i, 4));
        }
        for (int l = 1; l <= 5; l++) {
            for (int izz = 1; izz <= nr; izz++) {
                res(izz, 1, i, l) = res(izz, 1, i, l) + t(izz, 35+l) * t(izz, 41);
            }
        }
    }

    //
    // store finer-grid fluxes for enforcing conservation on coarser meshes
    //
    if (isf == 1) {
        if (jbctyp(1) == 21) {
            //
            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   installing j fluxes into qj0 for edge j=0 for block%3d", nbl);
            }
            //
            // left boundary
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                for (int l = 1; l <= 5; l++) {
                    int jk = (ipl-1)*jdim*kdim + 1 - jdim;
                    for (int k = 1; k <= kdim1; k++) {
                        jk = jk + jdim;
                        qj0(k, ii, l, 2) = t(jk, 30+l);
                    }
                }
            }
        }
        if (jbctyp(2) == 21) {
            //
            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   installing j fluxes into qj0 for edge j=jdim for block%3d", nbl);
            }
            //
            // right boundary
            //
            for (int ipl = 1; ipl <= npl; ipl++) {
                int ii = i + ipl - 1;
                for (int l = 1; l <= 5; l++) {
                    int jk = (ipl-1)*jdim*kdim;
                    for (int k = 1; k <= kdim1; k++) {
                        jk = jk + jdim;
                        qj0(k, ii, l, 4) = t(jk, 30+l);
                    }
                }
            }
        }
    }

    return;
}

} // namespace gfluxr_ns
