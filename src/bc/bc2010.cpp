// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc2010.h"
#include "runtime/common_blocks.h"
#include <cmath>
#include <cstring>
#include <algorithm>

namespace bc2010_ns {

void bc2010(int& jdim, int& kdim, int& idim,
            FortranArray4DRef<double> q,
            FortranArray4DRef<double> qj0,
            FortranArray4DRef<double> qk0,
            FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> sj,
            FortranArray4DRef<double> sk,
            FortranArray4DRef<double> si,
            FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck,
            FortranArray3DRef<double> bci,
            int& ista, int& iend, int& jsta, int& jend, int& ksta, int& kend,
            int& nface,
            FortranArray4DRef<double> tursav,
            FortranArray4DRef<double> tj0,
            FortranArray4DRef<double> tk0,
            FortranArray4DRef<double> ti0,
            FortranArray3DRef<double> vist3d,
            FortranArray4DRef<double> vj0,
            FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0,
            int& mdim, int& ndim,
            FortranArray4DRef<double> bcdata,
            char (&filname)[80],
            int& iuns,
            FortranArray1DRef<int> nou,
            FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim,
            int& myid, int& nummem)
{
    // COMMON block aliases
    float& gamma    = cmn_fluid.gamma;
    float& gm1      = cmn_fluid.gm1;
    float& gp1      = cmn_fluid.gp1;
    float& gm1g     = cmn_fluid.gm1g;
    float& gp1g     = cmn_fluid.gp1g;
    float& ggm1     = cmn_fluid.ggm1;
    float& p0       = cmn_ivals.p0;
    float* tur10    = cmn_ivals.tur10;  // 1-based in Fortran, access as tur10[l-1]
    int& level      = cmn_mgrd.level;
    int& lglobal    = cmn_mgrd.lglobal;
    int* ivisc      = cmn_reyue.ivisc;  // 0-based: ivisc[0],[1],[2]
    float& radtodeg = cmn_conversion.radtodeg;

    // Local dimension variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jend1 = jend - 1;
    int kend1 = kend - 1;
    int iend1 = iend - 1;

    // This bc makes use of only one plane of data
    int ip = 1;

    // Local working variables
    double pte, tte, alpe, betae;
    double unormi, umag_sq, umag, ai, rminus, xmach_local;
    double hti, qa, qb, qc, a_radical, a_plus, a_minus, a_new;
    double u_new, xmach_new, pressure, temperature;
    double ubar, t1;
    int itrflg1, ii, kk, jj, js, l;


    // **************************************************************************
    //  j=1 boundary             nozzle total BCs                   type 2010
    // **************************************************************************

    if (nface == 3) {

        // check to see if turbulence data is input (itrflg1 = 1) or
        // if freestream values are to be used (itrflg1 = 0)
        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        int j = 1;
        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;
            js = (i - ista) * (kend - ksta) + 1;

            for (int k = ksta; k <= kend1; k++) {
                kk = k - ksta + 1;

                pte   = bcdata(kk,ii,ip,1);
                tte   = bcdata(kk,ii,ip,2);
                alpe  = bcdata(kk,ii,ip,3);
                betae = bcdata(kk,ii,ip,4);

                pte   = pte * p0;
                alpe  = alpe / radtodeg;
                betae = betae / radtodeg;

                // The following method is used by FUN3D:
                unormi = q(j,k,i,2)*sj(1,k,i,1) +
                         q(j,k,i,3)*sj(1,k,i,2) +
                         q(j,k,i,4)*sj(1,k,i,3);
                umag_sq = q(j,k,i,2)*q(j,k,i,2) + q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*2.0;
                umag = std::sqrt(umag_sq);
                ai = std::sqrt((double)gamma * q(j,k,i,5) / q(j,k,i,1));
                rminus = std::abs(unormi) - (2.0*ai) / (double)gm1;
                // Note: xmach_local should be < 1 (not checking here)
                xmach_local = umag / ai;
                hti = (q(j,k,i,5)/q(j,k,i,1)) * ((double)gamma/(double)gm1) + 0.5*umag_sq;
                qa = 1.0 + 2.0/(double)gm1;
                qb = 2.0 * rminus;
                qc = 0.5*(double)gm1*rminus*rminus - (double)gm1*hti;
                a_radical = qb*qb - 4.0*qa*qc;
                if ((float)a_radical > 0.0f) {
                    a_plus  = -qb/(2.0*qa) + std::sqrt(a_radical)/(2.0*qa);
                    a_minus = -qb/(2.0*qa) - std::sqrt(a_radical)/(2.0*qa);
                } else if ((float)a_radical < 0.0f) {
                    a_plus  = ai;
                    a_minus = ai;
                } else {
                    a_plus  = -qb/(2.0*qa);
                    a_minus = -qb/(2.0*qa);
                }
                if ((float)a_plus > (float)a_minus) {
                    a_new = a_plus;
                } else {
                    a_new = a_minus;
                }
                u_new    = (2.0*a_new/(double)gm1) + rminus;
                xmach_new = u_new / a_new;
                pressure    = pte * std::pow(1.0 + 0.5*(double)gm1*xmach_new*xmach_new,
                                             -(double)gamma/((double)gamma - 1.0));
                temperature = tte * std::pow(pressure/pte, (double)gm1/(double)gamma);
                qj0(k,i,1,1) = (double)gamma * pressure / temperature;
                qj0(k,i,2,1) = u_new * std::cos(alpe) * std::cos(betae);
                qj0(k,i,3,1) = -u_new * std::sin(betae);
                qj0(k,i,4,1) = u_new * std::sin(alpe) * std::cos(betae);
                qj0(k,i,5,1) = pressure;
                qj0(k,i,1,2) = qj0(k,i,1,1);
                qj0(k,i,2,2) = qj0(k,i,2,1);
                qj0(k,i,3,2) = qj0(k,i,3,1);
                qj0(k,i,4,2) = qj0(k,i,4,1);
                qj0(k,i,5,2) = qj0(k,i,5,1);

                bcj(k,i,1) = 0.0;
            } // k loop
        } // i loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,1) = vist3d(1,k,i);
                    vj0(k,i,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            kk = k - ksta + 1;
                            ubar = -(qj0(k,i,2,1)*sj(1,k,i,1) +
                                     qj0(k,i,3,1)*sj(1,k,i,2) +
                                     qj0(k,i,4,1)*sj(1,k,i,3));
                            t1 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(kk,ii,ip,4+l);
                            if ((float)ubar < 0.0f) {
                                tj0(k,i,l,1) = t1;
                                tj0(k,i,l,2) = t1;
                            } else {
                                tj0(k,i,l,1) = tursav(1,k,i,l);
                                tj0(k,i,l,2) = tj0(k,i,l,1);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==3


    // **************************************************************************
    //  j=jdim boundary          nozzle total BCs                   type 2010
    // **************************************************************************

    if (nface == 4) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        int j = jdim1;
        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;
            js = (i - ista) * (kend - ksta) + 1;

            for (int k = ksta; k <= kend1; k++) {
                kk = k - ksta + 1;

                pte   = bcdata(kk,ii,ip,1);
                tte   = bcdata(kk,ii,ip,2);
                alpe  = bcdata(kk,ii,ip,3);
                betae = bcdata(kk,ii,ip,4);

                pte   = pte * p0;
                alpe  = alpe / radtodeg;
                betae = betae / radtodeg;

                // The following method is used by FUN3D:
                unormi = q(j,k,i,2)*sj(jdim,k,i,1) +
                         q(j,k,i,3)*sj(jdim,k,i,2) +
                         q(j,k,i,4)*sj(jdim,k,i,3);
                umag_sq = q(j,k,i,2)*q(j,k,i,2) + q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*2.0;
                umag = std::sqrt(umag_sq);
                ai = std::sqrt((double)gamma * q(j,k,i,5) / q(j,k,i,1));
                rminus = std::abs(unormi) - (2.0*ai) / (double)gm1;
                xmach_local = umag / ai;
                hti = (q(j,k,i,5)/q(j,k,i,1)) * ((double)gamma/(double)gm1) + 0.5*umag_sq;
                qa = 1.0 + 2.0/(double)gm1;
                qb = 2.0 * rminus;
                qc = 0.5*(double)gm1*rminus*rminus - (double)gm1*hti;
                a_radical = qb*qb - 4.0*qa*qc;
                if ((float)a_radical > 0.0f) {
                    a_plus  = -qb/(2.0*qa) + std::sqrt(a_radical)/(2.0*qa);
                    a_minus = -qb/(2.0*qa) - std::sqrt(a_radical)/(2.0*qa);
                } else if ((float)a_radical < 0.0f) {
                    a_plus  = ai;
                    a_minus = ai;
                } else {
                    a_plus  = -qb/(2.0*qa);
                    a_minus = -qb/(2.0*qa);
                }
                if ((float)a_plus > (float)a_minus) {
                    a_new = a_plus;
                } else {
                    a_new = a_minus;
                }
                u_new     = (2.0*a_new/(double)gm1) + rminus;
                xmach_new = u_new / a_new;
                pressure    = pte * std::pow(1.0 + 0.5*(double)gm1*xmach_new*xmach_new,
                                             -(double)gamma/((double)gamma - 1.0));
                temperature = tte * std::pow(pressure/pte, (double)gm1/(double)gamma);
                qj0(k,i,1,3) = (double)gamma * pressure / temperature;
                qj0(k,i,2,3) = u_new * std::cos(alpe) * std::cos(betae);
                qj0(k,i,3,3) = -u_new * std::sin(betae);
                qj0(k,i,4,3) = u_new * std::sin(alpe) * std::cos(betae);
                qj0(k,i,5,3) = pressure;
                qj0(k,i,1,4) = qj0(k,i,1,3);
                qj0(k,i,2,4) = qj0(k,i,2,3);
                qj0(k,i,3,4) = qj0(k,i,3,3);
                qj0(k,i,4,4) = qj0(k,i,4,3);
                qj0(k,i,5,4) = qj0(k,i,5,3);

                bcj(k,i,2) = 0.0;
            } // k loop
        } // i loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int k = ksta; k <= kend1; k++) {
                    vj0(k,i,1,3) = vist3d(jdim-1,k,i);
                    vj0(k,i,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        ii = i - ista + 1;
                        for (int k = ksta; k <= kend1; k++) {
                            kk = k - ksta + 1;
                            ubar = qj0(k,i,2,3)*sj(jdim,k,i,1) +
                                   qj0(k,i,3,3)*sj(jdim,k,i,2) +
                                   qj0(k,i,4,3)*sj(jdim,k,i,3);
                            t1 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(kk,ii,ip,4+l);
                            if ((float)ubar < 0.0f) {
                                tj0(k,i,l,3) = t1;
                                tj0(k,i,l,4) = t1;
                            } else {
                                tj0(k,i,l,3) = tursav(jdim-1,k,i,l);
                                tj0(k,i,l,4) = tj0(k,i,l,3);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==4


    // **************************************************************************
    //  k=1 boundary             nozzle total BCs                   type 2010
    // **************************************************************************

    if (nface == 5) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        int k = 1;
        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;
            js = (i - ista) * (jend - jsta) + 1;
            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;

                pte   = bcdata(jj,ii,ip,1);
                tte   = bcdata(jj,ii,ip,2);
                alpe  = bcdata(jj,ii,ip,3);
                betae = bcdata(jj,ii,ip,4);

                pte   = pte * p0;
                alpe  = alpe / radtodeg;
                betae = betae / radtodeg;

                // The following method is used by FUN3D:
                unormi = q(j,k,i,2)*sk(j,1,i,1) +
                         q(j,k,i,3)*sk(j,1,i,2) +
                         q(j,k,i,4)*sk(j,1,i,3);
                umag_sq = q(j,k,i,2)*q(j,k,i,2) + q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*2.0;
                umag = std::sqrt(umag_sq);
                ai = std::sqrt((double)gamma * q(j,k,i,5) / q(j,k,i,1));
                rminus = std::abs(unormi) - (2.0*ai) / (double)gm1;
                xmach_local = umag / ai;
                hti = (q(j,k,i,5)/q(j,k,i,1)) * ((double)gamma/(double)gm1) + 0.5*umag_sq;
                qa = 1.0 + 2.0/(double)gm1;
                qb = 2.0 * rminus;
                qc = 0.5*(double)gm1*rminus*rminus - (double)gm1*hti;
                a_radical = qb*qb - 4.0*qa*qc;
                if ((float)a_radical > 0.0f) {
                    a_plus  = -qb/(2.0*qa) + std::sqrt(a_radical)/(2.0*qa);
                    a_minus = -qb/(2.0*qa) - std::sqrt(a_radical)/(2.0*qa);
                } else if ((float)a_radical < 0.0f) {
                    a_plus  = ai;
                    a_minus = ai;
                } else {
                    a_plus  = -qb/(2.0*qa);
                    a_minus = -qb/(2.0*qa);
                }
                if ((float)a_plus > (float)a_minus) {
                    a_new = a_plus;
                } else {
                    a_new = a_minus;
                }
                u_new     = (2.0*a_new/(double)gm1) + rminus;
                xmach_new = u_new / a_new;
                pressure    = pte * std::pow(1.0 + 0.5*(double)gm1*xmach_new*xmach_new,
                                             -(double)gamma/((double)gamma - 1.0));
                temperature = tte * std::pow(pressure/pte, (double)gm1/(double)gamma);
                qk0(j,i,1,1) = (double)gamma * pressure / temperature;
                qk0(j,i,2,1) = u_new * std::cos(alpe) * std::cos(betae);
                qk0(j,i,3,1) = -u_new * std::sin(betae);
                qk0(j,i,4,1) = u_new * std::sin(alpe) * std::cos(betae);
                qk0(j,i,5,1) = pressure;
                qk0(j,i,1,2) = qk0(j,i,1,1);
                qk0(j,i,2,2) = qk0(j,i,2,1);
                qk0(j,i,3,2) = qk0(j,i,3,1);
                qk0(j,i,4,2) = qk0(j,i,4,1);
                qk0(j,i,5,2) = qk0(j,i,5,1);

                bck(j,i,1) = 0.0;
            } // j loop
        } // i loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,1) = vist3d(j,1,i);
                    vk0(j,i,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            jj = j - jsta + 1;
                            ubar = -(qk0(j,i,2,1)*sk(j,1,i,1) +
                                     qk0(j,i,3,1)*sk(j,1,i,2) +
                                     qk0(j,i,4,1)*sk(j,1,i,3));
                            t1 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,ii,ip,4+l);
                            if ((float)ubar < 0.0f) {
                                tk0(j,i,l,1) = t1;
                                tk0(j,i,l,2) = t1;
                            } else {
                                tk0(j,i,l,1) = tursav(j,1,i,l);
                                tk0(j,i,l,2) = tk0(j,i,l,1);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==5


    // **************************************************************************
    //  k=kdim boundary          nozzle total BCs                   type 2010
    // **************************************************************************

    if (nface == 6) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        int k = kdim1;
        for (int i = ista; i <= iend1; i++) {
            ii = i - ista + 1;
            js = (i - ista) * (jend - jsta) + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;

                pte   = bcdata(jj,ii,ip,1);
                tte   = bcdata(jj,ii,ip,2);
                alpe  = bcdata(jj,ii,ip,3);
                betae = bcdata(jj,ii,ip,4);

                pte   = pte * p0;
                alpe  = alpe / radtodeg;
                betae = betae / radtodeg;

                // The following method is used by FUN3D:
                unormi = q(j,k,i,2)*sk(j,kdim,i,1) +
                         q(j,k,i,3)*sk(j,kdim,i,2) +
                         q(j,k,i,4)*sk(j,kdim,i,3);
                umag_sq = q(j,k,i,2)*q(j,k,i,2) + q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*2.0;
                umag = std::sqrt(umag_sq);
                ai = std::sqrt((double)gamma * q(j,k,i,5) / q(j,k,i,1));
                rminus = std::abs(unormi) - (2.0*ai) / (double)gm1;
                xmach_local = umag / ai;
                hti = (q(j,k,i,5)/q(j,k,i,1)) * ((double)gamma/(double)gm1) + 0.5*umag_sq;
                qa = 1.0 + 2.0/(double)gm1;
                qb = 2.0 * rminus;
                qc = 0.5*(double)gm1*rminus*rminus - (double)gm1*hti;
                a_radical = qb*qb - 4.0*qa*qc;
                if ((float)a_radical > 0.0f) {
                    a_plus  = -qb/(2.0*qa) + std::sqrt(a_radical)/(2.0*qa);
                    a_minus = -qb/(2.0*qa) - std::sqrt(a_radical)/(2.0*qa);
                } else if ((float)a_radical < 0.0f) {
                    a_plus  = ai;
                    a_minus = ai;
                } else {
                    a_plus  = -qb/(2.0*qa);
                    a_minus = -qb/(2.0*qa);
                }
                if ((float)a_plus > (float)a_minus) {
                    a_new = a_plus;
                } else {
                    a_new = a_minus;
                }
                u_new     = (2.0*a_new/(double)gm1) + rminus;
                xmach_new = u_new / a_new;
                pressure    = pte * std::pow(1.0 + 0.5*(double)gm1*xmach_new*xmach_new,
                                             -(double)gamma/((double)gamma - 1.0));
                temperature = tte * std::pow(pressure/pte, (double)gm1/(double)gamma);
                qk0(j,i,1,3) = (double)gamma * pressure / temperature;
                qk0(j,i,2,3) = u_new * std::cos(alpe) * std::cos(betae);
                qk0(j,i,3,3) = -u_new * std::sin(betae);
                qk0(j,i,4,3) = u_new * std::sin(alpe) * std::cos(betae);
                qk0(j,i,5,3) = pressure;
                qk0(j,i,1,4) = qk0(j,i,1,3);
                qk0(j,i,2,4) = qk0(j,i,2,3);
                qk0(j,i,3,4) = qk0(j,i,3,3);
                qk0(j,i,4,4) = qk0(j,i,4,3);
                qk0(j,i,5,4) = qk0(j,i,5,3);

                bck(j,i,2) = 0.0;
            } // j loop
        } // i loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int i = ista; i <= iend1; i++) {
                for (int j = jsta; j <= jend1; j++) {
                    vk0(j,i,1,3) = vist3d(j,kdim-1,i);
                    vk0(j,i,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int i = ista; i <= iend1; i++) {
                        ii = i - ista + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            jj = j - jsta + 1;
                            ubar = qk0(j,i,2,3)*sk(j,kdim,i,1) +
                                   qk0(j,i,3,3)*sk(j,kdim,i,2) +
                                   qk0(j,i,4,3)*sk(j,kdim,i,3);
                            t1 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,ii,ip,4+l);
                            if ((float)ubar < 0.0f) {
                                tk0(j,i,l,3) = t1;
                                tk0(j,i,l,4) = t1;
                            } else {
                                tk0(j,i,l,3) = tursav(j,kdim-1,i,l);
                                tk0(j,i,l,4) = tk0(j,i,l,3);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==6


    // **************************************************************************
    //  i=1 boundary             nozzle total BCs                   type 2010
    // **************************************************************************

    if (nface == 1) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        int i = 1;
        for (int k = ksta; k <= kend1; k++) {
            kk = k - ksta + 1;
            js = (k - ksta) * (jend - jsta) + 1;

            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;

                pte   = bcdata(jj,kk,ip,1);
                tte   = bcdata(jj,kk,ip,2);
                alpe  = bcdata(jj,kk,ip,3);
                betae = bcdata(jj,kk,ip,4);

                pte   = pte * p0;
                alpe  = alpe / radtodeg;
                betae = betae / radtodeg;

                // The following method is used by FUN3D:
                unormi = q(j,k,i,2)*si(j,k,1,1) +
                         q(j,k,i,3)*si(j,k,1,2) +
                         q(j,k,i,4)*si(j,k,1,3);
                umag_sq = q(j,k,i,2)*q(j,k,i,2) + q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*2.0;
                umag = std::sqrt(umag_sq);
                ai = std::sqrt((double)gamma * q(j,k,i,5) / q(j,k,i,1));
                rminus = std::abs(unormi) - (2.0*ai) / (double)gm1;
                xmach_local = umag / ai;
                hti = (q(j,k,i,5)/q(j,k,i,1)) * ((double)gamma/(double)gm1) + 0.5*umag_sq;
                qa = 1.0 + 2.0/(double)gm1;
                qb = 2.0 * rminus;
                qc = 0.5*(double)gm1*rminus*rminus - (double)gm1*hti;
                a_radical = qb*qb - 4.0*qa*qc;
                if ((float)a_radical > 0.0f) {
                    a_plus  = -qb/(2.0*qa) + std::sqrt(a_radical)/(2.0*qa);
                    a_minus = -qb/(2.0*qa) - std::sqrt(a_radical)/(2.0*qa);
                } else if ((float)a_radical < 0.0f) {
                    a_plus  = ai;
                    a_minus = ai;
                } else {
                    a_plus  = -qb/(2.0*qa);
                    a_minus = -qb/(2.0*qa);
                }
                if ((float)a_plus > (float)a_minus) {
                    a_new = a_plus;
                } else {
                    a_new = a_minus;
                }
                u_new     = (2.0*a_new/(double)gm1) + rminus;
                xmach_new = u_new / a_new;
                pressure    = pte * std::pow(1.0 + 0.5*(double)gm1*xmach_new*xmach_new,
                                             -(double)gamma/((double)gamma - 1.0));
                temperature = tte * std::pow(pressure/pte, (double)gm1/(double)gamma);
                qi0(j,k,1,1) = (double)gamma * pressure / temperature;
                qi0(j,k,2,1) = u_new * std::cos(alpe) * std::cos(betae);
                qi0(j,k,3,1) = -u_new * std::sin(betae);
                qi0(j,k,4,1) = u_new * std::sin(alpe) * std::cos(betae);
                qi0(j,k,5,1) = pressure;
                qi0(j,k,1,2) = qi0(j,k,1,1);
                qi0(j,k,2,2) = qi0(j,k,2,1);
                qi0(j,k,3,2) = qi0(j,k,3,1);
                qi0(j,k,4,2) = qi0(j,k,4,1);
                qi0(j,k,5,2) = qi0(j,k,5,1);

                bci(j,k,1) = 0.0;
            } // j loop
        } // k loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,1) = vist3d(j,k,1);
                    vi0(j,k,1,2) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            jj = j - jsta + 1;
                            ubar = -(qi0(j,k,2,1)*si(j,k,1,1) +
                                     qi0(j,k,3,1)*si(j,k,1,2) +
                                     qi0(j,k,4,1)*si(j,k,1,3));
                            t1 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,kk,ip,4+l);
                            if ((float)ubar < 0.0f) {
                                ti0(j,k,l,1) = t1;
                                ti0(j,k,l,2) = t1;
                            } else {
                                ti0(j,k,l,1) = tursav(j,k,1,l);
                                ti0(j,k,l,2) = ti0(j,k,l,1);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==1


    // **************************************************************************
    //  i=idim boundary          nozzle total BCs                   type 2010
    // **************************************************************************

    if (nface == 2) {

        itrflg1 = 0;
        if ((float)bcdata(1,1,ip,5) > -1.e10f) itrflg1 = 1;

        int i = idim1;
        for (int k = ksta; k <= kend1; k++) {
            kk = k - ksta + 1;
            js = (k - ksta) * (jend - jsta) + 1;
            for (int j = jsta; j <= jend1; j++) {
                jj = j - jsta + 1;

                pte   = bcdata(jj,kk,ip,1);
                tte   = bcdata(jj,kk,ip,2);
                alpe  = bcdata(jj,kk,ip,3);
                betae = bcdata(jj,kk,ip,4);

                pte   = pte * p0;
                alpe  = alpe / radtodeg;
                betae = betae / radtodeg;

                // The following method is used by FUN3D:
                unormi = q(j,k,i,2)*si(j,k,idim,1) +
                         q(j,k,i,3)*si(j,k,idim,2) +
                         q(j,k,i,4)*si(j,k,idim,3);
                umag_sq = q(j,k,i,2)*q(j,k,i,2) + q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*2.0;
                umag = std::sqrt(umag_sq);
                ai = std::sqrt((double)gamma * q(j,k,i,5) / q(j,k,i,1));
                rminus = std::abs(unormi) - (2.0*ai) / (double)gm1;
                xmach_local = umag / ai;
                hti = (q(j,k,i,5)/q(j,k,i,1)) * ((double)gamma/(double)gm1) + 0.5*umag_sq;
                qa = 1.0 + 2.0/(double)gm1;
                qb = 2.0 * rminus;
                qc = 0.5*(double)gm1*rminus*rminus - (double)gm1*hti;
                a_radical = qb*qb - 4.0*qa*qc;
                if ((float)a_radical > 0.0f) {
                    a_plus  = -qb/(2.0*qa) + std::sqrt(a_radical)/(2.0*qa);
                    a_minus = -qb/(2.0*qa) - std::sqrt(a_radical)/(2.0*qa);
                } else if ((float)a_radical < 0.0f) {
                    a_plus  = ai;
                    a_minus = ai;
                } else {
                    a_plus  = -qb/(2.0*qa);
                    a_minus = -qb/(2.0*qa);
                }
                if ((float)a_plus > (float)a_minus) {
                    a_new = a_plus;
                } else {
                    a_new = a_minus;
                }
                u_new     = (2.0*a_new/(double)gm1) + rminus;
                xmach_new = u_new / a_new;
                pressure    = pte * std::pow(1.0 + 0.5*(double)gm1*xmach_new*xmach_new,
                                             -(double)gamma/((double)gamma - 1.0));
                temperature = tte * std::pow(pressure/pte, (double)gm1/(double)gamma);
                qi0(j,k,1,3) = (double)gamma * pressure / temperature;
                qi0(j,k,2,3) = u_new * std::cos(alpe) * std::cos(betae);
                qi0(j,k,3,3) = -u_new * std::sin(betae);
                qi0(j,k,4,3) = u_new * std::sin(alpe) * std::cos(betae);
                qi0(j,k,5,3) = pressure;
                qi0(j,k,1,4) = qi0(j,k,1,3);
                qi0(j,k,2,4) = qi0(j,k,2,3);
                qi0(j,k,3,4) = qi0(j,k,3,3);
                qi0(j,k,4,4) = qi0(j,k,4,3);
                qi0(j,k,5,4) = qi0(j,k,5,3);

                bci(j,k,2) = 0.0;
            } // j loop
        } // k loop

        if (ivisc[2] >= 2 || ivisc[1] >= 2 || ivisc[0] >= 2) {
            for (int k = ksta; k <= kend1; k++) {
                for (int j = jsta; j <= jend1; j++) {
                    vi0(j,k,1,3) = vist3d(j,k,idim-1);
                    vi0(j,k,1,4) = 0.0;
                }
            }
        }
        // only need to do advanced model turbulence B.C.s on finest grid
        if (level >= lglobal) {
            if (ivisc[2] >= 4 || ivisc[1] >= 4 || ivisc[0] >= 4) {
                for (l = 1; l <= nummem; l++) {
                    for (int k = ksta; k <= kend1; k++) {
                        kk = k - ksta + 1;
                        for (int j = jsta; j <= jend1; j++) {
                            jj = j - jsta + 1;
                            ubar = qi0(j,k,2,3)*si(j,k,idim,1) +
                                   qi0(j,k,3,3)*si(j,k,idim,2) +
                                   qi0(j,k,4,3)*si(j,k,idim,3);
                            t1 = (1 - itrflg1)*(double)tur10[l-1] + itrflg1*bcdata(jj,kk,ip,4+l);
                            if ((float)ubar < 0.0f) {
                                ti0(j,k,l,3) = t1;
                                ti0(j,k,l,4) = t1;
                            } else {
                                ti0(j,k,l,3) = tursav(j,k,idim-1,l);
                                ti0(j,k,l,4) = ti0(j,k,l,3);
                            }
                        }
                    }
                }
            }
        }

    } // end nface==2

} // end bc2010()



void bc(int& ntime, int& nbl,
        FortranArray2DRef<int> lw,
        FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w,
        int& mgwk,
        FortranArray1DRef<double> wk,
        int& nwork,
        double& cl,
        FortranArray1DRef<int> nou,
        FortranArray2DRef<char[120]> bou,
        int& nbuf,
        int& ibufdim,
        int& maxbl,
        int& maxgr,
        int& maxseg,
        FortranArray1DRef<int> itrans,
        FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm,
        FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg,
        FortranArray1DRef<int> nbci0,
        FortranArray1DRef<int> nbcj0,
        FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim,
        FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim,
        FortranArray4DRef<int> ibcinfo,
        FortranArray4DRef<int> jbcinfo,
        FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei,
        FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek,
        FortranArray3DRef<int> lwdat,
        int& myid,
        FortranArray1DRef<int> idimg,
        FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg,
        FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil,
        int& nummem)
{
    int iuns = 0;

    for (int ibl = 1; ibl <= nbl; ibl++) {
        int jdim = jdimg(ibl);
        int kdim = kdimg(ibl);
        int idim = idimg(ibl);

        // Build array views from w[] using lw offsets
        FortranArray4DRef<double> q(   &w(lw(1,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray4DRef<double> qj0( &w(lw(2,ibl)),  kdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qk0( &w(lw(3,ibl)),  jdim,   idim-1, 5,      4);
        FortranArray4DRef<double> qi0( &w(lw(4,ibl)),  jdim,   kdim,   5,      4);
        FortranArray4DRef<double> sj(  &w(lw(5,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> sk(  &w(lw(6,ibl)),  jdim,   kdim,   idim-1, 5);
        FortranArray4DRef<double> si(  &w(lw(7,ibl)),  jdim,   kdim,   idim,   5);
        FortranArray3DRef<double> bcj( &w(lw(8,ibl)),  kdim,   idim-1, 2);
        FortranArray3DRef<double> bck( &w(lw(9,ibl)),  jdim,   idim-1, 2);
        FortranArray3DRef<double> bci( &w(lw(10,ibl)), jdim,   kdim,   2);

        // Build turbulence/viscosity array views from wk[] using lw2 offsets
        FortranArray4DRef<double> tursav(&w(lw2(1,ibl)), jdim,   kdim,   idim,   nummem);
        FortranArray4DRef<double> tj0(   &w(lw2(2,ibl)), kdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> tk0(   &w(lw2(3,ibl)), jdim,   idim-1, nummem, 4);
        FortranArray4DRef<double> ti0(   &w(lw2(4,ibl)), jdim,   kdim,   nummem, 4);
        FortranArray3DRef<double> vist3d(&w(lw2(5,ibl)), jdim,   kdim,   idim);
        FortranArray4DRef<double> vj0(   &w(lw2(6,ibl)), kdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vk0(   &w(lw2(7,ibl)), jdim,   idim-1, 1,      4);
        FortranArray4DRef<double> vi0(   &w(lw2(8,ibl)), jdim,   kdim,   1,      4);

        // mdim and ndim for bcdata
        int mdim = nbcjdim(ibl);
        int ndim = nbcidim(ibl);

        // Loop over j-face BC segments
        for (int iseg = 1; iseg <= nbcj0(ibl); iseg++) {
            if (jbcinfo(ibl,iseg,1,1) == 2010) {
                int ista  = jbcinfo(ibl,iseg,1,2);
                int iend  = jbcinfo(ibl,iseg,1,3);
                int jsta  = jbcinfo(ibl,iseg,1,4);
                int jend  = jbcinfo(ibl,iseg,1,5);
                int ksta  = jbcinfo(ibl,iseg,1,6);
                int kend  = jbcinfo(ibl,iseg,1,7);
                int nface = jbcinfo(ibl,iseg,1,8);
                int lbcdata = lwdat(ibl,iseg,1);
                FortranArray4DRef<double> bcdata(&w(lbcdata), mdim, ndim, 2, 12);
                int ibcfil = bcfilej(ibl,iseg,1);
                char (&filname)[80] = bcfiles(ibcfil);
                bc2010(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // Loop over k-face BC segments
        mdim = nbckdim(ibl);
        ndim = nbcidim(ibl);
        for (int iseg = 1; iseg <= nbck0(ibl); iseg++) {
            if (kbcinfo(ibl,iseg,1,1) == 2010) {
                int ista  = kbcinfo(ibl,iseg,1,2);
                int iend  = kbcinfo(ibl,iseg,1,3);
                int jsta  = kbcinfo(ibl,iseg,1,4);
                int jend  = kbcinfo(ibl,iseg,1,5);
                int ksta  = kbcinfo(ibl,iseg,1,6);
                int kend  = kbcinfo(ibl,iseg,1,7);
                int nface = kbcinfo(ibl,iseg,1,8);
                int lbcdata = lwdat(ibl,iseg,2);
                FortranArray4DRef<double> bcdata(&w(lbcdata), mdim, ndim, 2, 12);
                int ibcfil = bcfilek(ibl,iseg,1);
                char (&filname)[80] = bcfiles(ibcfil);
                bc2010(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

        // Loop over i-face BC segments
        mdim = nbcjdim(ibl);
        ndim = nbckdim(ibl);
        for (int iseg = 1; iseg <= nbci0(ibl); iseg++) {
            if (ibcinfo(ibl,iseg,1,1) == 2010) {
                int ista  = ibcinfo(ibl,iseg,1,2);
                int iend  = ibcinfo(ibl,iseg,1,3);
                int jsta  = ibcinfo(ibl,iseg,1,4);
                int jend  = ibcinfo(ibl,iseg,1,5);
                int ksta  = ibcinfo(ibl,iseg,1,6);
                int kend  = ibcinfo(ibl,iseg,1,7);
                int nface = ibcinfo(ibl,iseg,1,8);
                int lbcdata = lwdat(ibl,iseg,3);
                FortranArray4DRef<double> bcdata(&w(lbcdata), mdim, ndim, 2, 12);
                int ibcfil = bcfilei(ibl,iseg,1);
                char (&filname)[80] = bcfiles(ibcfil);
                bc2010(jdim, kdim, idim, q, qj0, qk0, qi0, sj, sk, si,
                       bcj, bck, bci,
                       ista, iend, jsta, jend, ksta, kend, nface,
                       tursav, tj0, tk0, ti0, vist3d, vj0, vk0, vi0,
                       mdim, ndim, bcdata, filname, iuns,
                       nou, bou, nbuf, ibufdim, myid, nummem);
            }
        }

    } // end ibl loop

} // end bc()

} // namespace bc2010_ns
