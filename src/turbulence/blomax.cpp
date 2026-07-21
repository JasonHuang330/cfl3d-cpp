// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "blomax.h"
#include "q8smax.h"
#include "q8smin.h"
#include "termn8.h"
#include "ccomplex.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace blomax_ns {

void blomax(int& jdim, int& kdim, int& idim, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> vor, FortranArray3DRef<double> snj0, FortranArray3DRef<double> snk0, FortranArray3DRef<double> sni0, FortranArray3DRef<double> snjm, FortranArray3DRef<double> snkm, FortranArray3DRef<double> snim, FortranArray3DRef<double> vist3d, FortranArray3DRef<double> eoms, int& iprint, int& inmx, FortranArray1DRef<double> eomui, FortranArray1DRef<double> fbl, FortranArray1DRef<double> rhon, FortranArray1DRef<double> amun, FortranArray1DRef<double> vortn, FortranArray1DRef<double> disn, FortranArray1DRef<double> utot, FortranArray1DRef<double> eomuo, FortranArray1DRef<double> damp, int& nblt, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray3DRef<double> blank, int& iover, FortranArray3DRef<double> bci, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    // DATA statement constants
    const double aplus   = 26.0e0;
    const double ccp     = 1.6e0;
    const double ckleb   = 0.3e0;
    const double cwk     = 1.0e0;
    const double vk      = 0.4e0;
    const double clauser = 0.0168e0;

    // COMMON block aliases
    int32_t* ideg   = cmn_degshf.ideg;   // 1-based: ideg[1-1]=ideg[0], etc.
    float& gamma    = cmn_fluid.gamma;
    float& cbar     = cmn_fluid2.cbar;
    float& xmach    = cmn_info.xmach;
    int32_t* iwf    = cmn_wallfun.iwf;   // 1-based: iwf[1-1]=iwf[0]
    int32_t& ilamlo = cmn_lam.ilamlo;
    int32_t& ilamhi = cmn_lam.ilamhi;
    int32_t& jlamlo = cmn_lam.jlamlo;
    int32_t& jlamhi = cmn_lam.jlamhi;
    int32_t& klamlo = cmn_lam.klamlo;
    int32_t& klamhi = cmn_lam.klamhi;
    float& reue     = cmn_reyue.reue;
    float& tinf     = cmn_reyue.tinf;
    int32_t* ivisc  = cmn_reyue.ivisc;  // 1-based: ivisc[1-1]=ivisc[0]
    int32_t& isklton = cmn_sklton.isklton;
    int32_t& nfreeze = cmn_turbconv.nfreeze;
    int32_t& icyc    = cmn_mgrd.icyc;

    // Local variables
    int ilfreq1, ilfreq2, ipfreq1, ipfreq2;
    int jdim1, kdim1, idim1;
    double aplusi, xmi, ckout, c2b, c2bp;
    int iwall, ifit;
    int ibcjmin, ibcjmax, ibckmin, ibckmax, ibcimin, ibcimax;
    int i, j, k, in, nnn;
    int inswtch, infmax, in1, in2, inmax, inmax1;
    int kloop, kloop1, inmaxt, inmaxt1;
    int jloop, jloop1;
    int iloop, iloop1;
    int in21, ihead, nss;
    double t5, t6, amixl, yplusoy, yplus, term;
    double fblmax, ymax, udif, fwake, fwake2, utmax, utmin;
    double dfm, dfp, eomax, utw, s1, s2;
    int nsmurf;


    // Lines 101-111: Early return if icyc <= nfreeze
    if (icyc <= nfreeze) {
        if (isklton > 0) {
            nss = std::min(icyc <= 0 ? nfreeze : (int)std::min((double)icyc, (double)nfreeze), nfreeze);
            // Fortran: nss=min(ncyc,nfreeze)
            // ncyc is cmn_mgrd.ncyc
            nss = std::min(cmn_mgrd.ncyc, nfreeze);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " B-L turbulence model is frozen for %5d iterations or subits", nss);
        }
        return;
    }

    // Lines 135-142: print frequency setup
    ilfreq1 = 1;
    ilfreq2 = 5;
    ipfreq1 = 10;
    ipfreq2 = 10;
    if (idim == 2 || jdim == 2 || kdim == 2) {
        ilfreq2 = 1;
        ipfreq2 = 1;
    }

    // Lines 144-161: dimension and constant setup
    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    aplusi = 1.0e0 / aplus;
    xmi    = 1.0e0 / (double)xmach;
    ckout  = (double)reue * xmi * clauser * ccp;
    c2b    = (double)cbar / (double)tinf;
    c2bp   = c2b + 1.0e0;

    // iwall = 0 for yplus damping scaled with local values
    iwall = 0;
    // ifit  = 0 take fmax as largest discrete value of f(in)
    ifit  = 0;

    // Lines 163-186: determine which walls to use, based on bcj,bck,bci
    ibcjmin = 0;
    ibcjmax = 0;
    for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim; k++) {
            if ((float)bcj(k,i,1) > 0.f) ibcjmin = 1;
            if ((float)bcj(k,i,2) > 0.f) ibcjmax = 1;
        }
    }
    ibckmin = 0;
    ibckmax = 0;
    for (i = 1; i <= idim1; i++) {
        for (j = 1; j <= jdim; j++) {
            if ((float)bck(j,i,1) > 0.f) ibckmin = 1;
            if ((float)bck(j,i,2) > 0.f) ibckmax = 1;
        }
    }
    ibcimin = 0;
    ibcimax = 0;
    for (k = 1; k <= kdim; k++) {
        for (j = 1; j <= jdim; j++) {
            if ((float)bci(j,k,1) > 0.f) ibcimin = 1;
            if ((float)bci(j,k,2) > 0.f) ibcimax = 1;
        }
    }


    // Lines 187-794: ivisc(3)>1: k-wall turbulence
    if (ivisc[3-1] > 1) {
        // evaluate turbulent viscosity along normal to k=0 wall
        // defaults to using min face for distance if neither min nor max contains wall
        if (ibckmin == 1 || ibckmax == 0) {

            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " computing turb viscosity using Baldwin-Lomax, K=1   , block=%4d", nblt);
            }
            if (isklton > 0 && ideg[3-1] == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned on");
            }
            if (isklton > 0 && ideg[3-1] != 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned off");
            }

            if (inmx < kdim) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "error: inmx, kdim=%5d%5d", inmx, kdim);
                { int m1 = -1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
            }

            ihead = 0;
            // do 1000 i=1,idim1; do 1000 j=1,jdim1
            for (i = 1; i <= idim1; i++) {
            for (j = 1; j <= jdim1; j++) {

                // bound turbulent region
                kloop  = (int)(0.80 * kdim);
                kloop1 = kloop - 1;
                inmax  = kloop;
                inmax1 = inmax - 1;

                // bound search region for fmax
                in1 = (int)(inmax1 * 0.20) + 1;
                if (iwf[3-1] > 0) {
                    in1 = 3;
                }
                in2 = (int)(inmax1 * 0.80) + 1;

                // Chimera scheme modification
                if (iover == 1) {
                    if ((float)blank(j,in1,i) == 0.f) {
                        for (k = 1; k <= kdim1; k++) {
                            vist3d(j,k,i) = 0.0e0;
                        }
                        if (iprint >= 1) {
                            if (j == j/ilfreq1*ilfreq1 && i == i/ilfreq2*ilfreq2) {
                                if (ihead == 0) {
                                    nou(3) = std::min(nou(3)+1, ibufdim);
                                    std::snprintf(bou(nou(3),3), 120,
                                        " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                                }
                                nou(3) = std::min(nou(3)+1, ibufdim);
                                std::snprintf(bou(nou(3),3), 120,
                                    "%9d%5d     skipping eddy viscosity evaluation at this station due to hole at-near wall", j, i);
                                ihead = ihead + 1;
                            }
                        }
                        goto label1000;
                    }
                    in21 = in2;
                    for (in = in1; in <= in2; in++) {
                        if ((float)blank(j,in,i) == 0.f) goto label1524;
                        in21 = in;
                    }
                    label1524:;
                    in2 = in21;
                }


                // viscosity thru sutherland relation at cell centers
                for (k = 1; k <= kloop1; k++) {
                    t5         = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
                    t6         = std::sqrt(t5);
                    amun(k+1)  = c2bp * t5 * t6 / (c2b + t5);
                }
                for (k = 1; k <= kloop1; k++) {
                    rhon(k+1)  = q(j,k,i,1);
                    vortn(k+1) = vor(j,k,i);
                    { double tmp = snk0(j,k,i); disn(k+1) = std::abs(tmp); }
                    utot(k+1)  = std::sqrt( q(j,k,i,2)*q(j,k,i,2)
                                          + q(j,k,i,3)*q(j,k,i,3)
                                          + q(j,k,i,4)*q(j,k,i,4) );
                }

                rhon(1)  = rhon(2);
                amun(1)  = amun(2);
                disn(1)  = 0.0;
                utot(1)  = utot(2);
                vortn(1) = vortn(2);

                if ((float)bck(j,i,1) == 1.f) {
                    // modify wall values in no-slip region
                    t5 = (double)gamma * qk0(j,i,5,1) / qk0(j,i,1,1);
                    t6 = std::sqrt(t5);
                    amun(1)  = c2bp * t5 * t6 / (c2b + t5);
                    rhon(1)  = qk0(j,i,1,1);
                    utot(1)  = std::sqrt( qk0(j,i,2,1)*qk0(j,i,2,1)
                                        + qk0(j,i,3,1)*qk0(j,i,3,1)
                                        + qk0(j,i,4,1)*qk0(j,i,4,1) );
                    vortn(1) = (utot(2) - utot(1)) / disn(2);
                }

                // convenient groupings and initial values at k=1
                { double tmp = rhon(1)*vortn(1)*(double)reue*xmi/amun(1);
                  yplusoy = std::abs(tmp); }
                yplusoy = std::sqrt(yplusoy) * aplusi;
                yplus   = yplusoy * disn(2) * aplus * 2.0;

                if ((float)bck(j,i,1) == 0.f) {
                    // wake damping
                    term = 1.0 - std::exp(-50.0e0);
                    for (in = 2; in <= inmax; in++) {
                        damp(in) = term;
                    }
                } else {
                    // wall damping
                    for (in = 2; in <= inmax; in++) {
                        if (iwall == 0) {
                            term = yplusoy * disn(in) * std::sqrt(rhon(in)/rhon(1))
                                   * (amun(1)/amun(in));
                        } else {
                            term = yplusoy * disn(in);
                        }
                        damp(in) = 1.0e0 - std::exp(-term);
                    }
                }

                fbl(1)   = 0.0e0;
                eomui(1) = 0.0e0;
                for (in = 2; in <= inmax; in++) {
                    fbl(in)   = vortn(in) * disn(in) * damp(in);
                    amixl     = vk * disn(in) * damp(in);
                    eomui(in) = (double)reue * xmi * rhon(in) * vortn(in)
                                * amixl * amixl;
                }

                { FortranArray1DRef<double> utot_ref = utot;
                  utmax = q8smax_ns::q8smax(inmax, utot_ref);
                  utmin = q8smin_ns::q8smin(inmax, utot_ref); }
                if ((float)bck(j,i,1) == 1.f) utmin = utot(1);

                // locate max value and location of fbl
                fblmax = 1.0e-10;
                infmax = in1;
                for (in = in1; in <= in2; in++) {
                    if ((float)fbl(in) > (float)fblmax) {
                        fblmax = fbl(in);
                        infmax = in;
                    }
                    // first maximum (degani-schiff)
                    if (ideg[3-1] == 1 && (float)fbl(in) < 0.9e0*(float)fblmax) goto label9125;
                }
                label9125:;

                // curve fit of f near infmax to improve fmax location
                if (ifit > 0) {
                    if (infmax > in1 && infmax < in2) {
                        dfm = fbl(infmax) - fbl(infmax-1);
                        dfp = fbl(infmax) - fbl(infmax+1);
                        if ((float)fbl(infmax-1) < (float)fbl(infmax+1)) {
                            if ((float)dfm > 0.f) {
                                ymax = disn(infmax) + 0.5*(disn(infmax+1) - disn(infmax))*(1.0-dfp/dfm);
                            } else {
                                ymax = disn(infmax);
                            }
                        } else {
                            if ((float)dfp > 0.f) {
                                ymax = disn(infmax) - 0.5*(disn(infmax) - disn(infmax-1))*(1.0-dfm/dfp);
                            } else {
                                ymax = disn(infmax);
                            }
                        }
                    }
                } else {
                    ymax = disn(infmax);
                }

                { double tmp1e20 = 1.0e-20;
                  ymax = ccomplex_ns::ccmaxcr(ymax, tmp1e20); }
                udif   = utmax - utmin;
                fwake  = ymax * fblmax;
                fwake2 = cwk * udif * udif * ymax / fblmax;
                { fwake = std::min(fwake2, fwake); }

                for (in = 2; in <= inmax; in++) {
                    eomuo(in) = 1.0e0 / (1.0e0 + 5.5e0 * std::pow(ckleb*disn(in)/ymax, 6));
                    eomuo(in) = ckout * eomuo(in) * rhon(in) * fwake;
                }
                eomuo(1) = eomuo(2);

                for (in = 2; in <= in2; in++) {
                    if ((float)eomui(in) > (float)eomuo(in)) goto label9127;
                }
                label9127:
                inswtch = in - 1;
                for (in = inswtch; in <= inmax; in++) {
                    eomui(in) = eomuo(in);
                }

                // fill eomu array corresponding to k=0 wall
                for (k = 1; k <= kloop1; k++) {
                    vist3d(j,k,i) = eomui(k+1);
                }
                if (kloop1 < kdim1) {
                    for (k = kloop; k <= kdim1; k++) {
                        vist3d(j,k,i) = 0.0;
                    }
                }


                // iprint output for k=0 wall
                if (iprint >= 1) {
                    if (j == (j/ilfreq1)*ilfreq1 && i == (i/ilfreq2)*ilfreq2) {
                        // length scale parameters
                        if (ihead == 0) {
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                        }
                        ihead = ihead + 1;
                        { eomax = eomui(2);
                          for (int _in=2; _in<=inmax; _in++) if(eomui(_in)>eomax) eomax=eomui(_in); }
                        utmax = utmax / (double)xmach;
                        utmin = utmin / (double)xmach;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "%9d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%6d%6d%4d%4d%4d",
                            j, i, (float)yplus, (float)utmax,
                            (float)utmin, (float)fblmax, (float)ymax, (float)eomax,
                            inswtch, infmax, in1, in2, inmax);
                    }
                }

                // profile data
                if (iprint > 1) {
                    if (j == (j/ipfreq1)*ipfreq1 && i == (i/ipfreq2)*ipfreq2) {
                        ihead = 0;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "   j   i   k      disn      utot      eomu       fbl      amun    eomu-o      vort");
                        for (in = 1; in <= inmax; in++) {
                            utw = utot(in) / (double)xmach;
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                "%5d%5d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e",
                                j, i, in, (float)disn(in), (float)utw,
                                (float)eomui(in), (float)fbl(in),
                                (float)amun(in), (float)eomuo(in), (float)vortn(in));
                        }
                    }
                }

                label1000:;
            } // j loop
            } // i loop

        } // end if (ibckmin==1 || ibckmax==0)


        // evaluate turbulent viscosity along normal to k=kdim wall
        if (ibckmax == 1) {

            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " computing turb viscosity using Baldwin-Lomax, K=kdim, block=%4d", nblt);
            }
            if (isklton > 0 && ideg[3-1] == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned on");
            }
            if (isklton > 0 && ideg[3-1] != 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned off");
            }

            if (inmx < kdim) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "error: inmx, kdim=%5d%5d", inmx, kdim);
                { int m1 = -1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
            }

            ihead = 0;
            // do 10 i=1,idim1; do 10 j=1,jdim1
            for (i = 1; i <= idim1; i++) {
            for (j = 1; j <= jdim1; j++) {

                // bound turbulent region
                kloop  = (int)(0.80 * kdim);
                kloop1 = kloop - 1;
                inmaxt  = kloop;
                inmaxt1 = inmaxt - 1;

                kloop  = kdim1 - kloop + 1;
                kloop  = std::max(1, kloop);
                kloop1 = kloop + 1;
                inmax  = kloop;
                inmax1 = inmax + 1;

                // bound search region for fmax
                in1 = (int)(inmaxt1 * 0.20) + 1;
                if (iwf[3-1] > 0) {
                    in1 = 3;
                }
                in2 = (int)(inmaxt1 * 0.80) + 1;
                in1 = kdim1 - in1 + 1;
                in2 = kdim1 - in2 + 1;

                // Chimera scheme modification
                if (iover == 1) {
                    if ((float)blank(j,in1,i) == 0.f) {
                        for (k = 1; k <= kdim1; k++) {
                            vist3d(j,k,i) = 0.0e0;
                        }
                        if (iprint >= 1) {
                            if (j == j/ilfreq1*ilfreq1 && i == i/ilfreq2*ilfreq2) {
                                if (ihead == 0) {
                                    nou(3) = std::min(nou(3)+1, ibufdim);
                                    std::snprintf(bou(nou(3),3), 120,
                                        " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                                }
                                nou(3) = std::min(nou(3)+1, ibufdim);
                                std::snprintf(bou(nou(3),3), 120,
                                    "%9d%5d     skipping eddy viscosity evaluation at this station due to hole at-near wall", j, i);
                                ihead = ihead + 1;
                            }
                        }
                        goto label10;
                    }
                    in21 = in2;
                    for (in = in1; in >= in2; in--) {
                        if ((float)blank(j,in,i) == 0.f) goto label24;
                        in21 = in;
                    }
                    label24:;
                    in2 = in21;
                }

                // viscosity thru sutherland relation at cell centers
                for (k = kdim1; k >= kloop1; k--) {
                    t5       = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
                    t6       = std::sqrt(t5);
                    amun(k)  = c2bp * t5 * t6 / (c2b + t5);
                }
                for (k = kdim1; k >= kloop1; k--) {
                    rhon(k)  = q(j,k,i,1);
                    vortn(k) = vor(j,k,i);
                    { double tmp = snkm(j,k,i); disn(k) = std::abs(tmp); }
                    utot(k)  = std::sqrt( q(j,k,i,2)*q(j,k,i,2)
                                        + q(j,k,i,3)*q(j,k,i,3)
                                        + q(j,k,i,4)*q(j,k,i,4) );
                }

                rhon(kdim)  = rhon(kdim1);
                amun(kdim)  = amun(kdim1);
                disn(kdim)  = 0.0;
                utot(kdim)  = utot(kdim1);
                vortn(kdim) = vortn(kdim1);

                if ((float)bck(j,i,2) == 1.f) {
                    t5 = (double)gamma * qk0(j,i,5,3) / qk0(j,i,1,3);
                    t6 = std::sqrt(t5);
                    amun(kdim)  = c2bp * t5 * t6 / (c2b + t5);
                    rhon(kdim)  = qk0(j,i,1,3);
                    utot(kdim)  = std::sqrt( qk0(j,i,2,3)*qk0(j,i,2,3)
                                           + qk0(j,i,3,3)*qk0(j,i,3,3)
                                           + qk0(j,i,4,3)*qk0(j,i,4,3) );
                    vortn(kdim) = (utot(kdim1) - utot(kdim)) / disn(kdim1);
                }

                // convenient groupings and initial values at k=kdim
                { double tmp = rhon(kdim)*vortn(kdim)*(double)reue*xmi/amun(kdim);
                  yplusoy = std::abs(tmp); }
                yplusoy = std::sqrt(yplusoy) * aplusi;
                yplus   = yplusoy * disn(kdim1) * aplus * 2.0;

                if ((float)bck(j,i,2) == 0.f) {
                    // wake damping
                    term = 1.0 - std::exp(-50.0e0);
                    for (in = kdim1; in >= inmax1; in--) {
                        damp(in) = term;
                    }
                } else {
                    // wall damping
                    for (in = kdim1; in >= inmax1; in--) {
                        if (iwall == 0) {
                            term = yplusoy * disn(in) * std::sqrt(rhon(in)/rhon(kdim))
                                   * (amun(kdim)/amun(in));
                        } else {
                            term = yplusoy * disn(in);
                        }
                        damp(in) = 1.0e0 - std::exp(-term);
                    }
                }

                fbl(kdim)   = 0.0e0;
                eomui(kdim) = 0.0e0;
                for (in = kdim1; in >= inmax1; in--) {
                    fbl(in)   = vortn(in) * disn(in) * damp(in);
                    amixl     = vk * disn(in) * damp(in);
                    eomui(in) = (double)reue * xmi * rhon(in) * vortn(in)
                                * amixl * amixl;
                }

                { utmax = utot(inmax1); utmin = utot(inmax1);
                  for (int _in=inmax1; _in<=inmax1+inmaxt1-1; _in++) {
                      if(utot(_in)>utmax) utmax=utot(_in);
                      if(utot(_in)<utmin) utmin=utot(_in); } }
                if ((float)bck(j,i,2) == 1.f) utmin = utot(kdim);


                // locate max value and location of fbl (k=kdim wall)
                fblmax = 1.0e-10;
                infmax = in1 + 1;
                for (in = in1+1; in >= in2+1; in--) {
                    if ((float)fbl(in) > (float)fblmax) {
                        fblmax = fbl(in);
                        infmax = in;
                    }
                    // first maximum (degani-schiff)
                    if (ideg[3-1] == 1 && (float)fbl(in) < 0.9e0*(float)fblmax) goto label125;
                }
                label125:;

                // curve fit of f near infmax to improve fmax location
                if (ifit > 0) {
                    if (infmax < in1 && infmax > in2) {
                        dfm = fbl(infmax) - fbl(infmax+1);
                        dfp = fbl(infmax) - fbl(infmax-1);
                        if ((float)fbl(infmax+1) < (float)fbl(infmax-1)) {
                            if ((float)dfm > 0.f) {
                                ymax = disn(infmax) + 0.5*(disn(infmax-1) - disn(infmax))*(1.0-dfp/dfm);
                            } else {
                                ymax = disn(infmax);
                            }
                        } else {
                            if ((float)dfp > 0.f) {
                                ymax = disn(infmax) - 0.5*(disn(infmax) - disn(infmax+1))*(1.0-dfm/dfp);
                            } else {
                                ymax = disn(infmax);
                            }
                        }
                    }
                } else {
                    ymax = disn(infmax);
                }

                { double tmp1e20 = 1.0e-20;
                  ymax = ccomplex_ns::ccmaxcr(ymax, tmp1e20); }
                udif   = utmax - utmin;
                fwake  = ymax * fblmax;
                fwake2 = cwk * udif * udif * ymax / fblmax;
                { fwake = std::min(fwake2, fwake); }

                for (in = kdim1; in >= inmax1; in--) {
                    eomuo(in) = 1.0e0 / (1.0e0 + 5.5e0 * std::pow(ckleb*disn(in)/ymax, 6));
                    eomuo(in) = ckout * eomuo(in) * rhon(in) * fwake;
                }
                eomuo(kdim) = eomuo(kdim1);

                for (in = kdim1; in >= in2; in--) {
                    if ((float)eomui(in) > (float)eomuo(in)) goto label127;
                }
                label127:
                inswtch = in + 1;
                for (in = inswtch; in >= inmax1; in--) {
                    eomui(in) = eomuo(in);
                }

                // fill eomu array corresponding to k=kdim wall
                if (ibckmin == 1) {
                    // for k=1 & k=kdim walls
                    for (k = kdim1; k >= kloop1; k--) {
                        s1 = snk0(j,k,i) * snk0(j,k,i);
                        s2 = snkm(j,k,i) * snkm(j,k,i);
                        vist3d(j,k,i) = ( s1*eomui(k) + s2*vist3d(j,k,i) ) / ( s1 + s2 );
                    }
                } else {
                    // for k=kdim wall only
                    for (k = kdim1; k >= kloop1; k--) {
                        vist3d(j,k,i) = eomui(k);
                    }
                    if (kloop1 > 1) {
                        for (k = kloop; k >= 1; k--) {
                            vist3d(j,k,i) = 0.0;
                        }
                    }
                }

                // iprint output for k=kdim wall
                if (iprint >= 1) {
                    if (j == (j/ilfreq1)*ilfreq1 && i == (i/ilfreq2)*ilfreq2) {
                        if (ihead == 0) {
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                        }
                        ihead = ihead + 1;
                        { eomax = eomui(kdim1);
                          for (int _in=kdim1; _in>=kdim1-inmax1+1; _in--) if(eomui(_in)>eomax) eomax=eomui(_in); }
                        utmax = utmax / (double)xmach;
                        utmin = utmin / (double)xmach;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "%9d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%6d%6d%4d%4d%4d",
                            j, i, (float)yplus, (float)utmax,
                            (float)utmin, (float)fblmax, (float)ymax, (float)eomax,
                            inswtch, infmax, in1, in2, inmax);
                    }
                }

                // profile data
                if (iprint > 1) {
                    if (j == (j/ipfreq1)*ipfreq1 && i == (i/ipfreq2)*ipfreq2) {
                        ihead = 0;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "   j   i   k      disn      utot      eomu       fbl      amun    eomu-o      vort");
                        for (in = kdim1; in >= inmax1; in--) {
                            utw = utot(in) / (double)xmach;
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                "%5d%5d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e",
                                j, i, in, (float)disn(in), (float)utw,
                                (float)eomui(in), (float)fbl(in),
                                (float)amun(in), (float)eomuo(in), (float)vortn(in));
                        }
                    }
                }

                label10:;
            } // j loop
            } // i loop

        } // end if (ibckmax==1)

        // Lines 773-794: laminar region zeroing for ivisc(3)>1
        if (isklton > 0) {
            if (ilamlo == 0 || jlamlo == 0 || klamlo == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   block %4d in B-L turb model has no laminar regions", nblt);
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   block %4d in B-L turb model - laminar region is:", nblt);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   i=%5d to%5d, j=%5d to%5d, k=%5d to%5d",
                    ilamlo, ilamhi, jlamlo, jlamhi, klamlo, klamhi);
            }
        }
        for (i = ilamlo; i <= ilamhi-1; i++) {
        for (k = klamlo; k <= klamhi-1; k++) {
        for (j = jlamlo; j <= jlamhi-1; j++) {
            vist3d(j,k,i) = 0.0;
        }}}

    } // end if (ivisc[3-1] > 1)


    // Lines 796-1079: ivisc(2)>1: j-wall turbulence (j=0 face)
    if (ivisc[2-1] > 1) {
        // evaluate turbulent viscosity along normal to j=0 wall
        // defaults to using min face for distance if neither min nor max contains wall
        if (ibcjmin == 1 || ibcjmax == 0) {

            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " computing turb viscosity using Baldwin-Lomax, J=1   , block=%4d", nblt);
            }
            if (isklton > 0 && ideg[2-1] == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned on");
            }
            if (isklton > 0 && ideg[2-1] != 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned off");
            }

            if (inmx < jdim) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "error: inmx, jdim=%5d%5d", inmx, jdim);
                { int m1 = -1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
            }

            ihead = 0;
            // do 7000 i=1,idim1; do 7000 k=1,kdim1
            for (i = 1; i <= idim1; i++) {
            for (k = 1; k <= kdim1; k++) {

                // bound turbulent region
                int jloop_loc  = (int)(0.80 * jdim);
                int jloop1_loc = jloop_loc - 1;
                inmax  = jloop_loc;
                inmax1 = inmax - 1;

                // bound search region for fmax
                in1 = (int)(inmax1 * 0.20) + 1;
                if (iwf[2-1] > 0) { in1 = 3; }
                in2 = (int)(inmax1 * 0.80) + 1;

                // Chimera scheme modification
                if (iover == 1) {
                    if ((float)blank(in1,k,i) == 0.f) {
                        for (j = 1; j <= jdim1; j++) { eoms(j,k,i) = 0.0e0; }
                        if (iprint >= 1) {
                            if (k == k/ilfreq1*ilfreq1 && i == i/ilfreq2*ilfreq2) {
                                if (ihead == 0) {
                                    nou(3) = std::min(nou(3)+1, ibufdim);
                                    std::snprintf(bou(nou(3),3), 120,
                                        " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                                }
                                nou(3) = std::min(nou(3)+1, ibufdim);
                                std::snprintf(bou(nou(3),3), 120,
                                    "%9d%5d     skipping eddy viscosity evaluation at this station due to hole at-near wall", k, i);
                                ihead = ihead + 1;
                            }
                        }
                        goto label7000;
                    }
                    in21 = in2;
                    for (in = in1; in <= in2; in++) {
                        if ((float)blank(in,k,i) == 0.f) goto label1624;
                        in21 = in;
                    }
                    label1624:;
                    in2 = in21;
                }

                // viscosity thru sutherland relation at cell centers
                for (j = 1; j <= jloop1_loc; j++) {
                    t5 = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
                    t6 = std::sqrt(t5);
                    amun(j+1) = c2bp * t5 * t6 / (c2b + t5);
                }
                for (j = 1; j <= jloop1_loc; j++) {
                    rhon(j+1)  = q(j,k,i,1);
                    vortn(j+1) = vor(j,k,i);
                    { double tmp = snj0(j,k,i); disn(j+1) = std::abs(tmp); }
                    utot(j+1)  = std::sqrt( q(j,k,i,2)*q(j,k,i,2)
                                          + q(j,k,i,3)*q(j,k,i,3)
                                          + q(j,k,i,4)*q(j,k,i,4) );
                }
                rhon(1) = rhon(2); amun(1) = amun(2);
                disn(1) = 0.0; utot(1) = utot(2); vortn(1) = vortn(2);

                if ((float)bcj(k,i,1) == 1.f) {
                    t5 = (double)gamma * qj0(k,i,5,1) / qj0(k,i,1,1);
                    t6 = std::sqrt(t5);
                    amun(1)  = c2bp * t5 * t6 / (c2b + t5);
                    rhon(1)  = qj0(k,i,1,1);
                    utot(1)  = std::sqrt( qj0(k,i,2,1)*qj0(k,i,2,1)
                                        + qj0(k,i,3,1)*qj0(k,i,3,1)
                                        + qj0(k,i,4,1)*qj0(k,i,4,1) );
                    vortn(1) = (utot(2) - utot(1)) / disn(2);
                }

                { double tmp = rhon(1)*vortn(1)*(double)reue*xmi/amun(1);
                  yplusoy = std::abs(tmp); }
                yplusoy = std::sqrt(yplusoy) * aplusi;
                yplus   = yplusoy * disn(2) * aplus * 2.0;


                if ((float)bcj(k,i,1) == 0.f) {
                    term = 1.0 - std::exp(-50.0e0);
                    for (in = 2; in <= inmax; in++) { damp(in) = term; }
                } else {
                    for (in = 2; in <= inmax; in++) {
                        if (iwall == 0) {
                            term = yplusoy * disn(in) * std::sqrt(rhon(in)/rhon(1)) * (amun(1)/amun(in));
                        } else {
                            term = yplusoy * disn(in);
                        }
                        damp(in) = 1.0e0 - std::exp(-term);
                    }
                }

                fbl(1) = 0.0e0; eomui(1) = 0.0e0;
                for (in = 2; in <= inmax; in++) {
                    fbl(in)   = vortn(in) * disn(in) * damp(in);
                    amixl     = vk * disn(in) * damp(in);
                    eomui(in) = (double)reue * xmi * rhon(in) * vortn(in) * amixl * amixl;
                }

                { FortranArray1DRef<double> utot_ref = utot;
                  utmax = q8smax_ns::q8smax(inmax, utot_ref);
                  utmin = q8smin_ns::q8smin(inmax, utot_ref); }
                if ((float)bcj(k,i,1) == 1.f) utmin = utot(1);

                fblmax = 1.0e-10; infmax = in1;
                for (in = in1; in <= in2; in++) {
                    if ((float)fbl(in) > (float)fblmax) { fblmax = fbl(in); infmax = in; }
                    if (ideg[2-1] == 1 && (float)fbl(in) < 0.9e0*(float)fblmax) goto label6125;
                }
                label6125:;

                if (ifit > 0) {
                    if (infmax > in1 && infmax < in2) {
                        dfm = fbl(infmax) - fbl(infmax-1);
                        dfp = fbl(infmax) - fbl(infmax+1);
                        if ((float)fbl(infmax-1) < (float)fbl(infmax+1)) {
                            if ((float)dfm > 0.f) {
                                ymax = disn(infmax) + 0.5*(disn(infmax+1)-disn(infmax))*(1.0-dfp/dfm);
                            } else { ymax = disn(infmax); }
                        } else {
                            if ((float)dfp > 0.f) {
                                ymax = disn(infmax) - 0.5*(disn(infmax)-disn(infmax-1))*(1.0-dfm/dfp);
                            } else { ymax = disn(infmax); }
                        }
                    }
                } else { ymax = disn(infmax); }

                { double tmp1e20 = 1.0e-20; ymax = ccomplex_ns::ccmaxcr(ymax, tmp1e20); }
                udif = utmax - utmin;
                fwake = ymax * fblmax;
                fwake2 = cwk * udif * udif * ymax / fblmax;
                { fwake = std::min(fwake2, fwake); }

                for (in = 2; in <= inmax; in++) {
                    eomuo(in) = 1.0e0 / (1.0e0 + 5.5e0 * std::pow(ckleb*disn(in)/ymax, 6));
                    eomuo(in) = ckout * eomuo(in) * rhon(in) * fwake;
                }
                eomuo(1) = eomuo(2);

                for (in = 2; in <= in2; in++) {
                    if ((float)eomui(in) > (float)eomuo(in)) goto label6127;
                }
                label6127:
                inswtch = in - 1;
                for (in = inswtch; in <= inmax; in++) { eomui(in) = eomuo(in); }

                // fill eomu array corresponding to j=0 wall
                for (j = 1; j <= jloop1_loc; j++) { eoms(j,k,i) = eomui(j+1); }
                if (jloop1_loc < jdim1) {
                    for (j = jloop_loc; j <= jdim1; j++) { eoms(j,k,i) = 0.0; }
                }

                if (iprint >= 1) {
                    if (k == (k/ilfreq1)*ilfreq1 && i == (i/ilfreq2)*ilfreq2) {
                        if (ihead == 0) {
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                " j=0   k   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                        }
                        ihead = ihead + 1;
                        { eomax = eomui(2);
                          for (int _in=2; _in<=inmax; _in++) if(eomui(_in)>eomax) eomax=eomui(_in); }
                        utmax = utmax / (double)xmach; utmin = utmin / (double)xmach;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "%9d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%6d%6d%4d%4d%4d",
                            k, i, (float)yplus, (float)utmax,
                            (float)utmin, (float)fblmax, (float)ymax, (float)eomax,
                            inswtch, infmax, in1, in2, inmax);
                    }
                }
                if (iprint > 1) {
                    if (k == (k/ipfreq1)*ipfreq1 && i == (i/ipfreq2)*ipfreq2) {
                        ihead = 0;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "   k   i   j      disn      utot      eomu       fbl      amun    eomu-o      vort");
                        for (in = 1; in <= inmax; in++) {
                            utw = utot(in) / (double)xmach;
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                "%5d%5d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e",
                                k, i, in, (float)disn(in), (float)utw,
                                (float)eomui(in), (float)fbl(in),
                                (float)amun(in), (float)eomuo(in), (float)vortn(in));
                        }
                    }
                }

                label7000:;
            } // k loop
            } // i loop

        } // end if (ibcjmin==1 || ibcjmax==0)


        // evaluate turbulent viscosity along normal to j=jdim wall
        if (ibcjmax == 1) {

            if (isklton > 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    " computing turb viscosity using Baldwin-Lomax, J=jdim, block=%4d", nblt);
            }
            if (isklton > 0 && ideg[2-1] == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned on");
            }
            if (isklton > 0 && ideg[2-1] != 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned off");
            }

            if (inmx < jdim) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "error: inmx, jdim=%5d%5d", inmx, jdim);
                { int m1 = -1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
            }

            ihead = 0;
            // do 7610 i=1,idim1; do 7610 k=1,kdim1
            for (i = 1; i <= idim1; i++) {
            for (k = 1; k <= kdim1; k++) {

                int jloop_loc  = (int)(0.80 * jdim);
                int jloop1_loc = jloop_loc - 1;
                int inmaxt_loc  = jloop_loc;
                int inmaxt1_loc = inmaxt_loc - 1;

                jloop_loc  = jdim1 - jloop_loc + 1;
                jloop_loc  = std::max(1, jloop_loc);
                jloop1_loc = jloop_loc + 1;
                inmax  = jloop_loc;
                inmax1 = inmax + 1;

                in1 = (int)(inmaxt1_loc * 0.20) + 1;
                if (iwf[2-1] > 0) { in1 = 3; }
                in2 = (int)(inmaxt1_loc * 0.80) + 1;
                in1 = jdim1 - in1 + 1;
                in2 = jdim1 - in2 + 1;

                // Chimera scheme modification
                if (iover == 1) {
                    if ((float)blank(in1,k,i) == 0.f) {
                        for (j = 1; j <= jdim1; j++) { eoms(j,k,i) = 0.0e0; }
                        if (iprint >= 1) {
                            if (k == k/ilfreq1*ilfreq1 && i == i/ilfreq2*ilfreq2) {
                                if (ihead == 0) {
                                    nou(3) = std::min(nou(3)+1, ibufdim);
                                    std::snprintf(bou(nou(3),3), 120,
                                        " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                                }
                                nou(3) = std::min(nou(3)+1, ibufdim);
                                std::snprintf(bou(nou(3),3), 120,
                                    "%9d%5d     skipping eddy viscosity evaluation at this station due to hole at-near wall", k, i);
                                ihead = ihead + 1;
                            }
                        }
                        goto label7610;
                    }
                    in21 = in2;
                    for (in = in1; in >= in2; in--) {
                        if ((float)blank(in,k,i) == 0.f) goto label624;
                        in21 = in;
                    }
                    label624:;
                    in2 = in21;
                }

                // viscosity thru sutherland relation at cell centers
                for (j = jdim1; j >= jloop1_loc; j--) {
                    t5 = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
                    t6 = std::sqrt(t5);
                    amun(j) = c2bp * t5 * t6 / (c2b + t5);
                }
                for (j = jdim1; j >= jloop1_loc; j--) {
                    rhon(j)  = q(j,k,i,1);
                    vortn(j) = vor(j,k,i);
                    { double tmp = snjm(j,k,i); disn(j) = std::abs(tmp); }
                    utot(j)  = std::sqrt( q(j,k,i,2)*q(j,k,i,2)
                                        + q(j,k,i,3)*q(j,k,i,3)
                                        + q(j,k,i,4)*q(j,k,i,4) );
                }
                rhon(jdim) = rhon(jdim1); amun(jdim) = amun(jdim1);
                disn(jdim) = 0.0; utot(jdim) = utot(jdim1); vortn(jdim) = vortn(jdim1);

                if ((float)bcj(k,i,2) == 1.f) {
                    t5 = (double)gamma * qj0(k,i,5,3) / qj0(k,i,1,3);
                    t6 = std::sqrt(t5);
                    amun(jdim)  = c2bp * t5 * t6 / (c2b + t5);
                    rhon(jdim)  = qj0(k,i,1,3);
                    utot(jdim)  = std::sqrt( qj0(k,i,2,3)*qj0(k,i,2,3)
                                           + qj0(k,i,3,3)*qj0(k,i,3,3)
                                           + qj0(k,i,4,3)*qj0(k,i,4,3) );
                    vortn(jdim) = (utot(jdim1) - utot(jdim)) / disn(jdim1);
                }

                { double tmp = rhon(jdim)*vortn(jdim)*(double)reue*xmi/amun(jdim);
                  yplusoy = std::abs(tmp); }
                yplusoy = std::sqrt(yplusoy) * aplusi;
                yplus   = yplusoy * disn(jdim1) * aplus * 2.0;


                if ((float)bcj(k,i,2) == 0.f) {
                    term = 1.0 - std::exp(-50.0e0);
                    for (in = jdim1; in >= inmax1; in--) { damp(in) = term; }
                } else {
                    for (in = jdim1; in >= inmax1; in--) {
                        if (iwall == 0) {
                            term = yplusoy * disn(in) * std::sqrt(rhon(in)/rhon(jdim)) * (amun(jdim)/amun(in));
                        } else { term = yplusoy * disn(in); }
                        damp(in) = 1.0e0 - std::exp(-term);
                    }
                }

                fbl(jdim) = 0.0e0; eomui(jdim) = 0.0e0;
                for (in = jdim1; in >= inmax1; in--) {
                    fbl(in)   = vortn(in) * disn(in) * damp(in);
                    amixl     = vk * disn(in) * damp(in);
                    eomui(in) = (double)reue * xmi * rhon(in) * vortn(in) * amixl * amixl;
                }

                { utmax = utot(inmax1); utmin = utot(inmax1);
                  for (int _in=inmax1; _in<=inmax1+inmaxt1_loc-1; _in++) {
                      if(utot(_in)>utmax) utmax=utot(_in);
                      if(utot(_in)<utmin) utmin=utot(_in); } }
                if ((float)bcj(k,i,2) == 1.f) utmin = utot(jdim);

                fblmax = 1.0e-10; infmax = in1 + 1;
                for (in = in1+1; in >= in2+1; in--) {
                    if ((float)fbl(in) > (float)fblmax) { fblmax = fbl(in); infmax = in; }
                    if (ideg[2-1] == 1 && (float)fbl(in) < 0.9e0*(float)fblmax) goto label925;
                }
                label925:;

                if (ifit > 0) {
                    if (infmax < in1 && infmax > in2) {
                        dfm = fbl(infmax) - fbl(infmax+1);
                        dfp = fbl(infmax) - fbl(infmax-1);
                        if ((float)fbl(infmax+1) < (float)fbl(infmax-1)) {
                            if ((float)dfm > 0.f) {
                                ymax = disn(infmax) + 0.5*(disn(infmax-1)-disn(infmax))*(1.0-dfp/dfm);
                            } else { ymax = disn(infmax); }
                        } else {
                            if ((float)dfp > 0.f) {
                                ymax = disn(infmax) - 0.5*(disn(infmax)-disn(infmax+1))*(1.0-dfm/dfp);
                            } else { ymax = disn(infmax); }
                        }
                    }
                } else { ymax = disn(infmax); }

                { double tmp1e20 = 1.0e-20; ymax = ccomplex_ns::ccmaxcr(ymax, tmp1e20); }
                udif = utmax - utmin;
                fwake = ymax * fblmax;
                fwake2 = cwk * udif * udif * ymax / fblmax;
                { fwake = std::min(fwake2, fwake); }

                for (in = jdim1; in >= inmax1; in--) {
                    eomuo(in) = 1.0e0 / (1.0e0 + 5.5e0 * std::pow(ckleb*disn(in)/ymax, 6));
                    eomuo(in) = ckout * eomuo(in) * rhon(in) * fwake;
                }
                eomuo(jdim) = eomuo(jdim1);

                for (in = jdim1; in >= in2; in--) {
                    if ((float)eomui(in) > (float)eomuo(in)) goto label137;
                }
                label137:
                inswtch = in + 1;
                for (in = inswtch; in >= inmax1; in--) { eomui(in) = eomuo(in); }

                // fill eomu array corresponding to j=jdim wall
                if (ibcjmin == 1) {
                    // for j=1 & j=jdim walls
                    for (j = jdim1; j >= jloop1_loc; j--) {
                        s1 = snj0(j,k,i) * snj0(j,k,i);
                        s2 = snjm(j,k,i) * snjm(j,k,i);
                        eoms(j,k,i) = ( s1*eomui(j) + s2*eoms(j,k,i) ) / ( s1 + s2 );
                    }
                } else {
                    // for j=jdim wall only
                    for (j = jdim1; j >= jloop1_loc; j--) { eoms(j,k,i) = eomui(j); }
                    if (jloop1_loc > 1) {
                        for (j = jloop_loc; j >= 1; j--) { eoms(j,k,i) = 0.0; }
                    }
                }

                if (iprint >= 1) {
                    if (k == (k/ilfreq1)*ilfreq1 && i == (i/ilfreq2)*ilfreq2) {
                        if (ihead == 0) {
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                        }
                        ihead = ihead + 1;
                        { eomax = eomui(jdim1);
                          for (int _in=jdim1; _in>=jdim1-inmax1+1; _in--) if(eomui(_in)>eomax) eomax=eomui(_in); }
                        utmax = utmax / (double)xmach; utmin = utmin / (double)xmach;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "%9d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%6d%6d%4d%4d%4d",
                            k, i, (float)yplus, (float)utmax,
                            (float)utmin, (float)fblmax, (float)ymax, (float)eomax,
                            inswtch, infmax, in1, in2, inmax);
                    }
                }
                if (iprint > 1) {
                    if (k == (k/ipfreq1)*ipfreq1 && i == (i/ipfreq2)*ipfreq2) {
                        ihead = 0;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "   j   i   k      disn      utot      eomu       fbl      amun    eomu-o      vort");
                        for (in = jdim1; in >= inmax1; in--) {
                            utw = utot(in) / (double)xmach;
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                "%5d%5d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e",
                                k, i, in, (float)disn(in), (float)utw,
                                (float)eomui(in), (float)fbl(in),
                                (float)amun(in), (float)eomuo(in), (float)vortn(in));
                        }
                    }
                }

                label7610:;
            } // k loop
            } // i loop

        } // end if (ibcjmax==1)


        // Lines 1379-1463: form composite eddy-viscosity from k & j walls
        if (ivisc[3-1] > 1) {
            for (i = 1; i <= idim1; i++) {
            for (k = 1; k <= kdim1; k++) {
            for (j = 1; j <= jdim1; j++) {
                if (ibckmin == 0 && ibckmax == 1) {
                    s1 = snkm(j,k,i) * snkm(j,k,i);
                } else if (ibckmin == 1 && ibckmax == 1) {
                    { double a = snk0(j,k,i), b = snkm(j,k,i); s1 = std::min(a,b); }
                    s1 = s1 * s1;
                } else {
                    s1 = snk0(j,k,i) * snk0(j,k,i);
                }
                if (ibcjmin == 0 && ibcjmax == 1) {
                    s2 = snjm(j,k,i) * snjm(j,k,i);
                } else if (ibcjmin == 1 && ibcjmax == 1) {
                    { double a = snj0(j,k,i), b = snjm(j,k,i); s2 = std::min(a,b); }
                    s2 = s2 * s2;
                } else {
                    s2 = snj0(j,k,i) * snj0(j,k,i);
                }
                vist3d(j,k,i) = ( vist3d(j,k,i)*s2 + eoms(j,k,i)*s1 ) / ( s1 + s2 );
            }}}

            // smooth interior values
            nsmurf = 3;
            if (nsmurf == 0) return;

            for (nnn = 1; nnn <= nsmurf; nnn++) {
                for (i = 1; i <= idim1; i++) {
                for (k = 1; k <= kdim1; k++) {
                for (j = 1; j <= jdim1; j++) {
                    eoms(j,k,i) = vist3d(j,k,i);
                }}}
                for (i = 2; i <= idim1-1; i++) {
                for (k = 2; k <= kdim1-1; k++) {
                for (j = 2; j <= jdim1-1; j++) {
                    vist3d(j,k,i) = 0.125e0 * ( eoms(j+1,k-1,i-1) + eoms(j-1,k-1,i-1)
                                               + eoms(j+1,k+1,i-1) + eoms(j-1,k+1,i-1)
                                               + eoms(j+1,k-1,i+1) + eoms(j-1,k-1,i+1)
                                               + eoms(j+1,k+1,i+1) + eoms(j-1,k+1,i+1) );
                }}}
            }

        } else {
            // single turbulent wall - install vist3d array as eoms array
            for (i = 1; i <= idim1; i++) {
            for (k = 1; k <= kdim1; k++) {
            for (j = 1; j <= jdim1; j++) {
                vist3d(j,k,i) = eoms(j,k,i);
            }}}
        }

        // laminar region zeroing for ivisc(2)>1
        if (isklton > 0) {
            if (ilamlo == 0 || jlamlo == 0 || klamlo == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   block %4d in B-L turb model has no laminar regions", nblt);
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   block %4d in B-L turb model - laminar region is:", nblt);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120,
                    "   i=%5d to%5d, j=%5d to%5d, k=%5d to%5d",
                    ilamlo, ilamhi, jlamlo, jlamhi, klamlo, klamhi);
            }
        }
        for (i = ilamlo; i <= ilamhi-1; i++) {
        for (k = klamlo; k <= klamhi-1; k++) {
        for (j = jlamlo; j <= jlamhi-1; j++) {
            vist3d(j,k,i) = 0.0;
        }}}
        return;

    } // end if (ivisc[2-1] > 1)


    // Lines 1465-2132: ivisc(1)>1: i-wall turbulence
    if (ivisc[1-1] <= 1) return;

    // evaluate turbulent viscosity along normal to i=0 wall
    // defaults to using min face for distance if neither min nor max contains wall
    if (ibcimin == 1 || ibcimax == 0) {

        if (isklton > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " computing turb viscosity using Baldwin-Lomax, I=1   , block=%4d", nblt);
        }
        if (isklton > 0 && ideg[1-1] == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned on");
        }
        if (isklton > 0 && ideg[1-1] != 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned off");
        }

        ihead = 0;
        if (inmx < idim) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "error: inmx, idim=%5d%5d", inmx, idim);
            { int m1 = -1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
        }

        // do 2000 k=1,kdim1; do 2000 j=1,jdim1
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {

            iloop  = (int)(0.80 * idim);
            iloop1 = iloop - 1;
            inmax  = iloop;
            inmax1 = inmax - 1;

            in1 = (int)(inmax1 * 0.20) + 1;
            if (iwf[1-1] > 0) { in1 = 3; }
            in2 = (int)(inmax1 * 0.80) + 1;

            // Chimera scheme modification
            if (iover == 1) {
                if ((float)blank(j,k,in1) == 0.f) {
                    for (i = 1; i <= idim1; i++) { eoms(j,k,i) = 0.0e0; }
                    if (iprint >= 1) {
                        if (j == j/ilfreq1*ilfreq1 && k == k/ilfreq2*ilfreq2) {
                            if (ihead == 0) {
                                nou(3) = std::min(nou(3)+1, ibufdim);
                                std::snprintf(bou(nou(3),3), 120,
                                    " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                            }
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                "%9d%5d     skipping eddy viscosity evaluation at this station due to hole at-near wall", j, k);
                            ihead = ihead + 1;
                        }
                    }
                    goto label2000;
                }
                in21 = in2;
                for (in = in1; in <= in2; in++) {
                    if ((float)blank(j,k,in) == 0.f) goto label1724;
                    in21 = in;
                }
                label1724:;
                in2 = in21;
            }

            // viscosity thru sutherland relation at cell centers
            for (i = 1; i <= iloop1; i++) {
                t5 = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
                t6 = std::sqrt(t5);
                amun(i+1) = c2bp * t5 * t6 / (c2b + t5);
            }
            for (i = 1; i <= iloop1; i++) {
                rhon(i+1)  = q(j,k,i,1);
                vortn(i+1) = vor(j,k,i);
                { double tmp = sni0(j,k,i); disn(i+1) = std::abs(tmp); }
                utot(i+1)  = std::sqrt( q(j,k,i,2)*q(j,k,i,2)
                                      + q(j,k,i,3)*q(j,k,i,3)
                                      + q(j,k,i,4)*q(j,k,i,4) );
            }
            rhon(1) = rhon(2); amun(1) = amun(2);
            disn(1) = 0.0; utot(1) = utot(2); vortn(1) = vortn(2);

            if ((float)bci(j,k,1) == 1.f) {
                t5 = (double)gamma * qi0(j,k,5,1) / qi0(j,k,1,1);
                t6 = std::sqrt(t5);
                amun(1)  = c2bp * t5 * t6 / (c2b + t5);
                rhon(1)  = qi0(j,k,1,1);
                utot(1)  = std::sqrt( qi0(j,k,2,1)*qi0(j,k,2,1)
                                    + qi0(j,k,3,1)*qi0(j,k,3,1)
                                    + qi0(j,k,4,1)*qi0(j,k,4,1) );
                vortn(1) = (utot(2) - utot(1)) / disn(2);
            }

            { double tmp = rhon(1)*vortn(1)*(double)reue*xmi/amun(1);
              yplusoy = std::abs(tmp); }
            yplusoy = std::sqrt(yplusoy) * aplusi;
            yplus   = yplusoy * disn(2) * aplus * 2.0;


            if ((float)bci(j,k,1) == 0.f) {
                term = 1.0 - std::exp(-50.0e0);
                for (in = 2; in <= inmax; in++) { damp(in) = term; }
            } else {
                for (in = 2; in <= inmax; in++) {
                    if (iwall == 0) {
                        term = yplusoy * disn(in) * std::sqrt(rhon(in)/rhon(1)) * (amun(1)/amun(in));
                    } else { term = yplusoy * disn(in); }
                    damp(in) = 1.0e0 - std::exp(-term);
                }
            }

            fbl(1) = 0.0e0; eomui(1) = 0.0e0;
            for (in = 2; in <= inmax; in++) {
                fbl(in)   = vortn(in) * disn(in) * damp(in);
                amixl     = vk * disn(in) * damp(in);
                eomui(in) = (double)reue * xmi * rhon(in) * vortn(in) * amixl * amixl;
            }

            { FortranArray1DRef<double> utot_ref = utot;
              utmax = q8smax_ns::q8smax(inmax, utot_ref);
              utmin = q8smin_ns::q8smin(inmax, utot_ref); }
            if ((float)bci(j,k,1) == 1.f) utmin = utot(1);

            fblmax = 1.0e-10; infmax = in1;
            for (in = in1; in <= in2; in++) {
                if ((float)fbl(in) > (float)fblmax) { fblmax = fbl(in); infmax = in; }
                if (ideg[1-1] == 1 && (float)fbl(in) < 0.9e0*(float)fblmax) goto label8125;
            }
            label8125:;

            if (ifit > 0) {
                if (infmax > in1 && infmax < in2) {
                    dfm = fbl(infmax) - fbl(infmax-1);
                    dfp = fbl(infmax) - fbl(infmax+1);
                    if ((float)fbl(infmax-1) < (float)fbl(infmax+1)) {
                        if ((float)dfm > 0.f) {
                            ymax = disn(infmax) + 0.5*(disn(infmax+1)-disn(infmax))*(1.0-dfp/dfm);
                        } else { ymax = disn(infmax); }
                    } else {
                        if ((float)dfp > 0.f) {
                            ymax = disn(infmax) - 0.5*(disn(infmax)-disn(infmax-1))*(1.0-dfm/dfp);
                        } else { ymax = disn(infmax); }
                    }
                }
            } else { ymax = disn(infmax); }

            { double tmp1e20 = 1.0e-20; ymax = ccomplex_ns::ccmaxcr(ymax, tmp1e20); }
            udif = utmax - utmin;
            fwake = ymax * fblmax;
            fwake2 = cwk * udif * udif * ymax / fblmax;
            { fwake = std::min(fwake2, fwake); }

            for (in = 2; in <= inmax; in++) {
                eomuo(in) = 1.0e0 / (1.0e0 + 5.5e0 * std::pow(ckleb*disn(in)/ymax, 6));
                eomuo(in) = ckout * eomuo(in) * rhon(in) * fwake;
            }
            eomuo(1) = eomuo(2);

            for (in = 2; in <= in2; in++) {
                if ((float)eomui(in) > (float)eomuo(in)) goto label8127;
            }
            label8127:
            inswtch = in - 1;
            for (in = inswtch; in <= inmax; in++) { eomui(in) = eomuo(in); }

            // fill eomu array corresponding to i=0 wall
            for (i = 1; i <= iloop1; i++) { eoms(j,k,i) = eomui(i+1); }
            if (iloop1 < idim1) {
                for (i = iloop; i <= idim1; i++) { eoms(j,k,i) = 0.0; }
            }

            if (iprint >= 1) {
                if (j == (j/ilfreq1)*ilfreq1 && k == (k/ilfreq2)*ilfreq2) {
                    if (ihead == 0) {
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            " i=0   j   k     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max");
                    }
                    ihead = ihead + 1;
                    { eomax = eomui(2);
                      for (int _in=2; _in<=inmax; _in++) if(eomui(_in)>eomax) eomax=eomui(_in); }
                    utmax = utmax / (double)xmach; utmin = utmin / (double)xmach;
                    nou(3) = std::min(nou(3)+1, ibufdim);
                    std::snprintf(bou(nou(3),3), 120,
                        "%9d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%6d%6d%4d%4d%4d",
                        j, k, (float)yplus, (float)utmax,
                        (float)utmin, (float)fblmax, (float)ymax, (float)eomax,
                        inswtch, infmax, in1, in2, inmax);
                }
            }
            if (iprint > 1) {
                if (j == (j/ipfreq1)*ipfreq1 && k == (k/ipfreq2)*ipfreq2) {
                    ihead = 0;
                    nou(3) = std::min(nou(3)+1, ibufdim);
                    std::snprintf(bou(nou(3),3), 120,
                        "   j   k   i      disn      utot      eomu       fbl      amun    eomu-o");
                    for (in = 1; in <= inmax; in++) {
                        utw = utot(in) / (double)xmach;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "%5d%5d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e",
                            j, k, in, (float)disn(in), (float)utw,
                            (float)eomui(in), (float)fbl(in),
                            (float)amun(in), (float)eomuo(in));
                    }
                }
            }

            label2000:;
        } // j loop
        } // k loop

    } // end if (ibcimin==1 || ibcimax==0)


    // evaluate turbulent viscosity along normal to i=idim wall
    if (ibcimax == 1) {

        if (isklton > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                " computing turb viscosity using Baldwin-Lomax, I=idim, block=%4d", nblt);
        }
        if (isklton > 0 && ideg[1-1] == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned on");
        }
        if (isklton > 0 && ideg[1-1] != 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   Degani-Schiff option turned off");
        }

        if (inmx < idim) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "error: inmx, idim=%5d%5d", inmx, idim);
            { int m1 = -1; termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou); }
        }

        ihead = 0;
        // do 2007 k=1,kdim1; do 2007 j=1,jdim1
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {

            iloop  = (int)(0.80 * idim);
            iloop1 = iloop - 1;
            inmaxt  = iloop;
            inmaxt1 = inmaxt - 1;

            iloop  = idim1 - iloop + 1;
            iloop  = std::max(1, iloop);
            iloop1 = iloop + 1;
            inmax  = iloop;
            inmax1 = inmax + 1;

            in1 = (int)(inmaxt1 * 0.20) + 1;
            if (iwf[1-1] > 0) { in1 = 3; }
            in2 = (int)(inmaxt1 * 0.80) + 1;
            in1 = idim1 - in1 + 1;
            in2 = idim1 - in2 + 1;

            // Chimera scheme modification
            if (iover == 1) {
                if ((float)blank(j,k,in1) == 0.f) {
                    for (i = 1; i <= idim1; i++) { eoms(j,k,i) = 0.0e0; }
                    if (iprint >= 1) {
                        if (j == j/ilfreq1*ilfreq1 && k == k/ilfreq2*ilfreq2) {
                            if (ihead == 0) {
                                nou(3) = std::min(nou(3)+1, ibufdim);
                                std::snprintf(bou(nou(3),3), 120,
                                    " k=0   j   i     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max 1/2");
                            }
                            nou(3) = std::min(nou(3)+1, ibufdim);
                            std::snprintf(bou(nou(3),3), 120,
                                "%9d%5d     skipping eddy viscosity evaluation at this station due to hole at-near wall", j, k);
                            ihead = ihead + 1;
                        }
                    }
                    goto label2007;
                }
                in21 = in2;
                for (in = in1; in >= in2; in--) {
                    if ((float)blank(j,k,in) == 0.f) goto label6724;
                    in21 = in;
                }
                label6724:;
                in2 = in21;
            }

            // viscosity thru sutherland relation at cell centers
            for (i = idim1; i >= iloop1; i--) {
                t5 = (double)gamma * q(j,k,i,5) / q(j,k,i,1);
                t6 = std::sqrt(t5);
                amun(i) = c2bp * t5 * t6 / (c2b + t5);
            }
            for (i = idim1; i >= iloop1; i--) {
                rhon(i)  = q(j,k,i,1);
                vortn(i) = vor(j,k,i);
                { double tmp = snim(j,k,i); disn(i) = std::abs(tmp); }
                utot(i)  = std::sqrt( q(j,k,i,2)*q(j,k,i,2)
                                    + q(j,k,i,3)*q(j,k,i,3)
                                    + q(j,k,i,4)*q(j,k,i,4) );
            }
            rhon(idim) = rhon(idim1); amun(idim) = amun(idim1);
            disn(idim) = 0.0; utot(idim) = utot(idim1); vortn(idim) = vortn(idim1);

            if ((float)bci(j,k,2) == 1.f) {
                t5 = (double)gamma * qi0(j,k,5,3) / qi0(j,k,1,3);
                t6 = std::sqrt(t5);
                amun(idim)  = c2bp * t5 * t6 / (c2b + t5);
                rhon(idim)  = qi0(j,k,1,1);
                utot(idim)  = std::sqrt( qi0(j,k,2,3)*qi0(j,k,2,3)
                                       + qi0(j,k,3,3)*qi0(j,k,3,3)
                                       + qi0(j,k,4,3)*qi0(j,k,4,3) );
                vortn(idim) = (utot(idim1) - utot(idim)) / disn(idim1);
            }

            { double tmp = rhon(idim)*vortn(idim)*(double)reue*xmi/amun(idim);
              yplusoy = std::abs(tmp); }
            yplusoy = std::sqrt(yplusoy) * aplusi;
            yplus   = yplusoy * disn(idim1) * aplus * 2.0;


            if ((float)bci(j,k,2) == 0.f) {
                term = 1.0 - std::exp(-50.0e0);
                for (in = idim1; in >= inmax1; in--) { damp(in) = term; }
            } else {
                for (in = idim1; in >= inmax1; in--) {
                    if (iwall == 0) {
                        term = yplusoy * disn(in) * std::sqrt(rhon(in)/rhon(idim)) * (amun(idim)/amun(in));
                    } else { term = yplusoy * disn(in); }
                    damp(in) = 1.0e0 - std::exp(-term);
                }
            }

            fbl(idim) = 0.0e0; eomui(idim) = 0.0e0;
            for (in = idim1; in >= inmax1; in--) {
                fbl(in)   = vortn(in) * disn(in) * damp(in);
                amixl     = vk * disn(in) * damp(in);
                eomui(in) = (double)reue * xmi * rhon(in) * vortn(in) * amixl * amixl;
            }

            { utmax = utot(inmax1); utmin = utot(inmax1);
              for (int _in=inmax1; _in<=inmax1+inmaxt1-1; _in++) {
                  if(utot(_in)>utmax) utmax=utot(_in);
                  if(utot(_in)<utmin) utmin=utot(_in); } }
            if ((float)bci(j,k,2) == 1.f) utmin = utot(idim);

            fblmax = 1.0e-10; infmax = in1 + 1;
            for (in = in1+1; in >= in2+1; in--) {
                if ((float)fbl(in) > (float)fblmax) { fblmax = fbl(in); infmax = in; }
                if (ideg[1-1] == 1 && (float)fbl(in) < 0.9e0*(float)fblmax) goto label8175;
            }
            label8175:;

            if (ifit > 0) {
                if (infmax < in1 && infmax > in2) {
                    dfm = fbl(infmax) - fbl(infmax+1);
                    dfp = fbl(infmax) - fbl(infmax-1);
                    if ((float)fbl(infmax+1) < (float)fbl(infmax-1)) {
                        if ((float)dfm > 0.f) {
                            ymax = disn(infmax) + 0.5*(disn(infmax-1)-disn(infmax))*(1.0-dfp/dfm);
                        } else { ymax = disn(infmax); }
                    } else {
                        if ((float)dfp > 0.f) {
                            ymax = disn(infmax) - 0.5*(disn(infmax)-disn(infmax+1))*(1.0-dfm/dfp);
                        } else { ymax = disn(infmax); }
                    }
                }
            } else { ymax = disn(infmax); }

            { double tmp1e20 = 1.0e-20; ymax = ccomplex_ns::ccmaxcr(ymax, tmp1e20); }
            udif = utmax - utmin;
            fwake = ymax * fblmax;
            fwake2 = cwk * udif * udif * ymax / fblmax;
            { fwake = std::min(fwake2, fwake); }

            for (in = idim1; in >= inmax1; in--) {
                eomuo(in) = 1.0e0 / (1.0e0 + 5.5e0 * std::pow(ckleb*disn(in)/ymax, 6));
                eomuo(in) = ckout * eomuo(in) * rhon(in) * fwake;
            }
            eomuo(idim) = eomuo(idim1);

            for (in = idim1; in >= in2; in--) {
                if ((float)eomui(in) > (float)eomuo(in)) goto label8197;
            }
            label8197:
            inswtch = in + 1;
            for (in = inswtch; in >= inmax1; in--) { eomui(in) = eomuo(in); }

            // fill eomu array corresponding to i=idim wall
            if (ibcimin == 1) {
                // for i=1 & i=idim walls
                for (i = idim1; i >= iloop1; i--) {
                    s1 = sni0(j,k,i) * sni0(j,k,i);
                    s2 = snim(j,k,i) * snim(j,k,i);
                    eoms(j,k,i) = ( s1*eomui(k) + s2*eoms(j,k,i) ) / ( s1 + s2 );
                }
            } else {
                // for i=idim wall only
                for (i = idim1; i >= iloop1; i--) { eoms(j,k,i) = eomui(i); }
                if (iloop1 > 1) {
                    for (i = iloop; i >= 1; i--) { eoms(j,k,i) = 0.0; }
                }
            }

            if (iprint >= 1) {
                if (j == (j/ilfreq1)*ilfreq1 && k == (k/ilfreq2)*ilfreq2) {
                    if (ihead == 0) {
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            " i=0   j   k     yplus     utmax     ut(2)      fmax      ymax    emumax swtch  fmax in1 in2 max");
                    }
                    ihead = ihead + 1;
                    { eomax = eomui(idim1);
                      for (int _in=idim1; _in>=idim1-inmax1+1; _in--) if(eomui(_in)>eomax) eomax=eomui(_in); }
                    utmax = utmax / (double)xmach; utmin = utmin / (double)xmach;
                    nou(3) = std::min(nou(3)+1, ibufdim);
                    std::snprintf(bou(nou(3),3), 120,
                        "%9d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%6d%6d%4d%4d%4d",
                        j, k, (float)yplus, (float)utmax,
                        (float)utmin, (float)fblmax, (float)ymax, (float)eomax,
                        inswtch, infmax, in1, in2, inmax);
                }
            }
            if (iprint > 1) {
                if (j == (j/ipfreq1)*ipfreq1 && k == (k/ipfreq2)*ipfreq2) {
                    ihead = 0;
                    nou(3) = std::min(nou(3)+1, ibufdim);
                    std::snprintf(bou(nou(3),3), 120,
                        "   j   k   i      disn      utot      eomu       fbl      amun    eomu-o");
                    for (in = idim1; in >= inmax1; in--) {
                        utw = utot(in) / (double)xmach;
                        nou(3) = std::min(nou(3)+1, ibufdim);
                        std::snprintf(bou(nou(3),3), 120,
                            "%5d%5d%5d%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e",
                            j, k, in, (float)disn(in), (float)utw,
                            (float)eomui(in), (float)fbl(in),
                            (float)amun(in), (float)eomuo(in));
                    }
                }
            }

            label2007:;
        } // j loop
        } // k loop

    } // end if (ibcimax==1)


    // Lines 2049-2132: form composite eddy-viscosity from k & i walls (ivisc(3)>1)
    if (ivisc[3-1] > 1) {
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            if (ibckmin == 0 && ibckmax == 1) {
                s1 = snkm(j,k,i) * snkm(j,k,i);
            } else if (ibckmin == 1 && ibckmax == 1) {
                { double a = snk0(j,k,i), b = snkm(j,k,i); s1 = std::min(a,b); }
                s1 = s1 * s1;
            } else {
                s1 = snk0(j,k,i) * snk0(j,k,i);
            }
            if (ibcimin == 0 && ibcimax == 1) {
                s2 = snim(j,k,i) * snim(j,k,i);
            } else if (ibcimin == 1 && ibcimax == 1) {
                { double a = sni0(j,k,i), b = snim(j,k,i); s2 = std::min(a,b); }
                s2 = s2 * s2;
            } else {
                s2 = sni0(j,k,i) * sni0(j,k,i);
            }
            vist3d(j,k,i) = ( vist3d(j,k,i)*s2 + eoms(j,k,i)*s1 ) / ( s1 + s2 );
        }}}

        // smooth interior values
        nsmurf = 3;
        if (nsmurf == 0) return;

        for (nnn = 1; nnn <= nsmurf; nnn++) {
            for (i = 1; i <= idim1; i++) {
            for (k = 1; k <= kdim1; k++) {
            for (j = 1; j <= jdim1; j++) {
                eoms(j,k,i) = vist3d(j,k,i);
            }}}
            for (i = 2; i <= idim1-1; i++) {
            for (k = 2; k <= kdim1-1; k++) {
            for (j = 2; j <= jdim1-1; j++) {
                vist3d(j,k,i) = 0.125e0 * ( eoms(j+1,k-1,i-1) + eoms(j-1,k-1,i-1)
                                           + eoms(j+1,k+1,i-1) + eoms(j-1,k+1,i-1)
                                           + eoms(j+1,k-1,i+1) + eoms(j-1,k-1,i+1)
                                           + eoms(j+1,k+1,i+1) + eoms(j-1,k+1,i+1) );
            }}}
        }

    } else {
        // single turbulent wall - install vist3d array as eoms array
        for (i = 1; i <= idim1; i++) {
        for (k = 1; k <= kdim1; k++) {
        for (j = 1; j <= jdim1; j++) {
            vist3d(j,k,i) = eoms(j,k,i);
        }}}
    }

    // laminar region zeroing for ivisc(1)>1
    if (isklton > 0) {
        if (ilamlo == 0 || jlamlo == 0 || klamlo == 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "   block %4d in B-L turb model has no laminar regions", nblt);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "   block %4d in B-L turb model - laminar region is:", nblt);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120,
                "   i=%5d to%5d, j=%5d to%5d, k=%5d to%5d",
                ilamlo, ilamhi, jlamlo, jlamhi, klamlo, klamhi);
        }
    }
    for (i = ilamlo; i <= ilamhi-1; i++) {
    for (k = klamlo; k <= klamhi-1; k++) {
    for (j = jlamlo; j <= jlamhi-1; j++) {
        vist3d(j,k,i) = 0.0;
    }}}
    return;

} // end blomax

} // namespace blomax_ns
