// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "plot3t.h"
#include "runtime/fortran_io.h"
#include "delv.h"
#include "ccomplex.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace plot3t_ns {

void plot3t(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> q, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> xw, FortranArray3DRef<double> blank2, FortranArray3DRef<double> blank, FortranArray4DRef<double> xg, FortranArray3DRef<double> vist3d, int& iover, int& nblk, FortranArray1DRef<int> nmap, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> smin, FortranArray4DRef<double> ux, FortranArray4DRef<double> turre, FortranArray3DRef<double> vol, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> wt, FortranArray3DRef<double> cmuv, int& jdw, int& kdw, int& idw, int& nplots, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay, int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d, int& nblock, FortranArray1DRef<int> nblkpt, FortranArray3DRef<double> volj0, FortranArray3DRef<double> volk0, FortranArray3DRef<double> voli0, FortranArray1DRef<double> vormax, FortranArray1DRef<int> ivmax, FortranArray1DRef<int> jvmax, FortranArray1DRef<int> kvmax, int& nummem, int& ifuncdim)
{
    // COMMON block aliases
    int&     ibin        = cmn_bin.ibin;
    int&     iblnk       = cmn_bin.iblnk;
    float&   gamma_f     = cmn_fluid.gamma;
    float&   cbar        = cmn_fluid2.cbar;
    int&     i2d         = cmn_twod.i2d;
    float*   tur10       = cmn_ivals.tur10;   // 0-based
    int&     ilamlo      = cmn_lam.ilamlo;
    int&     ilamhi      = cmn_lam.ilamhi;
    int&     jlamlo      = cmn_lam.jlamlo;
    int&     jlamhi      = cmn_lam.jlamhi;
    int&     klamlo      = cmn_lam.klamlo;
    int&     klamhi      = cmn_lam.klamhi;
    float&   xmach       = cmn_info.xmach;
    float&   alpha       = cmn_info.alpha;
    int&     ivmx        = cmn_maxiv.ivmx;
    float&   reue        = cmn_reyue.reue;
    float&   tinf        = cmn_reyue.tinf;
    float&   time_f      = cmn_unst.time;
    int&     ialph       = cmn_igrdtyp.ialph;
    int&     lhdr        = cmn_moov.lhdr;
    float&   radtodeg    = cmn_conversion.radtodeg;
    int&     iturbprod   = cmn_turbconv.iturbprod;
    int&     ikoprod     = cmn_konew.ikoprod;
    int&     ifunct      = cmn_plot3dtyp.ifunct;
    int&     i_nonlin    = cmn_constit.i_nonlin;
    float&   c_nonlin_f  = cmn_constit.c_nonlin;
    float&   snonlin_lim_f = cmn_constit.snonlin_lim;
    int&     i_qcr2000   = cmn_constit.i_qcr2000;
    int&     i_qcr2013   = cmn_constit.i_qcr2013;
    int&     i_qcr2013v  = cmn_constit.i_qcr2013v;
    float&   c10_f       = cmn_easmv.c10;
    float&   c2_e_f      = cmn_easmv.c2;
    float&   c3_e_f      = cmn_easmv.c3;
    float&   c4_e_f      = cmn_easmv.c4;
    float&   c5_e_f      = cmn_easmv.c5;
    int&     ieasm_type  = cmn_easmv.ieasm_type;
    int&     ifort50write = cmn_writestuff.ifort50write;
    int&     j_ifort50write = cmn_writestuff.j_ifort50write;
    int&     i_ifort50write = cmn_writestuff.i_ifort50write;

    // local variables
    int    jdim1, kdim1, idim1;
    int    jw, kw, iw;
    int    i, j, k, l, m;
    int    inter, lmax1, nbl, lst, lcoord, lend;
    int    j21, j22, k21, k22, ll, mblk;
    int    n, it, ir, itime, iti;
    int    is, ie, js, je, ks, ke;
    int    is1, ie1, js1, je1, ks1, ke1;
    int    nblc, nblcc, ibembed;
    double re, ccr2;
    double alphaw, xmachw, alphww, reuew, timew;
    double cutoff;
    double s11, s22, s33, s11star, s22star, s33star;
    double s12, s13, s23, w12, w13, w23;
    double xis, wis, c1, gg;
    double zk;
    double t11, t22, t33, t12, t13, t23;
    double alpa1, alpa2;
    double omegatemp, eta, squig, factre;
    double al10, al1, al3, al4, eta1;
    double eta1_girimaji, c4new;
    double denom;
    double s11t, s22t, s33t;
    double o11, o22, o33, o12, o13, o23;
    double t11p, t22p, t33p, t12p, t13p, t23p;
    double xis2013;
    double pk;
    double uprime, vprime, wprime;
    double uuprime, vvprime, wwprime;
    double uvprime, uwprime, vwprime;
    double ske, prod, fd, povere, turbre, def, tint;
    double b11, b22, b33, b12, b13, b23;
    double xii, xiii, ff, qinvar;
    double uuu, vvv, www;
    double sigma;
    double r11slow, r22slow, r33slow, r12slow, r13slow, r23slow;
    double r11rapid, r22rapid, r33rapid, r12rapid, r13rapid, r23rapid;
    double r11, r22, r33, r12, r13, r23;
    double epsilon_loc, epsill;
    double c2b, c2bp, tt, fnu;
    double dist, velterm, rd;
    double temp;
    double yget, qset, fnuw, utau, ypl, upl, zkplus;
    double uuplus, vvplus, wwplus, uwplus, vwplus, uvplus;
    double uloc, ttx, eplus, cmuxx;
    double s1kbk1, s2kbk2, s3kbk3, s1kbk2, s1kbk3, s2kbk3;
    double b1ksk1, b2ksk2, b3ksk3, b1ksk2, b1ksk3, b2ksk3;
    double sklblk;
    double w11, w22, w33;
    double w1kbk1, w2kbk2, w3kbk3, w1kbk2, w1kbk3, w2kbk3;
    double b1kwk1, b2kwk2, b3kwk3, b1kwk2, b1kwk3, b2kwk3;
    double arg, f, c2star;
    // promote float COMMON fields to double for arithmetic
    double c10  = (double)c10_f;
    double c2_e = (double)c2_e_f;
    double c3_e = (double)c3_e_f;
    double c4_e = (double)c4_e_f;
    double c5_e = (double)c5_e_f;
    double c_nonlin   = (double)c_nonlin_f;
    double snonlin_lim = (double)snonlin_lim_f;



    jdim1 = jdim - 1;
    kdim1 = kdim - 1;
    idim1 = idim - 1;

    re   = (double)reue / (double)xmach;
    ccr2 = 2.5;

    // initialize xw and xg arrays
    jw = (j2 - j1) / j3 + 1;
    kw = (k2 - k1) / k3 + 1;
    iw = (i2 - i1) / i3 + 1;
    for (j = 1; j <= jw; j++) {
        for (k = 1; k <= kw; k++) {
            for (i = 1; i <= iw; i++) {
                for (l = 1; l <= ifuncdim; l++) {
                    xw(j,k,i,l) = 0.0;
                }
                for (l = 1; l <= 4; l++) {
                    xg(j,k,i,l) = 0.0;
                }
            }
        }
    }

    // assign single precision scalars
    alphaw = (double)radtodeg * ((double)alpha + thetay(nblk));
    xmachw = (double)xmach;
    alphww = alphaw;
    reuew  = (double)reue;
    timew  = (double)time_f;

    // get ux values
    delv_ns::delv(jdim, kdim, idim, q, sj, sk, si, vol, ux, wt,
                  blank, iover, qj0, qk0, qi0, bcj, bck, bci, nblk,
                  volj0, volk0, voli0, maxbl, vormax, ivmax, jvmax, kvmax);

    // average grid points to get cell centers and load into 1st 3 locations of xg
    iw = 0;
    for (i = i1; i <= i2; i += i3) {
        iw = iw + 1;
        kw = 0;
        for (k = k1; k <= k2; k += k3) {
            kw = kw + 1;
            jw = 0;
            for (j = j1; j <= j2; j += j3) {
                jw = jw + 1;
                xg(jw,kw,iw,1) = 0.125*(  x(j,k,i)     + x(j,k,i+1)
                                         + x(j+1,k,i)   + x(j+1,k,i+1)
                                         + x(j,k+1,i)   + x(j,k+1,i+1)
                                         + x(j+1,k+1,i) + x(j+1,k+1,i+1));
                xg(jw,kw,iw,2) = 0.125*(  y(j,k,i)     + y(j,k,i+1)
                                         + y(j+1,k,i)   + y(j+1,k,i+1)
                                         + y(j,k+1,i)   + y(j,k+1,i+1)
                                         + y(j+1,k+1,i) + y(j+1,k+1,i+1));
                xg(jw,kw,iw,3) = 0.125*(  z(j,k,i)     + z(j,k,i+1)
                                         + z(j+1,k,i)   + z(j+1,k,i+1)
                                         + z(j,k+1,i)   + z(j,k+1,i+1)
                                         + z(j+1,k+1,i) + z(j+1,k+1,i+1));
            }
        }
    }

    // set iblank (blank2) array
    // don't compute blanking info if iblnk = 0
    if (iblnk > 0) {
        // assign default iblank (blank2) array
        for (i = 1; i <= idim; i++)
            for (j = 1; j <= jdim; j++)
                for (k = 1; k <= kdim; k++)
                    blank2(j,k,i) = 1.0;

        // iblank (blank2) array for generalized patch interface boundaries
        if (std::abs(ninter) > 0) {
            for (inter = 1; inter <= std::abs(ninter); inter++) {
                lmax1  = iindex(inter,1);
                nbl    = iindex(inter,lmax1+2);
                if (nbl != nblk) continue;
                lst    = iindex(inter,2*lmax1+5);
                lcoord = iindex(inter,2*lmax1+3)/10;
                lend   = iindex(inter,2*lmax1+3) - lcoord*10;
                j21    = iindex(inter,2*lmax1+6);
                j22    = iindex(inter,2*lmax1+7);
                k21    = iindex(inter,2*lmax1+8);
                k22    = iindex(inter,2*lmax1+9);
                if (lcoord == 1) {
                    if (lend == 1) i = 1;
                    if (lend == 2) i = idim1;
                    for (j = j21; j <= j22-1; j++) {
                        for (k = k21; k <= k22-1; k++) {
                            ll   = lst + (j22-j21)*(k-k21) + (j-j21);
                            mblk = iindex(inter, nblkpt(ll)+1);
                            blank2(j,k,i) = -(double)nmap(mblk);
                        }
                    }
                }
                if (lcoord == 2) {
                    if (lend == 1) j = 1;
                    if (lend == 2) j = jdim1;
                    for (i = k21; i <= k22-1; i++) {
                        for (k = j21; k <= j22-1; k++) {
                            ll   = lst + (j22-j21)*(i-k21) + (k-j21);
                            mblk = iindex(inter, nblkpt(ll)+1);
                            blank2(j,k,i) = -(double)nmap(mblk);
                        }
                    }
                }
                if (lcoord == 3) {
                    if (lend == 1) k = 1;
                    if (lend == 2) k = kdim1;
                    for (i = k21; i <= k22-1; i++) {
                        for (j = j21; j <= j22-1; j++) {
                            ll   = lst + (j22-j21)*(i-k21) + (j-j21);
                            mblk = iindex(inter, nblkpt(ll)+1);
                            blank2(j,k,i) = -(double)nmap(mblk);
                        }
                    }
                }
            }
        }

        // iblank (blank2) array for 1:1 interface boundaries
        if (nbli > 0) {
            for (n = 1; n <= std::abs(nbli); n++) {
                if (nblon(n) >= 0) {
                    if (nblk == nblkk(1,n) || nblk == nblkk(2,n)) {
                        it = 1;
                        ir = 2;
                        if (nblk == nblkk(2,n)) {
                            it = 2;
                            ir = 1;
                        }
                        itime = 1;
                        if (nblkk(1,n) == nblkk(2,n)) itime = 2;
                        for (iti = 1; iti <= itime; iti++) {
                            if (iti > 1) {
                                it = 1;
                                ir = 2;
                            }
                            is = limblk(it,1,n);
                            ie = limblk(it,4,n);
                            js = limblk(it,2,n);
                            je = limblk(it,5,n);
                            ks = limblk(it,3,n);
                            ke = limblk(it,6,n);
                            is1 = std::min(is,ie);
                            ie1 = std::max(is,ie);
                            js1 = std::min(js,je);
                            je1 = std::max(js,je);
                            ks1 = std::min(ks,ke);
                            ke1 = std::max(ks,ke);
                            ie1 = std::min(ie1,idim1);
                            je1 = std::min(je1,jdim1);
                            ke1 = std::min(ke1,kdim1);
                            is1 = std::min(is1,idim1);
                            js1 = std::min(js1,jdim1);
                            ks1 = std::min(ks1,kdim1);
                            for (i = is1; i <= ie1; i++)
                                for (j = js1; j <= je1; j++)
                                    for (k = ks1; k <= ke1; k++)
                                        blank2(j,k,i) = -(double)nmap(nblkk(ir,n));
                        }
                    }
                }
            }
        }

        // iblank (blank2) array for embedded grids
        ibembed = 1;
        if (ibembed > 0) {
            for (nblc = 1; nblc <= nblock; nblc++) {
                if (nblk == nblc) continue;
                nblcc = nblcg(nblc);
                if (nblcc == nblk) {
                    js = jsg(nblc);
                    ks = ksg(nblc);
                    is = isg(nblc);
                    je = jeg(nblc);
                    ke = keg(nblc);
                    ie = ieg(nblc);
                    if (je > js) je = je - 1;
                    if (ke > ks) ke = ke - 1;
                    if (ie > is) ie = ie - 1;
                    for (i = is; i <= ie; i++)
                        for (j = js; j <= je; j++)
                            for (k = ks; k <= ke; k++)
                                blank2(j,k,i) = 0.0;
                }
            }
        }

        // iblank (blank2) array for overlapped grids
        if (iover == 1) {
            for (i = 1; i <= idim1; i++)
                for (j = 1; j <= jdim1; j++)
                    for (k = 1; k <= kdim1; k++)
                        if ((float)blank(j,k,i) == 0.0f) blank2(j,k,i) = 0.0;
        }
    } // end if (iblnk > 0)



    // plot3d data - write header to unit 11
    if (lhdr > 0 && myid == myhost) {
        if (i2d == 1) {
            if (iblnk == 0) {
                fortran_write_unit(11, "writing turb plot3d file for JDIM X KDIM =%5d x %5d grid\n", jdim, kdim);
                fortran_write_unit(11, "   plot3d files to be read with 2d/mgrid qualifiers\n");
            } else {
                fortran_write_unit(11, "writing turb plot3d file for JDIM X KDIM =%5d x %5d grid\n", jdim, kdim);
                fortran_write_unit(11, "   plot3d files to be read with 2d/mgrid/blank qualifiers\n");
            }
        } else {
            fortran_write_unit(11, " writing turb plot3d file for IDIM X JDIM X KDIM =%5d x %5d x %5d grid\n", idim, jdim, kdim);
            if (iblnk == 0) {
                fortran_write_unit(11, "   plot3d files to be read with /mgrid qualifiers\n");
            } else {
                fortran_write_unit(11, "   plot3d files to be read with /mgrid/blank qualifiers\n");
            }
        }
        if (ifunct == 0) {
            fortran_write_unit(11, "   plot3d file is Q file of turbulence quantities\n");
            if (i2d == 1) {
                fortran_write_unit(11, "  Defaults are production-term, uw, uu, ww for 2D\n");
            } else {
                fortran_write_unit(11, "   Defaults are production-term, uw, uu, ww, Sk/eps for 3D\n");
            }
        } else {
            fortran_write_unit(11, "   plot3d file is function file of%5d turbulence quantities\n", ifunct);
        }
    }

    // load blank2 into 4th location of the single precision xg array
    iw = 0;
    for (i = i1; i <= i2; i += i3) {
        iw = iw + 1;
        jw = 0;
        for (j = j1; j <= j2; j += j3) {
            jw = jw + 1;
            kw = 0;
            for (k = k1; k <= k2; k += k3) {
                kw = kw + 1;
                xg(jw,kw,iw,4) = blank2(j,k,i);
            }
        }
    }

    if (myid == myhost && iblnk > 0) {
        // ialph > 0 for a grid that was read in plot3d format with alpha measured
        // in the xy plane (TLNS3D convention)
        if (ialph != 0 && i2d != 1) {
            for (i = 1; i <= iw; i++) {
                for (j = 1; j <= jw; j++) {
                    for (k = 1; k <= kw; k++) {
                        temp        = xg(j,k,i,2);
                        xg(j,k,i,2) = xg(j,k,i,3);
                        xg(j,k,i,3) = -temp;
                    }
                }
            }
        }

        // output grid to unit 3
        if (ibin == 0) {
            if (i2d == 0) {
                if (iblnk == 0) {
                    // write(3,'(5e14.6)') (((xg(j,k,i,1),...),xg(2),...),xg(3),...)
                    {
                        FILE* f3 = fortran_get_unit(3);
                        int cnt = 0;
                        for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                            fprintf(f3, "%14.6E", (float)xg(j,k,i,1));
                            if (++cnt % 5 == 0) fprintf(f3, "\n");
                        }
                        for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                            fprintf(f3, "%14.6E", (float)xg(j,k,i,2));
                            if (++cnt % 5 == 0) fprintf(f3, "\n");
                        }
                        for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                            fprintf(f3, "%14.6E", (float)xg(j,k,i,3));
                            if (++cnt % 5 == 0) fprintf(f3, "\n");
                        }
                        if (cnt % 5 != 0) fprintf(f3, "\n");
                    }
                } else {
                    {
                        FILE* f3 = fortran_get_unit(3);
                        int cnt = 0;
                        for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                            fprintf(f3, "%14.6E", (float)xg(j,k,i,1));
                            if (++cnt % 5 == 0) fprintf(f3, "\n");
                        }
                        for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                            fprintf(f3, "%14.6E", (float)xg(j,k,i,2));
                            if (++cnt % 5 == 0) fprintf(f3, "\n");
                        }
                        for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                            fprintf(f3, "%14.6E", (float)xg(j,k,i,3));
                            if (++cnt % 5 == 0) fprintf(f3, "\n");
                        }
                        if (cnt % 5 != 0) fprintf(f3, "\n");
                        // write iblank as integers
                        int icnt = 0;
                        for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                            fprintf(f3, "%8d", (int)xg(j,k,i,4));
                            if (++icnt % 10 == 0) fprintf(f3, "\n");
                        }
                        if (icnt % 10 != 0) fprintf(f3, "\n");
                    }
                }
            } else {
                // i2d == 1: only x and z (indices 1 and 3)
                if (iblnk == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    int cnt = 0;
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        fprintf(f3, "%14.6E", (float)xg(j,k,i,1));
                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        fprintf(f3, "%14.6E", (float)xg(j,k,i,3));
                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                    }
                    if (cnt % 5 != 0) fprintf(f3, "\n");
                } else {
                    FILE* f3 = fortran_get_unit(3);
                    int cnt = 0;
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        fprintf(f3, "%14.6E", (float)xg(j,k,i,1));
                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        fprintf(f3, "%14.6E", (float)xg(j,k,i,3));
                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                    }
                    if (cnt % 5 != 0) fprintf(f3, "\n");
                    int icnt = 0;
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        fprintf(f3, "%8d", (int)xg(j,k,i,4));
                        if (++icnt % 10 == 0) fprintf(f3, "\n");
                    }
                    if (icnt % 10 != 0) fprintf(f3, "\n");
                }
            }
        } else {
            // binary output
            if (i2d == 0) {
                if (iblnk == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,1); fwrite(&v, sizeof(float), 1, f3);
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,2); fwrite(&v, sizeof(float), 1, f3);
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,3); fwrite(&v, sizeof(float), 1, f3);
                    }
                } else {
                    FILE* f3 = fortran_get_unit(3);
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,1); fwrite(&v, sizeof(float), 1, f3);
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,2); fwrite(&v, sizeof(float), 1, f3);
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,3); fwrite(&v, sizeof(float), 1, f3);
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        int iv = (int)xg(j,k,i,4); fwrite(&iv, sizeof(int), 1, f3);
                    }
                }
            } else {
                // i2d == 1
                if (iblnk == 0) {
                    FILE* f3 = fortran_get_unit(3);
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,1); fwrite(&v, sizeof(float), 1, f3);
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,3); fwrite(&v, sizeof(float), 1, f3);
                    }
                } else {
                    FILE* f3 = fortran_get_unit(3);
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,1); fwrite(&v, sizeof(float), 1, f3);
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        float v = (float)xg(j,k,i,3); fwrite(&v, sizeof(float), 1, f3);
                    }
                    for (k = 1; k <= kw; k++) for (j = 1; j <= jw; j++) for (i = 1; i <= iw; i++) {
                        int iv = (int)xg(j,k,i,4); fwrite(&iv, sizeof(int), 1, f3);
                    }
                }
            }
        }
    } // end if (myid == myhost && iblnk > 0)



    // turbulent quantities - initialize scalars
    uprime  = 0.0; vprime  = 0.0; wprime  = 0.0;
    uuprime = 0.0; vvprime = 0.0; wwprime = 0.0;
    uvprime = 0.0; uwprime = 0.0; vwprime = 0.0;
    ske     = 0.0; prod    = 0.0; fd      = 0.0;
    povere  = 0.0; turbre  = 0.0; def     = 0.0;
    tint    = 0.0;
    b11     = 0.0; b22     = 0.0; b33     = 0.0;
    b12     = 0.0; b13     = 0.0; b23     = 0.0;
    xii     = 0.0; xiii    = 0.0; ff      = 0.0;
    zk      = 0.0; qinvar  = 0.0;
    uuu     = 0.0; vvv     = 0.0; www     = 0.0;
    sigma   = 0.0;
    r11slow = 0.0; r22slow = 0.0; r33slow = 0.0;
    r12slow = 0.0; r13slow = 0.0; r23slow = 0.0;
    r11rapid= 0.0; r22rapid= 0.0; r33rapid= 0.0;
    r12rapid= 0.0; r13rapid= 0.0; r23rapid= 0.0;
    r11     = 0.0; r22     = 0.0; r33     = 0.0;
    r12     = 0.0; r13     = 0.0; r23     = 0.0;

    // optional printout to fort.50 header
    if (ifort50write == 1) {
        fortran_write_unit(50, "variables=\"y\",\"uuplus\",\"wwplus\",\"vvplus\",\"uwplus\","
                               "\"uvplus\",\"vwplus\",\"u_scaled\",\"kplus\",\"yplus\",\"uplus\","
                               "\"mu_turbulent/mu_lam_ref\",\"eplus\",\"cmu\"\n");
        fortran_write_unit(50, "#%5d  CFL3D scaled turb results along k-line at j,i=%4d%4d x=%16.6f\n",
                           k2-k1+1, j_ifort50write, i_ifort50write,
                           0.5*(x(j_ifort50write,1,i_ifort50write)+x(j_ifort50write+1,1,i_ifort50write)));
        fortran_write_unit(50, "# NOTE: eplus applies to k-eps models only\n");
    }

    // main loop over i, k, j
    iw = 0;
    for (i = i1; i <= i2; i += i3) {
        iw = iw + 1;
        kw = 0;
        for (k = k1; k <= k2; k += k3) {
            kw = kw + 1;
            jw = 0;
            for (j = j1; j <= j2; j += j3) {
                jw = jw + 1;

                // cutoff
                if (i >= ilamlo && i < ilamhi &&
                    j >= jlamlo && j < jlamhi &&
                    k >= klamlo && k < klamhi) {
                    cutoff = 0.0;
                } else {
                    cutoff = 1.0;
                }

                // Determine Sij and Wij values
                s11 = ux(j,k,i,1);
                s22 = ux(j,k,i,5);
                s33 = ux(j,k,i,9);
                s11star = s11 - ((s11+s22+s33)/3.0);
                s22star = s22 - ((s11+s22+s33)/3.0);
                s33star = s33 - ((s11+s22+s33)/3.0);
                s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4));
                s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7));
                s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8));
                w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4));
                w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7));
                w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8));
                // xis=SijSij, wis=WijWij
                xis = s11*s11 + s22*s22 + s33*s33 +
                      2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                wis = 2.0*w12*w12 + 2.0*w13*w13 + 2.0*w23*w23;
                c1  = 2.0*c10;
                gg  = 1.0/(c1/2.0 + c5_e - 1.0);

                if (ivmx < 6) {
                    // estimate for turb kinetic energy for models that do not compute it directly
                    zk = vist3d(j,k,i)/q(j,k,i,1)*std::sqrt(2.0*xis)/(re*0.31);
                } else {
                    // use exact computed turb kinetic energy
                    zk = turre(j,k,i,2);
                }
                // If using QCR2013, then k is already added, so do not add here
                if (i_qcr2013 == 1 || i_qcr2013v == 1) {
                    zk = 0.0;
                }



                // Stress tensor branches
                if ((ivmx <= 7 || ivmx == 10 || ivmx == 15 || ivmx == 16) && i_nonlin == 0) {
                    t11 = -2.0*vist3d(j,k,i)*s11star;
                    t22 = -2.0*vist3d(j,k,i)*s22star;
                    t33 = -2.0*vist3d(j,k,i)*s33star;
                    // Add 2/3*k to the diagonal stress terms for eddy viscosity models
                    t11 = t11 + 2.0*q(j,k,i,1)*zk*re/3.0;
                    t22 = t22 + 2.0*q(j,k,i,1)*zk*re/3.0;
                    t33 = t33 + 2.0*q(j,k,i,1)*zk*re/3.0;
                    t12 = -2.0*vist3d(j,k,i)*s12;
                    t13 = -2.0*vist3d(j,k,i)*s13;
                    t23 = -2.0*vist3d(j,k,i)*s23;
                    if (i_qcr2000 == 1 || i_qcr2013 == 1 || i_qcr2013v == 1) {
                        denom = std::sqrt(ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,2)*ux(j,k,i,2) +
                                          ux(j,k,i,3)*ux(j,k,i,3) + ux(j,k,i,4)*ux(j,k,i,4) +
                                          ux(j,k,i,5)*ux(j,k,i,5) + ux(j,k,i,6)*ux(j,k,i,6) +
                                          ux(j,k,i,7)*ux(j,k,i,7) + ux(j,k,i,8)*ux(j,k,i,8) +
                                          ux(j,k,i,9)*ux(j,k,i,9)) + 1.0e-20;
                        o11 = 0.0;
                        o22 = 0.0;
                        o33 = 0.0;
                        o12 = 2.0*w12/denom;
                        o13 = 2.0*w13/denom;
                        o23 = 2.0*w23/denom;
                        t11p = t11 - 0.3*(o11*t11+o12*t12+o13*t13 + o11*t11+o12*t12+o13*t13);
                        t22p = t22 - 0.3*(-o12*t12+o22*t22+o23*t23 + -o12*t12+o22*t22+o23*t23);
                        t33p = t33 - 0.3*(-o13*t13-o23*t23+o33*t33 + -o13*t13-o23*t23+o33*t33);
                        t12p = t12 - 0.3*(o11*t12+o12*t22+o13*t23 + -o12*t11+o22*t12+o23*t13);
                        t13p = t13 - 0.3*(o11*t13+o12*t23+o13*t33 + -o13*t11-o23*t12+o33*t13);
                        t23p = t23 - 0.3*(-o12*t13+o22*t23+o23*t33 + -o13*t12-o23*t22+o33*t23);
                        t11 = t11p; t22 = t22p; t33 = t33p;
                        t12 = t12p; t13 = t13p; t23 = t23p;
                    }
                    if (i_qcr2013 == 1) {
                        xis2013 = s11star*s11star + s22star*s22star + s33star*s33star +
                                  2.0*s12*s12 + 2.0*s13*s13 + 2.0*s23*s23;
                        t11 = t11 + ccr2*vist3d(j,k,i)*std::sqrt(2.0*xis2013);
                        t22 = t22 + ccr2*vist3d(j,k,i)*std::sqrt(2.0*xis2013);
                        t33 = t33 + ccr2*vist3d(j,k,i)*std::sqrt(2.0*xis2013);
                    }
                    if (i_qcr2013v == 1) {
                        t11 = t11 + ccr2*vist3d(j,k,i)*std::sqrt(2.0*wis);
                        t22 = t22 + ccr2*vist3d(j,k,i)*std::sqrt(2.0*wis);
                        t33 = t33 + ccr2*vist3d(j,k,i)*std::sqrt(2.0*wis);
                    }
                } else if (ivmx == 12) {
                    // nonlinear k-eps (ivmx==12)
                    alpa1 = (2.0-c4_e)/2.0*gg;
                    alpa2 = (2.0-c3_e)*gg;
                    {
                        double a1 = turre(j,k,i,1);
                        double a2 = (double)tur10[0];
                        omegatemp = ccomplex_ns::ccmax(a1, a2);
                    }
                    eta   = (2.0-c3_e)*(2.0-c3_e)*(gg*gg/4.0)*xis/(omegatemp*re)/(omegatemp*re);
                    squig = (2.0-c4_e)*(2.0-c4_e)*(gg*gg/4.0)*wis/(omegatemp*re)/(omegatemp*re);
                    { double a1=eta,   a2=10.0; eta   = ccomplex_ns::ccmincr(a1,a2); }
                    { double a1=squig, a2=10.0; squig = ccomplex_ns::ccmincr(a1,a2); }
                    factre = (3.0*(1.0+eta) + 0.2e-8*(eta*eta*eta + squig*squig*squig)) /
                             (3.0*(1.0+eta) +    0.2*(eta*eta*eta + squig*squig*squig));
                    t11 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*factre/omegatemp*(-s12*w12 - s13*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre/omegatemp*(s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re;
                    t22 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*factre/omegatemp*(s12*w12 - s23*w23)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre/omegatemp*(s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re;
                    t33 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*factre/omegatemp*(s23*w23 + s13*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre/omegatemp*(s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re;
                    t12 = -2.0*vist3d(j,k,i)*s12
                          - 2.0*alpa1*vist3d(j,k,i)*factre/omegatemp*(s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre/omegatemp*(s11*s12 + s12*s22 + s13*s23)/re;
                    t13 = -2.0*vist3d(j,k,i)*s13
                          - 2.0*alpa1*vist3d(j,k,i)*factre/omegatemp*(s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre/omegatemp*(s11*s13 + s12*s23 + s13*s33)/re;
                    t23 = -2.0*vist3d(j,k,i)*s23
                          - 2.0*alpa1*vist3d(j,k,i)*factre/omegatemp*(s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre/omegatemp*(s12*s13 + s22*s23 + s23*s33)/re;
                } else if (ivmx == 11) {
                    alpa1 = (2.0-c4_e)/2.0*gg;
                    alpa2 = (2.0-c3_e)*gg;
                    eta   = (2.0-c3_e)*(2.0-c3_e)*(gg*gg/4.0)*xis*
                            turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
                    squig = (2.0-c4_e)*(2.0-c4_e)*(gg*gg/4.0)*wis*
                            turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
                    { double a1=eta,   a2=10.0; eta   = ccomplex_ns::ccmincr(a1,a2); }
                    { double a1=squig, a2=10.0; squig = ccomplex_ns::ccmincr(a1,a2); }
                    factre = (3.0*(1.0+eta) + 0.2e-8*(eta*eta*eta + squig*squig*squig)) /
                             (3.0*(1.0+eta) +    0.2*(eta*eta*eta + squig*squig*squig));
                    t11 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re;
                    t22 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s12*w12 - s23*w23)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re;
                    t33 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s23*w23 + s13*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re;
                    t12 = -2.0*vist3d(j,k,i)*s12
                          - 2.0*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s12 + s12*s22 + s13*s23)/re;
                    t13 = -2.0*vist3d(j,k,i)*s13
                          - 2.0*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s13 + s12*s23 + s13*s33)/re;
                    t23 = -2.0*vist3d(j,k,i)*s23
                          - 2.0*alpa1*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
                          + 2.0*alpa2*vist3d(j,k,i)*factre*turre(j,k,i,2)/turre(j,k,i,1)*(s12*s13 + s22*s23 + s23*s33)/re;



                } else if (ivmx == 9 || ivmx == 13) {
                    al10 = 0.25*c10 - 1.0;
                    al1  = 3.8;
                    if (ieasm_type == 0) {
                        al10 = al10 + 1.8864;
                        al1  = al1  - 2.0;
                    }
                    al3  = 0.5*c3_e - 1.0;
                    al4  = 0.5*c4_e - 1.0;
                    eta1 = xis*turre(j,k,i,2)*turre(j,k,i,2)/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
                    alpa1 = -al4/(al10 - eta1*al1*cmuv(j,k,i));
                    alpa2 = -2.0*al3/(al10 - eta1*cmuv(j,k,i)*al1);
                    if (ivmx == 9) {
                        alpa1 = 0.0;
                        alpa2 = 0.0;
                    }
                    t11 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re;
                    t22 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s12*w12 - s23*w23)/re
                          + 2.0*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re;
                    t33 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s23*w23 + s13*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re;
                    t12 = -2.0*vist3d(j,k,i)*s12
                          - 2.0*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s12 + s12*s22 + s13*s23)/re;
                    t13 = -2.0*vist3d(j,k,i)*s13
                          - 2.0*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s11*s13 + s12*s23 + s13*s33)/re;
                    t23 = -2.0*vist3d(j,k,i)*s23
                          - 2.0*alpa1*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
                          + 2.0*alpa2*vist3d(j,k,i)*turre(j,k,i,2)/turre(j,k,i,1)*(s12*s13 + s22*s23 + s23*s33)/re;
                } else if (ivmx == 8 || ivmx == 14) {
                    al10 = 0.25*c10 - 1.0;
                    al1  = 3.8;
                    if (ieasm_type == 0 || ieasm_type == 3 || ieasm_type == 4) {
                        al10 = al10 + 1.8864;
                        al1  = al1  - 2.0;
                    }
                    al3  = 0.5*c3_e - 1.0;
                    al4  = 0.5*c4_e - 1.0;
                    // Girimaji JFM 2000 fix to c4
                    if (ieasm_type == 4) {
                        eta1_girimaji = xis/(xis + wis);
                        if ((float)eta1_girimaji < 0.5f) {
                            c4new = 2.0 - ((2.0-c4_e)*std::pow(eta1_girimaji/(1.0-eta1_girimaji), 0.75));
                        } else {
                            c4new = c4_e;
                        }
                        al4 = 0.5*c4new - 1.0;
                    }
                    eta1  = xis/(turre(j,k,i,1)*re)/(turre(j,k,i,1)*re);
                    alpa1 = -al4/(al10 - eta1*al1*cmuv(j,k,i));
                    alpa2 = -2.0*al3/(al10 - eta1*cmuv(j,k,i)*al1);
                    if (ivmx == 8) {
                        alpa1 = 0.0;
                        alpa2 = 0.0;
                    }
                    t11 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s11 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(-s12*w12 - s13*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re;
                    t22 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s22 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s12*w12 - s23*w23)/re
                          + 2.0*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re;
                    t33 = 2.0*q(j,k,i,1)*turre(j,k,i,2)*re/3.0
                          - 2.0*vist3d(j,k,i)*(s33 - 0.33333*(s11+s22+s33))
                          - 4.0*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s23*w23 + s13*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re;
                    t12 = -2.0*vist3d(j,k,i)*s12
                          - 2.0*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s11*s12 + s12*s22 + s13*s23)/re;
                    t13 = -2.0*vist3d(j,k,i)*s13
                          - 2.0*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
                          + 2.0*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s11*s13 + s12*s23 + s13*s33)/re;
                    t23 = -2.0*vist3d(j,k,i)*s23
                          - 2.0*alpa1*vist3d(j,k,i)/turre(j,k,i,1)*(s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
                          + 2.0*alpa2*vist3d(j,k,i)/turre(j,k,i,1)*(s12*s13 + s22*s23 + s23*s33)/re;
                } else if (i_nonlin != 0) {
                    denom = ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,2)*ux(j,k,i,2) + ux(j,k,i,3)*ux(j,k,i,3) +
                            ux(j,k,i,4)*ux(j,k,i,4) + ux(j,k,i,5)*ux(j,k,i,5) + ux(j,k,i,6)*ux(j,k,i,6) +
                            ux(j,k,i,7)*ux(j,k,i,7) + ux(j,k,i,8)*ux(j,k,i,8) + ux(j,k,i,9)*ux(j,k,i,9);
                    denom = std::sqrt(denom);
                    { double a1=denom, a2=snonlin_lim; denom = ccomplex_ns::ccmax(a1,a2); }
                    s11t = s11 - (s11+s22+s33)/3.0;
                    s22t = s22 - (s11+s22+s33)/3.0;
                    s33t = s33 - (s11+s22+s33)/3.0;
                    t11 = 2.0*q(j,k,i,1)*zk*re/3.0
                          - 2.0*vist3d(j,k,i)*s11t
                          - 8.0*c_nonlin*vist3d(j,k,i)/denom*(-s12*w12 - s13*w13);
                    t22 = 2.0*q(j,k,i,1)*zk*re/3.0
                          - 2.0*vist3d(j,k,i)*s22t
                          - 8.0*c_nonlin*vist3d(j,k,i)/denom*(s12*w12 - s23*w23);
                    t33 = 2.0*q(j,k,i,1)*zk*re/3.0
                          - 2.0*vist3d(j,k,i)*s33t
                          - 8.0*c_nonlin*vist3d(j,k,i)/denom*(s23*w23 + s13*w13);
                    t12 = -2.0*vist3d(j,k,i)*s12
                          - 4.0*c_nonlin*vist3d(j,k,i)/denom*(s11t*w12 - s22t*w12 - s13*w23 - s23*w13);
                    t13 = -2.0*vist3d(j,k,i)*s13
                          - 4.0*c_nonlin*vist3d(j,k,i)/denom*(s11t*w13 + s12*w23 - s23*w12 - s33t*w13);
                    t23 = -2.0*vist3d(j,k,i)*s23
                          - 4.0*c_nonlin*vist3d(j,k,i)/denom*(s12*w13 + s22t*w23 + s13*w12 - s33t*w23);
                }



                // Calculate 2-eqn-specific production term
                if (ivmx == 6 || ivmx == 7 || ivmx == 10 || ivmx == 16) {
                    if (ikoprod == 1) {
                        pk = 2.0*vist3d(j,k,i)*xis/(q(j,k,i,1)*re);
                    } else {
                        pk = 2.0*vist3d(j,k,i)*wis/(q(j,k,i,1)*re);
                    }
                } else if (ivmx == 8 || ivmx == 9 || ivmx == 13 || ivmx == 14) {
                    if (iturbprod == 1) {
                        pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
                               + t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
                               + t13*(ux(j,k,i,3)+ux(j,k,i,7))
                               + t23*(ux(j,k,i,6)+ux(j,k,i,8))) / (q(j,k,i,1)*re);
                    } else {
                        pk = 2.0*vist3d(j,k,i)*xis/(q(j,k,i,1)*re);
                    }
                } else if (ivmx == 11 || ivmx == 12) {
                    pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
                           + t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
                           + t13*(ux(j,k,i,3)+ux(j,k,i,7))
                           + t23*(ux(j,k,i,6)+ux(j,k,i,8))) / (q(j,k,i,1)*re);
                } else if (ivmx == 15) {
                    if (ikoprod == 2) {
                        pk = -(t11*ux(j,k,i,1) + t22*ux(j,k,i,5)
                               + t33*ux(j,k,i,9) + t12*(ux(j,k,i,2)+ux(j,k,i,4))
                               + t13*(ux(j,k,i,3)+ux(j,k,i,7))
                               + t23*(ux(j,k,i,6)+ux(j,k,i,8))) / (q(j,k,i,1)*re);
                    } else if (ikoprod == 1) {
                        pk = 2.0*vist3d(j,k,i)*xis/(q(j,k,i,1)*re);
                    } else {
                        pk = 2.0*vist3d(j,k,i)*wis/(q(j,k,i,1)*re);
                    }
                } else {
                    pk = 2.0*vist3d(j,k,i)*wis/(q(j,k,i,1)*re);
                }
                pk = pk * cutoff;

                // optional printout to fort.50
                if (ifort50write == 1) {
                    if (j == j_ifort50write && i == i_ifort50write) {
                        { double a1 = smin(j,k,i); yget = ccomplex_ns::ccabs(a1); }
                        qset = std::sqrt(q(j,1,i,2)*q(j,1,i,2) + q(j,1,i,3)*q(j,1,i,3) + q(j,1,i,4)*q(j,1,i,4));
                        tt   = (double)gamma_f*q(j,1,i,5)/q(j,1,i,1);
                        c2b  = 198.6/(double)tinf;
                        c2bp = c2b + 1.0;
                        fnuw = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                        { double a1 = smin(j,1,i); double absval = ccomplex_ns::ccabs(a1);
                          utau = std::sqrt(fnuw*qset/(absval*q(j,1,i,1)*re)); }
                        { double a1 = smin(j,k,i); double absval = ccomplex_ns::ccabs(a1);
                          ypl = re*q(j,1,i,1)*utau*absval/fnuw; }
                        upl  = std::sqrt(q(j,k,i,2)*q(j,k,i,2) + q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*q(j,k,i,4)) / utau;
                        zkplus = zk/(utau*utau);
                        uuplus = std::abs(t11)/(re*utau*utau);
                        vvplus = std::abs(t22)/(re*utau*utau);
                        wwplus = std::abs(t33)/(re*utau*utau);
                        uwplus = t13/(re*utau*utau);
                        vwplus = t23/(re*utau*utau);
                        uvplus = t12/(re*utau*utau);
                        uloc   = std::sqrt(q(j,k,i,2)*q(j,k,i,2) + q(j,k,i,3)*q(j,k,i,3) + q(j,k,i,4)*q(j,k,i,4))/(double)xmach;
                        ttx    = (double)gamma_f*q(j,k,i,5)/q(j,k,i,1);
                        fnu    = c2bp*ttx*std::sqrt(ttx)/(c2b+ttx);
                        if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                            eplus = turre(j,k,i,1)*fnu/(q(j,k,i,1)*utau*utau*utau*utau);
                        } else {
                            eplus = 0.0;
                        }
                        if (ivmx == 8 || ivmx == 9 || ivmx >= 13) {
                            cmuxx = -cmuv(j,k,i);
                        } else {
                            cmuxx = 0.09;
                        }
                        fortran_write_unit(50, "%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E%15.5E\n",
                                           yget, uuplus, wwplus, vvplus, uwplus,
                                           uvplus, vwplus, uloc, zkplus, ypl, upl, vist3d(j,k,i), eplus, cmuxx);
                    }
                }



                // u'/u_inf
                uprime = std::sqrt(std::abs(t11)/(q(j,k,i,1)*re))/(double)xmach;
                // v'/u_inf
                if (ialph == 0) {
                    vprime = std::sqrt(std::abs(t22)/(q(j,k,i,1)*re))/(double)xmach;
                } else {
                    vprime = std::sqrt(std::abs(t33)/(q(j,k,i,1)*re))/(double)xmach;
                }
                // w'/u_inf
                if (ialph == 0) {
                    wprime = std::sqrt(std::abs(t33)/(q(j,k,i,1)*re))/(double)xmach;
                } else {
                    wprime = std::sqrt(std::abs(t22)/(q(j,k,i,1)*re))/(double)xmach;
                }
                // u'u'/(u_inf)**2
                uuprime = t11/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                // v'v'/(u_inf)**2
                if (ialph == 0) {
                    vvprime = t22/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                } else {
                    vvprime = t33/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                }
                // w'w'/(u_inf)**2
                if (ialph == 0) {
                    wwprime = t33/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                } else {
                    wwprime = t22/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                }
                // u'v'/(u_inf)**2
                if (ialph == 0) {
                    uvprime = t12/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                } else {
                    uvprime = t13/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                }
                // u'w'/(u_inf)**2
                if (ialph == 0) {
                    uwprime = t13/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                } else {
                    uwprime = -t12/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                }
                // v'w'/(u_inf)**2
                if (ialph == 0) {
                    vwprime = t23/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                } else {
                    vwprime = -t23/(q(j,k,i,1)*re)/(double)xmach/(double)xmach;
                }
                // Laminar Viscosity
                c2b  = (double)cbar/(double)tinf;
                c2bp = c2b + 1.0;
                tt   = (double)gamma_f*q(j,k,i,5)/q(j,k,i,1);
                fnu  = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                // S*k/epsilon (first pass)
                if (ivmx == 6 || ivmx == 7) {
                    ske = std::sqrt(2.0*xis)/(turre(j,k,i,1)*re*0.09);
                } else if (ivmx == 8 || ivmx == 12 || ivmx == 14) {
                    ske = std::sqrt(2.0*xis)/(turre(j,k,i,1)*re);
                } else if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    ske = std::sqrt(2.0*xis)*turre(j,k,i,2)/(turre(j,k,i,1)*re);
                } else if (ivmx == 15) {
                    ske = std::sqrt(2.0*xis)*turre(j,k,i,2)*q(j,k,i,1)/(fnu*turre(j,k,i,1)*re);
                } else if (ivmx == 16) {
                    epsill = std::pow(0.09, 0.75)*std::pow(turre(j,k,i,2), 2.5)/turre(j,k,i,1);
                    ske    = std::sqrt(2.0*xis)*turre(j,k,i,2)/(epsill*re);
                }
                // P(k)*Lr/(u_inf)**3
                prod = pk/((double)xmach*(double)xmach*(double)xmach);
                // Fd blending function for DDES
                { double a1 = smin(j,k,i); dist = ccomplex_ns::ccabs(a1); }
                velterm = ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,2)*ux(j,k,i,2) + ux(j,k,i,3)*ux(j,k,i,3) +
                          ux(j,k,i,4)*ux(j,k,i,4) + ux(j,k,i,5)*ux(j,k,i,5) + ux(j,k,i,6)*ux(j,k,i,6) +
                          ux(j,k,i,7)*ux(j,k,i,7) + ux(j,k,i,8)*ux(j,k,i,8) + ux(j,k,i,9)*ux(j,k,i,9);
                rd = (vist3d(j,k,i)+fnu)/(q(j,k,i,1)*std::sqrt(velterm)*0.41*0.41*dist*dist*re);
                { double a1 = (8.0*rd)*(8.0*rd)*(8.0*rd); fd = 1.0 - ccomplex_ns::cctanh(a1); }
                // P(k)/epsilon
                if (ivmx == 6 || ivmx == 7) {
                    povere = pk/(turre(j,k,i,1)*turre(j,k,i,2)*re*0.09);
                } else if (ivmx == 8 || ivmx == 12 || ivmx == 14) {
                    povere = pk/(turre(j,k,i,1)*turre(j,k,i,2)*re);
                } else if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    povere = pk/(turre(j,k,i,1)*re);
                } else if (ivmx == 15) {
                    povere = pk*q(j,k,i,1)/(fnu*turre(j,k,i,1)*re);
                } else if (ivmx == 16) {
                    epsill = std::pow(0.09, 0.75)*std::pow(turre(j,k,i,2), 2.5)/turre(j,k,i,1);
                    povere = pk/(epsill*re);
                }
                // Turb Reynolds number
                if (ivmx == 6 || ivmx == 7) {
                    turbre = 4.0*q(j,k,i,1)*turre(j,k,i,2)/(9.0*fnu*turre(j,k,i,1)*0.09);
                } else if (ivmx == 8 || ivmx == 12 || ivmx == 14) {
                    turbre = 4.0*q(j,k,i,1)*turre(j,k,i,2)/(9.0*fnu*turre(j,k,i,1));
                } else if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    turbre = 4.0*q(j,k,i,1)*turre(j,k,i,2)*turre(j,k,i,2)/(9.0*fnu*turre(j,k,i,1));
                } else if (ivmx == 15) {
                    turbre = 4.0*q(j,k,i,1)*q(j,k,i,1)*turre(j,k,i,2)*turre(j,k,i,2)/(9.0*fnu*fnu*turre(j,k,i,1));
                } else if (ivmx == 16) {
                    epsill = std::pow(0.09, 0.75)*std::pow(turre(j,k,i,2), 2.5)/turre(j,k,i,1);
                    turbre = 4.0*q(j,k,i,1)*turre(j,k,i,2)*turre(j,k,i,2)/(9.0*fnu*epsill);
                }
                // Mean deformation
                def = (xis - wis)/(xis + wis + 1.0e-20);
                // Turbulence intensity
                if (ivmx >= 6) {
                    tint = std::sqrt(2.0*turre(j,k,i,2)/3.0)/(double)xmach;
                }



                // anisotropy tensor components
                b11 = t11/(2.0*q(j,k,i,1)*zk*re) - 0.33333333;
                if (ialph == 0) {
                    b22 = t22/(2.0*q(j,k,i,1)*zk*re) - 0.33333333;
                } else {
                    b22 = t33/(2.0*q(j,k,i,1)*zk*re) - 0.33333333;
                }
                if (ialph == 0) {
                    b33 = t33/(2.0*q(j,k,i,1)*zk*re) - 0.33333333;
                } else {
                    b33 = t22/(2.0*q(j,k,i,1)*zk*re) - 0.33333333;
                }
                if (ialph == 0) {
                    b12 = t12/(2.0*q(j,k,i,1)*zk*re);
                } else {
                    b12 = t13/(2.0*q(j,k,i,1)*zk*re);
                }
                if (ialph == 0) {
                    b13 = t13/(2.0*q(j,k,i,1)*zk*re);
                } else {
                    b13 = -t12/(2.0*q(j,k,i,1)*zk*re);
                }
                if (ialph == 0) {
                    b23 = t23/(2.0*q(j,k,i,1)*zk*re);
                } else {
                    b23 = -t23/(2.0*q(j,k,i,1)*zk*re);
                }
                // scalar invariant II
                if (i2d == 1) {
                    if (ialph == 0) {
                        xii = -0.5*(b11*b11 + b22*b22 + b33*b33 + 2.0*b13*b13);
                    } else {
                        xii = -0.5*(b11*b11 + b22*b22 + b33*b33 + 2.0*b12*b12);
                    }
                } else {
                    xii = -0.5*(b11*b11 + b22*b22 + b33*b33 + 2.0*b12*b12 + 2.0*b13*b13 + 2.0*b23*b23);
                }
                // scalar invariant III
                if (i2d == 1) {
                    if (ialph == 0) {
                        xiii = 0.333333*(b11*b11*b11 + b22*b22*b22 + b33*b33*b33 +
                                         3.0*b11*b13*b13 + 3.0*b33*b13*b13);
                    } else {
                        xiii = 0.333333*(b11*b11*b11 + b22*b22*b22 + b33*b33*b33 +
                                         3.0*b11*b12*b12 + 3.0*b22*b12*b12);
                    }
                } else {
                    xiii = 0.333333*(b11*b11*b11 + b22*b22*b22 + b33*b33*b33 +
                                     3.0*b11*b12*b12 + 3.0*b11*b13*b13 +
                                     3.0*b22*b12*b12 + 3.0*b22*b23*b23 +
                                     3.0*b33*b13*b13 + 3.0*b33*b23*b23 + 6.0*b12*b13*b23);
                }
                // scalar invariant F ("flatness parameter")
                ff = 1.0 + 9.0*xii + 27.0*xiii;
                // Q-criterion 2nd invariant of velocity gradient tensor
                if (i2d == 1) {
                    if (ialph == 0) {
                        qinvar = -0.5*(ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,9)*ux(j,k,i,9) +
                                       2.0*ux(j,k,i,3)*ux(j,k,i,7));
                    } else {
                        qinvar = -0.5*(ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,5)*ux(j,k,i,5) +
                                       2.0*ux(j,k,i,2)*ux(j,k,i,4));
                    }
                } else {
                    qinvar = -0.5*(ux(j,k,i,1)*ux(j,k,i,1) + ux(j,k,i,5)*ux(j,k,i,5) + ux(j,k,i,9)*ux(j,k,i,9) +
                                   2.0*ux(j,k,i,2)*ux(j,k,i,4) +
                                   2.0*ux(j,k,i,3)*ux(j,k,i,7) +
                                   2.0*ux(j,k,i,6)*ux(j,k,i,8));
                }
                // velocities by Uref
                uuu = q(j,k,i,2)/(double)xmach;
                if (ialph == 0) {
                    vvv = q(j,k,i,3)/(double)xmach;
                    www = q(j,k,i,4)/(double)xmach;
                } else {
                    vvv =  q(j,k,i,4)/(double)xmach;
                    www = -q(j,k,i,3)/(double)xmach;
                }
                // dissipation*Lr/(u_inf)**3
                if (ivmx == 6 || ivmx == 7) {
                    epsilon_loc = turre(j,k,i,1)*turre(j,k,i,2)*re*0.09/
                                  ((double)xmach*(double)xmach*(double)xmach);
                } else if (ivmx == 8 || ivmx == 12 || ivmx == 14) {
                    epsilon_loc = turre(j,k,i,1)*turre(j,k,i,2)*re/
                                  ((double)xmach*(double)xmach*(double)xmach);
                } else if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    epsilon_loc = turre(j,k,i,1)*re/((double)xmach*(double)xmach*(double)xmach);
                } else if (ivmx == 15) {
                    epsilon_loc = fnu*turre(j,k,i,1)*re/q(j,k,i,1)/
                                  ((double)xmach*(double)xmach*(double)xmach);
                } else if (ivmx == 16) {
                    epsill      = std::pow(0.09, 0.75)*std::pow(turre(j,k,i,2), 2.5)/turre(j,k,i,1);
                    epsilon_loc = epsill*re/((double)xmach*(double)xmach*(double)xmach);
                }
                // sigma is Sk/e, using Wallin-Johansson defn of S (second pass for ske)
                if (ivmx == 6 || ivmx == 7) {
                    ske = std::sqrt(xis/2.0)/(turre(j,k,i,1)*re*0.09);
                } else if (ivmx == 8 || ivmx == 12 || ivmx == 14) {
                    ske = std::sqrt(xis/2.0)/(turre(j,k,i,1)*re);
                } else if (ivmx == 9 || ivmx == 10 || ivmx == 11 || ivmx == 13) {
                    ske = std::sqrt(xis/2.0)*turre(j,k,i,2)/(turre(j,k,i,1)*re);
                } else if (ivmx == 15) {
                    ske = std::sqrt(xis/2.0)*turre(j,k,i,2)*q(j,k,i,1)/(fnu*turre(j,k,i,1)*re);
                } else if (ivmx == 16) {
                    epsill = std::pow(0.09, 0.75)*std::pow(turre(j,k,i,2), 2.5)/turre(j,k,i,1);
                    ske    = std::sqrt(xis/2.0)*turre(j,k,i,2)/(epsill*re);
                }



                // pressure-strain (only applies for EASMs)
                if (ivmx == 8 || ivmx == 9 || ivmx == 13 || ivmx == 14) {
                    s1kbk1 = s11*b11 + s12*b12 + s13*b13;
                    s2kbk2 = s12*b12 + s22*b22 + s23*b23;
                    s3kbk3 = s13*b13 + s23*b23 + s33*b33;
                    s1kbk2 = s11*b12 + s12*b22 + s13*b23;
                    s1kbk3 = s11*b13 + s12*b23 + s13*b33;
                    s2kbk3 = s12*b13 + s22*b23 + s23*b33;
                    b1ksk1 = b11*s11 + b12*s12 + b13*s13;
                    b2ksk2 = b12*s12 + b22*s22 + b23*s23;
                    b3ksk3 = b13*s13 + b23*s23 + b33*s33;
                    b1ksk2 = b11*s12 + b12*s22 + b13*s23;
                    b1ksk3 = b11*s13 + b12*s23 + b13*s33;
                    b2ksk3 = b12*s13 + b22*s23 + b23*s33;
                    sklblk = s11*b11 + s22*b22 + s33*b33 +
                             2.0*s12*b12 + 2.0*s13*b13 + 2.0*s23*b23;
                    w11 = 0.0;
                    w22 = 0.0;
                    w33 = 0.0;
                    w1kbk1 =  w11*b11 + w12*b12 + w13*b13;
                    w2kbk2 = -w12*b12 + w22*b22 + w23*b23;
                    w3kbk3 = -w13*b13 - w23*b23 + w33*b33;
                    w1kbk2 =  w11*b12 + w12*b22 + w13*b23;
                    w1kbk3 =  w11*b13 + w12*b23 + w13*b33;
                    w2kbk3 = -w12*b13 + w22*b23 + w23*b33;
                    b1kwk1 =  b11*w11 - b12*w12 - b13*w13;
                    b2kwk2 =  b12*w12 + b22*w22 - b23*w23;
                    b3kwk3 =  b13*w13 + b23*w23 + b33*w33;
                    b1kwk2 =  b11*w12 + b12*w22 - b13*w23;
                    b1kwk3 =  b11*w13 + b12*w23 + b13*w33;
                    b2kwk3 =  b12*w13 + b22*w23 + b23*w33;
                    arg = 4.0*sigma - 10.0;
                    if (ivmx == 14 && (ieasm_type == 3 || ieasm_type == 4)) {
                        { double a1 = arg; f = 0.5*(1.0 + ccomplex_ns::cctanh(a1)); }
                    } else {
                        f = 1.0;
                    }
                    c2star = 1.2 + (c2_e - 1.2)*f + 0.84*povere*(f - 1.0);
                    // pstrain*Lr/(u_inf)**3
                    r11slow = -(3.4*epsilon_loc + 1.8*prod)*b11;
                    r22slow = -(3.4*epsilon_loc + 1.8*prod)*b22;
                    r33slow = -(3.4*epsilon_loc + 1.8*prod)*b33;
                    r12slow = -(3.4*epsilon_loc + 1.8*prod)*b12;
                    r13slow = -(3.4*epsilon_loc + 1.8*prod)*b13;
                    r23slow = -(3.4*epsilon_loc + 1.8*prod)*b23;
                    r11rapid = (c2star*zk*s11 + zk*1.25*(s1kbk1+b1ksk1-0.66667*sklblk) +
                                zk*0.4*(w1kbk1-b1kwk1)) / ((double)xmach*(double)xmach*(double)xmach);
                    r22rapid = (c2star*zk*s22 + zk*1.25*(s2kbk2+b2ksk2-0.66667*sklblk) +
                                zk*0.4*(w2kbk2-b2kwk2)) / ((double)xmach*(double)xmach*(double)xmach);
                    r33rapid = (c2star*zk*s33 + zk*1.25*(s3kbk3+b3ksk3-0.66667*sklblk) +
                                zk*0.4*(w3kbk3-b3kwk3)) / ((double)xmach*(double)xmach*(double)xmach);
                    r12rapid = (c2star*zk*s12 + zk*1.25*(s1kbk2+b1ksk2) +
                                zk*0.4*(w1kbk2-b1kwk2)) / ((double)xmach*(double)xmach*(double)xmach);
                    r13rapid = (c2star*zk*s13 + zk*1.25*(s1kbk3+b1ksk3) +
                                zk*0.4*(w1kbk3-b1kwk3)) / ((double)xmach*(double)xmach*(double)xmach);
                    r23rapid = (c2star*zk*s23 + zk*1.25*(s2kbk3+b2ksk3) +
                                zk*0.4*(w2kbk3-b2kwk3)) / ((double)xmach*(double)xmach*(double)xmach);
                    r11 = r11slow + r11rapid;
                    r22 = r22slow + r22rapid;
                    r33 = r33slow + r33rapid;
                    r12 = r12slow + r12rapid;
                    r13 = r13slow + r13rapid;
                    r23 = r23slow + r23rapid;
                }

                // load appropriate data into single precision array
                xw(jw,kw,iw,1) = (double)(float)prod;
                xw(jw,kw,iw,2) = (double)(float)uwprime;
                xw(jw,kw,iw,3) = (double)(float)uuprime;
                xw(jw,kw,iw,4) = (double)(float)wwprime;
                xw(jw,kw,iw,5) = (double)(float)ske;
                if (ifunct >= 6)  xw(jw,kw,iw,6)  = (double)(float)turre(j,k,i,1);
                if (ifunct >= 7)  xw(jw,kw,iw,7)  = (double)(float)turre(j,k,i,2);
                if (ifunct >= 8)  xw(jw,kw,iw,8)  = (double)(float)vist3d(j,k,i);
                if (ifunct >= 9)  xw(jw,kw,iw,9)  = (double)(float)xii;
                if (ifunct >= 10) xw(jw,kw,iw,10) = (double)(float)xiii;
                if (ifunct >= 11) xw(jw,kw,iw,11) = (double)(float)povere;
                if (ifunct >= 12) xw(jw,kw,iw,12) = (double)(float)cmuv(j,k,i);
                if (ifunct >= 13) xw(jw,kw,iw,13) = (double)(float)b11;
                if (ifunct >= 14) xw(jw,kw,iw,14) = (double)(float)b22;
                if (ifunct >= 15) xw(jw,kw,iw,15) = (double)(float)b33;
                if (ifunct >= 16) xw(jw,kw,iw,16) = (double)(float)b13;
                if (ifunct >= 17) xw(jw,kw,iw,17) = (double)(float)uuu;
                if (ifunct >= 18) xw(jw,kw,iw,18) = (double)(float)www;
                if (ifunct >= 19) xw(jw,kw,iw,19) = (double)(float)r11;
                if (ifunct >= 20) xw(jw,kw,iw,20) = (double)(float)r33;
                if (ifunct >= 21) xw(jw,kw,iw,21) = (double)(float)r13;
                if (ifunct >= 22) xw(jw,kw,iw,22) = (double)(float)vvprime;
                if (ifunct >= 23) xw(jw,kw,iw,23) = (double)(float)uvprime;
                if (ifunct >= 24) xw(jw,kw,iw,24) = (double)(float)vwprime;
                if (ifunct >= 25) xw(jw,kw,iw,25) = (double)(float)turre(j,k,i,3);
                if (ifunct >= 26) xw(jw,kw,iw,26) = (double)(float)turre(j,k,i,4);
                if (ifunct >= 27) xw(jw,kw,iw,27) = (double)(float)ux(j,k,i,1);
                if (ifunct >= 28) xw(jw,kw,iw,28) = (double)(float)ux(j,k,i,3);
                if (ifunct >= 29) xw(jw,kw,iw,29) = (double)(float)ux(j,k,i,7);
                if (ifunct >= 30) xw(jw,kw,iw,30) = (double)(float)ux(j,k,i,9);
                if (ifunct >= 31) xw(jw,kw,iw,31) = (double)(float)smin(j,k,i);

            } // end j loop
        } // end k loop
    } // end i loop



    // output turb data in plot3d q file (unit 4)
    if (myid == myhost) {
        if (ibin == 0) {
            // formatted output
            if (ifunct == 0) {
                if (i2d == 0) {
                    // write(4,'(5e14.6)') xmachw,alphww,reuew,timew
                    {
                        FILE* f4 = fortran_get_unit(4);
                        fprintf(f4, "%14.6E%14.6E%14.6E%14.6E\n",
                                (float)xmachw, (float)alphww, (float)reuew, (float)timew);
                        // write(4,'(5e14.6)') ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
                        int cnt = 0;
                        for (m = 1; m <= 5; m++)
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f4, "%14.6E", (float)xw(j,k,i,m));
                                        if (++cnt % 5 == 0) fprintf(f4, "\n");
                                    }
                        if (cnt % 5 != 0) fprintf(f4, "\n");
                    }
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    fprintf(f4, "%14.6E%14.6E%14.6E%14.6E\n",
                            (float)xmachw, (float)alphww, (float)reuew, (float)timew);
                    int cnt = 0;
                    for (m = 1; m <= 4; m++)
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    fprintf(f4, "%14.6E", (float)xw(j,k,i,m));
                                    if (++cnt % 5 == 0) fprintf(f4, "\n");
                                }
                    if (cnt % 5 != 0) fprintf(f4, "\n");
                }
            } else {
                // ifunct != 0: write ifunct variables
                FILE* f4 = fortran_get_unit(4);
                int cnt = 0;
                for (m = 1; m <= ifunct; m++)
                    for (k = 1; k <= kw; k++)
                        for (j = 1; j <= jw; j++)
                            for (i = 1; i <= iw; i++) {
                                fprintf(f4, "%14.6E", (float)xw(j,k,i,m));
                                if (++cnt % 5 == 0) fprintf(f4, "\n");
                            }
                if (cnt % 5 != 0) fprintf(f4, "\n");
            }
        } else {
            // binary output
            if (ifunct == 0) {
                if (i2d == 0) {
                    FILE* f4 = fortran_get_unit(4);
                    float sv[4] = {(float)xmachw, (float)alphww, (float)reuew, (float)timew};
                    fwrite(sv, sizeof(float), 4, f4);
                    for (m = 1; m <= 5; m++)
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xw(j,k,i,m);
                                    fwrite(&v, sizeof(float), 1, f4);
                                }
                } else {
                    FILE* f4 = fortran_get_unit(4);
                    float sv[4] = {(float)xmachw, (float)alphww, (float)reuew, (float)timew};
                    fwrite(sv, sizeof(float), 4, f4);
                    for (m = 1; m <= 4; m++)
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xw(j,k,i,m);
                                    fwrite(&v, sizeof(float), 1, f4);
                                }
                }
            } else {
                FILE* f4 = fortran_get_unit(4);
                for (m = 1; m <= ifunct; m++)
                    for (k = 1; k <= kw; k++)
                        for (j = 1; j <= jw; j++)
                            for (i = 1; i <= iw; i++) {
                                float v = (float)xw(j,k,i,m);
                                fwrite(&v, sizeof(float), 1, f4);
                            }
            }
        }
    } // end if (myid == myhost)

    return;
}

} // namespace plot3t_ns
