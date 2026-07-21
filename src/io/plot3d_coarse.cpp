// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "plot3d_coarse.h"
#include "plot3d.h"
#include "cctogp.h"
#include "ccomplex.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace plot3d_coarse_ns {

void plot3d(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> xw, FortranArray4DRef<double> blank2, FortranArray3DRef<double> blank, FortranArray4DRef<double> xg, int& iflag, FortranArray3DRef<double> vist3d, int& iover, int& nblk, FortranArray1DRef<int> nmap, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& ifunc, int& iplot, int& jdw, int& kdw, int& idw, int& nplots, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay, int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d, int& nblock, FortranArray1DRef<int> nblkpt, FortranArray4DRef<double> xv, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, int& nset)
{
    plot3d_ns::plot3d(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
        q, qi0, qj0, qk0, x, y, z, xw, blank2, blank, xg, iflag, vist3d,
        iover, nblk, nmap, bcj, bck, bci, vj0, vk0, vi0, ifunc, iplot,
        jdw, kdw, idw, nplots, jdimg, kdimg, idimg, nblcg, jsg, ksg, isg,
        jeg, keg, ieg, ninter, iindex, intmax, nsub1, maxxe, nblkk, nbli,
        limblk, isva, nblon, mxbli, thetay, maxbl, maxgr, myid, myhost,
        mycomm, mblk2nd, inpl3d, nblock, nblkpt, xv, sj, sk, si, vol, nset);
}

void plot3d_coarse(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> xw, FortranArray4DRef<double> blank2, FortranArray3DRef<double> blank, FortranArray4DRef<double> xg, int& iflag, FortranArray3DRef<double> vist3d, int& iover, int& nblk, FortranArray1DRef<int> nmap, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, FortranArray4DRef<double> vi0, int& ifunc, int& iplot, int& jdw, int& kdw, int& idw, int& nplots, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay, int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d, int& nblock, FortranArray1DRef<int> nblkpt, FortranArray4DRef<double> xv, FortranArray4DRef<double> sj, FortranArray4DRef<double> sk, FortranArray4DRef<double> si, FortranArray3DRef<double> vol, int& nset)
{
    // COMMON block aliases
    int& ibin       = cmn_bin.ibin;
    int& iblnk      = cmn_bin.iblnk;
    float& gamma_f  = cmn_fluid.gamma;
    float& gm1_f    = cmn_fluid.gm1;
    float& gp1_f    = cmn_fluid.gp1;
    float& pr_f     = cmn_fluid2.pr;
    float& cbar_f   = cmn_fluid2.cbar;
    float& xmach_f  = cmn_info.xmach;
    float& alpha_f  = cmn_info.alpha;
    float& reue_f   = cmn_reyue.reue;
    float& tinf_f   = cmn_reyue.tinf;
    int*   ivisc    = cmn_reyue.ivisc;  // ivisc[3], 0-based C array
    int&   i2d      = cmn_twod.i2d;
    float& time_f   = cmn_unst.time;
    int*   iwf      = cmn_wallfun.iwf;  // iwf[3], 0-based C array
    int&   ialph    = cmn_igrdtyp.ialph;
    int&   lhdr     = cmn_moov.lhdr;
    int&   icallcrs = cmn_moovcrs2d.icallcrs;
    float& radtodeg_f = cmn_conversion.radtodeg;
    float& p0_f     = cmn_ivals.p0;

    // Promote COMMON floats to double for arithmetic
    double gamma    = (double)gamma_f;
    double gm1      = (double)gm1_f;
    double gp1      = (double)gp1_f;
    double pr       = (double)pr_f;
    double cbar     = (double)cbar_f;
    double xmach    = (double)xmach_f;
    double alpha    = (double)alpha_f;
    double reue     = (double)reue_f;
    double tinf     = (double)tinf_f;
    double time     = (double)time_f;
    double radtodeg = (double)radtodeg_f;
    double p0       = (double)p0_f;

    // Local variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;
    int jw, kw, iw;
    int j, k, i, l, m;
    double alphaw, xmachw, alphww, reuew, timew;
    double temp;
    int izero;
    int j11, j22, k11, k22, i11, i22;
    int jj, kk, ii;
    int nnn;
    int inter, lmax1, nbl, lst, lcoord, lend;
    int j21, k21, ie1, je1, ke1, is1, js1, ks1;
    int ll, mblk1, mblk2, mblk3, mblk4;
    int n, it, ir, itime, iti;
    int is, ie, js, je, ks, ke;
    int nblc, nblcc;
    int ibembed;
    int ldw, leq;
    double cpc;
    double term1, term2, term3;
    double q1, q2, q3, q4, q5, t1, xm1, edvis, pitot;
    int icp, isay;
    double id_d, id1_d;
    int id, id1, kd, kd1, kdx, jd, jd1, jdx, idx;
    double sgn;
    double q1k1, q2k1, q3k1, q4k1, q5k1, t1k1;
    double q1k2, q2k2, q3k2, q4k2, q5k2, t1k2;
    double q1j1, q2j1, q3j1, q4j1, q5j1, t1j1;
    double q1j2, q2j2, q3j2, q4j2, q5j2, t1j2;
    double q1i1, q2i1, q3i1, q4i1, q5i1, t1i1;
    double q1i2, q2i2, q3i2, q4i2, q5i2, t1i2;
    double dn, dx, dy, dz;
    double avgmut, emuka;
    double urel, vrel, wrel;
    double sk1, sk2, sk3, sj1, sj2, sj3, si1, si2, si3;
    double vnorm, upar, vpar, wpar;
    double Cf, cfx, cfy, cfz, Ch;
    double eps, tty1, t1m1;
    double pres1, temp1, yplus;
    int iold;
    int jv, kv, iv, jx, kx, ix;


    // initialize xw, xv and xg arrays
    jw = (j2-j1)/j3 + 1;
    kw = (k2-k1)/k3 + 1;
    iw = (i2-i1)/i3 + 1;
    for (j=1; j<=jw; j++) {
        for (k=1; k<=kw; k++) {
            for (i=1; i<=iw; i++) {
                for (l=1; l<=5; l++) {
                    xw(j,k,i,l) = 0.;
                    xv(j,k,i,l) = 0.;
                }
                for (l=1; l<=4; l++) {
                    xg(j,k,i,l) = 0.;
                }
            }
        }
    }

    // assign single precision scalars
    alphaw = radtodeg*(alpha+thetay(nblk));
    xmachw = xmach;
    alphww = alphaw;
    reuew  = reue;
    timew  = time;

    // load grid into first 3 locations of the single precision xg array
    iw = 0;
    for (i=i1; i<=i2; i+=i3) {
        iw = iw + 1;
        jw = 0;
        for (j=j1; j<=j2; j+=j3) {
            jw = jw + 1;
            kw = 0;
            for (k=k1; k<=k2; k+=k3) {
                kw = kw + 1;
                xg(jw,kw,iw,1) = x(j,k,i);
                xg(jw,kw,iw,2) = y(j,k,i);
                xg(jw,kw,iw,3) = z(j,k,i);
            }
        }
    }

    // set iblank (blank2) array
    if (iflag == 1 && iblnk > 0) {

        // assign default iblank (blank2) array
        for (i=1; i<=idim; i++)
            for (j=1; j<=jdim; j++)
                for (k=1; k<=kdim; k++) {
                    blank2(j,k,i,1) = 1.;
                    blank2(j,k,i,2) = 1.;
                }

        // zero out edges and corners if desired (izero > 0)
        izero = 0;

        if (izero > 0) {
            for (i=1; i<=idim; i++)
                for (j=1; j<=jdim; j+=jdim1)
                    for (k=1; k<=kdim; k+=kdim1)
                        blank2(j,k,i,1) = 0.;

            if (idim > 2) {
                for (i=1; i<=idim; i+=idim1) {
                    for (j=1; j<=jdim; j++)
                        blank2(j,1,i,1) = blank2(j,kdim,i,1) = 0.;
                    for (k=1; k<=kdim; k++)
                        blank2(1,k,i,1) = blank2(jdim,k,i,1) = 0.;
                }
            } else {
                for (j=1; j<=jdim; j++)
                    for (k=1; k<=kdim; k++)
                        blank2(j,k,idim,1) = 0.;
            }

            if (jdim == 2) {
                for (i=1; i<=idim; i++)
                    for (k=1; k<=kdim; k++)
                        blank2(jdim,k,i,1) = 0.;
            }

            if (kdim == 2) {
                for (i=1; i<=idim; i++)
                    for (k=1; k<=kdim; k++)
                        blank2(j,kdim,i,1) = 0.;
            }
        }

        // solid surface iblank (blank2) values (iblank=2 for solid surface)
        j11 = 1; j22 = jdim1;
        if (jdim == 2) { j11 = 1; j22 = 1; }
        k11 = 1; k22 = kdim1;
        if (kdim == 2) { k11 = 1; k22 = 1; }
        i11 = 1; i22 = idim1;
        if (idim == 2) { i11 = 1; i22 = 1; }

        i = 1;
        for (nnn=1; nnn<=2; nnn++) {
            for (j=j11; j<=j22; j++) {
                for (k=k11; k<=k22; k++) {
                    kk = std::min(k+1,k22);
                    jj = std::min(j+1,j22);
                    double a1=bci(j,k,nnn), a2=bci(jj,k,nnn), a3=bci(j,kk,nnn), a4=bci(jj,kk,nnn);
                    blank2(j,k,i,1) = 1. + ccomplex_ns::ccmax4(a1,a2,a3,a4);
                }
            }
            i = idim;
        }
        j = 1;
        for (nnn=1; nnn<=2; nnn++) {
            for (i=i11; i<=i22; i++) {
                for (k=k11; k<=k22; k++) {
                    kk = std::min(k+1,k22);
                    ii = std::min(i+1,i22);
                    double a1=bcj(k,i,nnn), a2=bcj(kk,i,nnn), a3=bcj(k,ii,nnn), a4=bcj(kk,ii,nnn);
                    blank2(j,k,i,1) = 1. + ccomplex_ns::ccmax4(a1,a2,a3,a4);
                }
            }
            j = jdim;
        }
        k = 1;
        for (nnn=1; nnn<=2; nnn++) {
            for (j=j11; j<=j22; j++) {
                for (i=i11; i<=i22; i++) {
                    jj = std::min(j+1,j22);
                    ii = std::min(i+1,i22);
                    double a1=bck(j,i,nnn), a2=bck(jj,i,nnn), a3=bck(j,ii,nnn), a4=bck(jj,ii,nnn);
                    blank2(j,k,i,1) = 1. + ccomplex_ns::ccmax4(a1,a2,a3,a4);
                }
            }
            k = kdim;
        }


        // iblank (blank2) array for patch interface boundaries
        if (std::abs(ninter) > 0) {
            for (inter=1; inter<=std::abs(ninter); inter++) {
                lmax1  = iindex(inter,1);
                nbl    = iindex(inter,lmax1+2);
                if (nbl != nblk) continue;
                lst    = iindex(inter,2*lmax1+5);
                lcoord = iindex(inter,2*lmax1+3)/10;
                lend   = iindex(inter,2*lmax1+3)-lcoord*10;
                j21    = iindex(inter,2*lmax1+6);
                j22    = iindex(inter,2*lmax1+7);
                k21    = iindex(inter,2*lmax1+8);
                k22    = iindex(inter,2*lmax1+9);

                if (lcoord == 1) {
                    if (lend == 1) i = 1;
                    if (lend == 2) i = idim;
                    if (jdim > 2 && kdim > 2) {
                        for (j=j21+1; j<=j22-1; j++) {
                            for (k=k21+1; k<=k22-1; k++) {
                                ll = lst + (j22-j21)*(k-1-k21) + (j-1-j21);
                                mblk1 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(k-k21-1) + (j-j21);
                                mblk2 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(k-k21) + (j-j21);
                                mblk3 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(k-k21) + (j-j21-1);
                                mblk4 = iindex(inter,nblkpt(ll)+1);
                                if (mblk1==mblk2 && mblk1==mblk3 && mblk1==mblk4)
                                    blank2(j,k,i,1) = -(float)nmap(mblk1);
                            }
                        }
                    } else if (jdim > 2) {
                        k = k21;
                        for (j=j21+1; j<=j22-1; j++) {
                            ll = lst + (j22-j21)*(k-k21) + (j-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22-j21)*(k-k21) + (j-j21-1);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1==mblk2) blank2(j,k,i,1) = -(float)nmap(mblk1);
                        }
                    } else if (kdim > 2) {
                        j = j21;
                        for (k=k21+1; k<=k22-1; k++) {
                            ll = lst + (j22-j21)*(k-k21-1) + (j-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22-j21)*(k-k21) + (j-j21);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1==mblk2) blank2(j,k,i,1) = -(float)nmap(mblk1);
                        }
                    }
                }

                if (lcoord == 2) {
                    if (lend == 1) j = 1;
                    if (lend == 2) j = jdim;
                    if (idim > 2 && kdim > 2) {
                        for (i=k21+1; i<=k22-1; i++) {
                            for (k=j21+1; k<=j22-1; k++) {
                                ll = lst + (j22-j21)*(i-1-k21) + (k-1-j21);
                                mblk1 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(i-k21-1) + (k-j21);
                                mblk2 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(i-k21) + (k-j21);
                                mblk3 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(i-k21) + (k-j21-1);
                                mblk4 = iindex(inter,nblkpt(ll)+1);
                                if (mblk1==mblk2 && mblk1==mblk3 && mblk1==mblk4)
                                    blank2(j,k,i,1) = -(float)nmap(mblk1);
                            }
                        }
                    } else if (idim > 2) {
                        k = j21;
                        for (i=k21+1; i<=k22-1; i++) {
                            ll = lst + (j22-j21)*(i-k21-1) + (k-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22-j21)*(i-k21) + (k-j21);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1==mblk2) blank2(j,k,i,1) = -(float)nmap(mblk1);
                        }
                    } else if (kdim > 2) {
                        i = k21;
                        for (k=j21+1; k<=j22-1; k++) {
                            ll = lst + (j22-j21)*(i-k21) + (k-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22-j21)*(i-k21) + (k-j21-1);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1==mblk2) blank2(j,k,i,1) = -(float)nmap(mblk1);
                        }
                    }
                }

                if (lcoord == 3) {
                    if (lend == 1) k = 1;
                    if (lend == 2) k = kdim;
                    if (idim > 2 && jdim > 2) {
                        for (i=k21+1; i<=k22-1; i++) {
                            for (j=j21+1; j<=j22-1; j++) {
                                ll = lst + (j22-j21)*(i-1-k21) + (j-1-j21);
                                mblk1 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(i-k21-1) + (j-j21);
                                mblk2 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(i-k21) + (j-j21);
                                mblk3 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22-j21)*(i-k21) + (j-j21-1);
                                mblk4 = iindex(inter,nblkpt(ll)+1);
                                if (mblk1==mblk2 && mblk1==mblk3 && mblk1==mblk4)
                                    blank2(j,k,i,1) = -(float)nmap(mblk1);
                            }
                        }
                    } else if (idim > 2) {
                        j = j21;
                        for (i=k21+1; i<=k22-1; i++) {
                            ll = lst + (j22-j21)*(i-k21-1) + (j-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22-j21)*(i-k21) + (j-j21);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1==mblk2) blank2(j,k,i,1) = -(float)nmap(mblk1);
                        }
                    } else if (jdim > 2) {
                        i = k21;
                        for (j=j21+1; j<=j22-1; j++) {
                            ll = lst + (j22-j21)*(i-k21) + (j-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22-j21)*(i-k21) + (j-j21-1);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1==mblk2) blank2(j,k,i,1) = -(float)nmap(mblk1);
                        }
                    }
                }
            } // end inter loop
        } // end if abs(ninter)>0


        // iblank (blank2) array for 1:1 interface boundaries
        if (nbli > 0) {
            for (n=1; n<=std::abs(nbli); n++) {
                if (nblon(n) >= 0) {
                    if (nblk==nblkk(1,n) || nblk==nblkk(2,n)) {
                        it = 1; ir = 2;
                        if (nblk == nblkk(2,n)) { it = 2; ir = 1; }

                        itime = 1;
                        if (nblkk(1,n) == nblkk(2,n)) itime = 2;
                        for (iti=1; iti<=itime; iti++) {
                            if (iti > 1) { it = 1; ir = 2; }

                            is = limblk(it,1,n);
                            ie = limblk(it,4,n);
                            js = limblk(it,2,n);
                            je = limblk(it,5,n);
                            ks = limblk(it,3,n);
                            ke = limblk(it,6,n);

                            if (isva(it,1,n)+isva(it,2,n) == 5) {
                                // i = constant interface
                                if (js > je) js = js+1;
                                if (js < je) je = je+1;
                                if (ks > ke) ks = ks+1;
                                if (ks < ke) ke = ke+1;
                                if (jdim == 2) { js = 1; je = 2; }
                                if (kdim == 2) { ks = 1; ke = 2; }
                            }
                            if (isva(it,1,n)+isva(it,2,n) == 4) {
                                // j = constant interface
                                if (is > ie) is = is+1;
                                if (is < ie) ie = ie+1;
                                if (ks > ke) ks = ks+1;
                                if (ks < ke) ke = ke+1;
                                if (idim == 2) { is = 1; ie = 2; }
                                if (kdim == 2) { ks = 1; ke = 2; }
                            }
                            if (isva(it,1,n)+isva(it,2,n) == 3) {
                                // k = constant interface
                                if (js > je) js = js+1;
                                if (js < je) je = je+1;
                                if (is > ie) is = is+1;
                                if (is < ie) ie = ie+1;
                                if (jdim == 2) { js = 1; je = 2; }
                                if (idim == 2) { is = 1; ie = 2; }
                            }

                            is1 = std::min(is,ie); ie1 = std::max(is,ie);
                            js1 = std::min(js,je); je1 = std::max(js,je);
                            ks1 = std::min(ks,ke); ke1 = std::max(ks,ke);

                            for (i=is1; i<=ie1; i++)
                                for (j=js1; j<=je1; j++)
                                    for (k=ks1; k<=ke1; k++)
                                        blank2(j,k,i,1) = -(float)nmap(nblkk(ir,n));
                        }
                    }
                }
            }
        }

        // iblank (blank2) array for embedded grids
        ibembed = 1;
        if (ibembed > 0) {
            for (nblc=1; nblc<=nblock; nblc++) {
                if (nblk == nblc) continue;
                nblcc = nblcg(nblc);
                if (nblcc == nblk) {
                    js = jsg(nblc);
                    if (js < jdimg(nblcc) && js > 1) js = js + 1;
                    ks = ksg(nblc);
                    if (ks < kdimg(nblcc) && ks > 1) ks = ks + 1;
                    is = isg(nblc);
                    if (is < idimg(nblcc) && is > 1) is = is + 1;
                    je = jeg(nblc);
                    if (je > 2 && je < jdimg(nblcc)) je = je - 1;
                    ke = keg(nblc);
                    if (ke > 2 && ke < kdimg(nblcc)) ke = ke - 1;
                    ie = ieg(nblc);
                    if (ie > 2 && ie < idimg(nblcc)) ie = ie - 1;
                    for (i=is; i<=ie; i++)
                        for (j=js; j<=je; j++)
                            for (k=ks; k<=ke; k++)
                                blank2(j,k,i,1) = 0.;
                }
            }
        }

        // iblank (blank2) array for overlapped grids
        if (iover == 1) {
            // interior of faces (i=1 and i=idim)
            // Correct loop: j from 2 to jdim1
            for (j=2; j<=jdim1; j++) {
                for (k=2; k<=kdim1; k++) {
                    double a1=blank(j,k,1), a2=blank(j-1,k,1), a3=blank(j,k-1,1), a4=blank(j-1,k-1,1);
                    blank2(j,k,1,2)    = ccomplex_ns::ccmin4(a1,a2,a3,a4);
                    double b1=blank(j,k,idim1), b2=blank(j-1,k,idim1), b3=blank(j,k-1,idim1), b4=blank(j-1,k-1,idim1);
                    blank2(j,k,idim,2) = ccomplex_ns::ccmin4(b1,b2,b3,b4);
                }
            }


            // edges and corners (i=1 and i=idim)
            for (m=1; m<=2; m++) {
                k  = 1; kk = 1;
                if (m == 2) { k = kdim; kk = kdim1; }
                for (j=2; j<=jdim1; j++) {
                    double a1=blank(j,kk,1), a2=blank(j-1,kk,1);
                    blank2(j,k,1,2)    = ccomplex_ns::ccmin(a1,a2);
                    double b1=blank(j,kk,idim1), b2=blank(j-1,kk,idim1);
                    blank2(j,k,idim,2) = ccomplex_ns::ccmin(b1,b2);
                }
            }
            for (m=1; m<=2; m++) {
                j  = 1; jj = 1;
                if (m == 2) { j = jdim; jj = jdim1; }
                for (k=2; k<=kdim1; k++) {
                    double a1=blank(jj,k,1), a2=blank(jj,k-1,1);
                    blank2(j,k,1,2)    = ccomplex_ns::ccmin(a1,a2);
                    double b1=blank(jj,k,idim1), b2=blank(jj,k-1,idim1);
                    blank2(j,k,idim,2) = ccomplex_ns::ccmin(b1,b2);
                }
            }
            blank2(jdim,kdim,1,2)    = blank2(jdim-1,kdim-1,1,2);
            blank2(1,1,1,2)          = blank2(2,2,1,2);
            blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim-1,idim,2);
            blank2(1,1,idim,2)       = blank2(2,2,idim,2);

            // interior of faces (j=1 and j=jdim)
            for (i=2; i<=idim1; i++) {
                for (k=2; k<=kdim1; k++) {
                    double a1=blank(1,k,i), a2=blank(1,k,i-1), a3=blank(1,k-1,i), a4=blank(1,k-1,i-1);
                    blank2(1,k,i,2)    = ccomplex_ns::ccmin4(a1,a2,a3,a4);
                    double b1=blank(jdim1,k,i), b2=blank(jdim1,k,i-1), b3=blank(jdim1,k-1,i), b4=blank(jdim1,k-1,i-1);
                    blank2(jdim,k,i,2) = ccomplex_ns::ccmin4(b1,b2,b3,b4);
                }
            }

            // edges and corners (j=1 and j=jdim)
            for (m=1; m<=2; m++) {
                k  = 1; kk = 1;
                if (m == 2) { k = kdim; kk = kdim1; }
                for (i=2; i<=idim1; i++) {
                    double a1=blank(1,kk,i), a2=blank(1,kk,i-1);
                    blank2(1,k,i,2)    = ccomplex_ns::ccmin(a1,a2);
                    double b1=blank(jdim1,kk,i), b2=blank(jdim1,kk,i-1);
                    blank2(jdim,k,i,2) = ccomplex_ns::ccmin(b1,b2);
                }
            }
            for (m=1; m<=2; m++) {
                i  = 1; ii = 1;
                if (m == 2) { i = idim; ii = idim1; }
                for (k=2; k<=kdim1; k++) {
                    double a1=blank(1,k,ii), a2=blank(1,k-1,ii);
                    blank2(1,k,i,2)    = ccomplex_ns::ccmin(a1,a2);
                    double b1=blank(jdim1,k,ii), b2=blank(jdim1,k-1,ii);
                    blank2(jdim,k,i,2) = ccomplex_ns::ccmin(b1,b2);
                }
            }
            blank2(1,kdim,idim,2)    = blank2(1,kdim-1,idim-1,2);
            blank2(1,1,1,2)          = blank2(1,2,2,2);
            blank2(jdim,kdim,idim,2) = blank2(jdim,kdim-1,idim-1,2);
            blank2(jdim,1,1,2)       = blank2(jdim,2,2,2);

            // interior of faces (k=1 and k=kdim)
            for (i=2; i<=idim1; i++) {
                for (j=2; j<=jdim1; j++) {
                    double a1=blank(j,1,i), a2=blank(j,1,i-1), a3=blank(j-1,1,i), a4=blank(j-1,1,i-1);
                    blank2(j,1,i,2)    = ccomplex_ns::ccmin4(a1,a2,a3,a4);
                    double b1=blank(j,kdim1,i), b2=blank(j,kdim1,i-1), b3=blank(j-1,kdim1,i), b4=blank(j-1,kdim1,i-1);
                    blank2(j,kdim,i,2) = ccomplex_ns::ccmin4(b1,b2,b3,b4);
                }
            }

            // edges and corners (k=1 and k=kdim)
            for (m=1; m<=2; m++) {
                j  = 1; jj = 1;
                if (m == 2) { j = jdim; jj = jdim1; }
                for (i=2; i<=idim1; i++) {
                    double a1=blank(jj,1,i), a2=blank(jj,1,i-1);
                    blank2(j,1,i,2)    = ccomplex_ns::ccmin(a1,a2);
                    double b1=blank(jj,kdim1,i), b2=blank(jj,kdim1,i-1);
                    blank2(j,kdim,i,2) = ccomplex_ns::ccmin(b1,b2);
                }
            }
            for (m=1; m<=2; m++) {
                i  = 1; ii = 1;
                if (m == 2) { i = idim; ii = idim1; }
                for (j=2; j<=jdim1; j++) {
                    double a1=blank(j,1,ii), a2=blank(j-1,1,ii);
                    blank2(j,1,i,2)    = ccomplex_ns::ccmin(a1,a2);
                    double b1=blank(j,kdim1,ii), b2=blank(j-1,kdim1,ii);
                    blank2(j,kdim,i,2) = ccomplex_ns::ccmin(b1,b2);
                }
            }
            blank2(jdim,1,idim,2)    = blank2(jdim-1,1,idim-1,2);
            blank2(1,1,1,2)          = blank2(2,1,2,2);
            blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim,idim-1,2);
            blank2(1,kdim,1,2)       = blank2(2,kdim,2,2);

            // interior cells
            for (i=2; i<=idim1; i++) {
                for (j=2; j<=jdim1; j++) {
                    for (k=2; k<=kdim1; k++) {
                        double a1=blank(j,k,i), a2=blank(j-1,k,i), a3=blank(j,k-1,i), a4=blank(j-1,k-1,i);
                        double a5=blank(j,k,i-1), a6=blank(j-1,k,i-1), a7=blank(j,k-1,i-1), a8=blank(j-1,k-1,i-1);
                        blank2(j,k,i,2) = ccomplex_ns::ccmin8(a1,a2,a3,a4,a5,a6,a7,a8);
                    }
                }
            }

            // combining topology and overset grid blankings
            for (i=1; i<=idim; i++)
                for (j=1; j<=jdim; j++)
                    for (k=1; k<=kdim; k++) {
                        double a1=blank2(j,k,i,1), a2=blank2(j,k,i,2);
                        blank2(j,k,i,1) = ccomplex_ns::ccmin(a1,a2);
                    }

        } // end if (iover==1)

    } // end if (iflag==1 && iblnk>0)


    if (iflag == 1) {
        // plot3d data

        if (lhdr > 0 && myid == myhost) {
            if (i2d == 1) {
                fortran_write_unit(11, "writing plot3d file for JDIM X KDIM =%5d x%5d grid\n", jdim, kdim);
            } else {
                fortran_write_unit(11, "writing plot3d file for IDIM X JDIM X KDIM =%5d x %5d x %5d grid\n", idim, jdim, kdim);
            }
            if (inpl3d(iplot,2) == 0) {
                fortran_write_unit(11, "   plot3dg file is an xyz file at grid points\n");
                fortran_write_unit(11, "   plot3dq file is a    q file at grid points\n");
            }
            if (ifunc == -5) {
                fortran_write_unit(11, "   plot3dg file is an xyz file at grid points\n");
                fortran_write_unit(11, "   plot3dq file is a function file of Cp at grid points\n");
            }
            if (ifunc == -6) {
                fortran_write_unit(11, "   plot3dg file is an xyz file at grid points\n");
                fortran_write_unit(11, "   plot3dq file is a function file of P/Pinf at grid points\n");
            }
            if (i2d == 1) {
                if (iblnk == 0) {
                    fortran_write_unit(11, "   plot3d files to be read with /mgrid/2d qualifiers\n");
                } else {
                    fortran_write_unit(11, "   plot3d files to be read with /mgrid/blank/2d qualifiers\n");
                }
            } else {
                if (iblnk == 0) {
                    fortran_write_unit(11, "   plot3d files to be read with /mgrid qualifiers\n");
                } else {
                    fortran_write_unit(11, "   plot3d files to be read with /mgrid/blank qualifiers\n");
                }
            }
        }

        // load blank2 into 4th location of the single precision xg array
        iw = 0;
        for (i=i1; i<=i2; i+=i3) {
            iw = iw + 1;
            jw = 0;
            for (j=j1; j<=j2; j+=j3) {
                jw = jw + 1;
                kw = 0;
                for (k=k1; k<=k2; k+=k3) {
                    kw = kw + 1;
                    xg(jw,kw,iw,4) = blank2(j,k,i,1);
                }
            }
        }

        if (myid == myhost && icallcrs == 0) {
            // ialph > 0 for a grid that was read in plot3d format with alpha measured
            // in the xy plane (TLNS3D convention)
            if (ialph != 0 && i2d != 1) {
                for (i=1; i<=iw; i++) {
                    for (j=1; j<=jw; j++) {
                        for (k=1; k<=kw; k++) {
                            temp        = xg(j,k,i,2);
                            xg(j,k,i,2) = xg(j,k,i,3);
                            xg(j,k,i,3) = -temp;
                        }
                    }
                }
            }

            // output grid
            FILE* f101 = fortran_get_unit(101);
            if (ibin == 0) {
                if (i2d == 0) {
                    if (iblnk == 0) {
                        // write(101,'(5e14.6)') (((xg(j,k,i,1),...),xg(2),...,xg(3)...)
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,1));
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,2));
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,3));
                        fprintf(f101, "\n");
                    } else {
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,1));
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,2));
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,3));
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)(int)xg(j,k,i,4));
                        fprintf(f101, "\n");
                    }
                } else {
                    if (iblnk == 0) {
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,1));
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,3));
                        fprintf(f101, "\n");
                    } else {
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,1));
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)xg(j,k,i,3));
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) fprintf(f101, "%14.6E", (float)(int)xg(j,k,i,4));
                        fprintf(f101, "\n");
                    }
                }
            } else {
                // binary
                if (i2d == 0) {
                    if (iblnk == 0) {
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,1); fwrite(&v,sizeof(float),1,f101); }
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,2); fwrite(&v,sizeof(float),1,f101); }
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,3); fwrite(&v,sizeof(float),1,f101); }
                    } else {
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,1); fwrite(&v,sizeof(float),1,f101); }
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,2); fwrite(&v,sizeof(float),1,f101); }
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,3); fwrite(&v,sizeof(float),1,f101); }
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { int v=(int)xg(j,k,i,4); fwrite(&v,sizeof(int),1,f101); }
                    }
                } else {
                    if (iblnk == 0) {
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,1); fwrite(&v,sizeof(float),1,f101); }
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,3); fwrite(&v,sizeof(float),1,f101); }
                    } else {
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,1); fwrite(&v,sizeof(float),1,f101); }
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { float v=(float)xg(j,k,i,3); fwrite(&v,sizeof(float),1,f101); }
                        for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++) { int v=(int)xg(j,k,i,4); fwrite(&v,sizeof(int),1,f101); }
                    }
                }
            }
        } // end if (myid==myhost && icallcrs==0)


        // determine q values at grid points and load into xv array
        jdw = (j2-j1)/j3 + 1;
        kdw = (k2-k1)/k3 + 1;
        idw = (i2-i1)/i3 + 1;
        ldw = 5;

        cctogp_ns::cctogp(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                          q, qi0, qj0, qk0, jdw, kdw, idw, xv, ldw);

        jw = (j2-j1)/j3 + 1;
        kw = (k2-k1)/k3 + 1;
        iw = (i2-i1)/i3 + 1;

        if (ifunc == 0) {
            // load solution (q) in xv array
            for (i=1; i<=iw; i++)
                for (k=1; k<=kw; k++)
                    for (j=1; j<=jw; j++) {
                        xv(j,k,i,5) = xv(j,k,i,5)/gm1 + 0.5*(xv(j,k,i,2)*xv(j,k,i,2)
                                    + xv(j,k,i,3)*xv(j,k,i,3) + xv(j,k,i,4)*xv(j,k,i,4))*xv(j,k,i,1);
                        xv(j,k,i,2) = xv(j,k,i,1)*xv(j,k,i,2);
                        xv(j,k,i,3) = xv(j,k,i,1)*xv(j,k,i,3);
                        xv(j,k,i,4) = xv(j,k,i,1)*xv(j,k,i,4);
                    }
        } else {
            // load desired function in xv array
            if (ifunc == -5) {
                cpc = 2.e0/(gamma*xmach*xmach);
                for (i=1; i<=iw; i++)
                    for (k=1; k<=kw; k++)
                        for (j=1; j<=jw; j++)
                            xv(j,k,i,1) = (xv(j,k,i,5)/p0-1)*cpc;
            }
            if (ifunc == -6) {
                for (i=1; i<=iw; i++)
                    for (k=1; k<=kw; k++)
                        for (j=1; j<=jw; j++)
                            xv(j,k,i,1) = xv(j,k,i,5)/p0;
            }
        }

        // load xv into single precision xw array
        leq = 5;
        if (ifunc != 0) leq = 1;

        for (l=1; l<=leq; l++)
            for (i=1; i<=iw; i++)
                for (k=1; k<=kw; k++)
                    for (j=1; j<=jw; j++)
                        xw(j,k,i,l) = xv(j,k,i,l);

        if (myid == myhost) {
            // output solution
            FILE* f92 = fortran_get_unit(92);
            if (ifunc == 0) {
                if (ibin == 0) {
                    fprintf(f92, "%14.6E%14.6E%14.6E%14.6E\n", (float)xmachw, (float)alphww, (float)reuew, (float)timew);
                    if (i2d == 0) {
                        if (ialph == 0) {
                            for (m=1; m<=5; m++)
                                for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                    fprintf(f92, "%14.6E", (float)xw(j,k,i,m));
                            fprintf(f92, "\n");
                        } else {
                            for (i=1; i<=iw; i++) for (j=1; j<=jw; j++) for (k=1; k<=kw; k++)
                                xw(j,k,i,3) = -xw(j,k,i,3);
                            for (m=1; m<=2; m++)
                                for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                    fprintf(f92, "%14.6E", (float)xw(j,k,i,m));
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                fprintf(f92, "%14.6E", (float)xw(j,k,i,4));
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                fprintf(f92, "%14.6E", (float)xw(j,k,i,3));
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                fprintf(f92, "%14.6E", (float)xw(j,k,i,5));
                            fprintf(f92, "\n");
                        }
                    } else {
                        for (m=1; m<=2; m++)
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                fprintf(f92, "%14.6E", (float)xw(j,k,i,m));
                        for (m=4; m<=5; m++)
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                fprintf(f92, "%14.6E", (float)xw(j,k,i,m));
                        fprintf(f92, "\n");
                    }
                } else {
                    // binary
                    { float v=(float)xmachw; fwrite(&v,sizeof(float),1,f92); }
                    { float v=(float)alphww;  fwrite(&v,sizeof(float),1,f92); }
                    { float v=(float)reuew;   fwrite(&v,sizeof(float),1,f92); }
                    { float v=(float)timew;   fwrite(&v,sizeof(float),1,f92); }
                    if (i2d == 0) {
                        if (ialph == 0) {
                            for (m=1; m<=5; m++)
                                for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                    { float v=(float)xw(j,k,i,m); fwrite(&v,sizeof(float),1,f92); }
                        } else {
                            for (i=1; i<=iw; i++) for (j=1; j<=jw; j++) for (k=1; k<=kw; k++)
                                xw(j,k,i,3) = -xw(j,k,i,3);
                            for (m=1; m<=2; m++)
                                for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                    { float v=(float)xw(j,k,i,m); fwrite(&v,sizeof(float),1,f92); }
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                { float v=(float)xw(j,k,i,4); fwrite(&v,sizeof(float),1,f92); }
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                { float v=(float)xw(j,k,i,3); fwrite(&v,sizeof(float),1,f92); }
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                { float v=(float)xw(j,k,i,5); fwrite(&v,sizeof(float),1,f92); }
                        }
                    } else {
                        for (m=1; m<=2; m++)
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                { float v=(float)xw(j,k,i,m); fwrite(&v,sizeof(float),1,f92); }
                        for (m=4; m<=5; m++)
                            for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                                { float v=(float)xw(j,k,i,m); fwrite(&v,sizeof(float),1,f92); }
                    }
                }
            } else {
                // function file
                if (ibin == 0) {
                    for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                        fprintf(f92, "%14.6E", (float)xw(j,k,i,1));
                    fprintf(f92, "\n");
                } else {
                    for (k=1; k<=kw; k++) for (j=1; j<=jw; j++) for (i=1; i<=iw; i++)
                        { float v=(float)xw(j,k,i,1); fwrite(&v,sizeof(float),1,f92); }
                }
            }
        } // end if (myid==myhost)

    } // end if (iflag==1)


    if (iflag == 2) {
        // print out data

        if (lhdr > 0 && myid == myhost) {
            fortran_write_unit(11, "writing printout file for IDIM X JDIM X KDIM =%5d x %5d x %5d grid\n", idim, jdim, kdim);
        }

        // determine q values at grid points and load into single precision array
        jdw = (j2-j1)/j3 + 1;
        kdw = (k2-k1)/k3 + 1;
        idw = (i2-i1)/i3 + 1;
        ldw = 5;
        cctogp_ns::cctogp(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                          q, qi0, qj0, qk0, jdw, kdw, idw, xv, ldw);
        for (l=1; l<=ldw; l++)
            for (i=1; i<=idw; i++)
                for (k=1; k<=kdw; k++)
                    for (j=1; j<=jdw; j++)
                        xw(j,k,i,l) = xv(j,k,i,l);

        // determine eddy viscosity values at grid points
        if (ivisc[2] > 1 || ivisc[1] > 1 || ivisc[0] > 1) {
            ldw = 1;
            // cctogp with vist3d (3D) as 4D dum, and xv(1,1,1,4) as output
            FortranArray4DRef<double> vist3d_4d(&vist3d(1,1,1), jdim, kdim, idim, 1);
            FortranArray4DRef<double> xv_4d(&xv(1,1,1,4), jdw, kdw, idw, 1);
            cctogp_ns::cctogp(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                              vist3d_4d, vi0, vj0, vk0, jdw, kdw, idw, xv_4d, ldw);
            for (i=1; i<=idw; i++)
                for (k=1; k<=kdw; k++)
                    for (j=1; j<=jdw; j++)
                        xg(j,k,i,4) = xv(j,k,i,4);
        }

        // old way for eddy viscosity (iold=0, kept for reference)
        iold = 0;
        if (iold > 0) {
            if (ivisc[2] > 1 || ivisc[1] > 1 || ivisc[0] > 1) {
                iw = 0;
                for (i=i1; i<=i2; i+=i3) {
                    iw = iw + 1;
                    jw = 0;
                    for (j=j1; j<=j2; j+=j3) {
                        jw = jw + 1;
                        kw = 0;
                        for (k=k1; k<=k2; k+=k3) {
                            kw = kw + 1;
                            jv = j; kv = k; iv = i;
                            jv = std::min(jdim-1, jv);
                            kv = std::min(kdim-1, kv);
                            iv = std::min(idim-1, iv);
                            jx = j-1; kx = k-1; ix = i-1;
                            jx = std::max(1, jx);
                            kx = std::max(1, kx);
                            ix = std::max(1, ix);
                            xg(jw,kw,iw,4) = 0.125*(vist3d(jv,kv,iv) + vist3d(jv,kv,ix)
                                           + vist3d(jx,kv,iv)   + vist3d(jx,kv,ix)
                                           + vist3d(jv,kx,iv)   + vist3d(jv,kx,ix)
                                           + vist3d(jx,kx,iv)   + vist3d(jx,kx,ix));
                        }
                    }
                }
            }
        }

        if (myid == myhost) {
            // icp=1 prints out cp's; icp=0 prints out pitot pressures
            icp = 1;

            FILE* f17 = fortran_get_unit(17);
            if (icp == 1) {
                fprintf(f17, "\n   I   J   K         X                Y                Z"
                        "              U/Uinf           V/Vinf           W/Winf"
                        "           P/Pinf           T/Tinf            MACH              cp"
                        "            tur. vis.\n");
            } else {
                fprintf(f17, "\n   I   J   K         X                Y                Z"
                        "              U/Uinf           V/Vinf           W/Winf"
                        "           P/Pinf           T/Tinf            MACH          pitotp"
                        "            tur. vis.\n");
            }

            term3 = 1./std::pow(1.+0.5*gm1*xmach*xmach, gamma/gm1);
            iw = 0;
            for (i=i1; i<=i2; i+=i3) {
                iw = iw + 1;
                kw = 0;
                for (k=k1; k<=k2; k+=k3) {
                    kw = kw + 1;
                    jw = 0;
                    for (j=j1; j<=j2; j+=j3) {
                        jw = jw + 1;
                        q1 = xw(jw,kw,iw,1);
                        q2 = xw(jw,kw,iw,2)/xmach;
                        q3 = xw(jw,kw,iw,3)/xmach;
                        q4 = xw(jw,kw,iw,4)/xmach;
                        q5 = gamma*xw(jw,kw,iw,5);
                        t1 = q5/q1;
                        xm1 = std::sqrt(xmach*xmach*(q2*q2+q3*q3+q4*q4)/t1);

                        // turbulent viscosity
                        edvis = xg(jw,kw,iw,4);

                        // cp or pitot pressure
                        if (icp == 1) {
                            pitot = 2.*(q5-1.)/(gamma*xmach*xmach);
                        } else {
                            if ((float)xm1 > 1.0f) {
                                term1 = std::pow(0.5*gp1*xm1*xm1, gamma/gm1);
                                term2 = std::pow(2.*gamma*xm1*xm1/gp1 - gm1/gp1, 1./gm1);
                                pitot = q5*term1*term3/term2;
                            } else {
                                term1 = std::pow(1.0+0.5*gm1*xm1*xm1, gamma/gm1);
                                pitot = q5*term1*term3;
                            }
                        }
                        if (ialph == 0) {
                            fprintf(f17, "%4d%4d%4d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                i, j, k,
                                (float)xg(jw,kw,iw,1), (float)xg(jw,kw,iw,2), (float)xg(jw,kw,iw,3),
                                (float)q2, (float)q3, (float)q4, (float)q5, (float)t1, (float)xm1,
                                (float)pitot, (float)edvis);
                        } else {
                            xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                            q3 = -q3;
                            fprintf(f17, "%4d%4d%4d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                i, j, k,
                                (float)xg(jw,kw,iw,1), (float)xg(jw,kw,iw,3), (float)xg(jw,kw,iw,2),
                                (float)q2, (float)q4, (float)q3, (float)q5, (float)t1, (float)xm1,
                                (float)pitot, (float)edvis);
                        }
                    }
                }
            }
        } // end if (myid==myhost)


        if (ivisc[2] != 0) {
            // output viscous-flow data on k=1 and/or k=kdim surfaces
            iw = 0;
            for (i=i1; i<=i2; i+=i3) {
                iw = iw + 1;
                id   = i;
                id1  = id-1;
                if (id1 < 1)    id1 = 1;
                if (id > idim1) id  = idim1;
                kw = 0;
                for (k=k1; k<=k2; k+=k3) {
                    kw = kw + 1;
                    jw = 0;
                    if (k != 1 && k != kdim) continue;
                    for (j=j1; j<=j2; j+=j3) {
                        jw = jw + 1;
                        if (j == 1 || j == jdim) continue;
                        if (k == 1) {
                            kd  = 1;
                            kd1 = kd+1;
                            kdx = 1;
                            m   = 2;
                            sgn = 1.;
                        } else {
                            kd  = kdim1;
                            kd1 = kd-1;
                            kdx = kdim;
                            m   = 4;
                            sgn = -1.;
                        }

                        // wall value - average in J and I
                        q1k1 =       .25*(qk0(j,id,1,m)  +qk0(j-1,id,1,m)
                                        + qk0(j,id1,1,m) +qk0(j-1,id1,1,m));
                        q2k1 =       .25*(qk0(j,id,2,m)  +qk0(j-1,id,2,m)
                                        + qk0(j,id1,2,m) +qk0(j-1,id1,2,m));
                        q3k1 =       .25*(qk0(j,id,3,m)  +qk0(j-1,id,3,m)
                                        + qk0(j,id1,3,m) +qk0(j-1,id1,3,m));
                        q4k1 =       .25*(qk0(j,id,4,m)  +qk0(j-1,id,4,m)
                                        + qk0(j,id1,4,m) +qk0(j-1,id1,4,m));
                        q5k1 = gamma*.25*(qk0(j,id,5,m)  +qk0(j-1,id,5,m)
                                        + qk0(j,id1,5,m) +qk0(j-1,id1,5,m));
                        t1k1 = q5k1/q1k1;

                        // first cell center location - average in J and I
                        q1k2 =       .25*(q(j,kd,id,1)   +q(j-1,kd,id,1)
                                        + q(j,kd,id1,1)  +q(j-1,kd,id1,1));
                        q2k2 =       .25*(q(j,kd,id,2)   +q(j-1,kd,id,2)
                                        + q(j,kd,id1,2)  +q(j-1,kd,id1,2));
                        q3k2 =       .25*(q(j,kd,id,3)   +q(j-1,kd,id,3)
                                        + q(j,kd,id1,3)  +q(j-1,kd,id1,3));
                        q4k2 =       .25*(q(j,kd,id,4)   +q(j-1,kd,id,4)
                                        + q(j,kd,id1,4)  +q(j-1,kd,id1,4));
                        q5k2 = gamma*.25*(q(j,kd,id,5)   +q(j-1,kd,id,5)
                                        + q(j,kd,id1,5)  +q(j-1,kd,id1,5));
                        t1k2 = q5k2/q1k2;

                        if (sk(j  ,kdx,id ,4) == 0. ||
                            sk(j-1,kdx,id ,4) == 0. ||
                            sk(j  ,kdx,id1,4) == 0. ||
                            sk(j-1,kdx,id1,4) == 0.) {
                            dx = x(j,kd+1,i)-x(j,kd,i);
                            dy = y(j,kd+1,i)-y(j,kd,i);
                            dz = z(j,kd+1,i)-z(j,kd,i);
                            dn = std::sqrt(dx*dx+dy*dy+dz*dz);
                        } else {
                            dn = (vol(j  ,kd,id )/sk(j  ,kdx,id ,4)+
                                  vol(j-1,kd,id )/sk(j-1,kdx,id ,4)+
                                  vol(j  ,kd,id1)/sk(j  ,kdx,id1,4)+
                                  vol(j-1,kd,id1)/sk(j-1,kdx,id1,4))/4.;
                        }

                        // Get turb viscosity at wall
                        if (iwf[2] == 1) {
                            avgmut = .25*(vk0(j,id,1,m)  +vk0(j-1,id,1,m)
                                        + vk0(j,id1,1,m) +vk0(j-1,id1,1,m));
                        } else {
                            avgmut = 0.;
                        }
                        emuka = std::pow(t1k1,1.5)*((1.0+cbar/tinf)/(t1k1+cbar/tinf));

                        // Use component of velocity parallel to wall
                        urel = q2k2-q2k1;
                        vrel = q3k2-q3k1;
                        wrel = q4k2-q4k1;
                        sk1  = (sk(j  ,kdx,id ,1)+sk(j-1,kdx,id ,1)+
                                sk(j  ,kdx,id1,1)+sk(j-1,kdx,id1,1))/4.;
                        sk2  = (sk(j  ,kdx,id ,2)+sk(j-1,kdx,id ,2)+
                                sk(j  ,kdx,id1,2)+sk(j-1,kdx,id1,2))/4.;
                        sk3  = (sk(j  ,kdx,id ,3)+sk(j-1,kdx,id ,3)+
                                sk(j  ,kdx,id1,3)+sk(j-1,kdx,id1,3))/4.;
                        vnorm = (urel*sk1 + vrel*sk2 + wrel*sk3)*sgn;
                        upar = urel-vnorm*sk1*sgn;
                        vpar = vrel-vnorm*sk2*sgn;
                        wpar = wrel-vnorm*sk3*sgn;
                        Cf   = 2.*std::sqrt(upar*upar + vpar*vpar + wpar*wpar);
                        Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn;
                        if ((float)q2k2 < 0.) Cf = -Cf;
                        // Cf vector
                        cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn;
                        cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn;
                        cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn;

                        eps  = 1.0e-6;
                        tty1 = 1.0 + gm1*0.5*xmach*xmach;
                        t1m1 = t1k1-tty1;
                        if (std::abs((float)t1m1) <= (float)eps) {
                            Ch = 999.99999;
                        } else {
                            Ch = 2.*(t1k2-t1k1);
                            Ch = (emuka+avgmut)/(reue*pr*(t1k1-tty1))*Ch/dn;
                        }

                        pres1 = q5k1;
                        temp1 = t1k1;
                        yplus = 0.;
                        if (std::abs((float)Cf) > 0.) {
                            double cfabs_arg = Cf*0.5/q1k1;
                            yplus = dn*reue*q1k1*std::sqrt(ccomplex_ns::ccabs(cfabs_arg))/emuka;
                        }

                        // store data in single precision arrays
                        xg(jw,kw,iw,1) = x(j,k,i);
                        xg(jw,kw,iw,2) = y(j,k,i);
                        xg(jw,kw,iw,3) = z(j,k,i);
                        xg(jw,kw,iw,4) = dn;

                        xw(jw,kw,iw,1) = pres1;
                        xw(jw,kw,iw,2) = temp1;
                        xw(jw,kw,iw,3) = Cf;
                        xw(jw,kw,iw,4) = Ch;
                        xw(jw,kw,iw,5) = yplus;
                        xw(jw,kw,iw,6) = cfx;
                        xw(jw,kw,iw,7) = cfy;
                        xw(jw,kw,iw,8) = cfz;
                    } // j loop
                } // k loop
            } // i loop

            // write printout data to file
            if (myid == myhost) {
                isay = 0;
                iw = 0;
                FILE* f17 = fortran_get_unit(17);
                for (i=i1; i<=i2; i+=i3) {
                    iw = iw + 1;
                    kw = 0;
                    for (k=k1; k<=k2; k+=k3) {
                        kw = kw + 1;
                        jw = 0;
                        if (k != 1 && k != kdim) continue;
                        for (j=j1; j<=j2; j+=j3) {
                            jw = jw + 1;
                            if (j == 1 || j == jdim) continue;
                            isay = isay+1;
                            if (isay == 1)
                                fprintf(f17, "\n   I   J   K         X                Y                Z"
                                        "               dn           P/Pinf           T/Tinf"
                                        "              Cf              Ch            yplus"
                                        "             Cfx             Cfy             Cfz\n");
                            if (ialph == 0) {
                                fprintf(f17, "%4d%4d%4d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                    i, j, k,
                                    (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,7),(float)xw(jw,kw,iw,8));
                            } else {
                                xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                                xw(jw,kw,iw,7) = -xw(jw,kw,iw,7);
                                fprintf(f17, "%4d%4d%4d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                    i, j, k,
                                    (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,8),(float)xw(jw,kw,iw,7));
                            }
                        }
                    }
                }
            } // end if (myid==myhost)
        } // end if (ivisc[2]!=0)


        if (ivisc[1] != 0) {
            // output viscous-flow data on j=1 and/or j=jdim surfaces
            iw = 0;
            for (i=i1; i<=i2; i+=i3) {
                iw = iw + 1;
                id   = i;
                id1  = id-1;
                if (id1 < 1)    id1 = 1;
                if (id > idim1) id  = idim1;
                jw = 0;
                for (j=j1; j<=j2; j+=j3) {
                    jw = jw + 1;
                    kw = 0;
                    if (j != 1 && j != jdim) continue;
                    for (k=k1; k<=k2; k+=k3) {
                        kw = kw + 1;
                        if (k == 1 || k == kdim) continue;
                        if (j == 1) {
                            jd  = 1;
                            jd1 = jd+1;
                            jdx = 1;
                            m   = 2;
                            sgn = 1.;
                        } else {
                            jd  = jdim1;
                            jd1 = jd-1;
                            jdx = jdim;
                            m   = 4;
                            sgn = -1.;
                        }

                        // wall value - average in K and I
                        q1j1 =       .25*(qj0(k,id,1,m)  +qj0(k-1,id,1,m)
                                        + qj0(k,id1,1,m) +qj0(k-1,id1,1,m));
                        q2j1 =       .25*(qj0(k,id,2,m)  +qj0(k-1,id,2,m)
                                        + qj0(k,id1,2,m) +qj0(k-1,id1,2,m));
                        q3j1 =       .25*(qj0(k,id,3,m)  +qj0(k-1,id,3,m)
                                        + qj0(k,id1,3,m) +qj0(k-1,id1,3,m));
                        q4j1 =       .25*(qj0(k,id,4,m)  +qj0(k-1,id,4,m)
                                        + qj0(k,id1,4,m) +qj0(k-1,id1,4,m));
                        q5j1 = gamma*.25*(qj0(k,id,5,m)  +qj0(k-1,id,5,m)
                                        + qj0(k,id1,5,m) +qj0(k-1,id1,5,m));
                        t1j1 = q5j1/q1j1;

                        // first cell center location - average in K and I
                        q1j2 =       .25*(q(jd,k,id,1)   +q(jd,k-1,id,1)
                                        + q(jd,k,id1,1)  +q(jd,k-1,id1,1));
                        q2j2 =       .25*(q(jd,k,id,2)   +q(jd,k-1,id,2)
                                        + q(jd,k,id1,2)  +q(jd,k-1,id1,2));
                        q3j2 =       .25*(q(jd,k,id,3)   +q(jd,k-1,id,3)
                                        + q(jd,k,id1,3)  +q(jd,k-1,id1,3));
                        q4j2 =       .25*(q(jd,k,id,4)   +q(jd,k-1,id,4)
                                        + q(jd,k,id1,4)  +q(jd,k-1,id1,4));
                        q5j2 = gamma*.25*(q(jd,k,id,5)   +q(jd,k-1,id,5)
                                        + q(jd,k,id1,5)  +q(jd,k-1,id1,5));
                        t1j2 = q5j2/q1j2;

                        if (sj(jdx,k  ,id ,4) == 0. ||
                            sj(jdx,k-1,id ,4) == 0. ||
                            sj(jdx,k  ,id1,4) == 0. ||
                            sj(jdx,k-1,id1,4) == 0.) {
                            dx = x(jd+1,k,i)-x(jd,k,i);
                            dy = y(jd+1,k,i)-y(jd,k,i);
                            dz = z(jd+1,k,i)-z(jd,k,i);
                            dn = std::sqrt(dx*dx+dy*dy+dz*dz);
                        } else {
                            dn = (vol(jd,k  ,id )/sj(jdx,k  ,id ,4)+
                                  vol(jd,k-1,id )/sj(jdx,k-1,id ,4)+
                                  vol(jd,k  ,id1)/sj(jdx,k  ,id1,4)+
                                  vol(jd,k-1,id1)/sj(jdx,k-1,id1,4))/4.;
                        }

                        // Get turb viscosity at wall
                        if (iwf[1] == 1) {
                            avgmut = .25*(vj0(k,id,1,m)  +vj0(k-1,id,1,m)
                                        + vj0(k,id1,1,m) +vj0(k-1,id1,1,m));
                        } else {
                            avgmut = 0.;
                        }
                        emuka = std::pow(t1j1,1.5)*((1.0+cbar/tinf)/(t1j1+cbar/tinf));

                        // Use component of velocity parallel to wall
                        urel = q2j2-q2j1;
                        vrel = q3j2-q3j1;
                        wrel = q4j2-q4j1;
                        sj1  = (sj(jdx,k  ,id ,1)+sj(jdx,k-1,id ,1)+
                                sj(jdx,k  ,id1,1)+sj(jdx,k-1,id1,1))/4.;
                        sj2  = (sj(jdx,k  ,id ,2)+sj(jdx,k-1,id ,2)+
                                sj(jdx,k  ,id1,2)+sj(jdx,k-1,id1,2))/4.;
                        sj3  = (sj(jdx,k  ,id ,3)+sj(jdx,k-1,id ,3)+
                                sj(jdx,k  ,id1,3)+sj(jdx,k-1,id1,3))/4.;
                        vnorm = (urel*sj1 + vrel*sj2 + wrel*sj3)*sgn;
                        upar = urel-vnorm*sj1*sgn;
                        vpar = vrel-vnorm*sj2*sgn;
                        wpar = wrel-vnorm*sj3*sgn;
                        Cf   = 2.*std::sqrt(upar*upar + vpar*vpar + wpar*wpar);
                        Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn;
                        if ((float)q2j2 < 0.) Cf = -Cf;
                        // Cf vector
                        cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn;
                        cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn;
                        cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn;

                        eps  = 1.0e-6;
                        tty1 = 1.0 + gm1*0.5*xmach*xmach;
                        t1m1 = t1j1-tty1;
                        if (std::abs((float)t1m1) <= (float)eps) {
                            Ch = 999.99999;
                        } else {
                            Ch = 2.*(t1j2-t1j1);
                            Ch = (emuka+avgmut)/(reue*pr*(t1j1-tty1))*Ch/dn;
                        }

                        pres1 = q5j1;
                        temp1 = t1j1;
                        yplus = 0.;
                        if (std::abs((float)Cf) > 0.) {
                            double cfabs_arg = Cf*0.5/q1j1;
                            yplus = dn*reue*q1j1*std::sqrt(ccomplex_ns::ccabs(cfabs_arg))/emuka;
                        }

                        // store data in single precision arrays
                        xg(jw,kw,iw,1) = x(j,k,i);
                        xg(jw,kw,iw,2) = y(j,k,i);
                        xg(jw,kw,iw,3) = z(j,k,i);
                        xg(jw,kw,iw,4) = dn;

                        xw(jw,kw,iw,1) = pres1;
                        xw(jw,kw,iw,2) = temp1;
                        xw(jw,kw,iw,3) = Cf;
                        xw(jw,kw,iw,4) = Ch;
                        xw(jw,kw,iw,5) = yplus;
                        xw(jw,kw,iw,6) = cfx;
                        xw(jw,kw,iw,7) = cfy;
                        xw(jw,kw,iw,8) = cfz;
                    } // k loop
                } // j loop
            } // i loop

            // write printout data to file
            if (myid == myhost) {
                isay = 0;
                iw = 0;
                FILE* f17 = fortran_get_unit(17);
                for (i=i1; i<=i2; i+=i3) {
                    iw = iw + 1;
                    jw = 0;
                    for (j=j1; j<=j2; j+=j3) {
                        jw = jw + 1;
                        kw = 0;
                        if (j != 1 && j != jdim) continue;
                        for (k=k1; k<=k2; k+=k3) {
                            kw = kw + 1;
                            if (k == 1 || k == kdim) continue;
                            isay = isay+1;
                            if (isay == 1)
                                fprintf(f17, "\n   I   J   K         X                Y                Z"
                                        "               dn           P/Pinf           T/Tinf"
                                        "              Cf              Ch            yplus"
                                        "             Cfx             Cfy             Cfz\n");
                            if (ialph == 0) {
                                fprintf(f17, "%4d%4d%4d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                    i, j, k,
                                    (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,7),(float)xw(jw,kw,iw,8));
                            } else {
                                xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                                xw(jw,kw,iw,7) = -xw(jw,kw,iw,7);
                                fprintf(f17, "%4d%4d%4d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                    i, j, k,
                                    (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,8),(float)xw(jw,kw,iw,7));
                            }
                        }
                    }
                }
            } // end if (myid==myhost)
        } // end if (ivisc[1]!=0)


        if (idim == 2) return;
        if (ivisc[0] != 0) {
            // output viscous-flow data on i=1 and/or i=idim surfaces
            jw = 0;
            for (j=j1; j<=j2; j+=j3) {
                jw = jw + 1;
                iw = 0;
                if (j == 1 || j == jdim) continue;
                for (i=i1; i<=i2; i+=i3) {
                    iw = iw + 1;
                    kw = 0;
                    if (i != 1 && i != idim) continue;
                    for (k=k1; k<=k2; k+=k3) {
                        kw = kw + 1;
                        if (k == 1 || k == kdim) continue;
                        if (i == 1) {
                            id  = 1;
                            id1 = id+1;
                            idx = 1;
                            m   = 2;
                            sgn = 1.;
                        } else {
                            id  = idim1;
                            id1 = id-1;
                            idx = idim;
                            m   = 4;
                            sgn = -1.;
                        }

                        // wall value - average in J and K
                        q1i1 =       .25*(qi0(j,k,1,m)   +qi0(j-1,k,1,m)
                                        + qi0(j,k-1,1,m) +qi0(j-1,k-1,1,m));
                        q2i1 =       .25*(qi0(j,k,2,m)   +qi0(j-1,k,2,m)
                                        + qi0(j,k-1,2,m) +qi0(j-1,k-1,2,m));
                        q3i1 =       .25*(qi0(j,k,3,m)   +qi0(j-1,k,3,m)
                                        + qi0(j,k-1,3,m) +qi0(j-1,k-1,3,m));
                        q4i1 =       .25*(qi0(j,k,4,m)   +qi0(j-1,k,4,m)
                                        + qi0(j,k-1,4,m) +qi0(j-1,k-1,4,m));
                        q5i1 = gamma*.25*(qi0(j,k,5,m)   +qi0(j-1,k,5,m)
                                        + qi0(j,k-1,5,m) +qi0(j-1,k-1,5,m));
                        t1i1 = q5i1/q1i1;

                        // first cell center location - average in J and K
                        q1i2 =       .25*(q(j,k,id,1)    +q(j,k-1,id,1)
                                        + q(j-1,k,id,1)  +q(j-1,k-1,id,1));
                        q2i2 =       .25*(q(j,k,id,2)    +q(j,k-1,id,2)
                                        + q(j-1,k,id,2)  +q(j-1,k-1,id,2));
                        q3i2 =       .25*(q(j,k,id,3)    +q(j,k-1,id,3)
                                        + q(j-1,k,id,3)  +q(j-1,k-1,id,3));
                        q4i2 =       .25*(q(j,k,id,4)    +q(j,k-1,id,4)
                                        + q(j-1,k,id,4)  +q(j-1,k-1,id,4));
                        q5i2 = gamma*.25*(q(j,k,id,5)    +q(j,k-1,id,5)
                                        + q(j-1,k,id,5)  +q(j-1,k-1,id,5));
                        t1i2 = q5i2/q1i2;

                        if (si(j  ,k  ,idx,4) == 0. ||
                            si(j  ,k-1,idx,4) == 0. ||
                            si(j-1,k  ,idx,4) == 0. ||
                            si(j-1,k-1,idx,4) == 0.) {
                            dx = x(j,k,id+1)-x(j,k,id);
                            dy = y(j,k,id+1)-y(j,k,id);
                            dz = z(j,k,id+1)-z(j,k,id);
                            dn = std::sqrt(dx*dx+dy*dy+dz*dz);
                        } else {
                            dn = (vol(j  ,k  ,id)/si(j  ,k  ,idx,4)+
                                  vol(j  ,k-1,id)/si(j  ,k-1,idx,4)+
                                  vol(j-1,k  ,id)/si(j-1,k  ,idx,4)+
                                  vol(j-1,k-1,id)/si(j-1,k-1,idx,4))/4.;
                        }

                        // Get turb viscosity at wall
                        if (iwf[0] == 1) {
                            avgmut = .25*(vi0(j,k,1,m)   +vi0(j-1,k,1,m)
                                        + vi0(j,k-1,1,m) +vi0(j-1,k-1,1,m));
                        } else {
                            avgmut = 0.;
                        }
                        emuka = std::pow(t1i1,1.5)*((1.0+cbar/tinf)/(t1i1+cbar/tinf));

                        // Use component of velocity parallel to wall
                        urel = q2i2-q2i1;
                        vrel = q3i2-q3i1;
                        wrel = q4i2-q4i1;
                        si1  = (si(j  ,k  ,idx,1)+si(j  ,k-1,idx,1)+
                                si(j-1,k  ,idx,1)+si(j-1,k-1,idx,1))/4.;
                        si2  = (si(j  ,k  ,idx,2)+si(j  ,k-1,idx,2)+
                                si(j-1,k  ,idx,2)+si(j-1,k-1,idx,2))/4.;
                        si3  = (si(j  ,k  ,idx,3)+si(j  ,k-1,idx,3)+
                                si(j-1,k  ,idx,3)+si(j-1,k-1,idx,3))/4.;
                        vnorm = (urel*si1 + vrel*si2 + wrel*si3)*sgn;
                        upar = urel-vnorm*si1*sgn;
                        vpar = vrel-vnorm*si2*sgn;
                        wpar = wrel-vnorm*si3*sgn;
                        Cf   = 2.*std::sqrt(upar*upar + vpar*vpar + wpar*wpar);
                        Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn;
                        if ((float)q2i2 < 0.) Cf = -Cf;
                        // Cf vector
                        cfx  = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn;
                        cfy  = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn;
                        cfz  = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn;

                        eps  = 1.0e-6;
                        tty1 = 1.0 + gm1*0.5*xmach*xmach;
                        t1m1 = t1i1-tty1;
                        if (std::abs((float)t1m1) <= (float)eps) {
                            Ch = 999.99999;
                        } else {
                            Ch = 2.*(t1i2-t1i1);
                            Ch = (emuka+avgmut)/(reue*pr*(t1i1-tty1))*Ch/dn;
                        }

                        pres1 = q5i1;
                        temp1 = t1i1;
                        yplus = 0.;
                        if (std::abs((float)Cf) > 0.) {
                            double cfabs_arg = Cf*0.5/q1i1;
                            yplus = dn*reue*q1i1*std::sqrt(ccomplex_ns::ccabs(cfabs_arg))/emuka;
                        }

                        // store data in single precision arrays
                        xg(jw,kw,iw,1) = x(j,k,i);
                        xg(jw,kw,iw,2) = y(j,k,i);
                        xg(jw,kw,iw,3) = z(j,k,i);
                        xg(jw,kw,iw,4) = dn;

                        xw(jw,kw,iw,1) = pres1;
                        xw(jw,kw,iw,2) = temp1;
                        xw(jw,kw,iw,3) = Cf;
                        xw(jw,kw,iw,4) = Ch;
                        xw(jw,kw,iw,5) = yplus;
                        xw(jw,kw,iw,6) = cfx;
                        xw(jw,kw,iw,7) = cfy;
                        xw(jw,kw,iw,8) = cfz;
                    } // k loop
                } // i loop
            } // j loop

            // write printout data to file
            if (myid == myhost) {
                isay = 0;
                jw = 0;
                FILE* f17 = fortran_get_unit(17);
                for (j=j1; j<=j2; j+=j3) {
                    jw = jw + 1;
                    iw = 0;
                    if (j == 1 || j == jdim) continue;
                    for (i=i1; i<=i2; i+=i3) {
                        iw = iw + 1;
                        kw = 0;
                        if (i != 1 && i != idim) continue;
                        for (k=k1; k<=k2; k+=k3) {
                            kw = kw + 1;
                            if (k == 1 || k == kdim) continue;
                            isay = isay+1;
                            if (isay == 1)
                                fprintf(f17, "\n   I   J   K         X                Y                Z"
                                        "               dn           P/Pinf           T/Tinf"
                                        "              Cf              Ch            yplus"
                                        "             Cfx             Cfy             Cfz\n");
                            if (ialph == 0) {
                                fprintf(f17, "%4d%4d%4d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                    i, j, k,
                                    (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,7),(float)xw(jw,kw,iw,8));
                            } else {
                                xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                                xw(jw,kw,iw,7) = -xw(jw,kw,iw,7);
                                fprintf(f17, "%4d%4d%4d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                    i, j, k,
                                    (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                    (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,8),(float)xw(jw,kw,iw,7));
                            }
                        }
                    }
                }
            } // end if (myid==myhost)
        } // end if (ivisc[0]!=0)

    } // end if (iflag==2)

    return;
}

} // namespace plot3d_coarse_ns
