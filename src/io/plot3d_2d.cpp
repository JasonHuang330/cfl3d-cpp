// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "plot3d_2d.h"
#include "plot3d.h"
#include "cctogp.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace plot3d_2d_ns {

// -----------------------------------------------------------------------
// plot3d: thin wrapper that delegates to plot3d_ns::plot3d
// -----------------------------------------------------------------------
void plot3d(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3,
            int& j1, int& j2, int& j3, int& k1, int& k2, int& k3,
            FortranArray4DRef<double> q, FortranArray4DRef<double> qi0,
            FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
            FortranArray3DRef<double> x, FortranArray3DRef<double> y,
            FortranArray3DRef<double> z, FortranArray4DRef<double> xw,
            FortranArray4DRef<double> blank2, FortranArray3DRef<double> blank,
            FortranArray4DRef<double> xg, int& iflag,
            FortranArray3DRef<double> vist3d, int& iover, int& nblk,
            FortranArray1DRef<int> nmap, FortranArray3DRef<double> bcj,
            FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
            FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
            FortranArray4DRef<double> vi0, int& ifunc, int& iplot,
            int& jdw, int& kdw, int& idw, int& nplots,
            FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
            FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg,
            FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg,
            FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg,
            FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
            int& ninter, FortranArray2DRef<int> iindex, int& intmax,
            int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk,
            int& nbli, FortranArray3DRef<int> limblk,
            FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon,
            int& mxbli, FortranArray1DRef<double> thetay,
            int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm,
            FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d,
            int& nblock, FortranArray1DRef<int> nblkpt,
            FortranArray4DRef<double> xv, FortranArray4DRef<double> sj,
            FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
            FortranArray3DRef<double> vol, int& nset)
{
    plot3d_ns::plot3d(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                      q, qi0, qj0, qk0, x, y, z, xw, blank2, blank, xg,
                      iflag, vist3d, iover, nblk, nmap, bcj, bck, bci,
                      vj0, vk0, vi0, ifunc, iplot, jdw, kdw, idw, nplots,
                      jdimg, kdimg, idimg, nblcg, jsg, ksg, isg, jeg, keg, ieg,
                      ninter, iindex, intmax, nsub1, maxxe, nblkk, nbli,
                      limblk, isva, nblon, mxbli, thetay, maxbl, maxgr,
                      myid, myhost, mycomm, mblk2nd, inpl3d, nblock, nblkpt,
                      xv, sj, sk, si, vol, nset);
}

// -----------------------------------------------------------------------
// plot3d_2d: write output file at grid points in PLOT3D format
// -----------------------------------------------------------------------
void plot3d_2d(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3,
               int& j1, int& j2, int& j3, int& k1, int& k2, int& k3,
               FortranArray4DRef<double> q, FortranArray4DRef<double> qi0,
               FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0,
               FortranArray3DRef<double> x, FortranArray3DRef<double> y,
               FortranArray3DRef<double> z, FortranArray4DRef<double> xw,
               FortranArray4DRef<double> blank2, FortranArray3DRef<double> blank,
               FortranArray4DRef<double> xg, int& iflag,
               FortranArray3DRef<double> vist3d, int& iover, int& nblk,
               FortranArray1DRef<int> nmap, FortranArray3DRef<double> bcj,
               FortranArray3DRef<double> bck, FortranArray3DRef<double> bci,
               FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0,
               FortranArray4DRef<double> vi0, int& ifunc, int& iplot,
               int& jdw, int& kdw, int& idw, int& nplots,
               FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
               FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg,
               FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg,
               FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg,
               FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg,
               int& ninter, FortranArray2DRef<int> iindex, int& intmax,
               int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk,
               int& nbli, FortranArray3DRef<int> limblk,
               FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon,
               int& mxbli, FortranArray1DRef<double> thetay,
               int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm,
               FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d,
               int& nblock, FortranArray1DRef<int> nblkpt,
               FortranArray4DRef<double> xv, FortranArray4DRef<double> sj,
               FortranArray4DRef<double> sk, FortranArray4DRef<double> si,
               FortranArray3DRef<double> vol, int& nset)
{
    // COMMON block aliases
    int&   ibin       = cmn_bin.ibin;
    int&   iblnk      = cmn_bin.iblnk;
    double gamma      = (double)cmn_fluid.gamma;
    double gm1        = (double)cmn_fluid.gm1;
    double gp1        = (double)cmn_fluid.gp1;
    double pr         = (double)cmn_fluid2.pr;
    double cbar       = (double)cmn_fluid2.cbar;
    double xmach      = (double)cmn_info.xmach;
    double alpha      = (double)cmn_info.alpha;
    double reue       = (double)cmn_reyue.reue;
    double tinf       = (double)cmn_reyue.tinf;
    int&   i2d        = cmn_twod.i2d;
    double time_v     = (double)cmn_unst.time;
    int&   ialph      = cmn_igrdtyp.ialph;
    int&   lhdr       = cmn_moov.lhdr;
    int&   icall2d    = cmn_moovcrs2d.icall2d;
    double radtodeg   = (double)cmn_conversion.radtodeg;
    double p0         = (double)cmn_ivals.p0;

    // Local variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;

    int jw = 0, kw = 0, iw = 0;
    int j = 0, k = 0, i = 0, l = 0, m = 0;
    double alphaw = 0., xmachw = 0., alphww = 0., reuew = 0., timew = 0.;
    double temp = 0.;
    int izero = 0;
    int j11 = 0, j22_loc = 0, k11 = 0, k22_loc = 0, i11 = 0, i22_loc = 0;
    int jj = 0, kk = 0, ii = 0;
    int nnn = 0;
    int inter = 0, lmax1 = 0, nbl = 0, lst = 0, lcoord = 0, lend = 0;
    int j21 = 0, k21 = 0;
    int ll = 0, mblk1 = 0, mblk2 = 0, mblk3 = 0, mblk4 = 0;
    int n = 0, it = 0, ir = 0, itime = 0, iti = 0;
    int is = 0, ie = 0, js = 0, je = 0, ks = 0, ke = 0;
    int is1 = 0, ie1 = 0, js1 = 0, je1 = 0, ks1 = 0, ke1 = 0;
    int nblc = 0, nblcc = 0;
    int ibembed = 0;
    double cpc = 0.;
    int leq = 0;
    // printout section
    int icp = 0;
    double term1 = 0., term2 = 0., term3 = 0., pitot = 0.;
    double q1 = 0., q2 = 0., q3 = 0., q4 = 0., q5 = 0., t1 = 0., xm1 = 0., edvis = 0.;
    int isay = 0;
    // viscous surface section
    int id = 0, id1 = 0, kd = 0, kd1 = 0, kdx = 0, jd = 0, jd1 = 0, jdx = 0, idx = 0;
    double q1k1 = 0., q2k1 = 0., q3k1 = 0., q4k1 = 0., q5k1 = 0., t1k1 = 0.;
    double q1k2 = 0., q2k2 = 0., q3k2 = 0., q4k2 = 0., q5k2 = 0., t1k2 = 0.;
    double q1j1 = 0., q2j1 = 0., q3j1 = 0., q4j1 = 0., q5j1 = 0., t1j1 = 0.;
    double q1j2 = 0., q2j2 = 0., q3j2 = 0., q4j2 = 0., q5j2 = 0., t1j2 = 0.;
    double q1i1 = 0., q2i1 = 0., q3i1 = 0., q4i1 = 0., q5i1 = 0., t1i1 = 0.;
    double q1i2 = 0., q2i2 = 0., q3i2 = 0., q4i2 = 0., q5i2 = 0., t1i2 = 0.;
    double dx = 0., dy = 0., dz = 0., dn = 0.;
    double avgmut = 0., emuka = 0.;
    double urel = 0., vrel = 0., wrel = 0.;
    double sk1 = 0., sk2 = 0., sk3 = 0., sj1 = 0., sj2 = 0., sj3 = 0.;
    double si1 = 0., si2 = 0., si3 = 0.;
    double vnorm = 0., upar = 0., vpar = 0., wpar = 0.;
    double Cf = 0., cfx = 0., cfy = 0., cfz = 0., Ch = 0.;
    double eps = 0., tty1 = 0., t1m1 = 0.;
    double pres1 = 0., temp1 = 0., yplus = 0.;
    double sgn = 0.;
    int iold = 0, jv = 0, kv = 0, iv = 0, jx = 0, kx = 0, ix = 0;
    int ldw = 0;

    // ---- body appended below ----


    // initialize xw, xv and xg arrays
    jw = (j2-j1)/j3 + 1;
    kw = (k2-k1)/k3 + 1;
    iw = (i2-i1)/i3 + 1;
    for (j = 1; j <= jw; j++) {
        for (k = 1; k <= kw; k++) {
            for (i = 1; i <= iw; i++) {
                for (l = 1; l <= 5; l++) {
                    xw(j,k,i,l) = 0.;
                    xv(j,k,i,l) = 0.;
                }
                for (l = 1; l <= 4; l++) {
                    xg(j,k,i,l) = 0.;
                }
            }
        }
    }

    // assign single precision scalars
    alphaw = radtodeg * (alpha + (double)thetay(nblk));
    xmachw = xmach;
    alphww = alphaw;
    reuew  = reue;
    timew  = time_v;

    // load grid into first 3 locations of the single precision xg array
    iw = 0;
    for (i = i1; i <= i2; i += i3) {
        iw = iw + 1;
        jw = 0;
        for (j = j1; j <= j2; j += j3) {
            jw = jw + 1;
            kw = 0;
            for (k = k1; k <= k2; k += k3) {
                kw = kw + 1;
                xg(jw,kw,iw,1) = x(j,k,i);
                xg(jw,kw,iw,2) = y(j,k,i);
                xg(jw,kw,iw,3) = z(j,k,i);
            }
        }
    }

    // set iblank (blank2) array
    // currently don't need blanking data for printout option (iflag=2),
    // so only compute blank2 array for plot3d output option
    // also, don't compute blanking info if iblnk = 0
    if (iflag == 1 && iblnk > 0) {

        // assign default iblank (blank2) array
        for (i = 1; i <= idim; i++) {
            for (j = 1; j <= jdim; j++) {
                for (k = 1; k <= kdim; k++) {
                    blank2(j,k,i,1) = 1.;
                    blank2(j,k,i,2) = 1.;
                }
            }
        }

        // zero out edges and corners if desired (izero > 0)
        izero = 0;
        if (izero > 0) {
            for (i = 1; i <= idim; i++) {
                for (j = 1; j <= jdim; j += jdim1) {
                    for (k = 1; k <= kdim; k += kdim1) {
                        blank2(j,k,i,1) = 0.;
                    }
                }
            }
            if (idim > 2) {
                for (i = 1; i <= idim; i += idim1) {
                    for (j = 1; j <= jdim; j++) {
                        blank2(j,1,i,1)    = 0.;
                        blank2(j,kdim,i,1) = 0.;
                    }
                    for (k = 1; k <= kdim; k++) {
                        blank2(1,k,i,1)    = 0.;
                        blank2(jdim,k,i,1) = 0.;
                    }
                }
            } else {
                for (j = 1; j <= jdim; j++) {
                    for (k = 1; k <= kdim; k++) {
                        blank2(j,k,idim,1) = 0.;
                    }
                }
            }
            if (jdim == 2) {
                for (i = 1; i <= idim; i++) {
                    for (k = 1; k <= kdim; k++) {
                        blank2(jdim,k,i,1) = 0.;
                    }
                }
            }
            if (kdim == 2) {
                for (i = 1; i <= idim; i++) {
                    for (k = 1; k <= kdim; k++) {
                        blank2(j,kdim,i,1) = 0.;
                    }
                }
            }
        } // end izero > 0

        // solid surface iblank (blank2) values (iblank=2 for solid surface)
        j11 = 1;
        j22_loc = jdim1;
        if (jdim == 2) { j11 = 1; j22_loc = 1; }
        k11 = 1;
        k22_loc = kdim1;
        if (kdim == 2) { k11 = 1; k22_loc = 1; }
        i11 = 1;
        i22_loc = idim1;
        if (idim == 2) { i11 = 1; i22_loc = 1; }

        i = 1;
        for (nnn = 1; nnn <= 2; nnn++) {
            for (j = j11; j <= j22_loc; j++) {
                for (k = k11; k <= k22_loc; k++) {
                    kk = std::min(k+1, k22_loc);
                    jj = std::min(j+1, j22_loc);
                    blank2(j,k,i,1) = 1. + std::max({(double)bci(j,k,nnn),
                                                      (double)bci(jj,k,nnn),
                                                      (double)bci(j,kk,nnn),
                                                      (double)bci(jj,kk,nnn)});
                }
            }
            i = idim;
        }
        j = 1;
        for (nnn = 1; nnn <= 2; nnn++) {
            for (i = i11; i <= i22_loc; i++) {
                for (k = k11; k <= k22_loc; k++) {
                    kk = std::min(k+1, k22_loc);
                    ii = std::min(i+1, i22_loc);
                    blank2(j,k,i,1) = 1. + std::max({(double)bcj(k,i,nnn),
                                                      (double)bcj(kk,i,nnn),
                                                      (double)bcj(k,ii,nnn),
                                                      (double)bcj(kk,ii,nnn)});
                }
            }
            j = jdim;
        }
        k = 1;
        for (nnn = 1; nnn <= 2; nnn++) {
            for (j = j11; j <= j22_loc; j++) {
                for (i = i11; i <= i22_loc; i++) {
                    jj = std::min(j+1, j22_loc);
                    ii = std::min(i+1, i22_loc);
                    blank2(j,k,i,1) = 1. + std::max({(double)bck(j,i,nnn),
                                                      (double)bck(jj,i,nnn),
                                                      (double)bck(j,ii,nnn),
                                                      (double)bck(jj,ii,nnn)});
                }
            }
            k = kdim;
        }



        // iblank (blank2) array for patch interface boundaries
        if (std::abs(ninter) > 0) {
            for (inter = 1; inter <= std::abs(ninter); inter++) {
                lmax1  = iindex(inter,1);
                nbl    = iindex(inter,lmax1+2);
                if (nbl != nblk) continue; // go to 1600
                lst    = iindex(inter,2*lmax1+5);
                lcoord = iindex(inter,2*lmax1+3)/10;
                lend   = iindex(inter,2*lmax1+3)-lcoord*10;
                j21    = iindex(inter,2*lmax1+6);
                j22_loc= iindex(inter,2*lmax1+7);
                k21    = iindex(inter,2*lmax1+8);
                k22_loc= iindex(inter,2*lmax1+9);

                if (lcoord == 1) {
                    if (lend == 1) i = 1;
                    if (lend == 2) i = idim;
                    if (jdim > 2 && kdim > 2) {
                        for (j = j21+1; j <= j22_loc-1; j++) {
                            for (k = k21+1; k <= k22_loc-1; k++) {
                                ll = lst + (j22_loc-j21)*(k-1-k21) + (j-1-j21);
                                mblk1 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(k-k21-1) + (j-j21);
                                mblk2 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(k-k21) + (j-j21);
                                mblk3 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(k-k21) + (j-j21-1);
                                mblk4 = iindex(inter,nblkpt(ll)+1);
                                if (mblk1==mblk2 && mblk1==mblk3 && mblk1==mblk4)
                                    blank2(j,k,i,1) = -(double)nmap(mblk1);
                            }
                        }
                    } else if (jdim > 2) {
                        k = k21;
                        for (j = j21+1; j <= j22_loc-1; j++) {
                            ll = lst + (j22_loc-j21)*(k-k21) + (j-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22_loc-j21)*(k-k21) + (j-j21-1);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1 == mblk2) blank2(j,k,i,1) = -(double)nmap(mblk1);
                        }
                    } else if (kdim > 2) {
                        j = j21;
                        for (k = k21+1; k <= k22_loc-1; k++) {
                            ll = lst + (j22_loc-j21)*(k-k21-1) + (j-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22_loc-j21)*(k-k21) + (j-j21);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1 == mblk2) blank2(j,k,i,1) = -(double)nmap(mblk1);
                        }
                    }
                }

                if (lcoord == 2) {
                    if (lend == 1) j = 1;
                    if (lend == 2) j = jdim;
                    if (idim > 2 && kdim > 2) {
                        for (i = k21+1; i <= k22_loc-1; i++) {
                            for (k = j21+1; k <= j22_loc-1; k++) {
                                ll = lst + (j22_loc-j21)*(i-1-k21) + (k-1-j21);
                                mblk1 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(i-k21-1) + (k-j21);
                                mblk2 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(i-k21) + (k-j21);
                                mblk3 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(i-k21) + (k-j21-1);
                                mblk4 = iindex(inter,nblkpt(ll)+1);
                                if (mblk1==mblk2 && mblk1==mblk3 && mblk1==mblk4)
                                    blank2(j,k,i,1) = -(double)nmap(mblk1);
                            }
                        }
                    } else if (idim > 2) {
                        k = j21;
                        for (i = k21+1; i <= k22_loc-1; i++) {
                            ll = lst + (j22_loc-j21)*(i-k21-1) + (k-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22_loc-j21)*(i-k21) + (k-j21);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1 == mblk2) blank2(j,k,i,1) = -(double)nmap(mblk1);
                        }
                    } else if (kdim > 2) {
                        i = k21;
                        for (k = j21+1; k <= j22_loc-1; k++) {
                            ll = lst + (j22_loc-j21)*(i-k21) + (k-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22_loc-j21)*(i-k21) + (k-j21-1);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1 == mblk2) blank2(j,k,i,1) = -(double)nmap(mblk1);
                        }
                    }
                }

                if (lcoord == 3) {
                    if (lend == 1) k = 1;
                    if (lend == 2) k = kdim;
                    if (idim > 2 && jdim > 2) {
                        for (i = k21+1; i <= k22_loc-1; i++) {
                            for (j = j21+1; j <= j22_loc-1; j++) {
                                ll = lst + (j22_loc-j21)*(i-1-k21) + (j-1-j21);
                                mblk1 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(i-k21-1) + (j-j21);
                                mblk2 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(i-k21) + (j-j21);
                                mblk3 = iindex(inter,nblkpt(ll)+1);
                                ll = lst + (j22_loc-j21)*(i-k21) + (j-j21-1);
                                mblk4 = iindex(inter,nblkpt(ll)+1);
                                if (mblk1==mblk2 && mblk1==mblk3 && mblk1==mblk4)
                                    blank2(j,k,i,1) = -(double)nmap(mblk1);
                            }
                        }
                    } else if (idim > 2) {
                        j = j21;
                        for (i = k21+1; i <= k22_loc-1; i++) {
                            ll = lst + (j22_loc-j21)*(i-k21-1) + (j-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22_loc-j21)*(i-k21) + (j-j21);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1 == mblk2) blank2(j,k,i,1) = -(double)nmap(mblk1);
                        }
                    } else if (jdim > 2) {
                        i = k21;
                        for (j = j21+1; j <= j22_loc-1; j++) {
                            ll = lst + (j22_loc-j21)*(i-k21) + (j-j21);
                            mblk1 = iindex(inter,nblkpt(ll)+1);
                            ll = lst + (j22_loc-j21)*(i-k21) + (j-j21-1);
                            mblk2 = iindex(inter,nblkpt(ll)+1);
                            if (mblk1 == mblk2) blank2(j,k,i,1) = -(double)nmap(mblk1);
                        }
                    }
                }
            } // end inter loop (1600)
        } // end abs(ninter) > 0



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
                        // allow for 1-1 blocking in same grid
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
                            // cell center indices ---> grid point indices
                            if (isva(it,1,n)+isva(it,2,n) == 5) {
                                // i = constant interface
                                if (js > je) js = js+1;
                                if (js < je) je = je+1;
                                if (ks > ke) ks = ks+1;
                                if (ks < ke) ke = ke+1;
                                // 2d cases
                                if (jdim == 2) { js = 1; je = 2; }
                                if (kdim == 2) { ks = 1; ke = 2; }
                            }
                            if (isva(it,1,n)+isva(it,2,n) == 4) {
                                // j = constant interface
                                if (is > ie) is = is+1;
                                if (is < ie) ie = ie+1;
                                if (ks > ke) ks = ks+1;
                                if (ks < ke) ke = ke+1;
                                // 2d cases
                                if (idim == 2) { is = 1; ie = 2; }
                                if (kdim == 2) { ks = 1; ke = 2; }
                            }
                            if (isva(it,1,n)+isva(it,2,n) == 3) {
                                // k = constant interface
                                if (js > je) js = js+1;
                                if (js < je) je = je+1;
                                if (is > ie) is = is+1;
                                if (is < ie) ie = ie+1;
                                // 2d cases
                                if (jdim == 2) { js = 1; je = 2; }
                                if (idim == 2) { is = 1; ie = 2; }
                            }
                            is1 = std::min(is,ie);
                            ie1 = std::max(is,ie);
                            js1 = std::min(js,je);
                            je1 = std::max(js,je);
                            ks1 = std::min(ks,ke);
                            ke1 = std::max(ks,ke);
                            for (i = is1; i <= ie1; i++) {
                                for (j = js1; j <= je1; j++) {
                                    for (k = ks1; k <= ke1; k++) {
                                        blank2(j,k,i,1) = -(double)nmap(nblkk(ir,n));
                                    }
                                }
                            }
                        } // iti loop (101)
                    }
                }
            } // n loop (100)
        } // nbli > 0

        // iblank (blank2) array for embedded grids
        ibembed = 1;
        if (ibembed > 0) {
            for (nblc = 1; nblc <= nblock; nblc++) {
                if (nblk == nblc) continue; // go to 7500
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
                    for (i = is; i <= ie; i++) {
                        for (j = js; j <= je; j++) {
                            for (k = ks; k <= ke; k++) {
                                blank2(j,k,i,1) = 0.;
                            }
                        }
                    }
                }
            } // nblc loop (7500)
        } // ibembed > 0



        // iblank (blank2) array for overlapped grids
        if (iover == 1) {
            // interior of faces (i=1 and i=idim)
            for (j = 2; j <= jdim1; j++) {
                for (k = 2; k <= kdim1; k++) {
                    blank2(j,k,1,2)    = std::min({(double)blank(j,k,1),
                                                   (double)blank(j-1,k,1),
                                                   (double)blank(j,k-1,1),
                                                   (double)blank(j-1,k-1,1)});
                    blank2(j,k,idim,2) = std::min({(double)blank(j,k,idim1),
                                                   (double)blank(j-1,k,idim1),
                                                   (double)blank(j,k-1,idim1),
                                                   (double)blank(j-1,k-1,idim1)});
                }
            }
            // edges and corners (i=1 and i=idim)
            for (m = 1; m <= 2; m++) {
                k  = 1;
                kk = 1;
                if (m == 2) { k = kdim; kk = kdim1; }
                for (j = 2; j <= jdim1; j++) {
                    blank2(j,k,1,2)    = std::min((double)blank(j,kk,1),
                                                   (double)blank(j-1,kk,1));
                    blank2(j,k,idim,2) = std::min((double)blank(j,kk,idim1),
                                                   (double)blank(j-1,kk,idim1));
                }
            }
            for (m = 1; m <= 2; m++) {
                j  = 1;
                jj = 1;
                if (m == 2) { j = jdim; jj = jdim1; }
                for (k = 2; k <= kdim1; k++) {
                    blank2(j,k,1,2)    = std::min((double)blank(jj,k,1),
                                                   (double)blank(jj,k-1,1));
                    blank2(j,k,idim,2) = std::min((double)blank(jj,k,idim1),
                                                   (double)blank(jj,k-1,idim1));
                }
            }
            blank2(jdim,kdim,1,2)    = blank2(jdim-1,kdim-1,1,2);
            blank2(1,1,1,2)          = blank2(2,2,1,2);
            blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim-1,idim,2);
            blank2(1,1,idim,2)       = blank2(2,2,idim,2);

            // interior of faces (j=1 and j=jdim)
            for (i = 2; i <= idim1; i++) {
                for (k = 2; k <= kdim1; k++) {
                    blank2(1,k,i,2)    = std::min({(double)blank(1,k,i),
                                                   (double)blank(1,k,i-1),
                                                   (double)blank(1,k-1,i),
                                                   (double)blank(1,k-1,i-1)});
                    blank2(jdim,k,i,2) = std::min({(double)blank(jdim1,k,i),
                                                   (double)blank(jdim1,k,i-1),
                                                   (double)blank(jdim1,k-1,i),
                                                   (double)blank(jdim1,k-1,i-1)});
                }
            }
            // edges and corners (j=1 and j=jdim)
            for (m = 1; m <= 2; m++) {
                k  = 1;
                kk = 1;
                if (m == 2) { k = kdim; kk = kdim1; }
                for (i = 2; i <= idim1; i++) {
                    blank2(1,k,i,2)    = std::min((double)blank(1,kk,i),
                                                   (double)blank(1,kk,i-1));
                    blank2(jdim,k,i,2) = std::min((double)blank(jdim1,kk,i),
                                                   (double)blank(jdim1,kk,i-1));
                }
            }
            for (m = 1; m <= 2; m++) {
                i  = 1;
                ii = 1;
                if (m == 2) { i = idim; ii = idim1; }
                for (k = 2; k <= kdim1; k++) {
                    blank2(1,k,i,2)    = std::min((double)blank(1,k,ii),
                                                   (double)blank(1,k-1,ii));
                    blank2(jdim,k,i,2) = std::min((double)blank(jdim1,k,ii),
                                                   (double)blank(jdim1,k-1,ii));
                }
            }
            blank2(1,kdim,idim,2)    = blank2(1,kdim-1,idim-1,2);
            blank2(1,1,1,2)          = blank2(1,2,2,2);
            blank2(jdim,kdim,idim,2) = blank2(jdim,kdim-1,idim-1,2);
            blank2(jdim,1,1,2)       = blank2(jdim,2,2,2);

            // interior of faces (k=1 and k=kdim)
            for (i = 2; i <= idim1; i++) {
                for (j = 2; j <= jdim1; j++) {
                    blank2(j,1,i,2)    = std::min({(double)blank(j,1,i),
                                                   (double)blank(j,1,i-1),
                                                   (double)blank(j-1,1,i),
                                                   (double)blank(j-1,1,i-1)});
                    blank2(j,kdim,i,2) = std::min({(double)blank(j,kdim1,i),
                                                   (double)blank(j,kdim1,i-1),
                                                   (double)blank(j-1,kdim1,i),
                                                   (double)blank(j-1,kdim1,i-1)});
                }
            }
            // edges and corners (k=1 and k=kdim)
            for (m = 1; m <= 2; m++) {
                j  = 1;
                jj = 1;
                if (m == 2) { j = jdim; jj = jdim1; }
                for (i = 2; i <= idim1; i++) {
                    blank2(j,1,i,2)    = std::min((double)blank(jj,1,i),
                                                   (double)blank(jj,1,i-1));
                    blank2(j,kdim,i,2) = std::min((double)blank(jj,kdim1,i),
                                                   (double)blank(jj,kdim1,i-1));
                }
            }
            for (m = 1; m <= 2; m++) {
                i  = 1;
                ii = 1;
                if (m == 2) { i = idim; ii = idim1; }
                for (j = 2; j <= jdim1; j++) {
                    blank2(j,1,i,2)    = std::min((double)blank(j,1,ii),
                                                   (double)blank(j-1,1,ii));
                    blank2(j,kdim,i,2) = std::min((double)blank(j,kdim1,ii),
                                                   (double)blank(j-1,kdim1,ii));
                }
            }
            blank2(jdim,1,idim,2)    = blank2(jdim-1,1,idim-1,2);
            blank2(1,1,1,2)          = blank2(2,1,2,2);
            blank2(jdim,kdim,idim,2) = blank2(jdim-1,kdim,idim-1,2);
            blank2(1,kdim,1,2)       = blank2(2,kdim,2,2);

            // interior cells
            for (i = 2; i <= idim1; i++) {
                for (j = 2; j <= jdim1; j++) {
                    for (k = 2; k <= kdim1; k++) {
                        blank2(j,k,i,2) = std::min({(double)blank(j,k,i),
                                                    (double)blank(j-1,k,i),
                                                    (double)blank(j,k-1,i),
                                                    (double)blank(j-1,k-1,i),
                                                    (double)blank(j,k,i-1),
                                                    (double)blank(j-1,k,i-1),
                                                    (double)blank(j,k-1,i-1),
                                                    (double)blank(j-1,k-1,i-1)});
                    }
                }
            }

            // combining topology and overset grid blankings
            for (i = 1; i <= idim; i++) {
                for (j = 1; j <= jdim; j++) {
                    for (k = 1; k <= kdim; k++) {
                        blank2(j,k,i,1) = std::min((double)blank2(j,k,i,1),
                                                    (double)blank2(j,k,i,2));
                    }
                }
            }
        } // iover == 1

    } // end iflag==1 && iblnk>0



    // *** plot3d data ***
    if (iflag == 1) {

        if (lhdr > 0 && myid == myhost) {
            FILE* f11 = fortran_get_unit(11);
            if (i2d == 1) {
                fprintf(f11, "writing plot3d file for JDIM X KDIM =%5d x%5d grid\n",
                        jdim, kdim);
            } else {
                fprintf(f11, "writing plot3d file for IDIM X JDIM X KDIM =%5d x%5d x%5d grid\n",
                        idim, jdim, kdim);
            }
            if (inpl3d(iplot,2) == 0) {
                fprintf(f11, "   plot3dg file is an xyz file at grid points\n");
                fprintf(f11, "   plot3dq file is a    q file at grid points\n");
            }
            if (ifunc == -5) {
                fprintf(f11, "   plot3dg file is an xyz file at grid points\n");
                fprintf(f11, "   plot3dq file is a function file of Cp at grid points\n");
            }
            if (ifunc == -6) {
                fprintf(f11, "   plot3dg file is an xyz file at grid points\n");
                fprintf(f11, "   plot3dq file is a function file of P/Pinf at grid points\n");
            }
            if (i2d == 1) {
                if (iblnk == 0) {
                    fprintf(f11, "   plot3d files to be read with /mgrid/2d qualifiers\n");
                } else {
                    fprintf(f11, "   plot3d files to be read with /mgrid/blank/2d qualifiers\n");
                }
            } else {
                if (iblnk == 0) {
                    fprintf(f11, "   plot3d files to be read with /mgrid qualifiers\n");
                } else {
                    fprintf(f11, "   plot3d files to be read with /mgrid/blank qualifiers\n");
                }
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
                    xg(jw,kw,iw,4) = blank2(j,k,i,1);
                }
            }
        }

        if (myid == myhost && icall2d == 0) {
            // ialph > 0: grid was read in plot3d format with alpha in xy plane
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

            // output grid to unit 93
            FILE* f93 = fortran_get_unit(93);
            if (ibin == 0) {
                // formatted output
                if (i2d == 0) {
                    if (iblnk == 0) {
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,1));
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,2));
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,3));
                        fprintf(f93,"\n");
                    } else {
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,1));
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,2));
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,3));
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)(int)xg(j,k,i,4));
                        fprintf(f93,"\n");
                    }
                } else {
                    if (iblnk == 0) {
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,1));
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,3));
                        fprintf(f93,"\n");
                    } else {
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,1));
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)xg(j,k,i,3));
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f93,"%14.6e",(float)(int)xg(j,k,i,4));
                        fprintf(f93,"\n");
                    }
                }
            } else {
                // binary output
                if (i2d == 0) {
                    if (iblnk == 0) {
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,1); fwrite(&v,sizeof(float),1,f93); }
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,2); fwrite(&v,sizeof(float),1,f93); }
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,3); fwrite(&v,sizeof(float),1,f93); }
                    } else {
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,1); fwrite(&v,sizeof(float),1,f93); }
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,2); fwrite(&v,sizeof(float),1,f93); }
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,3); fwrite(&v,sizeof(float),1,f93); }
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { int iv=(int)xg(j,k,i,4); fwrite(&iv,sizeof(int),1,f93); }
                    }
                } else {
                    if (iblnk == 0) {
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,1); fwrite(&v,sizeof(float),1,f93); }
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,3); fwrite(&v,sizeof(float),1,f93); }
                    } else {
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,1); fwrite(&v,sizeof(float),1,f93); }
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xg(j,k,i,3); fwrite(&v,sizeof(float),1,f93); }
                        for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { int iv=(int)xg(j,k,i,4); fwrite(&iv,sizeof(int),1,f93); }
                    }
                }
            }
        } // myid==myhost && icall2d==0



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
            for (i = 1; i <= iw; i++) {
                for (k = 1; k <= kw; k++) {
                    for (j = 1; j <= jw; j++) {
                        xv(j,k,i,5) = xv(j,k,i,5)/gm1 + 0.5*(xv(j,k,i,2)*xv(j,k,i,2)
                                    + xv(j,k,i,3)*xv(j,k,i,3) + xv(j,k,i,4)*xv(j,k,i,4))*xv(j,k,i,1);
                        xv(j,k,i,2) = xv(j,k,i,1)*xv(j,k,i,2);
                        xv(j,k,i,3) = xv(j,k,i,1)*xv(j,k,i,3);
                        xv(j,k,i,4) = xv(j,k,i,1)*xv(j,k,i,4);
                    }
                }
            }
        } else {
            // load desired function in xv array
            if (ifunc == -5) {
                cpc = 2.e0/(gamma*xmach*xmach);
                for (i = 1; i <= iw; i++) {
                    for (k = 1; k <= kw; k++) {
                        for (j = 1; j <= jw; j++) {
                            xv(j,k,i,1) = (xv(j,k,i,5)/p0 - 1.)*cpc;
                        }
                    }
                }
            }
            if (ifunc == -6) {
                for (i = 1; i <= iw; i++) {
                    for (k = 1; k <= kw; k++) {
                        for (j = 1; j <= jw; j++) {
                            xv(j,k,i,1) = xv(j,k,i,5)/p0;
                        }
                    }
                }
            }
        }

        // load xv into single precision xw array
        leq = 5;
        if (ifunc != 0) leq = 1;
        for (l = 1; l <= leq; l++) {
            for (i = 1; i <= iw; i++) {
                for (k = 1; k <= kw; k++) {
                    for (j = 1; j <= jw; j++) {
                        xw(j,k,i,l) = xv(j,k,i,l);
                    }
                }
            }
        }

        if (myid == myhost) {
            // output solution to unit 94
            FILE* f94 = fortran_get_unit(94);
            if (ifunc == 0) {
                if (ibin == 0) {
                    // formatted
                    fprintf(f94, "%14.6e%14.6e%14.6e%14.6e\n",
                            (float)xmachw, (float)alphww, (float)reuew, (float)timew);
                    if (i2d == 0) {
                        if (ialph == 0) {
                            for (m=1;m<=5;m++) for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f94,"%14.6e",(float)xw(j,k,i,m));
                            fprintf(f94,"\n");
                        } else {
                            for (i=1;i<=iw;i++) for (j=1;j<=jw;j++) for (k=1;k<=kw;k++) xw(j,k,i,3) = -xw(j,k,i,3);
                            for (m=1;m<=2;m++) for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f94,"%14.6e",(float)xw(j,k,i,m));
                            for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f94,"%14.6e",(float)xw(j,k,i,4));
                            for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f94,"%14.6e",(float)xw(j,k,i,3));
                            for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f94,"%14.6e",(float)xw(j,k,i,5));
                            fprintf(f94,"\n");
                        }
                    } else {
                        for (m=1;m<=2;m++) for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f94,"%14.6e",(float)xw(j,k,i,m));
                        for (m=4;m<=5;m++) for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f94,"%14.6e",(float)xw(j,k,i,m));
                        fprintf(f94,"\n");
                    }
                } else {
                    // binary
                    { float v=(float)xmachw; fwrite(&v,sizeof(float),1,f94); }
                    { float v=(float)alphww;  fwrite(&v,sizeof(float),1,f94); }
                    { float v=(float)reuew;   fwrite(&v,sizeof(float),1,f94); }
                    { float v=(float)timew;   fwrite(&v,sizeof(float),1,f94); }
                    if (i2d == 0) {
                        if (ialph == 0) {
                            for (m=1;m<=5;m++) for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xw(j,k,i,m); fwrite(&v,sizeof(float),1,f94); }
                        } else {
                            for (i=1;i<=iw;i++) for (j=1;j<=jw;j++) for (k=1;k<=kw;k++) xw(j,k,i,3) = -xw(j,k,i,3);
                            for (m=1;m<=2;m++) for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xw(j,k,i,m); fwrite(&v,sizeof(float),1,f94); }
                            for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xw(j,k,i,4); fwrite(&v,sizeof(float),1,f94); }
                            for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xw(j,k,i,3); fwrite(&v,sizeof(float),1,f94); }
                            for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xw(j,k,i,5); fwrite(&v,sizeof(float),1,f94); }
                        }
                    } else {
                        for (m=1;m<=2;m++) for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xw(j,k,i,m); fwrite(&v,sizeof(float),1,f94); }
                        for (m=4;m<=5;m++) for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xw(j,k,i,m); fwrite(&v,sizeof(float),1,f94); }
                    }
                }
            } else {
                // ifunc != 0: write only xw(1)
                if (ibin == 0) {
                    for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) fprintf(f94,"%14.6e",(float)xw(j,k,i,1));
                    fprintf(f94,"\n");
                } else {
                    for (k=1;k<=kw;k++) for (j=1;j<=jw;j++) for (i=1;i<=iw;i++) { float v=(float)xw(j,k,i,1); fwrite(&v,sizeof(float),1,f94); }
                }
            }
        } // myid==myhost

    } // end iflag==1



    // *** print out data ***
    if (iflag == 2) {

        if (lhdr > 0 && myid == myhost) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, "writing printout file for IDIM X JDIM X KDIM =%5d x%5d x%5d grid\n",
                    idim, jdim, kdim);
        }

        // determine q values at grid points and load into single precision array
        jdw = (j2-j1)/j3 + 1;
        kdw = (k2-k1)/k3 + 1;
        idw = (i2-i1)/i3 + 1;
        ldw = 5;
        cctogp_ns::cctogp(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                          q, qi0, qj0, qk0, jdw, kdw, idw, xv, ldw);
        for (l = 1; l <= ldw; l++) {
            for (i = 1; i <= idw; i++) {
                for (k = 1; k <= kdw; k++) {
                    for (j = 1; j <= jdw; j++) {
                        xw(j,k,i,l) = xv(j,k,i,l);
                    }
                }
            }
        }

        // determine eddy viscosity values at grid points
        if (cmn_reyue.ivisc[2] > 1 || cmn_reyue.ivisc[1] > 1 || cmn_reyue.ivisc[0] > 1) {
            ldw = 1;
            // Fortran: call cctogp(...,vist3d,vi0,vj0,vk0,...,xv(1,1,1,4),ldw)
            // vist3d(jdim,kdim,idim) passed as 4D dummy -> treat as (jdim,kdim,idim,1)
            // xv(1,1,1,4) is the start of the 4th "layer" of xv
            FortranArray4DRef<double> vist3d_4d(&vist3d(1,1,1), jdim, kdim, idim, 1);
            FortranArray4DRef<double> xv4(&xv(1,1,1,4), jdw, kdw, idw, 1);
            cctogp_ns::cctogp(jdim, kdim, idim, i1, i2, i3, j1, j2, j3, k1, k2, k3,
                              vist3d_4d, vi0, vj0, vk0, jdw, kdw, idw, xv4, ldw);
            for (i = 1; i <= idw; i++) {
                for (k = 1; k <= kdw; k++) {
                    for (j = 1; j <= jdw; j++) {
                        xg(j,k,i,4) = xv(j,k,i,4);
                    }
                }
            }
        }

        // old way for eddy viscosity (iold=0, so skipped)
        iold = 0;
        if (iold > 0) {
            if (cmn_reyue.ivisc[2] > 1 || cmn_reyue.ivisc[1] > 1 || cmn_reyue.ivisc[0] > 1) {
                iw = 0;
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    jw = 0;
                    for (j = j1; j <= j2; j += j3) {
                        jw = jw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jv = j;
                            kv = k;
                            iv = i;
                            jv = std::min(jdim-1, jv);
                            kv = std::min(kdim-1, kv);
                            iv = std::min(idim-1, iv);
                            jx = j-1;
                            kx = k-1;
                            ix = i-1;
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
                        "           P/Pinf           T/Tinf            MACH"
                        "              cp           tur. vis.\n");
            } else {
                fprintf(f17, "\n   I   J   K         X                Y                Z"
                        "              U/Uinf           V/Vinf           W/Winf"
                        "           P/Pinf           T/Tinf            MACH"
                        "            pitotp           tur. vis.\n");
            }

            term3 = 1./(std::pow(1.+0.5*gm1*xmach*xmach, gamma/gm1));
            iw = 0;
            for (i = i1; i <= i2; i += i3) {
                iw = iw + 1;
                kw = 0;
                for (k = k1; k <= k2; k += k3) {
                    kw = kw + 1;
                    jw = 0;
                    for (j = j1; j <= j2; j += j3) {
                        jw = jw + 1;
                        q1    = xw(jw,kw,iw,1);
                        q2    = xw(jw,kw,iw,2)/xmach;
                        q3    = xw(jw,kw,iw,3)/xmach;
                        q4    = xw(jw,kw,iw,4)/xmach;
                        q5    = gamma*xw(jw,kw,iw,5);
                        t1    = q5/q1;
                        xm1   = std::sqrt(xmach*xmach*(q2*q2+q3*q3+q4*q4)/t1);
                        edvis = xg(jw,kw,iw,4);
                        if (icp == 1) {
                            pitot = 2.*(q5-1.)/(gamma*xmach*xmach);
                        } else {
                            if (xm1 > 1.0) {
                                term1 = std::pow(0.5*gp1*xm1*xm1, gamma/gm1);
                                term2 = std::pow(2.*gamma*xm1*xm1/gp1 - gm1/gp1, 1./gm1);
                                pitot = q5*term1*term3/term2;
                            } else {
                                term1 = std::pow(1.0+0.5*gm1*xm1*xm1, gamma/gm1);
                                pitot = q5*term1*term3;
                            }
                        }
                        if (ialph == 0) {
                            fprintf(f17, "%4d%4d%4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                                    i,j,k,
                                    (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,3),
                                    (float)q2,(float)q3,(float)q4,(float)q5,(float)t1,(float)xm1,(float)pitot,(float)edvis);
                        } else {
                            xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                            q3 = -q3;
                            fprintf(f17, "%4d%4d%4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                                    i,j,k,
                                    (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,2),
                                    (float)q2,(float)q4,(float)q3,(float)q5,(float)t1,(float)xm1,(float)pitot,(float)edvis);
                        }
                    }
                }
            }
        } // myid==myhost



        // output viscous-flow data on k=1 and/or k=kdim surfaces
        if (cmn_reyue.ivisc[2] != 0) {
            iw = 0;
            for (i = i1; i <= i2; i += i3) {
                iw = iw + 1;
                id   = i;
                id1  = id-1;
                if (id1 < 1)    id1 = 1;
                if (id > idim1) id  = idim1;
                kw = 0;
                for (k = k1; k <= k2; k += k3) {
                    kw = kw + 1;
                    jw = 0;
                    if (k != 1 && k != kdim) continue; // go to 6001
                    for (j = j1; j <= j2; j += j3) {
                        jw = jw + 1;
                        if (j == 1 || j == jdim) continue; // go to 6000
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
                        if (sk(j  ,kdx,id ,4) == 0. || sk(j-1,kdx,id ,4) == 0. ||
                            sk(j  ,kdx,id1,4) == 0. || sk(j-1,kdx,id1,4) == 0.) {
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
                        if (cmn_wallfun.iwf[2] == 1) {
                            avgmut = .25*(vk0(j,id,1,m)  +vk0(j-1,id,1,m)
                                        + vk0(j,id1,1,m) +vk0(j-1,id1,1,m));
                        } else {
                            avgmut = 0.;
                        }
                        emuka = std::pow(t1k1,1.5)*((1.0+cbar/tinf)/(t1k1+cbar/tinf));
                        // velocity parallel to wall
                        urel = q2k2-q2k1;
                        vrel = q3k2-q3k1;
                        wrel = q4k2-q4k1;
                        sk1 = (sk(j  ,kdx,id ,1)+sk(j-1,kdx,id ,1)+
                               sk(j  ,kdx,id1,1)+sk(j-1,kdx,id1,1))/4.;
                        sk2 = (sk(j  ,kdx,id ,2)+sk(j-1,kdx,id ,2)+
                               sk(j  ,kdx,id1,2)+sk(j-1,kdx,id1,2))/4.;
                        sk3 = (sk(j  ,kdx,id ,3)+sk(j-1,kdx,id ,3)+
                               sk(j  ,kdx,id1,3)+sk(j-1,kdx,id1,3))/4.;
                        vnorm = (urel*sk1 + vrel*sk2 + wrel*sk3)*sgn;
                        upar = urel-vnorm*sk1*sgn;
                        vpar = vrel-vnorm*sk2*sgn;
                        wpar = wrel-vnorm*sk3*sgn;
                        Cf   = 2.*std::sqrt(upar*upar + vpar*vpar + wpar*wpar);
                        Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn;
                        if (q2k2 < 0.) Cf = -Cf;
                        cfx = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn;
                        cfy = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn;
                        cfz = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn;
                        eps  = 1.0e-6;
                        tty1 = 1.0 + gm1*0.5*xmach*xmach;
                        t1m1 = t1k1-tty1;
                        if (std::abs(t1m1) <= eps) {
                            Ch = 999.99999;
                        } else {
                            Ch = 2.*(t1k2-t1k1);
                            Ch = (emuka+avgmut)/(reue*pr*(t1k1-tty1))*Ch/dn;
                        }
                        pres1 = q5k1;
                        temp1 = t1k1;
                        yplus = 0.;
                        if (std::abs(Cf) > 0.) {
                            yplus = dn*reue*q1k1*std::sqrt(std::abs(Cf*0.5/q1k1))/emuka;
                        }
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
                    } // j loop (6000)
                } // k loop (6001)
            } // i loop (6002)



            // write printout data to file (k-surface)
            if (myid == myhost) {
                isay = 0;
                iw = 0;
                FILE* f17 = fortran_get_unit(17);
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    kw = 0;
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        jw = 0;
                        if (k != 1 && k != kdim) continue;
                        for (j = j1; j <= j2; j += j3) {
                            jw = jw + 1;
                            if (j == 1 || j == jdim) continue;
                            isay = isay+1;
                            if (isay == 1)
                                fprintf(f17, "\n   I   J   K         X                Y                Z"
                                        "               dn           P/Pinf           T/Tinf"
                                        "               Cf               Ch            yplus"
                                        "              Cfx              Cfy              Cfz\n");
                            if (ialph == 0) {
                                fprintf(f17, "%4d%4d%4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                                        i,j,k,
                                        (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,7),(float)xw(jw,kw,iw,8));
                            } else {
                                xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                                xw(jw,kw,iw,7) = -xw(jw,kw,iw,7);
                                fprintf(f17, "%4d%4d%4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                                        i,j,k,
                                        (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,8),(float)xw(jw,kw,iw,7));
                            }
                        }
                    }
                }
            } // myid==myhost
        } // ivisc(3) != 0

        // output viscous-flow data on j=1 and/or j=jdim surfaces
        if (cmn_reyue.ivisc[1] != 0) {
            iw = 0;
            for (i = i1; i <= i2; i += i3) {
                iw = iw + 1;
                id   = i;
                id1  = id-1;
                if (id1 < 1)    id1 = 1;
                if (id > idim1) id  = idim1;
                jw = 0;
                for (j = j1; j <= j2; j += j3) {
                    jw = jw + 1;
                    kw = 0;
                    if (j != 1 && j != jdim) continue; // go to 7001
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        if (k == 1 || k == kdim) continue; // go to 7000
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
                        if (sj(jdx,k  ,id ,4) == 0. || sj(jdx,k-1,id ,4) == 0. ||
                            sj(jdx,k  ,id1,4) == 0. || sj(jdx,k-1,id1,4) == 0.) {
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
                        if (cmn_wallfun.iwf[1] == 1) {
                            avgmut = .25*(vj0(k,id,1,m)  +vj0(k-1,id,1,m)
                                        + vj0(k,id1,1,m) +vj0(k-1,id1,1,m));
                        } else {
                            avgmut = 0.;
                        }
                        emuka = std::pow(t1j1,1.5)*((1.0+cbar/tinf)/(t1j1+cbar/tinf));
                        urel = q2j2-q2j1;
                        vrel = q3j2-q3j1;
                        wrel = q4j2-q4j1;
                        sj1 = (sj(jdx,k  ,id ,1)+sj(jdx,k-1,id ,1)+
                               sj(jdx,k  ,id1,1)+sj(jdx,k-1,id1,1))/4.;
                        sj2 = (sj(jdx,k  ,id ,2)+sj(jdx,k-1,id ,2)+
                               sj(jdx,k  ,id1,2)+sj(jdx,k-1,id1,2))/4.;
                        sj3 = (sj(jdx,k  ,id ,3)+sj(jdx,k-1,id ,3)+
                               sj(jdx,k  ,id1,3)+sj(jdx,k-1,id1,3))/4.;
                        vnorm = (urel*sj1 + vrel*sj2 + wrel*sj3)*sgn;
                        upar = urel-vnorm*sj1*sgn;
                        vpar = vrel-vnorm*sj2*sgn;
                        wpar = wrel-vnorm*sj3*sgn;
                        Cf   = 2.*std::sqrt(upar*upar + vpar*vpar + wpar*wpar);
                        Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn;
                        if (q2j2 < 0.) Cf = -Cf;
                        cfx = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn;
                        cfy = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn;
                        cfz = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn;
                        eps  = 1.0e-6;
                        tty1 = 1.0 + gm1*0.5*xmach*xmach;
                        t1m1 = t1j1-tty1;
                        if (std::abs(t1m1) <= eps) {
                            Ch = 999.99999;
                        } else {
                            Ch = 2.*(t1j2-t1j1);
                            Ch = (emuka+avgmut)/(reue*pr*(t1j1-tty1))*Ch/dn;
                        }
                        pres1 = q5j1;
                        temp1 = t1j1;
                        yplus = 0.;
                        if (std::abs(Cf) > 0.) {
                            yplus = dn*reue*q1j1*std::sqrt(std::abs(Cf*0.5/q1j1))/emuka;
                        }
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
                    } // k loop (7000)
                } // j loop (7001)
            } // i loop (7002)



            // write printout data to file (j-surface)
            if (myid == myhost) {
                isay = 0;
                iw = 0;
                FILE* f17 = fortran_get_unit(17);
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    jw = 0;
                    for (j = j1; j <= j2; j += j3) {
                        jw = jw + 1;
                        kw = 0;
                        if (j != 1 && j != jdim) continue;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            if (k == 1 || k == kdim) continue;
                            isay = isay+1;
                            if (isay == 1)
                                fprintf(f17, "\n   I   J   K         X                Y                Z"
                                        "               dn           P/Pinf           T/Tinf"
                                        "               Cf               Ch            yplus"
                                        "              Cfx              Cfy              Cfz\n");
                            if (ialph == 0) {
                                fprintf(f17, "%4d%4d%4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                                        i,j,k,
                                        (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,7),(float)xw(jw,kw,iw,8));
                            } else {
                                xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                                xw(jw,kw,iw,7) = -xw(jw,kw,iw,7);
                                fprintf(f17, "%4d%4d%4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                                        i,j,k,
                                        (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,8),(float)xw(jw,kw,iw,7));
                            }
                        }
                    }
                }
            } // myid==myhost
        } // ivisc(2) != 0

        // i-surface viscous data: return early if idim==2
        if (idim == 2) return;

        if (cmn_reyue.ivisc[0] != 0) {
            jw = 0;
            for (j = j1; j <= j2; j += j3) {
                jw = jw + 1;
                iw = 0;
                if (j == 1 || j == jdim) continue; // go to 8002
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    kw = 0;
                    if (i != 1 && i != idim) continue; // go to 8001
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        if (k == 1 || k == kdim) continue; // go to 8000
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
                        if (si(j  ,k  ,idx,4) == 0. || si(j  ,k-1,idx,4) == 0. ||
                            si(j-1,k  ,idx,4) == 0. || si(j-1,k-1,idx,4) == 0.) {
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
                        if (cmn_wallfun.iwf[0] == 1) {
                            avgmut = .25*(vi0(j,k,1,m)   +vi0(j-1,k,1,m)
                                        + vi0(j,k-1,1,m) +vi0(j-1,k-1,1,m));
                        } else {
                            avgmut = 0.;
                        }
                        emuka = std::pow(t1i1,1.5)*((1.0+cbar/tinf)/(t1i1+cbar/tinf));
                        urel = q2i2-q2i1;
                        vrel = q3i2-q3i1;
                        wrel = q4i2-q4i1;
                        si1 = (si(j  ,k  ,idx,1)+si(j  ,k-1,idx,1)+
                               si(j-1,k  ,idx,1)+si(j-1,k-1,idx,1))/4.;
                        si2 = (si(j  ,k  ,idx,2)+si(j  ,k-1,idx,2)+
                               si(j-1,k  ,idx,2)+si(j-1,k-1,idx,2))/4.;
                        si3 = (si(j  ,k  ,idx,3)+si(j  ,k-1,idx,3)+
                               si(j-1,k  ,idx,3)+si(j-1,k-1,idx,3))/4.;
                        vnorm = (urel*si1 + vrel*si2 + wrel*si3)*sgn;
                        upar = urel-vnorm*si1*sgn;
                        vpar = vrel-vnorm*si2*sgn;
                        wpar = wrel-vnorm*si3*sgn;
                        Cf   = 2.*std::sqrt(upar*upar + vpar*vpar + wpar*wpar);
                        Cf = 2.0*(emuka+avgmut)/(reue*xmach)*Cf/dn;
                        if (q2i2 < 0.) Cf = -Cf;
                        cfx = 4.*(emuka+avgmut)/(reue*xmach)*upar/dn;
                        cfy = 4.*(emuka+avgmut)/(reue*xmach)*vpar/dn;
                        cfz = 4.*(emuka+avgmut)/(reue*xmach)*wpar/dn;
                        eps  = 1.0e-6;
                        tty1 = 1.0 + gm1*0.5*xmach*xmach;
                        t1m1 = t1i1-tty1;
                        if (std::abs(t1m1) <= eps) {
                            Ch = 999.99999;
                        } else {
                            Ch = 2.*(t1i2-t1i1);
                            Ch = (emuka+avgmut)/(reue*pr*(t1i1-tty1))*Ch/dn;
                        }
                        pres1 = q5i1;
                        temp1 = t1i1;
                        yplus = 0.;
                        if (std::abs(Cf) > 0.) {
                            yplus = dn*reue*q1i1*std::sqrt(std::abs(Cf*0.5/q1i1))/emuka;
                        }
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
                    } // k loop (8000)
                } // i loop (8001)
            } // j loop (8002)



            // write printout data to file (i-surface)
            if (myid == myhost) {
                isay = 0;
                jw = 0;
                FILE* f17 = fortran_get_unit(17);
                for (j = j1; j <= j2; j += j3) {
                    jw = jw + 1;
                    iw = 0;
                    if (j == 1 || j == jdim) continue;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        if (i != 1 && i != idim) continue;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            if (k == 1 || k == kdim) continue;
                            isay = isay+1;
                            if (isay == 1)
                                fprintf(f17, "\n   I   J   K         X                Y                Z"
                                        "               dn           P/Pinf           T/Tinf"
                                        "               Cf               Ch            yplus"
                                        "              Cfx              Cfy              Cfz\n");
                            if (ialph == 0) {
                                fprintf(f17, "%4d%4d%4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                                        i,j,k,
                                        (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,7),(float)xw(jw,kw,iw,8));
                            } else {
                                xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                                xw(jw,kw,iw,7) = -xw(jw,kw,iw,7);
                                fprintf(f17, "%4d%4d%4d %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n",
                                        i,j,k,
                                        (float)xg(jw,kw,iw,1),(float)xg(jw,kw,iw,3),(float)xg(jw,kw,iw,2),(float)xg(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,1),(float)xw(jw,kw,iw,2),(float)xw(jw,kw,iw,3),(float)xw(jw,kw,iw,4),
                                        (float)xw(jw,kw,iw,5),(float)xw(jw,kw,iw,6),(float)xw(jw,kw,iw,8),(float)xw(jw,kw,iw,7));
                            }
                        }
                    }
                }
            } // myid==myhost
        } // ivisc(1) != 0

    } // end iflag==2

} // end plot3d_2d function

} // namespace plot3d_2d_ns
