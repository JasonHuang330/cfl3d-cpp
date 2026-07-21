// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "plot3c.h"
#include "runtime/fortran_io.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace plot3c_ns {

void plot3c(int& jdim, int& kdim, int& idim, int& i1, int& i2, int& i3, int& j1, int& j2, int& j3, int& k1, int& k2, int& k3, FortranArray4DRef<double> q, FortranArray4DRef<double> qi0, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> xw, FortranArray3DRef<double> blank2, FortranArray3DRef<double> blank, FortranArray4DRef<double> xg, int& iflag, FortranArray3DRef<double> vist3d, FortranArray4DRef<double> vi0, FortranArray4DRef<double> vj0, FortranArray4DRef<double> vk0, int& iover, int& nblk, FortranArray1DRef<int> nmap, FortranArray3DRef<double> smin, int& ifunc, int& iplot, int& jdw, int& kdw, int& idw, int& nplots, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nblcg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> isg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, int& ninter, FortranArray2DRef<int> iindex, int& intmax, int& nsub1, int& maxxe, FortranArray2DRef<int> nblkk, int& nbli, FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva, FortranArray1DRef<int> nblon, int& mxbli, FortranArray1DRef<double> thetay, int& maxbl, int& maxgr, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> inpl3d, int& nblock, FortranArray1DRef<int> nblkpt, int& ipsurf)
{
    // COMMON block aliases
    int&   ibin      = cmn_bin.ibin;
    int&   iblnk     = cmn_bin.iblnk;
    float& gamma_f   = cmn_fluid.gamma;
    float& gm1_f     = cmn_fluid.gm1;
    float& gp1_f     = cmn_fluid.gp1;
    float& xmach_f   = cmn_info.xmach;
    float& alpha_f   = cmn_info.alpha;
    float& reue_f    = cmn_reyue.reue;
    float& tinf_f    = cmn_reyue.tinf;
    int&   i2d       = cmn_twod.i2d;
    float& time_f    = cmn_unst.time;
    int&   ialph     = cmn_igrdtyp.ialph;
    int&   icall1    = cmn_moov.icall1;
    int&   lhdr      = cmn_moov.lhdr;
    float& radtodeg  = cmn_conversion.radtodeg;
    float& p0_f      = cmn_ivals.p0;
    float& cbar_f    = cmn_fluid2.cbar;

    // Promote to double for arithmetic
    double gamma    = (double)gamma_f;
    double gm1      = (double)gm1_f;
    double gp1      = (double)gp1_f;
    double xmach    = (double)xmach_f;
    double alpha    = (double)alpha_f;
    double reue     = (double)reue_f;
    double tinf     = (double)tinf_f;
    double time_v   = (double)time_f;
    double p0       = (double)p0_f;
    double cbar     = (double)cbar_f;

    // Local variables
    int jdim1 = jdim - 1;
    int kdim1 = kdim - 1;
    int idim1 = idim - 1;
    int jw, kw, iw;
    int i, j, k, l, n;
    double alphaw, xmachw, alphww, reuew, timew;
    double temp;
    int mm;
    double cpc;
    double vor, tt, c2b, c2bp, xmu, xmut, xnuwiggle, delta_v;
    double term1, term2, term3, pitot, q1, q2, q3, q4, q5, t1, xm1, edvis;
    int icp;
    int inter, lmax1, nbl, lst, lcoord, lend, j21, j22, k21, k22, ll, mblk;
    int it, ir, itime, iti, is, ie, js, je, ks, ke;
    int is1, ie1, js1, je1, ks1, ke1;
    int ibembed;
    int nblc, nblcc;


    // initialize xw and xg arrays
    jw = (j2-j1)/j3 + 1;
    kw = (k2-k1)/k3 + 1;
    iw = (i2-i1)/i3 + 1;
    for (j = 1; j <= jw; j++) {
        for (k = 1; k <= kw; k++) {
            for (i = 1; i <= iw; i++) {
                for (l = 1; l <= 5; l++) {
                    xw(j,k,i,l) = 0.;
                }
                for (l = 1; l <= 4; l++) {
                    xg(j,k,i,l) = 0.;
                }
            }
        }
    }

    // assign single precision scalars
    alphaw = (double)radtodeg * (alpha + (double)thetay(nblk));
    xmachw = xmach;
    alphww = alphaw;
    reuew  = reue;
    timew  = time_v;

    // average grid points to get cell centers and load into first 3
    // locations of the single precision xg array
    if (ipsurf > 0) {
        // surface face centers
        if (i1 == i2) {
            iw = 1;
            if (i1 == 1) {
                i = 1;
            } else {
                i = idim;
            }
            kw = 0;
            for (k = k1; k <= k2; k += k3) {
                kw = kw + 1;
                jw = 0;
                for (j = j1; j <= j2; j += j3) {
                    jw = jw + 1;
                    xg(jw,kw,1,1) = .25*(   x(j,k,i)   + x(j+1,k,i)
                                           + x(j,k+1,i) + x(j+1,k+1,i));
                    xg(jw,kw,1,2) = .25*(   y(j,k,i)   + y(j+1,k,i)
                                           + y(j,k+1,i) + y(j+1,k+1,i));
                    xg(jw,kw,1,3) = .25*(   z(j,k,i)   + z(j+1,k,i)
                                           + z(j,k+1,i) + z(j+1,k+1,i));
                }
            }
        } else if (j1 == j2) {
            jw = 1;
            if (j1 == 1) {
                j = 1;
            } else {
                j = jdim;
            }
            iw = 0;
            for (i = i1; i <= i2; i += i3) {
                iw = iw + 1;
                kw = 0;
                for (k = k1; k <= k2; k += k3) {
                    kw = kw + 1;
                    xg(1,kw,iw,1) = .25*(   x(j,k,i)   + x(j,k,i+1)
                                           + x(j,k+1,i) + x(j,k+1,i+1));
                    xg(1,kw,iw,2) = .25*(   y(j,k,i)   + y(j,k,i+1)
                                           + y(j,k+1,i) + y(j,k+1,i+1));
                    xg(1,kw,iw,3) = .25*(   z(j,k,i)   + z(j,k,i+1)
                                           + z(j,k+1,i) + z(j,k+1,i+1));
                }
            }
        } else if (k1 == k2) {
            kw = 1;
            if (k1 == 1) {
                k = 1;
            } else {
                k = kdim;
            }
            iw = 0;
            for (i = i1; i <= i2; i += i3) {
                iw = iw + 1;
                jw = 0;
                for (j = j1; j <= j2; j += j3) {
                    jw = jw + 1;
                    xg(jw,1,iw,1) = .25*(   x(j,k,i)   + x(j,k,i+1)
                                           + x(j+1,k,i) + x(j+1,k,i+1));
                    xg(jw,1,iw,2) = .25*(   y(j,k,i)   + y(j,k,i+1)
                                           + y(j+1,k,i) + y(j+1,k,i+1));
                    xg(jw,1,iw,3) = .25*(   z(j,k,i)   + z(j,k,i+1)
                                           + z(j+1,k,i) + z(j+1,k,i+1));
                }
            }
        }
    } else {
        // cell centers
        iw = 0;
        for (i = i1; i <= i2; i += i3) {
            iw = iw + 1;
            kw = 0;
            for (k = k1; k <= k2; k += k3) {
                kw = kw + 1;
                jw = 0;
                for (j = j1; j <= j2; j += j3) {
                    jw = jw + 1;
                    xg(jw,kw,iw,1) = .125*(   x(j,k,i)    +x(j,k,i+1)
                                             + x(j+1,k,i)  +x(j+1,k,i+1)
                                             + x(j,k+1,i)  +x(j,k+1,i+1)
                                             + x(j+1,k+1,i)+x(j+1,k+1,i+1));
                    xg(jw,kw,iw,2) = .125*(   y(j,k,i)    +y(j,k,i+1)
                                             + y(j+1,k,i)  +y(j+1,k,i+1)
                                             + y(j,k+1,i)  +y(j,k+1,i+1)
                                             + y(j+1,k+1,i)+y(j+1,k+1,i+1));
                    xg(jw,kw,iw,3) = .125*(   z(j,k,i)    +z(j,k,i+1)
                                             + z(j+1,k,i)  +z(j+1,k,i+1)
                                             + z(j,k+1,i)  +z(j,k+1,i+1)
                                             + z(j+1,k+1,i)+z(j+1,k+1,i+1));
                }
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
                    blank2(j,k,i) = 1.;
                }
            }
        }

        // iblank (blank2) array for generalized patch interface boundaries
        if (std::abs(ninter) > 0) {
            for (inter = 1; inter <= std::abs(ninter); inter++) {
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
                    if (lend == 2) i = idim1;
                    for (j = j21; j <= j22-1; j++) {
                        for (k = k21; k <= k22-1; k++) {
                            ll = lst + (j22-j21)*(k-k21) + (j-j21);
                            mblk = iindex(inter,nblkpt(ll)+1);
                            blank2(j,k,i) = -(double)(nmap(mblk));
                        }
                    }
                }
                if (lcoord == 2) {
                    if (lend == 1) j = 1;
                    if (lend == 2) j = jdim1;
                    for (i = k21; i <= k22-1; i++) {
                        for (k = j21; k <= j22-1; k++) {
                            ll = lst + (j22-j21)*(i-k21) + (k-j21);
                            mblk = iindex(inter,nblkpt(ll)+1);
                            blank2(j,k,i) = -(double)(nmap(mblk));
                        }
                    }
                }
                if (lcoord == 3) {
                    if (lend == 1) k = 1;
                    if (lend == 2) k = kdim1;
                    for (i = k21; i <= k22-1; i++) {
                        for (j = j21; j <= j22-1; j++) {
                            ll = lst + (j22-j21)*(i-k21) + (j-j21);
                            mblk = iindex(inter,nblkpt(ll)+1);
                            blank2(j,k,i) = -(double)(nmap(mblk));
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
                            for (i = is1; i <= ie1; i++) {
                                for (j = js1; j <= je1; j++) {
                                    for (k = ks1; k <= ke1; k++) {
                                        blank2(j,k,i) = -(double)(nmap(nblkk(ir,n)));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        // iblank (blank2) array for embedded grids - the underlying
        // coarse grid areas are blanked out if the parameter ibembed > 0
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
                    for (i = is; i <= ie; i++) {
                        for (j = js; j <= je; j++) {
                            for (k = ks; k <= ke; k++) {
                                blank2(j,k,i) = 0.;
                            }
                        }
                    }
                }
            }
        }

        // iblank (blank2) array for overlapped grids
        if (iover == 1) {
            for (i = 1; i <= idim1; i++) {
                for (j = 1; j <= jdim1; j++) {
                    for (k = 1; k <= kdim1; k++) {
                        if (blank(j,k,i) == 0.0) blank2(j,k,i) = 0.;
                    }
                }
            }
        }

    } // end if (iflag == 1 && iblnk > 0)

    if (iflag == 1) {
        // plot3d data

        if (lhdr > 0 && myid == myhost) {
            FILE* f11 = fortran_get_unit(11);
            if (i2d == 1) {
                fprintf(f11, "writing plot3d file for JDIM X KDIM =%5d x%5d grid\n", jdim, kdim);
            } else {
                fprintf(f11, "writing plot3d file for IDIM X JDIM X KDIM =%5d x %5d x %5d grid\n",
                        idim, jdim, kdim);
            }
            if (inpl3d(iplot,2) == 1) {
                if (ipsurf == 0) {
                    fprintf(f11, "   plot3dg file is an xyz file at cell centers\n");
                    fprintf(f11, "   plot3dq file is a    q file at cell centers\n");
                } else {
                    fprintf(f11, "   plot3dg file is an xyz file at cell face centers\n");
                    fprintf(f11, "   plot3dq file is a    q file at cell face centers\n");
                }
            }
            if (ifunc == 3) {
                if (ipsurf == 0) {
                    fprintf(f11, "   plot3dg file is an xyz file at cell centers\n");
                    fprintf(f11, "   plot3dq file is a function file of minimum distance at cell centers\n");
                } else {
                    fprintf(f11, "   plot3dg file is an xyz file at cell face centers\n");
                    fprintf(f11, "   plot3dq file is a function file of minimum distance at cell face centers\n");
                }
            }
            if (ifunc == 4) {
                if (ipsurf == 0) {
                    fprintf(f11, "    plot3dg file is an xyz file at cell centers\n");
                    fprintf(f11, "    plot3dq file is a function file of mu_t/mu_lam  at cell centers\n");
                } else {
                    fprintf(f11, "   plot3dg file is an xyz file at cell face centers\n");
                    fprintf(f11, "   plot3dq file is a function file of mu_t/mu_lam at cell face centers\n");
                }
            }
            if (ifunc == 5) {
                if (ipsurf == 0) {
                    fprintf(f11, "   plot3dg file is an xyz file at cell centers\n");
                    fprintf(f11, "   plot3dq file is a function file of Cp at cell centers\n");
                } else {
                    fprintf(f11, "   plot3dg file is an xyz file at cell face centers\n");
                    fprintf(f11, "   plot3dq file is a function file of Cp at cell face centers\n");
                }
            }
            if (ifunc == 6) {
                if (ipsurf == 0) {
                    fprintf(f11, "   plot3dg file is an xyz file at cell centers\n");
                    fprintf(f11, "   plot3dq file is a function file of P/Pinf at cell centers\n");
                } else {
                    fprintf(f11, "   plot3dg file is an xyz file at cell face centers\n");
                    fprintf(f11, "   plot3dq file is a function file of P/Pinf at cell face centers\n");
                }
            }
            if (ifunc == 7) {
                if (ipsurf == 0) {
                    fprintf(f11, "   plot3dg file is an xyz file at cell centers\n");
                    fprintf(f11, "   plot3dq file is a function file of turbulence index at the wall\n");
                } else {
                    fprintf(f11, "   plot3dg file is an xyz file at cell face centers\n");
                    fprintf(f11, "   plot3dq file is a function file of turbulence index at the wall\n");
                }
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
                    xg(jw,kw,iw,4) = blank2(j,k,i);
                }
            }
        }

        if (myid == myhost && icall1 == 0) {
            // ialph > 0 for a grid that was read in plot3d format with alpha measured
            //           in the xy plane (TLNS3D convention)
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

            // output grid
            FILE* f3 = fortran_get_unit(3);
            if (ibin == 0) {
                if (i2d == 0) {
                    if (iblnk == 0) {
                        {
                            int cnt = 0;
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,1));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,2));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,3));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            if (cnt % 5 != 0) fprintf(f3, "\n");
                        }
                    } else {
                        {
                            int cnt = 0;
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,1));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,2));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,3));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)((int)xg(j,k,i,4)));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            if (cnt % 5 != 0) fprintf(f3, "\n");
                        }
                    }
                } else {
                    // i2d != 0, formatted
                    if (iblnk == 0) {
                        {
                            int cnt = 0;
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,1));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,3));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            if (cnt % 5 != 0) fprintf(f3, "\n");
                        }
                    } else {
                        {
                            int cnt = 0;
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,1));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)xg(j,k,i,3));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        fprintf(f3, "%14.6E", (float)((int)xg(j,k,i,4)));
                                        if (++cnt % 5 == 0) fprintf(f3, "\n");
                                    }
                            if (cnt % 5 != 0) fprintf(f3, "\n");
                        }
                    }
                }
            } else {

                // ibin != 0 (binary)
                if (i2d == 0) {
                    if (iblnk == 0) {
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,1);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,2);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,3);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                    } else {
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,1);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,2);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,3);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    int iv = (int)xg(j,k,i,4);
                                    fwrite(&iv, sizeof(int), 1, f3);
                                }
                    }
                } else {
                    // i2d != 0, binary
                    if (iblnk == 0) {
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,1);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,3);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                    } else {
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,1);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    float v = (float)xg(j,k,i,3);
                                    fwrite(&v, sizeof(float), 1, f3);
                                }
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    int iv = (int)xg(j,k,i,4);
                                    fwrite(&iv, sizeof(int), 1, f3);
                                }
                    }
                }
            } // end ibin
        } // end if (myid == myhost && icall1 == 0)


        // load appropriate data into single precision array
        jw = (j2-j1)/j3 + 1;
        kw = (k2-k1)/k3 + 1;
        iw = (i2-i1)/i3 + 1;

        if (ifunc == 0) {
            // load solution (q) in xw array
            if (ipsurf == 0) {
                iw = 0;
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    kw = 0;
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        jw = 0;
                        for (j = j1; j <= j2; j += j3) {
                            jw = jw + 1;
                            xw(jw,kw,iw,1) = q(j,k,i,1);
                            xw(jw,kw,iw,5) = q(j,k,i,5)/gm1+0.5*(q(j,k,i,2)*q(j,k,i,2)
                                           + q(j,k,i,3)*q(j,k,i,3)+q(j,k,i,4)*q(j,k,i,4))*q(j,k,i,1);
                            xw(jw,kw,iw,2) = q(j,k,i,1)*q(j,k,i,2);
                            xw(jw,kw,iw,3) = q(j,k,i,1)*q(j,k,i,3);
                            xw(jw,kw,iw,4) = q(j,k,i,1)*q(j,k,i,4);
                        }
                    }
                }
            } else {
                if (i1 == i2) {
                    if (i1 == 1) {
                        mm = 2;
                    } else {
                        mm = 4;
                    }
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xw(jw,kw,iw,1) = qi0(j,k,1,mm);
                                xw(jw,kw,iw,5) = qi0(j,k,5,mm)/gm1+0.5*(qi0(j,k,2,mm)*qi0(j,k,2,mm)
                                               + qi0(j,k,3,mm)*qi0(j,k,3,mm)
                                               + qi0(j,k,4,mm)*qi0(j,k,4,mm))*qi0(j,k,1,mm);
                                xw(jw,kw,iw,2) = qi0(j,k,1,mm)*qi0(j,k,2,mm);
                                xw(jw,kw,iw,3) = qi0(j,k,1,mm)*qi0(j,k,3,mm);
                                xw(jw,kw,iw,4) = qi0(j,k,1,mm)*qi0(j,k,4,mm);
                            }
                        }
                    }
                } else if (j1 == j2) {
                    if (j1 == 1) {
                        mm = 2;
                    } else {
                        mm = 4;
                    }
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xw(jw,kw,iw,1) = qj0(k,i,1,mm);
                                xw(jw,kw,iw,5) = qj0(k,i,5,mm)/gm1+0.5*(qj0(k,i,2,mm)*qj0(k,i,2,mm)
                                               + qj0(k,i,3,mm)*qj0(k,i,3,mm)
                                               + qj0(k,i,4,mm)*qj0(k,i,4,mm))*qj0(k,i,1,mm);
                                xw(jw,kw,iw,2) = qj0(k,i,1,mm)*qj0(k,i,2,mm);
                                xw(jw,kw,iw,3) = qj0(k,i,1,mm)*qj0(k,i,3,mm);
                                xw(jw,kw,iw,4) = qj0(k,i,1,mm)*qj0(k,i,4,mm);
                            }
                        }
                    }
                } else if (k1 == k2) {
                    if (k1 == 1) {
                        mm = 2;
                    } else {
                        mm = 4;
                    }
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xw(jw,kw,iw,1) = qk0(j,i,1,mm);
                                xw(jw,kw,iw,5) = qk0(j,i,5,mm)/gm1+0.5*(qk0(j,i,2,mm)*qk0(j,i,2,mm)
                                               + qk0(j,i,3,mm)*qk0(j,i,3,mm)
                                               + qk0(j,i,4,mm)*qk0(j,i,4,mm))*qk0(j,i,1,mm);
                                xw(jw,kw,iw,2) = qk0(j,i,1,mm)*qk0(j,i,2,mm);
                                xw(jw,kw,iw,3) = qk0(j,i,1,mm)*qk0(j,i,3,mm);
                                xw(jw,kw,iw,4) = qk0(j,i,1,mm)*qk0(j,i,4,mm);
                            }
                        }
                    }
                }
            }
        } else {
            // load desired function in xw array


            // ifunc == 3: minimum distance
            if (ifunc == 3) {
                if (ipsurf == 0) {
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xw(jw,kw,iw,1) = (double)std::abs(smin(j,k,i));
                            }
                        }
                    }
                } else {
                    // at cell face centers in a wall, smin is by definition zero
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xw(jw,kw,iw,1) = 0.0;
                            }
                        }
                    }
                }
            }
            // ifunc == 4: turbulent viscosity
            if (ifunc == 4) {
                if (ipsurf == 0) {
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xw(jw,kw,iw,1) = vist3d(j,k,i);
                            }
                        }
                    }
                } else {
                    if (i1 == i2) {
                        if (i1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = vi0(j,k,1,mm);
                                }
                            }
                        }
                    } else if (j1 == j2) {
                        if (j1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = vj0(k,i,1,mm);
                                }
                            }
                        }
                    } else if (k1 == k2) {
                        if (k1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = vk0(j,i,1,mm);
                                }
                            }
                        }
                    }
                }
            }


            // ifunc == 5: pressure coefficient (cp)
            if (ifunc == 5) {
                cpc = 2.e0/(gamma*xmach*xmach);
                if (ipsurf == 0) {
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xw(jw,kw,iw,1) = (q(j,k,i,5)/p0-1)*cpc;
                            }
                        }
                    }
                } else {
                    if (i1 == i2) {
                        if (i1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = (qi0(j,k,5,mm)/p0-1)*cpc;
                                }
                            }
                        }
                    } else if (j1 == j2) {
                        if (j1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = (qj0(k,i,5,mm)/p0-1)*cpc;
                                }
                            }
                        }
                    } else if (k1 == k2) {
                        if (k1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = (qk0(j,i,5,mm)/p0-1)*cpc;
                                }
                            }
                        }
                    }
                }
            }
            // ifunc == 6: p/pinf
            if (ifunc == 6) {
                if (ipsurf == 0) {
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xw(jw,kw,iw,1) = q(j,k,i,5)/p0;
                            }
                        }
                    }
                } else {
                    if (i1 == i2) {
                        if (i1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = qi0(j,k,5,mm)/p0;
                                }
                            }
                        }
                    } else if (j1 == j2) {
                        if (j1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = qj0(k,i,5,mm)/p0;
                                }
                            }
                        }
                    } else if (k1 == k2) {
                        if (k1 == 1) { mm = 2; } else { mm = 4; }
                        iw = 0;
                        for (i = i1; i <= i2; i += i3) {
                            iw = iw + 1;
                            kw = 0;
                            for (k = k1; k <= k2; k += k3) {
                                kw = kw + 1;
                                jw = 0;
                                for (j = j1; j <= j2; j += j3) {
                                    jw = jw + 1;
                                    xw(jw,kw,iw,1) = qk0(j,i,5,mm)/p0;
                                }
                            }
                        }
                    }
                }
            }


            // ifunc == 7: turbulence index
            if (ifunc == 7) {
                iw = 0;
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    kw = 0;
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        jw = 0;
                        for (j = j1; j <= j2; j += j3) {
                            jw = jw + 1;
                            // vorticity at wall is approximate
                            vor = std::sqrt(q(j,k,i,2)*q(j,k,i,2)+q(j,k,i,3)*q(j,k,i,3)+q(j,k,i,4)*q(j,k,i,4))/
                                  std::abs(smin(j,k,i));
                            tt = gamma*q(j,k,i,5)/q(j,k,i,1);
                            c2b = cbar/tinf;
                            c2bp = c2b + 1.0;
                            xmu = c2bp*tt*std::sqrt(tt)/(c2b+tt);
                            // iterate to get nuwiggle
                            xmut = vist3d(j,k,i);
                            xnuwiggle = 500.;
                            for (n = 1; n <= 50; n++) {
                                delta_v = (xnuwiggle*xnuwiggle*xnuwiggle*xnuwiggle
                                          - xmut*xnuwiggle*xnuwiggle*xnuwiggle
                                          - xmut*(7.1*7.1*7.1)) /
                                          (4.*xnuwiggle*xnuwiggle*xnuwiggle
                                          - 3.*xmut*xnuwiggle*xnuwiggle);
                                xnuwiggle = xnuwiggle - delta_v;
                                // iteration may not be working
                                xnuwiggle = std::max(xnuwiggle, 1.e-20);
                                xnuwiggle = std::min(xnuwiggle, 5000.);
                            }
                            xw(jw,kw,iw,1) = xnuwiggle/(std::abs(smin(j,k,i)))*
                                std::sqrt(q(j,k,i,1)*xmach/reue)/(0.41*std::sqrt(xmu*vor));
                            // limit maximum to 1 (represents turbulent)
                            if (xw(jw,kw,iw,1) > 1.0) {
                                xw(jw,kw,iw,1) = 1.0;
                            }
                        }
                    }
                }
            }

        } // end else (ifunc != 0)

        // output solution
        if (myid == myhost) {
            FILE* f4 = fortran_get_unit(4);
            if (ifunc == 0) {
                if (ibin == 0) {
                    // write(4,'(5e14.6)') xmachw,alphww,reuew,timew
                    {
                        int cnt = 0;
                        fprintf(f4, "%14.6E", (float)xmachw); ++cnt;
                        fprintf(f4, "%14.6E", (float)alphww); ++cnt;
                        fprintf(f4, "%14.6E", (float)reuew); ++cnt;
                        fprintf(f4, "%14.6E", (float)timew); ++cnt;
                        if (cnt % 5 != 0) fprintf(f4, "\n");
                    }
                    if (i2d == 0) {
                        if (ialph == 0) {
                            // write(4,'(5e14.6)') ((((xw(j,k,i,m),i=1,iw),j=1,jw),k=1,kw),m=1,5)
                            {
                                int cnt = 0;
                                for (int m = 1; m <= 5; m++)
                                    for (k = 1; k <= kw; k++)
                                        for (j = 1; j <= jw; j++)
                                            for (i = 1; i <= iw; i++) {
                                                fprintf(f4, "%14.6E", (float)xw(j,k,i,m));
                                                if (++cnt % 5 == 0) fprintf(f4, "\n");
                                            }
                                if (cnt % 5 != 0) fprintf(f4, "\n");
                            }
                        } else {
                            for (i = 1; i <= iw; i++)
                                for (j = 1; j <= jw; j++)
                                    for (k = 1; k <= kw; k++)
                                        xw(j,k,i,3) = -xw(j,k,i,3);
                            {
                                int cnt = 0;
                                for (int m = 1; m <= 2; m++)
                                    for (k = 1; k <= kw; k++)
                                        for (j = 1; j <= jw; j++)
                                            for (i = 1; i <= iw; i++) {
                                                fprintf(f4, "%14.6E", (float)xw(j,k,i,m));
                                                if (++cnt % 5 == 0) fprintf(f4, "\n");
                                            }
                                for (k = 1; k <= kw; k++)
                                    for (j = 1; j <= jw; j++)
                                        for (i = 1; i <= iw; i++) {
                                            fprintf(f4, "%14.6E", (float)xw(j,k,i,4));
                                            if (++cnt % 5 == 0) fprintf(f4, "\n");
                                        }
                                for (k = 1; k <= kw; k++)
                                    for (j = 1; j <= jw; j++)
                                        for (i = 1; i <= iw; i++) {
                                            fprintf(f4, "%14.6E", (float)xw(j,k,i,3));
                                            if (++cnt % 5 == 0) fprintf(f4, "\n");
                                        }
                                for (k = 1; k <= kw; k++)
                                    for (j = 1; j <= jw; j++)
                                        for (i = 1; i <= iw; i++) {
                                            fprintf(f4, "%14.6E", (float)xw(j,k,i,5));
                                            if (++cnt % 5 == 0) fprintf(f4, "\n");
                                        }
                                if (cnt % 5 != 0) fprintf(f4, "\n");
                            }
                        }
                    } else {
                        // i2d != 0, formatted, ifunc==0
                        {
                            int cnt = 0;
                            for (int m = 1; m <= 2; m++)
                                for (k = 1; k <= kw; k++)
                                    for (j = 1; j <= jw; j++)
                                        for (i = 1; i <= iw; i++) {
                                            fprintf(f4, "%14.6E", (float)xw(j,k,i,m));
                                            if (++cnt % 5 == 0) fprintf(f4, "\n");
                                        }
                            for (int m = 4; m <= 5; m++)
                                for (k = 1; k <= kw; k++)
                                    for (j = 1; j <= jw; j++)
                                        for (i = 1; i <= iw; i++) {
                                            fprintf(f4, "%14.6E", (float)xw(j,k,i,m));
                                            if (++cnt % 5 == 0) fprintf(f4, "\n");
                                        }
                            if (cnt % 5 != 0) fprintf(f4, "\n");
                        }
                    }
                } else {

                    // ibin != 0, ifunc == 0
                    {
                        float fv;
                        fv = (float)xmachw; fwrite(&fv, sizeof(float), 1, f4);
                        fv = (float)alphww;  fwrite(&fv, sizeof(float), 1, f4);
                        fv = (float)reuew;  fwrite(&fv, sizeof(float), 1, f4);
                        fv = (float)timew;  fwrite(&fv, sizeof(float), 1, f4);
                    }
                    if (i2d == 0) {
                        if (ialph == 0) {
                            for (int m = 1; m <= 5; m++)
                                for (k = 1; k <= kw; k++)
                                    for (j = 1; j <= jw; j++)
                                        for (i = 1; i <= iw; i++) {
                                            float fv = (float)xw(j,k,i,m);
                                            fwrite(&fv, sizeof(float), 1, f4);
                                        }
                        } else {
                            for (i = 1; i <= iw; i++)
                                for (j = 1; j <= jw; j++)
                                    for (k = 1; k <= kw; k++)
                                        xw(j,k,i,3) = -xw(j,k,i,3);
                            for (int m = 1; m <= 2; m++)
                                for (k = 1; k <= kw; k++)
                                    for (j = 1; j <= jw; j++)
                                        for (i = 1; i <= iw; i++) {
                                            float fv = (float)xw(j,k,i,m);
                                            fwrite(&fv, sizeof(float), 1, f4);
                                        }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        float fv = (float)xw(j,k,i,4);
                                        fwrite(&fv, sizeof(float), 1, f4);
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        float fv = (float)xw(j,k,i,3);
                                        fwrite(&fv, sizeof(float), 1, f4);
                                    }
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        float fv = (float)xw(j,k,i,5);
                                        fwrite(&fv, sizeof(float), 1, f4);
                                    }
                        }
                    } else {
                        // i2d != 0, binary, ifunc==0
                        for (int m = 1; m <= 2; m++)
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        float fv = (float)xw(j,k,i,m);
                                        fwrite(&fv, sizeof(float), 1, f4);
                                    }
                        for (int m = 4; m <= 5; m++)
                            for (k = 1; k <= kw; k++)
                                for (j = 1; j <= jw; j++)
                                    for (i = 1; i <= iw; i++) {
                                        float fv = (float)xw(j,k,i,m);
                                        fwrite(&fv, sizeof(float), 1, f4);
                                    }
                    }
                }
            } else {
                // ifunc != 0: function file (only 1 variable)
                if (ibin == 0) {
                    {
                        int cnt = 0;
                        for (k = 1; k <= kw; k++)
                            for (j = 1; j <= jw; j++)
                                for (i = 1; i <= iw; i++) {
                                    fprintf(f4, "%14.6E", (float)xw(j,k,i,1));
                                    if (++cnt % 5 == 0) fprintf(f4, "\n");
                                }
                        if (cnt % 5 != 0) fprintf(f4, "\n");
                    }
                } else {
                    for (k = 1; k <= kw; k++)
                        for (j = 1; j <= jw; j++)
                            for (i = 1; i <= iw; i++) {
                                float fv = (float)xw(j,k,i,1);
                                fwrite(&fv, sizeof(float), 1, f4);
                            }
                }
            }
        } // end if (myid == myhost)

    } // end if (iflag == 1)


    if (iflag == 2) {
        // print out data
        if (lhdr > 0 && myid == myhost) {
            FILE* f11 = fortran_get_unit(11);
            fprintf(f11, "writing printout file for IDIM X JDIM X KDIM =%5d x %5d x %5d grid\n",
                    idim, jdim, kdim);
        }

        // load q and turb. viscosity into single precision arrays
        jdw = (j2-j1)/j3 + 1;
        kdw = (k2-k1)/k3 + 1;
        idw = (i2-i1)/i3 + 1;

        if (ipsurf == 0) {
            iw = 0;
            for (i = i1; i <= i2; i += i3) {
                iw = iw + 1;
                jw = 0;
                for (j = j1; j <= j2; j += j3) {
                    jw = jw + 1;
                    kw = 0;
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        for (l = 1; l <= 5; l++) {
                            xw(jw,kw,iw,l) = q(j,k,i,l);
                        }
                    }
                }
            }
        } else {
            if (i1 == i2) {
                if (i1 == 1) { mm = 2; } else { mm = 4; }
                iw = 0;
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    kw = 0;
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        jw = 0;
                        for (j = j1; j <= j2; j += j3) {
                            jw = jw + 1;
                            for (l = 1; l <= 5; l++) {
                                xw(jw,kw,iw,l) = qi0(j,k,l,mm);
                            }
                        }
                    }
                }
            } else if (j1 == j2) {
                if (j1 == 1) { mm = 2; } else { mm = 4; }
                iw = 0;
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    kw = 0;
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        jw = 0;
                        for (j = j1; j <= j2; j += j3) {
                            jw = jw + 1;
                            for (l = 1; l <= 5; l++) {
                                xw(jw,kw,iw,l) = qj0(k,i,l,mm);
                            }
                        }
                    }
                }
            } else if (k1 == k2) {
                if (k1 == 1) { mm = 2; } else { mm = 4; }
                iw = 0;
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    kw = 0;
                    for (k = k1; k <= k2; k += k3) {
                        kw = kw + 1;
                        jw = 0;
                        for (j = j1; j <= j2; j += j3) {
                            jw = jw + 1;
                            for (l = 1; l <= 5; l++) {
                                xw(jw,kw,iw,l) = qk0(j,i,l,mm);
                            }
                        }
                    }
                }
            }
        }

        // load turbulent viscosity into xg(jw,kw,iw,4)
        if (cmn_reyue.ivisc[2] > 1 || cmn_reyue.ivisc[1] > 1 || cmn_reyue.ivisc[0] > 1) {
            if (ipsurf == 0) {
                iw = 0;
                for (i = i1; i <= i2; i += i3) {
                    iw = iw + 1;
                    jw = 0;
                    for (j = j1; j <= j2; j += j3) {
                        jw = jw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            xg(jw,kw,iw,4) = vist3d(j,k,i);
                        }
                    }
                }
            } else {
                if (i1 == i2) {
                    if (i1 == 1) { mm = 2; } else { mm = 4; }
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xg(jw,kw,iw,4) = vi0(j,k,1,mm);
                            }
                        }
                    }
                } else if (j1 == j2) {
                    if (j1 == 1) { mm = 2; } else { mm = 4; }
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xg(jw,kw,iw,4) = vj0(k,i,1,mm);
                            }
                        }
                    }
                } else if (k1 == k2) {
                    if (k1 == 1) { mm = 2; } else { mm = 4; }
                    iw = 0;
                    for (i = i1; i <= i2; i += i3) {
                        iw = iw + 1;
                        kw = 0;
                        for (k = k1; k <= k2; k += k3) {
                            kw = kw + 1;
                            jw = 0;
                            for (j = j1; j <= j2; j += j3) {
                                jw = jw + 1;
                                xg(jw,kw,iw,4) = vk0(j,i,1,mm);
                            }
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
                        "              cp            tur. vis.\n");
            } else {
                fprintf(f17, "\n   I   J   K         X                Y                Z"
                        "              U/Uinf           V/Vinf           W/Winf"
                        "           P/Pinf           T/Tinf            MACH"
                        "            pitotp            tur. vis.\n");
            }

            term3 = 1./( std::pow(1.+0.5*gm1*xmach*xmach, gamma/gm1) );
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
                        // turbulent viscosity
                        edvis = xg(jw,kw,iw,4);
                        // cp or pitot pressure
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
                            // write(17,29)i,j,k,xg(jw,kw,iw,1),...
                            // format(3i5,11(1x,e16.9))
                            fprintf(f17, "%5d%5d%5d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                    i, j, k,
                                    (float)xg(jw,kw,iw,1),
                                    (float)xg(jw,kw,iw,2),
                                    (float)xg(jw,kw,iw,3),
                                    (float)q2, (float)q3, (float)q4,
                                    (float)q5, (float)t1, (float)xm1,
                                    (float)pitot, (float)edvis);
                        } else {
                            xg(jw,kw,iw,2) = -xg(jw,kw,iw,2);
                            q3 = -q3;
                            fprintf(f17, "%5d%5d%5d %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                                    i, j, k,
                                    (float)xg(jw,kw,iw,1),
                                    (float)xg(jw,kw,iw,3),
                                    (float)xg(jw,kw,iw,2),
                                    (float)q2, (float)q4, (float)q3,
                                    (float)q5, (float)t1, (float)xm1,
                                    (float)pitot, (float)edvis);
                        }
                    }
                }
            }
        } // end if (myid == myhost)

    } // end if (iflag == 2)

} // end plot3c

} // namespace plot3c_ns
