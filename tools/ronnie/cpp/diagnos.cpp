// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// diagnos.cpp — faithful C++ translation of diagnos.F (CFL3D ronnie).
// Perform diagnostic checks on interpolation coefficients (generalized
// coordinates) found via the search and inversion routines.
//
// Sequential real build: real(x)->x, float(i)->f_float(i), int(x)->f_int(x)
// (truncates toward zero), abs on real -> f_abs, abs on int -> std::abs,
// ccabs(x)->ccabs (=fabs), sqrt->std::sqrt.
//
// Commons:
//   /igrdtyp/ ip3dgrd,ialph -> igrdtyp_.ip3dgrd, igrdtyp_.ialph
//   /tol/ epsc,epsc0,epsreen,epscoll -> tol_.epsc (slot0). c0tol=epsc.
//
// The C++ signature (matching invert.cpp's forward declaration) renames the
// Fortran array formals x2int,y2int,z2int,x2fit,y2fit,z2fit to
// sxie,seta,sxie2,seta2,xie2s,eta2s respectively; inside the body we #define
// the Fortran names back onto the C++ params.
//
// Array layouts (column-major, 1-based), d1=jdim1, d2=kdim1:
//   x1,y1,z1,x1mid,y1mid,z1mid,x1mide,y1mide,z1mide : (jdim1,kdim1,msub1)
//   x2int,y2int : (jdim1,kdim1,msub1) ; z2int : (jdim1,kdim1,msub2)
//   x2fit : (jdim1,kdim1,msub2) ; y2fit,z2fit : (jdim1,kdim1)  [2-D]
//   x2,y2,z2 : (jdim1,kdim1,msub2)
//   xie2,eta2,mblkpt : (npt) ; jjmax1,kkmax1 : (msub1) ; jjmax2,kkmax2 : (msub2)
//   xif1,xif2,etf1,etf2 : (msub1)
#include "ron_common.h"

// forward decls (mirror the Fortran argument lists; copied from shear.cpp)
void extra(int jdim,int kdim,int nsub,int l,
           double* x,double* y,double* z,int jcell,int kcell,int jcl,int jcr,
           double& x5,double& y5,double& z5,int& icase,int ifit);
void extrae(int jdim,int kdim,int nsub,int l,
            double* x,double* y,double* z,int jcell,int kcell,int kcl,int kcr,
            double& x7,double& y7,double& z7,int& icase,int ifit);
void direct(double x5,double x6,double x7,double x8,
            double y5,double y6,double y7,double y8,
            double z5,double z6,double z7,double z8,
            double& a1,double& a2,double& a3,int& imaxa,
            int* nou,FStr* bou,int nbuf,int ibufdim);

void diagnos(int istop,int iout,int igap,int jdim1,int kdim1,int msub1,
             int msub2,int* jjmax1,int* kkmax1,int lmax1,double* x1,double* y1,
             double* z1,double* x1mid,double* y1mid,double* z1mid,
             double* x1mide,double* y1mide,double* z1mide,double* sxie,
             double* seta,double* sxie2,double* seta2,double* xie2s,
             double* eta2s,int* jjmax2,int* kkmax2,double* x2,double* y2,
             double* z2,double* xie2,double* eta2,int* mblkpt,int icheck,
             int intmx,double xc,double yc,double zc,int ifit,int j21,int j22,
             int k21,int k22,int npt,int ic0,int iorph,int* xif1,int* xif2,
             int* etf1,int* etf2,int itoss0,int iself,int* nou,FStr* bou,
             int nbuf,int ibufdim,int myid)
{
    (void)intmx; (void)npt;
    const int d1 = jdim1, d2 = kdim1;
    (void)msub1; (void)msub2;

    // 3-D arrays (jdim1,kdim1,*)
    #define X1(i,j,k)     x1    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y1(i,j,k)     y1    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z1(i,j,k)     z1    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define X1MID(i,j,k)  x1mid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y1MID(i,j,k)  y1mid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z1MID(i,j,k)  z1mid [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define X1MIDE(i,j,k) x1mide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y1MIDE(i,j,k) y1mide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z1MIDE(i,j,k) z1mide[((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define X2(i,j,k)     x2    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y2(i,j,k)     y2    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z2(i,j,k)     z2    [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    // x2int/y2int -> sxie/seta (msub1) ; z2int -> sxie2 (msub2)
    #define X2INT(i,j,k)  sxie  [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y2INT(i,j,k)  seta  [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Z2INT(i,j,k)  sxie2 [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    // x2fit -> seta2 (msub2, 3-D) ; y2fit,z2fit -> xie2s,eta2s (2-D)
    #define X2FIT(i,j,k)  seta2 [((i)-1)+((j)-1)*d1+((k)-1)*d1*d2]
    #define Y2FIT(i,j)    xie2s [((i)-1)+((j)-1)*d1]
    #define Z2FIT(i,j)    eta2s [((i)-1)+((j)-1)*d1]
    // 1-D arrays (1-based)
    #define XIE2(i)   xie2[(i)-1]
    #define ETA2(i)   eta2[(i)-1]
    #define MBLKPT(i) mblkpt[(i)-1]
    #define JJMAX1(i) jjmax1[(i)-1]
    #define KKMAX1(i) kkmax1[(i)-1]
    #define JJMAX2(i) jjmax2[(i)-1]
    #define KKMAX2(i) kkmax2[(i)-1]
    #define XIF1(i)   xif1[(i)-1]
    #define XIF2(i)   xif2[(i)-1]
    #define ETF1(i)   etf1[(i)-1]
    #define ETF2(i)   etf2[(i)-1]

    char line[256];

    int jmax2 = JJMAX2(1) - 1;
    int kmax2 = KKMAX2(1) - 1;

    // -- FORMAT-1222..1334 literals (used repeatedly) --
    const char* F1222 = "            possible causes:";
    const char* F1223 = "              ...mmcxie and/or mmceta incorrectly set";
    const char* F1224 = "              ...\"from\" block input data incorrectly set";
    const char* F1225 = "              ...small gaps between \"from\" blocks";
    const char* F1226 = "              ...\"to\" and \"from\" blocks do not lie on the same surface";
    const char* F1227 = "              ...factj or factk too small";
    const char* F1228 = "              ...mmcxie incorrectly set";
    const char* F1229 = "              ...mmceta incorrectly set";
    const char* F1330 = "              ...interface not really C-0, check grid and/or C-0 flag";
    const char* F1331 = "              ...search routine error...uh-oh";
    const char* F1332 = "              ...these may be the orphan points anticipated by setting iorph > 0";
    const char* F1333 = "              ...branch cut on the \"from\" side may falsely trigger this message";
    const char* F1334 = "                 and the results may in fact be OK - verify with patch_p3d.xx file";

    int j,k,l,ll,l1,jc,kc;

    // diagnostics for case where search/inversion routines have failed
    if (istop == 1) {
        bou_put(nou,bou,ibufdim,4,F1222);
        bou_put(nou,bou,ibufdim,4,F1224);
        bou_put(nou,bou,ibufdim,4,F1226);
        bou_put(nou,bou,ibufdim,4,F1227);
        goto L552;
    }

    // diagnostics for completed generalized coordinate interpolation

    // list-directed write: leading space + '    beginning diagnostic checks'
    bou_put(nou,bou,ibufdim,4,"     beginning diagnostic checks");

    // ****** first check ******
    // check final cell center generalized coordinates to be sure that xie,eta
    // of the "to" grid lies inside the domain of the from grid(s)
    {
        int neta = 0;
        int nxie = 0;

        for (k = k21; k <= k22-1; ++k) {
            for (j = j21; j <= j22-1; ++j) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                l1 = MBLKPT(ll);
                if (l1 > 0) {
                    if (ETA2(ll) < 1. || ETA2(ll) > (double)(KKMAX1(l1)-2)) {
                        std::snprintf(line,sizeof(line),
                            "interp. no.: %4d, \"to\" cell center j,k = %4d%4d found in \"from\" block %4d",
                            icheck,j,k,l1);
                        bou_put(nou,bou,ibufdim,2,line);
                        std::snprintf(line,sizeof(line),"   with eta = %s",
                            tc::fmtE(ETA2(ll),14,7).c_str());
                        bou_put(nou,bou,ibufdim,2,line);
                        std::snprintf(line,sizeof(line),
                            "   in this \"from\" block, legal range of eta is 1 to %4d",
                            KKMAX1(l1)-2);
                        bou_put(nou,bou,ibufdim,2,line);
                        neta = neta + 1;
                    }
                }
            }
        }

        for (j = j21; j <= j22-1; ++j) {
            for (k = k21; k <= k22-1; ++k) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                l1 = MBLKPT(ll);
                if (l1 > 0) {
                    if (XIE2(ll) < 1. || XIE2(ll) > (double)(JJMAX1(l1)-2)) {
                        std::snprintf(line,sizeof(line),
                            "interp. no.: %4d, \"to\" cell center j,k = %4d%4d found in \"from\" block %4d",
                            icheck,j,k,l1);
                        bou_put(nou,bou,ibufdim,2,line);
                        std::snprintf(line,sizeof(line),"   with xie = %s",
                            tc::fmtE(XIE2(ll),14,7).c_str());
                        bou_put(nou,bou,ibufdim,2,line);
                        std::snprintf(line,sizeof(line),
                            "   in this \"from\" block, legal range of xie is 1 to %4d",
                            JJMAX1(l1)-2);
                        bou_put(nou,bou,ibufdim,2,line);
                        nxie = nxie + 1;
                    }
                }
            }
        }

        // format 710
        std::snprintf(line,sizeof(line),
            "     %4d points were found outside the legal range of xie",nxie);
        bou_put(nou,bou,ibufdim,4,line);
        // format 711
        std::snprintf(line,sizeof(line),
            "     %4d points were found outside the legal range of eta",neta);
        bou_put(nou,bou,ibufdim,4,line);

        if (neta+nxie > 0) {
            // list-directed
            bou_put(nou,bou,ibufdim,4,
                "         check fort.9 to make sure values are not far outside legal range");
            bou_put(nou,bou,ibufdim,4,F1222);
            if (igap > 0) {
                bou_put(nou,bou,ibufdim,4,F1225);
            } else {
                if (nxie > 0) bou_put(nou,bou,ibufdim,4,F1228);
                if (neta > 0) bou_put(nou,bou,ibufdim,4,F1229);
            }
            if (iorph > 0) {
                bou_put(nou,bou,ibufdim,4,F1332);
            }
            iout = 1;
        }
    }

    // ****** second check (grids with orphan points only) ******
    // count number of orphans
    if (iorph > 0) {
        int norph = 0;
        for (j = j21; j <= j22-1; ++j) {
            for (k = k21; k <= k22-1; ++k) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                l1 = MBLKPT(ll);
                if (l1 == 0) norph = norph + 1;
            }
        }
        // format 721
        std::snprintf(line,sizeof(line),
            "     %4d points were flagged as orphans",norph);
        bou_put(nou,bou,ibufdim,4,line);
    }

    // ****** third check (grids that self-connect with a branch cut) ******
    // check for points that get interpolated from themselves (branch cut error)
    if (iself > 0) {
        int ierr = 0;
        for (j = j21; j <= j22-1; ++j) {
            for (k = k21; k <= k22-1; ++k) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                l1 = MBLKPT(ll);
                jc = f_int(XIE2(ll));
                kc = f_int(ETA2(ll));
                if (jc == j && kc == k) {
                    ierr = ierr + 1;
                    std::snprintf(line,sizeof(line),
                        "branch cut problem at j,k =%4d%4d",j,k);
                    bou_put(nou,bou,ibufdim,2,line);
                }
            }
        }
        if (ierr > 0) {
            // format 821
            std::snprintf(line,sizeof(line),
                "     %4d points are interpolated from themselves",ierr);
            bou_put(nou,bou,ibufdim,4,line);
            // format 822
            bou_put(nou,bou,ibufdim,4,
                "          search routine failure due to branch cut - DO NOT run flow solver");
        }
    }

    // ****** fourth check (C-0 continuous grids only) ******
    // for C-0 grids, check to make sure all interpolation coefficients are
    // xie=m+0.5, eta=n+0.5, to within the tolerence set on convergence for the
    // generalized coordinates
    if (ic0 > 0) {
        double c0tol  = tol_.epsc;
        int    nc00   = (j22-j21)*(k22-k21);
        int    nc0    = nc00;
        double errmx  = 0.;
        int    jerrmx = 1;
        int    kerrmx = 1;

        for (j = j21; j <= j22-1; ++j) {
            for (k = k21; k <= k22-1; ++k) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                l1 = MBLKPT(ll); (void)l1;
                jc = f_int(XIE2(ll));
                kc = f_int(ETA2(ll));
                double xiec = ccabs(XIE2(ll)-jc-.5);
                double etac = ccabs(ETA2(ll)-kc-.5);
                if (xiec > errmx) {
                    jerrmx = j;
                    kerrmx = k;
                    errmx  = xiec;
                }
                if (etac > errmx) {
                    jerrmx = j;
                    kerrmx = k;
                    errmx  = etac;
                }
                if (xiec > c0tol || etac > c0tol) {
                    nc0 = nc0 - 1;
                }
            }
        }
        // format 722
        std::snprintf(line,sizeof(line),
            "     %4d points  flagged as C-0 have |xie-.5| or |eta-.5| >%s",
            nc00-nc0, tc::fmtE(c0tol,10,3).c_str());
        bou_put(nou,bou,ibufdim,4,line);
        // format 723
        std::snprintf(line,sizeof(line),
            "          the max. error %s occurs at j,k = %4d,%4d",
            tc::fmtE(errmx,10,3).c_str(),jerrmx,kerrmx);
        bou_put(nou,bou,ibufdim,4,line);
        if (nc00-nc0 != 0) {
            iout = 1;
            bou_put(nou,bou,ibufdim,4,F1222);
            bou_put(nou,bou,ibufdim,4,F1330);
            bou_put(nou,bou,ibufdim,4,F1224);
            bou_put(nou,bou,ibufdim,4,F1331);
        } else {
            for (j = j21; j <= j22-1; ++j) {
                for (k = k21; k <= k22-1; ++k) {
                    ll = (j22-j21)*(k-k21) + (j-j21+1);
                    jc = f_int(XIE2(ll));
                    kc = f_int(ETA2(ll));
                    XIE2(ll) = jc + 0.5;
                    ETA2(ll) = kc + 0.5;
                }
            }
            // format 720
            bou_put(nou,bou,ibufdim,4,
                "          C-0 check acceptable, resetting xie and eta to exact values (0.5,0.5)");
        }
    }

    // ****** fifth check ******
    // check jacobians of grid in generalized coordinates for anomolies
    // (negative areas). not applicable in 2D cases or cases with orphan points!
    {
        int jseg = j22 - j21;
        int kseg = k22 - k21;
        if (jseg > 1 && kseg > 1 && iorph <= 0) {

            int nct  = 0;
            int iflg = 0;
            for (l = 1; l <= lmax1; ++l) {

                // find j11,k11, the first "to" cell located in "from" block l.
                int j11 = -1;
                int k11 = -1;
                {
                    bool found = false;
                    for (j = j21; j <= j22-2 && !found; ++j) {
                        for (k = k21; k <= k22-2; ++k) {
                            ll         = (j22-j21)*(k-k21)   + (j-j21+1);
                            int ll2    = (j22-j21)*(k-k21)   + (j+1-j21+1);
                            int ll3    = (j22-j21)*(k+1-k21) + (j+1-j21+1);
                            int ll4    = (j22-j21)*(k+1-k21) + (j-j21+1);
                            if (MBLKPT(ll)==l && MBLKPT(ll2)==l &&
                                MBLKPT(ll3)==l && MBLKPT(ll4)==l) {
                                int jtest = JJMAX1(MBLKPT(ll))/2;
                                int ktest = KKMAX1(MBLKPT(ll))/2;
                                int jdif24 = std::abs(f_int(XIE2(ll2) - XIE2(ll4)));
                                int kdif24 = std::abs(f_int(ETA2(ll2) - ETA2(ll4)));
                                int jdif31 = std::abs(f_int(XIE2(ll3) - XIE2(ll)));
                                int kdif31 = std::abs(f_int(ETA2(ll3) - ETA2(ll)));
                                int ibrnch = 0;
                                if (jdif24>jtest || kdif24>ktest ||
                                    jdif31>jtest || kdif31>ktest) ibrnch = 1;
                                if (ibrnch == 0) {
                                    j11 = j;
                                    k11 = k;
                                    found = true;
                                    break;   // go to 1756
                                }
                            }
                        }
                    }
                }
                // 1756 continue
                if (j11 < 1 || k11 < 1) continue;   // go to 1750
                iflg = 1;
                j  = j11;
                k  = k11;
                ll        = (j22-j21)*(k-k21)   + (j-j21+1);
                int ll2   = (j22-j21)*(k-k21)   + (j+1-j21+1);
                int ll3   = (j22-j21)*(k+1-k21) + (j+1-j21+1);
                int ll4   = (j22-j21)*(k+1-k21) + (j-j21+1);
                double axie = XIE2(ll3) - XIE2(ll);
                double aeta = ETA2(ll3) - ETA2(ll);
                double bxie = XIE2(ll4) - XIE2(ll2);
                double beta = ETA2(ll4) - ETA2(ll2);
                double bj0  = axie*beta-bxie*aeta;

                for (j = j21; j <= j22-2; ++j) {
                    for (k = k21; k <= k22-2; ++k) {
                        ll      = (j22-j21)*(k-k21)   + (j-j21+1);
                        int m2  = (j22-j21)*(k-k21)   + (j+1-j21+1);
                        int m3  = (j22-j21)*(k+1-k21) + (j+1-j21+1);
                        int m4  = (j22-j21)*(k+1-k21) + (j-j21+1);
                        if (MBLKPT(ll)==l && MBLKPT(m2)==l &&
                            MBLKPT(m3)==l && MBLKPT(m4)==l) {
                            int jtest = JJMAX1(MBLKPT(ll))/2;
                            int ktest = KKMAX1(MBLKPT(ll))/2;
                            int jdif24 = std::abs(f_int(XIE2(m2) - XIE2(m4)));
                            int kdif24 = std::abs(f_int(ETA2(m2) - ETA2(m4)));
                            int jdif31 = std::abs(f_int(XIE2(m3) - XIE2(ll)));
                            int kdif31 = std::abs(f_int(ETA2(m3) - ETA2(ll)));
                            int ibrnch = 0;
                            if (jdif24>jtest || kdif24>ktest ||
                                jdif31>jtest || kdif31>ktest) ibrnch = 1;
                            if (ibrnch == 0) {
                                double axie2 = XIE2(m3) - XIE2(ll);
                                double aeta2 = ETA2(m3) - ETA2(ll);
                                double bxie2 = XIE2(m4) - XIE2(m2);
                                double beta2 = ETA2(m4) - ETA2(m2);
                                double bj    = axie2*beta2-bxie2*aeta2;
                                if (bj*bj0 <= 0.) {
                                    nct = nct + 1;
                                    std::snprintf(line,sizeof(line),
                                        "non-unique point at j,k = %5d,%5d",j,k);
                                    bou_put(nou,bou,ibufdim,2,line);
                                    std::snprintf(line,sizeof(line),
                                        "xie2(j,k),eta2(j,k)         = %s,%s",
                                        f_fw(XIE2(ll),11,5).c_str(),
                                        f_fw(ETA2(ll),11,5).c_str());
                                    bou_put(nou,bou,ibufdim,2,line);
                                    std::snprintf(line,sizeof(line),
                                        "xie2(j+1,k),eta2(j+1,k)     = %s,%s",
                                        f_fw(XIE2(m2),11,5).c_str(),
                                        f_fw(ETA2(m2),11,5).c_str());
                                    bou_put(nou,bou,ibufdim,2,line);
                                    std::snprintf(line,sizeof(line),
                                        "xie2(j+1,k+1),eta2(j+1,k+1) = %s,%s",
                                        f_fw(XIE2(m3),11,5).c_str(),
                                        f_fw(ETA2(m3),11,5).c_str());
                                    bou_put(nou,bou,ibufdim,2,line);
                                    std::snprintf(line,sizeof(line),
                                        "xie2(j,k+1),eta2(j,k+1)     = %s,%s",
                                        f_fw(XIE2(m4),11,5).c_str(),
                                        f_fw(ETA2(m4),11,5).c_str());
                                    bou_put(nou,bou,ibufdim,2,line);
                                }
                            }
                        }
                    }
                }
            } // 1750

            if (iflg == 0) {
                // list-directed
                bou_put(nou,bou,ibufdim,4,
                    "   caution...no consecutive \"to\" cells found in any of the \"from\" blocks");
            }
            // format 724 (note the format has a trailing space, trimmed by outbuf)
            std::snprintf(line,sizeof(line),
                "     %4d non-unique points were found in the generalized-coordinate ",nct);
            bou_put(nou,bou,ibufdim,4,line);
            // format 725
            bou_put(nou,bou,ibufdim,4,
                "       mapping between the \"to\" and \"from\" grids");
            if (nct > 0) {
                bou_put(nou,bou,ibufdim,4,
                    "         Caution: this may indicate a serious problem");
                bou_put(nou,bou,ibufdim,4,F1222);
                bou_put(nou,bou,ibufdim,4,F1223);
                bou_put(nou,bou,ibufdim,4,F1226);
                bou_put(nou,bou,ibufdim,4,F1333);
                bou_put(nou,bou,ibufdim,4,F1334);
                iout = 1;
            }
        }
    }

    // ****** sixth check ******
    // compute cell centers of "to" grid from the interpolation coefficients
    // and compare with a direct calculation from the "to" grid. both are done
    // consistant with the input value of ifit.
    {
        int    nbpt   = 0;
        double percnt = 0.10;
        double errmax = 0.;
        int    jerrmx = 1;
        int    kerrmx = 1;

        for (j = j21; j <= j22-1; ++j) {
            for (k = k21; k <= k22-1; ++k) {
                ll = (j22-j21)*(k-k21) + (j-j21+1);
                l1 = MBLKPT(ll);

                // skip if orphan point
                if (l1 == 0) continue;   // go to 9113

                jc = f_int(XIE2(ll));
                kc = f_int(ETA2(ll));

                // add 1 to account for expanded "from" grids
                jc = jc + 1;
                kc = kc + 1;
                double xiejk = XIE2(ll) + 1;
                double etajk = ETA2(ll) + 1;
                double xiec  = f_float(jc);
                double etac  = f_float(kc);

                double x11 = X1(jc,kc,l1);
                double y11 = Y1(jc,kc,l1);
                double z11 = Z1(jc,kc,l1);
                double x12 = X1(jc+1,kc,l1);
                double y12 = Y1(jc+1,kc,l1);
                double z12 = Z1(jc+1,kc,l1);
                double x13 = X1(jc+1,kc+1,l1);
                double y13 = Y1(jc+1,kc+1,l1);
                double z13 = Z1(jc+1,kc+1,l1);
                double x14 = X1(jc,kc+1,l1);
                double y14 = Y1(jc,kc+1,l1);
                double z14 = Z1(jc,kc+1,l1);
                double x15 = X1MID(jc,kc,l1);
                double y15 = Y1MID(jc,kc,l1);
                double z15 = Z1MID(jc,kc,l1);
                double x16 = X1MID(jc,kc+1,l1);
                double y16 = Y1MID(jc,kc+1,l1);
                double z16 = Z1MID(jc,kc+1,l1);
                double x17 = X1MIDE(jc,kc,l1);
                double y17 = Y1MIDE(jc,kc,l1);
                double z17 = Z1MIDE(jc,kc,l1);
                double x18 = X1MIDE(jc+1,kc,l1);
                double y18 = Y1MIDE(jc+1,kc,l1);
                double z18 = Z1MIDE(jc+1,kc,l1);

                double dx2 = x12 - x11, dy2 = y12 - y11, dz2 = z12 - z11;
                double dx3 = x13 - x11, dy3 = y13 - y11, dz3 = z13 - z11;
                double dx4 = x14 - x11, dy4 = y14 - y11, dz4 = z14 - z11;
                double dx5 = x15 - x11, dy5 = y15 - y11, dz5 = z15 - z11;
                double dx6 = x16 - x11, dy6 = y16 - y11, dz6 = z16 - z11;
                double dx7 = x17 - x11, dy7 = y17 - y11, dz7 = z17 - z11;
                double dx8 = x18 - x11, dy8 = y18 - y11, dz8 = z18 - z11;

                double a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0;
                double b2=0,b3=0,b4=0,b5=0,b6=0,b7=0,b8=0;
                double c2=0,c3=0,c4=0,c5=0,c6=0,c7=0,c8=0;
                double df1,df2,df3;

                if (ifit == 1) {
                    // bi-linear fit
                    a2 = dx2; a3 = dx4; a4 = dx3 - a2 - a3;
                    a5 = 0.; a6 = 0.; a7 = 0.; a8 = 0.;
                    b2 = dy2; b3 = dy4; b4 = dy3 - b2 - b3;
                    b5 = 0.; b6 = 0.; b7 = 0.; b8 = 0.;
                    c2 = dz2; c3 = dz4; c4 = dz3 - c2 - c3;
                    c5 = 0.; c6 = 0.; c7 = 0.; c8 = 0.;
                }

                if (ifit == 2) {
                    // (degenerate) bi-quadratic fit
                    a2  = -dx2 + 4.*dx5;
                    a3  = -dx4 + 4.*dx7;
                    a5  = 2.*dx2 - 4.*dx5;
                    a7  = 2.*dx4 - 4.*dx7;
                    df1 = dx3 - a2 - a3 - a5 - a7;
                    df2 = dx6 - .5*a2 - a3 - .25*a5 - a7;
                    df3 = dx8 - a2 - .5*a3 - a5 - .25*a7;
                    a4  = -3.*df1 + 4.*df2 + 4.*df3;
                    a6  = 2.*df1 - 4.*df2;
                    a8  = 2.*df1 - 4.*df3;
                    b2  = -dy2 + 4.*dy5;
                    b3  = -dy4 + 4.*dy7;
                    b5  = 2.*dy2 - 4.*dy5;
                    b7  = 2.*dy4 - 4.*dy7;
                    df1 = dy3 - b2 - b3 - b5 - b7;
                    df2 = dy6 - .5*b2 - b3 - .25*b5 - b7;
                    df3 = dy8 - b2 - .5*b3 - b5 - .25*b7;
                    b4  = -3.*df1 + 4.*df2 + 4.*df3;
                    b6  = 2.*df1 - 4.*df2;
                    b8  = 2.*df1 - 4.*df3;
                    c2  = -dz2 + 4.*dz5;
                    c3  = -dz4 + 4.*dz7;
                    c5  = 2.*dz2 - 4.*dz5;
                    c7  = 2.*dz4 - 4.*dz7;
                    df1 = dz3 - c2 - c3 - c5 - c7;
                    df2 = dz6 - .5*c2 - c3 - .25*c5 - c7;
                    df3 = dz8 - c2 - .5*c3 - c5 - .25*c7;
                    c4  = -3.*df1 + 4.*df2 + 4.*df3;
                    c6  = 2.*df1 - 4.*df2;
                    c8  = 2.*df1 - 4.*df3;
                }

                if (ifit == 3) {
                    // quadratic fit in xie, linear fit in eta
                    a3  = dx4; b3  = dy4; c3  = dz4;
                    a2  = -dx2 + 4.*dx5;
                    b2  = -dy2 + 4.*dy5;
                    c2  = -dz2 + 4.*dz5;
                    a5  = 2.*dx2  - 4.*dx5;
                    b5  = 2.*dy2  - 4.*dy5;
                    c5  = 2.*dz2  - 4.*dz5;
                    df1 = x13 - x12 - a3;
                    df2 = x16 - x15 - a3;
                    a4  =   -df1 + 4.*df2;
                    a6  = 2.*df1 - 4.*df2;
                    df1 = y13 - y12 - b3;
                    df2 = y16 - y15 - b3;
                    b4  =   -df1 + 4.*df2;
                    b6  = 2.*df1 - 4.*df2;
                    df1 = z13 - z12 - c3;
                    df2 = z16 - z15 - c3;
                    c4  =   -df1 + 4.*df2;
                    c6  = 2.*df1 - 4.*df2;
                    a7  = 0.; a8  = 0.;
                    b7  = 0.; b8  = 0.;
                    c7  = 0.; c8  = 0.;
                }

                if (ifit == 4) {
                    // linear fit in xie, quadratic fit in eta
                    a2  = dx2;
                    a3  = -dx4 + 4.*dx7;
                    a7  = 2.*dx4 -4.*dx7;
                    df1 = dx3 - a2 - a3 - a7;
                    df2 = dx8 - a2 -.5*a3 - .25*a7;
                    a4  = -df1 + 4.*df2;
                    a8  = 2.*df1 - 4.*df2;
                    a5  = 0.; a6  = 0.;
                    b2  = dy2;
                    b3  = -dy4 + 4.*dy7;
                    b7  = 2.*dy4 -4.*dy7;
                    df1 = dy3 - b2 - b3 - b7;
                    df2 = dy8 - b2 -.5*b3 - .25*b7;
                    b4  = -df1 + 4.*df2;
                    b8  = 2.*df1 - 4.*df2;
                    b5  = 0.; b6  = 0.;
                    c2  = dz2;
                    c3  = -dz4 + 4.*dz7;
                    c7  = 2.*dz4 -4.*dz7;
                    df1 = dz3 - c2 - c3 - c7;
                    df2 = dz8 - c2 -.5*c3 - .25*c7;
                    c4  = -df1 + 4.*df2;
                    c8  = 2.*df1 - 4.*df2;
                    c5  = 0.; c6  = 0.;
                }

                double xie = xiejk - xiec;
                double eta = etajk - etac;

                X2INT(j,k,1) = x11 + a3*eta + eta*( a7*eta + a8*xie*eta )
                    + xie*( a2 + a4*eta + a5*xie + a6*xie*eta );
                Y2INT(j,k,1) = y11 + b3*eta + eta*( b7*eta + b8*xie*eta )
                    + xie*( b2 + b4*eta + b5*xie + b6*xie*eta );
                Z2INT(j,k,1) = z11 + c3*eta + eta*( c7*eta + c8*xie*eta )
                    + xie*( c2 + c4*eta + c5*xie + c6*xie*eta );

                // compute center of "to" cell directly from "to" grid, consistent
                // with ifit
                int jl = 1;
                int jr = jmax2;
                int kl = 1;
                int kr = kmax2;
                int kcall = k+1;
                int icase;
                double x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;
                extra(jdim1,kdim1,msub2,1,x2,y2,z2,j,kcall,jl,jr,x6,y6,z6,icase,ifit);
                extra(jdim1,kdim1,msub2,1,x2,y2,z2,j,k,jl,jr,x5,y5,z5,icase,ifit);
                extrae(jdim1,kdim1,msub2,1,x2,y2,z2,j,k,kl,kr,x7,y7,z7,icase,ifit);
                int jcall = j+1;
                extrae(jdim1,kdim1,msub2,1,x2,y2,z2,jcall,k,kl,kr,x8,y8,z8,icase,ifit);

                double x2c=0,y2c=0,z2c=0;
                if (ifit == 1) {
                    // bi-linear
                    x2c = 0.25*( X2(j,k,1) + X2(j+1,k,1)
                          + X2(j+1,k+1,1) + X2(j,k+1,1) );
                    y2c = 0.25*( Y2(j,k,1) + Y2(j+1,k,1)
                          + Y2(j+1,k+1,1) + Y2(j,k+1,1) );
                    z2c = 0.25*( Z2(j,k,1) + Z2(j+1,k,1)
                          + Z2(j+1,k+1,1) + Z2(j,k+1,1) );
                }
                if (ifit == 2) {
                    // bi-quadratic
                    x2c = 0.5* ( x5 + x6 + x7 + x8 )
                        -0.25*( X2(j,k,1)     + X2(j+1,k,1)
                        +       X2(j+1,k+1,1) + X2(j,k+1,1) );
                    y2c = 0.5* ( y5 + y6 + y7 + y8 )
                        -0.25*( Y2(j,k,1)    + Y2(j+1,k,1)
                        +      Y2(j+1,k+1,1) + Y2(j,k+1,1) );
                    z2c = 0.5* ( z5 + z6 + z7 + z8 )
                        -0.25*( Z2(j,k,1)     + Z2(j+1,k,1)
                        +       Z2(j+1,k+1,1) + Z2(j,k+1,1) );
                }
                if (ifit == 3) {
                    // quadratic in xie, linear in eta
                    x2c = .5*(x5 + x6);
                    y2c = .5*(y5 + y6);
                    z2c = .5*(z5 + z6);
                }
                if (ifit == 4) {
                    // linear in xie, quadratic in eta
                    x2c = .5*(x7 + x8);
                    y2c = .5*(y7 + y8);
                    z2c = .5*(z7 + z8);
                }

                X2FIT(j,k,1) = x2c;
                Y2FIT(j,k)   = y2c;
                Z2FIT(j,k)   = z2c;

                // compare "to" cell center from generalized-coord interpolation
                // with "to" cell center from applying ifit to grid points.
                double x00 = X2(j,k,1);
                double y00 = Y2(j,k,1);
                double z00 = Z2(j,k,1);
                double xscal = 0.;
                double yscal = 0.;
                double zscal = 0.;
                for (int jj = 1; jj <= 2; ++jj) {
                    for (int kk = 1; kk <= 2; ++kk) {
                        double xsc = ccabs(X2(j+jj-1,k+kk-1,1)-x00);
                        if (xsc > xscal) xscal = xsc;
                        double ysc = ccabs(Y2(j+jj-1,k+kk-1,1)-y00);
                        if (ysc > yscal) yscal = ysc;
                        double zsc = ccabs(Z2(j+jj-1,k+kk-1,1)-z00);
                        if (zsc > zscal) zscal = zsc;
                    }
                }
                if (xscal == 0) xscal = -1.;
                if (yscal == 0) yscal = -1.;
                if (zscal == 0) zscal = -1.;

                // check errors only in transverse directions - not projected
                double xerr = 0;
                double yerr = 0.;
                double zerr = 0.;
                int itoss;
                if (itoss0 == 0) {
                    double a1,ta2,ta3;
                    direct(x15,x16,x17,x18,y15,y16,y17,y18,z15,z16,z17,z18,
                           a1,ta2,ta3,itoss,nou,bou,nbuf,ibufdim);
                } else {
                    itoss = itoss0;
                }
                if (itoss == 1) {
                    yerr = ccabs(Y2INT(j,k,1) - y2c)/yscal;
                    zerr = ccabs(Z2INT(j,k,1) - z2c)/zscal;
                }
                if (itoss == 2) {
                    xerr = ccabs(X2INT(j,k,1) - x2c)/xscal;
                    zerr = ccabs(Z2INT(j,k,1) - z2c)/zscal;
                }
                if (itoss == 3) {
                    xerr = ccabs(X2INT(j,k,1) - x2c)/xscal;
                    yerr = ccabs(Y2INT(j,k,1) - y2c)/yscal;
                }
                double error = xerr;
                if (yerr > error) error = yerr;
                if (zerr > error) error = zerr;
                if (error > errmax) {
                    errmax = error;
                    jerrmx = j;
                    kerrmx = k;
                }
                if (xerr > percnt || yerr > percnt || zerr > percnt) {
                    nbpt = nbpt + 1;

                    std::snprintf(line,sizeof(line),
                        " interpolation no. %5d for to cell j,k %5d,%5d",icheck,j,k);
                    bou_put(nou,bou,ibufdim,2,line);
                    std::snprintf(line,sizeof(line),
                        "    xint,yint,zint= %s,%s,%s",
                        f_fw(X2INT(j,k,1),11,5).c_str(),
                        f_fw(Y2INT(j,k,1),11,5).c_str(),
                        f_fw(Z2INT(j,k,1),11,5).c_str());
                    bou_put(nou,bou,ibufdim,2,line);
                    std::snprintf(line,sizeof(line),
                        "    x,y,z= %s,%s,%s",
                        f_fw(x2c,11,5).c_str(),
                        f_fw(y2c,11,5).c_str(),
                        f_fw(z2c,11,5).c_str());
                    bou_put(nou,bou,ibufdim,2,line);
                    std::snprintf(line,sizeof(line),
                        "    xscal,yscal,zscal= %s,%s,%s",
                        f_fw(xscal,11,5).c_str(),
                        f_fw(yscal,11,5).c_str(),
                        f_fw(zscal,11,5).c_str());
                    bou_put(nou,bou,ibufdim,2,line);
                }
            } // 9113 (k)
        } // 9113 (j)

        // format 726
        std::snprintf(line,sizeof(line),
            "     %4d interpolated cell centers differ by more than %3d percent",
            nbpt, f_int(percnt*100));
        bou_put(nou,bou,ibufdim,4,line);
        // format 7261
        bou_put(nou,bou,ibufdim,4,
            "          from those obtained directly from the grid points");
        // format 727
        std::snprintf(line,sizeof(line),
            "          the max. difference %s percent occurs at j,k = %4d,%4d",
            tc::fmtE(errmax*100.,10,3).c_str(),jerrmx,kerrmx);
        bou_put(nou,bou,ibufdim,4,line);
        if (nbpt > 0) {
            bou_put(nou,bou,ibufdim,4,
                "           note: differences between interpolated and directly obtained cell");
            bou_put(nou,bou,ibufdim,4,
                "           centers often (correctly) arise when boundaries are being rendered");
            bou_put(nou,bou,ibufdim,4,
                "           coincident. however, if differences occur when boundaries are");
            bou_put(nou,bou,ibufdim,4,
                "           not rendered coincident, or if the other diagnostic counts are");
            bou_put(nou,bou,ibufdim,4,
                "           non-zero, a problem in the input or the grid is likely");
        }
    }

    // ****** seventh check ******
    // check normals on (cell center) grids as calculated 1) directly from the
    // "to" grid, and 2) via generalized coordinate interpolation.
    // not applicable in 2D cases or for cases with orphans!
    {
        int jseg = j22 - j21;
        int kseg = k22 - k21;
        if (jseg > 1 && kseg > 1 && iorph <= 0) {

            double epsn   = .1;
            int    nnorm  = 0;
            double prdsum = 0.;
            double dnmx   = 0.;
            int    jdnmx  = 1;
            int    kdnmx  = 1;

            for (j = j21; j <= j22-2; ++j) {
                for (k = k21; k <= k22-2; ++k) {

                    // normals at "to" cell centers directly from "to" grid points
                    double ax = X2FIT(j+1,k+1,1) - X2FIT(j,k,1);
                    double ay = Y2FIT(j+1,k+1) - Y2FIT(j,k);
                    double az = Z2FIT(j+1,k+1) - Z2FIT(j,k);
                    double bx = X2FIT(j,k+1,1) - X2FIT(j+1,k,1);
                    double by = Y2FIT(j,k+1) - Y2FIT(j+1,k);
                    double bz = Z2FIT(j,k+1) - Z2FIT(j+1,k);
                    double dnx = ay*bz - az*by;
                    double dny = az*bx - ax*bz;
                    double dnz = ax*by - ay*bx;
                    double d = std::sqrt(dnx*dnx +dny*dny +dnz*dnz);
                    if (d <= 0.) d = 1.;
                    double xnfit = dnx/d;
                    double ynfit = dny/d;
                    double znfit = dnz/d;

                    // normals at "to" cell centers from "to" cell centers
                    ax = X2INT(j+1,k,1) - X2INT(j,k,1);
                    ay = Y2INT(j+1,k,1) - Y2INT(j,k,1);
                    az = Z2INT(j+1,k,1) - Z2INT(j,k,1);
                    bx = X2INT(j,k+1,1) - X2INT(j,k,1);
                    by = Y2INT(j,k+1,1) - Y2INT(j,k,1);
                    bz = Z2INT(j,k+1,1) - Z2INT(j,k,1);
                    dnx = ay*bz - az*by;
                    dny = az*bx - ax*bz;
                    dnz = ax*by - ay*bx;
                    d = std::sqrt(dnx*dnx +dny*dny +dnz*dnz);
                    if (d <= 0) d = 1.;
                    double xnint = dnx/d;
                    double ynint = dny/d;
                    double znint = dnz/d;

                    // inner product of the two unit normals
                    double prod = xnfit*xnint + ynfit*ynint + znfit*znint;
                    prdsum = prdsum + prod;
                    if (prod < 1.-epsn || prod > 1.+epsn) {
                        nnorm = nnorm + 1;
                    }
                    double dn = prod - 1.;
                    if (f_abs(dn) > f_abs(dnmx)) {
                        dnmx  = dn;
                        jdnmx = j;
                        kdnmx = k;
                    }
                }
            }

            // format 728
            std::snprintf(line,sizeof(line),
                "     %4d normals at interpolated cell centers differ by more than %3d percent",
                nnorm, f_int(epsn*100));
            bou_put(nou,bou,ibufdim,4,line);
            // list-directed
            bou_put(nou,bou,ibufdim,4,
                "         from the normals at cell centers obtained directly from grid points");
            // format 729
            std::snprintf(line,sizeof(line),
                "       the max. difference %s percent occurs at j,k = %4d,%4d",
                tc::fmtE(dnmx*100.,10,3).c_str(),jdnmx,kdnmx);
            bou_put(nou,bou,ibufdim,4,line);
            // format 730
            std::snprintf(line,sizeof(line),
                "       the average inner product over the interface is %s",
                tc::fmtE(prdsum/(double)(j22-j21-1)/(double)(k22-k21-1),10,3).c_str());
            bou_put(nou,bou,ibufdim,4,line);

            if (nnorm > 0) {
                bou_put(nou,bou,ibufdim,4,F1222);
                bou_put(nou,bou,ibufdim,4,F1223);
                bou_put(nou,bou,ibufdim,4,F1226);
            }
        }
    }

    // ****** eighth check ******
    // Following output for visual check of the interpolation procedure.
L552:
    if (iout > 0) {
        int iunit = 40;
        (void)iunit;

        char titl[16];
        int len1;
        if (icheck > 99) {
            len1 = 13;
            std::snprintf(titl,sizeof(titl),"patch_p3d.%3d",icheck);
        } else if (icheck > 9) {
            len1 = 12;
            std::snprintf(titl,sizeof(titl),"patch_p3d.%2d",icheck);
        } else {
            len1 = 11;
            std::snprintf(titl,sizeof(titl),"patch_p3d.%1d",icheck);
        }
        // pad to 14 chars (titlptchgrd is character*14)
        for (int i = len1; i < 14; ++i) titl[i] = ' ';
        titl[14] = '\0';
        std::string fname(titl, len1);

        FILE* f = std::fopen(fname.c_str(), "w");
        if (f) {
            if (iorph <= 0) {
                std::snprintf(line,sizeof(line),
                    "   plot3d (/mg/for) file for this patch interface written to%s",
                    titl);   // a14
                bou_put(nou,bou,ibufdim,4,line);
            } else {
                std::snprintf(line,sizeof(line),
                    "   plot3d (/mg/for/blank) file for this patch interface written to%s",
                    titl);
                bou_put(nou,bou,ibufdim,4,line);
            }

            // helper for list-directed real
            auto gl = [&](double v){ std::fprintf(f,"%s", f_glist(v).c_str()); };

            if (iorph <= 0) std::fprintf(f,"%s\n", f_ild(lmax1+2).c_str());
            else            std::fprintf(f,"%s\n", f_ild(lmax1+3).c_str());

            int iwhole = 0;
            int maxx   = 0;
            if (iwhole == 0) {
                for (l = 1; l <= lmax1; ++l) {
                    if (XIF1(l) == 1)         XIF1(l) = XIF1(l) + 1;
                    if (ETF1(l) == 1)         ETF1(l) = ETF1(l) + 1;
                    if (XIF2(l) == JJMAX1(l)) XIF2(l) = XIF2(l) - 1;
                    if (ETF2(l) == KKMAX1(l)) ETF2(l) = ETF2(l) - 1;
                }
            }
            if (istop != 1) {
                if (iorph <= 0) {
                    for (l = 1; l <= lmax1; ++l) {
                        std::fprintf(f,"%s%s%s",
                            f_ild(XIF2(l)-XIF1(l)+1-maxx).c_str(),
                            f_ild(ETF2(l)-ETF1(l)+1-maxx).c_str(),
                            f_ild(1).c_str());
                    }
                    std::fprintf(f,"%s%s%s%s%s%s\n",
                        f_ild(j22-j21+1).c_str(),f_ild(k22-k21+1).c_str(),f_ild(1).c_str(),
                        f_ild(j22-j21).c_str(),f_ild(k22-k21).c_str(),f_ild(1).c_str());
                } else {
                    for (l = 1; l <= lmax1; ++l) {
                        std::fprintf(f,"%s%s%s",
                            f_ild(XIF2(l)-XIF1(l)+1-maxx).c_str(),
                            f_ild(ETF2(l)-ETF1(l)+1-maxx).c_str(),
                            f_ild(1).c_str());
                    }
                    std::fprintf(f,"%s%s%s%s%s%s%s%s%s\n",
                        f_ild(j22-j21+1).c_str(),f_ild(k22-k21+1).c_str(),f_ild(1).c_str(),
                        f_ild(j22-j21).c_str(),f_ild(k22-k21).c_str(),f_ild(1).c_str(),
                        f_ild(j22-j21).c_str(),f_ild(k22-k21).c_str(),f_ild(1).c_str());
                }
            } else {
                for (l = 1; l <= lmax1; ++l) {
                    std::fprintf(f,"%s%s%s",
                        f_ild(XIF2(l)-XIF1(l)+1-maxx).c_str(),
                        f_ild(ETF2(l)-ETF1(l)+1-maxx).c_str(),
                        f_ild(1).c_str());
                }
                std::fprintf(f,"%s%s%s%s%s%s\n",
                    f_ild(j22-j21+1).c_str(),f_ild(k22-k21+1).c_str(),f_ild(1).c_str(),
                    f_ild(1).c_str(),f_ild(1).c_str(),f_ild(1).c_str());
            }

            for (l = 1; l <= lmax1; ++l) {
                int js = XIF1(l);
                int ks = ETF1(l);
                int je = XIF2(l);
                int ke = ETF2(l);
                if (iorph <= 0) {
                    if (igrdtyp_.ialph == 0) {
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(X1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(Y1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(Z1(j,k,l));
                    } else {
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(X1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(Z1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(-Y1(j,k,l));
                    }
                } else {
                    if (igrdtyp_.ialph == 0) {
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(X1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(Y1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(Z1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) std::fprintf(f,"%s",f_ild(1).c_str());
                    } else {
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(X1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(Z1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) gl(-Y1(j,k,l));
                        for (k=ks;k<=ke;++k) for (j=js;j<=je;++j) std::fprintf(f,"%s",f_ild(1).c_str());
                    }
                }
                std::fprintf(f,"\n");
            }
            // "to" grid
            if (iorph <= 0) {
                if (igrdtyp_.ialph == 0) {
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(X2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(Y2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(Z2(j,k,1));
                } else {
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(X2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(Z2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(-Y2(j,k,1));
                }
            } else {
                if (igrdtyp_.ialph == 0) {
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(X2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(Y2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(Z2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) std::fprintf(f,"%s",f_ild(1).c_str());
                } else {
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(X2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(Z2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) gl(-Y2(j,k,1));
                    for (k=k21;k<=k22;++k) for (j=j21;j<=j22;++j) std::fprintf(f,"%s",f_ild(1).c_str());
                }
            }
            std::fprintf(f,"\n");

            if (istop == 1) {
                if (igrdtyp_.ialph == 0) {
                    gl(xc); gl(yc); gl(zc);
                } else {
                    gl(xc); gl(zc); gl(-yc);
                }
                std::fprintf(f,"\n");
            } else {
                if (iorph <= 0) {
                    if (igrdtyp_.ialph == 0) {
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(X2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Y2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Z2INT(j,k,1));
                    } else {
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(X2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Z2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(-Y2INT(j,k,1));
                    }
                    std::fprintf(f,"\n");
                } else {
                    for (j=j21;j<=j22-1;++j) for (k=k21;k<=k22-1;++k) {
                        ll = (j22-j21)*(k-k21) + (j-j21+1);
                        l1 = MBLKPT(ll);
                        X2FIT(j,k,1) = 1.;
                        if (l1 == 0) X2FIT(j,k,1) = 0.;
                    }
                    if (igrdtyp_.ialph == 0) {
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(X2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Y2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Z2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) std::fprintf(f,"%s",f_ild(f_int(X2FIT(j,k,1))).c_str());
                    } else {
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(X2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Z2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(-Y2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) std::fprintf(f,"%s",f_ild(f_int(X2FIT(j,k,1))).c_str());
                    }
                    std::fprintf(f,"\n");
                    for (j=j21;j<=j22-1;++j) for (k=k21;k<=k22-1;++k) {
                        ll = (j22-j21)*(k-k21) + (j-j21+1);
                        l1 = MBLKPT(ll);
                        X2FIT(j,k,1) = 0.;
                        if (l1 == 0) {
                            X2FIT(j,k,1) = 1.;
                            X2INT(j,k,1) = 0.25*( X2(j,k,1) + X2(j+1,k,1)
                                 + X2(j+1,k+1,1) + X2(j,k+1,1) );
                            Y2INT(j,k,1) = 0.25*( Y2(j,k,1) + Y2(j+1,k,1)
                                 + Y2(j+1,k+1,1) + Y2(j,k+1,1) );
                            Z2INT(j,k,1) = 0.25*( Z2(j,k,1) + Z2(j+1,k,1)
                                 + Z2(j+1,k+1,1) + Z2(j,k+1,1) );
                        }
                    }
                    if (igrdtyp_.ialph == 0) {
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(X2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Y2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Z2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) std::fprintf(f,"%s",f_ild(f_int(X2FIT(j,k,1))).c_str());
                    } else {
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(X2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(Z2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) gl(-Y2INT(j,k,1));
                        for (k=k21;k<=k22-1;++k) for (j=j21;j<=j22-1;++j) std::fprintf(f,"%s",f_ild(f_int(X2FIT(j,k,1))).c_str());
                    }
                    std::fprintf(f,"\n");
                }
            }
            std::fclose(f);
        }
        if (istop > 0) termn8(myid,-1,ibufdim,nbuf,bou,nou);
    }

    (void)xc; (void)yc; (void)zc;
    return;

    #undef X1
    #undef Y1
    #undef Z1
    #undef X1MID
    #undef Y1MID
    #undef Z1MID
    #undef X1MIDE
    #undef Y1MIDE
    #undef Z1MIDE
    #undef X2
    #undef Y2
    #undef Z2
    #undef X2INT
    #undef Y2INT
    #undef Z2INT
    #undef X2FIT
    #undef Y2FIT
    #undef Z2FIT
    #undef XIE2
    #undef ETA2
    #undef MBLKPT
    #undef JJMAX1
    #undef KKMAX1
    #undef JJMAX2
    #undef KKMAX2
    #undef XIF1
    #undef XIF2
    #undef ETF1
    #undef ETF2
}
