// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "deform.h"
#include "runtime/fortran_io.h"
#include "umalloc.h"
#include "arclen.h"
#include "tfiedge.h"
#include "tfiface.h"
#include "tfivol.h"
#include "bc_delt.h"
#include "updatedg.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>


namespace deform_ns {

void deform_el(FortranArray3DRef<int> islavept, int& nslave, int& nmaster,
               FortranArray1DRef<double> u,
               FortranArray1DRef<double> xst, FortranArray1DRef<double> yst, FortranArray1DRef<double> zst,
               int& nt, int& myhost, int& mycomm, int& myid, int& nnodes,
               FortranArray1DRef<int> mblk2nd, FortranArray2DRef<int> nblelst,
               int& maxbl, int& iseqr)
{
    // COMMON blocks
    int& i2d     = cmn_twod.i2d;
    int& iexp    = cmn_zero.iexp;
    int& meshdef = cmn_deformz.meshdef;
    int& isktyp  = cmn_deformz.isktyp;
    int& nsprgit = cmn_deformz.nsprgit;
    int& negvol  = cmn_deformz.negvol;
    int& ndgrd   = cmn_deformz.ndgrd;
    int& ndwrt   = cmn_deformz.ndwrt;
    int& ntstep  = cmn_unst.ntstep;

    int stats = 0;
    int memuse = 0;

    // Allocatable arrays
    FortranArray1D<double> sa;
    FortranArray2D<double> stiffl;
    FortranArray1D<int>    ija;
    FortranArray1D<double> b;
    FortranArray2D<double> ripm;
    FortranArray1D<double> ei, ej, ek, gij, gjk, gik;
    FortranArray1D<double> xix, xiy, xiz;
    FortranArray1D<double> etax, etay, etaz;
    FortranArray1D<double> zetax, zetay, zetaz;
    FortranArray1D<double> ooj;
    FortranArray1D<double> volij;
    FortranArray1D<double> volik;
    FortranArray2D<double> tinv, cjm1, cjp1, cim1, cip1, ckm1, ckp1;
    FortranArray2D<double> c, t, ul;
    FortranArray1D<int>    n1;

    sa.allocate(245*nslave+2);
    { int nw=245*nslave+2, iz=0; char txt[1]={'s'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    stiffl.allocate(24,24);
    { int nw=576, iz=0; char txt[1]={'s'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ija.allocate(245*nslave+2);
    { int nw=245*nslave+2, iz=1; char txt[1]={'i'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    b.allocate(3*nslave);
    { int nw=3*nslave, iz=0; char txt[1]={'b'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ripm.allocate(7,nslave);
    { int nw=7*nslave, iz=0; char txt[1]={'r'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ei.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'e'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ej.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'e'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ek.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'e'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    gjk.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'g'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    gij.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'g'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    gik.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'g'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    xix.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'x'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    xiy.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'x'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    xiz.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'x'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    etax.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'e'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    etay.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'e'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    etaz.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'e'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    zetax.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'z'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    zetay.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'z'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    zetaz.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'z'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ooj.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'o'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    volij.allocate(3*nslave);
    { int nw=3*nslave, iz=0; char txt[1]={'v'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    volik.allocate(nslave);
    { int nw=nslave, iz=0; char txt[1]={'v'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    tinv.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'t'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    cjm1.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'c'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    cjp1.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'c'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    cim1.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'c'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    cip1.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'c'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ckm1.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'c'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ckp1.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'c'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    c.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'c'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    t.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'t'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ul.allocate(6,6);
    { int nw=36, iz=0; char txt[1]={'u'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    n1.allocate(20);
    { int nw=20, iz=1; char txt[1]={'n'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }

    double eps = 100.0 * std::pow(10.0, (double)(-iexp));

    // Initialize ija
    int ipt = 3*nslave+2;
    ija(1) = ipt;
    for (int n = 1; n <= nslave; n++) {
        ija(3*(n-1)+2) = ija(3*(n-1)+1) + islavept(n,10,iseqr);
        ija(3*(n-1)+3) = ija(3*(n-1)+2) + islavept(n,10,iseqr);
        ija(3*n    +1) = ija(3*(n-1)+3) + islavept(n,10,iseqr);
    }
    for (int n = 3*nslave+2; n <= 245*nslave+2; n++) {
        ija(n) = 0;
    }

    if (nt == 1 && meshdef == 1) {
        for (int n = 1; n <= nslave; n++) {
            int iimax = islavept(n,11,iseqr);
            if (iimax > 1) {
                for (int ii2 = 2; ii2 <= iimax; ii2++) {
                    int n2 = islavept(n,12+ii2-2,iseqr);
                    double test = (xst(n2)-xst(n))*(xst(n2)-xst(n))
                                 +(yst(n2)-yst(n))*(yst(n2)-yst(n))
                                 +(zst(n2)-zst(n))*(zst(n2)-zst(n));
                    test = std::sqrt(test);
                    if (test > 100.0*eps) {
                        FILE* f11 = fortran_get_unit(11);
                        fprintf(f11, " WARNING: Macro-elements %6d and %6d"
                                " at 1-1 blocking interface have geometric mismatch = %16.8e\n",
                                n, n2, test);
                    }
                }
            }
        }
    }

    double gini = 0.0, eini = 0.0, arg1 = 0.0, arg2 = 0.0;
    int ndim = 0;

    if (isktyp > 0) {
        gini = 10.0;
        eini = 10.0;
        int mxiter1 = 1;

        FortranArray1DRef<double> b_ref = b.ref();
        elrhs(b_ref, u, islavept, nslave, nmaster, iseqr);

        for (int itr1 = 1; itr1 <= mxiter1; itr1++) {
            for (int n = 1; n <= 245*nslave+2; n++) {
                sa(n) = 0.0;
            }

            if (meshdef == 1) {
                if (myid == myhost) {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " Calculating macroelement moduli for %8d element nodes\n", nslave);
                }
            }

            hookefe(ei.ref(), ej.ref(), ek.ref(), gij.ref(), gjk.ref(), gik.ref(),
                    xst, yst, zst, volij.ref(), volik.ref(), eps, eini, gini,
                    arg1, arg2, islavept, nslave, nmaster, nt, nnodes,
                    myhost, myid, mycomm, mblk2nd, maxbl, iseqr);

            if (meshdef == 1) {
                if (myid == myhost) {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " Calculating macroelement metrics \n");
                }
            }

            elmetricsfe(xix.ref(), xiy.ref(), xiz.ref(),
                        etax.ref(), etay.ref(), etaz.ref(),
                        zetax.ref(), zetay.ref(), zetaz.ref(),
                        ooj.ref(), xst, yst, zst, eps,
                        islavept, nslave, nmaster, nnodes,
                        myhost, myid, mycomm, mblk2nd, maxbl, iseqr);

            if (meshdef == 1) {
                if (myid == myhost) {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " Calculating and assembling macroelement stiffness  matrix \n");
                }
            }

            updatedg_ns::elglobfe(sa.ref(), xst, yst, zst,
                     xix.ref(), xiy.ref(), xiz.ref(),
                     etax.ref(), etay.ref(), etaz.ref(),
                     zetax.ref(), zetay.ref(), zetaz.ref(),
                     ei.ref(), ej.ref(), ek.ref(),
                     gij.ref(), gjk.ref(), gik.ref(),
                     ooj.ref(), eps, stiffl.ref(),
                     islavept, ija.ref(),
                     nslave, nmaster, nnodes,
                     myhost, myid, mycomm,
                     mblk2nd, nblelst, maxbl, iseqr);

            if (meshdef == 1) {
                if (myid == myhost) {
                    FILE* f11 = fortran_get_unit(11);
                    fprintf(f11, " Solving the macroelement system \n");
                }
            }

            int ncount = 245*nslave+2;
            ndim = 3*nslave;
            double tol = 1.e-6;
            int itmax_gs = 9*nslave;
            double err = 0.0;
            gaussseidel(ncount, ndim, b_ref, u, sa.ref(), ija.ref(),
                        tol, itmax_gs, err, volij.ref(), myid, myhost);

            coincdef(volij.ref(), u, islavept, nslave, nmaster, iseqr, n1.ref());
        }
    } else {
        if (meshdef == 1) {
            if (myid == myhost) {
                FILE* f11 = fortran_get_unit(11);
                fprintf(f11, " Solving nodal displacements using exponential decay  \n");
            }
        }

        expdecay(xst, yst, zst, u, volij.ref(), ei.ref(), ripm.ref(), eps,
                 islavept, nslave, nmaster, ndim, nt,
                 nnodes, myhost, myid, mycomm, mblk2nd, maxbl, iseqr);
        coincdef(volij.ref(), u, islavept, nslave, nmaster, iseqr, n1.ref());
    }

    if (i2d != 0) {
        for (int n = 1; n <= nslave; n++) {
            if (islavept(n,6,iseqr) == n) {
                int nip1 = islavept(n,7,iseqr);
                u(3*(n-1)+1)    = 0.5*(u(3*(n-1)+1) + u(3*(nip1-1)+1));
                u(3*(nip1-1)+1) = u(3*(n-1)+1);
                u(3*(nip1-1)+2) = 0.0;
                u(3*(n-1)+2)    = 0.0;
                u(3*(n-1)+3)    = 0.5*(u(3*(n-1)+3) + u(3*(nip1-1)+3));
                u(3*(nip1-1)+3) = u(3*(n-1)+3);
            }
        }
    }

    // release memory (automatic via FortranArray destructors)
}


void coincdef(FortranArray1DRef<double> volij, FortranArray1DRef<double> u,
              FortranArray3DRef<int> islavept, int& nslave, int& nmaster,
              int& iseqr, FortranArray1DRef<int> n1)
{
    // volij = 0.
    for (int i = volij.lbound(1); i <= volij.ubound(1); i++) volij(i) = 0.0;

    for (int n = 1; n <= nslave; n++) {
        int iimax = islavept(n,11,iseqr);
        if (iimax > 1) {
            n1(1) = n;
            if (islavept(n,8,iseqr) != 0) {
                int ii3 = 0;
                for (int ii2 = 2; ii2 <= iimax; ii2++) {
                    n1(ii2) = islavept(n,12+ii2-2,iseqr);
                    if (islavept(n1(ii2),8,iseqr) == 0) ii3 = ii2;
                }
                if (ii3 == 0) {
                    double rooiim = 1.0 / (double)iimax;
                    for (int ii2 = 1; ii2 <= iimax; ii2++) {
                        for (int j = 1; j <= 3; j++) {
                            volij(3*(n-1)+j) = volij(3*(n-1)+j) +
                                               u(3*(n1(ii2)-1)+j)*rooiim;
                        }
                    }
                } else {
                    for (int j = 1; j <= 3; j++) {
                        volij(3*(n-1)+j) = u(3*(n1(ii3)-1)+j);
                    }
                }
            } else {
                for (int j = 1; j <= 3; j++) {
                    volij(3*(n-1)+j) = u(3*(n-1)+j);
                }
            }
        }
    }
    for (int n = 1; n <= nslave; n++) {
        int iimax = islavept(n,11,iseqr);
        if (iimax > 1) {
            n1(1) = n;
            for (int ii2 = 2; ii2 <= iimax; ii2++) {
                n1(ii2) = islavept(n,12+ii2-2,iseqr);
            }
            for (int ii2 = 1; ii2 <= iimax; ii2++) {
                for (int j = 1; j <= 3; j++) {
                    u(3*(n1(ii2)-1)+j) = volij(3*(n-1)+j);
                }
            }
        }
    }
}

void gaussseidel(int& nc, int& n, FortranArray1DRef<double> b,
                 FortranArray1DRef<double> x, FortranArray1DRef<double> sa,
                 FortranArray1DRef<int> ija, double& tol, int& itmax,
                 double& err, FortranArray1DRef<double> xt,
                 int& myid, int& myhost)
{
    const double EPS_gs = 1.0e-12;
    int& meshdef = cmn_deformz.meshdef;
    float& alpha1_f = cmn_deformz.alpha1;

    int it1 = 0;
    int itmx1 = 200;
    err = 1.0e+10;
    int nresmx = 1;
    int it = 1;
    for (it = 1; it <= itmax; it++) {
        double sum = 0.0;
        double resid = 0.0;
        for (int i = 1; i <= n; i++) {
            double ut = b(i)/sa(i);
            for (int j = ija(i); j <= ija(i+1)-1; j++) {
                ut = ut - sa(j)*x(ija(j))/sa(i);
            }
            xt(i) = (double)alpha1_f*ut + (1.0-(double)alpha1_f)*x(i);
            sum = sum + (xt(i)-x(i))*(xt(i)-x(i));
            double resid1 = std::sqrt((xt(i)-x(i))*(xt(i)-x(i)));
            if (resid1 > resid) {
                resid = resid1;
                nresmx = i;
            }
            x(i) = xt(i);
        }
        double errold = err;
        err = std::sqrt(sum);
        nresmx = 1+(nresmx-1)/3;
        if (meshdef == 1 && myid == myhost) {
            FILE* fu = fortran_get_unit(1000+myid);
            fprintf(fu, "%8d%16.8e%8d%16.8e\n", it, err, nresmx, resid);
        }
        if (err < tol) goto label2000;
        if (errold < err) {
            it1 = it1 + 1;
            if (it1 > itmx1) {
                FILE* f1000 = fortran_get_unit(1000);
                fprintf(f1000, " Stopping, Gauss-Seidel scheme not converging: \n"
                               " Reduce BETA1 or change Macro-element definition\n");
                exit(1);
            }
        } else {
            it1 = 0;
        }
    }
    label2000:;
    if (myid == myhost) {
        FILE* fu = fortran_get_unit(1000+myid);
        fprintf(fu, "%8d%16.8e\n", it, err);
    }
}



void linbcg(int& nc, int& n, FortranArray1DRef<double> b,
            FortranArray1DRef<double> x, FortranArray1DRef<double> sa,
            FortranArray1DRef<int> ija, int& itol, double& tol,
            int& itmax, int& iter, double& err, double& time,
            double& eps, int& myid)
{
    const double EPS_lcg = 1.0e-12;
    // local arrays
    FortranArray1D<double> p(n), pp(n), r(n), rr(n), z(n), zz(n);
    int itrnsp0 = 0, itrnsp1 = 1;

    iter = 0;
    double test = 0.0;
    for (int j = 1; j <= n; j++) {
        if (test < std::abs(b(j))) test = std::abs(b(j));
    }
    if (test <= eps) return;

    atimes(n, nc, x, r, sa, ija, itrnsp0, myid);
    for (int j = 1; j <= n; j++) {
        r(j)  = b(j) - r(j);
        rr(j) = r(j);
    }
    atimes(n, nc, r, rr, sa, ija, itrnsp0, myid);

    double znrm = 1.0;
    double bnrm = 0.0;
    if (itol == 1) {
        bnrm = snrm(n, b, itol);
    } else if (itol == 2) {
        asolve(n, nc, b, z, sa, itrnsp0);
        bnrm = snrm(n, z, itol);
    } else if (itol == 3 || itol == 4) {
        asolve(n, nc, b, z, sa, itrnsp0);
        bnrm = snrm(n, z, itol);
        asolve(n, nc, r, z, sa, itrnsp0);
        znrm = snrm(n, z, itol);
    } else {
        FILE* fu = fortran_get_unit(1000+myid);
        fprintf(fu, "illegal itol in linbcg\n");
        exit(1);
    }
    double bkden = 0.0;
    if (bnrm < EPS_lcg) {
        err  = 0.0;
        znrm = 0.0;
        goto label200;
    }
    asolve(n, nc, r, z, sa, itrnsp0);

    label100:
    if (iter <= itmax) {
        iter = iter + 1;
        asolve(n, nc, rr, zz, sa, itrnsp1);
        double bknum = 0.0;
        for (int j = 1; j <= n; j++) {
            bknum = bknum + z(j)*rr(j);
        }
        if (iter == 1) {
            for (int j = 1; j <= n; j++) {
                p(j)  = z(j);
                pp(j) = zz(j);
            }
        } else {
            double bk = bknum/bkden;
            for (int j = 1; j <= n; j++) {
                p(j)  = bk*p(j)  + z(j);
                pp(j) = bk*pp(j) + zz(j);
            }
        }
        bkden = bknum;
        atimes(n, nc, p, z, sa, ija, itrnsp0, myid);
        double akden = 0.0;
        for (int j = 1; j <= n; j++) {
            akden = akden + z(j)*pp(j);
        }
        double ak = bknum/akden;
        atimes(n, nc, pp, zz, sa, ija, itrnsp1, myid);

        for (int j = 1; j <= n; j++) {
            x(j)  = x(j)  + ak*p(j);
            r(j)  = r(j)  - ak*z(j);
            rr(j) = rr(j) - ak*zz(j);
        }

        asolve(n, nc, r, z, sa, itrnsp0);
        double zm1nrm = 0.0;
        double dxnrm  = 0.0;
        double xnrm   = 0.0;
        if ((double)itol == 1.0) {
            if (bnrm > EPS_lcg) {
                err = snrm(n, r, itol)/bnrm;
            } else {
                err = 0.0;
            }
            FILE* fu = fortran_get_unit(1000+myid);
            fprintf(fu, "1%8d  %12.5e %12.5e\n", iter, std::log10(err), znrm);
        } else if (itol == 2) {
            if (bnrm > EPS_lcg) {
                err = snrm(n, z, itol)/bnrm;
            } else {
                err = 0.0;
            }
            FILE* fu = fortran_get_unit(1000+myid);
            fprintf(fu, "1%8d  %12.5e %12.5e\n", iter, std::log10(err), znrm);
        } else if (itol == 3 || itol == 4) {
            zm1nrm = znrm;
            znrm   = snrm(n, z, itol);
            if (std::abs(zm1nrm-znrm) > EPS_lcg*znrm) {
                dxnrm = std::abs(ak)*snrm(n, p, itol);
                if (std::abs(zm1nrm-znrm) > EPS_lcg) {
                    err = znrm/std::abs(zm1nrm-znrm)*dxnrm;
                } else {
                    err = 0.0;
                }
            } else {
                if (bnrm > EPS_lcg) {
                    err = znrm/bnrm;
                } else {
                    err  = 0.0;
                    znrm = 0.0;
                }
                goto label100;
            }
            xnrm = snrm(n, x, itol);
            if (err <= 0.5*xnrm) {
                if (xnrm > EPS_lcg) {
                    err = err/xnrm;
                } else {
                    err = 0.0;
                }
                FILE* fu = fortran_get_unit(1000+myid);
                fprintf(fu, "4%8d  %12.5e %12.5e\n", iter, std::log10(err), znrm);
            } else {
                if (bnrm > EPS_lcg) {
                    err = znrm/bnrm;
                } else {
                    err  = 0.0;
                    znrm = 0.0;
                }
                goto label100;
            }
        }
        if (err > tol) goto label100;
    }
    label200:;
}



void atimes(int& n, int& nc, FortranArray1DRef<double> x,
            FortranArray1DRef<double> r, FortranArray1DRef<double> sa,
            FortranArray1DRef<int> ija, int& itrnsp, int& myid)
{
    if (itrnsp == 0) {
        dsprsax(sa, ija, x, r, n, nc, myid);
    } else {
        dsprstx(sa, ija, x, r, n, nc, myid);
    }
}

void dsprstx(FortranArray1DRef<double> sa, FortranArray1DRef<int> ija,
             FortranArray1DRef<double> x, FortranArray1DRef<double> b,
             int& n, int& nc, int& myid)
{
    if (ija(1) != n+2) {
        FILE* fu = fortran_get_unit(1000+myid);
        fprintf(fu, "mismatched vector and matrix in sprstx%8d%8d%8d\n",
                ija(1), n+2, nc);
        exit(1);
    }
    for (int i = 1; i <= n; i++) {
        b(i) = sa(i)*x(i);
    }
    for (int i = 1; i <= n; i++) {
        for (int k = ija(i); k <= ija(i+1)-1; k++) {
            int j = ija(k);
            b(j) = b(j) + sa(k)*x(i);
        }
    }
}

void dprec(FortranArray1DRef<double> sa, FortranArray1DRef<int> ija,
           FortranArray1DRef<double> b, int& n, int& nc)
{
    // Note: Fortran uses myid without declaration (implicit integer = 0)
    int myid_local = 0;
    if (ija(1) != n+2) {
        FILE* fu = fortran_get_unit(1000+myid_local);
        fprintf(fu, "mismatched vector and matrix in sprsax%8d%8d%8d\n",
                ija(1), n+2, nc);
        exit(1);
    }
    for (int i = 1; i <= n; i++) {
        b(i) = b(i)/sa(i);
        for (int k = ija(i); k <= ija(i+1)-1; k++) {
            sa(k) = sa(k)/sa(i);
        }
        sa(i) = 1.0;
    }
}

void dsprsax(FortranArray1DRef<double> sa, FortranArray1DRef<int> ija,
             FortranArray1DRef<double> x, FortranArray1DRef<double> b,
             int& n, int& nc, int& myid)
{
    if (ija(1) != n+2) {
        FILE* fu = fortran_get_unit(1000+myid);
        fprintf(fu, "mismatched vector and matrix in sprsax%8d%8d%8d\n",
                ija(1), n+2, nc);
        exit(1);
    }
    for (int i = 1; i <= n; i++) {
        b(i) = sa(i)*x(i);
        for (int k = ija(i); k <= ija(i+1)-1; k++) {
            b(i) = b(i) + sa(k)*x(ija(k));
        }
    }
}

void asolve(int& n, int& nc, FortranArray1DRef<double> b,
            FortranArray1DRef<double> x, FortranArray1DRef<double> sa,
            int& itrnsp)
{
    for (int i = 1; i <= n; i++) {
        x(i) = b(i)/sa(i);
    }
}

double snrm(int& n, FortranArray1DRef<double> sx, int& itol)
{
    double result = 0.0;
    if (itol <= 3) {
        for (int i = 1; i <= n; i++) {
            result = result + sx(i)*sx(i);
        }
        result = std::sqrt(result);
    } else {
        int isamax = 1;
        for (int i = 1; i <= n; i++) {
            if (std::abs(sx(i)) > std::abs(sx(isamax))) isamax = i;
        }
        result = std::abs(sx(isamax));
    }
    return result;
}



void invdet(FortranArray2DRef<double> UL, FortranArray2DRef<double> COMIN,
            int& N, double& DTNRM, double& DETM)
{
    // Local C(N,N) and J(550)
    FortranArray2D<double> C(N, N);
    int J[551]; // J(550) 1-based: use J[1..550]

    for (int K = 1; K <= N; K++) {
        for (int I = 1; I <= N; I++) {
            C(I,K) = COMIN(I,K);
        }
    }
    double PD = 1.0;
    for (int L = 1; L <= N; L++) {
        double DD = 0.0;
        for (int K = 1; K <= N; K++) {
            DD = DD + C(L,K)*C(L,K);
        }
        DD = std::sqrt(DD);
        PD = PD*DD;
    }
    DETM = 1.0;
    for (int L = 1; L <= N; L++) {
        J[L+20] = L;
    }
    for (int L = 1; L <= N; L++) {
        double CC = 0.0;
        int M = L;
        for (int K = L; K <= N; K++) {
            if ((std::abs(CC) - std::abs(C(L,K))) >= 0.0) continue; // goto 135
            // label 126:
            M = K;
            CC = C(L,K);
        }
        // label 127:
        if (L != M) {
            // label 128:
            int K = J[M+20];
            J[M+20] = J[L+20];
            J[L+20] = K;
            for (int Kk = 1; Kk <= N; Kk++) {
                double S = C(Kk,L);
                C(Kk,L) = C(Kk,M);
                C(Kk,M) = S;
            }
        }
        // label 138:
        C(L,L) = 1.0;
        DETM = DETM*CC;
        for (int Mm = 1; Mm <= N; Mm++) {
            C(L,Mm) = C(L,Mm)/CC;
        }
        for (int Mm = 1; Mm <= N; Mm++) {
            if (L == Mm) continue; // goto 142
            // label 129:
            CC = C(Mm,L);
            if (CC == 0.0) continue; // goto 142
            // label 130:
            C(Mm,L) = 0.0;
            for (int K = 1; K <= N; K++) {
                C(Mm,K) = C(Mm,K) - CC*C(L,K);
            }
        }
    }
    for (int L = 1; L <= N; L++) {
        if (J[L+20] == L) continue; // goto 143
        // label 131:
        int M = L;
        // label 132:
        label132:
        M = M + 1;
        if (J[M+20] == L) goto label133;
        // label 136:
        if (N > M) goto label132;
        label133:
        J[M+20] = J[L+20];
        for (int K = 1; K <= N; K++) {
            double CC = C(L,K);
            C(L,K) = C(M,K);
            C(M,K) = CC;
        }
        J[L+20] = L;
    }
    DETM  = std::abs(DETM);
    DTNRM = DETM/PD;
    for (int k = 1; k <= N; k++) {
        for (int I = 1; I <= N; I++) {
            UL(I,k) = C(I,k);
        }
    }
}



void expdecay(FortranArray1DRef<double> xst, FortranArray1DRef<double> yst,
              FortranArray1DRef<double> zst, FortranArray1DRef<double> u,
              FortranArray1DRef<double> ut, FortranArray1DRef<double> r,
              FortranArray2DRef<double> ripm, double& eps,
              FortranArray3DRef<int> islavept, int& nslave, int& nmaster,
              int& ndim, int& nt, int& nnodes, int& myhost, int& myid,
              int& mycomm, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iseqr)
{
    float& beta2_f  = cmn_deformz.beta2;
    float& alpha2_f = cmn_deformz.alpha2;
    int&   nsprgit  = cmn_deformz.nsprgit;
    float& rmin_f   = cmn_sgminmax.rmin;
    float& rmax_f   = cmn_sgminmax.rmax;

    double relax = 1.0;

    // local fixed arrays
    double arg[8], coef[9], rsn[8];  // 1-based: use [1..7]
    int nipm[8], nsf[8], n3[8];

    // Find closest surface point for each control point
    for (int n1 = 1; n1 <= nslave; n1++) {
        r(n1) = 1.0e6;
        if (islavept(n1,8,iseqr) != 0) {
            for (int n2 = 1; n2 <= nslave; n2++) {
                if (islavept(n2,8,iseqr) == 0) {
                    double test = (xst(n2)-xst(n1))*(xst(n2)-xst(n1))
                                 +(yst(n2)-yst(n1))*(yst(n2)-yst(n1))
                                 +(zst(n2)-zst(n1))*(zst(n2)-zst(n1));
                    test = std::sqrt(test);
                    if (r(n1) > test) {
                        r(n1) = test;
                        islavept(n1,21,iseqr) = n2;
                    }
                }
            }
        } else {
            r(n1) = 0.0;
            islavept(n1,21,iseqr) = n1;
        }
    }
    if (nt == 1) {
        rmax_f = 0.0f;
        rmin_f = 1.0e+30f;
        for (int n = 1; n <= nslave; n++) {
            double rave = r(n);
            if (rave > (double)rmax_f) rmax_f = (float)rave;
            if (rave < (double)rmin_f && rave > eps) rmin_f = (float)rave;
        }
    }
    double rmin = (double)rmin_f;
    double rmax = (double)rmax_f;
    double eps1 = rmin/10000.0;

    for (int n = 1; n <= nslave; n++) {
        nsf[1] = islavept(n,8,iseqr);
        if (nsf[1] != 0) {
            nipm[1] = n;
            nipm[2] = islavept(n,2,iseqr);
            nipm[3] = islavept(n,3,iseqr);
            nipm[4] = islavept(n,4,iseqr);
            nipm[5] = islavept(n,5,iseqr);
            nipm[6] = islavept(n,6,iseqr);
            nipm[7] = islavept(n,7,iseqr);
            for (int ii = 1; ii <= 7; ii++) {
                nsf[ii] = islavept(nipm[ii],8,iseqr);
                if (n != nipm[ii]) {
                    ripm(ii,n) = (xst(n)-xst(nipm[ii]))*(xst(n)-xst(nipm[ii]))
                                +(yst(n)-yst(nipm[ii]))*(yst(n)-yst(nipm[ii]))
                                +(zst(n)-zst(nipm[ii]))*(zst(n)-zst(nipm[ii]))
                                +eps1;
                    ripm(ii,n) = std::sqrt(ripm(ii,n));
                } else {
                    ripm(ii,n) = 1.0e+30;
                }
            }
            for (int ii = 1; ii <= 8; ii++) coef[ii] = 0.0;
            coef[8] = 1.0;
            coef[1] = 1.0;
            for (int ii = 2; ii <= 7; ii++) {
                if (n != nipm[ii]) {
                    coef[ii] = std::exp(-100.0*ripm(ii,n)/rmax);
                    coef[8]  = coef[8] + coef[ii];
                }
            }
            n3[1] = islavept(n,21,iseqr);
            rsn[1] = r(n);
            for (int ii = 2; ii <= 7; ii++) {
                n3[ii]  = islavept(nipm[ii],21,iseqr);
                rsn[ii] = r(nipm[ii]);
            }
            int n2 = n3[1];
            arg[1] = -(double)beta2_f*(rsn[1]/rmax - (double)alpha2_f);
            if (arg[1] > 0.0) arg[1] = 0.0;
            coef[1] = coef[1]/coef[8];
            for (int ii = 2; ii <= 7; ii++) {
                if (n != nipm[ii]) {
                    n2 = n3[ii];
                    arg[ii] = -(double)beta2_f*(rsn[ii]/rmax - (double)alpha2_f);
                    if (arg[ii] > 0.0) arg[ii] = 0.0;
                    coef[ii] = coef[ii]/coef[8];
                } else {
                    coef[ii] = 0.0;
                    arg[ii]  = 0.0;
                }
            }
            for (int j = 1; j <= 3; j++) {
                n2 = n3[1];
                u(3*(n-1)+j) = u(3*(n2-1)+j)*std::exp(arg[1]);
            }
        }
    }

    for (int it = 1; it <= nsprgit; it++) {
        for (int n = 1; n <= nslave; n++) {
            nipm[2] = islavept(n,2,iseqr);
            nipm[3] = islavept(n,3,iseqr);
            nipm[4] = islavept(n,4,iseqr);
            nipm[5] = islavept(n,5,iseqr);
            nipm[6] = islavept(n,6,iseqr);
            nipm[7] = islavept(n,7,iseqr);
            nsf[1]  = islavept(n,8,iseqr);
            if (nsf[1] != 0) {
                double coef1 = 0.0;
                for (int ii = 2; ii <= 7; ii++) {
                    coef1 = coef1 + 1.0/ripm(ii,n);
                }
                for (int j = 1; j <= 3; j++) {
                    ut(3*(n-1)+j) = 0.0;
                    for (int ii = 2; ii <= 7; ii++) {
                        ut(3*(n-1)+j) = ut(3*(n-1)+j)
                            + u(3*(nipm[ii]-1)+j)/ripm(ii,n)/coef1;
                    }
                }
            }
        }
        for (int n = 1; n <= nslave; n++) {
            nipm[2] = islavept(n,2,iseqr);
            nipm[3] = islavept(n,3,iseqr);
            nipm[4] = islavept(n,4,iseqr);
            nipm[5] = islavept(n,5,iseqr);
            nipm[6] = islavept(n,6,iseqr);
            nipm[7] = islavept(n,7,iseqr);
            nsf[1]  = islavept(n,8,iseqr);
            if (nsf[1] != 0) {
                u(3*(n-1)+1) = relax*ut(3*(n-1)+1) + (1.0-relax)*u(3*(n-1)+1);
                u(3*(n-1)+2) = relax*ut(3*(n-1)+2) + (1.0-relax)*u(3*(n-1)+2);
                u(3*(n-1)+3) = relax*ut(3*(n-1)+3) + (1.0-relax)*u(3*(n-1)+3);
            }
        }
    }
}



void hookefe(FortranArray1DRef<double> ei, FortranArray1DRef<double> ej,
             FortranArray1DRef<double> ek, FortranArray1DRef<double> gij,
             FortranArray1DRef<double> gjk, FortranArray1DRef<double> gik,
             FortranArray1DRef<double> xst, FortranArray1DRef<double> yst,
             FortranArray1DRef<double> zst, FortranArray1DRef<double> volij,
             FortranArray1DRef<double> volik, double& eps,
             double& eini, double& gini, double& arg1, double& arg2,
             FortranArray3DRef<int> islavept, int& nslave, int& nmaster,
             int& nt, int& nnodes, int& myhost, int& myid, int& mycomm,
             FortranArray1DRef<int> mblk2nd, int& maxbl, int& iseqr)
{
    float& rmin_f = cmn_sgminmax.rmin;
    float& rmax_f = cmn_sgminmax.rmax;
    float& beta1_f = cmn_deformz.beta1;

    // local arrays (0-based in Fortran: r1(0:8), ni1(0:8), ni0(0:8), test1(0:8))
    FortranArray1D<double> r(nslave);
    double r1[9], test1[9];
    int ni1[9], ni0[9];

    rmax_f = 0.0f;
    rmin_f = 1.0e+30f;

    for (int n1 = 1; n1 <= nslave; n1++) {
        r(n1) = 1.0e6;
        islavept(n1,21,iseqr) = 0;
        ni0[0] = n1;
        ni0[1] = islavept(n1    ,3,iseqr);
        ni0[2] = islavept(n1    ,5,iseqr);
        ni0[3] = islavept(ni0[1],5,iseqr);
        ni0[4] = islavept(n1    ,7,iseqr);
        ni0[5] = islavept(ni0[2],7,iseqr);
        ni0[6] = islavept(ni0[1],7,iseqr);
        ni0[7] = islavept(ni0[6],5,iseqr);
        if (ni0[1] != n1 && ni0[2] != n1 && ni0[4] != n1) {
            double xcent = 0.0, ycent = 0.0, zcent = 0.0;
            for (int ii = 0; ii <= 7; ii++) {
                xcent = xcent + xst(ni0[ii])/8.0;
                ycent = ycent + yst(ni0[ii])/8.0;
                zcent = zcent + zst(ni0[ii])/8.0;
            }
            for (int n2 = 1; n2 <= nslave; n2++) {
                if (islavept(n2,8,iseqr) == 0) {
                    ni1[0] = n2;
                    ni1[1] = islavept(n2    ,3,iseqr);
                    ni1[2] = islavept(n2    ,5,iseqr);
                    ni1[3] = islavept(ni1[1],5,iseqr);
                    ni1[4] = islavept(n2    ,7,iseqr);
                    ni1[5] = islavept(ni1[2],7,iseqr);
                    ni1[6] = islavept(ni1[1],7,iseqr);
                    ni1[7] = islavept(ni1[6],5,iseqr);
                    double test = (xst(n2)-xcent)*(xst(n2)-xcent)
                                 +(yst(n2)-ycent)*(yst(n2)-ycent)
                                 +(zst(n2)-zcent)*(zst(n2)-zcent);
                    test = std::sqrt(test);
                    if (r(n1) > test) {
                        r(n1) = test;
                        islavept(n1,21,iseqr) = n2;
                    }
                    if (islavept(ni1[1],8,iseqr) == 0 && islavept(ni1[2],8,iseqr) == 0) {
                        double xsct = (xst(n2)+xst(ni1[1])+xst(ni1[2])+xst(ni1[3]))/4.0;
                        double ysct = (yst(n2)+yst(ni1[1])+yst(ni1[2])+yst(ni1[3]))/4.0;
                        double zsct = (zst(n2)+zst(ni1[1])+zst(ni1[2])+zst(ni1[3]))/4.0;
                        test = (xsct-xcent)*(xsct-xcent)
                              +(ysct-ycent)*(ysct-ycent)
                              +(zsct-zcent)*(zsct-zcent);
                        test = std::sqrt(test);
                        if (r(n1) > test) {
                            r(n1) = test;
                            islavept(n1,21,iseqr) = n2;
                        }
                    }
                    if (islavept(ni1[2],8,iseqr) == 0 && islavept(ni1[4],8,iseqr) == 0) {
                        double xsct = (xst(n2)+xst(ni1[2])+xst(ni1[4])+xst(ni1[5]))/4.0;
                        double ysct = (yst(n2)+yst(ni1[2])+yst(ni1[4])+yst(ni1[5]))/4.0;
                        double zsct = (zst(n2)+zst(ni1[2])+zst(ni1[4])+zst(ni1[5]))/4.0;
                        test = (xsct-xcent)*(xsct-xcent)
                              +(ysct-ycent)*(ysct-ycent)
                              +(zsct-zcent)*(zsct-zcent);
                        test = std::sqrt(test);
                        if (r(n1) > test) {
                            r(n1) = test;
                            islavept(n1,21,iseqr) = n2;
                        }
                    }
                    if (islavept(ni1[1],8,iseqr) == 0 && islavept(ni1[4],8,iseqr) == 0) {
                        double xsct = (xst(n2)+xst(ni1[1])+xst(ni1[4])+xst(ni1[6]))/4.0;
                        double ysct = (yst(n2)+yst(ni1[1])+yst(ni1[4])+yst(ni1[6]))/4.0;
                        double zsct = (zst(n2)+zst(ni1[1])+zst(ni1[4])+zst(ni1[6]))/4.0;
                        test = (xsct-xcent)*(xsct-xcent)
                              +(ysct-ycent)*(ysct-ycent)
                              +(zsct-zcent)*(zsct-zcent);
                        test = std::sqrt(test);
                        if (r(n1) > test) {
                            r(n1) = test;
                            islavept(n1,21,iseqr) = n2;
                        }
                    }
                }
            }
        }
        if (islavept(n1,21,iseqr) == 0) {
            r(n1) = 0.0;
            islavept(n1,21,iseqr) = n1;
        }
    }
    rmax_f = 0.0f;
    rmin_f = 1.0e+30f;
    for (int n1 = 1; n1 <= nslave; n1++) {
        if (r(n1) > (double)rmax_f) rmax_f = (float)r(n1);
        if (r(n1) < (double)rmin_f && r(n1) > eps) rmin_f = (float)r(n1);
    }
    double rmax = (double)rmax_f;
    double rmin = (double)rmin_f;
    double beta1 = (double)beta1_f;

    double e0 = 500000.0;
    double r0 = (rmax*std::log(1.0-eini/e0)+beta1*rmin)
               /(     std::log(1.0-eini/e0)+beta1     );

    for (int n = 1; n <= nslave; n++) {
        ni0[1] = islavept(n,3,iseqr);
        ni0[2] = islavept(n,5,iseqr);
        ni0[4] = islavept(n,7,iseqr);
        if (ni0[1] != n && ni0[2] != n && ni0[4] != n) {
            double rave  = r(n);
            double denom = (1.0 - std::exp(-beta1*(rave-r0)/(rmax-r0)));
            if (denom < eps) denom = eps;
            double fact = 1.0/denom;
            ej(n)  = eini*fact;
            ek(n)  = eini*fact;
            ei(n)  = eini*fact;
            gij(n) = gini*fact;
            gik(n) = gini*fact;
            gjk(n) = gini*fact;
        }
    }
}



void elmetricsfe(FortranArray1DRef<double> xix, FortranArray1DRef<double> xiy,
                 FortranArray1DRef<double> xiz, FortranArray1DRef<double> etax,
                 FortranArray1DRef<double> etay, FortranArray1DRef<double> etaz,
                 FortranArray1DRef<double> zetax, FortranArray1DRef<double> zetay,
                 FortranArray1DRef<double> zetaz, FortranArray1DRef<double> ooj,
                 FortranArray1DRef<double> xs, FortranArray1DRef<double> ys,
                 FortranArray1DRef<double> zs, double& eps,
                 FortranArray3DRef<int> islavept, int& nslave, int& nmaster,
                 int& nnodes, int& myhost, int& myid, int& mycomm,
                 FortranArray1DRef<int> mblk2nd, int& maxbl, int& iseqr)
{
    for (int n = 1; n <= nslave; n++) {
        int njp    = islavept(n,3,iseqr);
        int nkp    = islavept(n,5,iseqr);
        int njpkp  = islavept(nkp,3,iseqr);
        int nip    = islavept(n,7,iseqr);
        int nipkp  = islavept(nkp,7,iseqr);
        int njpip  = islavept(njp,7,iseqr);
        int njpipkp= islavept(njpip,5,iseqr);
        if (njp != n && nip != n && nkp != n) {
            double xxi   = 0.25*(xs(njp)  -xs(n)  +xs(njpkp)  -xs(nkp)
                                +xs(njpip)-xs(nip)+xs(njpipkp)-xs(nipkp));
            double xeta  = 0.25*(xs(nkp)  -xs(n)  +xs(njpkp)  -xs(njp)
                                +xs(nipkp)-xs(nip)+xs(njpipkp)-xs(njpip));
            double xzeta = 0.25*(xs(nip)  -xs(n)  +xs(nipkp)  -xs(nkp)
                                +xs(njpip)-xs(njp)+xs(njpipkp)-xs(njpkp));
            double yxi   = 0.25*(ys(njp)  -ys(n)  +ys(njpkp)  -ys(nkp)
                                +ys(njpip)-ys(nip)+ys(njpipkp)-ys(nipkp));
            double yeta  = 0.25*(ys(nkp)  -ys(n)  +ys(njpkp)  -ys(njp)
                                +ys(nipkp)-ys(nip)+ys(njpipkp)-ys(njpip));
            double yzeta = 0.25*(ys(nip)  -ys(n)  +ys(nipkp)  -ys(nkp)
                                +ys(njpip)-ys(njp)+ys(njpipkp)-ys(njpkp));
            double zxi   = 0.25*(zs(njp)  -zs(n)  +zs(njpkp)  -zs(nkp)
                                +zs(njpip)-zs(nip)+zs(njpipkp)-zs(nipkp));
            double zeta  = 0.25*(zs(nkp)  -zs(n)  +zs(njpkp)  -zs(njp)
                                +zs(nipkp)-zs(nip)+zs(njpipkp)-zs(njpip));
            double zzeta = 0.25*(zs(nip)  -zs(n)  +zs(nipkp)  -zs(nkp)
                                +zs(njpip)-zs(njp)+zs(njpipkp)-zs(njpkp));
            ooj(n)   = 1.0/(  xxi*(yeta*zzeta-yzeta*zeta)
                             -xeta*(yxi *zzeta-yzeta*zxi )
                            +xzeta*(yxi *zeta -yeta *zxi ));
            xix(n)   = (yeta*zzeta-yzeta*zeta)*ooj(n);
            xiy(n)   = (xzeta*zeta-xeta*zzeta)*ooj(n);
            xiz(n)   = (xeta*yzeta-xzeta*yeta)*ooj(n);
            etax(n)  = (yzeta*zxi -yxi *zzeta)*ooj(n);
            etay(n)  = (xxi*zzeta -xzeta* zxi)*ooj(n);
            etaz(n)  = (xzeta*yxi -xxi *yzeta)*ooj(n);
            zetax(n) = (yxi *zeta -yeta*  zxi)*ooj(n);
            zetay(n) = (xeta*zxi  -xxi * zeta)*ooj(n);
            zetaz(n) = (xxi *yeta -xeta*  yxi)*ooj(n);
        } else {
            ooj(n)   = 0.0;
            xix(n)   = 0.0;
            xiy(n)   = 0.0;
            xiz(n)   = 0.0;
            etax(n)  = 0.0;
            etay(n)  = 0.0;
            etaz(n)  = 0.0;
            zetax(n) = 0.0;
            zetay(n) = 0.0;
            zetaz(n) = 0.0;
        }
    }
}

void elrhs(FortranArray1DRef<double> b, FortranArray1DRef<double> u,
           FortranArray3DRef<int> islavept, int& nslave, int& nmaster,
           int& iseqr)
{
    for (int n = 1; n <= nslave; n++) {
        int ii4 = 1;
        int ni3 = 0;
        if (islavept(n,8,iseqr) == 0) ii4 = 0;
        if (ii4 != 0) {
            int iimax = islavept(n,11,iseqr);
            for (int ii2 = 2; ii2 <= iimax; ii2++) {
                ni3 = islavept(n,12+ii2-2,iseqr);
                if (islavept(ni3,8,iseqr) == 0) {
                    ii4 = 0;
                    goto label1500;
                }
            }
        }
        label1500:;

        if (ii4 == 0 && ni3 == 0) {
            b(3*(n-1)+1) = u(3*(n-1)+1);
            b(3*(n-1)+2) = u(3*(n-1)+2);
            b(3*(n-1)+3) = u(3*(n-1)+3);
        } else if (ii4 == 0 && ni3 != 0) {
            b(3*(n-1)+1) = u(3*(ni3-1)+1);
            b(3*(n-1)+2) = u(3*(ni3-1)+2);
            b(3*(n-1)+3) = u(3*(ni3-1)+3);
        } else {
            b(3*(n-1)+1) = 0.0;
            b(3*(n-1)+2) = 0.0;
            b(3*(n-1)+3) = 0.0;
        }
    }
}



void deform(int& nbl, int& idim, int& jdim, int& kdim,
            FortranArray1DRef<double> x, FortranArray1DRef<double> y, FortranArray1DRef<double> z,
            FortranArray1DRef<double> xnm2, FortranArray1DRef<double> ynm2, FortranArray1DRef<double> znm2,
            FortranArray1DRef<double> xnm1, FortranArray1DRef<double> ynm1, FortranArray1DRef<double> znm1,
            FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk, FortranArray4DRef<double> delti,
            FortranArray1DRef<double> u, FortranArray4DRef<double> vel,
            FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
            FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
            FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
            int& maxbl, double& dt,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid,
            FortranArray1DRef<int> idefrm,
            FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim,
            FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim,
            FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim,
            FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
            int& maxseg, FortranArray1DRef<double> wk, int& nsurf,
            int& irst, int& iflag,
            FortranArray3DRef<int> islavept, int& nslave,
            FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip, FortranArray2DRef<int> kskip,
            FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg,
            FortranArray2DRef<int> iaesurf, int& maxsegdg, int& nmaster,
            int& iseq, FortranArray1DRef<int> iskmax, FortranArray1DRef<int> jskmax,
            FortranArray1DRef<int> kskmax, int& nt)
{
    int& i2d     = cmn_twod.i2d;
    int& meshdef = cmn_deformz.meshdef;
    int& isktyp  = cmn_deformz.isktyp;
    int& ntstep  = cmn_unst.ntstep;
    int& ita     = cmn_unst.ita;
    float& time_f = cmn_unst.time;
    int& movie   = cmn_moov.movie;

    int stats = 0;
    int memuse = 0;

    // Allocatable arrays
    FortranArray3D<double> arci, arcj, arck;
    FortranArray1D<double> dx, dy, dz;
    FortranArray3D<double> dx1, dx2, dx3, dx4, dx5, dx6, dx7;
    FortranArray3D<double> dy1, dy2, dy3, dy4, dy5, dy6, dy7;
    FortranArray3D<double> dz1, dz2, dz3, dz4, dz5, dz6, dz7;
    FortranArray1D<int>    ibl;

    arci.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'a'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    arcj.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'a'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    arck.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'a'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dx.allocate(jdim*kdim*idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dx1.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dx2.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dx3.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dx4.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dx5.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dx6.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dx7.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy.allocate(jdim*kdim*idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy1.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy2.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy3.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy4.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy5.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy6.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy7.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz.allocate(jdim*kdim*idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz1.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz2.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz3.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz4.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz5.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz6.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz7.allocate(jdim,kdim,idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    ibl.allocate(nsurf);
    { int nw=nsurf, iz=1; char txt[1]={'i'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }

    // arclen call - need 3D views of x,y,z
    {
        FortranArray3DRef<double> x3d(&x(1), jdim, kdim, idim);
        FortranArray3DRef<double> y3d(&y(1), jdim, kdim, idim);
        FortranArray3DRef<double> z3d(&z(1), jdim, kdim, idim);
        arclen_ns::arclen(idim, jdim, kdim, arci.ref(), arcj.ref(), arck.ref(),
                          x3d, y3d, z3d,
                          nbl, nou, bou, nbuf, ibufdim, myid);
    }

    // initialize deltas
    for (int j = 1; j <= idim*jdim*kdim; j++) {
        dx(j) = 0.0;
        dy(j) = 0.0;
        dz(j) = 0.0;
    }
    for (int n = 1; n <= nslave; n++) {
        int nbl1 = islavept(n,9,iseq);
        int ll   = islavept(n,1,iseq);
        if (nbl1 == nbl) {
            dx(ll+1) = u(3*(n-1)+1);
            dy(ll+1) = u(3*(n-1)+2);
            dz(ll+1) = u(3*(n-1)+3);
        }
    }



    // meshdef output
    if (meshdef == 1) {
        int inc = ntstep;
        if (movie != 0) {
            double rinc1 = (double)nt/(double)movie;
            int inc1 = nt/movie;
            if ((double)inc1 == rinc1) inc = nt;
        }
        if (nt == inc) {
            FILE* fu = fortran_get_unit(4000+myid);
            fprintf(fu, "ZONE I = %8d  J = %8d K = %8d T=\"Block %5d Time Step = %6d\"\n",
                    jskmax(nbl), kskmax(nbl), iskmax(nbl), nbl, nt);
            for (int n = 1; n <= nslave; n++) {
                int nbl1 = islavept(n,9,iseq);
                int ll   = islavept(n,1,iseq);
                int n2   = islavept(n,21,iseq);
                if (nbl1 == nbl) {
                    fprintf(fu, " %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e%8d%8d\n",
                            x(ll+1)+u(3*(n-1)+1), y(ll+1)+u(3*(n-1)+2),
                            z(ll+1)+u(3*(n-1)+3), u(3*(n-1)+1),
                            u(3*(n-1)+2), u(3*(n-1)+3), n, n2);
                }
            }
        }
    }

    // iskp, jskp, kskp
    int iskp = 1, jskp = 1, kskp = 1;
    if (std::abs(isktyp) == 1) {
        iskp = iskip(nbl,1);
        jskp = jskip(nbl,1);
        kskp = kskip(nbl,1);
    }

    // Build 3D views of dx,dy,dz for passing to subroutines
    FortranArray3DRef<double> dx3d(&dx(1), jdim, kdim, idim);
    FortranArray3DRef<double> dy3d(&dy(1), jdim, kdim, idim);
    FortranArray3DRef<double> dz3d(&dz(1), jdim, kdim, idim);

    // i=1 to idim subfaces (edges)
    if (std::abs(isktyp) == 1) {
        for (int i = 1; i <= idim; i += iskp) {
            int js = 1, je = jdim, ks = 1, ke = kdim;
            for (int j = js; j <= je-jskp; j += jskp) {
                for (int k = ks; k <= ke; k += kskp) {
                    int j2 = j+jskp;
                    tfiedge_ns::tfiedge(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                        i,i,j,j2,k,k,
                                        arci.ref(),arcj.ref(),arck.ref(),
                                        nou,bou,nbuf,ibufdim,myid,nbl);
                }
            }
            for (int j = js; j <= je; j += jskp) {
                for (int k = ks; k <= ke-kskp; k += kskp) {
                    int k2 = k+kskp;
                    tfiedge_ns::tfiedge(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                        i,i,j,j,k,k2,
                                        arci.ref(),arcj.ref(),arck.ref(),
                                        nou,bou,nbuf,ibufdim,myid,nbl);
                }
            }
        }
    } else {
        for (int i1 = 1; i1 <= iskmax(nbl); i1++) {
            int i = iskip(nbl,i1);
            for (int j1 = 1; j1 <= jskmax(nbl)-1; j1++) {
                int j = jskip(nbl,j1);
                jskp = jskip(nbl,j1+1);
                for (int k1 = 1; k1 <= kskmax(nbl); k1++) {
                    int k = kskip(nbl,k1);
                    tfiedge_ns::tfiedge(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                        i,i,j,jskp,k,k,
                                        arci.ref(),arcj.ref(),arck.ref(),
                                        nou,bou,nbuf,ibufdim,myid,nbl);
                }
            }
            for (int j1 = 1; j1 <= jskmax(nbl); j1++) {
                int j = jskip(nbl,j1);
                for (int k1 = 1; k1 <= kskmax(nbl)-1; k1++) {
                    int k = kskip(nbl,k1);
                    kskp = kskip(nbl,k1+1);
                    tfiedge_ns::tfiedge(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                        i,i,j,j,k,kskp,
                                        arci.ref(),arcj.ref(),arck.ref(),
                                        nou,bou,nbuf,ibufdim,myid,nbl);
                }
            }
        }
    }

    // preserve deltas on edges of solid surfaces
    if (idefrm(nbl) < 999) {
        bc_delt_ns::bc_delt(nbl,dx3d,dy3d,dz3d,deltj,deltk,delti,
                            jcsi,jcsf,kcsi,kcsf,icsi,icsf,
                            jdim,kdim,idim,maxbl,maxsegdg,nsegdfrm);
    }



    // i=1 subfaces
    for (int nseg = 1; nseg <= nbci0(nbl); nseg++) {
        int ii = 1;
        int ibctyp = ibcinfo(nbl,nseg,1,ii);
        if (std::abs(ibctyp) != 2004 && std::abs(ibctyp) != 2014 &&
            std::abs(ibctyp) != 2024 && std::abs(ibctyp) != 2034 &&
            std::abs(ibctyp) != 2016 && std::abs(ibctyp) != 1005 &&
            std::abs(ibctyp) != 1006) {
            int js = ibcinfo(nbl,nseg,2,ii);
            int je = ibcinfo(nbl,nseg,3,ii);
            int ks = ibcinfo(nbl,nseg,4,ii);
            int ke = ibcinfo(nbl,nseg,5,ii);
            if (std::abs(isktyp) == 1) {
                for (int j = js; j <= je-jskp; j += jskp) {
                    for (int k = ks; k <= ke-kskp; k += kskp) {
                        int i1v=1, i2v=1, j2=j+jskp, k2=k+kskp;
                        tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            dx1.ref(),dy1.ref(),dz1.ref(),
                                            dx2.ref(),dy2.ref(),dz2.ref(),
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            i1v,i2v,j,j2,k,k2,
                                            nou,bou,nbuf,ibufdim,myid);
                    }
                }
            } else {
                for (int j1 = 1; j1 <= jskmax(nbl)-1; j1++) {
                    int j = jskip(nbl,j1);
                    jskp = jskip(nbl,j1+1);
                    for (int k1 = 1; k1 <= kskmax(nbl)-1; k1++) {
                        int k = kskip(nbl,k1);
                        kskp = kskip(nbl,k1+1);
                        if ((j >= js && j < je) && (k >= ks && k < ke)) {
                            int i1v=1, i2v=1;
                            tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                dx1.ref(),dy1.ref(),dz1.ref(),
                                                dx2.ref(),dy2.ref(),dz2.ref(),
                                                arci.ref(),arcj.ref(),arck.ref(),
                                                i1v,i2v,j,jskp,k,kskp,
                                                nou,bou,nbuf,ibufdim,myid);
                        }
                    }
                }
            }
        }
    }

    // i=idim subfaces
    for (int nseg = 1; nseg <= nbcidim(nbl); nseg++) {
        int ii = 2;
        int ibctyp = ibcinfo(nbl,nseg,1,ii);
        if (std::abs(ibctyp) != 2004 && std::abs(ibctyp) != 2014 &&
            std::abs(ibctyp) != 2024 && std::abs(ibctyp) != 2034 &&
            std::abs(ibctyp) != 2016 && std::abs(ibctyp) != 1005 &&
            std::abs(ibctyp) != 1006) {
            int js = ibcinfo(nbl,nseg,2,ii);
            int je = ibcinfo(nbl,nseg,3,ii);
            int ks = ibcinfo(nbl,nseg,4,ii);
            int ke = ibcinfo(nbl,nseg,5,ii);
            if (std::abs(isktyp) == 1) {
                for (int j = js; j <= je-jskp; j += jskp) {
                    for (int k = ks; k <= ke-kskp; k += kskp) {
                        int j2=j+jskp, k2=k+kskp;
                        tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            dx1.ref(),dy1.ref(),dz1.ref(),
                                            dx2.ref(),dy2.ref(),dz2.ref(),
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            idim,idim,j,j2,k,k2,
                                            nou,bou,nbuf,ibufdim,myid);
                    }
                }
            } else {
                for (int j1 = 1; j1 <= jskmax(nbl)-1; j1++) {
                    int j = jskip(nbl,j1);
                    jskp = jskip(nbl,j1+1);
                    for (int k1 = 1; k1 <= kskmax(nbl)-1; k1++) {
                        int k = kskip(nbl,k1);
                        kskp = kskip(nbl,k1+1);
                        if ((j >= js && j < je) && (k >= ks && k < ke)) {
                            tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                dx1.ref(),dy1.ref(),dz1.ref(),
                                                dx2.ref(),dy2.ref(),dz2.ref(),
                                                arci.ref(),arcj.ref(),arck.ref(),
                                                idim,idim,j,jskp,k,kskp,
                                                nou,bou,nbuf,ibufdim,myid);
                        }
                    }
                }
            }
        }
    }

    // Intermediate i subfaces
    if (idim > 2) {
        if (std::abs(isktyp) == 1) {
            for (int i = 1+iskp; i <= idim-iskp; i++) {
                int js = 1, je = jdim, ks = 1, ke = kdim;
                for (int j = js; j <= je-jskp; j += jskp) {
                    for (int k = ks; k <= ke-kskp; k += kskp) {
                        int j2=j+jskp, k2=k+kskp;
                        tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            dx1.ref(),dy1.ref(),dz1.ref(),
                                            dx2.ref(),dy2.ref(),dz2.ref(),
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            i,i,j,j2,k,k2,
                                            nou,bou,nbuf,ibufdim,myid);
                    }
                }
            }
        } else {
            for (int i1 = 2; i1 <= iskmax(nbl)-1; i1++) {
                int i = iskip(nbl,i1);
                for (int j1 = 1; j1 <= jskmax(nbl)-1; j1++) {
                    int j = jskip(nbl,j1);
                    jskp = jskip(nbl,j1+1);
                    for (int k1 = 1; k1 <= kskmax(nbl)-1; k1++) {
                        int k = kskip(nbl,k1);
                        kskp = kskip(nbl,k1+1);
                        tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            dx1.ref(),dy1.ref(),dz1.ref(),
                                            dx2.ref(),dy2.ref(),dz2.ref(),
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            i,i,j,jskp,k,kskp,
                                            nou,bou,nbuf,ibufdim,myid);
                    }
                }
            }
        }
    }



    if (i2d == 0) {
        // j=1 to jdim subfaces (edges)
        if (std::abs(isktyp) == 1) {
            for (int j = 1; j <= jdim; j += jskp) {
                int is = 1, ie = idim, ks = 1, ke = kdim;
                for (int i = is; i <= ie-iskp; i += iskp) {
                    for (int k = ks; k <= ke; k += kskp) {
                        int i2 = i+iskp;
                        tfiedge_ns::tfiedge(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            i,i2,j,j,k,k,
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            nou,bou,nbuf,ibufdim,myid,nbl);
                    }
                }
            }
        } else {
            for (int j1 = 1; j1 <= jskmax(nbl); j1++) {
                int j = jskip(nbl,j1);
                for (int i1 = 1; i1 <= iskmax(nbl)-1; i1++) {
                    int i = iskip(nbl,i1);
                    iskp = iskip(nbl,i1+1);
                    for (int k1 = 1; k1 <= kskmax(nbl); k1++) {
                        int k = kskip(nbl,k1);
                        tfiedge_ns::tfiedge(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            i,iskp,j,j,k,k,
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            nou,bou,nbuf,ibufdim,myid,nbl);
                    }
                }
            }
        }

        // preserve deltas on edges of solid surfaces
        if (idefrm(nbl) < 999) {
            bc_delt_ns::bc_delt(nbl,dx3d,dy3d,dz3d,deltj,deltk,delti,
                                jcsi,jcsf,kcsi,kcsf,icsi,icsf,
                                jdim,kdim,idim,maxbl,maxsegdg,nsegdfrm);
        }

        // j=1 subfaces
        for (int nseg = 1; nseg <= nbcj0(nbl); nseg++) {
            int jj = 1;
            int jbctyp = jbcinfo(nbl,nseg,1,jj);
            if (std::abs(jbctyp) != 2004 && std::abs(jbctyp) != 2014 &&
                std::abs(jbctyp) != 2024 && std::abs(jbctyp) != 2034 &&
                std::abs(jbctyp) != 2016 && std::abs(jbctyp) != 1005 &&
                std::abs(jbctyp) != 1006) {
                int is = jbcinfo(nbl,nseg,2,jj);
                int ie = jbcinfo(nbl,nseg,3,jj);
                int ks = jbcinfo(nbl,nseg,4,jj);
                int ke = jbcinfo(nbl,nseg,5,jj);
                if (std::abs(isktyp) == 1) {
                    for (int i = is; i <= ie-iskp; i += iskp) {
                        for (int k = ks; k <= ke-kskp; k += kskp) {
                            int j1v=1, j2v=1, i2=i+iskp, k2=k+kskp;
                            tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                dx1.ref(),dy1.ref(),dz1.ref(),
                                                dx2.ref(),dy2.ref(),dz2.ref(),
                                                arci.ref(),arcj.ref(),arck.ref(),
                                                i,i2,j1v,j2v,k,k2,
                                                nou,bou,nbuf,ibufdim,myid);
                        }
                    }
                } else {
                    for (int i1 = 1; i1 <= iskmax(nbl)-1; i1++) {
                        int i = iskip(nbl,i1);
                        iskp = iskip(nbl,i1+1);
                        for (int k1 = 1; k1 <= kskmax(nbl)-1; k1++) {
                            int k = kskip(nbl,k1);
                            kskp = kskip(nbl,k1+1);
                            if ((i >= is && i < ie) && (k >= ks && k < ke)) {
                                int j1v=1, j2v=1;
                                tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                    dx1.ref(),dy1.ref(),dz1.ref(),
                                                    dx2.ref(),dy2.ref(),dz2.ref(),
                                                    arci.ref(),arcj.ref(),arck.ref(),
                                                    i,iskp,j1v,j2v,k,kskp,
                                                    nou,bou,nbuf,ibufdim,myid);
                            }
                        }
                    }
                }
            }
        }

        // j=jdim subfaces
        // j=jdim subfaces
        for (int nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
            int jj = 2;
            int jbctyp = jbcinfo(nbl,nseg,1,jj);
            if (std::abs(jbctyp) != 2004 && std::abs(jbctyp) != 2014 &&
                std::abs(jbctyp) != 2024 && std::abs(jbctyp) != 2034 &&
                std::abs(jbctyp) != 2016 && std::abs(jbctyp) != 1005 &&
                std::abs(jbctyp) != 1006) {
                int is = jbcinfo(nbl,nseg,2,jj);
                int ie = jbcinfo(nbl,nseg,3,jj);
                int ks = jbcinfo(nbl,nseg,4,jj);
                int ke = jbcinfo(nbl,nseg,5,jj);
                if (std::abs(isktyp) == 1) {
                    for (int i = is; i <= ie-iskp; i += iskp) {
                        for (int k = ks; k <= ke-kskp; k += kskp) {
                            int i2=i+iskp, k2=k+kskp;
                            tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                dx1.ref(),dy1.ref(),dz1.ref(),
                                                dx2.ref(),dy2.ref(),dz2.ref(),
                                                arci.ref(),arcj.ref(),arck.ref(),
                                                i,i2,jdim,jdim,k,k2,
                                                nou,bou,nbuf,ibufdim,myid);
                        }
                    }
                } else {
                    for (int i1 = 1; i1 <= iskmax(nbl)-1; i1++) {
                        int i = iskip(nbl,i1);
                        iskp = iskip(nbl,i1+1);
                        for (int k1 = 1; k1 <= kskmax(nbl)-1; k1++) {
                            int k = kskip(nbl,k1);
                            kskp = kskip(nbl,k1+1);
                            if ((i >= is && i < ie) && (k >= ks && k < ke)) {
                                tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                    dx1.ref(),dy1.ref(),dz1.ref(),
                                                    dx2.ref(),dy2.ref(),dz2.ref(),
                                                    arci.ref(),arcj.ref(),arck.ref(),
                                                    i,iskp,jdim,jdim,k,kskp,
                                                    nou,bou,nbuf,ibufdim,myid);
                            }
                        }
                    }
                }
            }
        }



        // Intermediate j subfaces
        if (std::abs(isktyp) == 1) {
            for (int j = 1+jskp; j <= jdim-jskp; j++) {
                int is = 1, ie = idim, ks = 1, ke = kdim;
                for (int i = is; i <= ie-iskp; i += iskp) {
                    for (int k = ks; k <= ke-kskp; k += kskp) {
                        int i2=i+iskp, k2=k+kskp;
                        tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            dx1.ref(),dy1.ref(),dz1.ref(),
                                            dx2.ref(),dy2.ref(),dz2.ref(),
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            i,i2,j,j,k,k2,
                                            nou,bou,nbuf,ibufdim,myid);
                    }
                }
            }
        } else {
            for (int j1 = 2; j1 <= jskmax(nbl)-1; j1++) {
                int j = jskip(nbl,j1);
                for (int i1 = 1; i1 <= iskmax(nbl)-1; i1++) {
                    int i = iskip(nbl,i1);
                    iskp = iskip(nbl,i1+1);
                    for (int k1 = 1; k1 <= kskmax(nbl)-1; k1++) {
                        int k = kskip(nbl,k1);
                        kskp = kskip(nbl,k1+1);
                        tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            dx1.ref(),dy1.ref(),dz1.ref(),
                                            dx2.ref(),dy2.ref(),dz2.ref(),
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            i,iskp,j,j,k,kskp,
                                            nou,bou,nbuf,ibufdim,myid);
                    }
                }
            }
        }

        // k=1 subfaces
        for (int nseg = 1; nseg <= nbck0(nbl); nseg++) {
            int kk = 1;
            int kbctyp = kbcinfo(nbl,nseg,1,kk);
            if (std::abs(kbctyp) != 2004 && std::abs(kbctyp) != 2014 &&
                std::abs(kbctyp) != 2024 && std::abs(kbctyp) != 2034 &&
                std::abs(kbctyp) != 2016 && std::abs(kbctyp) != 1005 &&
                std::abs(kbctyp) != 1006) {
                int is = kbcinfo(nbl,nseg,2,kk);
                int ie = kbcinfo(nbl,nseg,3,kk);
                int js = kbcinfo(nbl,nseg,4,kk);
                int je = kbcinfo(nbl,nseg,5,kk);
                if (std::abs(isktyp) == 1) {
                    for (int i = is; i <= ie-iskp; i += iskp) {
                        for (int j = js; j <= je-jskp; j += jskp) {
                            int k1v=1, k2v=1, i2=i+iskp, j2=j+jskp;
                            tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                dx1.ref(),dy1.ref(),dz1.ref(),
                                                dx2.ref(),dy2.ref(),dz2.ref(),
                                                arci.ref(),arcj.ref(),arck.ref(),
                                                i,i2,j,j2,k1v,k2v,
                                                nou,bou,nbuf,ibufdim,myid);
                        }
                    }
                } else {
                    for (int i1 = 1; i1 <= iskmax(nbl)-1; i1++) {
                        int i = iskip(nbl,i1);
                        iskp = iskip(nbl,i1+1);
                        for (int j1 = 1; j1 <= jskmax(nbl)-1; j1++) {
                            int j = jskip(nbl,j1);
                            jskp = jskip(nbl,j1+1);
                            if ((i >= is && i < ie) && (j >= js && j < je)) {
                                int k1v=1, k2v=1;
                                tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                    dx1.ref(),dy1.ref(),dz1.ref(),
                                                    dx2.ref(),dy2.ref(),dz2.ref(),
                                                    arci.ref(),arcj.ref(),arck.ref(),
                                                    i,iskp,j,jskp,k1v,k2v,
                                                    nou,bou,nbuf,ibufdim,myid);
                            }
                        }
                    }
                }
            }
        }

        // k=kdim subfaces
        for (int nseg = 1; nseg <= nbckdim(nbl); nseg++) {
            int kk = 2;
            int kbctyp = kbcinfo(nbl,nseg,1,kk);
            if (std::abs(kbctyp) != 2004 && std::abs(kbctyp) != 2014 &&
                std::abs(kbctyp) != 2024 && std::abs(kbctyp) != 2034 &&
                std::abs(kbctyp) != 2016 && std::abs(kbctyp) != 1005 &&
                std::abs(kbctyp) != 1006) {
                int is = kbcinfo(nbl,nseg,2,kk);
                int ie = kbcinfo(nbl,nseg,3,kk);
                int js = kbcinfo(nbl,nseg,4,kk);
                int je = kbcinfo(nbl,nseg,5,kk);
                if (std::abs(isktyp) == 1) {
                    for (int i = is; i <= ie-iskp; i += iskp) {
                        for (int j = js; j <= je-jskp; j += jskp) {
                            int i2=i+iskp, j2=j+jskp;
                            tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                dx1.ref(),dy1.ref(),dz1.ref(),
                                                dx2.ref(),dy2.ref(),dz2.ref(),
                                                arci.ref(),arcj.ref(),arck.ref(),
                                                i,i2,j,j2,kdim,kdim,
                                                nou,bou,nbuf,ibufdim,myid);
                        }
                    }
                } else {
                    for (int i1 = 1; i1 <= iskmax(nbl)-1; i1++) {
                        int i = iskip(nbl,i1);
                        iskp = iskip(nbl,i1+1);
                        for (int j1 = 1; j1 <= jskmax(nbl)-1; j1++) {
                            int j = jskip(nbl,j1);
                            jskp = jskip(nbl,j1+1);
                            if ((i >= is && i < ie) && (j >= js && j < je)) {
                                tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                                    dx1.ref(),dy1.ref(),dz1.ref(),
                                                    dx2.ref(),dy2.ref(),dz2.ref(),
                                                    arci.ref(),arcj.ref(),arck.ref(),
                                                    i,iskp,j,jskp,kdim,kdim,
                                                    nou,bou,nbuf,ibufdim,myid);
                            }
                        }
                    }
                }
            }
        }



        // Intermediate k subfaces
        if (std::abs(isktyp) == 1) {
            for (int k = 1+kskp; k <= kdim-kskp; k++) {
                int is = 1, ie = idim, js = 1, je = jdim;
                for (int i = is; i <= ie-iskp; i += iskp) {
                    for (int j = js; j <= je-jskp; j += jskp) {
                        int i2=i+iskp, j2=j+jskp;
                        tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            dx1.ref(),dy1.ref(),dz1.ref(),
                                            dx2.ref(),dy2.ref(),dz2.ref(),
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            i,i2,j,j2,k,k,
                                            nou,bou,nbuf,ibufdim,myid);
                    }
                }
            }
        } else {
            for (int k1 = 2; k1 <= kskmax(nbl)-1; k1++) {
                int k = kskip(nbl,k1);
                for (int i1 = 1; i1 <= iskmax(nbl)-1; i1++) {
                    int i = iskip(nbl,i1);
                    iskp = iskip(nbl,i1+1);
                    for (int j1 = 1; j1 <= jskmax(nbl)-1; j1++) {
                        int j = jskip(nbl,j1);
                        jskp = jskip(nbl,j1+1);
                        tfiface_ns::tfiface(idim,jdim,kdim,dx3d,dy3d,dz3d,
                                            dx1.ref(),dy1.ref(),dz1.ref(),
                                            dx2.ref(),dy2.ref(),dz2.ref(),
                                            arci.ref(),arcj.ref(),arck.ref(),
                                            i,iskp,j,jskp,k,k,
                                            nou,bou,nbuf,ibufdim,myid);
                    }
                }
            }
        }

        // TFI to get deltas throughout the volume
        tfivol_ns::tfivol(idim,jdim,kdim,iskp,jskp,kskp,
                          iskmax,jskmax,kskmax,iskip,jskip,kskip,
                          isktyp,dx3d,dy3d,dz3d,
                          dx1.ref(),dy1.ref(),dz1.ref(),
                          dx2.ref(),dy2.ref(),dz2.ref(),
                          dx3.ref(),dy3.ref(),dz3.ref(),
                          dx4.ref(),dy4.ref(),dz4.ref(),
                          dx5.ref(),dy5.ref(),dz5.ref(),
                          dx6.ref(),dy6.ref(),dz6.ref(),
                          dx7.ref(),dy7.ref(),dz7.ref(),
                          arci.ref(),arcj.ref(),arck.ref(),
                          nou,bou,nbuf,ibufdim,myid,maxbl,nbl);

    } // end if (i2d == 0)

    // compute grid velocities before overwriting current grid
    if (irst == 0) {
        double time_d = (double)time_f;
        if ((float)time_d <= (float)dt || std::abs(ita) == 1) {
            // first order backward derivatives
            double fact = 1.0/dt;
            for (int i = 1; i <= idim; i++) {
                for (int k = 1; k <= kdim; k++) {
                    for (int j = 1; j <= jdim; j++) {
                        int jj = j + (k-1)*jdim + (i-1)*jdim*kdim;
                        vel(j,k,i,1) = fact*((dx(jj)+x(jj)) - xnm1(jj));
                        vel(j,k,i,2) = fact*((dy(jj)+y(jj)) - ynm1(jj));
                        vel(j,k,i,3) = fact*((dz(jj)+z(jj)) - znm1(jj));
                    }
                }
            }
        } else {
            // second order backward derivatives
            double fact = 1.0/(2.0*dt);
            for (int i = 1; i <= idim; i++) {
                for (int k = 1; k <= kdim; k++) {
                    for (int j = 1; j <= jdim; j++) {
                        int jj = j + (k-1)*jdim + (i-1)*jdim*kdim;
                        vel(j,k,i,1) = fact*(3.0*(dx(jj)+x(jj)) - 4.0*xnm1(jj)+xnm2(jj));
                        vel(j,k,i,2) = fact*(3.0*(dy(jj)+y(jj)) - 4.0*ynm1(jj)+ynm2(jj));
                        vel(j,k,i,3) = fact*(3.0*(dz(jj)+z(jj)) - 4.0*znm1(jj)+znm2(jj));
                    }
                }
            }
        }
    }

    // for second order case, store current--->old before updating grid
    if (std::abs(ita) > 1) {
        for (int i = 1; i <= jdim*kdim*idim; i++) {
            xnm2(i) = xnm1(i);
            ynm2(i) = ynm1(i);
            znm2(i) = znm1(i);
        }
    }

    // add deltas to the current mesh to get the new one
    for (int i = 1; i <= jdim*kdim*idim; i++) {
        x(i) = x(i) + dx(i);
        y(i) = y(i) + dy(i);
        z(i) = z(i) + dz(i);
    }
    // release memory (automatic via FortranArray destructors)
}



void deform_surf(int& nbl, int& idim, int& jdim, int& kdim,
                 FortranArray4DRef<double> deltj, FortranArray4DRef<double> deltk,
                 FortranArray4DRef<double> delti,
                 FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
                 FortranArray2DRef<int> icsi, FortranArray2DRef<int> icsf,
                 FortranArray2DRef<int> jcsi, FortranArray2DRef<int> jcsf,
                 FortranArray2DRef<int> kcsi, FortranArray2DRef<int> kcsf,
                 int& maxbl, int& mseq, double& time, double& dt, int& ita,
                 FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
                 int& nbuf, int& ibufdim, int& myid,
                 FortranArray1DRef<int> idefrm,
                 FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim,
                 FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim,
                 FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim,
                 FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
                 FortranArray4DRef<int> kbcinfo,
                 int& maxseg, FortranArray1DRef<double> wk,
                 FortranArray1DRef<double> u, int& nsurf,
                 int& irst, int& iflag,
                 FortranArray3DRef<int> islavept, int& nslave,
                 FortranArray1DRef<int> nsegdfrm, FortranArray2DRef<int> idfrmseg,
                 FortranArray2DRef<int> iaesurf, int& maxsegdg, int& nmaster,
                 int& iseq)
{
    int& i2d = cmn_twod.i2d;

    int stats = 0;
    int memuse = 0;

    FortranArray1D<double> dx, dy, dz;
    dx.allocate(jdim*kdim*idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dy.allocate(jdim*kdim*idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }
    dz.allocate(jdim*kdim*idim);
    { int nw=jdim*kdim*idim, iz=0; char txt[1]={'d'}; umalloc_ns::umalloc(nw,iz,txt,memuse,stats); }

    // initialize deltas
    for (int i = 1; i <= jdim*kdim*idim; i++) {
        dx(i) = 0.0;
        dy(i) = 0.0;
        dz(i) = 0.0;
    }

    // preserve deltas on edges of solid surfaces
    if (idefrm(nbl) < 999) {
        FortranArray3DRef<double> dx3d(&dx(1), jdim, kdim, idim);
        FortranArray3DRef<double> dy3d(&dy(1), jdim, kdim, idim);
        FortranArray3DRef<double> dz3d(&dz(1), jdim, kdim, idim);
        bc_delt_ns::bc_delt(nbl,dx3d,dy3d,dz3d,deltj,deltk,delti,
                            jcsi,jcsf,kcsi,kcsf,icsi,icsf,
                            jdim,kdim,idim,maxbl,maxsegdg,nsegdfrm);

        for (int n = 1; n <= nslave; n++) {
            int nbl1 = islavept(n,9,iseq);
            int ll   = islavept(n,1,iseq);
            int isrf = islavept(n,8,iseq);
            if (nbl1 == nbl && isrf == 0) {
                u(3*(n-1)+1) = dx(ll+1);
                u(3*(n-1)+2) = dy(ll+1);
                u(3*(n-1)+3) = dz(ll+1);
            }
        }
    }
    // release memory (automatic via FortranArray destructors)
}


} // namespace deform_ns
