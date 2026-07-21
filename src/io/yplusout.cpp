// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "yplusout.h"
#include "runtime/fortran_io.h"
#include "lead.h"
#include "bc.h"
#include "qface.h"
#include "calyplus.h"
#include <cmath>
#include <cstdio>

namespace yplusout_ns {

void yplusout(int& iseq, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2, FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork, double& cl, int& maxbl, int& maxgr, int& maxseg, int& nblock, FortranArray3DRef<int> lwdat, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek, FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm, FortranArray1DRef<int> nblg, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray1DRef<int> iadvance, FortranArray1DRef<int> iovrlp, int& myid, int& myhost, int& mycomm, FortranArray1DRef<int> mblk2nd, FortranArray1DRef<int> idf, FortranArray1DRef<int> jdf, FortranArray1DRef<int> kdf, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, FortranArray1DRef<char[80]> bcfiles, int& mxbcfil, int& nummem)
{
    // COMMON block aliases
    int32_t& jdim    = cmn_ginfo.jdim;
    int32_t& kdim    = cmn_ginfo.kdim;
    int32_t& idim    = cmn_ginfo.idim;
    int32_t& lq      = cmn_ginfo.lq;
    int32_t& lqj0    = cmn_ginfo.lqj0;
    int32_t& lqk0    = cmn_ginfo.lqk0;
    int32_t& lqi0    = cmn_ginfo.lqi0;
    int32_t& lsj     = cmn_ginfo.lsj;
    int32_t& lsk     = cmn_ginfo.lsk;
    int32_t& lsi     = cmn_ginfo.lsi;
    int32_t& lvol    = cmn_ginfo.lvol;
    int32_t& lx      = cmn_ginfo.lx;
    int32_t& ly      = cmn_ginfo.ly;
    int32_t& lz      = cmn_ginfo.lz;
    int32_t& lvis    = cmn_ginfo.lvis;
    int32_t& lblk    = cmn_ginfo.lblk;
    int32_t& lvj0    = cmn_ginfo.lvj0;
    int32_t& lvk0    = cmn_ginfo.lvk0;
    int32_t& lvi0    = cmn_ginfo.lvi0;
    int32_t& lbcj    = cmn_ginfo.lbcj;
    int32_t& lbck    = cmn_ginfo.lbck;
    int32_t& lbci    = cmn_ginfo.lbci;

    int32_t& levt    = cmn_mgrd.levt;
    int32_t& lglobal = cmn_mgrd.lglobal;

    int32_t& ivmx    = cmn_maxiv.ivmx;
    int32_t& i2d     = cmn_twod.i2d;
    int32_t& ialph   = cmn_igrdtyp.ialph;

    int32_t* ivisc   = cmn_reyue.ivisc;  // 1-based: ivisc[0]=ivisc(1), etc.

    float& ypsumb    = cmn_ypinfo.ypsumb;
    float& ypsumsqb  = cmn_ypinfo.ypsumsqb;
    float& ypmaxb    = cmn_ypinfo.ypmaxb;
    float& ypminb    = cmn_ypinfo.ypminb;
    float& dnmaxb    = cmn_ypinfo.dnmaxb;
    float& dnminb    = cmn_ypinfo.dnminb;
    float& ypchk     = cmn_ypinfo.ypchk;
    int32_t& nptsb   = cmn_ypinfo.nptsb;
    int32_t& jypmaxb = cmn_ypinfo.jypmaxb;
    int32_t& kypmaxb = cmn_ypinfo.kypmaxb;
    int32_t& iypmaxb = cmn_ypinfo.iypmaxb;
    int32_t& jypminb = cmn_ypinfo.jypminb;
    int32_t& kypminb = cmn_ypinfo.kypminb;
    int32_t& iypminb = cmn_ypinfo.iypminb;
    int32_t& jdnmaxb = cmn_ypinfo.jdnmaxb;
    int32_t& kdnmaxb = cmn_ypinfo.kdnmaxb;
    int32_t& idnmaxb = cmn_ypinfo.idnmaxb;
    int32_t& jdnminb = cmn_ypinfo.jdnminb;
    int32_t& kdnminb = cmn_ypinfo.kdnminb;
    int32_t& idnminb = cmn_ypinfo.idnminb;
    int32_t& nypchkb = cmn_ypinfo.nypchkb;

    int32_t* ncyc1   = cmn_info.ncyc1;  // 1-based array in struct

    float& edvislim  = cmn_turbconv.edvislim;
    int32_t& iwarneddy = cmn_turbconv.iwarneddy;

    // Local variables
    int ifunc, iunit;
    float ypsumt, ypsumsqt, ypmaxt, ypmint, dnmaxt, dnmint;
    int iypmaxt, jypmaxt, kypmaxt;
    int iypmint, jypmint, kypmint;
    int idnmaxt, jdnmaxt, kdnmaxt;
    int idnmint, jdnmint, kdnmint;
    int nypchkt, nptst;
    int nbypmaxt, nbypmint, nbdnmaxt, nbdnmint;
    int iflag;
    int nbl, ns, m, nbc, nseg;
    int ivv1, ivv2, ivv;
    int levg;
    int ncount;
    int ktest, jtest, itest;
    int ldim;
    int lypj, lypk, lypi, lblj_local, lblk_local, lbli;
    int ldnj, ldnk, ldni, lvtj, lvtk, lvti;
    float ypavt, ypstdt;



    // check for number of points with yplus > ypchk
    ypchk = 5.f;

    // ifunc - flag to generate plot3d function files
    ifunc = 0;
    iunit  = 28;

    ypsumt   = 0.f;
    ypsumsqt = 0.f;
    ypmaxt   = 0.f;
    iypmaxt  = 0;
    jypmaxt  = 0;
    kypmaxt  = 0;
    ypmint   = 1.e9f;
    iypmint  = 0;
    jypmint  = 0;
    kypmint  = 0;
    dnmaxt   = 0.f;
    idnmaxt  = 0;
    jdnmaxt  = 0;
    kdnmaxt  = 0;
    dnmint   = 1.e9f;
    idnmint  = 0;
    jdnmint  = 0;
    kdnmint  = 0;
    nypchkt  = 0;
    nptst    = 0;
    nbypmaxt = 0;
    nbypmint = 0;
    nbdnmaxt = 0;
    nbdnmint = 0;
    iflag    = 0;

    if (myid == myhost) {

        if (ifunc > 0) {
            if (i2d == 0) {
                fortran_open_unit(iunit,   "surf_xyz.fmt",  "w");
                fortran_open_unit(iunit+1, "surf_y+.fmt",   "w");
                fortran_open_unit(iunit+2, "surf_dn.fmt",   "w");
                fortran_open_unit(iunit+3, "surf_vist.fmt", "w");
            } else {
                fortran_open_unit(iunit, "surf_y+_2d.fmt", "w");
            }
            int nvar = 1;

            // count number of surfaces to be output to function file
            ncount = 0;
            for (nbl = 1; nbl <= nblock; nbl++) {
                lead_ns::lead(nbl, lw, lw2, maxbl);

                ivv1 = 0;
                if (ivisc[2] > 1 || ivisc[1] > 1 || ivisc[0] > 1) ivv1 = 1;

                ivv2 = 0;
                for (ns = 1; ns <= nbci0(nbl); ns++) {
                    if (std::abs(ibcinfo(nbl,ns,1,1)) == 2004 ||
                        std::abs(ibcinfo(nbl,ns,1,1)) == 2014 ||
                        std::abs(ibcinfo(nbl,ns,1,1)) == 2024 ||
                        std::abs(ibcinfo(nbl,ns,1,1)) == 2034 ||
                        std::abs(ibcinfo(nbl,ns,1,1)) == 2016) ivv2 = ivv2 + 1;
                }
                for (ns = 1; ns <= nbcidim(nbl); ns++) {
                    if (std::abs(ibcinfo(nbl,ns,1,2)) == 2004 ||
                        std::abs(ibcinfo(nbl,ns,1,2)) == 2014 ||
                        std::abs(ibcinfo(nbl,ns,1,2)) == 2024 ||
                        std::abs(ibcinfo(nbl,ns,1,2)) == 2034 ||
                        std::abs(ibcinfo(nbl,ns,1,2)) == 2016) ivv2 = ivv2 + 1;
                }
                for (ns = 1; ns <= nbcj0(nbl); ns++) {
                    if (std::abs(jbcinfo(nbl,ns,1,1)) == 2004 ||
                        std::abs(jbcinfo(nbl,ns,1,1)) == 2014 ||
                        std::abs(jbcinfo(nbl,ns,1,1)) == 2024 ||
                        std::abs(jbcinfo(nbl,ns,1,1)) == 2034 ||
                        std::abs(jbcinfo(nbl,ns,1,1)) == 2016) ivv2 = ivv2 + 1;
                }
                for (ns = 1; ns <= nbcjdim(nbl); ns++) {
                    if (std::abs(jbcinfo(nbl,ns,1,2)) == 2004 ||
                        std::abs(jbcinfo(nbl,ns,1,2)) == 2014 ||
                        std::abs(jbcinfo(nbl,ns,1,2)) == 2024 ||
                        std::abs(jbcinfo(nbl,ns,1,2)) == 2034 ||
                        std::abs(jbcinfo(nbl,ns,1,2)) == 2016) ivv2 = ivv2 + 1;
                }
                for (ns = 1; ns <= nbck0(nbl); ns++) {
                    if (std::abs(kbcinfo(nbl,ns,1,1)) == 2004 ||
                        std::abs(kbcinfo(nbl,ns,1,1)) == 2014 ||
                        std::abs(kbcinfo(nbl,ns,1,1)) == 2024 ||
                        std::abs(kbcinfo(nbl,ns,1,1)) == 2034 ||
                        std::abs(kbcinfo(nbl,ns,1,1)) == 2016) ivv2 = ivv2 + 1;
                }
                for (ns = 1; ns <= nbckdim(nbl); ns++) {
                    if (std::abs(kbcinfo(nbl,ns,1,2)) == 2004 ||
                        std::abs(kbcinfo(nbl,ns,1,2)) == 2014 ||
                        std::abs(kbcinfo(nbl,ns,1,2)) == 2024 ||
                        std::abs(kbcinfo(nbl,ns,1,2)) == 2034 ||
                        std::abs(kbcinfo(nbl,ns,1,2)) == 2016) ivv2 = ivv2 + 1;
                }

                ivv = 0;
                if (ivv1 > 0 && ivv2 > 0) ivv = 1;

                levg = levelg(nbl);



                if (levg >= lglobal && levg <= levt
                    && ncyc1[iseq-1] > 0 && ivv > 0) {

                    if (ivisc[2] > 1) {
                        for (m = 1; m <= 2; m++) {
                            if (m == 1) {
                                nbc = nbck0(nbl);
                            } else {
                                nbc = nbckdim(nbl);
                            }
                            ktest = 0;
                            for (nseg = 1; nseg <= nbc; nseg++) {
                                if (std::abs(kbcinfo(nbl,nseg,1,m)) == 2004 ||
                                    std::abs(kbcinfo(nbl,nseg,1,m)) == 2014 ||
                                    std::abs(kbcinfo(nbl,nseg,1,m)) == 2024 ||
                                    std::abs(kbcinfo(nbl,nseg,1,m)) == 2034 ||
                                    std::abs(kbcinfo(nbl,nseg,1,m)) == 2016) ktest = 1;
                            }
                            if (ktest > 0) {
                                ncount = ncount + 1;
                                idf(ncount) = idim;
                                if (i2d == 1) idf(ncount) = 1;
                                jdf(ncount) = jdim;
                                kdf(ncount) = 1;
                            }
                        }
                    }
                    if (ivisc[1] > 1) {
                        for (m = 1; m <= 2; m++) {
                            if (m == 1) {
                                nbc = nbcj0(nbl);
                            } else {
                                nbc = nbcjdim(nbl);
                            }
                            jtest = 0;
                            for (nseg = 1; nseg <= nbc; nseg++) {
                                if (std::abs(jbcinfo(nbl,nseg,1,m)) == 2004 ||
                                    std::abs(jbcinfo(nbl,nseg,1,m)) == 2014 ||
                                    std::abs(jbcinfo(nbl,nseg,1,m)) == 2024 ||
                                    std::abs(jbcinfo(nbl,nseg,1,m)) == 2034 ||
                                    std::abs(jbcinfo(nbl,nseg,1,m)) == 2016) jtest = 1;
                            }
                            if (jtest > 0) {
                                ncount = ncount + 1;
                                idf(ncount) = idim;
                                if (i2d == 1) idf(ncount) = 1;
                                jdf(ncount) = 1;
                                kdf(ncount) = kdim;
                            }
                        }
                    }
                    if (ivisc[0] > 1) {
                        for (m = 1; m <= 2; m++) {
                            if (m == 1) {
                                nbc = nbci0(nbl);
                            } else {
                                nbc = nbcidim(nbl);
                            }
                            itest = 0;
                            for (nseg = 1; nseg <= nbc; nseg++) {
                                if (std::abs(ibcinfo(nbl,nseg,1,m)) == 2004 ||
                                    std::abs(ibcinfo(nbl,nseg,1,m)) == 2014 ||
                                    std::abs(ibcinfo(nbl,nseg,1,m)) == 2024 ||
                                    std::abs(ibcinfo(nbl,nseg,1,m)) == 2034 ||
                                    std::abs(ibcinfo(nbl,nseg,1,m)) == 2016) itest = 1;
                            }
                            if (itest > 0) {
                                ncount = ncount + 1;
                                idf(ncount) = 1;
                                jdf(ncount) = jdim;
                                kdf(ncount) = kdim;
                            }
                        }
                    }

                } // end if levg...

            } // end do 600 nbl

            if (ncount > 0) {
                if (i2d == 0) {
                    FILE* fu = fortran_get_unit(iunit);
                    fprintf(fu, " %d\n", ncount);
                    for (int l = 1; l <= ncount; l++) {
                        fprintf(fu, " %d %d %d", idf(l), jdf(l), kdf(l));
                    }
                    fprintf(fu, "\n");
                    FILE* fu1 = fortran_get_unit(iunit+1);
                    fprintf(fu1, " %d\n", ncount);
                    for (int l = 1; l <= ncount; l++) {
                        fprintf(fu1, " %d %d %d %d", idf(l), jdf(l), kdf(l), nvar);
                    }
                    fprintf(fu1, "\n");
                    FILE* fu2 = fortran_get_unit(iunit+2);
                    fprintf(fu2, " %d\n", ncount);
                    for (int l = 1; l <= ncount; l++) {
                        fprintf(fu2, " %d %d %d %d", idf(l), jdf(l), kdf(l), nvar);
                    }
                    fprintf(fu2, "\n");
                    FILE* fu3 = fortran_get_unit(iunit+3);
                    fprintf(fu3, " %d\n", ncount);
                    for (int l = 1; l <= ncount; l++) {
                        fprintf(fu3, " %d %d %d %d", idf(l), jdf(l), kdf(l), nvar);
                    }
                    fprintf(fu3, "\n");
                } else {
                    FILE* fu = fortran_get_unit(iunit);
                    fprintf(fu, " title = \"surface y+ data\"\n");
                    if (ialph == 0) {
                        fprintf(fu, " variables = x, z, y+, dn, turb_visc\n");
                    } else {
                        fprintf(fu, " variables = x, y, y+, dn, turb_visc\n");
                    }
                }
            }

        } // end if (ifunc > 0)

    } // end if (myid == myhost) [first block]



    // calculate y+ data
    for (nbl = 1; nbl <= nblock; nbl++) {

        if (myid == mblk2nd(nbl) || myid == myhost) {

            lead_ns::lead(nbl, lw, lw2, maxbl);

            // determine if yplus statistics should be output: block should
            // be turbulent *and* have at least one solid surface
            ivv1 = 0;
            if (ivisc[2] > 1 || ivisc[1] > 1 || ivisc[0] > 1) ivv1 = 1;

            ivv2 = 0;
            for (ns = 1; ns <= nbci0(nbl); ns++) {
                if (std::abs(ibcinfo(nbl,ns,1,1)) == 2004 ||
                    std::abs(ibcinfo(nbl,ns,1,1)) == 2014 ||
                    std::abs(ibcinfo(nbl,ns,1,1)) == 2024 ||
                    std::abs(ibcinfo(nbl,ns,1,1)) == 2034 ||
                    std::abs(ibcinfo(nbl,ns,1,1)) == 2016) ivv2 = ivv2 + 1;
            }
            for (ns = 1; ns <= nbcidim(nbl); ns++) {
                if (std::abs(ibcinfo(nbl,ns,1,2)) == 2004 ||
                    std::abs(ibcinfo(nbl,ns,1,2)) == 2014 ||
                    std::abs(ibcinfo(nbl,ns,1,2)) == 2024 ||
                    std::abs(ibcinfo(nbl,ns,1,2)) == 2034 ||
                    std::abs(ibcinfo(nbl,ns,1,2)) == 2016) ivv2 = ivv2 + 1;
            }
            for (ns = 1; ns <= nbcj0(nbl); ns++) {
                if (std::abs(jbcinfo(nbl,ns,1,1)) == 2004 ||
                    std::abs(jbcinfo(nbl,ns,1,1)) == 2014 ||
                    std::abs(jbcinfo(nbl,ns,1,1)) == 2024 ||
                    std::abs(jbcinfo(nbl,ns,1,1)) == 2034 ||
                    std::abs(jbcinfo(nbl,ns,1,1)) == 2016) ivv2 = ivv2 + 1;
            }
            for (ns = 1; ns <= nbcjdim(nbl); ns++) {
                if (std::abs(jbcinfo(nbl,ns,1,2)) == 2004 ||
                    std::abs(jbcinfo(nbl,ns,1,2)) == 2014 ||
                    std::abs(jbcinfo(nbl,ns,1,2)) == 2024 ||
                    std::abs(jbcinfo(nbl,ns,1,2)) == 2034 ||
                    std::abs(jbcinfo(nbl,ns,1,2)) == 2016) ivv2 = ivv2 + 1;
            }
            for (ns = 1; ns <= nbck0(nbl); ns++) {
                if (std::abs(kbcinfo(nbl,ns,1,1)) == 2004 ||
                    std::abs(kbcinfo(nbl,ns,1,1)) == 2014 ||
                    std::abs(kbcinfo(nbl,ns,1,1)) == 2024 ||
                    std::abs(kbcinfo(nbl,ns,1,1)) == 2034 ||
                    std::abs(kbcinfo(nbl,ns,1,1)) == 2016) ivv2 = ivv2 + 1;
            }
            for (ns = 1; ns <= nbckdim(nbl); ns++) {
                if (std::abs(kbcinfo(nbl,ns,1,2)) == 2004 ||
                    std::abs(kbcinfo(nbl,ns,1,2)) == 2014 ||
                    std::abs(kbcinfo(nbl,ns,1,2)) == 2024 ||
                    std::abs(kbcinfo(nbl,ns,1,2)) == 2034 ||
                    std::abs(kbcinfo(nbl,ns,1,2)) == 2016) ivv2 = ivv2 + 1;
            }

            ivv = 0;
            if (ivv1 > 0 && ivv2 > 0) ivv = 1;

            levg = levelg(nbl);

            // output statistics on global and embedded levels
            if (levg >= lglobal && levg <= levt
                && ncyc1[iseq-1] > 0 && ivv > 0) {

                if (myid == myhost) {
                    fortran_write_unit(11, "\n %s%6d%s%6d%s\n",
                        "YPLUS STATISTICS (endpts not included) - BLOCK",
                        nbl, " (GRID", igridg(nbl), ")");
                }



                if (mblk2nd(nbl) == myid) {
                    if (iadvance(nbl) >= 0) {
                        int ntime_bc = 1;
                        bc_ns::bc(ntime_bc, nbl, lw, lw2, w, mgwk, wk, nwork, cl,
                                  nou, bou, nbuf, ibufdim, maxbl, maxgr, maxseg,
                                  itrans, irotat, idefrm, igridg, nblg,
                                  nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
                                  ibcinfo, jbcinfo, kbcinfo,
                                  bcfilei, bcfilej, bcfilek, lwdat, myid,
                                  idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
                    }
                    ldim = 5;
                    {
                        FortranArray4DRef<double> w_lq(   &w(lq),    jdim, kdim, idim, ldim);
                        FortranArray4DRef<double> w_lqj0( &w(lqj0),  jdim, kdim, 5,    ldim);
                        FortranArray4DRef<double> w_lqk0( &w(lqk0),  jdim, kdim, 5,    ldim);
                        FortranArray4DRef<double> w_lqi0( &w(lqi0),  jdim, kdim, 5,    ldim);
                        FortranArray3DRef<double> w_lbcj( &w(lbcj),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lbck( &w(lbck),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lbci( &w(lbci),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lblk( &w(lblk),  jdim, kdim, idim);
                        qface_ns::qface(jdim, kdim, idim,
                                        w_lq, w_lqj0, w_lqk0, w_lqi0,
                                        w_lbcj, w_lbck, w_lbci, w_lblk, ldim);
                    }
                    if (ivmx >= 2) {
                        ldim = 1;
                        FortranArray4DRef<double> w_lvis( &w(lvis),  jdim, kdim, idim, ldim);
                        FortranArray4DRef<double> w_lvj0( &w(lvj0),  jdim, kdim, 5,    ldim);
                        FortranArray4DRef<double> w_lvk0( &w(lvk0),  jdim, kdim, 5,    ldim);
                        FortranArray4DRef<double> w_lvi0( &w(lvi0),  jdim, kdim, 5,    ldim);
                        FortranArray3DRef<double> w_lbcj( &w(lbcj),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lbck( &w(lbck),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lbci( &w(lbci),  jdim, kdim, idim);
                        FortranArray3DRef<double> w_lblk( &w(lblk),  jdim, kdim, idim);
                        qface_ns::qface(jdim, kdim, idim,
                                        w_lvis, w_lvj0, w_lvk0, w_lvi0,
                                        w_lbcj, w_lbck, w_lbci, w_lblk, ldim);
                    }
                }

                lypj = 1;
                lypk = lypj + kdim*idim*2;
                lypi = lypk + jdim*idim*2;
                lblj_local = lypi + jdim*kdim*2;
                lblk_local = lblj_local + kdim*idim*2;
                lbli = lblk_local + jdim*idim*2;
                ldnj = lbli + jdim*kdim*2;
                ldnk = ldnj + kdim*idim*2;
                ldni = ldnk + jdim*idim*2;
                lvtj = ldni + jdim*kdim*2;
                lvtk = lvtj + kdim*idim*2;
                lvti = lvtk + jdim*idim*2;



                {
                    FortranArray4DRef<double> w_lq(   &w(lq),   jdim, kdim, idim, 5);
                    FortranArray4DRef<double> w_lqi0( &w(lqi0), jdim, kdim, 5,    5);
                    FortranArray4DRef<double> w_lqj0( &w(lqj0), jdim, kdim, 5,    5);
                    FortranArray4DRef<double> w_lqk0( &w(lqk0), jdim, kdim, 5,    5);
                    FortranArray3DRef<double> w_lx(   &w(lx),   jdim, kdim, idim);
                    FortranArray3DRef<double> w_ly(   &w(ly),   jdim, kdim, idim);
                    FortranArray3DRef<double> w_lz(   &w(lz),   jdim, kdim, idim);
                    FortranArray3DRef<double> w_lvis( &w(lvis), jdim, kdim, idim);
                    FortranArray3DRef<double> w_lbcj( &w(lbcj), jdim, kdim, idim);
                    FortranArray3DRef<double> w_lbck( &w(lbck), jdim, kdim, idim);
                    FortranArray3DRef<double> w_lbci( &w(lbci), jdim, kdim, idim);
                    FortranArray4DRef<double> w_lsj(  &w(lsj),  jdim, kdim, idim, 5);
                    FortranArray4DRef<double> w_lsk(  &w(lsk),  jdim, kdim, idim, 5);
                    FortranArray4DRef<double> w_lsi(  &w(lsi),  jdim, kdim, idim, 5);
                    FortranArray3DRef<double> wk_lypj(&wk(lypj), kdim, idim, 2);
                    FortranArray3DRef<double> wk_lypk(&wk(lypk), jdim, idim, 2);
                    FortranArray3DRef<double> wk_lypi(&wk(lypi), jdim, kdim, 2);
                    FortranArray3DRef<double> wk_lblj(&wk(lblj_local), kdim, idim, 2);
                    FortranArray3DRef<double> wk_lblk(&wk(lblk_local), jdim, idim, 2);
                    FortranArray3DRef<double> wk_lbli(&wk(lbli), jdim, kdim, 2);
                    FortranArray3DRef<double> wk_ldnj(&wk(ldnj), kdim, idim, 2);
                    FortranArray3DRef<double> wk_ldnk(&wk(ldnk), jdim, idim, 2);
                    FortranArray3DRef<double> wk_ldni(&wk(ldni), jdim, kdim, 2);
                    FortranArray3DRef<double> wk_lvtj(&wk(lvtj), kdim, idim, 2);
                    FortranArray3DRef<double> wk_lvtk(&wk(lvtk), jdim, idim, 2);
                    FortranArray3DRef<double> wk_lvti(&wk(lvti), jdim, kdim, 2);
                    FortranArray4DRef<double> w_lvj0( &w(lvj0), jdim, kdim, 5,    1);
                    FortranArray4DRef<double> w_lvk0( &w(lvk0), jdim, kdim, 5,    1);
                    FortranArray4DRef<double> w_lvi0( &w(lvi0), jdim, kdim, 5,    1);
                    FortranArray3DRef<double> w_lvol( &w(lvol), jdim, kdim, idim);
                    calyplus_ns::calyplus(jdim, kdim, idim, nbl,
                                          w_lq, w_lqi0, w_lqj0, w_lqk0,
                                          w_lx, w_ly, w_lz, w_lvis,
                                          iovrlp(nbl),
                                          w_lbcj, w_lbck, w_lbci,
                                          w_lsj, w_lsk, w_lsi,
                                          wk_lypj, wk_lypk, wk_lypi,
                                          wk_lblj, wk_lblk, wk_lbli,
                                          wk_ldnj, wk_ldnk, wk_ldni,
                                          wk_lvtj, wk_lvtk, wk_lvti,
                                          ifunc, iunit,
                                          w_lvj0, w_lvk0, w_lvi0,
                                          maxbl, maxseg,
                                          nbci0, nbcj0, nbck0,
                                          nbcidim, nbcjdim, nbckdim,
                                          myid, myhost, mycomm, mblk2nd,
                                          ibcinfo, jbcinfo, kbcinfo,
                                          w_lvol);
                }



                // save yplus statistics on global level for final summary
                if (levg == lglobal) {
                    iflag = 1;
                    nptst    = nptst    + nptsb;
                    ypsumt   = ypsumt   + ypsumb;
                    ypsumsqt = ypsumsqt + ypsumsqb;
                    nypchkt  = nypchkt  + nypchkb;
                    if ((float)ypmaxb > (float)ypmaxt) {
                        ypmaxt   = ypmaxb;
                        iypmaxt  = iypmaxb;
                        jypmaxt  = jypmaxb;
                        kypmaxt  = kypmaxb;
                        nbypmaxt = nbl;
                    }
                    if ((float)ypminb < (float)ypmint) {
                        ypmint   = ypminb;
                        iypmint  = iypminb;
                        jypmint  = jypminb;
                        kypmint  = kypminb;
                        nbypmint = nbl;
                    }
                    if ((float)dnmaxb > (float)dnmaxt) {
                        dnmaxt   = dnmaxb;
                        idnmaxt  = idnmaxb;
                        jdnmaxt  = jdnmaxb;
                        kdnmaxt  = kdnmaxb;
                        nbdnmaxt = nbl;
                    }
                    if ((float)dnminb < (float)dnmint) {
                        dnmint   = dnminb;
                        idnmint  = idnminb;
                        jdnmint  = jdnminb;
                        kdnmint  = kdnminb;
                        nbdnmint = nbl;
                    }
                }

            } // end if levg...

        } // end if (myid == mblk2nd(nbl) || myid == myhost)

    } // end do 10 nbl



    if (myid == myhost) {
        if (iflag > 0) {
            ypavt  = ypsumt / (float)nptst;
            ypstdt = std::sqrt(ypsumsqt / (float)(nptst - 1));
            // 101: YPLUS STATISTICS (endpts not included) - ALL GLOBAL BLOCKS
            fortran_write_unit(11, "\n %s%s\n",
                "YPLUS STATISTICS (endpts not included)",
                " - ALL GLOBAL BLOCKS");
            // 102: header line Y+ MAX
            fortran_write_unit(11, "%s%s\n",
                "     Y+ MAX    ILOC    JLOC    KLOC",
                "     BLOCK      GRID");
            // 107: e10.3, 4x, i4, 4x, i4, 4x, i4, 4x, i6, 4x, i6
            fortran_write_unit(11, " %10.3E    %4d    %4d    %4d    %6d    %6d\n",
                (float)ypmaxt, iypmaxt, jypmaxt, kypmaxt, nbypmaxt, igridg(nbypmaxt));
            // 103: header line Y+ MIN
            fortran_write_unit(11, "%s%s\n",
                "     Y+ MIN    ILOC    JLOC    KLOC",
                "     BLOCK      GRID");
            // 107
            fortran_write_unit(11, " %10.3E    %4d    %4d    %4d    %6d    %6d\n",
                (float)ypmint, iypmint, jypmint, kypmint, nbypmint, igridg(nbypmint));
            // 104: header line DN MAX
            fortran_write_unit(11, "%s%s\n",
                "     DN MAX    ILOC    JLOC    KLOC",
                "     BLOCK      GRID");
            // 107
            fortran_write_unit(11, " %10.3E    %4d    %4d    %4d    %6d    %6d\n",
                (float)dnmaxt, idnmaxt, jdnmaxt, kdnmaxt, nbdnmaxt, igridg(nbdnmaxt));
            // 105: header line DN MIN
            fortran_write_unit(11, "%s%s\n",
                "     DN MIN    ILOC    JLOC    KLOC",
                "     BLOCK      GRID");
            // 107
            fortran_write_unit(11, " %10.3E    %4d    %4d    %4d    %6d    %6d\n",
                (float)dnmint, idnmint, jdnmint, kdnmint, nbdnmint, igridg(nbdnmint));
            // 106: Y+ AVG, Y+ STD DEV, NY+ >, NPTS
            fortran_write_unit(11, "     Y+ AVG    Y+ STD DEV       NY+ >%2d   NPTS\n",
                (int)(float)ypchk);
            // 108: e10.3, 4x, e10.3, 8x, i6, 1x, i6
            fortran_write_unit(11, " %10.3E    %10.3E        %6d %6d\n",
                (float)ypavt, (float)ypstdt, nypchkt, nptst);
            if ((float)ypavt > 2.5f) {
                // 114
                fortran_write_unit(11, "\n %s%s\n%s%s\n%s%s\n%s\n",
                    "WARNING: avg y+ is > 2.5.  It is recommended that",
                    " you revise the grid",
                    "to have smaller min spacing at walls.  (If you",
                    " are employing WALL FUNCTIONS,",
                    "large avg y+ values are acceptable, but wall",
                    " functions are ad hoc",
                    "and not recommended for general use.)");
            }
            if (iwarneddy > 0) {
                // 113
                fortran_write_unit(11, "\n %s%s\n%s %10.3E%s%s\n%s\n",
                    "WARNING: final value of eddy viscosity has been",
                    " limited at one or more",
                    "points to ", (float)edvislim, ".  Unless this run is far from",
                    " convergence, increase",
                    "edvislim via keyword input for subsequent runs");
            }
            if (ifunc > 0) {
                if (i2d == 0) {
                    // 109
                    fortran_write_unit(11, "\n %s\n %s\n %s\n %s\n %s\n",
                        "writing function file surf_y+.fmt (y+)",
                        "writing function file surf_dn.fmt (normal spacing)",
                        "writing function file surf_vist.fmt (turb. visc.)",
                        "surf_xyz.fmt is the corresponding surface grid",
                        "(use multi/formatted/blank as FAST read options)");
                }
                if (i2d == 1) {
                    if (ialph == 0) {
                        // 111
                        fortran_write_unit(11, "\n %s\n",
                            "writing 2d data file surf_y+_2d.fmt (x,z,y+,dn,turb. visc.)");
                    } else {
                        // 112
                        fortran_write_unit(11, "\n %s\n",
                            "writing 2d data file surf_y+_2d.fmt (x,y,y+,dn,turb. visc.)");
                    }
                }
            } else {
                // 110
                fortran_write_unit(11, "\n %s%s\n%s\n",
                    "set ifunc = 1 in subroutine yplusout and recompile",
                    " if function files of",
                    "               y+, dn, and turb visc are desired");
            }
        }
    }

} // end yplusout

} // namespace yplusout_ns
