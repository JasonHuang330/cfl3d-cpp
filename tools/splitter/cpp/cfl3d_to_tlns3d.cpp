// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
#include "cfl3d_to_tlns3d.h"
#include "runtime/fortran_io.h"
#include "shortinp.h"
#include "readkey.h"
#include "termn8.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

namespace cfl3d_to_tlns3d_ns {

void cfl3d_to_tlns3d(int& iver, int& ipatch, int& iptyp, int& ipar, int& nnodes, int& isd, int& mbloc, int& msegn, int& msegt, int& mxbli, FortranArray1DRef<int> il, FortranArray1DRef<int> jl, FortranArray1DRef<int> kl, FortranArray2DRef<double> rkap0g, FortranArray1DRef<int> levelg, FortranArray1DRef<int> igridg, FortranArray2DRef<int> iflimg, FortranArray2DRef<int> ifdsg, FortranArray2DRef<int> iviscg, FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg, FortranArray2DRef<int> idiagg, FortranArray1DRef<int> nblcg, FortranArray2DRef<int> idegg, FortranArray1DRef<int> jsg, FortranArray1DRef<int> ksg, FortranArray1DRef<int> jeg, FortranArray1DRef<int> keg, FortranArray1DRef<int> ieg, FortranArray2DRef<int> mit, FortranArray1DRef<int> ilamlog, FortranArray1DRef<int> ilamhig, FortranArray1DRef<int> jlamlog, FortranArray1DRef<int> jlamhig, FortranArray1DRef<int> klamlog, FortranArray1DRef<int> klamhig, FortranArray2DRef<int> iredundant, FortranArray1DRef<int> nseg1, FortranArray1DRef<int> nseg, FortranArray1DRef<int> nsilo, FortranArray1DRef<int> nsihi, FortranArray1DRef<int> nsjlo, FortranArray1DRef<int> nsjhi, FortranArray1DRef<int> nsklo, FortranArray1DRef<int> nskhi, FortranArray3DRef<int> nb1, FortranArray3DRef<int> ne1, FortranArray3DRef<int> nb2, FortranArray3DRef<int> ne2, FortranArray3DRef<int> ibct, FortranArray4DRef<double> bcval, FortranArray1DRef<int> iovrlp, FortranArray3DRef<int> ifoflg, FortranArray3DRef<int> ndat, FortranArray3DRef<double> xmap, FortranArray3DRef<int> imap, FortranArray2DRef<int> ivisb, FortranArray2DRef<double> twotref, FortranArray1DRef<int> itrb1, FortranArray1DRef<int> itrb2, FortranArray1DRef<int> jtrb1, FortranArray1DRef<int> jtrb2, FortranArray1DRef<int> ktrb1, FortranArray1DRef<int> ktrb2, FortranArray1DRef<int> iturbb, FortranArray1DRef<int> ibif1, FortranArray1DRef<int> ibif2, FortranArray1DRef<int> nbi1, FortranArray1DRef<int> nei1, FortranArray1DRef<int> nbj1, FortranArray1DRef<int> nej1, FortranArray1DRef<int> nbk1, FortranArray1DRef<int> nek1, FortranArray1DRef<int> nbi2, FortranArray1DRef<int> nei2, FortranArray1DRef<int> nbj2, FortranArray1DRef<int> nej2, FortranArray1DRef<int> nbk2, FortranArray1DRef<int> nek2, FortranArray1DRef<int> nd11, FortranArray1DRef<int> nd21, FortranArray1DRef<int> nd12, FortranArray1DRef<int> nd22, FortranArray1DRef<int> ifsor, FortranArray1DRef<int> nb1s, FortranArray1DRef<int> ne1s, FortranArray1DRef<int> nb2s, FortranArray1DRef<int> ne2s, int& nrotat, int& ntrans, double& tlref, double& rlref, int& ioflag, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou, int& nbuf, int& ibufdim, char (&cflout)[80])
{
    // Local variables
    char gridindum[80];
    char sdgridindum[80];
    int iunit5sav;
    int ititr, myid;
    int ibloc, iseg, iface;
    int nn, ib, ns, iii, nnsg;
    int ibl, ifc, isg, nsf;
    int nptpif, iptpif;
    int npatif, ipif;
    int ninter, ncs;
    int ilamflg, iembflg, isum;
    float fsum;
    int ibcfilen;
    int ndata;
    int ib_dummy, is_dummy;
    int idum, iplt3dtyp_local;
    int nface1, nface2;
    int icyc1, icyc2;
    int irev1, irev2;
    int ift1, ift2;
    int idel11, idel12, idel21, idel22;
    int idebug;
    int is, is1, is2;
    int im11, im12, im13, im14, im15, im16, im17, im18, im19, im110, im111, im112;
    int im21, im22, im23, im24, im25, im26, im27, im28, im29, im210, im211, im212;
    int ilamlo, ilamhi, jlamlo, jlamhi, klamlo, klamhi;
    int inewg, igridc, is_emb, js_emb, ks_emb, ie_emb, je_emb, ke_emb;
    int idi, idj, idk, ifi, ifj, ifk;
    int ifdsi, ifdsj, ifdsk;
    float rkapi, rkapj, rkapk;
    int irkapi, irkapj, irkapk;
    int ncg_loc, iem, iadv, ifor;
    int itr, iv;
    int m, i, n, mm;
    int mseq_loc, mgflag_loc, iconsf_loc, mtt_loc, ngam_loc;
    int issc_loc, issr_loc;
    int iipv_loc;
    int ifoo;
    int iif, iftar, nb1t, ne1t, nb2t, ne2t, nfsor_loc;
    int iblt, iijkt, imnmxt, ift;
    int nbeg1t, nend1t, nbeg2t, nend2t;
    int ibls, iijks, imnmxs, ifs;
    int nbeg1s, nend1s, nbeg2s, nend2s;
    int it;
    int itrace;
    int intr;
    int iiseg;
    int ifchk1, ifchk2, ifchk3, ifchk4, ifchk5, ifchk6;
    int ifill;
    int ii;
    int igd;
    int itrans_loc;
    float rfreqt, xmag, ymag, zmag;
    float dxmax, dymax, dzmax;
    int irotat;
    float rfreqr, thxmag, thymag, thzmag;
    float xorig, yorig, zorig;
    float thxmax, thymax, thzmax;
    int idir;
    int nbloc_local;
    // zero-initialize myid (single process context)
    myid = 0;
    int zero_val = 0;
    int neg3 = -3;
    int neg7 = -7;



    // initial some parameters that need to be set if input version differs
    // from output version
    cmn_info3.ialph  = 1;
    cmn_info3.isnd   = 0;
    cmn_info3.ihstry = 0;
    cmn_info2.cfltau = 5.0f;
    cmn_info2.c2spe  = 0.0f;
    nnodes = 1;

    // zero out xmap and ivisb arrays
    for (nn = 1; nn <= msegt; nn++) {
        for (ib = 1; ib <= mbloc; ib++) {
            ivisb(nn, ib) = 0;
            for (ns = 1; ns <= msegn; ns++) {
                xmap(nn, ns, ib) = 0.0;
            }
        }
    }

    // read cfl3d input file (unit 10)
    fortran_read_list(10);
    {
        FILE* f10 = fortran_get_unit(10);
        fgets(gridindum, sizeof(gridindum), f10);
        // strip newline
        int glen = (int)strlen(gridindum);
        if (glen > 0 && gridindum[glen-1] == '\n') gridindum[glen-1] = '\0';
    }
    {
        FILE* f10 = fortran_get_unit(10);
        char linebuf[256];
        // read plt3dg
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.plt3dg, linebuf, 80);
        }
        // read plt3dq
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.plt3dq, linebuf, 80);
        }
        // read output
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.output, linebuf, 80);
        }
        // read resid
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.resid, linebuf, 80);
        }
        // read turbres
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.turbres, linebuf, 80);
        }
        // read blomx
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.blomx, linebuf, 80);
        }
        // read output2
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.output2, linebuf, 80);
        }
        // read printout
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.printout, linebuf, 80);
        }
        // read pplunge
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.pplunge, linebuf, 80);
        }
        // read ovrlap
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.ovrlap, linebuf, 80);
        }
        // read patch
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.patch, linebuf, 80);
        }
        // read restrt
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            int l = (int)strlen(linebuf);
            if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
            strncpy(cmn_cflfiles.restrt, linebuf, 80);
        }
    }

    // check for keyword input
    iunit5sav      = cmn_unit5.iunit5;
    cmn_unit5.iunit5 = 10;
    {
        int iunit11_val = -99;
        int ierrflg_val = -1;
        readkey_ns::readkey(ititr, myid, ibufdim, nbuf, bou, nou, iunit11_val, ierrflg_val);
    }
    cmn_unit5.iunit5 = iunit5sav;

    // read title if not keyword input
    if (ititr == 0) {
        // read(10,10)(title(i),i=1,20)  format(20a4)
        FILE* f10 = fortran_get_unit(10);
        char linebuf[256];
        if (fgets(linebuf, sizeof(linebuf), f10)) {
            // 20a4: read 20 groups of 4 chars into title array (stored as float)
            // title is float[20] in common block, but used as character storage
            // Copy raw bytes into title array
            memset(cmn_info.title, 0, sizeof(cmn_info.title));
            int copylen = (int)strlen(linebuf);
            if (copylen > 80) copylen = 80;
            memcpy(cmn_info.title, linebuf, copylen);
        }
    }

    // read blank line
    fortran_read_list(10);

    // read flow parameters based on version
    if (iver == 4) {
        if (ipar == 0) {
            fortran_read_list(10, &cmn_info.xmach, &cmn_info.alpha, &cmn_info.beta,
                              &cmn_info2.reue, &cmn_info2.tinf, &cmn_info3.isnd, &cmn_info2.c2spe);
        } else {
            fortran_read_list(10, &cmn_info.xmach, &cmn_info.alpha, &cmn_info.beta,
                              &cmn_info2.reue, &cmn_info2.tinf, &cmn_info3.isnd, &cmn_info2.c2spe, &nnodes);
        }
    } else if (iver >= 5) {
        fortran_read_list(10, &cmn_info.xmach, &cmn_info.alpha, &cmn_info.beta,
                          &cmn_info2.reue, &cmn_info2.tinf, &cmn_info3.ialph, &cmn_info3.ihstry);
    }

    // handle Cl-specified option and preconditioning
    cmn_alphait.ialphit = 0;
    cmn_precon1.iprecon = 0;
    if (cmn_info.xmach < 0.0f) {
        if (iver == 4) {
            cmn_alphait.ialphit = 1;
        } else {
            cmn_precon1.iprecon = 1;
        }
        cmn_info.xmach = fabsf(cmn_info.xmach);
    }

    if (iver == 4) {
        fortran_read_list(10);
        if (cmn_alphait.ialphit == 0) {
            fortran_read_list(10, &cmn_info2.sref, &cmn_info2.cref, &cmn_info2.bref,
                              &cmn_info2.xmc, &cmn_info2.ymc, &cmn_info2.zmc);
        } else {
            fortran_read_list(10, &cmn_info2.sref, &cmn_info2.cref, &cmn_info2.bref,
                              &cmn_info2.xmc, &cmn_info2.ymc, &cmn_info2.zmc,
                              &cmn_alphait.cltarg, &cmn_alphait.resupdt);
            cmn_alphait.resupdt = -fabsf(cmn_alphait.resupdt);
        }
        fortran_read_list(10);
        fortran_read_list(10, &cmn_info.dt, &cmn_info3.irest, &cmn_info.iflagts,
                          &cmn_info.fmax, &cmn_info3.iunst, &cmn_info2.rfreq,
                          &cmn_info2.alphau, &cmn_info2.cloc);
    } else {
        fortran_read_list(10);
        fortran_read_list(10, &cmn_info2.sref, &cmn_info2.cref, &cmn_info2.bref,
                          &cmn_info2.xmc, &cmn_info2.ymc, &cmn_info2.zmc);
        fortran_read_list(10);
        fortran_read_list(10, &cmn_info.dt, &cmn_info3.irest, &cmn_info.iflagts,
                          &cmn_info.fmax, &cmn_info3.iunst, &cmn_info2.cfltau);
    }



    // read in number of blocs
    fortran_read_list(10);
    fortran_read_list(10, &cmn_grdinfo.nbloc, &cmn_info3.nplot3d, &cmn_info3.nprint,
                      &cmn_info3.nwrest, &cmn_info3.ichk, &cmn_info3.i2d,
                      &cmn_info3.ntstep, &cmn_info3.ita);
    cmn_info3.ip3dgrd = 0;
    if (cmn_grdinfo.nbloc < 0) cmn_info3.ip3dgrd = 1;
    cmn_grdinfo.nbloc = abs(cmn_grdinfo.nbloc);
    if (cmn_grdinfo.nbloc > mbloc) {
        printf("stopping: parameter mbloc must be at least %d\n", cmn_grdinfo.nbloc);
        termn8_ns::termn8(zero_val, neg3, ibufdim, nbuf, bou, nou);
    }

    // read block viscosity/turbulence data
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        // do 100
        fortran_read_list(10, &ncg_loc, &iem, &iadv, &ifor,
                          &ivisb(1, ibloc), &ivisb(2, ibloc), &ivisb(3, ibloc));
        ivisb(4, ibloc) = ncg_loc;
        ivisb(5, ibloc) = iem;
        ivisb(6, ibloc) = iadv;
        ivisb(7, ibloc) = ifor;

        itr = 0;
        for (iv = 1; iv <= 3; iv++) {
            itr = std::max(itr, abs(ivisb(iv, ibloc)));
        }
        iturbb(ibloc) = 0;
        if (itr >= 2) {
            iturbb(ibloc) = 1;
        }

        // itrb1, etc for tlns3d only - used as zero-range marker in map file
        itrb1(ibloc) = 0;
        itrb2(ibloc) = 0;
        jtrb1(ibloc) = 0;
        jtrb2(ibloc) = 0;
        ktrb1(ibloc) = 0;
        ktrb2(ibloc) = 0;
    } // end do 100

    // read in size of each block
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        // do 110
        fortran_read_list(10, &il(ibloc), &jl(ibloc), &kl(ibloc));
    } // end do 110

    // ilamlo, etc. NOT PRESERVED during splitting...issue a warning
    ilamflg = 0;
    isum    = 0;

    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        // do 120
        fortran_read_list(10, &ilamlo, &ilamhi, &jlamlo, &jlamhi, &klamlo, &klamhi);
        ivisb(8,  ibloc) = ilamlo;
        ivisb(9,  ibloc) = ilamhi;
        ivisb(10, ibloc) = jlamlo;
        ivisb(11, ibloc) = jlamhi;
        ivisb(12, ibloc) = klamlo;
        ivisb(13, ibloc) = klamhi;
        isum = isum + ilamlo + ilamhi + jlamlo + jlamhi + klamlo + klamhi;
    } // end do 120
    if (isum > 0 && strncmp(cflout, "null", 4) != 0) {
        printf("WARNING: laminar regions specified...these are NOT preserved during splitting\n");
        printf("\n");
    }

    // embedded grid data NOT PRESERVED during splitting...issue a warning
    iembflg = 0;
    isum    = 0;

    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        // do 121
        fortran_read_list(10, &inewg, &igridc, &is_emb, &js_emb, &ks_emb,
                          &ie_emb, &je_emb, &ke_emb);
        ivisb(14, ibloc) = inewg;
        ivisb(15, ibloc) = igridc;
        ivisb(16, ibloc) = is_emb;
        ivisb(17, ibloc) = js_emb;
        ivisb(18, ibloc) = ks_emb;
        ivisb(19, ibloc) = ie_emb;
        ivisb(20, ibloc) = je_emb;
        ivisb(21, ibloc) = ke_emb;
        fsum = (float)(isum + inewg + igridc + is_emb + js_emb + ks_emb + ie_emb + je_emb + ke_emb);
    } // end do 121
    if (isum > 0 && strncmp(cflout, "null", 4) != 0) {
        printf("WARNING: embeded grid data specified...these are NOT preserved during splitting\n");
        printf("\n");
    }

    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        // do 122
        fortran_read_list(10, &idi, &idj, &idk, &ifi, &ifj, &ifk);
        ivisb(22, ibloc) = idi;
        ivisb(23, ibloc) = idj;
        ivisb(24, ibloc) = idk;
        ivisb(25, ibloc) = ifi;
        ivisb(26, ibloc) = ifj;
        ivisb(27, ibloc) = ifk;
    } // end do 122

    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        // do 123
        fortran_read_list(10, &ifdsi, &ifdsj, &ifdsk, &rkapi, &rkapj, &rkapk);
        ivisb(28, ibloc) = ifdsi;
        ivisb(29, ibloc) = ifdsj;
        ivisb(30, ibloc) = ifdsk;
        // integer values to represent real input values for rkap
        if (rkapi < -0.9f  && rkapi > -1.1f)  irkapi = -1;
        if (rkapi <  0.1f  && rkapi > -0.1f)  irkapi =  0;
        if (rkapi <  1.1f  && rkapi >  0.9f)  irkapi =  1;
        if (rkapi <  0.40f && rkapi >  0.25f) irkapi =  3;
        if (rkapj < -0.9f  && rkapj > -1.1f)  irkapj = -1;
        if (rkapj <  0.1f  && rkapj > -0.1f)  irkapj =  0;
        if (rkapj <  1.1f  && rkapj >  0.9f)  irkapj =  1;
        if (rkapj <  0.40f && rkapj >  0.25f) irkapj =  3;
        if (rkapk < -0.9f  && rkapk > -1.1f)  irkapk = -1;
        if (rkapk <  0.1f  && rkapk > -0.1f)  irkapk =  0;
        if (rkapk <  1.1f  && rkapk >  0.9f)  irkapk =  1;
        if (rkapk <  0.40f && rkapk >  0.25f) irkapk =  3;
        ivisb(31, ibloc) = irkapi;
        ivisb(32, ibloc) = irkapj;
        ivisb(33, ibloc) = irkapk;
    } // end do 123



    // initialize bcval and ifoflg array
    for (ibloc = 1; ibloc <= mbloc; ibloc++) {
        for (iface = 1; iface <= 6; iface++) {
            for (iseg = 1; iseg <= msegn; iseg++) {
                for (int ll = 1; ll <= 7; ll++) {
                    bcval(ibloc, iface, iseg, ll) = 0.0;
                }
                ifoflg(ibloc, iface, iseg) = 1;
            }
        }
    }

    // initialize bc data file counter
    ibcfilen = 0;

    // read in number of segments on each face of each block
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        // do 130
        fortran_read_list(10, &ifoo, &nsilo(ibloc), &nsihi(ibloc), &nsjlo(ibloc),
                          &nsjhi(ibloc), &nsklo(ibloc), &nskhi(ibloc), &iovrlp(ibloc));
        ivisb(34, ibloc) = iovrlp(ibloc);
    } // end do 130

    // set up non-interface segments for i=imin faces
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        for (iseg = 1; iseg <= nsilo(ibloc); iseg++) {
            // do 140
            iface = 1;
            fortran_read_list(10, &ib_dummy, &is_dummy, &ibct(ibloc, iface, iseg),
                              &nb1(ibloc, iface, iseg), &ne1(ibloc, iface, iseg),
                              &nb2(ibloc, iface, iseg), &ne2(ibloc, iface, iseg), &ndata);
            shortinp_ns::shortinp(nb1(ibloc, iface, iseg), ne1(ibloc, iface, iseg),
                                  nb2(ibloc, iface, iseg), ne2(ibloc, iface, iseg),
                                  il(ibloc), jl(ibloc), kl(ibloc), iface);
            if (ibct(ibloc, iface, iseg) == 2005 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: periodic bc data specified...these are NOT preserved during splitting\n");
            } else if (ibct(ibloc, iface, iseg) == 2006 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: radial equilibrium bc data specified...these MAY NOT be preserved during splitting\n");
            }
            ndat(ibloc, iface, iseg) = ndata;
            ndat(ibloc, iface, iseg) = ndata;
            if (abs(ndata) > 0) {
                fortran_read_list(10);
                if (ndata > 0) {
                    // implied-DO: read all ndata values from ONE line
                    FILE* f10 = fortran_get_unit(10);
                    char linebuf[256];
                    if (fgets(linebuf, sizeof(linebuf), f10)) {
                        int pos = 0, tmp_pos = 0;
                        for (int ll = 1; ll <= ndata; ll++) {
                            double tmp_val = 0.0;
                            sscanf(linebuf + pos, " %lf%n", &tmp_val, &tmp_pos);
                            bcval(ibloc, iface, iseg, ll) = tmp_val;
                            pos += tmp_pos;
                        }
                    }
                } else {
                    // can't actually deal with file name yet, so just assign a number
                    fortran_read_list(10);
                    ibcfilen = ibcfilen + 1;
                    bcval(ibloc, iface, iseg, 1) = (double)ibcfilen;
                }
            }
            if (is_dummy < 0) ifoflg(ibloc, iface, iseg) = 0;
        }
    } // end do 140

    // set up non-interface segments for i=imax faces
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        for (iseg = 1; iseg <= nsihi(ibloc); iseg++) {
            // do 145
            iface = 2;
            fortran_read_list(10, &ib_dummy, &is_dummy, &ibct(ibloc, iface, iseg),
                              &nb1(ibloc, iface, iseg), &ne1(ibloc, iface, iseg),
                              &nb2(ibloc, iface, iseg), &ne2(ibloc, iface, iseg), &ndata);
            shortinp_ns::shortinp(nb1(ibloc, iface, iseg), ne1(ibloc, iface, iseg),
                                  nb2(ibloc, iface, iseg), ne2(ibloc, iface, iseg),
                                  il(ibloc), jl(ibloc), kl(ibloc), iface);
            if (ibct(ibloc, iface, iseg) == 2005 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: periodic bc data specified...these are NOT preserved during splitting\n");
            } else if (ibct(ibloc, iface, iseg) == 2006 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: radial equilibrium bc data specified...these MAY NOT be preserved during splitting\n");
            }
            ndat(ibloc, iface, iseg) = ndata;
            if (abs(ndata) > 0) {
                fortran_read_list(10);
                if (ndata > 0) {
                    FILE* f10 = fortran_get_unit(10);
                    char linebuf[256];
                    if (fgets(linebuf, sizeof(linebuf), f10)) {
                        int pos = 0, tmp_pos = 0;
                        for (int ll = 1; ll <= ndata; ll++) {
                            double tmp_val = 0.0;
                            sscanf(linebuf + pos, " %lf%n", &tmp_val, &tmp_pos);
                            bcval(ibloc, iface, iseg, ll) = tmp_val;
                            pos += tmp_pos;
                        }
                    }
                } else {
                    fortran_read_list(10);
                    ibcfilen = ibcfilen + 1;
                    bcval(ibloc, iface, iseg, 1) = (double)ibcfilen;
                }
            }
            if (is_dummy < 0) ifoflg(ibloc, iface, iseg) = 0;
        }
    } // end do 145

    // set up non-interface segments for j=jmin faces
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        for (iseg = 1; iseg <= nsjlo(ibloc); iseg++) {
            // do 150 - note: nb2/ne2 read before nb1/ne1 for j-faces
            iface = 3;
            fortran_read_list(10, &ib_dummy, &is_dummy, &ibct(ibloc, iface, iseg),
                              &nb2(ibloc, iface, iseg), &ne2(ibloc, iface, iseg),
                              &nb1(ibloc, iface, iseg), &ne1(ibloc, iface, iseg), &ndata);
            shortinp_ns::shortinp(nb1(ibloc, iface, iseg), ne1(ibloc, iface, iseg),
                                  nb2(ibloc, iface, iseg), ne2(ibloc, iface, iseg),
                                  il(ibloc), jl(ibloc), kl(ibloc), iface);
            if (ibct(ibloc, iface, iseg) == 2005 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: periodic bc data specified...these are NOT preserved during splitting\n");
            } else if (ibct(ibloc, iface, iseg) == 2006 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: radial equilibrium bc data specified...these MAY NOT be preserved during splitting\n");
            }
            ndat(ibloc, iface, iseg) = ndata;
            if (abs(ndata) > 0) {
                fortran_read_list(10);
                if (ndata > 0) {
                    FILE* f10 = fortran_get_unit(10);
                    char linebuf[256];
                    if (fgets(linebuf, sizeof(linebuf), f10)) {
                        int pos = 0, tmp_pos = 0;
                        for (int ll = 1; ll <= ndata; ll++) {
                            double tmp_val = 0.0;
                            sscanf(linebuf + pos, " %lf%n", &tmp_val, &tmp_pos);
                            bcval(ibloc, iface, iseg, ll) = tmp_val;
                            pos += tmp_pos;
                        }
                    }
                } else {
                    ibcfilen = ibcfilen + 1;
                    fortran_read_list(10);
                    bcval(ibloc, iface, iseg, 1) = (double)ibcfilen;
                }
            }
            if (is_dummy < 0) ifoflg(ibloc, iface, iseg) = 0;
        }
    } // end do 150



    // set up non-interface segments for j=jmax faces
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        for (iseg = 1; iseg <= nsjhi(ibloc); iseg++) {
            // do 155 - note: nb2/ne2 read before nb1/ne1 for j-faces
            iface = 4;
            fortran_read_list(10, &ib_dummy, &is_dummy, &ibct(ibloc, iface, iseg),
                              &nb2(ibloc, iface, iseg), &ne2(ibloc, iface, iseg),
                              &nb1(ibloc, iface, iseg), &ne1(ibloc, iface, iseg), &ndata);
            shortinp_ns::shortinp(nb1(ibloc, iface, iseg), ne1(ibloc, iface, iseg),
                                  nb2(ibloc, iface, iseg), ne2(ibloc, iface, iseg),
                                  il(ibloc), jl(ibloc), kl(ibloc), iface);
            if (ibct(ibloc, iface, iseg) == 2005 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: periodic bc data specified...these are NOT preserved during splitting\n");
            } else if (ibct(ibloc, iface, iseg) == 2006 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: radial equilibrium bc data specified...these MAY NOT be preserved during splitting\n");
            }
            ndat(ibloc, iface, iseg) = ndata;
            if (abs(ndata) > 0) {
                fortran_read_list(10);
                if (ndata > 0) {
                    FILE* f10 = fortran_get_unit(10);
                    char linebuf[256];
                    if (fgets(linebuf, sizeof(linebuf), f10)) {
                        int pos = 0, tmp_pos = 0;
                        for (int ll = 1; ll <= ndata; ll++) {
                            double tmp_val = 0.0;
                            sscanf(linebuf + pos, " %lf%n", &tmp_val, &tmp_pos);
                            bcval(ibloc, iface, iseg, ll) = tmp_val;
                            pos += tmp_pos;
                        }
                    }
                } else {
                    fortran_read_list(10);
                    ibcfilen = ibcfilen + 1;
                    bcval(ibloc, iface, iseg, 1) = (double)ibcfilen;
                }
            }
            if (is_dummy < 0) ifoflg(ibloc, iface, iseg) = 0;
        }
    } // end do 155

    // set up non-interface segments for k=kmin faces
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        for (iseg = 1; iseg <= nsklo(ibloc); iseg++) {
            // do 160
            iface = 5;
            fortran_read_list(10, &ib_dummy, &is_dummy, &ibct(ibloc, iface, iseg),
                              &nb1(ibloc, iface, iseg), &ne1(ibloc, iface, iseg),
                              &nb2(ibloc, iface, iseg), &ne2(ibloc, iface, iseg), &ndata);
            shortinp_ns::shortinp(nb1(ibloc, iface, iseg), ne1(ibloc, iface, iseg),
                                  nb2(ibloc, iface, iseg), ne2(ibloc, iface, iseg),
                                  il(ibloc), jl(ibloc), kl(ibloc), iface);
            if (ibct(ibloc, iface, iseg) == 2005 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: periodic bc data specified...these are NOT preserved during splitting\n");
            } else if (ibct(ibloc, iface, iseg) == 2006 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: radial equilibrium bc data specified...these MAY NOT be preserved during splitting\n");
            }
            ndat(ibloc, iface, iseg) = ndata;
            if (abs(ndata) > 0) {
                fortran_read_list(10);
                if (ndata > 0) {
                    FILE* f10 = fortran_get_unit(10);
                    char linebuf[256];
                    if (fgets(linebuf, sizeof(linebuf), f10)) {
                        int pos = 0, tmp_pos = 0;
                        for (int ll = 1; ll <= ndata; ll++) {
                            double tmp_val = 0.0;
                            sscanf(linebuf + pos, " %lf%n", &tmp_val, &tmp_pos);
                            bcval(ibloc, iface, iseg, ll) = tmp_val;
                            pos += tmp_pos;
                        }
                    }
                } else {
                    fortran_read_list(10);
                    ibcfilen = ibcfilen + 1;
                    bcval(ibloc, iface, iseg, 1) = (double)ibcfilen;
                }
            }
            if (is_dummy < 0) ifoflg(ibloc, iface, iseg) = 0;
        }
    } // end do 160

    // set up non-interface segments for k=kmax faces
    fortran_read_list(10);
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        for (iseg = 1; iseg <= nskhi(ibloc); iseg++) {
            // do 165
            iface = 6;
            fortran_read_list(10, &ib_dummy, &is_dummy, &ibct(ibloc, iface, iseg),
                              &nb1(ibloc, iface, iseg), &ne1(ibloc, iface, iseg),
                              &nb2(ibloc, iface, iseg), &ne2(ibloc, iface, iseg), &ndata);
            shortinp_ns::shortinp(nb1(ibloc, iface, iseg), ne1(ibloc, iface, iseg),
                                  nb2(ibloc, iface, iseg), ne2(ibloc, iface, iseg),
                                  il(ibloc), jl(ibloc), kl(ibloc), iface);
            if (ibct(ibloc, iface, iseg) == 2005 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: periodic bc data specified...these are NOT preserved during splitting\n");
            } else if (ibct(ibloc, iface, iseg) == 2006 && strncmp(cflout, "null", 4) != 0) {
                printf("WARNING: radial equilibrium bc data specified...these MAY NOT be preserved during splitting\n");
            }
            ndat(ibloc, iface, iseg) = ndata;
            if (abs(ndata) > 0) {
                fortran_read_list(10);
                if (ndata > 0) {
                    FILE* f10 = fortran_get_unit(10);
                    char linebuf[256];
                    if (fgets(linebuf, sizeof(linebuf), f10)) {
                        int pos = 0, tmp_pos = 0;
                        for (int ll = 1; ll <= ndata; ll++) {
                            double tmp_val = 0.0;
                            sscanf(linebuf + pos, " %lf%n", &tmp_val, &tmp_pos);
                            bcval(ibloc, iface, iseg, ll) = tmp_val;
                            pos += tmp_pos;
                        }
                    }
                } else {
                    fortran_read_list(10);
                    ibcfilen = ibcfilen + 1;
                    bcval(ibloc, iface, iseg, 1) = (double)ibcfilen;
                }
            }
            if (is_dummy < 0) ifoflg(ibloc, iface, iseg) = 0;
        }
    } // end do 165

    // initialize imap array
    for (iii = 1; iii <= msegt; iii++) {
        for (nnsg = 1; nnsg <= msegn; nnsg++) {
            for (ibl = 1; ibl <= mbloc; ibl++) {
                imap(iii, nnsg, ibl) = 0;
            }
        }
    }

    // set up imap array using the info that we have so far
    for (ibl = 1; ibl <= cmn_grdinfo.nbloc; ibl++) {
        // do 200
        nnsg = 0;
        for (ifc = 1; ifc <= 6; ifc++) {
            // do 201
            if (ifc == 1) nsf = nsilo(ibl);
            if (ifc == 2) nsf = nsihi(ibl);
            if (ifc == 3) nsf = nsjlo(ibl);
            if (ifc == 4) nsf = nsjhi(ibl);
            if (ifc == 5) nsf = nsklo(ibl);
            if (ifc == 6) nsf = nskhi(ibl);
            for (isg = 1; isg <= nsf; isg++) {
                // do 202
                if (ibct(ibl, ifc, isg) != 0) {
                    nnsg      = nnsg + 1;
                    nseg(ibl) = nnsg;
                    if (nseg(ibl) > msegn) {
                        printf("stopping: parameter msegn must be at least %d\n", nseg(ibl));
                        termn8_ns::termn8(zero_val, neg7, ibufdim, nbuf, bou, nou);
                    }

                    imap(1,  nnsg, ibl) = ibct(ibl, ifc, isg);
                    imap(2,  nnsg, ibl) = ifc;
                    imap(3,  nnsg, ibl) = nb1(ibl, ifc, isg);
                    imap(4,  nnsg, ibl) = ne1(ibl, ifc, isg);
                    imap(5,  nnsg, ibl) = nb2(ibl, ifc, isg);
                    imap(6,  nnsg, ibl) = ne2(ibl, ifc, isg);
                    imap(7,  nnsg, ibl) = 0;
                    imap(8,  nnsg, ibl) = 0;
                    if (imap(1, nnsg, ibl) == 4) {
                        idir = 2;
                        if (ifc <= 2) idir = 1;
                        if (ifc >= 5) idir = 3;
                        imap(8, nnsg, ibl) = idir;
                    }
                    imap(9,  nnsg, ibl) = 0;
                    imap(10, nnsg, ibl) = 0;
                    imap(11, nnsg, ibl) = 0;
                    imap(12, nnsg, ibl) = 0;

                    imap(21, nnsg, ibl) = ndat(ibl, ifc, isg);
                    imap(22, nnsg, ibl) = ifoflg(ibl, ifc, isg);
                    imap(23, nnsg, ibl) = iovrlp(ibl);
                    xmap(1,  nnsg, ibl) = bcval(ibl, ifc, isg, 1);
                    xmap(2,  nnsg, ibl) = bcval(ibl, ifc, isg, 2);
                    xmap(3,  nnsg, ibl) = bcval(ibl, ifc, isg, 3);
                    xmap(4,  nnsg, ibl) = bcval(ibl, ifc, isg, 4);
                    xmap(5,  nnsg, ibl) = bcval(ibl, ifc, isg, 5);
                    xmap(6,  nnsg, ibl) = bcval(ibl, ifc, isg, 6);
                    xmap(7,  nnsg, ibl) = bcval(ibl, ifc, isg, 7);
                }
            } // end do 202
        } // end do 201
    } // end do 200



    // read multigrid/sequencing data
    fortran_read_list(10);
    if (iver == 4) {
        fortran_read_list(10, &cmn_info.mseq, &cmn_info.mgflag, &cmn_info.iconsf,
                          &cmn_info3.mtt, &cmn_info.ngam);
    } else {
        if (cmn_precon1.iprecon == 0) {
            fortran_read_list(10, &cmn_info.mseq, &cmn_info.mgflag, &cmn_info.iconsf,
                              &cmn_info3.mtt, &cmn_info.ngam);
        } else {
            fortran_read_list(10, &cmn_info.mseq, &cmn_info.mgflag, &cmn_info.iconsf,
                              &cmn_info3.mtt, &cmn_info.ngam,
                              &cmn_precon1.cprec, &cmn_precon1.uref, &cmn_precon1.avn);
        }
    }
    fortran_read_list(10);
    fortran_read_list(10, &cmn_mgv.issc,
                      &cmn_mgv.epsssc[0], &cmn_mgv.epsssc[1], &cmn_mgv.epsssc[2],
                      &cmn_mgv.issr,
                      &cmn_mgv.epsssr[0], &cmn_mgv.epsssr[1], &cmn_mgv.epsssr[2]);
    fortran_read_list(10);
    for (m = 1; m <= cmn_info.mseq; m++) {
        fortran_read_list(10, &cmn_info.ncyc1[m-1], &cmn_info3.mglev[m-1],
                          &cmn_info3.nem[m-1], &cmn_info.nitfo1[m-1]);
    }
    fortran_read_list(10);
    for (m = 1; m <= cmn_info.mseq; m++) {
        int nmitl = cmn_info3.mglev[m-1] + cmn_info3.nem[m-1];
        // read (mitl(i,m),i=1,(mglev(m)+nem(m))) as list-directed on one line
        // mitl is stored column-major: mitl(i,m) -> mitl[(m-1)*5 + (i-1)]
        // This is an implied-DO: all values on one line
        {
            FILE* f10 = fortran_get_unit(10);
            char linebuf[256];
            if (fgets(linebuf, sizeof(linebuf), f10)) {
                int pos = 0, tmp_pos = 0;
                for (i = 1; i <= nmitl; i++) {
                    sscanf(linebuf + pos, " %d%n", &cmn_info3.mitl[(m-1)*5 + (i-1)], &tmp_pos);
                    pos += tmp_pos;
                }
            }
        }
    }
    fortran_read_list(10);
    fortran_read_list(10);

    // set up point-to-point interface segments
    fortran_read_list(10, &nptpif);

    if (iver == 4) {
        fortran_read_list(10);
        for (iptpif = 1; iptpif <= nptpif; iptpif++) {
            fortran_read_list(10, &ifoo, &ibif1(iptpif), &ibif2(iptpif));
        }
        fortran_read_list(10);
        for (iptpif = 1; iptpif <= nptpif; iptpif++) {
            fortran_read_list(10, &nbi1(iptpif), &nbj1(iptpif), &nbk1(iptpif),
                              &nei1(iptpif), &nej1(iptpif), &nek1(iptpif),
                              &nd11(iptpif), &nd21(iptpif));
        }
        fortran_read_list(10);
        for (iptpif = 1; iptpif <= nptpif; iptpif++) {
            fortran_read_list(10, &nbi2(iptpif), &nbj2(iptpif), &nbk2(iptpif),
                              &nei2(iptpif), &nej2(iptpif), &nek2(iptpif),
                              &nd12(iptpif), &nd22(iptpif));
        }
    } else {
        fortran_read_list(10);
        for (iptpif = 1; iptpif <= nptpif; iptpif++) {
            fortran_read_list(10, &ifoo, &ibif1(iptpif),
                              &nbi1(iptpif), &nbj1(iptpif), &nbk1(iptpif),
                              &nei1(iptpif), &nej1(iptpif), &nek1(iptpif),
                              &nd11(iptpif), &nd21(iptpif));
        }
        fortran_read_list(10);
        for (iptpif = 1; iptpif <= nptpif; iptpif++) {
            fortran_read_list(10, &ifoo, &ibif2(iptpif),
                              &nbi2(iptpif), &nbj2(iptpif), &nbk2(iptpif),
                              &nei2(iptpif), &nej2(iptpif), &nek2(iptpif),
                              &nd12(iptpif), &nd22(iptpif));
        }
    }



    // do 230: process point-to-point interfaces
    for (iptpif = 1; iptpif <= nptpif; iptpif++) {

        // determine face type for target and source interfaces
        if (nbi1(iptpif) == nei1(iptpif)) {
            nface1 = 1;
            if (nbi1(iptpif) != 1) nface1 = 2;
        } else if (nbj1(iptpif) == nej1(iptpif)) {
            nface1 = 3;
            if (nbj1(iptpif) != 1) nface1 = 4;
        } else {
            nface1 = 5;
            if (nbk1(iptpif) != 1) nface1 = 6;
        }

        if (nbi2(iptpif) == nei2(iptpif)) {
            nface2 = 1;
            if (nbi2(iptpif) != 1) nface2 = 2;
        } else if (nbj2(iptpif) == nej2(iptpif)) {
            nface2 = 3;
            if (nbj2(iptpif) != 1) nface2 = 4;
        } else {
            nface2 = 5;
            if (nbk2(iptpif) != 1) nface2 = 6;
        }

        // block 1 side (1=target, 2=source)
        ibl       = ibif1(iptpif);
        nseg(ibl) = nseg(ibl) + 1;
        if (nseg(ibl) > msegn) {
            printf("stopping: parameter msegn must be at least %d\n", nseg(ibl));
            termn8_ns::termn8(zero_val, neg7, ibufdim, nbuf, bou, nou);
        }
        isg = nseg(ibl);

        imap(1,  isg, ibl) = 1;
        imap(2,  isg, ibl) = nface1;
        imap(7,  isg, ibl) = ibif2(iptpif);
        imap(8,  isg, ibl) = nface2;
        imap(22, isg, ibl) = 1;
        imap(23, isg, ibl) = iovrlp(ibl);

        // determine whether nfaces should be negative
        icyc1 = 1;
        icyc2 = 1;

        if (nface1 <= 2) {
            if (nd11(iptpif) != 2) icyc1 = -1;
        } else if (nface1 >= 5) {
            if (nd11(iptpif) != 1) icyc1 = -1;
        } else {
            if (nd11(iptpif) != 3) icyc1 = -1;
        }

        if (nface2 <= 2) {
            if (nd12(iptpif) != 2) icyc2 = -1;
        } else if (nface2 >= 5) {
            if (nd12(iptpif) != 1) icyc2 = -1;
        } else {
            if (nd12(iptpif) != 3) icyc2 = -1;
        }

        if (icyc1 * icyc2 < 0) imap(8, isg, ibl) = -imap(8, isg, ibl);

        irev1 = 0;
        irev2 = 0;

        // i-face (target)
        if (nface1 <= 2) {
            imap(3, isg, ibl) = nbj1(iptpif);
            imap(4, isg, ibl) = nej1(iptpif);
            imap(5, isg, ibl) = nbk1(iptpif);
            imap(6, isg, ibl) = nek1(iptpif);
            if (nej1(iptpif) < nbj1(iptpif)) {
                irev1 = 1;
                imap(3, isg, ibl) = nej1(iptpif);
                imap(4, isg, ibl) = nbj1(iptpif);
            }
            if (nek1(iptpif) < nbk1(iptpif)) {
                irev2 = 1;
                imap(5, isg, ibl) = nek1(iptpif);
                imap(6, isg, ibl) = nbk1(iptpif);
            }
        // k-face (target)
        } else if (nface1 >= 5) {
            imap(3, isg, ibl) = nbi1(iptpif);
            imap(4, isg, ibl) = nei1(iptpif);
            imap(5, isg, ibl) = nbj1(iptpif);
            imap(6, isg, ibl) = nej1(iptpif);
            if (nei1(iptpif) < nbi1(iptpif)) {
                irev1 = 1;
                imap(3, isg, ibl) = nei1(iptpif);
                imap(4, isg, ibl) = nbi1(iptpif);
            }
            if (nej1(iptpif) < nbj1(iptpif)) {
                irev2 = 1;
                imap(5, isg, ibl) = nej1(iptpif);
                imap(6, isg, ibl) = nbj1(iptpif);
            }
        // j-face (target)
        } else {
            imap(3, isg, ibl) = nbk1(iptpif);
            imap(4, isg, ibl) = nek1(iptpif);
            imap(5, isg, ibl) = nbi1(iptpif);
            imap(6, isg, ibl) = nei1(iptpif);
            if (nek1(iptpif) < nbk1(iptpif)) {
                irev1 = 1;
                imap(3, isg, ibl) = nek1(iptpif);
                imap(4, isg, ibl) = nbk1(iptpif);
            }
            if (nei1(iptpif) < nbi1(iptpif)) {
                irev2 = 1;
                imap(5, isg, ibl) = nei1(iptpif);
                imap(6, isg, ibl) = nbi1(iptpif);
            }
        }



        // i-face (source) for block 1 side
        if (nface2 <= 2) {
            imap(9,  isg, ibl) = nbj2(iptpif);
            imap(10, isg, ibl) = nej2(iptpif);
            imap(11, isg, ibl) = nbk2(iptpif);
            imap(12, isg, ibl) = nek2(iptpif);
            if (irev1 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(9,  isg, ibl) = nej2(iptpif);
                    imap(10, isg, ibl) = nbj2(iptpif);
                } else {
                    imap(11, isg, ibl) = nek2(iptpif);
                    imap(12, isg, ibl) = nbk2(iptpif);
                }
            }
            if (irev2 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(11, isg, ibl) = nek2(iptpif);
                    imap(12, isg, ibl) = nbk2(iptpif);
                } else {
                    imap(9,  isg, ibl) = nej2(iptpif);
                    imap(10, isg, ibl) = nbj2(iptpif);
                }
            }
        // k-face (source) for block 1 side
        } else if (nface2 >= 5) {
            imap(9,  isg, ibl) = nbi2(iptpif);
            imap(10, isg, ibl) = nei2(iptpif);
            imap(11, isg, ibl) = nbj2(iptpif);
            imap(12, isg, ibl) = nej2(iptpif);
            if (irev1 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(9,  isg, ibl) = nei2(iptpif);
                    imap(10, isg, ibl) = nbi2(iptpif);
                } else {
                    imap(11, isg, ibl) = nej2(iptpif);
                    imap(12, isg, ibl) = nbj2(iptpif);
                }
            }
            if (irev2 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(11, isg, ibl) = nej2(iptpif);
                    imap(12, isg, ibl) = nbj2(iptpif);
                } else {
                    imap(9,  isg, ibl) = nei2(iptpif);
                    imap(10, isg, ibl) = nbi2(iptpif);
                }
            }
        // j-face (source) for block 1 side
        } else {
            imap(9,  isg, ibl) = nbk2(iptpif);
            imap(10, isg, ibl) = nek2(iptpif);
            imap(11, isg, ibl) = nbi2(iptpif);
            imap(12, isg, ibl) = nei2(iptpif);
            if (irev1 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(9,  isg, ibl) = nek2(iptpif);
                    imap(10, isg, ibl) = nbk2(iptpif);
                } else {
                    imap(11, isg, ibl) = nei2(iptpif);
                    imap(12, isg, ibl) = nbi2(iptpif);
                }
            }
            if (irev2 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(11, isg, ibl) = nei2(iptpif);
                    imap(12, isg, ibl) = nbi2(iptpif);
                } else {
                    imap(9,  isg, ibl) = nek2(iptpif);
                    imap(10, isg, ibl) = nbk2(iptpif);
                }
            }
        }

        // set imap(1,isg,ibl) = 0 for 1-1 segments that map to themselves
        if (ibif1(iptpif) == ibif2(iptpif)) {
            ift1   = abs(imap(2,  isg, ibl));
            ift2   = abs(imap(8,  isg, ibl));
            idel11 = abs(imap(3,  isg, ibl) - imap(4,  isg, ibl));
            idel12 = abs(imap(5,  isg, ibl) - imap(6,  isg, ibl));
            idel21 = abs(imap(9,  isg, ibl) - imap(10, isg, ibl));
            idel22 = abs(imap(11, isg, ibl) - imap(12, isg, ibl));
            if ((ift1 == ift2) && (idel11 == idel21) && (idel12 == idel22))
                imap(1, isg, ibl) = 0;
        }



        // block 2 side (1=source, 2=target)
        ibl       = ibif2(iptpif);
        nseg(ibl) = nseg(ibl) + 1;
        if (nseg(ibl) > msegn) {
            printf("stopping: parameter msegn must be at least %d\n", nseg(ibl));
            termn8_ns::termn8(zero_val, neg7, ibufdim, nbuf, bou, nou);
        }
        isg = nseg(ibl);

        imap(1,  isg, ibl) = 1;
        imap(2,  isg, ibl) = nface2;
        imap(7,  isg, ibl) = ibif1(iptpif);
        imap(8,  isg, ibl) = nface1;
        imap(22, isg, ibl) = 1;
        imap(23, isg, ibl) = iovrlp(ibl);

        // determine whether nfaces should be negative (block 2 side)
        icyc1 = 1;
        icyc2 = 1;

        if (nface2 <= 2) {
            if (nd12(iptpif) != 2) icyc2 = -1;
        } else if (nface2 >= 5) {
            if (nd12(iptpif) != 1) icyc2 = -1;
        } else {
            if (nd12(iptpif) != 3) icyc2 = -1;
        }

        if (nface1 <= 2) {
            if (nd11(iptpif) != 2) icyc1 = -1;
        } else if (nface1 >= 5) {
            if (nd11(iptpif) != 1) icyc1 = -1;
        } else {
            if (nd11(iptpif) != 3) icyc1 = -1;
        }

        if (icyc1 * icyc2 < 0) imap(8, isg, ibl) = -imap(8, isg, ibl);

        irev1 = 0;
        irev2 = 0;

        // i-face (target) for block 2 side
        if (nface2 <= 2) {
            imap(3, isg, ibl) = nbj2(iptpif);
            imap(4, isg, ibl) = nej2(iptpif);
            imap(5, isg, ibl) = nbk2(iptpif);
            imap(6, isg, ibl) = nek2(iptpif);
            if (nej2(iptpif) < nbj2(iptpif)) {
                irev1 = 1;
                imap(3, isg, ibl) = nej2(iptpif);
                imap(4, isg, ibl) = nbj2(iptpif);
            }
            if (nek2(iptpif) < nbk2(iptpif)) {
                irev2 = 1;
                imap(5, isg, ibl) = nek2(iptpif);
                imap(6, isg, ibl) = nbk2(iptpif);
            }
        // k-face (target) for block 2 side
        } else if (nface2 >= 5) {
            imap(3, isg, ibl) = nbi2(iptpif);
            imap(4, isg, ibl) = nei2(iptpif);
            imap(5, isg, ibl) = nbj2(iptpif);
            imap(6, isg, ibl) = nej2(iptpif);
            if (nei2(iptpif) < nbi2(iptpif)) {
                irev1 = 1;
                imap(3, isg, ibl) = nei2(iptpif);
                imap(4, isg, ibl) = nbi2(iptpif);
            }
            if (nej2(iptpif) < nbj2(iptpif)) {
                irev2 = 1;
                imap(5, isg, ibl) = nej2(iptpif);
                imap(6, isg, ibl) = nbj2(iptpif);
            }
        // j-face (target) for block 2 side
        } else {
            imap(3, isg, ibl) = nbk2(iptpif);
            imap(4, isg, ibl) = nek2(iptpif);
            imap(5, isg, ibl) = nbi2(iptpif);
            imap(6, isg, ibl) = nei2(iptpif);
            if (nek2(iptpif) < nbk2(iptpif)) {
                irev1 = 1;
                imap(3, isg, ibl) = nek2(iptpif);
                imap(4, isg, ibl) = nbk2(iptpif);
            }
            if (nei2(iptpif) < nbi2(iptpif)) {
                irev2 = 1;
                imap(5, isg, ibl) = nei2(iptpif);
                imap(6, isg, ibl) = nbi2(iptpif);
            }
        }



        // i-face (source) for block 2 side
        if (nface1 <= 2) {
            imap(9,  isg, ibl) = nbj1(iptpif);
            imap(10, isg, ibl) = nej1(iptpif);
            imap(11, isg, ibl) = nbk1(iptpif);
            imap(12, isg, ibl) = nek1(iptpif);
            if (irev1 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(9,  isg, ibl) = nej1(iptpif);
                    imap(10, isg, ibl) = nbj1(iptpif);
                } else {
                    imap(11, isg, ibl) = nek1(iptpif);
                    imap(12, isg, ibl) = nbk1(iptpif);
                }
            }
            if (irev2 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(11, isg, ibl) = nek1(iptpif);
                    imap(12, isg, ibl) = nbk1(iptpif);
                } else {
                    imap(9,  isg, ibl) = nej1(iptpif);
                    imap(10, isg, ibl) = nbj1(iptpif);
                }
            }
        // k-face (source) for block 2 side
        } else if (nface1 >= 5) {
            imap(9,  isg, ibl) = nbi1(iptpif);
            imap(10, isg, ibl) = nei1(iptpif);
            imap(11, isg, ibl) = nbj1(iptpif);
            imap(12, isg, ibl) = nej1(iptpif);
            if (irev1 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(9,  isg, ibl) = nei1(iptpif);
                    imap(10, isg, ibl) = nbi1(iptpif);
                } else {
                    imap(11, isg, ibl) = nej1(iptpif);
                    imap(12, isg, ibl) = nbj1(iptpif);
                }
            }
            if (irev2 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(11, isg, ibl) = nej1(iptpif);
                    imap(12, isg, ibl) = nbj1(iptpif);
                } else {
                    imap(9,  isg, ibl) = nei1(iptpif);
                    imap(10, isg, ibl) = nbi1(iptpif);
                }
            }
        // j-face (source) for block 2 side
        } else {
            imap(9,  isg, ibl) = nbk1(iptpif);
            imap(10, isg, ibl) = nek1(iptpif);
            imap(11, isg, ibl) = nbi1(iptpif);
            imap(12, isg, ibl) = nei1(iptpif);
            if (irev1 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(9,  isg, ibl) = nek1(iptpif);
                    imap(10, isg, ibl) = nbk1(iptpif);
                } else {
                    imap(11, isg, ibl) = nei1(iptpif);
                    imap(12, isg, ibl) = nbi1(iptpif);
                }
            }
            if (irev2 == 1) {
                if (imap(8, isg, ibl) > 0) {
                    imap(11, isg, ibl) = nei1(iptpif);
                    imap(12, isg, ibl) = nbi1(iptpif);
                } else {
                    imap(9,  isg, ibl) = nek1(iptpif);
                    imap(10, isg, ibl) = nbk1(iptpif);
                }
            }
        }

        // set imap(1,isg,ibl) = 0 for 1-1 segments that map to themselves (block 2)
        if (ibif1(iptpif) == ibif2(iptpif)) {
            ift1   = abs(imap(2,  isg, ibl));
            ift2   = abs(imap(8,  isg, ibl));
            idel11 = abs(imap(3,  isg, ibl) - imap(4,  isg, ibl));
            idel12 = abs(imap(5,  isg, ibl) - imap(6,  isg, ibl));
            idel21 = abs(imap(9,  isg, ibl) - imap(10, isg, ibl));
            idel22 = abs(imap(11, isg, ibl) - imap(12, isg, ibl));
            if ((ift1 == ift2) && (idel11 == idel21) && (idel12 == idel22))
                imap(1, isg, ibl) = 0;
        }

    } // end do 230



    // eliminate double-counted 1-1 interfaces
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        // identify those segments that are redundant, and flag them with iredundant() = 1
        idebug = 0;
        for (is = 1; is <= msegn; is++) {
            iredundant(is, ibloc) = 0;
        }
        for (is1 = 1; is1 <= nseg(ibloc); is1++) {
            im11  = imap(1,  is1, ibloc);
            im12  = imap(2,  is1, ibloc);
            im13  = imap(3,  is1, ibloc);
            im14  = imap(4,  is1, ibloc);
            im15  = imap(5,  is1, ibloc);
            im16  = imap(6,  is1, ibloc);
            im17  = imap(7,  is1, ibloc);
            im18  = imap(8,  is1, ibloc);
            im19  = imap(9,  is1, ibloc);
            im110 = imap(10, is1, ibloc);
            im111 = imap(11, is1, ibloc);
            im112 = imap(12, is1, ibloc);
            for (is2 = 1; is2 <= nseg(ibloc); is2++) {
                if (iredundant(is1, ibloc) == 0) {
                    im21  = imap(1,  is2, ibloc);
                    im22  = imap(2,  is2, ibloc);
                    im23  = imap(3,  is2, ibloc);
                    im24  = imap(4,  is2, ibloc);
                    im25  = imap(5,  is2, ibloc);
                    im26  = imap(6,  is2, ibloc);
                    im27  = imap(7,  is2, ibloc);
                    im28  = imap(8,  is2, ibloc);
                    im29  = imap(9,  is2, ibloc);
                    im210 = imap(10, is2, ibloc);
                    im211 = imap(11, is2, ibloc);
                    im212 = imap(12, is2, ibloc);
                    if (is2 != is1) {
                        if ((im11 == im21)   && (im12 == im22)   &&
                            (im12 == im22)   && (im13 == im23)   &&
                            (im13 == im23)   && (im14 == im24)   &&
                            (im15 == im25)   && (im16 == im26)   &&
                            (im17 == im27)   && (im18 == im28)   &&
                            (im19 == im29)   && (im110 == im210) &&
                            (im111 == im211) && (im112 == im212)) {
                            iredundant(is2, ibloc) = 1;
                            if (idebug > 0) {
                                printf("block %d seg %d is redundant to %d\n", ibloc, is1, is2);
                                printf("iredundant(is1) = %d redundant(is2) = %d\n",
                                       iredundant(is1, ibloc), iredundant(is2, ibloc));
                            }
                        }
                    }
                }
            }
        }
    }



    // read remainder of input file
    if (iver == 4) {
        // patch data
        fortran_read_list(10);
        fortran_read_list(10);
        fortran_read_list(10, &ninter);

        // plot3d data
        fortran_read_list(10);
        fortran_read_list(10);
        if (abs(cmn_info3.nplot3d) > 0) {
            for (mm = 1; mm <= abs(cmn_info3.nplot3d); mm++) {
                fortran_read_list(10, &idum, &iplt3dtyp_local);
            }
        }

        // movie data
        fortran_read_list(10);
        if (abs(cmn_info3.nplot3d) > 0) {
            fortran_read_list(10, &idum);
        }

        // print out data
        fortran_read_list(10);
        fortran_read_list(10);
        if (abs(cmn_info3.nprint) > 0) {
            for (mm = 1; mm <= abs(cmn_info3.nprint); mm++) {
                fortran_read_list(10);
            }
        }

        // sensitivity data - with end-of-file handling
        {
            FILE* f10 = fortran_get_unit(10);
            char linebuf[256];
            // read(10,*,end=8989)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_8989;
            // read(10,*,end=8989)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_8989;
            // read(10,*,end=8989) ndv,isdform
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_8989;
            sscanf(linebuf, "%d %d", &cmn_info3.ndv, &cmn_info3.isdform);
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_8989;
            // read(10,'(a60)') sdgridindum
            if (fgets(sdgridindum, sizeof(sdgridindum), f10) == nullptr) goto label_8989;
            {
                int l = (int)strlen(sdgridindum);
                if (l > 0 && sdgridindum[l-1] == '\n') sdgridindum[l-1] = '\0';
            }
            // read(10,'(a60)') dovrlap
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_8989;
            {
                int l = (int)strlen(linebuf);
                if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
                strncpy(cmn_cflfiles.dovrlap, linebuf, 80);
            }
            // read(10,'(a60)') dpatch
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_8989;
            {
                int l = (int)strlen(linebuf);
                if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
                strncpy(cmn_cflfiles.dpatch, linebuf, 80);
            }
            // read(10,'(a60)') dresid
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_8989;
            {
                int l = (int)strlen(linebuf);
                if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
                strncpy(cmn_cflfiles.dresid, linebuf, 80);
            }
        }
        label_8989:;

    } else {
        // patch data
        fortran_read_list(10);
        fortran_read_list(10);
        fortran_read_list(10, &ninter);

        // plot3d data
        fortran_read_list(10);
        fortran_read_list(10);
        if (abs(cmn_info3.nplot3d) > 0) {
            cmn_info3.nplot3d = abs(cmn_info3.nplot3d);
            for (mm = 1; mm <= cmn_info3.nplot3d; mm++) {
                fortran_read_list(10);
            }
        }

        // movie data
        fortran_read_list(10);
        fortran_read_list(10, &idum);

        // print out data
        fortran_read_list(10);
        fortran_read_list(10);
        if (abs(cmn_info3.nprint) > 0) {
            cmn_info3.nprint = abs(cmn_info3.nprint);
            for (mm = 1; mm <= cmn_info3.nprint; mm++) {
                fortran_read_list(10);
            }
        }

        // control surfaces
        fortran_read_list(10);
        fortran_read_list(10);
        fortran_read_list(10, &ncs);
        fortran_read_list(10);
        if (ncs > 0) {
            for (mm = 1; mm <= ncs; mm++) {
                fortran_read_list(10);
            }
        }

        // dynamic mesh data - with end-of-file handling
        {
            FILE* f10 = fortran_get_unit(10);
            char linebuf[256];
            // read(10,*,end=999)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            // read(10,*) ntrans
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            sscanf(linebuf, "%d", &ntrans);
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            // if (ntrans.gt.0) read(10,*) tlref
            if (ntrans > 0) {
                if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
                sscanf(linebuf, "%lf", &tlref);
            }
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            if (ntrans > 0) {
                for (n = 1; n <= ntrans; n++) {
                    if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
                    sscanf(linebuf, "%d %d %f %f %f %f", &igd, &itrans_loc, &rfreqt, &xmag, &ymag, &zmag);
                    xmap(8,  1, igd) = (double)itrans_loc;
                    xmap(9,  1, igd) = (double)rfreqt;
                    xmap(10, 1, igd) = (double)xmag;
                    xmap(11, 1, igd) = (double)ymag;
                    xmap(12, 1, igd) = (double)zmag;
                }
            }
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            if (ntrans > 0) {
                for (n = 1; n <= ntrans; n++) {
                    if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
                    sscanf(linebuf, "%d %f %f %f", &igd, &dxmax, &dymax, &dzmax);
                    xmap(13, 1, igd) = (double)dxmax;
                    xmap(14, 1, igd) = (double)dymax;
                    xmap(15, 1, igd) = (double)dzmax;
                }
            }
        }



        // rotation data
        {
            FILE* f10 = fortran_get_unit(10);
            char linebuf[256];
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            // read(10,*) nrotat
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            sscanf(linebuf, "%d", &nrotat);
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            // if (nrotat.gt.0) read(10,*) rlref
            if (nrotat > 0) {
                if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
                sscanf(linebuf, "%lf", &rlref);
            }
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            if (nrotat > 0) {
                for (n = 1; n <= nrotat; n++) {
                    if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
                    sscanf(linebuf, "%d %d %f %f %f %f %f %f %f",
                           &igd, &irotat, &rfreqr, &thxmag, &thymag, &thzmag,
                           &xorig, &yorig, &zorig);
                    for (iseg = 1; iseg <= msegn; iseg++) {
                        xmap(16, iseg, igd) = (double)irotat;
                        xmap(17, iseg, igd) = (double)rfreqr;
                        xmap(18, iseg, igd) = (double)thxmag;
                        xmap(19, iseg, igd) = (double)thymag;
                        xmap(20, iseg, igd) = (double)thzmag;
                        xmap(21, iseg, igd) = (double)xorig;
                        xmap(22, iseg, igd) = (double)yorig;
                        xmap(23, iseg, igd) = (double)zorig;
                    }
                }
            }
            // read(10,*)
            if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
            if (nrotat > 0) {
                for (n = 1; n <= nrotat; n++) {
                    if (fgets(linebuf, sizeof(linebuf), f10) == nullptr) goto label_999;
                    sscanf(linebuf, "%d %f %f %f", &igd, &thxmax, &thymax, &thzmax);
                    for (iseg = 1; iseg <= msegn; iseg++) {
                        xmap(24, iseg, igd) = (double)thxmax;
                        xmap(25, iseg, igd) = (double)thymax;
                        xmap(26, iseg, igd) = (double)thzmax;
                    }
                }
            }
        }
        // ignore dynamic mesh data for now!
        label_999:;

    } // end if (iver == 4) else



    // read in patched grid info from ronnie input file, if applicable
    if (ipatch != 0) {

        fortran_read_list(15);
        fortran_read_list(15);
        {
            FILE* f15 = fortran_get_unit(15);
            char linebuf[256];
            if (fgets(linebuf, sizeof(linebuf), f15)) {
                int l = (int)strlen(linebuf);
                if (l > 0 && linebuf[l-1] == '\n') linebuf[l-1] = '\0';
                strncpy(cmn_ronfiles.rout, linebuf, 80);
            }
        }
        fortran_read_list(15);
        if (iptyp > 0) {
            fortran_read_list(15);
            fortran_read_list(15, &ioflag, &itrace);
        } else {
            ioflag = 0;
        }
        // read(15,10) (titleron(i),i=1,20)  format(20a4)
        {
            FILE* f15 = fortran_get_unit(15);
            char linebuf[256];
            if (fgets(linebuf, sizeof(linebuf), f15)) {
                memset(cmn_ron.titleron, 0, sizeof(cmn_ron.titleron));
                int copylen = (int)strlen(linebuf);
                if (copylen > 80) copylen = 80;
                memcpy(cmn_ron.titleron, linebuf, copylen);
            }
        }
        fortran_read_list(15);
        fortran_read_list(15);
        fortran_read_list(15);

        // do 300
        for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
            fortran_read_list(15);
        }

        fortran_read_list(15);
        fortran_read_list(15, &npatif);
        fortran_read_list(15);

        // do 310
        for (intr = 1; intr <= npatif; intr++) {
            fortran_read_list(15);
        }

        fortran_read_list(15);
        if (ioflag != 0) fortran_read_list(15);

        // do 350
        for (ipif = 1; ipif <= npatif; ipif++) {
            if (ioflag == 0) {
                // type 1 ronnie input
                // read(15,*) iif,iftar,nb1t,ne1t,nb2t,ne2t,nfsor,(ifsor(is),is=1,nfsor)
                // This reads variable number of items; use raw file read
                {
                    FILE* f15 = fortran_get_unit(15);
                    char linebuf[1024];
                    if (fgets(linebuf, sizeof(linebuf), f15)) {
                        // parse: iif iftar nb1t ne1t nb2t ne2t nfsor ifsor(1)...ifsor(nfsor)
                        int pos = 0;
                        int nread;
                        nread = sscanf(linebuf + pos, "%d %d %d %d %d %d %d%n",
                                       &iif, &iftar, &nb1t, &ne1t, &nb2t, &ne2t, &nfsor_loc, &pos);
                        for (int is_idx = 1; is_idx <= nfsor_loc; is_idx++) {
                            int tmp_pos = 0;
                            sscanf(linebuf + pos, " %d%n", &ifsor(is_idx), &tmp_pos);
                            pos += tmp_pos;
                        }
                    }
                }

                // parse target face value & set up index ranges
                iblt   = iftar / 100;
                iftar  = iftar - iblt * 100;
                iijkt  = iftar / 10;
                imnmxt = iftar - iijkt * 10;
                if (iijkt == 1) {
                    ift = 1;
                    if (imnmxt == 2) ift = 2;
                } else if (iijkt == 2) {
                    ift = 3;
                    if (imnmxt == 2) ift = 4;
                } else {
                    ift = 5;
                    if (imnmxt == 2) ift = 6;
                }

                if (ift <= 2) {
                    // i-face
                    nbeg1t = nb1t;
                    nend1t = ne1t;
                    nbeg2t = nb2t;
                    nend2t = ne2t;
                    if (nbeg1t == 0) nbeg1t = 1;
                    if (nend1t == 0) nend1t = jl(iblt);
                    if (nbeg2t == 0) nbeg2t = 1;
                    if (nend2t == 0) nend2t = kl(iblt);
                } else if (ift >= 5) {
                    // k-face
                    nbeg1t = nb2t;
                    nend1t = ne2t;
                    nbeg2t = nb1t;
                    nend2t = ne1t;
                    if (nbeg1t == 0) nbeg1t = 1;
                    if (nend1t == 0) nend1t = il(iblt);
                    if (nbeg2t == 0) nbeg2t = 1;
                    if (nend2t == 0) nend2t = jl(iblt);
                } else {
                    // j-face
                    nbeg1t = nb1t;
                    nend1t = ne1t;
                    nbeg2t = nb2t;
                    nend2t = ne2t;
                    if (nbeg1t == 0) nbeg1t = 1;
                    if (nend1t == 0) nend1t = kl(iblt);
                    if (nbeg2t == 0) nbeg2t = 1;
                    if (nend2t == 0) nend2t = il(iblt);
                }

                if (nend1t < nbeg1t) { ifoo = nend1t; nend1t = nbeg1t; nbeg1t = ifoo; }
                if (nend2t < nbeg2t) { ifoo = nend2t; nend2t = nbeg2t; nbeg2t = ifoo; }

                // parse source face values and set up index ranges (do 320)
                for (i = 1; i <= nfsor_loc; i++) {
                    ibls     = ifsor(i) / 100;
                    ifsor(i) = ifsor(i) - ibls * 100;
                    iijks    = ifsor(i) / 10;
                    imnmxs   = ifsor(i) - iijks * 10;
                    if (iijks == 1) {
                        ifs = 1;
                        if (imnmxs == 2) ifs = 2;
                    } else if (iijks == 2) {
                        ifs = 3;
                        if (imnmxs == 2) ifs = 4;
                    } else {
                        ifs = 5;
                        if (imnmxs == 2) ifs = 6;
                    }

                    if (ifs <= 2) {
                        nbeg1s = 1; nend1s = jl(ibls); nbeg2s = 1; nend2s = kl(ibls);
                    } else if (ifs >= 5) {
                        nbeg1s = 1; nend1s = il(ibls); nbeg2s = 1; nend2s = jl(ibls);
                    } else {
                        nbeg1s = 1; nend1s = kl(ibls); nbeg2s = 1; nend2s = il(ibls);
                    }

                    nseg(iblt) = nseg(iblt) + 1;
                    it = nseg(iblt);
                    if (nseg(iblt) > msegn) {
                        printf("stopping: parameter msegn must be at least %d\n", nseg(iblt));
                        termn8_ns::termn8(zero_val, neg7, ibufdim, nbuf, bou, nou);
                    }
                    imap(1,  it, iblt) = -1;
                    imap(2,  it, iblt) = ift;
                    imap(3,  it, iblt) = nbeg1t;
                    imap(4,  it, iblt) = nend1t;
                    imap(5,  it, iblt) = nbeg2t;
                    imap(6,  it, iblt) = nend2t;
                    imap(7,  it, iblt) = ibls;
                    imap(8,  it, iblt) = ifs;
                    imap(9,  it, iblt) = nbeg1s;
                    imap(10, it, iblt) = nend1s;
                    imap(11, it, iblt) = nbeg2s;
                    imap(12, it, iblt) = nend2s;
                    imap(21, it, iblt) = 0;
                    imap(22, it, iblt) = 1;
                    imap(23, it, iblt) = iovrlp(iblt);
                } // end do 320

            } else {
                // type 2 ronnie input
                // read(15,*) iif,iftar,nb1t,ne1t,nb2t,ne2t,nfsor
                fortran_read_list(15, &iif, &iftar, &nb1t, &ne1t, &nb2t, &ne2t, &nfsor_loc);
                // do 330
                for (int is_idx = 1; is_idx <= nfsor_loc; is_idx++) {
                    fortran_read_list(15, &ifsor(is_idx), &nb1s(is_idx), &ne1s(is_idx),
                                     &nb2s(is_idx), &ne2s(is_idx));
                }



                // parse target face value & set up index ranges (type 2)
                iblt   = iftar / 100;
                iftar  = iftar - iblt * 100;
                iijkt  = iftar / 10;
                imnmxt = iftar - iijkt * 10;
                if (iijkt == 1) {
                    ift = 1;
                    if (imnmxt == 2) ift = 2;
                } else if (iijkt == 2) {
                    ift = 3;
                    if (imnmxt == 2) ift = 4;
                } else {
                    ift = 5;
                    if (imnmxt == 2) ift = 6;
                }

                if (ift <= 2) {
                    nbeg1t = nb1t; nend1t = ne1t; nbeg2t = nb2t; nend2t = ne2t;
                    if (nbeg1t == 0) nbeg1t = 1;
                    if (nend1t == 0) nend1t = jl(iblt);
                    if (nbeg2t == 0) nbeg2t = 1;
                    if (nend2t == 0) nend2t = kl(iblt);
                } else if (ift >= 5) {
                    nbeg1t = nb2t; nend1t = ne2t; nbeg2t = nb1t; nend2t = ne1t;
                    if (nbeg1t == 0) nbeg1t = 1;
                    if (nend1t == 0) nend1t = il(iblt);
                    if (nbeg2t == 0) nbeg2t = 1;
                    if (nend2t == 0) nend2t = jl(iblt);
                } else {
                    nbeg1t = nb1t; nend1t = ne1t; nbeg2t = nb2t; nend2t = ne2t;
                    if (nbeg1t == 0) nbeg1t = 1;
                    if (nend1t == 0) nend1t = kl(iblt);
                    if (nbeg2t == 0) nbeg2t = 1;
                    if (nend2t == 0) nend2t = il(iblt);
                }

                if (nend1t < nbeg1t) { ifoo = nend1t; nend1t = nbeg1t; nbeg1t = ifoo; }
                if (nend2t < nbeg2t) { ifoo = nend2t; nend2t = nbeg2t; nbeg2t = ifoo; }

                // parse source face values and set up index ranges (do 340)
                for (i = 1; i <= nfsor_loc; i++) {
                    ibls     = ifsor(i) / 100;
                    ifsor(i) = ifsor(i) - ibls * 100;
                    iijks    = ifsor(i) / 10;
                    imnmxs   = ifsor(i) - iijks * 10;
                    if (iijks == 1) {
                        ifs = 1;
                        if (imnmxs == 2) ifs = 2;
                    } else if (iijks == 2) {
                        ifs = 3;
                        if (imnmxs == 2) ifs = 4;
                    } else {
                        ifs = 5;
                        if (imnmxs == 2) ifs = 6;
                    }

                    if (ifs <= 2) {
                        nbeg1s = nb1s(i); nend1s = ne1s(i); nbeg2s = nb2s(i); nend2s = ne2s(i);
                        if (nbeg1s == 0) nbeg1s = 1;
                        if (nend1s == 0) nend1s = jl(ibls);
                        if (nbeg2s == 0) nbeg2s = 1;
                        if (nend2s == 0) nend2s = kl(ibls);
                    } else if (ifs >= 5) {
                        nbeg1s = nb2s(i); nend1s = ne2s(i); nbeg2s = nb1s(i); nend2s = ne1s(i);
                        if (nbeg1s == 0) nbeg1s = 1;
                        if (nend1s == 0) nend1s = il(ibls);
                        if (nbeg2s == 0) nbeg2s = 1;
                        if (nend2s == 0) nend2s = jl(ibls);
                    } else {
                        nbeg1s = nb1s(i); nend1s = ne1s(i); nbeg2s = nb2s(i); nend2s = ne2s(i);
                        if (nbeg1s == 0) nbeg1s = 1;
                        // Note: Fortran has bug here: "if(nbeg1s .eq. 0) nend1s = kl(ibls)"
                        // and "if(nend2s .eq. 0) nbeg2s = 1" - faithfully translated:
                        if (nbeg1s == 0) nend1s = kl(ibls);
                        if (nend2s == 0) nbeg2s = 1;
                        if (nend2s == 0) nend2s = il(ibls);
                    }

                    if (nend1s < nbeg1s) { ifoo = nend1s; nend1s = nbeg1s; nbeg1s = ifoo; }
                    if (nend2s < nbeg2s) { ifoo = nend2s; nend2s = nbeg2s; nbeg2s = ifoo; }

                    nseg(iblt) = nseg(iblt) + 1;
                    if (nseg(iblt) > msegn) {
                        printf("stopping: parameter msegn must be at least %d\n", nseg(iblt));
                        termn8_ns::termn8(zero_val, neg7, ibufdim, nbuf, bou, nou);
                    }
                    it = nseg(iblt);
                    imap(1,  it, iblt) = -1;
                    imap(2,  it, iblt) = ift;
                    imap(3,  it, iblt) = nbeg1t;
                    imap(4,  it, iblt) = nend1t;
                    imap(5,  it, iblt) = nbeg2t;
                    imap(6,  it, iblt) = nend2t;
                    imap(7,  it, iblt) = ibls;
                    imap(8,  it, iblt) = ifs;
                    imap(9,  it, iblt) = nbeg1s;
                    imap(10, it, iblt) = nend1s;
                    imap(11, it, iblt) = nbeg2s;
                    imap(12, it, iblt) = nend2s;
                    imap(21, it, iblt) = 0;
                    imap(22, it, iblt) = 1;
                    imap(23, it, iblt) = iovrlp(iblt);
                } // end do 340

            } // end if (ioflag == 0) else

        } // end do 350

    } // end if (ipatch != 0)



    // set parameters for baldwin-lomax model: imap(13-20,it,iblt)
    // (these are set non-zero only if the wall is viscous)
    for (int iblt = 1; iblt <= cmn_grdinfo.nbloc; iblt++) {
        for (it = 1; it <= nseg(iblt); it++) {

            imap(13, it, iblt) = 0;
            imap(14, it, iblt) = 0;
            imap(15, it, iblt) = 0;
            imap(16, it, iblt) = 0;
            imap(17, it, iblt) = 0;
            imap(18, it, iblt) = 0;
            imap(19, it, iblt) = 0;
            imap(20, it, iblt) = 0;

            idir = 2;
            if (imap(2, it, iblt) <= 2) idir = 1;
            if (imap(2, it, iblt) >= 5) idir = 3;

            if (imap(1, it, iblt) == 3 &&
                (abs(ivisb(idir, iblt)) >= 2 ||
                 abs(ivisb(idir, iblt)) <= 3)) {
                imap(13, it, iblt) = 1;
                imap(14, it, iblt) = imap(3, it, iblt);
                imap(15, it, iblt) = imap(4, it, iblt);
                imap(16, it, iblt) = imap(5, it, iblt);
                imap(17, it, iblt) = imap(6, it, iblt);
                // ranges for normal direction
                if (imap(2, it, iblt) == 1) {
                    imap(18, it, iblt) = 3;
                    imap(19, it, iblt) = il(iblt);
                    imap(20, it, iblt) = il(iblt);
                } else if (imap(2, it, iblt) == 2) {
                    imap(18, it, iblt) = il(iblt) - 2;
                    imap(19, it, iblt) = 1;
                    imap(20, it, iblt) = 1;
                } else if (imap(2, it, iblt) == 3) {
                    imap(18, it, iblt) = 3;
                    imap(19, it, iblt) = jl(iblt);
                    imap(20, it, iblt) = jl(iblt);
                } else if (imap(2, it, iblt) == 4) {
                    imap(18, it, iblt) = jl(iblt) - 2;
                    imap(19, it, iblt) = 1;
                    imap(20, it, iblt) = 1;
                } else if (imap(2, it, iblt) == 5) {
                    imap(18, it, iblt) = 3;
                    imap(19, it, iblt) = kl(iblt);
                    imap(20, it, iblt) = kl(iblt);
                } else if (imap(2, it, iblt) == 6) {
                    imap(18, it, iblt) = kl(iblt) - 2;
                    imap(19, it, iblt) = 1;
                    imap(20, it, iblt) = 1;
                }
            }

            // default to single tref
            twotref(it, iblt) = 0.0;

        }
    }



    // check for block faces that have not been set
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        ifchk1 = 0; ifchk2 = 0; ifchk3 = 0;
        ifchk4 = 0; ifchk5 = 0; ifchk6 = 0;
        for (iseg = 1; iseg <= nseg(ibloc); iseg++) {
            if (imap(2, iseg, ibloc) == 1) ifchk1 = 1;
            if (imap(2, iseg, ibloc) == 2) ifchk2 = 1;
            if (imap(2, iseg, ibloc) == 3) ifchk3 = 1;
            if (imap(2, iseg, ibloc) == 4) ifchk4 = 1;
            if (imap(2, iseg, ibloc) == 5) ifchk5 = 1;
            if (imap(2, iseg, ibloc) == 6) ifchk6 = 1;
        }
        if (ifchk1 == 0) printf("WARNING: block %d face imin, has no bc data\n", ibloc);
        if (ifchk2 == 0) printf("WARNING: block %d face imax, has no bc data\n", ibloc);
        if (ifchk3 == 0) printf("WARNING: block %d face jmin, has no bc data\n", ibloc);
        if (ifchk4 == 0) printf("WARNING: block %d face jmax, has no bc data\n", ibloc);
        if (ifchk5 == 0) printf("WARNING: block %d face kmin, has no bc data\n", ibloc);
        if (ifchk6 == 0) printf("WARNING: block %d face kmax, has no bc data\n", ibloc);

        ifill = 0;
        if (ifchk1 == 0 || ifchk2 == 0 || ifchk3 == 0 ||
            ifchk4 == 0 || ifchk5 == 0 || ifchk6 == 0) {
            printf(" will fill this missing data with extrapolation\n");
            ifill = 1;
        }

        // fill in missing faces with extrapolation
        if (ifchk1 == 0 && ifill == 1) {
            nseg(ibloc) = nseg(ibloc) + 1;
            imap(1, nseg(ibloc), ibloc) = 1002;
            imap(2, nseg(ibloc), ibloc) = 1;
            imap(3, nseg(ibloc), ibloc) = 1;
            imap(4, nseg(ibloc), ibloc) = jl(ibloc);
            imap(5, nseg(ibloc), ibloc) = 1;
            imap(6, nseg(ibloc), ibloc) = kl(ibloc);
            for (ii = 7; ii <= 20; ii++) imap(ii, nseg(ibloc), ibloc) = 0;
        }
        if (ifchk2 == 0 && ifill == 1) {
            nseg(ibloc) = nseg(ibloc) + 1;
            imap(1, nseg(ibloc), ibloc) = 1002;
            imap(2, nseg(ibloc), ibloc) = 2;
            imap(3, nseg(ibloc), ibloc) = 1;
            imap(4, nseg(ibloc), ibloc) = jl(ibloc);
            imap(5, nseg(ibloc), ibloc) = 1;
            imap(6, nseg(ibloc), ibloc) = kl(ibloc);
            for (ii = 7; ii <= 20; ii++) imap(ii, nseg(ibloc), ibloc) = 0;
        }
        if (ifchk3 == 0 && ifill == 1) {
            nseg(ibloc) = nseg(ibloc) + 1;
            imap(1, nseg(ibloc), ibloc) = 1002;
            imap(2, nseg(ibloc), ibloc) = 3;
            imap(3, nseg(ibloc), ibloc) = 1;
            imap(4, nseg(ibloc), ibloc) = kl(ibloc);
            imap(5, nseg(ibloc), ibloc) = 1;
            imap(6, nseg(ibloc), ibloc) = il(ibloc);
            for (ii = 7; ii <= 20; ii++) imap(ii, nseg(ibloc), ibloc) = 0;
        }
        if (ifchk4 == 0 && ifill == 1) {
            nseg(ibloc) = nseg(ibloc) + 1;
            imap(1, nseg(ibloc), ibloc) = 1002;
            imap(2, nseg(ibloc), ibloc) = 4;
            imap(3, nseg(ibloc), ibloc) = 1;
            imap(4, nseg(ibloc), ibloc) = kl(ibloc);
            imap(5, nseg(ibloc), ibloc) = 1;
            imap(6, nseg(ibloc), ibloc) = il(ibloc);
            for (ii = 7; ii <= 20; ii++) imap(ii, nseg(ibloc), ibloc) = 0;
        }
        if (ifchk5 == 0 && ifill == 1) {
            nseg(ibloc) = nseg(ibloc) + 1;
            imap(1, nseg(ibloc), ibloc) = 1002;
            imap(2, nseg(ibloc), ibloc) = 5;
            imap(3, nseg(ibloc), ibloc) = 1;
            imap(4, nseg(ibloc), ibloc) = il(ibloc);
            imap(5, nseg(ibloc), ibloc) = 1;
            imap(6, nseg(ibloc), ibloc) = jl(ibloc);
            for (ii = 7; ii <= 20; ii++) imap(ii, nseg(ibloc), ibloc) = 0;
        }
        if (ifchk6 == 0 && ifill == 1) {
            nseg(ibloc) = nseg(ibloc) + 1;
            imap(1, nseg(ibloc), ibloc) = 1002;
            imap(2, nseg(ibloc), ibloc) = 6;
            imap(3, nseg(ibloc), ibloc) = 1;
            imap(4, nseg(ibloc), ibloc) = il(ibloc);
            imap(5, nseg(ibloc), ibloc) = 1;
            imap(6, nseg(ibloc), ibloc) = jl(ibloc);
            for (ii = 7; ii <= 20; ii++) imap(ii, nseg(ibloc), ibloc) = 0;
        }
    } // end do ibloc



    // determine nseg1(ibloc) that doesn't include redundant 1-1 segments
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {
        nseg1(ibloc) = nseg(ibloc);
        for (is = 1; is <= nseg(ibloc); is++) {
            if (iredundant(is, ibloc) == 1) {
                nseg1(ibloc) = nseg1(ibloc) - 1;
            }
        }
    }

    // write out map file (unit 20)
    fortran_write_unit(20, "  nbloc\n");
    fortran_write_unit(20, "%5d\n", cmn_grdinfo.nbloc);
    fortran_write_unit(20, " \n");

    // do 600
    for (ibloc = 1; ibloc <= cmn_grdinfo.nbloc; ibloc++) {

        fortran_write_unit(20, "  nseg  ivisb1 ivisb2 ivisb3 iturbb itrb1 itrb2 jtrb1 jtrb2 ktrb1 ktrb2\n");
        fortran_write_unit(20, "%7d%7d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                           nseg1(ibloc), ivisb(1, ibloc), ivisb(2, ibloc), ivisb(3, ibloc),
                           iturbb(ibloc), itrb1(ibloc), itrb2(ibloc), jtrb1(ibloc),
                           jtrb2(ibloc), ktrb1(ibloc), ktrb2(ibloc));
        fortran_write_unit(20, "begin non-standard map data\n");
        fortran_write_unit(20, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                           ivisb(4,  ibloc), ivisb(5,  ibloc), ivisb(6,  ibloc),
                           ivisb(7,  ibloc), ivisb(8,  ibloc), ivisb(9,  ibloc),
                           ivisb(10, ibloc), ivisb(11, ibloc), ivisb(12, ibloc),
                           ivisb(13, ibloc), ivisb(14, ibloc), ivisb(15, ibloc));
        fortran_write_unit(20, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                           ivisb(16, ibloc), ivisb(17, ibloc), ivisb(18, ibloc),
                           ivisb(19, ibloc), ivisb(20, ibloc), ivisb(21, ibloc),
                           ivisb(22, ibloc), ivisb(23, ibloc), ivisb(24, ibloc),
                           ivisb(25, ibloc), ivisb(26, ibloc), ivisb(27, ibloc));
        fortran_write_unit(20, "%6d%6d%6d%6d%6d%6d%6d\n",
                           ivisb(28, ibloc), ivisb(29, ibloc), ivisb(30, ibloc),
                           ivisb(31, ibloc), ivisb(32, ibloc), ivisb(33, ibloc),
                           ivisb(34, ibloc));
        fortran_write_unit(20, "end non-standard map data\n");
        fortran_write_unit(20, "  nbt  nst   typ  nft 1tmn 1tmx 2tmn 2tmx  nbs  nfs 1smn 1smx 2smn 2smx\n");
        fortran_write_unit(20, "               ntrb 1beg 1end 2beg 2end fmx1 fmx2 nvis twotrf\n");
        fortran_write_unit(20, " \n");

        // write out imap array
        iiseg = 0;
        for (iseg = 1; iseg <= nseg(ibloc); iseg++) {
            if (iredundant(iseg, ibloc) == 0) {
                iiseg = iiseg + 1;
                fortran_write_unit(20, "%5d%5d%6d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d\n",
                                   ibloc, iiseg,
                                   imap(1,  iseg, ibloc), imap(2,  iseg, ibloc),
                                   imap(3,  iseg, ibloc), imap(4,  iseg, ibloc),
                                   imap(5,  iseg, ibloc), imap(6,  iseg, ibloc),
                                   imap(7,  iseg, ibloc), imap(8,  iseg, ibloc),
                                   imap(9,  iseg, ibloc), imap(10, iseg, ibloc),
                                   imap(11, iseg, ibloc), imap(12, iseg, ibloc));
                fortran_write_unit(20, "               %5d%5d%5d%5d%5d%5d%5d%5d%10.4f\n",
                                   imap(13, iseg, ibloc), imap(14, iseg, ibloc),
                                   imap(15, iseg, ibloc), imap(16, iseg, ibloc),
                                   imap(17, iseg, ibloc), imap(18, iseg, ibloc),
                                   imap(19, iseg, ibloc), imap(20, iseg, ibloc),
                                   twotref(iseg, ibloc));
                fortran_write_unit(20, "                   0.0   0.0  0.0  1.0  0.0\n");
                fortran_write_unit(20, "               begin non-standard map data\n");
                fortran_write_unit(20, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                                   xmap(1,  iseg, ibloc), xmap(2,  iseg, ibloc),
                                   xmap(3,  iseg, ibloc), xmap(4,  iseg, ibloc),
                                   xmap(5,  iseg, ibloc));
                fortran_write_unit(20, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                                   xmap(6,  iseg, ibloc), xmap(7,  iseg, ibloc),
                                   xmap(8,  iseg, ibloc), xmap(9,  iseg, ibloc),
                                   xmap(10, iseg, ibloc));
                fortran_write_unit(20, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                                   xmap(11, iseg, ibloc), xmap(12, iseg, ibloc),
                                   xmap(13, iseg, ibloc), xmap(14, iseg, ibloc),
                                   xmap(15, iseg, ibloc));
                fortran_write_unit(20, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                                   xmap(16, iseg, ibloc), xmap(17, iseg, ibloc),
                                   xmap(18, iseg, ibloc), xmap(19, iseg, ibloc),
                                   xmap(20, iseg, ibloc));
                fortran_write_unit(20, "               %10.4f%10.4f%10.4f%10.4f%10.4f\n",
                                   xmap(21, iseg, ibloc), xmap(22, iseg, ibloc),
                                   xmap(23, iseg, ibloc), xmap(24, iseg, ibloc),
                                   xmap(25, iseg, ibloc));
                fortran_write_unit(20, "               %10.4f\n",
                                   xmap(26, iseg, ibloc));
                fortran_write_unit(20, "               %5d%5d%5d\n",
                                   imap(21, iseg, ibloc), imap(22, iseg, ibloc),
                                   imap(23, iseg, ibloc));
                fortran_write_unit(20, "               end non-standard map data\n");
                fortran_write_unit(20, " \n");
            }
        } // end do 600 iseg

    } // end do 600 ibloc

    // rewind unsplit input file so it can be read again in order
    // to set keyword data for the split input file.
    rewind(fortran_get_unit(10));

} // end cfl3d_to_tlns3d

} // namespace cfl3d_to_tlns3d_ns
