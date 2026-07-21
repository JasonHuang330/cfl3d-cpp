// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include <cstdio>
#include "bc.h"
#include "bc1000.h"
#include "bc1001.h"
#include "bc1002.h"
#include "bc1003.h"
#include "bc1005.h"
#include "bc1008.h"
#include "bc1011.h"
#include "bc1012.h"
#include "bc1013.h"
#include "bc2002.h"
#include "bc2003.h"
#include "bc2004.h"
#include "bc2006.h"
#include "bc2007.h"
#include "bc2008.h"
#include "bc2009.h"
#include "bc2010.h"
#include "bc2016.h"
#include "bc2019.h"
#include "bc2026.h"
#include "bc2102.h"
#include "bc9999.h"
#include "bcchk.h"
#include "bcnonin.h"
#include "chksym.h"
#include "chkrap.h"
#include "termn8.h"
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"
#include <cstring>
#include <cstdlib>
#include <algorithm>

namespace bc_ns {

void bc(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
        FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk, int& nwork,
        double& cl, FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
        int& nbuf, int& ibufdim, int& maxbl, int& maxgr, int& maxseg,
        FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
        FortranArray1DRef<int> idefrm, FortranArray1DRef<int> igridg,
        FortranArray1DRef<int> nblg, FortranArray1DRef<int> nbci0,
        FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0,
        FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim,
        FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo,
        FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo,
        FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
        FortranArray3DRef<int> bcfilek, FortranArray3DRef<int> lwdat,
        int& myid, FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
        FortranArray1DRef<int> kdimg, FortranArray1DRef<char[80]> bcfiles,
        int& mxbcfil, int& nummem)
{
    // COMMON block aliases
    int& jdim     = cmn_ginfo.jdim;
    int& kdim     = cmn_ginfo.kdim;
    int& idim     = cmn_ginfo.idim;
    int& lq       = cmn_ginfo.lq;
    int& lqj0     = cmn_ginfo.lqj0;
    int& lqk0     = cmn_ginfo.lqk0;
    int& lqi0     = cmn_ginfo.lqi0;
    int& lsj      = cmn_ginfo.lsj;
    int& lsk      = cmn_ginfo.lsk;
    int& lsi      = cmn_ginfo.lsi;
    int& lbcj     = cmn_ginfo.lbcj;
    int& lbck     = cmn_ginfo.lbck;
    int& lbci     = cmn_ginfo.lbci;
    int& lxtbj    = cmn_ginfo.lxtbj;
    int& lxtbk    = cmn_ginfo.lxtbk;
    int& lxtbi    = cmn_ginfo.lxtbi;
    int& latbj    = cmn_ginfo.latbj;
    int& latbk    = cmn_ginfo.latbk;
    int& latbi    = cmn_ginfo.latbi;
    int& lxib     = cmn_ginfo.lxib;
    int& ltj0     = cmn_ginfo.ltj0;
    int& ltk0     = cmn_ginfo.ltk0;
    int& lti0     = cmn_ginfo.lti0;
    int& lvis     = cmn_ginfo.lvis;
    int& lvj0     = cmn_ginfo.lvj0;
    int& lvk0     = cmn_ginfo.lvk0;
    int& lvi0     = cmn_ginfo.lvi0;
    int& lsnk0    = cmn_ginfo.lsnk0;
    int& lsni0    = cmn_ginfo.lsni0;
    int& lnbl     = cmn_ginfo.lnbl;
    int& lxkb     = cmn_ginfo.lxkb;
    int& lx       = cmn_ginfo.lx;
    int& ly       = cmn_ginfo.ly;
    int& lz       = cmn_ginfo.lz;
    int& lblk     = cmn_ginfo.lblk;
    int& isklton  = cmn_sklton.isklton;
    int& isym     = cmn_chkbc.isym;
    int& jsym     = cmn_chkbc.jsym;
    int& ksym     = cmn_chkbc.ksym;
    int& iwrap    = cmn_chkbc.iwrap;
    int& jwrap    = cmn_chkbc.jwrap;
    int& kwrap    = cmn_chkbc.kwrap;
    int& noninflag = cmn_noninertial.noninflag;

    // Local variables
    char filname[80];
    int isklt_sav, iuns, igrid, jvdim, ixwk;
    int nseg, ista, iend, jsta, jend, ksta, kend;
    int ldata, mdim, ndim;
    int iflg, jvchk;
    int ibcflg, istop;
    int lsta, lend, ll;
    int ngh, nblnum, ichk, ngnew;

    // temporarily set isklton = 0 to suppress output messages...
    // the messages are now output from bc_info
    isklt_sav = isklton;
    isklton   = 0;

    iuns = std::max({irotat(nbl), itrans(nbl), idefrm(nbl)});


    if (isklton > 0) {
        if (iuns == 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " boundary conditions for block%6d %s(grid %6d)",
                          nbl, "", igridg(nbl));
        }
        if (iuns > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " ");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " boundary conditions for block%6d %s(grid %6d) - dynamic mesh",
                          nbl, "", igridg(nbl));
        }
    }

    igrid = igridg(nbl);
    jvdim = std::max({jdim*kdim, jdim*idim, kdim*idim});
    ixwk  = 1;

    // check available storage in work array
    if (isklt_sav > 0) {
        // check storage for bc1003
        iflg = 0;
        for (nseg = 1; nseg <= nbci0(nbl); nseg++)
            if (ibcinfo(nbl,nseg,1,1) == 1003) iflg = 1;
        for (nseg = 1; nseg <= nbcidim(nbl); nseg++)
            if (ibcinfo(nbl,nseg,1,2) == 1003) iflg = 1;
        for (nseg = 1; nseg <= nbcj0(nbl); nseg++)
            if (jbcinfo(nbl,nseg,1,1) == 1003) iflg = 1;
        for (nseg = 1; nseg <= nbcjdim(nbl); nseg++)
            if (jbcinfo(nbl,nseg,1,2) == 1003) iflg = 1;
        for (nseg = 1; nseg <= nbck0(nbl); nseg++)
            if (kbcinfo(nbl,nseg,1,1) == 1003) iflg = 1;
        for (nseg = 1; nseg <= nbckdim(nbl); nseg++)
            if (kbcinfo(nbl,nseg,1,2) == 1003) iflg = 1;
        if (iflg > 0) {
            jvchk = ixwk + jvdim*23 - 1;
            if (jvchk > nwork) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " stopping in bc...insufficient wk storage for bctype 1003");
                int ierrflg = -1;
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }

        // check storage for bc2003
        iflg = 0;
        for (nseg = 1; nseg <= nbci0(nbl); nseg++)
            if (ibcinfo(nbl,nseg,1,1) == 2003) iflg = 1;
        for (nseg = 1; nseg <= nbcidim(nbl); nseg++)
            if (ibcinfo(nbl,nseg,1,2) == 2003) iflg = 1;
        for (nseg = 1; nseg <= nbcj0(nbl); nseg++)
            if (jbcinfo(nbl,nseg,1,1) == 2003) iflg = 1;
        for (nseg = 1; nseg <= nbcjdim(nbl); nseg++)
            if (jbcinfo(nbl,nseg,1,2) == 2003) iflg = 1;
        for (nseg = 1; nseg <= nbck0(nbl); nseg++)
            if (kbcinfo(nbl,nseg,1,1) == 2003) iflg = 1;
        for (nseg = 1; nseg <= nbckdim(nbl); nseg++)
            if (kbcinfo(nbl,nseg,1,2) == 2003) iflg = 1;
        if (iflg > 0) {
            jvchk = ixwk + jvdim*25 - 1;
            if (jvchk > nwork) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, " stopping in bc...insufficient wk storage for bctype 2003");
                int ierrflg = -1;
                termn8_ns::termn8(myid, ierrflg, ibufdim, nbuf, bou, nou);
            }
        }
    }



    // Initialize density and pressure boundary conditions for comparison
    // later in boundary condition check.
    ibcflg = 0;
    if (isklton > 0 && ntime == 1) {
        istop = 1;
        {
            int jd2 = jdim, kd2 = kdim, id2 = idim;
            FortranArray4DRef<double> q_ref(&w(lq), jd2, kd2, id2, 5);
            FortranArray4DRef<double> qi0_ref(&w(lqi0), jd2, kd2, 5, 4);
            FortranArray4DRef<double> qj0_ref(&w(lqj0), kd2, id2-1, 5, 4);
            FortranArray4DRef<double> qk0_ref(&w(lqk0), jd2, id2-1, 5, 4);
            FortranArray3DRef<double> blank_ref(&w(lblk), jd2, kd2, id2);
            bcchk_ns::bcchk(id2, jd2, kd2, q_ref, qi0_ref, qj0_ref, qk0_ref,
                            blank_ref, ibcflg, nbl, nou, bou, nbuf, ibufdim,
                            myid, istop, igridg, maxbl);
        }
    }

    // Initialize boundary condition types to all cell-center type.
    lsta = lbcj;
    lend = lbci + jdim*kdim*2 - 1;
    for (ll = lsta; ll <= lend; ll++)
        w(ll) = 0.0;

    // Update physical boundary conditions.
    // ********************
    // i=0 boundary
    // ********************
    for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
        ista = 1;
        iend = 1;
        jsta = ibcinfo(nbl,nseg,2,1);
        jend = ibcinfo(nbl,nseg,3,1);
        ksta = ibcinfo(nbl,nseg,4,1);
        kend = ibcinfo(nbl,nseg,5,1);
        ldata = lwdat(nbl,nseg,1);
        mdim = jend - jsta;
        ndim = kend - ksta;
        std::memcpy(filname, bcfiles(bcfilei(nbl,nseg,1)), 80);

        int jd = jdim, kd = kdim, id = idim;
        FortranArray4DRef<double> q_r(&w(lq), jd, kd, id, 5);
        FortranArray4DRef<double> qj0_r(&w(lqj0), kd, id-1, 5, 4);
        FortranArray4DRef<double> qk0_r(&w(lqk0), jd, id-1, 5, 4);
        FortranArray4DRef<double> qi0_r(&w(lqi0), jd, kd, 5, 4);
        FortranArray4DRef<double> sj_r(&w(lsj), jd, kd, id-1, 5);
        FortranArray4DRef<double> sk_r(&w(lsk), jd, kd, id-1, 5);
        FortranArray4DRef<double> si_r(&w(lsi), jd, kd, id, 5);
        FortranArray3DRef<double> bcj_r(&w(lbcj), kd, id-1, 2);
        FortranArray3DRef<double> bck_r(&w(lbck), jd, id-1, 2);
        FortranArray3DRef<double> bci_r(&w(lbci), jd, kd, 2);
        FortranArray4DRef<double> xtbj_r(&w(lxtbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> xtbk_r(&w(lxtbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> xtbi_r(&w(lxtbi), jd, kd, 3, 2);
        FortranArray4DRef<double> atbj_r(&w(latbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> atbk_r(&w(latbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> atbi_r(&w(latbi), jd, kd, 3, 2);
        FortranArray4DRef<double> tursav_r(&w(lxib), jd, kd, id, nummem);
        FortranArray4DRef<double> tj0_r(&w(ltj0), kd, id-1, nummem, 4);
        FortranArray4DRef<double> tk0_r(&w(ltk0), jd, id-1, nummem, 4);
        FortranArray4DRef<double> ti0_r(&w(lti0), jd, kd, nummem, 4);
        FortranArray3DRef<double> vist3d_r(&w(lvis), jd, kd, id);
        FortranArray4DRef<double> vj0_r(&w(lvj0), kd, id-1, 1, 4);
        FortranArray4DRef<double> vk0_r(&w(lvk0), jd, id-1, 1, 4);
        FortranArray4DRef<double> vi0_r(&w(lvi0), jd, kd, 1, 4);
        FortranArray3DRef<double> smin_r(&w(lsnk0), jd-1, kd-1, id-1);
        FortranArray3DRef<double> snj0_r(&w(lsni0), kd, id-1, 2);
        FortranArray3DRef<double> snk0_r(&w(lsnk0), jd, id-1, 2);
        FortranArray3DRef<double> sni0_r(&w(lsni0), jd, kd, 2);
        FortranArray3DRef<double> snjm_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> snkm_r(&w(lxkb), jd, kd, id);
        FortranArray3DRef<double> snim_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> x_r(&w(lx), jd, kd, id);
        FortranArray3DRef<double> y_r(&w(ly), jd, kd, id);
        FortranArray3DRef<double> z_r(&w(lz), jd, kd, id);
        FortranArray4DRef<double> bcdata_r(&w(ldata), mdim, ndim, 2, 12);
        FortranArray2DRef<double> wkt_r(&wk(ixwk), jvdim, 23);
        int nface_val = 1;

        if (ibcinfo(nbl,nseg,1,1) == 9999)
            bc9999_ns::bc9999(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, y_r, z_r, smin_r);

        if (ibcinfo(nbl,nseg,1,1) == 1000)
            bc1000_ns::bc1000(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (ibcinfo(nbl,nseg,1,1) == 1000 && noninflag > 0)
            bcnonin_ns::bcnonin(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                iuns, nou, bou, nbuf, ibufdim, x_r, y_r, z_r, nbl);

        if (ibcinfo(nbl,nseg,1,1) == 1001)
            bc1001_ns::bc1001(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, z_r);

        if (ibcinfo(nbl,nseg,1,1) == 1002)
            bc1002_ns::bc1002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);
        // (i=0 loop continues below)

        if (ibcinfo(nbl,nseg,1,1) == 1003)
            bc1003_ns::bc1003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, x_r, z_r, cl,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,1) == 1005) {
            int inormmom = 0;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (ibcinfo(nbl,nseg,1,1) == 1006) {
            int inormmom = 1;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (ibcinfo(nbl,nseg,1,1) == 1008)
            bc1008_ns::bc1008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (ibcinfo(nbl,nseg,1,1) == 1011) {
            if (isklt_sav > 0)
                chksym_ns::chksym(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    isym, jsym, ksym, nou, bou, nbuf, ibufdim, myid);
            bc1011_ns::bc1011(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, isym, jsym, ksym, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (ibcinfo(nbl,nseg,1,1) == 1012) {
            if (isklt_sav > 0)
                chkrap_ns::chkrap(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    iwrap, jwrap, kwrap, nou, bou, nbuf, ibufdim, myid);
            bc1012_ns::bc1012(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iwrap, jwrap, kwrap, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (ibcinfo(nbl,nseg,1,1) == 1013)
            bc1013_ns::bc1013(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (ibcinfo(nbl,nseg,1,1) == 2002)
            bc2002_ns::bc2002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,1) == 2102) {
            int iflag = 0;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (ibcinfo(nbl,nseg,1,1) == 2103) {
            int iflag = 1;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (ibcinfo(nbl,nseg,1,1) == 2003)
            bc2003_ns::bc2003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,1) == 2009)
            bc2009_ns::bc2009(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,1) == 2010)
            bc2010_ns::bc2010(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,1) == 2019)
            bc2019_ns::bc2019(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);
        // (i=0 loop continues)

        if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2004 ||
            std::abs(ibcinfo(nbl,nseg,1,1)) == 2014) {
            int irelv = ibcinfo(nbl,nseg,1,1);
            int iuse3 = 0;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2024) {
            int irelv = ibcinfo(nbl,nseg,1,1);
            int iuse3 = 2;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2034) {
            int irelv = ibcinfo(nbl,nseg,1,1);
            int iuse3 = 1;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (ibcinfo(nbl,nseg,1,1) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ichk = 1;
            ngnew = 0;
            if (isklt_sav > 0 && ntime == 1) ichk = 0;
            if (nblnum > 0) {
                ngnew = nblg(nblnum) + (nbl - nblg(ngh));
                if (ngnew < nbl) ichk = 1;
            }
            int jdnew = jdimg(ngnew), kdnew = kdimg(ngnew), idnew = idimg(ngnew);
            FortranArray4DRef<double> qj0c_r(&w(lw(2,ngnew)), kdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qk0c_r(&w(lw(3,ngnew)), jdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qi0c_r(&w(lw(4,ngnew)), jdnew, kdnew, 5, 4);
            bc2006_ns::bc2006(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                x_r, y_r, z_r, ngnew, jdnew, kdnew, idnew,
                qj0c_r, qk0c_r, qi0c_r, ichk,
                nou, bou, nbuf, ibufdim, myid, nummem);
        }

        if (ibcinfo(nbl,nseg,1,1) == 2007)
            bc2007_ns::bc2007(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,1) == 2008 ||
            ibcinfo(nbl,nseg,1,1) == 2018 ||
            ibcinfo(nbl,nseg,1,1) == 2028 ||
            ibcinfo(nbl,nseg,1,1) == 2038) {
            int iflag = ibcinfo(nbl,nseg,1,1);
            bc2008_ns::bc2008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }

        if (ibcinfo(nbl,nseg,1,1) == 2016) {
            int irelv = ibcinfo(nbl,nseg,1,1);
            bc2016_ns::bc2016(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }

        if (ibcinfo(nbl,nseg,1,1) == 2026) {
            int irelv = ibcinfo(nbl,nseg,1,1);
            bc2026_ns::bc2026(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }
    } // end do 802 (i=0 boundary)


    // ********************
    // i=idim boundary
    // ********************
    for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
        ista = idim;
        iend = idim;
        jsta = ibcinfo(nbl,nseg,2,2);
        jend = ibcinfo(nbl,nseg,3,2);
        ksta = ibcinfo(nbl,nseg,4,2);
        kend = ibcinfo(nbl,nseg,5,2);
        ldata = lwdat(nbl,nseg,2);
        mdim = jend - jsta;
        ndim = kend - ksta;
        std::memcpy(filname, bcfiles(bcfilei(nbl,nseg,2)), 80);

        int jd = jdim, kd = kdim, id = idim;
        FortranArray4DRef<double> q_r(&w(lq), jd, kd, id, 5);
        FortranArray4DRef<double> qj0_r(&w(lqj0), kd, id-1, 5, 4);
        FortranArray4DRef<double> qk0_r(&w(lqk0), jd, id-1, 5, 4);
        FortranArray4DRef<double> qi0_r(&w(lqi0), jd, kd, 5, 4);
        FortranArray4DRef<double> sj_r(&w(lsj), jd, kd, id-1, 5);
        FortranArray4DRef<double> sk_r(&w(lsk), jd, kd, id-1, 5);
        FortranArray4DRef<double> si_r(&w(lsi), jd, kd, id, 5);
        FortranArray3DRef<double> bcj_r(&w(lbcj), kd, id-1, 2);
        FortranArray3DRef<double> bck_r(&w(lbck), jd, id-1, 2);
        FortranArray3DRef<double> bci_r(&w(lbci), jd, kd, 2);
        FortranArray4DRef<double> xtbj_r(&w(lxtbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> xtbk_r(&w(lxtbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> xtbi_r(&w(lxtbi), jd, kd, 3, 2);
        FortranArray4DRef<double> atbj_r(&w(latbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> atbk_r(&w(latbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> atbi_r(&w(latbi), jd, kd, 3, 2);
        FortranArray4DRef<double> tursav_r(&w(lxib), jd, kd, id, nummem);
        FortranArray4DRef<double> tj0_r(&w(ltj0), kd, id-1, nummem, 4);
        FortranArray4DRef<double> tk0_r(&w(ltk0), jd, id-1, nummem, 4);
        FortranArray4DRef<double> ti0_r(&w(lti0), jd, kd, nummem, 4);
        FortranArray3DRef<double> vist3d_r(&w(lvis), jd, kd, id);
        FortranArray4DRef<double> vj0_r(&w(lvj0), kd, id-1, 1, 4);
        FortranArray4DRef<double> vk0_r(&w(lvk0), jd, id-1, 1, 4);
        FortranArray4DRef<double> vi0_r(&w(lvi0), jd, kd, 1, 4);
        FortranArray3DRef<double> smin_r(&w(lsnk0), jd-1, kd-1, id-1);
        FortranArray3DRef<double> snj0_r(&w(lsni0), kd, id-1, 2);
        FortranArray3DRef<double> snk0_r(&w(lsnk0), jd, id-1, 2);
        FortranArray3DRef<double> sni0_r(&w(lsni0), jd, kd, 2);
        FortranArray3DRef<double> snjm_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> snkm_r(&w(lxkb), jd, kd, id);
        FortranArray3DRef<double> snim_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> x_r(&w(lx), jd, kd, id);
        FortranArray3DRef<double> y_r(&w(ly), jd, kd, id);
        FortranArray3DRef<double> z_r(&w(lz), jd, kd, id);
        FortranArray4DRef<double> bcdata_r(&w(ldata), mdim, ndim, 2, 12);
        FortranArray2DRef<double> wkt_r(&wk(ixwk), jvdim, 23);
        int nface_val = 2;

        if (ibcinfo(nbl,nseg,1,2) == 9999)
            bc9999_ns::bc9999(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, y_r, z_r, smin_r);

        if (ibcinfo(nbl,nseg,1,2) == 1000)
            bc1000_ns::bc1000(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (ibcinfo(nbl,nseg,1,2) == 1000 && noninflag > 0)
            bcnonin_ns::bcnonin(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                iuns, nou, bou, nbuf, ibufdim, x_r, y_r, z_r, nbl);

        if (ibcinfo(nbl,nseg,1,2) == 1001)
            bc1001_ns::bc1001(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, z_r);

        if (ibcinfo(nbl,nseg,1,2) == 1002)
            bc1002_ns::bc1002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (ibcinfo(nbl,nseg,1,2) == 1003)
            bc1003_ns::bc1003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, x_r, z_r, cl,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,2) == 1005) {
            int inormmom = 0;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (ibcinfo(nbl,nseg,1,2) == 1006) {
            int inormmom = 1;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (ibcinfo(nbl,nseg,1,2) == 1008)
            bc1008_ns::bc1008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);
        // (i=idim loop continues)

        if (ibcinfo(nbl,nseg,1,2) == 1011) {
            if (isklt_sav > 0)
                chksym_ns::chksym(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    isym, jsym, ksym, nou, bou, nbuf, ibufdim, myid);
            bc1011_ns::bc1011(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, isym, jsym, ksym, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (ibcinfo(nbl,nseg,1,2) == 1012) {
            if (isklt_sav > 0)
                chkrap_ns::chkrap(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    iwrap, jwrap, kwrap, nou, bou, nbuf, ibufdim, myid);
            bc1012_ns::bc1012(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iwrap, jwrap, kwrap, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (ibcinfo(nbl,nseg,1,2) == 1013)
            bc1013_ns::bc1013(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (ibcinfo(nbl,nseg,1,2) == 2002)
            bc2002_ns::bc2002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,2) == 2102) {
            int iflag = 0;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (ibcinfo(nbl,nseg,1,2) == 2103) {
            int iflag = 1;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (ibcinfo(nbl,nseg,1,2) == 2003)
            bc2003_ns::bc2003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,2) == 2009)
            bc2009_ns::bc2009(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,2) == 2010)
            bc2010_ns::bc2010(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,2) == 2019)
            bc2019_ns::bc2019(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);
        // (i=idim loop continues)

        if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2004 ||
            std::abs(ibcinfo(nbl,nseg,1,2)) == 2014) {
            int irelv = ibcinfo(nbl,nseg,1,2);
            int iuse3 = 0;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2024) {
            int irelv = ibcinfo(nbl,nseg,1,2);
            int iuse3 = 2;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2034) {
            int irelv = ibcinfo(nbl,nseg,1,2);
            int iuse3 = 1;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (ibcinfo(nbl,nseg,1,2) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ichk = 1;
            ngnew = 0;
            if (isklt_sav > 0 && ntime == 1) ichk = 0;
            if (nblnum > 0) {
                ngnew = nblg(nblnum) + (nbl - nblg(ngh));
                if (ngnew < nbl) ichk = 1;
            }
            int jdnew = jdimg(ngnew), kdnew = kdimg(ngnew), idnew = idimg(ngnew);
            FortranArray4DRef<double> qj0c_r(&w(lw(2,ngnew)), kdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qk0c_r(&w(lw(3,ngnew)), jdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qi0c_r(&w(lw(4,ngnew)), jdnew, kdnew, 5, 4);
            bc2006_ns::bc2006(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                x_r, y_r, z_r, ngnew, jdnew, kdnew, idnew,
                qj0c_r, qk0c_r, qi0c_r, ichk,
                nou, bou, nbuf, ibufdim, myid, nummem);
        }

        if (ibcinfo(nbl,nseg,1,2) == 2007)
            bc2007_ns::bc2007(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (ibcinfo(nbl,nseg,1,2) == 2008 ||
            ibcinfo(nbl,nseg,1,2) == 2018 ||
            ibcinfo(nbl,nseg,1,2) == 2028 ||
            ibcinfo(nbl,nseg,1,2) == 2038) {
            int iflag = ibcinfo(nbl,nseg,1,2);
            bc2008_ns::bc2008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }

        if (ibcinfo(nbl,nseg,1,2) == 2016) {
            int irelv = ibcinfo(nbl,nseg,1,2);
            bc2016_ns::bc2016(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }

        if (ibcinfo(nbl,nseg,1,2) == 2026) {
            int irelv = ibcinfo(nbl,nseg,1,2);
            bc2026_ns::bc2026(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }
    } // end do 803 (i=idim boundary)


    // ********************
    // j=0 boundary
    // ********************
    for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
        jsta = 1;
        jend = 1;
        ista = jbcinfo(nbl,nseg,2,1);
        iend = jbcinfo(nbl,nseg,3,1);
        ksta = jbcinfo(nbl,nseg,4,1);
        kend = jbcinfo(nbl,nseg,5,1);
        ldata = lwdat(nbl,nseg,3);
        mdim = kend - ksta;
        ndim = iend - ista;
        std::memcpy(filname, bcfiles(bcfilej(nbl,nseg,1)), 80);

        int jd = jdim, kd = kdim, id = idim;
        FortranArray4DRef<double> q_r(&w(lq), jd, kd, id, 5);
        FortranArray4DRef<double> qj0_r(&w(lqj0), kd, id-1, 5, 4);
        FortranArray4DRef<double> qk0_r(&w(lqk0), jd, id-1, 5, 4);
        FortranArray4DRef<double> qi0_r(&w(lqi0), jd, kd, 5, 4);
        FortranArray4DRef<double> sj_r(&w(lsj), jd, kd, id-1, 5);
        FortranArray4DRef<double> sk_r(&w(lsk), jd, kd, id-1, 5);
        FortranArray4DRef<double> si_r(&w(lsi), jd, kd, id, 5);
        FortranArray3DRef<double> bcj_r(&w(lbcj), kd, id-1, 2);
        FortranArray3DRef<double> bck_r(&w(lbck), jd, id-1, 2);
        FortranArray3DRef<double> bci_r(&w(lbci), jd, kd, 2);
        FortranArray4DRef<double> xtbj_r(&w(lxtbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> xtbk_r(&w(lxtbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> xtbi_r(&w(lxtbi), jd, kd, 3, 2);
        FortranArray4DRef<double> atbj_r(&w(latbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> atbk_r(&w(latbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> atbi_r(&w(latbi), jd, kd, 3, 2);
        FortranArray4DRef<double> tursav_r(&w(lxib), jd, kd, id, nummem);
        FortranArray4DRef<double> tj0_r(&w(ltj0), kd, id-1, nummem, 4);
        FortranArray4DRef<double> tk0_r(&w(ltk0), jd, id-1, nummem, 4);
        FortranArray4DRef<double> ti0_r(&w(lti0), jd, kd, nummem, 4);
        FortranArray3DRef<double> vist3d_r(&w(lvis), jd, kd, id);
        FortranArray4DRef<double> vj0_r(&w(lvj0), kd, id-1, 1, 4);
        FortranArray4DRef<double> vk0_r(&w(lvk0), jd, id-1, 1, 4);
        FortranArray4DRef<double> vi0_r(&w(lvi0), jd, kd, 1, 4);
        FortranArray3DRef<double> smin_r(&w(lsnk0), jd-1, kd-1, id-1);
        FortranArray3DRef<double> snj0_r(&w(lsni0), kd, id-1, 2);
        FortranArray3DRef<double> snk0_r(&w(lsnk0), jd, id-1, 2);
        FortranArray3DRef<double> sni0_r(&w(lsni0), jd, kd, 2);
        FortranArray3DRef<double> snjm_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> snkm_r(&w(lxkb), jd, kd, id);
        FortranArray3DRef<double> snim_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> x_r(&w(lx), jd, kd, id);
        FortranArray3DRef<double> y_r(&w(ly), jd, kd, id);
        FortranArray3DRef<double> z_r(&w(lz), jd, kd, id);
        FortranArray4DRef<double> bcdata_r(&w(ldata), mdim, ndim, 2, 12);
        FortranArray2DRef<double> wkt_r(&wk(ixwk), jvdim, 23);
        int nface_val = 3;

        if (jbcinfo(nbl,nseg,1,1) == 9999)
            bc9999_ns::bc9999(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, y_r, z_r, smin_r);

        if (jbcinfo(nbl,nseg,1,1) == 1000)
            bc1000_ns::bc1000(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (jbcinfo(nbl,nseg,1,1) == 1000 && noninflag > 0)
            bcnonin_ns::bcnonin(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                iuns, nou, bou, nbuf, ibufdim, x_r, y_r, z_r, nbl);

        if (jbcinfo(nbl,nseg,1,1) == 1001)
            bc1001_ns::bc1001(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, z_r);

        if (jbcinfo(nbl,nseg,1,1) == 1002)
            bc1002_ns::bc1002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (jbcinfo(nbl,nseg,1,1) == 1003)
            bc1003_ns::bc1003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, x_r, z_r, cl,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,1) == 1005) {
            int inormmom = 0;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (jbcinfo(nbl,nseg,1,1) == 1006) {
            int inormmom = 1;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (jbcinfo(nbl,nseg,1,1) == 1008)
            bc1008_ns::bc1008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);
        // (j=0 loop continues)

        if (jbcinfo(nbl,nseg,1,1) == 1011) {
            if (isklt_sav > 0)
                chksym_ns::chksym(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    isym, jsym, ksym, nou, bou, nbuf, ibufdim, myid);
            bc1011_ns::bc1011(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, isym, jsym, ksym, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (jbcinfo(nbl,nseg,1,1) == 1012) {
            if (isklt_sav > 0)
                chkrap_ns::chkrap(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    iwrap, jwrap, kwrap, nou, bou, nbuf, ibufdim, myid);
            bc1012_ns::bc1012(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iwrap, jwrap, kwrap, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (jbcinfo(nbl,nseg,1,1) == 1013)
            bc1013_ns::bc1013(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (jbcinfo(nbl,nseg,1,1) == 2002)
            bc2002_ns::bc2002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,1) == 2102) {
            int iflag = 0;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (jbcinfo(nbl,nseg,1,1) == 2103) {
            int iflag = 1;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (jbcinfo(nbl,nseg,1,1) == 2003)
            bc2003_ns::bc2003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,1) == 2009)
            bc2009_ns::bc2009(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,1) == 2010)
            bc2010_ns::bc2010(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,1) == 2019)
            bc2019_ns::bc2019(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2004 ||
            std::abs(jbcinfo(nbl,nseg,1,1)) == 2014) {
            int irelv = jbcinfo(nbl,nseg,1,1);
            int iuse3 = 0;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }
        // (j=0 loop continues)

        if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2024) {
            int irelv = jbcinfo(nbl,nseg,1,1);
            int iuse3 = 2;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2034) {
            int irelv = jbcinfo(nbl,nseg,1,1);
            int iuse3 = 1;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (jbcinfo(nbl,nseg,1,1) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ichk = 1;
            ngnew = 0;
            if (isklt_sav > 0 && ntime == 1) ichk = 0;
            if (nblnum > 0) {
                ngnew = nblg(nblnum) + (nbl - nblg(ngh));
                if (ngnew < nbl) ichk = 1;
            }
            int jdnew = jdimg(ngnew), kdnew = kdimg(ngnew), idnew = idimg(ngnew);
            FortranArray4DRef<double> qj0c_r(&w(lw(2,ngnew)), kdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qk0c_r(&w(lw(3,ngnew)), jdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qi0c_r(&w(lw(4,ngnew)), jdnew, kdnew, 5, 4);
            bc2006_ns::bc2006(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                x_r, y_r, z_r, ngnew, jdnew, kdnew, idnew,
                qj0c_r, qk0c_r, qi0c_r, ichk,
                nou, bou, nbuf, ibufdim, myid, nummem);
        }

        if (jbcinfo(nbl,nseg,1,1) == 2007)
            bc2007_ns::bc2007(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,1) == 2008 ||
            jbcinfo(nbl,nseg,1,1) == 2018 ||
            jbcinfo(nbl,nseg,1,1) == 2028 ||
            jbcinfo(nbl,nseg,1,1) == 2038) {
            int iflag = jbcinfo(nbl,nseg,1,1);
            bc2008_ns::bc2008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }

        if (jbcinfo(nbl,nseg,1,1) == 2016) {
            int irelv = jbcinfo(nbl,nseg,1,1);
            bc2016_ns::bc2016(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }

        if (jbcinfo(nbl,nseg,1,1) == 2026) {
            int irelv = jbcinfo(nbl,nseg,1,1);
            bc2026_ns::bc2026(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }
    } // end do 804 (j=0 boundary)


    // ********************
    // j=jdim boundary
    // ********************
    for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
        jsta = jdim;
        jend = jdim;
        ista = jbcinfo(nbl,nseg,2,2);
        iend = jbcinfo(nbl,nseg,3,2);
        ksta = jbcinfo(nbl,nseg,4,2);
        kend = jbcinfo(nbl,nseg,5,2);
        ldata = lwdat(nbl,nseg,4);
        mdim = kend - ksta;
        ndim = iend - ista;
        std::memcpy(filname, bcfiles(bcfilej(nbl,nseg,2)), 80);

        int jd = jdim, kd = kdim, id = idim;
        FortranArray4DRef<double> q_r(&w(lq), jd, kd, id, 5);
        FortranArray4DRef<double> qj0_r(&w(lqj0), kd, id-1, 5, 4);
        FortranArray4DRef<double> qk0_r(&w(lqk0), jd, id-1, 5, 4);
        FortranArray4DRef<double> qi0_r(&w(lqi0), jd, kd, 5, 4);
        FortranArray4DRef<double> sj_r(&w(lsj), jd, kd, id-1, 5);
        FortranArray4DRef<double> sk_r(&w(lsk), jd, kd, id-1, 5);
        FortranArray4DRef<double> si_r(&w(lsi), jd, kd, id, 5);
        FortranArray3DRef<double> bcj_r(&w(lbcj), kd, id-1, 2);
        FortranArray3DRef<double> bck_r(&w(lbck), jd, id-1, 2);
        FortranArray3DRef<double> bci_r(&w(lbci), jd, kd, 2);
        FortranArray4DRef<double> xtbj_r(&w(lxtbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> xtbk_r(&w(lxtbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> xtbi_r(&w(lxtbi), jd, kd, 3, 2);
        FortranArray4DRef<double> atbj_r(&w(latbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> atbk_r(&w(latbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> atbi_r(&w(latbi), jd, kd, 3, 2);
        FortranArray4DRef<double> tursav_r(&w(lxib), jd, kd, id, nummem);
        FortranArray4DRef<double> tj0_r(&w(ltj0), kd, id-1, nummem, 4);
        FortranArray4DRef<double> tk0_r(&w(ltk0), jd, id-1, nummem, 4);
        FortranArray4DRef<double> ti0_r(&w(lti0), jd, kd, nummem, 4);
        FortranArray3DRef<double> vist3d_r(&w(lvis), jd, kd, id);
        FortranArray4DRef<double> vj0_r(&w(lvj0), kd, id-1, 1, 4);
        FortranArray4DRef<double> vk0_r(&w(lvk0), jd, id-1, 1, 4);
        FortranArray4DRef<double> vi0_r(&w(lvi0), jd, kd, 1, 4);
        FortranArray3DRef<double> smin_r(&w(lsnk0), jd-1, kd-1, id-1);
        FortranArray3DRef<double> snj0_r(&w(lsni0), kd, id-1, 2);
        FortranArray3DRef<double> snk0_r(&w(lsnk0), jd, id-1, 2);
        FortranArray3DRef<double> sni0_r(&w(lsni0), jd, kd, 2);
        FortranArray3DRef<double> snjm_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> snkm_r(&w(lxkb), jd, kd, id);
        FortranArray3DRef<double> snim_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> x_r(&w(lx), jd, kd, id);
        FortranArray3DRef<double> y_r(&w(ly), jd, kd, id);
        FortranArray3DRef<double> z_r(&w(lz), jd, kd, id);
        FortranArray4DRef<double> bcdata_r(&w(ldata), mdim, ndim, 2, 12);
        FortranArray2DRef<double> wkt_r(&wk(ixwk), jvdim, 23);
        int nface_val = 4;

        if (jbcinfo(nbl,nseg,1,2) == 9999)
            bc9999_ns::bc9999(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, y_r, z_r, smin_r);

        if (jbcinfo(nbl,nseg,1,2) == 1000)
            bc1000_ns::bc1000(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (jbcinfo(nbl,nseg,1,2) == 1000 && noninflag > 0)
            bcnonin_ns::bcnonin(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                iuns, nou, bou, nbuf, ibufdim, x_r, y_r, z_r, nbl);

        if (jbcinfo(nbl,nseg,1,2) == 1001)
            bc1001_ns::bc1001(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, z_r);

        if (jbcinfo(nbl,nseg,1,2) == 1002)
            bc1002_ns::bc1002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (jbcinfo(nbl,nseg,1,2) == 1003)
            bc1003_ns::bc1003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, x_r, z_r, cl,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,2) == 1005) {
            int inormmom = 0;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (jbcinfo(nbl,nseg,1,2) == 1006) {
            int inormmom = 1;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (jbcinfo(nbl,nseg,1,2) == 1008)
            bc1008_ns::bc1008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);
        // (j=jdim loop continues)

        if (jbcinfo(nbl,nseg,1,2) == 1011) {
            if (isklt_sav > 0)
                chksym_ns::chksym(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    isym, jsym, ksym, nou, bou, nbuf, ibufdim, myid);
            bc1011_ns::bc1011(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, isym, jsym, ksym, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (jbcinfo(nbl,nseg,1,2) == 1012) {
            if (isklt_sav > 0)
                chkrap_ns::chkrap(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    iwrap, jwrap, kwrap, nou, bou, nbuf, ibufdim, myid);
            bc1012_ns::bc1012(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iwrap, jwrap, kwrap, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (jbcinfo(nbl,nseg,1,2) == 1013)
            bc1013_ns::bc1013(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (jbcinfo(nbl,nseg,1,2) == 2002)
            bc2002_ns::bc2002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,2) == 2102) {
            int iflag = 0;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (jbcinfo(nbl,nseg,1,2) == 2103) {
            int iflag = 1;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (jbcinfo(nbl,nseg,1,2) == 2003)
            bc2003_ns::bc2003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,2) == 2009)
            bc2009_ns::bc2009(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,2) == 2010)
            bc2010_ns::bc2010(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,2) == 2019)
            bc2019_ns::bc2019(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2004 ||
            std::abs(jbcinfo(nbl,nseg,1,2)) == 2014) {
            int irelv = jbcinfo(nbl,nseg,1,2);
            int iuse3 = 0;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }
        // (j=jdim loop continues)

        if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2024) {
            int irelv = jbcinfo(nbl,nseg,1,2);
            int iuse3 = 2;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2034) {
            int irelv = jbcinfo(nbl,nseg,1,2);
            int iuse3 = 1;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (jbcinfo(nbl,nseg,1,2) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ichk = 1;
            ngnew = 0;
            if (isklt_sav > 0 && ntime == 1) ichk = 0;
            if (nblnum > 0) {
                ngnew = nblg(nblnum) + (nbl - nblg(ngh));
                if (ngnew < nbl) ichk = 1;
            }
            int jdnew = jdimg(ngnew), kdnew = kdimg(ngnew), idnew = idimg(ngnew);
            FortranArray4DRef<double> qj0c_r(&w(lw(2,ngnew)), kdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qk0c_r(&w(lw(3,ngnew)), jdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qi0c_r(&w(lw(4,ngnew)), jdnew, kdnew, 5, 4);
            bc2006_ns::bc2006(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                x_r, y_r, z_r, ngnew, jdnew, kdnew, idnew,
                qj0c_r, qk0c_r, qi0c_r, ichk,
                nou, bou, nbuf, ibufdim, myid, nummem);
        }

        if (jbcinfo(nbl,nseg,1,2) == 2007)
            bc2007_ns::bc2007(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (jbcinfo(nbl,nseg,1,2) == 2008 ||
            jbcinfo(nbl,nseg,1,2) == 2018 ||
            jbcinfo(nbl,nseg,1,2) == 2028 ||
            jbcinfo(nbl,nseg,1,2) == 2038) {
            int iflag = jbcinfo(nbl,nseg,1,2);
            bc2008_ns::bc2008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }

        if (jbcinfo(nbl,nseg,1,2) == 2016) {
            int irelv = jbcinfo(nbl,nseg,1,2);
            bc2016_ns::bc2016(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }

        if (jbcinfo(nbl,nseg,1,2) == 2026) {
            int irelv = jbcinfo(nbl,nseg,1,2);
            bc2026_ns::bc2026(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }
    } // end do 805 (j=jdim boundary)


    // ********************
    // k=0 boundary
    // ********************
    for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
        ksta = 1;
        kend = 1;
        ista = kbcinfo(nbl,nseg,2,1);
        iend = kbcinfo(nbl,nseg,3,1);
        jsta = kbcinfo(nbl,nseg,4,1);
        jend = kbcinfo(nbl,nseg,5,1);
        ldata = lwdat(nbl,nseg,5);
        mdim = jend - jsta;
        ndim = iend - ista;
        std::memcpy(filname, bcfiles(bcfilek(nbl,nseg,1)), 80);

        int jd = jdim, kd = kdim, id = idim;
        FortranArray4DRef<double> q_r(&w(lq), jd, kd, id, 5);
        FortranArray4DRef<double> qj0_r(&w(lqj0), kd, id-1, 5, 4);
        FortranArray4DRef<double> qk0_r(&w(lqk0), jd, id-1, 5, 4);
        FortranArray4DRef<double> qi0_r(&w(lqi0), jd, kd, 5, 4);
        FortranArray4DRef<double> sj_r(&w(lsj), jd, kd, id-1, 5);
        FortranArray4DRef<double> sk_r(&w(lsk), jd, kd, id-1, 5);
        FortranArray4DRef<double> si_r(&w(lsi), jd, kd, id, 5);
        FortranArray3DRef<double> bcj_r(&w(lbcj), kd, id-1, 2);
        FortranArray3DRef<double> bck_r(&w(lbck), jd, id-1, 2);
        FortranArray3DRef<double> bci_r(&w(lbci), jd, kd, 2);
        FortranArray4DRef<double> xtbj_r(&w(lxtbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> xtbk_r(&w(lxtbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> xtbi_r(&w(lxtbi), jd, kd, 3, 2);
        FortranArray4DRef<double> atbj_r(&w(latbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> atbk_r(&w(latbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> atbi_r(&w(latbi), jd, kd, 3, 2);
        FortranArray4DRef<double> tursav_r(&w(lxib), jd, kd, id, nummem);
        FortranArray4DRef<double> tj0_r(&w(ltj0), kd, id-1, nummem, 4);
        FortranArray4DRef<double> tk0_r(&w(ltk0), jd, id-1, nummem, 4);
        FortranArray4DRef<double> ti0_r(&w(lti0), jd, kd, nummem, 4);
        FortranArray3DRef<double> vist3d_r(&w(lvis), jd, kd, id);
        FortranArray4DRef<double> vj0_r(&w(lvj0), kd, id-1, 1, 4);
        FortranArray4DRef<double> vk0_r(&w(lvk0), jd, id-1, 1, 4);
        FortranArray4DRef<double> vi0_r(&w(lvi0), jd, kd, 1, 4);
        FortranArray3DRef<double> smin_r(&w(lsnk0), jd-1, kd-1, id-1);
        FortranArray3DRef<double> snj0_r(&w(lsni0), kd, id-1, 2);
        FortranArray3DRef<double> snk0_r(&w(lsnk0), jd, id-1, 2);
        FortranArray3DRef<double> sni0_r(&w(lsni0), jd, kd, 2);
        FortranArray3DRef<double> snjm_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> snkm_r(&w(lxkb), jd, kd, id);
        FortranArray3DRef<double> snim_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> x_r(&w(lx), jd, kd, id);
        FortranArray3DRef<double> y_r(&w(ly), jd, kd, id);
        FortranArray3DRef<double> z_r(&w(lz), jd, kd, id);
        FortranArray4DRef<double> bcdata_r(&w(ldata), mdim, ndim, 2, 12);
        FortranArray2DRef<double> wkt_r(&wk(ixwk), jvdim, 23);
        int nface_val = 5;

        if (kbcinfo(nbl,nseg,1,1) == 9999)
            bc9999_ns::bc9999(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, y_r, z_r, smin_r);

        if (kbcinfo(nbl,nseg,1,1) == 1000)
            bc1000_ns::bc1000(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (kbcinfo(nbl,nseg,1,1) == 1000 && noninflag > 0)
            bcnonin_ns::bcnonin(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                iuns, nou, bou, nbuf, ibufdim, x_r, y_r, z_r, nbl);

        if (kbcinfo(nbl,nseg,1,1) == 1001)
            bc1001_ns::bc1001(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, z_r);

        if (kbcinfo(nbl,nseg,1,1) == 1002)
            bc1002_ns::bc1002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (kbcinfo(nbl,nseg,1,1) == 1003)
            bc1003_ns::bc1003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, x_r, z_r, cl,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,1) == 1005) {
            int inormmom = 0;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (kbcinfo(nbl,nseg,1,1) == 1006) {
            int inormmom = 1;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (kbcinfo(nbl,nseg,1,1) == 1008)
            bc1008_ns::bc1008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);
        // (k=0 loop continues)

        if (kbcinfo(nbl,nseg,1,1) == 1011) {
            if (isklt_sav > 0)
                chksym_ns::chksym(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    isym, jsym, ksym, nou, bou, nbuf, ibufdim, myid);
            bc1011_ns::bc1011(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, isym, jsym, ksym, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (kbcinfo(nbl,nseg,1,1) == 1012) {
            if (isklt_sav > 0)
                chkrap_ns::chkrap(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    iwrap, jwrap, kwrap, nou, bou, nbuf, ibufdim, myid);
            bc1012_ns::bc1012(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iwrap, jwrap, kwrap, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (kbcinfo(nbl,nseg,1,1) == 1013)
            bc1013_ns::bc1013(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (kbcinfo(nbl,nseg,1,1) == 2002)
            bc2002_ns::bc2002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,1) == 2102) {
            int iflag = 0;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (kbcinfo(nbl,nseg,1,1) == 2103) {
            int iflag = 1;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (kbcinfo(nbl,nseg,1,1) == 2003)
            bc2003_ns::bc2003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,1) == 2009)
            bc2009_ns::bc2009(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,1) == 2010)
            bc2010_ns::bc2010(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,1) == 2019)
            bc2019_ns::bc2019(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2004 ||
            std::abs(kbcinfo(nbl,nseg,1,1)) == 2014) {
            int irelv = kbcinfo(nbl,nseg,1,1);
            int iuse3 = 0;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }
        // (k=0 loop continues)

        if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2024) {
            int irelv = kbcinfo(nbl,nseg,1,1);
            int iuse3 = 2;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2034) {
            int irelv = kbcinfo(nbl,nseg,1,1);
            int iuse3 = 1;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (kbcinfo(nbl,nseg,1,1) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ichk = 1;
            ngnew = 0;
            if (isklt_sav > 0 && ntime == 1) ichk = 0;
            if (nblnum > 0) {
                ngnew = nblg(nblnum) + (nbl - nblg(ngh));
                if (ngnew < nbl) ichk = 1;
            }
            int jdnew = jdimg(ngnew), kdnew = kdimg(ngnew), idnew = idimg(ngnew);
            FortranArray4DRef<double> qj0c_r(&w(lw(2,ngnew)), kdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qk0c_r(&w(lw(3,ngnew)), jdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qi0c_r(&w(lw(4,ngnew)), jdnew, kdnew, 5, 4);
            bc2006_ns::bc2006(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                x_r, y_r, z_r, ngnew, jdnew, kdnew, idnew,
                qj0c_r, qk0c_r, qi0c_r, ichk,
                nou, bou, nbuf, ibufdim, myid, nummem);
        }

        if (kbcinfo(nbl,nseg,1,1) == 2007)
            bc2007_ns::bc2007(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,1) == 2008 ||
            kbcinfo(nbl,nseg,1,1) == 2018 ||
            kbcinfo(nbl,nseg,1,1) == 2028 ||
            kbcinfo(nbl,nseg,1,1) == 2038) {
            int iflag = kbcinfo(nbl,nseg,1,1);
            bc2008_ns::bc2008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }

        if (kbcinfo(nbl,nseg,1,1) == 2016) {
            int irelv = kbcinfo(nbl,nseg,1,1);
            bc2016_ns::bc2016(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }

        if (kbcinfo(nbl,nseg,1,1) == 2026) {
            int irelv = kbcinfo(nbl,nseg,1,1);
            bc2026_ns::bc2026(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }
    } // end do 806 (k=0 boundary)


    // ********************
    // k=kdim boundary
    // ********************
    for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
        ksta = kdim;
        kend = kdim;
        ista = kbcinfo(nbl,nseg,2,2);
        iend = kbcinfo(nbl,nseg,3,2);
        jsta = kbcinfo(nbl,nseg,4,2);
        jend = kbcinfo(nbl,nseg,5,2);
        ldata = lwdat(nbl,nseg,6);
        mdim = jend - jsta;
        ndim = iend - ista;
        std::memcpy(filname, bcfiles(bcfilek(nbl,nseg,2)), 80);

        int jd = jdim, kd = kdim, id = idim;
        FortranArray4DRef<double> q_r(&w(lq), jd, kd, id, 5);
        FortranArray4DRef<double> qj0_r(&w(lqj0), kd, id-1, 5, 4);
        FortranArray4DRef<double> qk0_r(&w(lqk0), jd, id-1, 5, 4);
        FortranArray4DRef<double> qi0_r(&w(lqi0), jd, kd, 5, 4);
        FortranArray4DRef<double> sj_r(&w(lsj), jd, kd, id-1, 5);
        FortranArray4DRef<double> sk_r(&w(lsk), jd, kd, id-1, 5);
        FortranArray4DRef<double> si_r(&w(lsi), jd, kd, id, 5);
        FortranArray3DRef<double> bcj_r(&w(lbcj), kd, id-1, 2);
        FortranArray3DRef<double> bck_r(&w(lbck), jd, id-1, 2);
        FortranArray3DRef<double> bci_r(&w(lbci), jd, kd, 2);
        FortranArray4DRef<double> xtbj_r(&w(lxtbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> xtbk_r(&w(lxtbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> xtbi_r(&w(lxtbi), jd, kd, 3, 2);
        FortranArray4DRef<double> atbj_r(&w(latbj), kd, id-1, 3, 2);
        FortranArray4DRef<double> atbk_r(&w(latbk), jd, id-1, 3, 2);
        FortranArray4DRef<double> atbi_r(&w(latbi), jd, kd, 3, 2);
        FortranArray4DRef<double> tursav_r(&w(lxib), jd, kd, id, nummem);
        FortranArray4DRef<double> tj0_r(&w(ltj0), kd, id-1, nummem, 4);
        FortranArray4DRef<double> tk0_r(&w(ltk0), jd, id-1, nummem, 4);
        FortranArray4DRef<double> ti0_r(&w(lti0), jd, kd, nummem, 4);
        FortranArray3DRef<double> vist3d_r(&w(lvis), jd, kd, id);
        FortranArray4DRef<double> vj0_r(&w(lvj0), kd, id-1, 1, 4);
        FortranArray4DRef<double> vk0_r(&w(lvk0), jd, id-1, 1, 4);
        FortranArray4DRef<double> vi0_r(&w(lvi0), jd, kd, 1, 4);
        FortranArray3DRef<double> smin_r(&w(lsnk0), jd-1, kd-1, id-1);
        FortranArray3DRef<double> snj0_r(&w(lsni0), kd, id-1, 2);
        FortranArray3DRef<double> snk0_r(&w(lsnk0), jd, id-1, 2);
        FortranArray3DRef<double> sni0_r(&w(lsni0), jd, kd, 2);
        FortranArray3DRef<double> snjm_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> snkm_r(&w(lxkb), jd, kd, id);
        FortranArray3DRef<double> snim_r(&w(lnbl), jd, kd, id);
        FortranArray3DRef<double> x_r(&w(lx), jd, kd, id);
        FortranArray3DRef<double> y_r(&w(ly), jd, kd, id);
        FortranArray3DRef<double> z_r(&w(lz), jd, kd, id);
        FortranArray4DRef<double> bcdata_r(&w(ldata), mdim, ndim, 2, 12);
        FortranArray2DRef<double> wkt_r(&wk(ixwk), jvdim, 23);
        int nface_val = 6;

        if (kbcinfo(nbl,nseg,1,2) == 9999)
            bc9999_ns::bc9999(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, y_r, z_r, smin_r);

        if (kbcinfo(nbl,nseg,1,2) == 1000)
            bc1000_ns::bc1000(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (kbcinfo(nbl,nseg,1,2) == 1000 && noninflag > 0)
            bcnonin_ns::bcnonin(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                iuns, nou, bou, nbuf, ibufdim, x_r, y_r, z_r, nbl);

        if (kbcinfo(nbl,nseg,1,2) == 1001)
            bc1001_ns::bc1001(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem, x_r, z_r);

        if (kbcinfo(nbl,nseg,1,2) == 1002)
            bc1002_ns::bc1002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (kbcinfo(nbl,nseg,1,2) == 1003)
            bc1003_ns::bc1003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, x_r, z_r, cl,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,2) == 1005) {
            int inormmom = 0;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (kbcinfo(nbl,nseg,1,2) == 1006) {
            int inormmom = 1;
            bc1005_ns::bc1005(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, inormmom, nou, bou, nbuf, ibufdim,
                nummem);
        }

        if (kbcinfo(nbl,nseg,1,2) == 1008)
            bc1008_ns::bc1008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);
        // (k=kdim loop continues)

        if (kbcinfo(nbl,nseg,1,2) == 1011) {
            if (isklt_sav > 0)
                chksym_ns::chksym(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    isym, jsym, ksym, nou, bou, nbuf, ibufdim, myid);
            bc1011_ns::bc1011(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, isym, jsym, ksym, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (kbcinfo(nbl,nseg,1,2) == 1012) {
            if (isklt_sav > 0)
                chkrap_ns::chkrap(nbl, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, nface_val, idim, jdim, kdim,
                    maxbl, maxseg, ista, iend, jsta, jend, ksta, kend,
                    iwrap, jwrap, kwrap, nou, bou, nbuf, ibufdim, myid);
            bc1012_ns::bc1012(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iwrap, jwrap, kwrap, iuns,
                nou, bou, nbuf, ibufdim, nummem);
        }

        if (kbcinfo(nbl,nseg,1,2) == 1013)
            bc1013_ns::bc1013(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, iuns, nou, bou, nbuf, ibufdim,
                nummem);

        if (kbcinfo(nbl,nseg,1,2) == 2002)
            bc2002_ns::bc2002(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,2) == 2102) {
            int iflag = 0;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (kbcinfo(nbl,nseg,1,2) == 2103) {
            int iflag = 1;
            bc2102_ns::bc2102(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nbl, nou, bou, nbuf, ibufdim, myid, nummem, iflag);
        }

        if (kbcinfo(nbl,nseg,1,2) == 2003)
            bc2003_ns::bc2003(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, wkt_r, jvdim, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,2) == 2009)
            bc2009_ns::bc2009(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,2) == 2010)
            bc2010_ns::bc2010(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,2) == 2019)
            bc2019_ns::bc2019(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2004 ||
            std::abs(kbcinfo(nbl,nseg,1,2)) == 2014) {
            int irelv = kbcinfo(nbl,nseg,1,2);
            int iuse3 = 0;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }
        // (k=kdim loop continues)

        if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2024) {
            int irelv = kbcinfo(nbl,nseg,1,2);
            int iuse3 = 2;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2034) {
            int irelv = kbcinfo(nbl,nseg,1,2);
            int iuse3 = 1;
            bc2004_ns::bc2004(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem, iuse3);
        }

        if (kbcinfo(nbl,nseg,1,2) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ichk = 1;
            ngnew = 0;
            if (isklt_sav > 0 && ntime == 1) ichk = 0;
            if (nblnum > 0) {
                ngnew = nblg(nblnum) + (nbl - nblg(ngh));
                if (ngnew < nbl) ichk = 1;
            }
            int jdnew = jdimg(ngnew), kdnew = kdimg(ngnew), idnew = idimg(ngnew);
            FortranArray4DRef<double> qj0c_r(&w(lw(2,ngnew)), kdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qk0c_r(&w(lw(3,ngnew)), jdnew, idnew-1, 5, 4);
            FortranArray4DRef<double> qi0c_r(&w(lw(4,ngnew)), jdnew, kdnew, 5, 4);
            bc2006_ns::bc2006(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                x_r, y_r, z_r, ngnew, jdnew, kdnew, idnew,
                qj0c_r, qk0c_r, qi0c_r, ichk,
                nou, bou, nbuf, ibufdim, myid, nummem);
        }

        if (kbcinfo(nbl,nseg,1,2) == 2007)
            bc2007_ns::bc2007(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, nummem);

        if (kbcinfo(nbl,nseg,1,2) == 2008 ||
            kbcinfo(nbl,nseg,1,2) == 2018 ||
            kbcinfo(nbl,nseg,1,2) == 2028 ||
            kbcinfo(nbl,nseg,1,2) == 2038) {
            int iflag = kbcinfo(nbl,nseg,1,2);
            bc2008_ns::bc2008(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r, filname, iuns,
                nou, bou, nbuf, ibufdim, myid, iflag, nummem);
        }

        if (kbcinfo(nbl,nseg,1,2) == 2016) {
            int irelv = kbcinfo(nbl,nseg,1,2);
            bc2016_ns::bc2016(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }

        if (kbcinfo(nbl,nseg,1,2) == 2026) {
            int irelv = kbcinfo(nbl,nseg,1,2);
            bc2026_ns::bc2026(jd, kd, id, q_r, qj0_r, qk0_r, qi0_r,
                sj_r, sk_r, si_r, bcj_r, bck_r, bci_r,
                xtbj_r, xtbk_r, xtbi_r, atbj_r, atbk_r, atbi_r,
                ista, iend, jsta, jend, ksta, kend, nface_val,
                tursav_r, tj0_r, tk0_r, ti0_r,
                smin_r, vist3d_r, vj0_r, vk0_r, vi0_r, mdim, ndim, bcdata_r,
                filname, iuns, irelv, snj0_r, snk0_r, sni0_r,
                ntime, snjm_r, snkm_r, snim_r, nou, bou, nbuf, ibufdim, myid,
                nummem);
        }
    } // end do 807 (k=kdim boundary)

    // Restore isklton
    isklton = isklt_sav;

} // end bc

} // namespace bc_ns
