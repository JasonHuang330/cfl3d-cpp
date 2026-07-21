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
#include "bc_xmera.h"
#include "bc.h"
#include "avghole.h"
#include "runtime/fortran_array.h"
#include "runtime/common_blocks.h"

namespace bc_xmera_ns {

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
    bc_ns::bc(ntime, nbl, lw, lw2, w, mgwk, wk, nwork, cl, nou, bou, nbuf, ibufdim,
              maxbl, maxgr, maxseg, itrans, irotat, idefrm, igridg, nblg,
              nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
              ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek,
              lwdat, myid, idimg, jdimg, kdimg, bcfiles, mxbcfil, nummem);
}

void bc_xmera(int& ntime, int& nbl, FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
              FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
              int& nwork, int& maxbl, int& iitot, FortranArray2DRef<int> iviscg,
              FortranArray1DRef<int> iovrlp, FortranArray1DRef<int> lbg,
              FortranArray2DRef<int> ibpntsg, FortranArray3DRef<double> qb,
              FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg,
              FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg,
              FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
              int& nbuf, int& ibufdim, int& int_updt, int& nummem)
{
    // COMMON /ginfo/
    int& jdim  = cmn_ginfo.jdim;
    int& kdim  = cmn_ginfo.kdim;
    int& idim  = cmn_ginfo.idim;
    int& lq    = cmn_ginfo.lq;
    int& lqj0  = cmn_ginfo.lqj0;
    int& lqk0  = cmn_ginfo.lqk0;
    int& lqi0  = cmn_ginfo.lqi0;
    int& lblk  = cmn_ginfo.lblk;
    int& lbcj  = cmn_ginfo.lbcj;
    int& lbck  = cmn_ginfo.lbck;
    int& lbci  = cmn_ginfo.lbci;
    int& lvis  = cmn_ginfo.lvis;
    int& lvj0  = cmn_ginfo.lvj0;
    int& lvk0  = cmn_ginfo.lvk0;
    int& lvi0  = cmn_ginfo.lvi0;
    int& lxib  = cmn_ginfo.lxib;
    int& ltj0  = cmn_ginfo.ltj0;
    int& ltk0  = cmn_ginfo.ltk0;
    int& lti0  = cmn_ginfo.lti0;

    // chimera grid boundaries
    if (iovrlp(nbl) == 1 && ntime > 0) {
        int ldim = 5;
        {
            FortranArray4DRef<double> q_view(&w(lq), jdim, kdim, idim, ldim);
            FortranArray3DRef<double> blank_view(&w(lblk), jdim, kdim, idim);
            avghole_ns::avghole(q_view, blank_view, jdim, kdim, idim, nbl, ldim, int_updt);
        }
        {
            FortranArray4DRef<double> q_view(&w(lq), jdim, kdim, idim, ldim);
            FortranArray4DRef<double> qj0_view(&w(lqj0), kdim, idim-1, ldim, 4);
            FortranArray4DRef<double> qk0_view(&w(lqk0), jdim, idim-1, ldim, 4);
            FortranArray4DRef<double> qi0_view(&w(lqi0), jdim, kdim, ldim, 4);
            FortranArray3DRef<double> bcj_view(&w(lbcj), kdim, idim-1, 2);
            FortranArray3DRef<double> bck_view(&w(lbck), jdim, idim-1, 2);
            FortranArray3DRef<double> bci_view(&w(lbci), jdim, kdim, 2);
            avghole_ns::xupdt(q_view, qj0_view, qk0_view, qi0_view,
                              jdim, kdim, idim, nbl, ldim,
                              bcj_view, bck_view, bci_view,
                              maxbl, iitot, iibg, kkbg, jjbg, ibcg,
                              lbg, ibpntsg, qb, nou, bou, nbuf, ibufdim, int_updt);
        }

        // update turbulence quantities
        if (iviscg(nbl, 1) >= 2 || iviscg(nbl, 2) >= 2 || iviscg(nbl, 3) >= 2) {
            int ldim2 = 1;
            {
                FortranArray4DRef<double> vis_view(&w(lvis), jdim, kdim, idim, ldim2);
                FortranArray3DRef<double> blank_view(&w(lblk), jdim, kdim, idim);
                avghole_ns::avghole(vis_view, blank_view, jdim, kdim, idim, nbl, ldim2, int_updt);
            }
            {
                FortranArray4DRef<double> vis_view(&w(lvis), jdim, kdim, idim, ldim2);
                FortranArray4DRef<double> vj0_view(&w(lvj0), kdim, idim-1, ldim2, 4);
                FortranArray4DRef<double> vk0_view(&w(lvk0), jdim, idim-1, ldim2, 4);
                FortranArray4DRef<double> vi0_view(&w(lvi0), jdim, kdim, ldim2, 4);
                FortranArray3DRef<double> bcj_view(&w(lbcj), kdim, idim-1, 2);
                FortranArray3DRef<double> bck_view(&w(lbck), jdim, idim-1, 2);
                FortranArray3DRef<double> bci_view(&w(lbci), jdim, kdim, 2);
                avghole_ns::xupdt(vis_view, vj0_view, vk0_view, vi0_view,
                                  jdim, kdim, idim, nbl, ldim2,
                                  bcj_view, bck_view, bci_view,
                                  maxbl, iitot, iibg, kkbg, jjbg, ibcg,
                                  lbg, ibpntsg, qb, nou, bou, nbuf, ibufdim, int_updt);
            }
        }
        if (iviscg(nbl, 1) >= 4 || iviscg(nbl, 2) >= 4 || iviscg(nbl, 3) >= 4) {
            int ldim3 = nummem;
            {
                FortranArray4DRef<double> xib_view(&w(lxib), jdim, kdim, idim, ldim3);
                FortranArray3DRef<double> blank_view(&w(lblk), jdim, kdim, idim);
                avghole_ns::avghole(xib_view, blank_view, jdim, kdim, idim, nbl, ldim3, int_updt);
            }
            {
                FortranArray4DRef<double> xib_view(&w(lxib), jdim, kdim, idim, ldim3);
                FortranArray4DRef<double> tj0_view(&w(ltj0), kdim, idim-1, ldim3, 4);
                FortranArray4DRef<double> tk0_view(&w(ltk0), jdim, idim-1, ldim3, 4);
                FortranArray4DRef<double> ti0_view(&w(lti0), jdim, kdim, ldim3, 4);
                FortranArray3DRef<double> bcj_view(&w(lbcj), kdim, idim-1, 2);
                FortranArray3DRef<double> bck_view(&w(lbck), jdim, idim-1, 2);
                FortranArray3DRef<double> bci_view(&w(lbci), jdim, kdim, 2);
                avghole_ns::xupdt(xib_view, tj0_view, tk0_view, ti0_view,
                                  jdim, kdim, idim, nbl, ldim3,
                                  bcj_view, bck_view, bci_view,
                                  maxbl, iitot, iibg, kkbg, jjbg, ibcg,
                                  lbg, ibpntsg, qb, nou, bou, nbuf, ibufdim, int_updt);
            }
        }
    }
    return;
}

} // namespace bc_xmera_ns
