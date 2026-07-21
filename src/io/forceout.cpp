// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "forceout.h"
#include "force.h"
#include "runtime/fortran_io.h"
#include <cstdio>

namespace forceout_ns {

void force(int& jdim, int& kdim, int& idim, FortranArray3DRef<double> x, FortranArray3DRef<double> y, FortranArray3DRef<double> z, FortranArray4DRef<double> sk, FortranArray4DRef<double> sj, FortranArray4DRef<double> si, double& cl, double& cd, double& cz, double& cy, double& cx, double& cmy, double& cmx, double& cmz, double& chd, double& swet, int& i00, FortranArray3DRef<double> ub, FortranArray3DRef<double> vb, FortranArray3DRef<double> wb, FortranArray3DRef<double> vmuk, FortranArray3DRef<double> vmuj, FortranArray3DRef<double> vmui, FortranArray3DRef<double> vol, int& ifo, int& jfo, int& kfo, FortranArray3DRef<double> bcj, FortranArray3DRef<double> bck, FortranArray3DRef<double> bci, FortranArray3DRef<double> blank, int& nbl, FortranArray4DRef<double> xtbj, FortranArray4DRef<double> xtbk, FortranArray4DRef<double> xtbi, int& iuns, FortranArray4DRef<double> qj0, FortranArray4DRef<double> qk0, FortranArray4DRef<double> qi0, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, int& nn, int& maxseg)
{
    force_ns::force(jdim, kdim, idim, x, y, z, sk, sj, si, cl, cd, cz, cy, cx, cmy, cmx, cmz, chd, swet, i00, ub, vb, wb, vmuk, vmuj, vmui, vol, ifo, jfo, kfo, bcj, bck, bci, blank, nbl, xtbj, xtbk, xtbi, iuns, qj0, qk0, qi0, nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim, ibcinfo, jbcinfo, kbcinfo, nn, maxseg);
}

void forceout(int& iseq, int& maxbl, int& maxgr, int& maxseg, int& nblock, FortranArray1DRef<int> iforce, FortranArray1DRef<int> igridg, FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim, FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim, FortranArray1DRef<int> levelg, FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo, FortranArray4DRef<int> kbcinfo, FortranArray1DRef<double> swett, FortranArray1DRef<double> clt, FortranArray1DRef<double> cdt, FortranArray1DRef<double> cxt, FortranArray1DRef<double> cyt, FortranArray1DRef<double> czt, FortranArray1DRef<double> cmxt, FortranArray1DRef<double> cmyt, FortranArray1DRef<double> cmzt, FortranArray1DRef<double> cdpt, FortranArray1DRef<double> cdvt, FortranArray1DRef<double> swetw, FortranArray1DRef<double> clw, FortranArray1DRef<double> cdw, FortranArray1DRef<double> cxw, FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw, FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw, FortranArray1DRef<double> cmzw, FortranArray1DRef<double> cdpw, FortranArray1DRef<double> cdvw, int& ncycmax, int& myhost, int& myid, int& mycomm, FortranArray1DRef<int> mblk2nd)
{
    // Access COMMON blocks
    int32_t* ncyc1  = cmn_info.ncyc1;   // 1-based array (0-based in C)
    int32_t& lglobal = cmn_mgrd.lglobal;

    // Local variables
    int ifor, levg, ifo, jfo, kfo;
    int jsta, jend, ksta, kend, ista, iend;
    int nseg, n;

    // individual block summaries
    if (myid == myhost) {
        fortran_write_unit(11, "\n***** FORCE AND MOMENT SUMMARIES - FINAL TIME STEP/MULTIGRID CYCLE *****\n");
    }

    for (n = 1; n <= nblock; n++) {
        ifor = iforce(n);
        levg = levelg(n);
        if (levg == lglobal && ncyc1[iseq-1] > 0 && ifor > 0) {
            if (myid == myhost) {
                ifo = ifor / 100;
                jfo = (ifor - ifo * 100) / 10;
                kfo = (ifor - ifo * 100 - jfo * 10);
                fortran_write_unit(11, "\n SUMMARY OF FORCES AND MOMENTS - BLOCK%6d (GRID%6d)\n\n", n, igridg(n));
                if (ifo == 1) {
                    fortran_write_unit(11, " forces computed on i=1 surface segment(s):\n");
                    for (nseg = 1; nseg <= nbci0(n); nseg++) {
                        jsta = ibcinfo(n, nseg, 2, 1);
                        jend = ibcinfo(n, nseg, 3, 1);
                        ksta = ibcinfo(n, nseg, 4, 1);
                        kend = ibcinfo(n, nseg, 5, 1);
                        if (ibcinfo(n, nseg, 6, 1) > 0)
                            fortran_write_unit(11, "    segment%5d  j=%5d,%5d  k=%5d,%5d\n", nseg, jsta, jend, ksta, kend);
                    }
                }
                if (ifo == 2) {
                    fortran_write_unit(11, " forces computed on i=idim surface segment(s):\n");
                    for (nseg = 1; nseg <= nbcidim(n); nseg++) {
                        jsta = ibcinfo(n, nseg, 2, 2);
                        jend = ibcinfo(n, nseg, 3, 2);
                        ksta = ibcinfo(n, nseg, 4, 2);
                        kend = ibcinfo(n, nseg, 5, 2);
                        if (ibcinfo(n, nseg, 6, 2) > 0)
                            fortran_write_unit(11, "    segment%5d  j=%5d,%5d  k=%5d,%5d\n", nseg, jsta, jend, ksta, kend);
                    }
                }
                if (ifo == 3) {
                    fortran_write_unit(11, " forces computed on i=1 surface segment(s):\n");
                    for (nseg = 1; nseg <= nbci0(n); nseg++) {
                        jsta = ibcinfo(n, nseg, 2, 1);
                        jend = ibcinfo(n, nseg, 3, 1);
                        ksta = ibcinfo(n, nseg, 4, 1);
                        kend = ibcinfo(n, nseg, 5, 1);
                        if (ibcinfo(n, nseg, 6, 1) > 0)
                            fortran_write_unit(11, "    segment%5d  j=%5d,%5d  k=%5d,%5d\n", nseg, jsta, jend, ksta, kend);
                    }
                    fortran_write_unit(11, " forces computed on i=idim surface segment(s):\n");
                    for (nseg = 1; nseg <= nbcidim(n); nseg++) {
                        jsta = ibcinfo(n, nseg, 2, 2);
                        jend = ibcinfo(n, nseg, 3, 2);
                        ksta = ibcinfo(n, nseg, 4, 2);
                        kend = ibcinfo(n, nseg, 5, 2);
                        if (ibcinfo(n, nseg, 6, 2) > 0)
                            fortran_write_unit(11, "    segment%5d  j=%5d,%5d  k=%5d,%5d\n", nseg, jsta, jend, ksta, kend);
                    }
                }
                if (jfo == 1) {
                    fortran_write_unit(11, " forces computed on j=1 surface segment(s):\n");
                    for (nseg = 1; nseg <= nbcj0(n); nseg++) {
                        ista = jbcinfo(n, nseg, 2, 1);
                        iend = jbcinfo(n, nseg, 3, 1);
                        ksta = jbcinfo(n, nseg, 4, 1);
                        kend = jbcinfo(n, nseg, 5, 1);
                        if (jbcinfo(n, nseg, 6, 1) > 0)
                            fortran_write_unit(11, "    segment%5d  i=%5d,%5d  k=%5d,%5d\n", nseg, ista, iend, ksta, kend);
                    }
                }
                if (jfo == 2) {
                    fortran_write_unit(11, " forces computed on j=jdim surface segment(s):\n");
                    for (nseg = 1; nseg <= nbcjdim(n); nseg++) {
                        ista = jbcinfo(n, nseg, 2, 2);
                        iend = jbcinfo(n, nseg, 3, 2);
                        ksta = jbcinfo(n, nseg, 4, 2);
                        kend = jbcinfo(n, nseg, 5, 2);
                        if (jbcinfo(n, nseg, 6, 2) > 0)
                            fortran_write_unit(11, "    segment%5d  i=%5d,%5d  k=%5d,%5d\n", nseg, ista, iend, ksta, kend);
                    }
                }
                if (jfo == 3) {
                    fortran_write_unit(11, " forces computed on j=1 surface segment(s):\n");
                    for (nseg = 1; nseg <= nbcj0(n); nseg++) {
                        ista = jbcinfo(n, nseg, 2, 1);
                        iend = jbcinfo(n, nseg, 3, 1);
                        ksta = jbcinfo(n, nseg, 4, 1);
                        kend = jbcinfo(n, nseg, 5, 1);
                        if (jbcinfo(n, nseg, 6, 1) > 0)
                            fortran_write_unit(11, "    segment%5d  i=%5d,%5d  k=%5d,%5d\n", nseg, ista, iend, ksta, kend);
                    }
                    fortran_write_unit(11, " forces computed on j=jdim surface segment(s):\n");
                    for (nseg = 1; nseg <= nbcjdim(n); nseg++) {
                        ista = jbcinfo(n, nseg, 2, 2);
                        iend = jbcinfo(n, nseg, 3, 2);
                        ksta = jbcinfo(n, nseg, 4, 2);
                        kend = jbcinfo(n, nseg, 5, 2);
                        if (jbcinfo(n, nseg, 6, 2) > 0)
                            fortran_write_unit(11, "    segment%5d  i=%5d,%5d  k=%5d,%5d\n", nseg, ista, iend, ksta, kend);
                    }
                }
                if (kfo == 1) {
                    fortran_write_unit(11, " forces computed on k=1 surface segment(s):\n");
                    for (nseg = 1; nseg <= nbck0(n); nseg++) {
                        ista = kbcinfo(n, nseg, 2, 1);
                        iend = kbcinfo(n, nseg, 3, 1);
                        jsta = kbcinfo(n, nseg, 4, 1);
                        jend = kbcinfo(n, nseg, 5, 1);
                        if (kbcinfo(n, nseg, 6, 1) > 0)
                            fortran_write_unit(11, "    segment%5d  i=%5d,%5d  j=%5d,%5d\n", nseg, ista, iend, jsta, jend);
                    }
                }
                if (kfo == 2) {
                    fortran_write_unit(11, " forces computed on k=kdim surface segment(s):\n");
                    for (nseg = 1; nseg <= nbckdim(n); nseg++) {
                        ista = kbcinfo(n, nseg, 2, 2);
                        iend = kbcinfo(n, nseg, 3, 2);
                        jsta = kbcinfo(n, nseg, 4, 2);
                        jend = kbcinfo(n, nseg, 5, 2);
                        if (kbcinfo(n, nseg, 6, 2) > 0)
                            fortran_write_unit(11, "    segment%5d  i=%5d,%5d  j=%5d,%5d\n", nseg, ista, iend, jsta, jend);
                    }
                }
                if (kfo == 3) {
                    fortran_write_unit(11, " forces computed on k=1 surface segment(s):\n");
                    for (nseg = 1; nseg <= nbck0(n); nseg++) {
                        ista = kbcinfo(n, nseg, 2, 1);
                        iend = kbcinfo(n, nseg, 3, 1);
                        jsta = kbcinfo(n, nseg, 4, 1);
                        jend = kbcinfo(n, nseg, 5, 1);
                        if (kbcinfo(n, nseg, 6, 1) > 0)
                            fortran_write_unit(11, "    segment%5d  i=%5d,%5d  j=%5d,%5d\n", nseg, ista, iend, jsta, jend);
                    }
                    fortran_write_unit(11, " forces computed on k=kdim surface segment(s):\n");
                    for (nseg = 1; nseg <= nbckdim(n); nseg++) {
                        ista = kbcinfo(n, nseg, 2, 2);
                        iend = kbcinfo(n, nseg, 3, 2);
                        jsta = kbcinfo(n, nseg, 4, 2);
                        jend = kbcinfo(n, nseg, 5, 2);
                        if (kbcinfo(n, nseg, 6, 2) > 0)
                            fortran_write_unit(11, "    segment%5d  i=%5d,%5d  j=%5d,%5d\n", nseg, ista, iend, jsta, jend);
                    }
                }
            } // end if (myid == myhost)
            fortran_write_unit(11, "\n         CL-b              CD-b             CDp-b             CDv-b\n");
            fortran_write_unit(11, " %18.11e %18.11e %18.11e %18.11e\n",
                               (float)clt(n), (float)cdt(n), (float)cdpt(n), (float)cdvt(n));
            fortran_write_unit(11, "         CZ-b              CY-b              CX-b          wetted area\n");
            fortran_write_unit(11, " %18.11e %18.11e %18.11e %18.11e\n",
                               (float)czt(n), (float)cyt(n), (float)cxt(n), (float)swett(n));
            fortran_write_unit(11, "         CMY-b             CMX-b             CMZ-b\n");
            fortran_write_unit(11, " %18.11e %18.11e %18.11e\n",
                               (float)cmyt(n), (float)cmxt(n), (float)cmzt(n));
        } // end if (levg == lglobal ...)
    } // end do n=1,nblock

    // summary of global blocks
    int nres = cmn_info.nres;
    fortran_write_unit(11, "\n SUMMARY OF FORCES AND MOMENTS - ALL GLOBAL BLOCKS\n");
    fortran_write_unit(11, "\n          CL                CD               CDp               CDv\n");
    fortran_write_unit(11, " %18.11e %18.11e %18.11e %18.11e\n",
                       (float)clw(nres), (float)cdw(nres), (float)cdpw(nres), (float)cdvw(nres));
    fortran_write_unit(11, "          CZ                CY               CX            wetted area\n");
    fortran_write_unit(11, " %18.11e %18.11e %18.11e %18.11e\n",
                       (float)czw(nres), (float)cyw(nres), (float)cxw(nres), (float)swetw(nres));
    fortran_write_unit(11, "          CMY               CMX               CMZ\n");
    fortran_write_unit(11, " %18.11e %18.11e %18.11e\n",
                       (float)cmyw(nres), (float)cmxw(nres), (float)cmzw(nres));
    fortran_write_unit(11, "\n");

    return;
}

} // namespace forceout_ns
