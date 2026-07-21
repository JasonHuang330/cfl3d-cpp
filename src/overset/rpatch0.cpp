// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "rpatch0.h"
#include "rpatch.h"
#include "setup.h"
#include "runtime/fortran_io.h"
#include "runtime/common_blocks.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

namespace rpatch0_ns {

// ---------------------------------------------------------------------------
// rpatch0_ns::rpatch — delegates to rpatch_ns::rpatch
// ---------------------------------------------------------------------------
void rpatch(int& maxbl, int& maxxe, int& intmax, int& nsub1,
            FortranArray2DRef<double> windex, int& ninter,
            FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt,
            FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy,
            FortranArray2DRef<double> dthetzz, FortranArray1DRef<int> jdimg,
            FortranArray1DRef<int> kdimg, FortranArray1DRef<int> idimg,
            FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
            int& nbuf, int& ibufdim, int& myid)
{
    rpatch_ns::rpatch(maxbl, maxxe, intmax, nsub1, windex, ninter, iindex,
                      nblkpt, dthetxx, dthetyy, dthetzz, jdimg, kdimg, idimg,
                      nou, bou, nbuf, ibufdim, myid);
}

// ---------------------------------------------------------------------------
// rpatch0 — Read in generalized-coordinate interpolation data from a file;
//           only the data needed to evaluate sizing requirements are stored.
// ---------------------------------------------------------------------------
void rpatch0(int& intmx, int& msub1, FortranArray2DRef<int> iindex, int& ninter)
{
    // local variables
    // Fortran REAL under -fdefault-real-8 is 8 bytes; reading it as float (4)
    // under-consumes the record and misaligns every subsequent read → garbage
    // len for interfaces past the first → lmaxxe<0.
    double dum;
    int nfb, idum1, idum2, idum3, nbl, lst, len;
    int ndum, nnn;

    // Aliases for COMMON /params/
    int& lintmax  = cmn_params.lintmax;
    int& lnsub1   = cmn_params.lnsub1;
    int& lmaxxe   = cmn_params.lmaxxe;

    FILE* fp = fortran_get_unit(22);

    // Helper lambdas for reading Fortran unformatted sequential records
    // Each record: [4-byte length][data][4-byte length]
    auto read_rec_begin = [&](int32_t& marker) {
        fread(&marker, sizeof(int32_t), 1, fp);
    };
    auto read_rec_end = [&](int32_t& marker) {
        fread(&marker, sizeof(int32_t), 1, fp);
    };

    // read(22) ninter
    {
        int32_t m1, m2;
        read_rec_begin(m1);
        fread(&ninter, sizeof(int), 1, fp);
        read_rec_end(m2);
    }

    lintmax = ninter;
    lnsub1  = 0;

    len = 0;  // initialise so lmaxxe = len is valid even if abs(ninter)==0

    for (int n = 1; n <= std::abs(ninter); n++) {

        // read(22) iindex(n,1)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, 1), sizeof(int), 1, fp);
            read_rec_end(m2);
        }
        nfb = iindex(n, 1);
        lnsub1 = std::max(lnsub1, nfb);

        // read(22) (iindex(n,1+ll),ll=1,nfb)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            for (int ll = 1; ll <= nfb; ll++) {
                fread(&iindex(n, 1 + ll), sizeof(int), 1, fp);
            }
            read_rec_end(m2);
        }

        // read(22) iindex(n,nfb+2)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, nfb + 2), sizeof(int), 1, fp);
            read_rec_end(m2);
        }
        idum1 = iindex(n, nfb + 2);

        // read(22) (iindex(n,nfb+2+ll),ll=1,nfb)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            for (int ll = 1; ll <= nfb; ll++) {
                fread(&iindex(n, nfb + 2 + ll), sizeof(int), 1, fp);
            }
            read_rec_end(m2);
        }

        // read(22) iindex(n,2*nfb+3)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, 2 * nfb + 3), sizeof(int), 1, fp);
            read_rec_end(m2);
        }

        // read(22) iindex(n,2*nfb+4)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, 2 * nfb + 4), sizeof(int), 1, fp);
            read_rec_end(m2);
        }
        idum2 = iindex(n, 2 * nfb + 4);

        // read(22) iindex(n,2*nfb+5)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, 2 * nfb + 5), sizeof(int), 1, fp);
            read_rec_end(m2);
        }
        idum3 = iindex(n, 2 * nfb + 5);

        // read(22) iindex(n,2*nfb+6)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, 2 * nfb + 6), sizeof(int), 1, fp);
            read_rec_end(m2);
        }

        // read(22) iindex(n,2*nfb+7)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, 2 * nfb + 7), sizeof(int), 1, fp);
            read_rec_end(m2);
        }

        // read(22) iindex(n,2*nfb+8)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, 2 * nfb + 8), sizeof(int), 1, fp);
            read_rec_end(m2);
        }

        // read(22) iindex(n,2*nfb+9)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            fread(&iindex(n, 2 * nfb + 9), sizeof(int), 1, fp);
            read_rec_end(m2);
        }

        nbl = idum1;
        lst = idum3;
        len = lst + idum2 - 1;

        // read(22) (ndum,nnn=lst,len)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            for (nnn = lst; nnn <= len; nnn++) {
                fread(&ndum, sizeof(int), 1, fp);
            }
            read_rec_end(m2);
        }

        // read(22) ((dum,nnn=lst,len),ll=1,2)
        {
            int32_t m1, m2;
            read_rec_begin(m1);
            for (int ll = 1; ll <= 2; ll++) {
                for (nnn = lst; nnn <= len; nnn++) {
                    fread(&dum, sizeof(double), 1, fp);
                }
            }
            read_rec_end(m2);
        }

    } // end do 1500

    // don't need to read angular displacements for sizing purposes
    lmaxxe = len;

    // set ninter = -ninter to allow full read of
    // the patch file in subroutine setup
    ninter = -ninter;
}

// ---------------------------------------------------------------------------
// rpatch0_ns::setup — delegates to setup_ns::setup
// ---------------------------------------------------------------------------
void setup(FortranArray2DRef<int> lw, FortranArray2DRef<int> lw2,
           FortranArray1DRef<double> w, int& mgwk, FortranArray1DRef<double> wk,
           int& nwork, FortranArray1DRef<int> iwk, int& iwork,
           FortranArray1DRef<int> itest, FortranArray1DRef<int> jtest,
           FortranArray1DRef<int> ktest, int& maxbl, int& mxbli, int& maxgr,
           int& maxseg, int& nsub1, int& maxxe, int& intmax, int& iitot,
           int& ncycmax, FortranArray3DRef<int> lwdat, FortranArray1DRef<int> lig,
           FortranArray1DRef<int> lbg, FortranArray1DRef<int> iovrlp,
           FortranArray3DRef<double> qb, int& nblock, FortranArray2DRef<int> iviscg,
           FortranArray1DRef<int> jdimg, FortranArray1DRef<int> kdimg,
           FortranArray1DRef<int> idimg, FortranArray1DRef<double> utrans,
           FortranArray1DRef<double> vtrans, FortranArray1DRef<double> wtrans,
           FortranArray1DRef<double> omegax, FortranArray1DRef<double> omegay,
           FortranArray1DRef<double> omegaz, FortranArray1DRef<double> xorig,
           FortranArray1DRef<double> yorig, FortranArray1DRef<double> zorig,
           FortranArray1DRef<double> dxmx, FortranArray1DRef<double> dymx,
           FortranArray1DRef<double> dzmx, FortranArray1DRef<double> dthxmx,
           FortranArray1DRef<double> dthymx, FortranArray1DRef<double> dthzmx,
           FortranArray1DRef<double> thetax, FortranArray1DRef<double> thetay,
           FortranArray1DRef<double> thetaz, FortranArray1DRef<double> rfreqt,
           FortranArray1DRef<double> rfreqr, FortranArray1DRef<double> xorig0,
           FortranArray1DRef<double> yorig0, FortranArray1DRef<double> zorig0,
           FortranArray1DRef<double> time2, FortranArray1DRef<double> thetaxl,
           FortranArray1DRef<double> thetayl, FortranArray1DRef<double> thetazl,
           FortranArray1DRef<int> itrans, FortranArray1DRef<int> irotat,
           FortranArray1DRef<int> idefrm, FortranArray4DRef<double> bcvali,
           FortranArray4DRef<double> bcvalj, FortranArray4DRef<double> bcvalk,
           FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcidim,
           FortranArray1DRef<int> nbcj0, FortranArray1DRef<int> nbcjdim,
           FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbckdim,
           FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
           FortranArray4DRef<int> kbcinfo, FortranArray3DRef<int> bcfilei,
           FortranArray3DRef<int> bcfilej, FortranArray3DRef<int> bcfilek,
           int& ngrid, FortranArray1DRef<int> ncgg, FortranArray1DRef<int> nblg,
           FortranArray1DRef<int> iemg, FortranArray1DRef<int> inewgg,
           FortranArray1DRef<double> rms, FortranArray1DRef<double> clw,
           FortranArray1DRef<double> cdw, FortranArray1DRef<double> cdpw,
           FortranArray1DRef<double> cdvw, FortranArray1DRef<double> cxw,
           FortranArray1DRef<double> cyw, FortranArray1DRef<double> czw,
           FortranArray1DRef<double> cmxw, FortranArray1DRef<double> cmyw,
           FortranArray1DRef<double> cmzw, int& n_clcd,
           FortranArray3DRef<double> clcd, int& nblocks_clcd,
           FortranArray2DRef<int> blocks_clcd, FortranArray1DRef<double> chdw,
           FortranArray1DRef<double> swetw, FortranArray1DRef<double> fmdotw,
           FortranArray1DRef<double> cfttotw, FortranArray1DRef<double> cftmomw,
           FortranArray1DRef<double> cftpw, FortranArray1DRef<double> cftvw,
           FortranArray2DRef<double> rmstr, FortranArray2DRef<int> nneg,
           int& ntr, FortranArray2DRef<double> windex, int& ninter,
           FortranArray2DRef<int> iindex, FortranArray1DRef<int> nblkpt,
           FortranArray2DRef<double> dthetxx, FortranArray2DRef<double> dthetyy,
           FortranArray2DRef<double> dthetzz, FortranArray1DRef<int> iibg,
           FortranArray1DRef<int> kkbg, FortranArray1DRef<int> jjbg,
           FortranArray1DRef<int> ibcg, FortranArray1DRef<double> dxintg,
           FortranArray1DRef<double> dyintg, FortranArray1DRef<double> dzintg,
           FortranArray1DRef<int> iiig, FortranArray1DRef<int> jjig,
           FortranArray1DRef<int> kkig, FortranArray2DRef<int> ibpntsg,
           FortranArray1DRef<int> iipntsg, FortranArray1DRef<int> mblk2nd,
           FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
           int& nbuf, int& ibufdim, FortranArray1DRef<int> ireq_qb,
           FortranArray1DRef<int> igridg, FortranArray1DRef<char[80]> bcfiles,
           int& mxbcfil, FortranArray2DRef<double> utrnsae,
           FortranArray2DRef<double> vtrnsae, FortranArray2DRef<double> wtrnsae,
           FortranArray2DRef<double> omgxae, FortranArray2DRef<double> omgyae,
           FortranArray2DRef<double> omgzae, FortranArray2DRef<double> xorgae,
           FortranArray2DRef<double> yorgae, FortranArray2DRef<double> zorgae,
           FortranArray2DRef<double> thtxae, FortranArray2DRef<double> thtyae,
           FortranArray2DRef<double> thtzae, FortranArray2DRef<double> rfrqtae,
           FortranArray2DRef<double> rfrqrae, FortranArray2DRef<int> icsi,
           FortranArray2DRef<int> icsf, FortranArray2DRef<int> jcsi,
           FortranArray2DRef<int> jcsf, FortranArray2DRef<int> kcsi,
           FortranArray2DRef<int> kcsf, FortranArray2DRef<double> freq,
           FortranArray2DRef<double> gmass, FortranArray2DRef<double> damp,
           FortranArray2DRef<double> x0, FortranArray2DRef<double> gf0,
           int& nmds, int& maxaes, FortranArray2DRef<double> aesrfdat,
           FortranArray3DRef<double> perturb, FortranArray3DRef<int> islavept,
           int& nslave, FortranArray2DRef<int> iskip, FortranArray2DRef<int> jskip,
           FortranArray2DRef<int> kskip, FortranArray3DRef<double> bmat,
           FortranArray3DRef<double> stm, FortranArray3DRef<double> stmi,
           FortranArray2DRef<double> gforcn, FortranArray2DRef<double> gforcnm,
           FortranArray2DRef<double> xxn, FortranArray1DRef<int> nsegdfrm,
           FortranArray2DRef<int> idfrmseg, FortranArray2DRef<int> iaesurf,
           int& maxsegdg, int& nmaster, FortranArray4DRef<double> aehist,
           FortranArray1DRef<double> timekeep, FortranArray2DRef<int> inpl3d,
           int& nplots, int& nplot3d, FortranArray1DRef<int> levelg,
           FortranArray1DRef<int> iadvance, FortranArray2DRef<double> xs,
           FortranArray2DRef<double> gforcs, FortranArray2DRef<double> xorgae0,
           FortranArray2DRef<double> yorgae0, FortranArray2DRef<double> zorgae0,
           FortranArray2DRef<int> icouple, int& lfgm, FortranArray2DRef<int> nblk,
           FortranArray3DRef<int> limblk, FortranArray3DRef<int> isva,
           FortranArray2DRef<int> nblelst, FortranArray1DRef<int> iskmax,
           FortranArray1DRef<int> jskmax, FortranArray1DRef<int> kskmax,
           FortranArray1DRef<double> ue, FortranArray1DRef<int> irdrea,
           int& nbli, int& nummem)
{
    setup_ns::setup(lw, lw2, w, mgwk, wk, nwork, iwk, iwork, itest, jtest,
                    ktest, maxbl, mxbli, maxgr, maxseg, nsub1, maxxe, intmax,
                    iitot, ncycmax, lwdat, lig, lbg, iovrlp, qb, nblock,
                    iviscg, jdimg, kdimg, idimg, utrans, vtrans, wtrans,
                    omegax, omegay, omegaz, xorig, yorig, zorig, dxmx, dymx,
                    dzmx, dthxmx, dthymx, dthzmx, thetax, thetay, thetaz,
                    rfreqt, rfreqr, xorig0, yorig0, zorig0, time2, thetaxl,
                    thetayl, thetazl, itrans, irotat, idefrm, bcvali, bcvalj,
                    bcvalk, nbci0, nbcidim, nbcj0, nbcjdim, nbck0, nbckdim,
                    ibcinfo, jbcinfo, kbcinfo, bcfilei, bcfilej, bcfilek,
                    ngrid, ncgg, nblg, iemg, inewgg, rms, clw, cdw, cdpw,
                    cdvw, cxw, cyw, czw, cmxw, cmyw, cmzw, n_clcd, clcd,
                    nblocks_clcd, blocks_clcd, chdw, swetw, fmdotw, cfttotw,
                    cftmomw, cftpw, cftvw, rmstr, nneg, ntr, windex, ninter,
                    iindex, nblkpt, dthetxx, dthetyy, dthetzz, iibg, kkbg,
                    jjbg, ibcg, dxintg, dyintg, dzintg, iiig, jjig, kkig,
                    ibpntsg, iipntsg, mblk2nd, nou, bou, nbuf, ibufdim,
                    ireq_qb, igridg, bcfiles, mxbcfil, utrnsae, vtrnsae,
                    wtrnsae, omgxae, omgyae, omgzae, xorgae, yorgae, zorgae,
                    thtxae, thtyae, thtzae, rfrqtae, rfrqrae, icsi, icsf,
                    jcsi, jcsf, kcsi, kcsf, freq, gmass, damp, x0, gf0,
                    nmds, maxaes, aesrfdat, perturb, islavept, nslave, iskip,
                    jskip, kskip, bmat, stm, stmi, gforcn, gforcnm, xxn,
                    nsegdfrm, idfrmseg, iaesurf, maxsegdg, nmaster, aehist,
                    timekeep, inpl3d, nplots, nplot3d, levelg, iadvance, xs,
                    gforcs, xorgae0, yorgae0, zorgae0, icouple, lfgm, nblk,
                    limblk, isva, nblelst, iskmax, jskmax, kskmax, ue, irdrea,
                    nbli, nummem);
}

} // namespace rpatch0_ns
