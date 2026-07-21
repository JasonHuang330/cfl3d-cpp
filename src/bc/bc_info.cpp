// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// This file is part of cfl3d-cpp, a C++ translation (MODIFIED) of NASA CFL3D
// (https://github.com/nasa/CFL3D). Independent, unofficial; not endorsed by
// or affiliated with NASA.
#include "bc_info.h"
#include "termn8.h"
#include "runtime/fortran_io.h"
#include <cstring>
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace bc_info_ns {

// ============================================================
// bc() - dispatcher that calls bc_info()
// ============================================================
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
    // mblk2nd is stored in wk starting at nummem (integer data in double array)
    int* mblk2nd_ptr = reinterpret_cast<int*>(&wk(nummem));
    FortranArray1D<int> mblk2nd_arr(maxbl);
    for (int i = 1; i <= maxbl; i++) {
        mblk2nd_arr(i) = mblk2nd_ptr[i-1];
    }
    FortranArray1DRef<int> mblk2nd_ref = mblk2nd_arr.ref();

    // Extract scalar parameters from lw/lw2
    int iitot  = lw(1,nbl);
    int intmax = lw(2,nbl);
    int nsub1  = lw(3,nbl);
    int mxbli  = lw(4,nbl);
    int lbcprd = lw(5,nbl);
    int nbli   = lw2(1,nbl);
    int ninter = lw2(2,nbl);

    double epsc0  = 1.e-3;
    double epsrot = 1.e-3;

    // Build array refs from wk for the chimera/patch/period arrays
    // These are stored as integer data in wk (Fortran equivalence trick)
    // Layout (Fortran column-major, integers packed 2 per double):
    // iibg(iitot), kkbg(iitot), jjbg(iitot), ibcg(iitot)
    // ibpntsg(maxbl,4), iipntsg(maxbl), lig(maxbl), lbg(maxbl)
    // isav_blk(2*mxbli,17), isav_prd(lbcprd,12)
    // iemg(maxgr)
    // geom_miss(2*mxbli), period_miss(lbcprd)  [doubles]
    // iindex(intmax,6*nsub1+9)
    // iovrlp(maxbl)

    int* iwk = reinterpret_cast<int*>(&wk(1));

    int off = 0;
    int liibg    = off; off += iitot;
    int lkkbg    = off; off += iitot;
    int ljjbg    = off; off += iitot;
    int libcg    = off; off += iitot;
    int libpntsg = off; off += maxbl*4;
    int liipntsg = off; off += maxbl;
    int llig     = off; off += maxbl;
    int llbg     = off; off += maxbl;
    int lisav_blk= off; off += 2*mxbli*17;
    int lisav_prd= off; off += lbcprd*12;
    int liemg    = off; off += maxgr;
    // geom_miss and period_miss are double arrays
    // off integers = off/2 doubles (rounded up to next double boundary)
    int off_dbl = (off + 1) / 2;
    int lgeom_dbl   = off_dbl;
    int lperiod_dbl = lgeom_dbl + 2*mxbli;
    // iindex is integer, stored after period_miss
    int liindex = (lperiod_dbl + lbcprd) * 2;
    int liovrlp = liindex + intmax*(6*nsub1+9);

    FortranArray1D<int> iibg_arr(iitot), kkbg_arr(iitot), jjbg_arr(iitot), ibcg_arr(iitot);
    for (int i=0; i<iitot; i++) {
        iibg_arr(i+1) = iwk[liibg+i];
        kkbg_arr(i+1) = iwk[lkkbg+i];
        jjbg_arr(i+1) = iwk[ljjbg+i];
        ibcg_arr(i+1) = iwk[libcg+i];
    }

    FortranArray2D<int> ibpntsg_arr(maxbl, 4);
    for (int j=0; j<4; j++)
        for (int i=0; i<maxbl; i++)
            ibpntsg_arr(i+1,j+1) = iwk[libpntsg+j*maxbl+i];

    FortranArray1D<int> iipntsg_arr(maxbl), lig_arr(maxbl), lbg_arr(maxbl), iovrlp_arr(maxbl);
    for (int i=0; i<maxbl; i++) {
        iipntsg_arr(i+1) = iwk[liipntsg+i];
        lig_arr(i+1)     = iwk[llig+i];
        lbg_arr(i+1)     = iwk[llbg+i];
        iovrlp_arr(i+1)  = iwk[liovrlp+i];
    }

    FortranArray2D<int> isav_blk_arr(2*mxbli, 17);
    for (int j=0; j<17; j++)
        for (int i=0; i<2*mxbli; i++)
            isav_blk_arr(i+1,j+1) = iwk[lisav_blk+j*2*mxbli+i];

    FortranArray2D<int> isav_prd_arr(lbcprd, 12);
    for (int j=0; j<12; j++)
        for (int i=0; i<lbcprd; i++)
            isav_prd_arr(i+1,j+1) = iwk[lisav_prd+j*lbcprd+i];

    FortranArray1D<int> iemg_arr(maxgr);
    for (int i=0; i<maxgr; i++)
        iemg_arr(i+1) = iwk[liemg+i];

    FortranArray1D<double> geom_miss_arr(2*mxbli), period_miss_arr(lbcprd);
    for (int i=0; i<2*mxbli; i++)
        geom_miss_arr(i+1) = wk(lgeom_dbl+1+i);
    for (int i=0; i<lbcprd; i++)
        period_miss_arr(i+1) = wk(lperiod_dbl+1+i);

    FortranArray2D<int> iindex_arr(intmax, 6*nsub1+9);
    for (int j=0; j<6*nsub1+9; j++)
        for (int i=0; i<intmax; i++)
            iindex_arr(i+1,j+1) = iwk[liindex+j*intmax+i];

    FortranArray1DRef<int> iibg_ref    = iibg_arr.ref();
    FortranArray1DRef<int> kkbg_ref    = kkbg_arr.ref();
    FortranArray1DRef<int> jjbg_ref    = jjbg_arr.ref();
    FortranArray1DRef<int> ibcg_ref    = ibcg_arr.ref();
    FortranArray2DRef<int> ibpntsg_ref = ibpntsg_arr.ref();
    FortranArray1DRef<int> iipntsg_ref = iipntsg_arr.ref();
    FortranArray1DRef<int> lig_ref     = lig_arr.ref();
    FortranArray1DRef<int> lbg_ref     = lbg_arr.ref();
    FortranArray2DRef<int> isav_blk_ref= isav_blk_arr.ref();
    FortranArray2DRef<int> isav_prd_ref= isav_prd_arr.ref();
    FortranArray1DRef<int> iemg_ref    = iemg_arr.ref();
    FortranArray1DRef<double> geom_miss_ref   = geom_miss_arr.ref();
    FortranArray1DRef<double> period_miss_ref = period_miss_arr.ref();
    FortranArray2DRef<int> iindex_ref  = iindex_arr.ref();
    FortranArray1DRef<int> iovrlp_ref  = iovrlp_arr.ref();

    bc_info(ntime, nbl, w, mgwk, ibcinfo, jbcinfo, kbcinfo,
            nbci0, nbcj0, nbck0, nbcidim, nbcjdim, nbckdim,
            bcfilei, bcfilej, bcfilek, igridg, itrans, irotat, idefrm,
            idimg, jdimg, kdimg, nblg,
            iibg_ref, kkbg_ref, jjbg_ref, ibcg_ref,
            ibpntsg_ref, iipntsg_ref, lig_ref, lbg_ref,
            isav_blk_ref, isav_prd_ref, iemg_ref,
            nbli, geom_miss_ref, period_miss_ref,
            ninter, iindex_ref,
            nou, bou, nbuf, myid, mblk2nd_ref, ibufdim, maxbl, maxseg,
            iitot, intmax, nsub1, mxbli, lbcprd, epsc0, epsrot,
            lwdat, maxgr, iovrlp_ref, bcfiles, mxbcfil);
}


// ============================================================
// bc_info() - print boundary condition information
// ============================================================
void bc_info(int& ntime, int& nbl, FortranArray1DRef<double> w, int& mgwk,
             FortranArray4DRef<int> ibcinfo, FortranArray4DRef<int> jbcinfo,
             FortranArray4DRef<int> kbcinfo,
             FortranArray1DRef<int> nbci0, FortranArray1DRef<int> nbcj0,
             FortranArray1DRef<int> nbck0, FortranArray1DRef<int> nbcidim,
             FortranArray1DRef<int> nbcjdim, FortranArray1DRef<int> nbckdim,
             FortranArray3DRef<int> bcfilei, FortranArray3DRef<int> bcfilej,
             FortranArray3DRef<int> bcfilek,
             FortranArray1DRef<int> igridg, FortranArray1DRef<int> itrans,
             FortranArray1DRef<int> irotat, FortranArray1DRef<int> idefrm,
             FortranArray1DRef<int> idimg, FortranArray1DRef<int> jdimg,
             FortranArray1DRef<int> kdimg, FortranArray1DRef<int> nblg,
             FortranArray1DRef<int> iibg, FortranArray1DRef<int> kkbg,
             FortranArray1DRef<int> jjbg, FortranArray1DRef<int> ibcg,
             FortranArray2DRef<int> ibpntsg, FortranArray1DRef<int> iipntsg,
             FortranArray1DRef<int> lig, FortranArray1DRef<int> lbg,
             FortranArray2DRef<int> isav_blk, FortranArray2DRef<int> isav_prd,
             FortranArray1DRef<int> iemg,
             int& nbli, FortranArray1DRef<double> geom_miss,
             FortranArray1DRef<double> period_miss,
             int& ninter, FortranArray2DRef<int> iindex,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& myid, FortranArray1DRef<int> mblk2nd,
             int& ibufdim, int& maxbl, int& maxseg,
             int& iitot, int& intmax, int& nsub1, int& mxbli, int& lbcprd,
             double& epsc0, double& epsrot,
             FortranArray3DRef<int> lwdat, int& maxgr,
             FortranArray1DRef<int> iovrlp,
             FortranArray1DRef<char[80]> bcfiles, int& mxbcfil)
{
    char filname[80];
    int& jdim = cmn_ginfo.jdim;
    int& kdim = cmn_ginfo.kdim;
    int& idim = cmn_ginfo.idim;
    int& isklton = cmn_sklton.isklton;
    int& level = cmn_mgrd.level;
    int* is_blk = cmn_is_blockbc.is_blk;
    int* ie_blk = cmn_is_blockbc.ie_blk;
    int* is_prd = cmn_is_perbc.is_prd;
    int* ie_prd = cmn_is_perbc.ie_prd;

    int ista, iend, jsta, jend, ksta, kend;
    int ldata, mdim, ndim;
    int nseg, ngh, nblnum, ngnew, jdimp, kdimp, idimp, nblp;
    int nface, lcntp, nbl_chk, nface_chk, nseg_chk;
    int iuns;

    if (isklton != 1) return;

    iuns = std::max({irotat(nbl), itrans(nbl), idefrm(nbl)});

    if (iuns == 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " boundary conditions for block%6d (grid %6d)", nbl, igridg(nbl));
    }
    if (iuns > 0) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " ");
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " boundary conditions for block%6d (grid %6d) - dynamic mesh", nbl, igridg(nbl));
    }

    // i=0 boundary (face 1)
    for (nseg = 1; nseg <= nbci0(nbl); nseg++) {
        ista = 1; iend = 1;
        jsta = ibcinfo(nbl,nseg,2,1); jend = ibcinfo(nbl,nseg,3,1);
        ksta = ibcinfo(nbl,nseg,4,1); kend = ibcinfo(nbl,nseg,5,1);
        ldata = lwdat(nbl,nseg,1);
        mdim = jend-jsta; ndim = kend-ksta;
        std::memcpy(filname, bcfiles(bcfilei(nbl,nseg,1)), 80);
        nface = 1;

        if (ibcinfo(nbl,nseg,1,1) == 9999)
            out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,1) == 1000)
            out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,1) == 1001)
            out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,1) == 1002)
            out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,1) == 1003)
            out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,myid);
        if (ibcinfo(nbl,nseg,1,1) == 1005) { int inm=0;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (ibcinfo(nbl,nseg,1,1) == 1006) { int inm=1;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (ibcinfo(nbl,nseg,1,1) == 1008)
            out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,1) == 1011)
            out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,1) == 1012)
            out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,1) == 1013)
            out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);

        if (ibcinfo(nbl,nseg,1,1) == 2002) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2003) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2009) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=0;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2010) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=1;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2019) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2004) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2014) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2024) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(ibcinfo(nbl,nseg,1,1)) == 2034) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2005) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            jdimp = jdimg(ngnew); kdimp = kdimg(ngnew); idimp = idimg(ngnew);
            nblp = ngnew;
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,nbl,jdimp,kdimp,idimp,nblp,filname,myid,mblk2nd,maxbl);
            for (lcntp = is_prd[level-1]; lcntp <= ie_prd[level-1]; lcntp++) {
                nbl_chk   = isav_prd(lcntp,1);
                nface_chk = isav_prd(lcntp,2);
                nseg_chk  = isav_prd(lcntp,11);
                if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                    if ((float)period_miss(lcntp) > (float)epsrot) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "           blocking check....periodic mismatch = %14.7E", (float)period_miss(lcntp));
                    }
                }
            }
        }
        if (ibcinfo(nbl,nseg,1,1) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = 0;
            if (nblnum > 0) ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,ngnew,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2007) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2008) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2016) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=5;
            out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2026) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2018) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2028) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2038) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2102) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=0;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        if (ibcinfo(nbl,nseg,1,1) == 2103) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=1;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        // block interface / chimera / embedding checks for face 1
        for (int lcnt = is_blk[level-1]; lcnt <= ie_blk[level-1]; lcnt++) {
            nbl_chk   = isav_blk(lcnt,1);
            nface_chk = isav_blk(lcnt,2);
            nseg_chk  = isav_blk(lcnt,17);
            if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                if ((float)geom_miss(lcnt) > (float)epsc0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "           blocking check....geometry mismatch = %14.7E", (float)geom_miss(lcnt));
                }
            }
        }
    } // end do i=0 face


    // i=idim boundary (face 2)
    for (nseg = 1; nseg <= nbcidim(nbl); nseg++) {
        ista = idim; iend = idim;
        jsta = ibcinfo(nbl,nseg,2,2); jend = ibcinfo(nbl,nseg,3,2);
        ksta = ibcinfo(nbl,nseg,4,2); kend = ibcinfo(nbl,nseg,5,2);
        ldata = lwdat(nbl,nseg,2);
        mdim = jend-jsta; ndim = kend-ksta;
        std::memcpy(filname, bcfiles(bcfilei(nbl,nseg,2)), 80);
        nface = 2;

        if (ibcinfo(nbl,nseg,1,2) == 9999)
            out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,2) == 1000)
            out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,2) == 1001)
            out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,2) == 1002)
            out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,2) == 1003)
            out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,myid);
        if (ibcinfo(nbl,nseg,1,2) == 1005) { int inm=0;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (ibcinfo(nbl,nseg,1,2) == 1006) { int inm=1;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (ibcinfo(nbl,nseg,1,2) == 1008)
            out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,2) == 1011)
            out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,2) == 1012)
            out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,2) == 1013)
            out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (ibcinfo(nbl,nseg,1,2) == 2002) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2003) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2009) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=0;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2010) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=1;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2019) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2004) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2014) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2024) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(ibcinfo(nbl,nseg,1,2)) == 2034) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2005) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            jdimp = jdimg(ngnew); kdimp = kdimg(ngnew); idimp = idimg(ngnew);
            nblp = ngnew;
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,nbl,jdimp,kdimp,idimp,nblp,filname,myid,mblk2nd,maxbl);
            for (lcntp = is_prd[level-1]; lcntp <= ie_prd[level-1]; lcntp++) {
                nbl_chk   = isav_prd(lcntp,1);
                nface_chk = isav_prd(lcntp,2);
                nseg_chk  = isav_prd(lcntp,11);
                if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                    if ((float)period_miss(lcntp) > (float)epsrot) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "           blocking check....periodic mismatch = %14.7E", (float)period_miss(lcntp));
                    }
                }
            }
        }
        if (ibcinfo(nbl,nseg,1,2) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = 0;
            if (nblnum > 0) ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,ngnew,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2007) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2008) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2016) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=5;
            out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2026) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2018) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2028) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2038) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2102) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=0;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        if (ibcinfo(nbl,nseg,1,2) == 2103) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=1;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        for (int lcnt = is_blk[level-1]; lcnt <= ie_blk[level-1]; lcnt++) {
            nbl_chk   = isav_blk(lcnt,1);
            nface_chk = isav_blk(lcnt,2);
            nseg_chk  = isav_blk(lcnt,17);
            if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                if ((float)geom_miss(lcnt) > (float)epsc0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "           blocking check....geometry mismatch = %14.7E", (float)geom_miss(lcnt));
                }
            }
        }
    } // end do i=idim face


    // j=0 boundary (face 3)
    for (nseg = 1; nseg <= nbcj0(nbl); nseg++) {
        jsta = 1; jend = 1;
        ista = jbcinfo(nbl,nseg,2,1); iend = jbcinfo(nbl,nseg,3,1);
        ksta = jbcinfo(nbl,nseg,4,1); kend = jbcinfo(nbl,nseg,5,1);
        ldata = lwdat(nbl,nseg,3);
        mdim = iend-ista; ndim = kend-ksta;
        std::memcpy(filname, bcfiles(bcfilej(nbl,nseg,1)), 80);
        nface = 3;

        if (jbcinfo(nbl,nseg,1,1) == 9999)
            out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,1) == 1000)
            out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,1) == 1001)
            out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,1) == 1002)
            out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,1) == 1003)
            out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,myid);
        if (jbcinfo(nbl,nseg,1,1) == 1005) { int inm=0;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (jbcinfo(nbl,nseg,1,1) == 1006) { int inm=1;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (jbcinfo(nbl,nseg,1,1) == 1008)
            out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,1) == 1011)
            out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,1) == 1012)
            out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,1) == 1013)
            out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,1) == 2002) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2003) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2009) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=0;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2010) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=1;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2019) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2004) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2014) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2024) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(jbcinfo(nbl,nseg,1,1)) == 2034) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2005) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            jdimp = jdimg(ngnew); kdimp = kdimg(ngnew); idimp = idimg(ngnew);
            nblp = ngnew;
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,nbl,jdimp,kdimp,idimp,nblp,filname,myid,mblk2nd,maxbl);
            for (lcntp = is_prd[level-1]; lcntp <= ie_prd[level-1]; lcntp++) {
                nbl_chk   = isav_prd(lcntp,1);
                nface_chk = isav_prd(lcntp,2);
                nseg_chk  = isav_prd(lcntp,11);
                if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                    if ((float)period_miss(lcntp) > (float)epsrot) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "           blocking check....periodic mismatch = %14.7E", (float)period_miss(lcntp));
                    }
                }
            }
        }
        if (jbcinfo(nbl,nseg,1,1) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = 0;
            if (nblnum > 0) ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,ngnew,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2007) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2008) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2016) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=5;
            out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2026) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2018) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2028) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2038) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2102) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=0;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        if (jbcinfo(nbl,nseg,1,1) == 2103) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=1;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        for (int lcnt = is_blk[level-1]; lcnt <= ie_blk[level-1]; lcnt++) {
            nbl_chk   = isav_blk(lcnt,1);
            nface_chk = isav_blk(lcnt,2);
            nseg_chk  = isav_blk(lcnt,17);
            if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                if ((float)geom_miss(lcnt) > (float)epsc0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "           blocking check....geometry mismatch = %14.7E", (float)geom_miss(lcnt));
                }
            }
        }
    } // end do j=0 face


    // j=jdim boundary (face 4)
    for (nseg = 1; nseg <= nbcjdim(nbl); nseg++) {
        jsta = jdim; jend = jdim;
        ista = jbcinfo(nbl,nseg,2,2); iend = jbcinfo(nbl,nseg,3,2);
        ksta = jbcinfo(nbl,nseg,4,2); kend = jbcinfo(nbl,nseg,5,2);
        ldata = lwdat(nbl,nseg,4);
        mdim = iend-ista; ndim = kend-ksta;
        std::memcpy(filname, bcfiles(bcfilej(nbl,nseg,2)), 80);
        nface = 4;

        if (jbcinfo(nbl,nseg,1,2) == 9999)
            out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,2) == 1000)
            out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,2) == 1001)
            out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,2) == 1002)
            out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,2) == 1003)
            out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,myid);
        if (jbcinfo(nbl,nseg,1,2) == 1005) { int inm=0;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (jbcinfo(nbl,nseg,1,2) == 1006) { int inm=1;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (jbcinfo(nbl,nseg,1,2) == 1008)
            out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,2) == 1011)
            out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,2) == 1012)
            out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,2) == 1013)
            out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (jbcinfo(nbl,nseg,1,2) == 2002) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2003) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2009) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=0;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2010) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=1;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2019) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2004) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2014) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2024) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(jbcinfo(nbl,nseg,1,2)) == 2034) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2005) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            jdimp = jdimg(ngnew); kdimp = kdimg(ngnew); idimp = idimg(ngnew);
            nblp = ngnew;
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,nbl,jdimp,kdimp,idimp,nblp,filname,myid,mblk2nd,maxbl);
            for (lcntp = is_prd[level-1]; lcntp <= ie_prd[level-1]; lcntp++) {
                nbl_chk   = isav_prd(lcntp,1);
                nface_chk = isav_prd(lcntp,2);
                nseg_chk  = isav_prd(lcntp,11);
                if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                    if ((float)period_miss(lcntp) > (float)epsrot) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "           blocking check....periodic mismatch = %14.7E", (float)period_miss(lcntp));
                    }
                }
            }
        }
        if (jbcinfo(nbl,nseg,1,2) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = 0;
            if (nblnum > 0) ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,ngnew,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2007) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2008) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2016) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=5;
            out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2026) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2018) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2028) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2038) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2102) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=0;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        if (jbcinfo(nbl,nseg,1,2) == 2103) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=1;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        for (int lcnt = is_blk[level-1]; lcnt <= ie_blk[level-1]; lcnt++) {
            nbl_chk   = isav_blk(lcnt,1);
            nface_chk = isav_blk(lcnt,2);
            nseg_chk  = isav_blk(lcnt,17);
            if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                if ((float)geom_miss(lcnt) > (float)epsc0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "           blocking check....geometry mismatch = %14.7E", (float)geom_miss(lcnt));
                }
            }
        }
    } // end do j=jdim face


    // k=0 boundary (face 5)
    for (nseg = 1; nseg <= nbck0(nbl); nseg++) {
        ksta = 1; kend = 1;
        ista = kbcinfo(nbl,nseg,2,1); iend = kbcinfo(nbl,nseg,3,1);
        jsta = kbcinfo(nbl,nseg,4,1); jend = kbcinfo(nbl,nseg,5,1);
        ldata = lwdat(nbl,nseg,5);
        mdim = iend-ista; ndim = jend-jsta;
        std::memcpy(filname, bcfiles(bcfilek(nbl,nseg,1)), 80);
        nface = 5;

        if (kbcinfo(nbl,nseg,1,1) == 9999)
            out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,1) == 1000)
            out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,1) == 1001)
            out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,1) == 1002)
            out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,1) == 1003)
            out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,myid);
        if (kbcinfo(nbl,nseg,1,1) == 1005) { int inm=0;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (kbcinfo(nbl,nseg,1,1) == 1006) { int inm=1;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (kbcinfo(nbl,nseg,1,1) == 1008)
            out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,1) == 1011)
            out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,1) == 1012)
            out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,1) == 1013)
            out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,1) == 2002) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2003) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2009) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=0;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2010) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=1;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2019) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2004) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2014) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2024) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(kbcinfo(nbl,nseg,1,1)) == 2034) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2005) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            jdimp = jdimg(ngnew); kdimp = kdimg(ngnew); idimp = idimg(ngnew);
            nblp = ngnew;
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,nbl,jdimp,kdimp,idimp,nblp,filname,myid,mblk2nd,maxbl);
            for (lcntp = is_prd[level-1]; lcntp <= ie_prd[level-1]; lcntp++) {
                nbl_chk   = isav_prd(lcntp,1);
                nface_chk = isav_prd(lcntp,2);
                nseg_chk  = isav_prd(lcntp,11);
                if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                    if ((float)period_miss(lcntp) > (float)epsrot) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "           blocking check....periodic mismatch = %14.7E", (float)period_miss(lcntp));
                    }
                }
            }
        }
        if (kbcinfo(nbl,nseg,1,1) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = 0;
            if (nblnum > 0) ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,ngnew,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2007) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2008) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2016) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=5;
            out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2026) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2018) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2028) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2038) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2102) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=0;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        if (kbcinfo(nbl,nseg,1,1) == 2103) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=1;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        for (int lcnt = is_blk[level-1]; lcnt <= ie_blk[level-1]; lcnt++) {
            nbl_chk   = isav_blk(lcnt,1);
            nface_chk = isav_blk(lcnt,2);
            nseg_chk  = isav_blk(lcnt,17);
            if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                if ((float)geom_miss(lcnt) > (float)epsc0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "           blocking check....geometry mismatch = %14.7E", (float)geom_miss(lcnt));
                }
            }
        }
    } // end do k=0 face


    // k=kdim boundary (face 6)
    for (nseg = 1; nseg <= nbckdim(nbl); nseg++) {
        ksta = kdim; kend = kdim;
        ista = kbcinfo(nbl,nseg,2,2); iend = kbcinfo(nbl,nseg,3,2);
        jsta = kbcinfo(nbl,nseg,4,2); jend = kbcinfo(nbl,nseg,5,2);
        ldata = lwdat(nbl,nseg,6);
        mdim = iend-ista; ndim = jend-jsta;
        std::memcpy(filname, bcfiles(bcfilek(nbl,nseg,2)), 80);
        nface = 6;

        if (kbcinfo(nbl,nseg,1,2) == 9999)
            out9999(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,2) == 1000)
            out1000(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,2) == 1001)
            out1001(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,2) == 1002)
            out1002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,2) == 1003)
            out1003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,myid);
        if (kbcinfo(nbl,nseg,1,2) == 1005) { int inm=0;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (kbcinfo(nbl,nseg,1,2) == 1006) { int inm=1;
            out1005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,inm); }
        if (kbcinfo(nbl,nseg,1,2) == 1008)
            out1008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,2) == 1011)
            out1011(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,2) == 1012)
            out1012(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,2) == 1013)
            out1013(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim);
        if (kbcinfo(nbl,nseg,1,2) == 2002) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2002(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2003) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2003(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2009) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=0;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2010) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflgg=1;
            out2009(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflgg);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2019) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2019(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2004) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2014) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2024) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (std::abs(kbcinfo(nbl,nseg,1,2)) == 2034) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2004(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2005) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            jdimp = jdimg(ngnew); kdimp = kdimg(ngnew); idimp = idimg(ngnew);
            nblp = ngnew;
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2005(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,nbl,jdimp,kdimp,idimp,nblp,filname,myid,mblk2nd,maxbl);
            for (lcntp = is_prd[level-1]; lcntp <= ie_prd[level-1]; lcntp++) {
                nbl_chk   = isav_prd(lcntp,1);
                nface_chk = isav_prd(lcntp,2);
                nseg_chk  = isav_prd(lcntp,11);
                if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                    if ((float)period_miss(lcntp) > (float)epsrot) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, "           blocking check....periodic mismatch = %14.7E", (float)period_miss(lcntp));
                    }
                }
            }
        }
        if (kbcinfo(nbl,nseg,1,2) == 2006) {
            ngh = igridg(nbl);
            nblnum = (int)w(ldata);
            ngnew = 0;
            if (nblnum > 0) ngnew = nblg(nblnum) + (nbl - nblg(ngh));
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2006(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,ngnew,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2007) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2007(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2008) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=0;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2016) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=5;
            out2016(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2026) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            out2026(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2018) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=1;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2028) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=2;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2038) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iout=3;
            out2008(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iout);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2102) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=0;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        if (kbcinfo(nbl,nseg,1,2) == 2103) {
            FortranArray4DRef<double> bcdata(reinterpret_cast<double*>(&w(ldata)), mdim, ndim, 2, 12);
            int iflag=1;
            out2102(jdim,kdim,idim,ista,iend,jsta,jend,ksta,kend,nface,nou,bou,nbuf,ibufdim,mdim,ndim,bcdata,filname,myid,mblk2nd,maxbl,iflag);
        }
        for (int lcnt = is_blk[level-1]; lcnt <= ie_blk[level-1]; lcnt++) {
            nbl_chk   = isav_blk(lcnt,1);
            nface_chk = isav_blk(lcnt,2);
            nseg_chk  = isav_blk(lcnt,17);
            if (nbl == nbl_chk && nface == nface_chk && nseg == nseg_chk) {
                if ((float)geom_miss(lcnt) > (float)epsc0) {
                    nou(1) = std::min(nou(1)+1, ibufdim);
                    std::snprintf(bou(nou(1),1), 120, "           blocking check....geometry mismatch = %14.7E", (float)geom_miss(lcnt));
                }
            }
        }
    } // end do k=kdim face
} // end bc_info


// ============================================================
// out1000 - freestream type 1000
// ============================================================
void out1000(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    int& noninflag = cmn_noninertial.noninflag;
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  freestream                     type 1000  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        if (noninflag > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   j=   1  NONINERTIAL freestream         type 1000  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        }
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  freestream                     type 1000  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        if (noninflag > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   j=jdim  NONINERTIAL freestream         type 1000  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        }
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  freestream                     type 1000  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        if (noninflag > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   k=   1  NONINERTIAL freestream         type 1000  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        }
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  freestream                     type 1000  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        if (noninflag > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   k=kdim  NONINERTIAL freestream         type 1000  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        }
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  freestream                     type 1000  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        if (noninflag > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   i=   1  NONINERTIAL freestream         type 1000  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        }
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  freestream                     type 1000  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        if (noninflag > 0) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   i=idim  NONINERTIAL freestream         type 1000  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        }
    }
}

// ============================================================
// out1001 - symmetry plane type 1001
// ============================================================
void out1001(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  symmetry plane                 type 1001  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  symmetry plane                 type 1001  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  symmetry plane                 type 1001  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  symmetry plane                 type 1001  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  symmetry plane                 type 1001  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  symmetry plane                 type 1001  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
}

// ============================================================
// out1002 - 1-pt extrapolation type 1002
// ============================================================
void out1002(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  1-pt extrapolation             type 1002  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  1-pt extrapolation             type 1002  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  1-pt extrapolation             type 1002  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  1-pt extrapolation             type 1002  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  1-pt extrapolation             type 1002  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  1-pt extrapolation             type 1002  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
}


// ============================================================
// out1003 - characteristic inflow/outflow type 1003
// ============================================================
void out1003(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& myid)
{
    int& iipv = cmn_info.iipv;
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  characteristic inflow/outflow  type 1003  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        if (iipv == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           With farfield point vortex correction");
        }
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  characteristic inflow/outflow  type 1003  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        if (iipv == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           With farfield point vortex correction");
        }
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  characteristic inflow/outflow  type 1003  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        if (iipv == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           With farfield point vortex correction");
        }
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  characteristic inflow/outflow  type 1003  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        if (iipv == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           With farfield point vortex correction");
        }
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  characteristic inflow/outflow  type 1003  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        if (iipv == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " i=1 bc1003 not appropriate for i2d=-1 (farfield point vortex correction)");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  characteristic inflow/outflow  type 1003  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        if (iipv == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " i=1 bc1003 not appropriate for i2d=-1 (farfield point vortex correction)");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }
}

// ============================================================
// out1005 - inviscid surface type 1005/1006
// ============================================================
void out1005(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& inormmom)
{
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        if (inormmom == 0)
            std::snprintf(bou(nou(1),1), 120, "   j=   1  inviscid surface               type 1005  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        else
            std::snprintf(bou(nou(1),1), 120, "   j=   1  inviscid surface w norm mom eq type 1006  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        if (inormmom == 0)
            std::snprintf(bou(nou(1),1), 120, "   j=jdim  inviscid surface               type 1005  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        else
            std::snprintf(bou(nou(1),1), 120, "   j=jdim  inviscid surface w norm mom eq type 1006  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        if (inormmom == 0)
            std::snprintf(bou(nou(1),1), 120, "   k=   1  inviscid surface               type 1005  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        else
            std::snprintf(bou(nou(1),1), 120, "   k=   1  inviscid surface w norm mom eq type 1006  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        if (inormmom == 0)
            std::snprintf(bou(nou(1),1), 120, "   k=kdim  inviscid surface               type 1005  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        else
            std::snprintf(bou(nou(1),1), 120, "   k=kdim  inviscid surface w norm mom eq type 1006  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        if (inormmom == 0)
            std::snprintf(bou(nou(1),1), 120, "   i=   1  inviscid surface               type 1005  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        else
            std::snprintf(bou(nou(1),1), 120, "   i=   1  inviscid surface w norm mom eq type 1006  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        if (inormmom == 0)
            std::snprintf(bou(nou(1),1), 120, "   i=idim  inviscid surface               type 1005  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        else
            std::snprintf(bou(nou(1),1), 120, "   i=idim  inviscid surface w norm mom eq type 1006  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
}

// ============================================================
// out1008 - tunnel inflow type 1008
// ============================================================
void out1008(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  tunnel inflow                  type 1008  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  tunnel inflow                  type 1008  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  tunnel inflow                  type 1008  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  tunnel inflow                  type 1008  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  tunnel inflow                  type 1008  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  tunnel inflow                  type 1008  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
}


// ============================================================
// out1011 - singular axis half plane type 1011
// ============================================================
void out1011(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  singular axis - half plane     type 1011  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  singular axis - half plane     type 1011  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  singular axis - half plane     type 1011  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  singular axis - half plane     type 1011  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  singular axis - half plane     type 1011  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  singular axis - half plane     type 1011  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
}

// ============================================================
// out1012 - singular axis full plane type 1012
// ============================================================
void out1012(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  singular axis - full plane     type 1012  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  singular axis - full plane     type 1012  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  singular axis - full plane     type 1012  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  singular axis - full plane     type 1012  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  singular axis - full plane     type 1012  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  singular axis - full plane     type 1012  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
}

// ============================================================
// out1013 - singular axis extrapolation type 1013
// ============================================================
void out1013(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  singular axis - extrapolation  type 1013  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  singular axis - extrapolation  type 1013  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  singular axis - extrapolation  type 1013  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  singular axis - extrapolation  type 1013  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  singular axis - extrapolation  type 1013  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  singular axis - extrapolation  type 1013  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
}

// ============================================================
// out2002 - set p, extrapolate rho,u,v,w type 2002
// ============================================================
void out2002(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int iflg;
    float pratio;
    char filname_str[61];

    auto do_check = [&]() {
        iflg = 0;
        for (int ipp=1; ipp<=2; ipp++)
            for (int d1=1; d1<=mdim; d1++)
                for (int d2=1; d2<=ndim; d2++)
                    if ((float)bcdata(d1,d2,ipp,1) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2002:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  pressure data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    };
    auto do_print_file = [&]() {
        if (strncmp(filname, "null", 4) == 0) {
            pratio = (float)bcdata(1,1,1,1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           P/Pinf     = %8.4f", pratio);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           P/Pinf set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    if (nface == 3) {
        do_check();
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  set p, extrapolate rho,u,v,w   type 2002  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        do_print_file();
    }
    if (nface == 4) {
        do_check();
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  set p, extrapolate rho,u,v,w   type 2002  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        do_print_file();
    }
    if (nface == 5) {
        do_check();
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  set p, extrapolate rho,u,v,w   type 2002  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        do_print_file();
    }
    if (nface == 6) {
        do_check();
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  set p, extrapolate rho,u,v,w   type 2002  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        do_print_file();
    }
    if (nface == 1) {
        do_check();
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  set p, extrapolate rho,u,v,w   type 2002  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        do_print_file();
    }
    if (nface == 2) {
        do_check();
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  set p, extrapolate rho,u,v,w   type 2002  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        do_print_file();
    }
}


// ============================================================
// out2003 - engine inflow type 2003
// ============================================================
void out2003(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int& ivmx = cmn_maxiv.ivmx;
    int iflg, itrflg1, itrflg2;
    float xme, pte, tte, alpe, betae, t1, t2;
    char filname_str[61];

    auto do_check_and_print = [&](const char* hdr) {
        int ipp = 1;
        itrflg1 = 0; itrflg2 = 0;
        if ((float)bcdata(1,1,ipp,6) >= 0.f) itrflg1 = 1;
        if ((float)bcdata(1,1,ipp,7) >= 0.f) itrflg2 = 1;
        iflg = 0;
        for (int l=1; l<=5; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2003:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  engine inflow data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (itrflg1 > 0 || itrflg2 > 0) {
            if (ivmx > 5 && itrflg1*itrflg2 == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  stopping in bc2003...must set both turbulence");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  quantities for 2-eq turb. models...only one is set");
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
        }
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "%s", hdr);
        if (strncmp(filname, "null", 4) == 0) {
            xme   = (float)bcdata(1,1,1,1);
            pte   = (float)bcdata(1,1,1,2);
            tte   = (float)bcdata(1,1,1,3);
            alpe  = (float)bcdata(1,1,1,4);
            betae = (float)bcdata(1,1,1,5);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           inlet mach     = %8.4f", xme);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           Pt/Pinf        = %8.4f", pte);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           Tt/Tinf        = %8.4f", tte);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           alphe          = %8.4f", alpe);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           betae          = %8.4f", betae);
            if (itrflg1 > 0) {
                t1 = (float)bcdata(1,1,1,6);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb1 (nondim) = %12.4E", t1);
            }
            if (itrflg2 > 0) {
                t2 = (float)bcdata(1,1,1,7);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb2 (nondim) = %12.4E", t2);
            }
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           inflow data set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    char hdr[120];
    if (nface == 3) {
        std::snprintf(hdr, 120, "   j=   1  engine inflow                  type 2003  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        do_check_and_print(hdr);
    }
    if (nface == 4) {
        std::snprintf(hdr, 120, "   j=jdim  engine inflow                  type 2003  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
        do_check_and_print(hdr);
    }
    if (nface == 5) {
        std::snprintf(hdr, 120, "   k=   1  engine inflow                  type 2003  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        do_check_and_print(hdr);
    }
    if (nface == 6) {
        std::snprintf(hdr, 120, "   k=kdim  engine inflow                  type 2003  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
        do_check_and_print(hdr);
    }
    if (nface == 1) {
        std::snprintf(hdr, 120, "   i=   1  engine inflow                  type 2003  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        do_check_and_print(hdr);
    }
    if (nface == 2) {
        std::snprintf(hdr, 120, "   i=idim  engine inflow                  type 2003  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
        do_check_and_print(hdr);
    }
}


// ============================================================
// out2004 - viscous wall type 2004/2014/2024/2034
// ============================================================
void out2004(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iout)
{
    int* iwf = cmn_wallfun.iwf;
    float& uub_2034 = cmn_wallvel_2034.uub_2034;
    float& vvb_2034 = cmn_wallvel_2034.vvb_2034;
    float& wwb_2034 = cmn_wallvel_2034.wwb_2034;
    float& vvb_xaxisrot2034 = cmn_wallvel_2034.vvb_xaxisrot2034;

    int iextra = 0;
    if (iout == 1 || iout == 2) iextra = 1;

    int iend1 = iend-1, jend1 = jend-1, kend1 = kend-1;
    int iflg;
    float ctemp, cq, stopindex, gammawall;
    char filname_str[61];

    // Helper: print the face header based on nface and iout
    auto print_hdr = [&](int a1, int a2, int a3, int a4, bool is_ij) {
        // is_ij: true if indices are i,j; false if i,k or j,k
        const char* types[4] = {"2004","2014","2024","2034"};
        const char* descs[4] = {"viscous wall                  ","viscous laminar wall          ","viscous laminar wall          ","viscous wall + imposed vel    "};
        const char* faces[6] = {"j=   1","j=jdim","k=   1","k=kdim","i=   1","i=idim"};
        int fi = nface-1;
        nou(1) = std::min(nou(1)+1, ibufdim);
        if (nface == 3 || nface == 4) {
            std::snprintf(bou(nou(1),1), 120, "   %s  %s type %s  i=%5d,%5d  k=%5d,%5d",
                faces[fi], descs[iout], types[iout], a1,a2,a3,a4);
        } else if (nface == 5 || nface == 6) {
            std::snprintf(bou(nou(1),1), 120, "   %s  %s type %s  i=%5d,%5d  j=%5d,%5d",
                faces[fi], descs[iout], types[iout], a1,a2,a3,a4);
        } else {
            std::snprintf(bou(nou(1),1), 120, "   %s  %s type %s  j=%5d,%5d  k=%5d,%5d",
                faces[fi], descs[iout], types[iout], a1,a2,a3,a4);
        }
        if (iout == 3) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           Imposed u,v,w,vrot_about_x_axis = %10.5f%10.5f%10.5f%10.5f",
                uub_2034, vvb_2034, wwb_2034, vvb_xaxisrot2034);
        }
    };

    auto do_face = [&](int a1, int a2, int a3, int a4) {
        iflg = 0;
        for (int l=1; l<=2+iextra; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            if (iout == 0) std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2004:");
            else if (iout == 1) std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2014:");
            else if (iout == 2) std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2024:");
            else std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2034:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  viscous wall data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (strncmp(filname, "null", 4) == 0) {
            ctemp = (float)bcdata(1,1,1,1);
            cq    = (float)bcdata(1,1,1,2);
            if (iout == 1) stopindex = (float)bcdata(1,1,1,3);
            if (iout == 2) gammawall = (float)bcdata(1,1,1,3);
            print_hdr(a1,a2,a3,a4,false);
            if ((float)ctemp > 0.f) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           Tw/Tinf    = %8.4f", ctemp);
            } else if (ctemp == 0.f) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           Tw/Tinf    = adiabatic wall");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           Tw/Tinf    = stagnation");
            }
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           C_q        = %8.4f", cq);
            if (iout == 1) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           stopindex  = %8d", (int)stopindex);
            }
            if (iout == 2) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           intermtbc  = %8.4f", gammawall);
            }
        } else {
            print_hdr(a1,a2,a3,a4,false);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           Twtype, C_q set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    if (nface == 3) { do_face(ista,iend,ksta,kend); if (iwf[2-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           Wall function employed in j-direction"); } }
    if (nface == 4) { do_face(ista,iend,ksta,kend); if (iwf[2-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           Wall function employed in j-direction"); } }
    if (nface == 5) { do_face(ista,iend,jsta,jend); if (iwf[3-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           Wall function employed in k-direction"); } }
    if (nface == 6) { do_face(ista,iend,jsta,jend); if (iwf[3-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           Wall function employed in k-direction"); } }
    if (nface == 1) { do_face(jsta,jend,ksta,kend); if (iwf[1-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           Wall function employed in i-direction"); } }
    if (nface == 2) { do_face(jsta,jend,ksta,kend); if (iwf[1-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           Wall function employed in i-direction"); } }
}


// ============================================================
// out2005 - periodic BC type 2005
// ============================================================
void out2005(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, int& nbl,
             int& jdimp, int& kdimp, int& idimp, int& nblp,
             char (&filname)[80], int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int iflg, ip;
    float dthtx, dthty, dthtz;
    char filname_str[61];

    auto do_check = [&]() {
        iflg = 0;
        for (int l=1; l<=4; l++)
            for (int ip2=1; ip2<=2; ip2++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip2,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2005:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  periodic BC data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    };

    auto do_dtht_check = [&]() {
        ip = 1;
        for (int l=2; l<=4; l++)
            for (int d1=1; d1<=mdim; d1++)
                for (int d2=1; d2<=ndim; d2++)
                    if (bcdata(d1,d2,ip,l) != bcdata(1,1,ip,l)) {
                        nou(1) = std::min(nou(1)+1, ibufdim);
                        std::snprintf(bou(nou(1),1), 120, " dthtx,y,z currently must be const over whole face");
                        int m1 = -1;
                        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
                    }
        if (bcdata(1,1,ip,2) != 0. && (bcdata(1,1,ip,3) != 0. || bcdata(1,1,ip,4) != 0.)) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " 2 of the 3 dtht values currently must = 0");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (bcdata(1,1,ip,3) != 0. && (bcdata(1,1,ip,2) != 0. || bcdata(1,1,ip,4) != 0.)) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " 2 of the 3 dtht values currently must = 0");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (bcdata(1,1,ip,4) != 0. && (bcdata(1,1,ip,2) != 0. || bcdata(1,1,ip,3) != 0.)) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " 2 of the 3 dtht values currently must = 0");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    };

    if (nface == 3) {
        do_check();
        if (strncmp(filname, "null", 4) == 0) {
            dthtx = (float)bcdata(1,1,1,2); dthty = (float)bcdata(1,1,1,3); dthtz = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   j=   1  periodic - angular rotation    type 2005  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodic with j=jdim, block%5d", nblp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetax    = %8.4f", dthtx);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetay    = %8.4f", dthty);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetaz    = %8.4f", dthtz);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   j=   1  periodic - angular rotation    type 2005  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodicity set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
        do_dtht_check();
        if (kdim != kdimp || idim != idimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " idim,kdim,idimp,kdimp=%5d%5d%5d%5d", idim,kdim,idimp,kdimp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Block you are periodic with currently must be of same 2 dimensions");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " (and orientation) on the face as the BC block");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (nbl != nblp && jdimp == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Periodic block cannot be dimension 2 in j");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  (on ANY multigrid level)");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }
    if (nface == 4) {
        do_check();
        if (strncmp(filname, "null", 4) == 0) {
            dthtx = (float)bcdata(1,1,1,2); dthty = (float)bcdata(1,1,1,3); dthtz = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   j=jdim  periodic - angular rotation    type 2005  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodic with j=   1, block%5d", nblp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetax    = %8.4f", dthtx);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetay    = %8.4f", dthty);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetaz    = %8.4f", dthtz);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   j=jdim  periodic - angular rotation    type 2005  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodicity set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
        do_dtht_check();
        if (kdim != kdimp || idim != idimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " idim,kdim,idimp,kdimp=%5d%5d%5d%5d", idim,kdim,idimp,kdimp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Block you are periodic with currently must be of same 2 dimensions");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " (and orientation) on the face as the BC block");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (nbl != nblp && jdimp == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Periodic block cannot be dimension 2 in j");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  (on ANY multigrid level)");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }

    if (nface == 5) {
        do_check();
        if (strncmp(filname, "null", 4) == 0) {
            dthtx = (float)bcdata(1,1,1,2); dthty = (float)bcdata(1,1,1,3); dthtz = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   k=   1  periodic - angular rotation    type 2005  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodic with k=kdim, block%5d", nblp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetax    = %8.4f", dthtx);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetay    = %8.4f", dthty);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetaz    = %8.4f", dthtz);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   k=   1  periodic - angular rotation    type 2005  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodicity set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
        do_dtht_check();
        if (jdim != jdimp || idim != idimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " idim,jdim,idimp,jdimp=%5d%5d%5d%5d", idim,jdim,idimp,jdimp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Block you are periodic with currently must be of same 2 dimensions");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " (and orientation) on the face as the BC block");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (nbl != nblp && kdimp == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Periodic block cannot be dimension 2 in k");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  (on ANY multigrid level)");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }
    if (nface == 6) {
        do_check();
        if (strncmp(filname, "null", 4) == 0) {
            dthtx = (float)bcdata(1,1,1,2); dthty = (float)bcdata(1,1,1,3); dthtz = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   k=kdim  periodic - angular rotation    type 2005  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodic with k=   1, block%5d", nblp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetax    = %8.4f", dthtx);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetay    = %8.4f", dthty);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetaz    = %8.4f", dthtz);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   k=kdim  periodic - angular rotation    type 2005  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodicity set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
        do_dtht_check();
        if (jdim != jdimp || idim != idimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " idim,jdim,idimp,jdimp=%5d%5d%5d%5d", idim,jdim,idimp,jdimp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Block you are periodic with currently must be of same 2 dimensions");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " (and orientation) on the face as the BC block");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (nbl != nblp && kdimp == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Periodic block cannot be dimension 2 in k");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  (on ANY multigrid level)");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }
    if (nface == 1) {
        do_check();
        if (strncmp(filname, "null", 4) == 0) {
            dthtx = (float)bcdata(1,1,1,2); dthty = (float)bcdata(1,1,1,3); dthtz = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   i=   1  periodic - angular rotation    type 2005  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodic with i=idim, block%5d", nblp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetax    = %8.4f", dthtx);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetay    = %8.4f", dthty);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetaz    = %8.4f", dthtz);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   i=   1  periodic - angular rotation    type 2005  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodicity set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
        do_dtht_check();
        if (jdim != jdimp || kdim != kdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " jdim,kdim,jdimp,kdimp=%5d%5d%5d%5d", jdim,kdim,jdimp,kdimp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Block you are periodic with currently must be of same 2 dimensions");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " (and orientation) on the face as the BC block");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (nbl != nblp && idimp == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Periodic block cannot be dimension 2 in i");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  (on ANY multigrid level)");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }
    if (nface == 2) {
        do_check();
        if (strncmp(filname, "null", 4) == 0) {
            dthtx = (float)bcdata(1,1,1,2); dthty = (float)bcdata(1,1,1,3); dthtz = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   i=idim  periodic - angular rotation    type 2005  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodic with i=   1, block%5d", nblp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetax    = %8.4f", dthtx);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetay    = %8.4f", dthty);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           dthetaz    = %8.4f", dthtz);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "   i=idim  periodic - angular rotation    type 2005  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           periodicity set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
        do_dtht_check();
        if (jdim != jdimp || kdim != kdimp) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " jdim,kdim,jdimp,kdimp=%5d%5d%5d%5d", jdim,kdim,jdimp,kdimp);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Block you are periodic with currently must be of same 2 dimensions");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " (and orientation) on the face as the BC block");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (nbl != nblp && idimp == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, " Periodic block cannot be dimension 2 in i");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  (on ANY multigrid level)");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    }
}


// ============================================================
// out2006 - radial equilibrium type 2006
// ============================================================
void out2006(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& nblc, int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int iflg, lijk, ldir;
    float pratio;
    char filname_str[61];

    if ((float)bcdata(1,1,1,3) == 0.f) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, " stopping...bcdata(3) = 0, should be +/-1, +/-2, or +/-3");
        int m1 = -1;
        termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
    }
    lijk = std::abs((int)(float)bcdata(1,1,1,3));
    ldir = (int)(float)bcdata(1,1,1,3) / lijk;

    auto do_check = [&](int nlayers) {
        iflg = 0;
        for (int l=1; l<=nlayers; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++) {
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
                        if (std::abs((int)(float)bcdata(d1,d2,ip,3)) > 3) iflg = 1;
                    }
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2006:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  radial equilib. press. data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    };

    auto do_lijk_dir = [&]() {
        if (lijk == 1 && ldir == +1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           integrated in +i direction"); }
        if (lijk == 1 && ldir == -1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           integrated in -i direction"); }
        if (lijk == 2 && ldir == +1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           integrated in +j direction"); }
        if (lijk == 2 && ldir == -1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           integrated in -j direction"); }
        if (lijk == 3 && ldir == +1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           integrated in +k direction"); }
        if (lijk == 3 && ldir == -1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"           integrated in -k direction"); }
    };

    auto do_print = [&](const char* hdr1, const char* hdr2, int a1, int a2, int a3, int a4) {
        if (strncmp(filname, "null", 4) == 0) {
            if ((float)bcdata(1,1,1,1) <= 0.f) {
                pratio = (float)bcdata(1,1,1,2);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "%s", hdr1);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           P/Pinf     = %8.4f", pratio);
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "%s", hdr2);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           continued from block%5d", nblc);
            }
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s", hdr1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           P/Pinf set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
        do_lijk_dir();
    };

    char hdr1[120], hdr2[120];
    if (nface == 3) {
        if (lijk != 1 && lijk != 3) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120," stopping...abs(bcdata(3)) must be 1 or 3 on a j=constant face");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
        do_check(4);
        std::snprintf(hdr1,120,"   j=   1  set p via radial equilibrium   type 2006  i=%5d,%5d  k=%5d,%5d",ista,iend,ksta,kend);
        std::snprintf(hdr2,120,"   j=   1  set p via radial equilibrium   type 2006  i=%5d,%5d  k=%5d,%5d",ista,iend,ksta,kend);
        do_print(hdr1,hdr2,ista,iend,ksta,kend);
        if (nblc != 0 && (ista != 1 || iend != idim)) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120,"  stopping...bc segment must span the entire block face in the radial direction");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
    }
    if (nface == 4) {
        if (lijk != 1 && lijk != 3) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120," stopping...abs(bcdata(3)) must be 1 or 3 on a j=constant face");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
        do_check(2);
        std::snprintf(hdr1,120,"   j=jdim  set p via radial equilibrium   type 2006  i=%5d,%5d  k=%5d,%5d",ista,iend,ksta,kend);
        std::snprintf(hdr2,120,"   j=jdim  set p via radial equilibrium   type 2006  i=%5d,%5d  k=%5d,%5d",ista,iend,ksta,kend);
        do_print(hdr1,hdr2,ista,iend,ksta,kend);
        if (nblc != 0 && (ista != 1 || iend != idim)) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120,"  stopping...bc segment must span the entire block face in the radial direction");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
    }
    if (nface == 5) {
        if (lijk != 1 && lijk != 2) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120," stopping...abs(bcdata(3)) must be 1 or 2 on a k=constant face");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
        do_check(4);
        std::snprintf(hdr1,120,"   k=   1  set p via radial equilibrium   type 2006  i=%5d,%5d  j=%5d,%5d",ista,iend,jsta,jend);
        std::snprintf(hdr2,120,"   k=   1  set p via radial equilibrium   type 2006  i=%5d,%5d  j=%5d,%5d",ista,iend,jsta,jend);
        do_print(hdr1,hdr2,ista,iend,jsta,jend);
        if (nblc != 0 && (ista != 1 || iend != idim)) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120,"  stopping...bc segment must span the entire block face in the radial direction");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
    }
    if (nface == 6) {
        if (lijk != 1 && lijk != 2) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120," stopping...abs(bcdata(3)) must be 1 or 2 on a k=constant face");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
        do_check(4);
        std::snprintf(hdr1,120,"   k=kdim  set p via radial equilibrium   type 2006  i=%5d,%5d  j=%5d,%5d",ista,iend,jsta,jend);
        std::snprintf(hdr2,120,"   k=kdim  set p via radial equilibrium   type 2006  i=%5d,%5d  j=%5d,%5d",ista,iend,jsta,jend);
        do_print(hdr1,hdr2,ista,iend,jsta,jend);
        if (nblc != 0 && (ista != 1 || iend != idim)) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120,"  stopping...bc segment must span the entire block face in the radial direction");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
    }
    if (nface == 1) {
        if (lijk != 2 && lijk != 3) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120," stopping...abs(bcdata(3)) must be 2 or 3 on an i=constant face");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
        do_check(4);
        std::snprintf(hdr1,120,"   i=   1  set p via radial equilibrium   type 2006  j=%5d,%5d  k=%5d,%5d",jsta,jend,ksta,kend);
        std::snprintf(hdr2,120,"   i=   1  set p via radial equilibrium   type 2006  j=%5d,%5d  k=%5d,%5d",jsta,jend,ksta,kend);
        do_print(hdr1,hdr2,jsta,jend,ksta,kend);
        if (nblc != 0 && (jsta != 1 || jend != jdim)) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120,"  stopping...bc segment must span the entire block face in the radial direction");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
    }
    if (nface == 2) {
        if (lijk != 2 && lijk != 3) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120," stopping...abs(bcdata(3)) must be 2 or 3 on an i=constant face");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
        do_check(4);
        std::snprintf(hdr1,120,"   i=idim  set p via radial equilibrium   type 2006  j=%5d,%5d  k=%5d,%5d",jsta,jend,ksta,kend);
        std::snprintf(hdr2,120,"   i=idim  set p via radial equilibrium   type 2006  j=%5d,%5d  k=%5d,%5d",jsta,jend,ksta,kend);
        do_print(hdr1,hdr2,jsta,jend,ksta,kend);
        if (nblc != 0 && (jsta != 1 || jend != jdim)) {
            nou(1)=std::min(nou(1)+1,ibufdim);
            std::snprintf(bou(nou(1),1),120,"  stopping...bc segment must span the entire block face in the radial direction");
            int m1=-1; termn8_ns::termn8(myid,m1,ibufdim,nbuf,bou,nou);
        }
    }
}


// ============================================================
// out2007 - primitive variable inflow type 2007
// ============================================================
void out2007(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int& ivmx = cmn_maxiv.ivmx;
    int iflg, itrflg1, itrflg2;
    float rho1, u1, v1, w1, p1, t1, t2;
    char filname_str[61];

    auto do_check_and_print = [&](const char* hdr_null, const char* hdr_file) {
        int ipp = 1;
        itrflg1 = 0; itrflg2 = 0;
        if ((float)bcdata(1,1,ipp,6) >= 0.f) itrflg1 = 1;
        if ((float)bcdata(1,1,ipp,7) >= 0.f) itrflg2 = 1;
        iflg = 0;
        for (int l=1; l<=5; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2007:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  primative variable data incorrectly set");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (itrflg1 > 0 || itrflg2 > 0) {
            if (ivmx > 5 && itrflg1*itrflg2 == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  stopping in bc2007...must set both turbulence");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  quantities for 2-eq turb. models...only one is set");
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
        }
        if (strncmp(filname, "null", 4) == 0) {
            rho1 = (float)bcdata(1,1,1,1);
            u1   = (float)bcdata(1,1,1,2);
            v1   = (float)bcdata(1,1,1,3);
            w1   = (float)bcdata(1,1,1,4);
            p1   = (float)bcdata(1,1,1,5);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s", hdr_null);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           rho/rhoinf  = %8.4f", rho1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           u/ainf      = %8.4f", u1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           v/ainf      = %8.4f", v1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           w/ainf      = %8.4f", w1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           p/pinf      = %8.4f", p1);
            if (itrflg1 > 0) {
                t1 = (float)bcdata(1,1,1,6);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb1 (nondim) = %12.4E", t1);
            }
            if (itrflg2 > 0) {
                t2 = (float)bcdata(1,1,1,7);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb2 (nondim) = %12.4E", t2);
            }
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s", hdr_file);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           inflow data set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    char hdr1[120], hdr2[120];
    if (nface == 3) {
        std::snprintf(hdr1,120,"   j=   1  set rho,u,v,w,p (primative)    type 2007  i=%5d,%5d  k=%5d,%5d",ista,iend,ksta,kend);
        std::snprintf(hdr2,120,"   j=   1  set rho,u,v,w,p (primative)    type 2007  i=%5d,%5d  k=%5d,%5d",ista,iend,ksta,kend);
        do_check_and_print(hdr1,hdr2);
    }
    if (nface == 4) {
        std::snprintf(hdr1,120,"   j=jdim  set rho,u,v,w,p (primative)    type 2007  i=%5d,%5d  k=%5d,%5d",ista,iend,ksta,kend);
        std::snprintf(hdr2,120,"   j=jdim  set rho,u,v,w,p (primative)    type 2007  i=%5d,%5d  k=%5d,%5d",ista,iend,ksta,kend);
        do_check_and_print(hdr1,hdr2);
    }
    if (nface == 5) {
        std::snprintf(hdr1,120,"   k=   1  set rho,u,v,w,p (primative)    type 2007  i=%5d,%5d  j=%5d,%5d",ista,iend,jsta,jend);
        std::snprintf(hdr2,120,"   k=   1  set rho,u,v,w,p (primative)    type 2007  i=%5d,%5d  j=%5d,%5d",ista,iend,jsta,jend);
        do_check_and_print(hdr1,hdr2);
    }
    if (nface == 6) {
        std::snprintf(hdr1,120,"   k=kdim  set rho,u,v,w,p (primative)    type 2007  i=%5d,%5d  j=%5d,%5d",ista,iend,jsta,jend);
        std::snprintf(hdr2,120,"   k=kdim  set rho,u,v,w,p (primative)    type 2007  i=%5d,%5d  j=%5d,%5d",ista,iend,jsta,jend);
        do_check_and_print(hdr1,hdr2);
    }
    if (nface == 1) {
        std::snprintf(hdr1,120,"   i=   1  set rho,u,v,w,p (primative)    type 2007  j=%5d,%5d  k=%5d,%5d",jsta,jend,ksta,kend);
        std::snprintf(hdr2,120,"   i=   1  set rho,u,v,w,p (primative)    type 2007  j=%5d,%5d  k=%5d,%5d",jsta,jend,ksta,kend);
        do_check_and_print(hdr1,hdr2);
    }
    if (nface == 2) {
        std::snprintf(hdr1,120,"   i=idim  set rho,u,v,w,p (primative)    type 2007  j=%5d,%5d  k=%5d,%5d",jsta,jend,ksta,kend);
        std::snprintf(hdr2,120,"   i=idim  set rho,u,v,w,p (primative)    type 2007  j=%5d,%5d  k=%5d,%5d",jsta,jend,ksta,kend);
        do_check_and_print(hdr1,hdr2);
    }
}


// ============================================================
// out2008 - set rho,u,v,w, extrapolate p type 2008/2018/2028/2038
// ============================================================
void out2008(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iout)
{
    int iflg, itrflg1, itrflg2;
    float val1, u1, v1, w1, t1, t2;
    char filname_str[61];

    // Format strings for each iout and face
    // iout=0: type 2008, iout=1: type 2018, iout=2: type 2028, iout=3: type 2038
    static const char* hdr_j1[4] = {
        "   j=   1  set rho,u,v,w, extrapolate p   type 2008",
        "   j=   1  set T,rhou,rhov,rhow, extrap p type 2018",
        "   j=   1 set f,rhou,rhov,rhow, ext rho,p type 2028",
        "   j=   1  set rho,u,v,w+rndm, extrap p   type 2038"
    };
    static const char* hdr_jd[4] = {
        "   j=jdim  set rho,u,v,w, extrapolate p   type 2008",
        "   j=jdim  set T,rhou,rhov,rhow, extrap p type 2018",
        "   j=jdim set f,rhou,rhov,rhow, ext rho,p type 2028",
        "   j=jdim  set rho,u,v,w+rndm, extrap p   type 2038"
    };
    static const char* hdr_k1[4] = {
        "   k=   1  set rho,u,v,w, extrapolate p   type 2008",
        "   k=   1  set T,rhou,rhov,rhow, extrap p type 2018",
        "   k=   1 set f,rhou,rhov,rhow, ext rho,p type 2028",
        "   k=   1  set rho,u,v,w+rndm, extrap p   type 2038"
    };
    static const char* hdr_kd[4] = {
        "   k=kdim  set rho,u,v,w, extrapolate p   type 2008",
        "   k=kdim  set T,rhou,rhov,rhow, extrap p type 2018",
        "   k=kdim set f,rhou,rhov,rhow, ext rho,p type 2028",
        "   k=kdim  set rho,u,v,w+rndm, extrap p   type 2038"
    };
    static const char* hdr_i1[4] = {
        "   i=   1  set rho,u,v,w, extrapolate p   type 2008",
        "   i=   1  set T,rhou,rhov,rhow, extrap p type 2018",
        "   i=   1 set f,rhou,rhov,rhow, ext rho,p type 2028",
        "   i=   1  set rho,u,v,w+rndm, extrap p   type 2038"
    };
    static const char* hdr_id[4] = {
        "   i=idim  set rho,u,v,w, extrapolate p   type 2008",
        "   i=idim  set T,rhou,rhov,rhow, extrap p type 2018",
        "   i=idim set f,rhou,rhov,rhow, ext rho,p type 2028",
        "   i=idim  set rho,u,v,w+rndm, extrap p   type 2038"
    };
    static const char* val_lbl[4] = {
        "rho/rhoinf     = ",
        "T/Tinf         = ",
        "freq*lref/aref    = ",
        "rho/rhoinf     = "
    };
    static const char* u_lbl[4] = {
        "u/ainf         = ",
        "rhou/rhoinfainf= ",
        "rhoumax/rhoinfainf= ",
        "u/ainf         = "
    };
    static const char* v_lbl[4] = {
        "v/ainf         = ",
        "rhov/rhoinfainf= ",
        "rhovmax/rhoinfainf= ",
        "v/ainf         = "
    };
    static const char* w_lbl[4] = {
        "w/ainf         = ",
        "rhow/rhoinfainf= ",
        "rhowmax/rhoinfainf= ",
        "w/ainf         = "
    };

    auto do_check = [&]() {
        iflg = 0;
        for (int l=1; l<=4; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            if (iout == 0) std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2008:");
            else if (iout == 1) std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2018:");
            else if (iout == 2) std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2028:");
            else std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2038:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  variable data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        itrflg1 = 0; itrflg2 = 0;
        if ((float)bcdata(1,1,1,5) >= 0.f) itrflg1 = 1;
        if ((float)bcdata(1,1,1,6) >= 0.f) itrflg2 = 1;
    };

    auto do_print = [&](const char* hdr) {
        if (strncmp(filname, "null", 4) == 0) {
            val1 = (float)bcdata(1,1,1,1);
            u1   = (float)bcdata(1,1,1,2);
            v1   = (float)bcdata(1,1,1,3);
            w1   = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s", hdr);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           %s%8.4f", val_lbl[iout], val1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           %s%8.4f", u_lbl[iout], u1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           %s%8.4f", v_lbl[iout], v1);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           %s%8.4f", w_lbl[iout], w1);
            if (itrflg1 > 0) {
                t1 = (float)bcdata(1,1,1,5);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb1 (nondim) = %12.4E", t1);
            }
            if (itrflg2 > 0) {
                t2 = (float)bcdata(1,1,1,6);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb2 (nondim) = %12.4E", t2);
            }
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s", hdr);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           inflow data set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    char hdr[120];
    if (nface == 3) {
        do_check();
        std::snprintf(hdr,120,"%s  i=%5d,%5d  k=%5d,%5d",hdr_j1[iout],ista,iend,ksta,kend);
        do_print(hdr);
    }
    if (nface == 4) {
        do_check();
        std::snprintf(hdr,120,"%s  i=%5d,%5d  k=%5d,%5d",hdr_jd[iout],ista,iend,ksta,kend);
        do_print(hdr);
    }
    if (nface == 5) {
        do_check();
        std::snprintf(hdr,120,"%s  i=%5d,%5d  j=%5d,%5d",hdr_k1[iout],ista,iend,jsta,jend);
        do_print(hdr);
    }
    if (nface == 6) {
        do_check();
        std::snprintf(hdr,120,"%s  i=%5d,%5d  j=%5d,%5d",hdr_kd[iout],ista,iend,jsta,jend);
        do_print(hdr);
    }
    if (nface == 1) {
        do_check();
        std::snprintf(hdr,120,"%s  j=%5d,%5d  k=%5d,%5d",hdr_i1[iout],jsta,jend,ksta,kend);
        do_print(hdr);
    }
    if (nface == 2) {
        do_check();
        std::snprintf(hdr,120,"%s  j=%5d,%5d  k=%5d,%5d",hdr_id[iout],jsta,jend,ksta,kend);
        do_print(hdr);
    }
}


// ============================================================
// out2009 - nozzle total BC type 2009/2010
// ============================================================
void out2009(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iflgg)
{
    int& ivmx = cmn_maxiv.ivmx;
    int iflg, itrflg1, itrflg2;
    float pte, tte, alpe, betae, t1, t2;
    char filname_str[61];

    // type 2009 headers (iflgg==0) and type 2010 headers (iflgg==1)
    static const char* hdr09[6] = {
        "   j=   1  nozzle total BC                type 2009",
        "   j=jdim  nozzle total BC                type 2009",
        "   k=   1  nozzle total BC                type 2009",
        "   k=kdim  nozzle total BC                type 2009",
        "   i=   1  nozzle total BC                type 2009",
        "   i=idim  nozzle total BC                type 2009"
    };
    static const char* hdr10[6] = {
        "   j=   1  nozzle total BC                type 2010",
        "   j=jdim  nozzle total BC                type 2010",
        "   k=   1  nozzle total BC                type 2010",
        "   k=kdim  nozzle total BC                type 2010",
        "   i=   1  nozzle total BC                type 2010",
        "   i=idim  nozzle total BC                type 2010"
    };

    auto do_check = [&]() {
        itrflg1 = 0; itrflg2 = 0;
        if ((float)bcdata(1,1,1,5) >= 0.f) itrflg1 = 1;
        if ((float)bcdata(1,1,1,6) >= 0.f) itrflg2 = 1;
        iflg = 0;
        for (int l=1; l<=4; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2009/2010:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  nozzle total BC data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (itrflg1 > 0 || itrflg2 > 0) {
            if (ivmx > 5 && itrflg1*itrflg2 == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  stopping in bc2009/2010...must set both turbulence");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  quantities for 2-eq turb. models...only one is set");
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
        }
    };

    auto do_print = [&](const char* hdr, int a1, int a2, int a3, int a4) {
        if (strncmp(filname, "null", 4) == 0) {
            pte   = (float)bcdata(1,1,1,1);
            tte   = (float)bcdata(1,1,1,2);
            alpe  = (float)bcdata(1,1,1,3);
            betae = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d",
                hdr, (nface<=4?'i':'j'), a1, a2, (nface<=4?'k':'k'), a3, a4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           Pt/Pinf        = %8.4f", pte);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           Tt/Tinf        = %8.4f", tte);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           alphe          = %8.4f", alpe);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           betae          = %8.4f", betae);
            if (itrflg1 > 0) {
                t1 = (float)bcdata(1,1,1,5);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb1 (nondim) = %12.4E", t1);
            }
            if (itrflg2 > 0) {
                t2 = (float)bcdata(1,1,1,6);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb2 (nondim) = %12.4E", t2);
            }
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d",
                hdr, (nface<=4?'i':'j'), a1, a2, (nface<=4?'k':'k'), a3, a4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           inflow data set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    const char** hdrs = (iflgg == 0) ? hdr09 : hdr10;
    if (nface == 3) { do_check(); do_print(hdrs[0], ista,iend,ksta,kend); }
    if (nface == 4) { do_check(); do_print(hdrs[1], ista,iend,ksta,kend); }
    if (nface == 5) { do_check(); do_print(hdrs[2], ista,iend,jsta,jend); }
    if (nface == 6) { do_check(); do_print(hdrs[3], ista,iend,jsta,jend); }
    if (nface == 1) { do_check(); do_print(hdrs[4], jsta,jend,ksta,kend); }
    if (nface == 2) { do_check(); do_print(hdrs[5], jsta,jend,ksta,kend); }
}


// ============================================================
// out2019 - OVERFLOW-type nozzle total BC type 2019
// ============================================================
void out2019(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int& ivmx = cmn_maxiv.ivmx;
    int iflg, itrflg1, itrflg2;
    float pte, tte, t1, t2;
    char filname_str[61];

    static const char* hdrs[6] = {
        "   j=   1  OVERFLOW-type nozzle total BC  type 2019",
        "   j=jdim  OVERFLOW-type nozzle total BC  type 2019",
        "   k=   1  OVERFLOW-type nozzle total BC  type 2019",
        "   k=kdim  OVERFLOW-type nozzle total BC  type 2019",
        "   i=   1  OVERFLOW-type nozzle total BC  type 2019",
        "   i=idim  OVERFLOW-type nozzle total BC  type 2019"
    };

    auto do_check = [&]() {
        itrflg1 = 0; itrflg2 = 0;
        if ((float)bcdata(1,1,1,3) >= 0.f) itrflg1 = 1;
        if ((float)bcdata(1,1,1,4) >= 0.f) itrflg2 = 1;
        iflg = 0;
        for (int l=1; l<=2; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2019:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  nozzle total BC data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (itrflg1 > 0 || itrflg2 > 0) {
            if (ivmx > 5 && itrflg1*itrflg2 == 0) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  stopping in bc2019...must set both turbulence");
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "  quantities for 2-eq turb. models...only one is set");
                int m1 = -1;
                termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
            }
        }
    };

    auto do_print = [&](const char* hdr, int a1, int a2, int a3, int a4, char c1, char c2) {
        if (strncmp(filname, "null", 4) == 0) {
            pte = (float)bcdata(1,1,1,1);
            tte = (float)bcdata(1,1,1,2);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d", hdr, c1, a1, a2, c2, a3, a4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "            Pt/Ptinf        = %8.4f", pte);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "            Tt/Ttinf        = %8.4f", tte);
            if (itrflg1 > 0) {
                t1 = (float)bcdata(1,1,1,3);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb1 (nondim) = %12.4E", t1);
            }
            if (itrflg2 > 0) {
                t2 = (float)bcdata(1,1,1,4);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb2 (nondim) = %12.4E", t2);
            }
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d", hdr, c1, a1, a2, c2, a3, a4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           inflow data set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    if (nface == 3) { do_check(); do_print(hdrs[0], ista,iend,ksta,kend,'i','k'); }
    if (nface == 4) { do_check(); do_print(hdrs[1], ista,iend,ksta,kend,'i','k'); }
    if (nface == 5) { do_check(); do_print(hdrs[2], ista,iend,jsta,jend,'i','j'); }
    if (nface == 6) { do_check(); do_print(hdrs[3], ista,iend,jsta,jend,'i','j'); }
    if (nface == 1) { do_check(); do_print(hdrs[4], jsta,jend,ksta,kend,'j','k'); }
    if (nface == 2) { do_check(); do_print(hdrs[5], jsta,jend,ksta,kend,'j','k'); }
}


// ============================================================
// out2016 - viscous wall type 2016 (with synthetic jet)
// ============================================================
void out2016(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iout)
{
    int* iwf = cmn_wallfun.iwf;
    float dt = cmn_info.dt;

    int iflg;
    float ctemp, cq, cqu, sjetx, sjety, sjetz, rfreq;
    char filname_str[61];

    static const char* hdrs[6] = {
        "   j=   1  viscous wall                   type 2016",
        "   j=jdim  viscous wall                   type 2016",
        "   k=   1  viscous wall                   type 2016",
        "   k=kdim  viscous wall                   type 2016",
        "   i=   1  viscous wall                   type 2016",
        "   i=idim  viscous wall                   type 2016"
    };

    auto do_check = [&]() {
        iflg = 0;
        for (int l=1; l<=2+iout; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++) {
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
                        if (bcdata(d1,d2,ip,l) != 0.0 && (float)dt < 0.f && l == 3) iflg = 2;
                    }
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2016:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  viscous wall data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        if (iflg == 2) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2016:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  unsteady jet not allowed with dt < 0");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    };

    auto do_print = [&](const char* hdr, int a1, int a2, int a3, int a4) {
        if (strncmp(filname, "null", 4) == 0) {
            ctemp = (float)bcdata(1,1,1,1);
            cq    = (float)bcdata(1,1,1,2);
            cqu   = (float)bcdata(1,1,1,3);
            sjetx = (float)bcdata(1,1,1,4);
            sjety = (float)bcdata(1,1,1,5);
            sjetz = (float)bcdata(1,1,1,6);
            rfreq = (float)bcdata(1,1,1,7);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d",
                hdr, (nface<=4?'i':'j'), a1, a2, (nface<=4?'k':'k'), a3, a4);
            if (ctemp > 0.f) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           Tw/Tinf    = %8.4f", ctemp);
            } else if (ctemp == 0.f) {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           Tw/Tinf    = adiabatic wall");
            } else {
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           Tw/Tinf    = stagnation");
            }
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           cq         = %8.4f", cq);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           cqu        = %8.4f", cqu);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           sjetx      = %8.4f", sjetx);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           sjetx      = %8.4f", sjety);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           sjetx      = %8.4f", sjetz);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           rfreq      = %8.4f", rfreq);
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d",
                hdr, (nface<=4?'i':'j'), a1, a2, (nface<=4?'k':'k'), a3, a4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           Twtype, C_q set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    if (nface == 3) {
        do_check();
        do_print(hdrs[0], ista,iend,ksta,kend);
        if (iwf[2-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"          Wall function employed in j-direction"); }
    }
    if (nface == 4) {
        do_check();
        do_print(hdrs[1], ista,iend,ksta,kend);
        if (iwf[2-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"          Wall function employed in j-direction"); }
    }
    if (nface == 5) {
        do_check();
        do_print(hdrs[2], ista,iend,jsta,jend);
        if (iwf[3-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"          Wall function employed in k-direction"); }
    }
    if (nface == 6) {
        do_check();
        do_print(hdrs[3], ista,iend,jsta,jend);
        if (iwf[3-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"          Wall function employed in k-direction"); }
    }
    if (nface == 1) {
        do_check();
        do_print(hdrs[4], jsta,jend,ksta,kend);
        if (iwf[1-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"          Wall function employed in i-direction"); }
    }
    if (nface == 2) {
        do_check();
        do_print(hdrs[5], jsta,jend,ksta,kend);
        if (iwf[1-1] == 1) { nou(1)=std::min(nou(1)+1,ibufdim); std::snprintf(bou(nou(1),1),120,"          Wall function employed in i-direction"); }
    }
}


// ============================================================
// out2026 - sweeping jet type 2026
// ============================================================
void out2026(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl)
{
    int iflg, itrflg1, itrflg2;
    float vmag, rfreq, sideangj, sxa, sya, sza, sxb, syb, szb, t1, t2;
    char filname_str[61];

    static const char* hdrs[6] = {
        "   j=   1  sweeping jet   type 2026",
        "   j=jdim  sweeping jet   type 2026",
        "   k=   1  sweeping jet   type 2026",
        "   k=kdim  sweeping jet   type 2026",
        "   i=   1  sweeping jet   type 2026",
        "   i=idim  sweeping jet   type 2026"
    };

    auto do_check = [&]() {
        itrflg1 = 0; itrflg2 = 0;
        if ((float)bcdata(1,1,1,10) >= 0.f) itrflg1 = 1;
        if ((float)bcdata(1,1,1,11) >= 0.f) itrflg2 = 1;
        iflg = 0;
        for (int l=1; l<=9; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2026:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  variable data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
        iflg = 0;
        for (int ip=1; ip<=2; ip++)
            for (int d1=1; d1<=mdim; d1++)
                for (int d2=1; d2<=ndim; d2++) {
                    if ((float)bcdata(d1,d2,ip,4) == 0.f && (float)bcdata(d1,d2,ip,5) == 0.f && (float)bcdata(d1,d2,ip,6) == 0.f) iflg = 1;
                    if ((float)bcdata(d1,d2,ip,7) == 0.f && (float)bcdata(d1,d2,ip,8) == 0.f && (float)bcdata(d1,d2,ip,9) == 0.f) iflg = 1;
                }
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2026:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  normal vectors cannot be zero");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    };

    auto do_print = [&](const char* hdr, int a1, int a2, int a3, int a4, char c1, char c2) {
        if (strncmp(filname, "null", 4) == 0) {
            vmag     = (float)bcdata(1,1,1,1);
            rfreq    = (float)bcdata(1,1,1,2);
            sideangj = (float)bcdata(1,1,1,3);
            sxa      = (float)bcdata(1,1,1,4);
            sya      = (float)bcdata(1,1,1,5);
            sza      = (float)bcdata(1,1,1,6);
            sxb      = (float)bcdata(1,1,1,7);
            syb      = (float)bcdata(1,1,1,8);
            szb      = (float)bcdata(1,1,1,9);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d", hdr, c1, a1, a2, c2, a3, a4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           vmag/ainf      = %8.4f", vmag);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           rfreq*lref/ainf= %8.4f", rfreq);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           sideangj (deg) = %8.4f", sideangj);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           sxa,sya,sza    = %8.4f%8.4f%8.4f", sxa, sya, sza);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           sxb,syb,szb    = %8.4f%8.4f%8.4f", sxb, syb, szb);
            if (itrflg1 > 0) {
                t1 = (float)bcdata(1,1,1,10);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb1 (nondim) = %12.4E", t1);
            }
            if (itrflg2 > 0) {
                t2 = (float)bcdata(1,1,1,11);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           turb2 (nondim) = %12.4E", t2);
            }
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d", hdr, c1, a1, a2, c2, a3, a4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           inflow data set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    if (nface == 3) { do_check(); do_print(hdrs[0], ista,iend,ksta,kend,'i','k'); }
    if (nface == 4) { do_check(); do_print(hdrs[1], ista,iend,ksta,kend,'i','k'); }
    if (nface == 5) { do_check(); do_print(hdrs[2], ista,iend,jsta,jend,'i','j'); }
    if (nface == 6) { do_check(); do_print(hdrs[3], ista,iend,jsta,jend,'i','j'); }
    if (nface == 1) { do_check(); do_print(hdrs[4], jsta,jend,ksta,kend,'j','k'); }
    if (nface == 2) { do_check(); do_print(hdrs[5], jsta,jend,ksta,kend,'j','k'); }
}


// ============================================================
// out2102 - time-varying pressure type 2102/2103
// ============================================================
void out2102(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim, int& mdim, int& ndim,
             FortranArray4DRef<double> bcdata, char (&filname)[80],
             int& myid, FortranArray1DRef<int> mblk2nd, int& maxbl, int& iflag)
{
    float twopi = (float)(8.0 * std::atan(1.0));
    int l_set = (iflag == 0) ? 4 : 5;
    int iflg;
    float pratio, deltap, rfreqp, lref, phioff;
    char filname_str[61];

    // type 2102 headers (iflag==0) and type 2103 headers (iflag==1)
    static const char* hdr02_null[6] = {
        "   j=   1  time-varying pressure          type 2102",
        "   j=jdim  time-varying pressure          type 2102",
        "   k=   1  time-varying pressure          type 2102",
        "   k=kdim  time-varying pressure          type 2102",
        "   i=   1  time-varying pressure          type 2102",
        "   i=idim  time-varying pressure          type 2102"
    };
    static const char* hdr03_null[6] = {
        "   j=   1  time-varying pressure          type 2103",
        "   j=jdim  time-varying pressure          type 2103",
        "   k=   1  time-varying pressure          type 2103",
        "   k=kdim  time-varying pressure          type 2103",
        "   i=   1  time-varying pressure          type 2103",
        "   i=idim  time-varying pressure          type 2103"
    };

    auto do_check = [&]() {
        iflg = 0;
        for (int l=1; l<=l_set; l++)
            for (int ip=1; ip<=2; ip++)
                for (int d1=1; d1<=mdim; d1++)
                    for (int d2=1; d2<=ndim; d2++)
                        if ((float)bcdata(d1,d2,ip,l) < -1.e10f) iflg = 1;
        if (iflg == 1) {
            nou(1) = std::min(nou(1)+1, ibufdim);
            if (iflag == 0) std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2102:");
            else            std::snprintf(bou(nou(1),1), 120, "  Stopping in bc2103:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "  pressure data incorrectly set ");
            int m1 = -1;
            termn8_ns::termn8(myid, m1, ibufdim, nbuf, bou, nou);
        }
    };

    auto do_print = [&](const char* hdr_null, int a1, int a2, int a3, int a4, char c1, char c2) {
        if (strncmp(filname, "null", 4) == 0) {
            pratio = (float)bcdata(1,1,1,1);
            deltap = (float)bcdata(1,1,1,2);
            rfreqp = (float)bcdata(1,1,1,3);
            lref   = (float)bcdata(1,1,1,4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s  %c=%5d,%5d  %c=%5d,%5d", hdr_null, c1, a1, a2, c2, a3, a4);
            nou(1) = std::min(nou(1)+1, ibufdim);
            if (iflag == 0)
                std::snprintf(bou(nou(1),1), 120, "           P(t)/Pinf  = P1 + eps*sin(k*t)");
            else
                std::snprintf(bou(nou(1),1), 120, "           P(t)/Pinf  = P1 + eps*sin(k*t+phioff)");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           P1         = %8.4f", pratio);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           eps        = %8.4f", deltap);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           k          = %8.4f", twopi*rfreqp/lref);
            if (iflag == 1) {
                phioff = (float)bcdata(1,1,1,5);
                nou(1) = std::min(nou(1)+1, ibufdim);
                std::snprintf(bou(nou(1),1), 120, "           phioff(deg)= %8.4f", phioff);
            }
        } else {
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "%s", hdr_null);
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(bou(nou(1),1), 120, "           P(t)/Pinf set from file:");
            nou(1) = std::min(nou(1)+1, ibufdim);
            std::snprintf(filname_str, 61, "%.60s", filname);
            std::snprintf(bou(nou(1),1), 120, "           %s", filname_str);
        }
    };

    const char** hdrs = (iflag == 0) ? hdr02_null : hdr03_null;
    if (nface == 3) { do_check(); do_print(hdrs[0], ista,iend,ksta,kend,'i','k'); }
    if (nface == 4) { do_check(); do_print(hdrs[1], ista,iend,ksta,kend,'i','k'); }
    if (nface == 5) { do_check(); do_print(hdrs[2], ista,iend,jsta,jend,'i','j'); }
    if (nface == 6) { do_check(); do_print(hdrs[3], ista,iend,jsta,jend,'i','j'); }
    if (nface == 1) { do_check(); do_print(hdrs[4], jsta,jend,ksta,kend,'j','k'); }
    if (nface == 2) { do_check(); do_print(hdrs[5], jsta,jend,ksta,kend,'j','k'); }
}

// ============================================================
// out9999 - exact soln (MMS) type 9999
// ============================================================
void out9999(int& jdim, int& kdim, int& idim, int& ista, int& iend,
             int& jsta, int& jend, int& ksta, int& kend, int& nface,
             FortranArray1DRef<int> nou, FortranArray2DRef<char[120]> bou,
             int& nbuf, int& ibufdim)
{
    if (nface == 3) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=   1  exact soln (MMS)               type 9999  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 4) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   j=jdim  exact soln (MMS)               type 9999  i=%5d,%5d  k=%5d,%5d", ista,iend,ksta,kend);
    }
    if (nface == 5) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=   1  exact soln (MMS)               type 9999  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 6) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   k=kdim  exact soln (MMS)               type 9999  i=%5d,%5d  j=%5d,%5d", ista,iend,jsta,jend);
    }
    if (nface == 1) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=   1  exact soln (MMS)               type 9999  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
    if (nface == 2) {
        nou(1) = std::min(nou(1)+1, ibufdim);
        std::snprintf(bou(nou(1),1), 120, "   i=idim  exact soln (MMS)               type 9999  j=%5d,%5d  k=%5d,%5d", jsta,jend,ksta,kend);
    }
}

} // namespace bc_info_ns
