// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// clusterE.cpp — maggie "Cluster E: diagnostics + output + stencil checks"
// Faithful F77 -> C++ port. Byte-exact I/O (unit 2 binary, unit 6/8 text).
// Built with: g++ -std=c++17 -O2 -ffp-contract=off -fdefault-real-8 semantics.
//
// SIGNATURE: void diagnos();
// SIGNATURE: void output();
// SIGNATURE: void chkxyz();
// SIGNATURE: void chkstn(int* iblank);
// SIGNATURE: void ibchk(int* ib,int* iblank,int ic,int is,int* ji,int* ki,int* li,int m1,int jd,int kd,int ld,
//                       double* xb,double* yb,double* zb,int* jinew,int* kinew,int* linew,double* x,double* y,double* z);
//
// ---------------------------------------------------------------------------
// mag_protos.h vs Fortran MISMATCHES (matched proto as instructed, flagged here):
//
//  * cindex : proto  void cindex(int* ibc,int ibpnts,int iipnts,int* jb,int* kb,int* lb,int mesh)
//             Fortran: call cindex( ibc, ibpnts(m), iipnts(m), jb,kb,lb, m )
//             In Fortran ibpnts(m) and iipnts(m) are OUTPUT scalars (set inside
//             cindex, then read back by output & written to unit 2). The proto
//             passes them BY VALUE (plain int), so the outputs CANNOT flow back.
//             To compile against the proto I call cindex(ibc, ibpnts[m], iipnts[m],
//             ...) passing them by value — but ibpnts[m]/iipnts[m] then stay
//             whatever cindex could not set. This is almost certainly a broken
//             proto.  >>> RECONCILE REQUIRED: proto should be
//                 void cindex(int* ibc,int& ibpnts,int& iipnts,...)  <<<
//             With the by-value proto the unit-2 record 2 fields ibpnts/iipnts
//             are UNRELIABLE. Flagging strongly; parent must fix the proto+def.
//
//  * reorder: proto  void reorder(int idimA,int ibpnts,int intpts,int* jb,int* kb,int* lb,int* iord,int jd,int kd,int ld)
//             Fortran: dimension intpts(4); call reorder(idim,ibpnts(m),intpts,...)
//             `intpts` is an ARRAY(4) filled by reorder and then written to unit 2
//             (4 ints in record 2). The proto types it as a by-value `int`, so it
//             cannot return the 4 values.  >>> RECONCILE REQUIRED: proto should be
//                 void reorder(int idimA,int ibpnts,int* intpts,...)  <<<
//             To compile against the current proto I pass intpts[1] by value; the
//             filled array is therefore NOT populated. Unit-2 record 2 intpts(1..4)
//             is UNRELIABLE until proto+def take `int* intpts`. Flagging strongly.
//
//  * bound  : proto  void bound(...,double* x,double* y,double* z,int ife)
//             Fortran: dimension ife(6); call bound(...,x,y,z,ife)  (ife all 0)
//             `ife` is an ARRAY(6) in Fortran; proto types it by-value int. Here
//             ife is always all-zero on entry and (in these callers) its returned
//             contents are unused, so passing ife[1] (==0) by value is harmless
//             for THIS cluster.  >>> RECONCILE (proto should be `int* ife`). <<<
//
//  * pltpts : Fortran name `pltpts`; proto is `pltpts_`. Called as pltpts_ below.
// ---------------------------------------------------------------------------

#include "mag_common.h"

// ===========================================================================
//  subroutine diagnos
// ===========================================================================
void diagnos()
{
    std::fprintf(funit(6), " \n");
    std::fprintf(funit(6),
        "   *** beginning diagnostic checks of interpolation/extrapolation stencils ***\n");

    // check for hole points or interpolated points as part of interp stencil
    chkstn(book3.iblank);

    // check for cell-center coordinate mismatch
    chkxyz();

    // output plot3d diagnostic file
    if (diagno.iplt3d > 0) {
        std::fprintf(funit(6), " \n");
        std::fprintf(funit(6),
            "  plot3d diagnostic file (mg/blank) being written to %-60s\n",
            files.plt3d);
        std::fprintf(funit(6), "         the file has the following structure:\n");

        int nmesh3 = conec1.nmesh * 3;
        int mm = 0;
        for (int m = 1; m <= nmesh3; m += 3) {
            mm = mm + 1;
            std::fprintf(funit(6), " \n");
            // 101 format(8x,'plot3d grid ',i2,'.....mesh ',i2,' field points')
            std::fprintf(funit(6),
                "        plot3d grid %2d.....mesh %2d field points\n", m, mm);
            // 102
            std::fprintf(funit(6),
                "        plot3d grid %2d.....mesh %2d hole points\n", m + 1, mm);
            // 103
            std::fprintf(funit(6),
                "        plot3d grid %2d.....mesh %2d fringe points\n", m + 2, mm);
        }
        if (pltpt.nbnd > 0) {
            // 104 format(/,8x,'plot3d grid ',i2,'.....interpolated boundary points (all meshes)')
            std::fprintf(funit(6), "\n");
            std::fprintf(funit(6),
                "        plot3d grid %2d.....interpolated boundary points (all meshes)\n",
                nmesh3 + 1);
            // 1012 format(22x,'.....must plot as POINTS in plot3d')
            std::fprintf(funit(6),
                "                      .....must plot as POINTS in plot3d\n");
        }
        if (pltpt.norph > 0) {
            // 105 format(/,8x,...'.....extrapolated points (orphans) (all meshes)')
            std::fprintf(funit(6), "\n");
            if (pltpt.nbnd > 0)
                std::fprintf(funit(6),
                    "        plot3d grid %2d.....extrapolated points (orphans) (all meshes)\n",
                    nmesh3 + 2);
            else
                std::fprintf(funit(6),
                    "        plot3d grid %2d.....extrapolated points (orphans) (all meshes)\n",
                    nmesh3 + 1);
            // 1012
            std::fprintf(funit(6),
                "                      .....must plot as POINTS in plot3d\n");
        }
        if (pltpt.nill > 0) {
            // 106 format(/,8x,...'.....points with illegal stencils (all meshes)')
            std::fprintf(funit(6), "\n");
            if (pltpt.nbnd > 0 && pltpt.norph > 0)
                std::fprintf(funit(6),
                    "        plot3d grid %2d.....points with illegal stencils (all meshes)\n",
                    nmesh3 + 3);
            else if (pltpt.nbnd > 0 || pltpt.norph > 0)
                std::fprintf(funit(6),
                    "        plot3d grid %2d.....points with illegal stencils (all meshes)\n",
                    nmesh3 + 2);
            else
                std::fprintf(funit(6),
                    "        plot3d grid %2d.....points with illegal stencils (all meshes)\n",
                    nmesh3 + 1);
            // 1013 format(22x,'.....(does not reflect any illegal stencils',
            //             /,22x,'.....corrected by zeroth order data transfer)')
            std::fprintf(funit(6),
                "                      .....(does not reflect any illegal stencils\n");
            std::fprintf(funit(6),
                "                      .....corrected by zeroth order data transfer)\n");
            // 1012
            std::fprintf(funit(6),
                "                      .....must plot as POINTS in plot3d\n");
        }

        for (int m = 1; m <= conec1.nmesh; ++m) {
            int jd = grid2.mjmax[m];
            int kd = grid2.mkmax[m];
            int ld = grid2.mlmax[m];
            getgrd(m, grid1.x, grid1.y, grid1.z,
                   image.jimage, image.kimage, image.limage, jd, kd, ld);
            getibl(m, book3.iblank, jd, kd, ld);
            pltpts_(conec1.nmesh, m, jd, kd, ld, grid1.x, grid1.y, grid1.z,
                    image.jimage, book3.iblank);
        }
    }
}

// ===========================================================================
//  subroutine output
//     writes summary to unit 6 and interpolation/connection data to unit 2
// ===========================================================================
void output()
{
    // local arrays (Fortran: dimension ibc(idim),ibpnts(mdim),iipnts(mdim))
    //                        dimension intpts(4),iord(idim)
    static int    ibc[idim + 1];
    static int    ibpnts[mdim + 1];
    static int    iipnts[mdim + 1];
    int           intpts[4 + 1];        // 1-based intpts(4)
    static int    iord[idim + 1];

    // --- count hole/fringe points ---
    int iholpt = 0;
    int ifrgpt = 0;
    int npnts  = 0;

    for (int m = 1; m <= conec1.nmesh; ++m) {
        int jd = grid2.mjmax[m];
        int kd = grid2.mkmax[m];
        int ld = grid2.mlmax[m];
        npnts = npnts + jd * kd * ld;
        getibl(m, book3.iblank, jd, kd, ld);

        for (int l = 1; l <= ld; ++l)
        for (int k = 1; k <= kd; ++k)
        for (int j = 1; j <= jd; ++j) {
            int i = j + (k - 1) * jd + (l - 1) * jd * kd;
            if (book3.iblank[i] == 0) iholpt = iholpt + 1;
            if (book3.iblank[i] <  0) ifrgpt = ifrgpt + 1;
        }
    }

    // list-directed write(6,*): leading blank, then fields.
    std::fprintf(funit(6), " \n");
    std::fprintf(funit(6), "   *** summary of maggie preprocessor ***\n");
    std::fprintf(funit(6), " \n");
    // write(6,*) '    there are ', nmesh,' grids in the composite mesh'
    std::fprintf(funit(6), "     there are %12d grids in the composite mesh\n",
                 conec1.nmesh);
    // write(6,*) '    with a total of ',npnts,' points'
    std::fprintf(funit(6), "     with a total of %12d points\n", npnts);
    // write(6,*) '    of which: ',iholpt,' are hole points'
    std::fprintf(funit(6), "     of which: %12d are hole points\n", iholpt);
    // write(6,*) '              ',ifrgpt,' are fringe points'
    std::fprintf(funit(6), "               %12d are fringe points\n", ifrgpt);

    // --- set iblank for flow solver (remove connection info) ---
    for (int m = 1; m <= conec1.nmesh; ++m) {
        int jd = grid2.mjmax[m];
        int kd = grid2.mkmax[m];
        int ld = grid2.mlmax[m];
        getibl(m, book3.iblank, jd, kd, ld);
        // output iblank WITH connection info for checkmag2
        putibl2(m, book3.iblank, jd, kd, ld);
        iblcon(book3.iblank, jd, kd, ld);
        putibl(m, book3.iblank, jd, kd, ld);
    }

    int iitot = 0;
    for (int m = 1; m <= conec1.nmesh; ++m)
        iitot = iitot + book2.ibpts[m];

    // --- set up cross index array + pointers, write unit 2 ---
    int iisptr = 0;   // (declared to persist across meshes as in Fortran)
    int iieptr = 0;
    for (int m = 1; m <= conec1.nmesh; ++m) {

        std::fprintf(funit(6), " \n");
        // write(6,*) '      writing connection data for CFL3D to ',' unit 2 for mesh ',m
        std::fprintf(funit(6),
            "       writing connection data for CFL3D to  unit 2 for mesh %12d\n",
            m);

        if (m == 1) {
            iisptr = 1;
            iieptr = book2.ibpts[m];
        } else {
            iisptr = iisptr + book2.ibpts[m - 1];
            iieptr = iieptr + book2.ibpts[m];
        }

        // call cindex( ibc,ibpnts(m),iipnts(m),jb,kb,lb,m )
        // proto passes ibpnts/iipnts BY VALUE (FLAGGED at top): outputs cannot
        // return through this call. Matching the proto verbatim.
        cindex(ibc, ibpnts[m], iipnts[m],
               intrp2.jb, intrp2.kb, intrp2.lb, m);

        // write(6,*) '        ibpnts,iipnts,iieptr,iisptr ',ibpnts(m),iipnts(m),iieptr,iisptr
        std::fprintf(funit(6),
            "         ibpnts,iipnts,iieptr,iisptr %12d%12d%12d%12d\n",
            ibpnts[m], iipnts[m], iieptr, iisptr);

        // call reorder(idim,ibpnts(m),intpts,jb,kb,lb,iord,mjmax(m),mkmax(m),mlmax(m))
        // intpts is ARRAY(4) in Fortran; proto is by-value int (FLAGGED). Matching
        // proto verbatim by passing intpts[1]; the array is NOT filled by this call.
        reorder(idim, ibpnts[m], intpts,
                intrp2.jb, intrp2.kb, intrp2.lb, iord,
                grid2.mjmax[m], grid2.mkmax[m], grid2.mlmax[m]);

        // write(6,*) '        intpts = ',intpts   (list-directed: 4 ints)
        std::fprintf(funit(6), "         intpts = %12d%12d%12d%12d\n",
                     intpts[1], intpts[2], intpts[3], intpts[4]);

        // -------------------------------------------------------------------
        //  UNIT 2 RECORDS (ovrlp.bin), unformatted sequential.
        //  -fdefault-real-8: REAL -> 8-byte double; INTEGER -> 4-byte int.
        //  Each Fortran write(2) is ONE record.
        // -------------------------------------------------------------------

        // write(2) mjmax(m),mkmax(m),mlmax(m)          [record 1: 3 ints]
        {
            tc::RecordWriter rw;
            rw.put_i(grid2.mjmax[m]);
            rw.put_i(grid2.mkmax[m]);
            rw.put_i(grid2.mlmax[m]);
            rw.flush(funit(2));
        }

        // write(2) ibpnts(m),intpts,iipnts(m),iieptr,iisptr
        //   intpts is intpts(4) -> expands to 4 ints.
        //   [record 2: ibpnts, intpts(1..4), iipnts, iieptr, iisptr = 8 ints]
        {
            tc::RecordWriter rw;
            rw.put_i(ibpnts[m]);
            rw.put_i(intpts[1]);
            rw.put_i(intpts[2]);
            rw.put_i(intpts[3]);
            rw.put_i(intpts[4]);
            rw.put_i(iipnts[m]);
            rw.put_i(iieptr);
            rw.put_i(iisptr);
            rw.flush(funit(2));
        }

        // call getint( m,ji,ki,li,jbpt,kbpt,lbpt,dxint,dyint,dzint )
        getint(m, book2.ji, book2.ki, book2.li,
               book2.jbpt, book2.kbpt, book2.lbpt,
               book2.dxint, book2.dyint, book2.dzint);

        // write(2) ( ji(i),ki(i),li(i),dxint(i),dyint(i),dzint(i), i=1,iipnts(m) )
        //   [record 3: per i -> 3 ints + 3 doubles, contiguous within one record]
        {
            tc::RecordWriter rw;
            for (int i = 1; i <= iipnts[m]; ++i) {
                rw.put_i(book2.ji[i]);
                rw.put_i(book2.ki[i]);
                rw.put_i(book2.li[i]);
                rw.put_d(book2.dxint[i]);
                rw.put_d(book2.dyint[i]);
                rw.put_d(book2.dzint[i]);
            }
            rw.flush(funit(2));
        }

        // write(2) ( jb(iord(i)),kb(iord(i)),lb(iord(i)),ibc(iord(i)), i=1,ibpnts(m) )
        //   [record 4: per i -> 4 ints, contiguous within one record]
        {
            tc::RecordWriter rw;
            for (int i = 1; i <= ibpnts[m]; ++i) {
                int io = iord[i];
                rw.put_i(intrp2.jb[io]);
                rw.put_i(intrp2.kb[io]);
                rw.put_i(intrp2.lb[io]);
                rw.put_i(ibc[io]);
            }
            rw.flush(funit(2));
        }

        int jd = grid2.mjmax[m];
        int kd = grid2.mkmax[m];
        int ld = grid2.mlmax[m];

        // add iblank array to file 2, as required by CFL3D
        getibl(m, book3.iblank, jd, kd, ld);
        wiblnk(book3.iblank, jd, kd, ld);   // wiblnk itself writes record(s) to unit 2
    }

    // --- minimum parameter sizes ---
    int iidmax = 0;
    int ibdmax = 0;
    for (int m = 1; m <= conec1.nmesh; ++m) {
        iidmax = iidmax + iipnts[m];
        ibdmax = ibdmax + ibpnts[m];
    }
    int iitotl = std::max(iidmax, ibdmax);

    std::fprintf(funit(6), " \n");
    // write(6,*) '    minimum dimension for parameter IITOT',' in CFL3D (module RHS): ',iitotl
    std::fprintf(funit(6),
        "     minimum dimension for parameter IITOT in CFL3D (module RHS): %12d\n",
        iitotl);

    (void)iitot;   // computed as in Fortran; not otherwise used

    // --- free up scratch temp_* files (open then delete) ---
    char titl[21];
    for (int pass = 0; pass < 4; ++pass) {
        const char* prefix = (pass == 0) ? "temp_grd." :
                             (pass == 1) ? "temp_cen." :
                             (pass == 2) ? "temp_ibl." : "temp_int.";
        for (int m = 1; m <= conec1.nmesh; ++m) {
            int len;
            if (m > 99) {
                len = 12;
                std::snprintf(titl, sizeof(titl), "%s%3d", prefix, m);
            } else if (m > 9) {
                len = 11;
                std::snprintf(titl, sizeof(titl), "%s%2d", prefix, m);
            } else {
                len = 10;
                std::snprintf(titl, sizeof(titl), "%s%1d", prefix, m);
            }
            for (int i = len; i < 20; ++i) titl[i] = ' ';
            titl[len] = '\0';   // filename = titl(1:len)
            int iunit = 30;
            fopen_unit(iunit, titl, "wb+");   // open status='unknown'
            fclose_unit(iunit);               // close status='delete'
            std::remove(titl);
        }
    }
}

// ===========================================================================
//  subroutine chkxyz
// ===========================================================================
void chkxyz()
{
    View2<int> IPNTR(book1.ipntr, mdim, mhldim * mdim);
    View2<int> NPNTR(book1.npntr, mdim, mhldim * mdim);
    View2<int> MHBS (book1.mhbs,  mdim, mdim);
    View2<int> MOBS (book1.mobs,  mdim, mdim);

    int ife[6 + 1];
    for (int ll = 1; ll <= 6; ++ll) ife[ll] = 0;

    std::fprintf(funit(6), " \n");
    // write(6,*) '    checking for mismatch between cell',' centers of fringe and boundary pts'
    std::fprintf(funit(6),
        "     checking for mismatch between cell centers of fringe and boundary pts\n");
    // write(6,*) '    as calculated from 1) input meshes and 2)',' data transfer stencils'
    std::fprintf(funit(6),
        "     as calculated from 1) input meshes and 2) data transfer stencils\n");

    double sigdif = .1;

    for (int m = 1; m <= conec1.nmesh; ++m) {

        int itotm = 0;
        int itotn = 0;
        double errmax = 0.;
        int jerrmx = 1;
        int kerrmx = 1;
        int lerrmx = 1;
        int nsig = 0;
        int nh = conec1.nhole[m];
        int no = conec1.noutr[m];
        int jmax = grid2.mjmax[m];
        int kmax = grid2.mkmax[m];
        int lmax = grid2.mlmax[m];

        int jinc = 1, kinc = 1, linc = 1;
        if (jmax < 2) jinc = 0;
        if (kmax < 2) kinc = 0;
        if (lmax < 2) linc = 0;

        // --- check fringe points ---
        if (nh > 0) {
            for (int n = 1; n <= nh; ++n) {
                int nserch = conec1.MHOLE(m, n);
                for (int nn = 1; nn <= nserch; ++nn) {
                    int m1 = conec1.LHOLE(m, n, nn);
                    int jmax1 = grid2.mjmax[m1];
                    int kmax1 = grid2.mkmax[m1];
                    int lmax1 = grid2.mlmax[m1];
                    getint(m1, book2.ji, book2.ki, book2.li,
                           book2.jbpt, book2.kbpt, book2.lbpt,
                           book2.dxint, book2.dyint, book2.dzint);
                    int iset = MHBS(m, m1);
                    int is = IPNTR(m1, iset);
                    int ie = NPNTR(m1, iset);
                    getgrd(m, grid1.x, grid1.y, grid1.z,
                           image.jimage, image.kimage, image.limage,
                           jmax, kmax, lmax);

                    for (int i = is; i <= ie; ++i) {
                        double x0, y0, z0;
                        // proto bound: last arg int ife (FLAGGED: Fortran ife(6) array)
                        bound(book2.jbpt[i], book2.kbpt[i], book2.lbpt[i],
                              x0, y0, z0, jmax, kmax, lmax,
                              grid1.x, grid1.y, grid1.z, ife);

                        int j00 = book2.jbpt[i];
                        int k00 = book2.kbpt[i];
                        int l00 = book2.lbpt[i];
                        if (j00 >= jmax) j00 = jmax - 1;
                        if (k00 >= kmax) k00 = kmax - 1;
                        if (l00 >= lmax) l00 = lmax - 1;
                        if (j00 <= 0) j00 = 1;
                        if (k00 <= 0) k00 = 1;
                        if (l00 <= 0) l00 = 1;

                        double xscal = 0., yscal = 0., zscal = 0.;
                        for (int j = j00; j <= j00 + jinc; ++j)
                        for (int k = k00; k <= k00 + kinc; ++k)
                        for (int l = l00; l <= l00 + linc; ++l) {
                            int ii = j + (k - 1) * jmax + (l - 1) * jmax * kmax;
                            double dx = std::abs(grid1.x[ii] - x0);
                            double dy = std::abs(grid1.y[ii] - y0);
                            double dz = std::abs(grid1.z[ii] - z0);
                            if (dx > xscal) xscal = dx;
                            if (dy > yscal) yscal = dy;
                            if (dz > zscal) zscal = dz;
                        }
                        if (xscal <= 0.) xscal = 1.;
                        if (yscal <= 0.) yscal = 1.;
                        if (zscal <= 0.) zscal = 1.;

                        itotn = itotn + 1;
                        chkpt.xcc0[itotn]   = x0;
                        chkpt.ycc0[itotn]   = y0;
                        chkpt.zcc0[itotn]   = z0;
                        chkpt.xscale[itotn] = xscal;
                        chkpt.yscale[itotn] = yscal;
                        chkpt.zscale[itotn] = zscal;
                    } // 141

                    getgrd(m1, grid1.x, grid1.y, grid1.z,
                           image.jimage, image.kimage, image.limage,
                           jmax1, kmax1, lmax1);

                    for (int i = is; i <= ie; ++i) {
                        double xintrp, yintrp, zintrp;
                        interp(jmax1, kmax1, lmax1, grid1.x, grid1.y, grid1.z,
                               book2.ji[i], book2.ki[i], book2.li[i],
                               book2.dxint[i], book2.dyint[i], book2.dzint[i],
                               xintrp, yintrp, zintrp);

                        itotm = itotm + 1;
                        double x0    = chkpt.xcc0[itotm];
                        double y0    = chkpt.ycc0[itotm];
                        double z0    = chkpt.zcc0[itotm];
                        double xscal = chkpt.xscale[itotm];
                        double yscal = chkpt.yscale[itotm];
                        double zscal = chkpt.zscale[itotm];
                        chkpt.xccint[itotm] = xintrp;
                        chkpt.yccint[itotm] = yintrp;
                        chkpt.zccint[itotm] = zintrp;

                        double xerr = std::abs(xintrp - x0) / xscal;
                        double yerr = std::abs(yintrp - y0) / yscal;
                        double zerr = std::abs(zintrp - z0) / zscal;

                        double error = xerr;
                        if (yerr > error) error = yerr;
                        if (zerr > error) error = zerr;
                        if (error > errmax) {
                            errmax = error;
                            jerrmx = book2.jbpt[i];
                            kerrmx = book2.kbpt[i];
                            lerrmx = book2.lbpt[i];
                        }

                        if (xerr > sigdif || yerr > sigdif || zerr > sigdif) {
                            nsig = nsig + 1;
                            std::fprintf(funit(8),
                                "interpolated coordinates of fringe pt j,k,i = %12d%12d%12d mesh = %12d\n",
                                book2.jbpt[i], book2.kbpt[i], book2.lbpt[i], m);
                            std::fprintf(funit(8),
                                "differ by %s percent from the actual coordinates:\n",
                                tc::fmtE(100 * error, 15, 7).c_str());
                            std::fprintf(funit(8),
                                "        x0,y0,z0 = %s%s%s\n",
                                tc::fmtE(x0, 15, 7).c_str(),
                                tc::fmtE(y0, 15, 7).c_str(),
                                tc::fmtE(z0, 15, 7).c_str());
                            std::fprintf(funit(8),
                                "  xint,yint,zint = %s%s%s\n",
                                tc::fmtE(xintrp, 15, 7).c_str(),
                                tc::fmtE(yintrp, 15, 7).c_str(),
                                tc::fmtE(zintrp, 15, 7).c_str());
                            std::fprintf(funit(8),
                                "  xsca,ysca,zsca = %s%s%s\n",
                                tc::fmtE(xscal, 15, 7).c_str(),
                                tc::fmtE(yscal, 15, 7).c_str(),
                                tc::fmtE(zscal, 15, 7).c_str());
                            std::fprintf(funit(8),
                                "  target cell: j,k,l = %12d%12d%12d mesh = %12d\n",
                                book2.ji[i], book2.ki[i], book2.li[i], m1);
                            std::fprintf(funit(8),
                                "  xie,eta,zeta = %s%s%s\n",
                                tc::fmtE(book2.dxint[i], 15, 7).c_str(),
                                tc::fmtE(book2.dyint[i], 15, 7).c_str(),
                                tc::fmtE(book2.dzint[i], 15, 7).c_str());
                        }
                    } // 13
                }
            } // 14
        }

        // --- check boundary points ---
        if (no > 0) {
            for (int n = 1; n <= no; ++n) {
                int nnn = 1;
                int nserch = conec1.MOUTR(m, nnn);
                for (int nn = 1; nn <= nserch; ++nn) {
                    nnn = 1;
                    int m1 = conec1.LOUTR(m, nnn, nn);
                    int jmax1 = grid2.mjmax[m1];
                    int kmax1 = grid2.mkmax[m1];
                    int lmax1 = grid2.mlmax[m1];
                    getint(m1, book2.ji, book2.ki, book2.li,
                           book2.jbpt, book2.kbpt, book2.lbpt,
                           book2.dxint, book2.dyint, book2.dzint);
                    int iset = MOBS(m, m1);
                    int is = IPNTR(m1, iset);
                    int ie = NPNTR(m1, iset);
                    getgrd(m, grid1.x, grid1.y, grid1.z,
                           image.jimage, image.kimage, image.limage,
                           jmax, kmax, lmax);

                    for (int i = is; i <= ie; ++i) {
                        double x0, y0, z0;
                        bound(book2.jbpt[i], book2.kbpt[i], book2.lbpt[i],
                              x0, y0, z0, jmax, kmax, lmax,
                              grid1.x, grid1.y, grid1.z, ife);

                        int j00 = book2.jbpt[i];
                        int k00 = book2.kbpt[i];
                        int l00 = book2.lbpt[i];
                        if (j00 >= jmax) j00 = jmax - 1;
                        if (k00 >= kmax) k00 = kmax - 1;
                        if (l00 >= lmax) l00 = lmax - 1;
                        if (j00 <= 0) j00 = 1;
                        if (k00 <= 0) k00 = 1;
                        if (l00 <= 0) l00 = 1;

                        double xscal = 0., yscal = 0., zscal = 0.;
                        for (int j = j00; j <= j00 + jinc; ++j)
                        for (int k = k00; k <= k00 + kinc; ++k)
                        for (int l = l00; l <= l00 + linc; ++l) {
                            int ii = j + (k - 1) * jmax + (l - 1) * jmax * kmax;
                            double dx = std::abs(grid1.x[ii] - x0);
                            double dy = std::abs(grid1.y[ii] - y0);
                            double dz = std::abs(grid1.z[ii] - z0);
                            if (dx > xscal) xscal = dx;
                            if (dy > yscal) yscal = dy;
                            if (dz > zscal) zscal = dz;
                        }
                        if (xscal <= 0.) xscal = 1.;
                        if (yscal <= 0.) yscal = 1.;
                        if (zscal <= 0.) zscal = 1.;

                        itotn = itotn + 1;
                        chkpt.xcc0[itotn]   = x0;
                        chkpt.ycc0[itotn]   = y0;
                        chkpt.zcc0[itotn]   = z0;
                        chkpt.xscale[itotn] = xscal;
                        chkpt.yscale[itotn] = yscal;
                        chkpt.zscale[itotn] = zscal;
                    } // 161

                    getgrd(m1, grid1.x, grid1.y, grid1.z,
                           image.jimage, image.kimage, image.limage,
                           jmax1, kmax1, lmax1);

                    for (int i = is; i <= ie; ++i) {
                        double xintrp, yintrp, zintrp;
                        interp(jmax1, kmax1, lmax1, grid1.x, grid1.y, grid1.z,
                               book2.ji[i], book2.ki[i], book2.li[i],
                               book2.dxint[i], book2.dyint[i], book2.dzint[i],
                               xintrp, yintrp, zintrp);

                        itotm = itotm + 1;
                        double x0    = chkpt.xcc0[itotm];
                        double y0    = chkpt.ycc0[itotm];
                        double z0    = chkpt.zcc0[itotm];
                        double xscal = chkpt.xscale[itotm];
                        double yscal = chkpt.yscale[itotm];
                        double zscal = chkpt.zscale[itotm];
                        chkpt.xccint[itotm] = xintrp;
                        chkpt.yccint[itotm] = yintrp;
                        chkpt.zccint[itotm] = zintrp;

                        double xerr = std::abs(xintrp - x0) / xscal;
                        double yerr = std::abs(yintrp - y0) / yscal;
                        double zerr = std::abs(zintrp - z0) / zscal;

                        double error = xerr;
                        if (yerr > error) error = yerr;
                        if (zerr > error) error = zerr;
                        if (error > errmax) {
                            errmax = error;
                            jerrmx = book2.jbpt[i];
                            kerrmx = book2.kbpt[i];
                            lerrmx = book2.lbpt[i];
                        }

                        if (xerr > sigdif || yerr > sigdif || zerr > sigdif) {
                            nsig = nsig + 1;
                            std::fprintf(funit(8),
                                "interpolated coordinates of boundary pt j,k,i = %12d%12d%12d mesh = %12d\n",
                                book2.jbpt[i], book2.kbpt[i], book2.lbpt[i], m);
                            std::fprintf(funit(8),
                                "differ by %s percent from the actual coordinates:\n",
                                tc::fmtE(100 * error, 15, 7).c_str());
                            std::fprintf(funit(8),
                                "        x0,y0,z0 = %s%s%s\n",
                                tc::fmtE(x0, 15, 7).c_str(),
                                tc::fmtE(y0, 15, 7).c_str(),
                                tc::fmtE(z0, 15, 7).c_str());
                            std::fprintf(funit(8),
                                "  xint,yint,zint = %s%s%s\n",
                                tc::fmtE(xintrp, 15, 7).c_str(),
                                tc::fmtE(yintrp, 15, 7).c_str(),
                                tc::fmtE(zintrp, 15, 7).c_str());
                            std::fprintf(funit(8),
                                "  xsca,ysca,zsca = %s%s%s\n",
                                tc::fmtE(xscal, 15, 7).c_str(),
                                tc::fmtE(yscal, 15, 7).c_str(),
                                tc::fmtE(zscal, 15, 7).c_str());
                            std::fprintf(funit(8),
                                "  target cell: j,k,l = %12d%12d%12d mesh = %12d\n",
                                book2.ji[i], book2.ki[i], book2.li[i], m1);
                            std::fprintf(funit(8),
                                "  xie,eta,zeta = %s%s%s\n",
                                tc::fmtE(book2.dxint[i], 15, 7).c_str(),
                                tc::fmtE(book2.dyint[i], 15, 7).c_str(),
                                tc::fmtE(book2.dzint[i], 15, 7).c_str());
                        }
                    } // 15
                }
            } // 16
        }

        if (nh + no > 0) {
            std::fprintf(funit(6), " \n");
            // write(6,*) '      ',nsig,' cell centers of mesh ',m,' calculated using the transfer stencils'
            std::fprintf(funit(6),
                "       %12d cell centers of mesh %12d calculated using the transfer stencils\n",
                nsig, m);
            // write(6,*) '      differ by more than ',int(sigdif*100),' percent from the actual cell centers'
            std::fprintf(funit(6),
                "       differ by more than %12d percent from the actual cell centers\n",
                (int)(sigdif * 100));
            // write(6,*) '      the maximum deviation ',100*errmax,' percent occurs at ',jerrmx,kerrmx,lerrmx
            std::fprintf(funit(6),
                "       the maximum deviation %s percent occurs at %12d%12d%12d\n",
                tc::fmtE(100 * errmax, 15, 7).c_str(), jerrmx, kerrmx, lerrmx);
        }
    } // 1000
}

// ===========================================================================
//  subroutine chkstn( iblank )
// ===========================================================================
void chkstn(int* iblank)
{
    const int mdim2 = mdim * mdim;

    View2<int> IPNTR(book1.ipntr, mdim, mhldim * mdim);
    View2<int> NPNTR(book1.npntr, mdim, mhldim * mdim);
    View2<int> MHBS (book1.mhbs,  mdim, mdim);
    View2<int> MOBS (book1.mobs,  mdim, mdim);

    bool iholcor;

    int list[mdim + 1];
    // dimension mhhtrn(mdim,mdim) etc; data ...*0
    static int mhhtrn[mdim * mdim + 1];
    static int mhotrn[mdim * mdim + 1];
    static int mohtrn[mdim * mdim + 1];
    static int mootrn[mdim * mdim + 1];
    static int jinew[idim + 1], kinew[idim + 1], linew[idim + 1];
    int ife[6 + 1];

    // data mhhtrn/mdim2*0/,... : zero-initialise (also matches per-call reset
    // semantics — Fortran DATA is one-time init, but these are only accumulated
    // within this single call before being printed, then never read again.)
    for (int i = 1; i <= mdim2; ++i) { mhhtrn[i] = 0; mhotrn[i] = 0; mohtrn[i] = 0; mootrn[i] = 0; }
    View2<int> MHHTRN(mhhtrn, mdim, mdim);
    View2<int> MHOTRN(mhotrn, mdim, mdim);
    View2<int> MOHTRN(mohtrn, mdim, mdim);
    View2<int> MOOTRN(mootrn, mdim, mdim);

    iholcor = true;

    int iflag1 = 0;
    int iflag2 = 0;

    for (int ll = 1; ll <= 6; ++ll) ife[ll] = 0;

    for (int m = 1; m <= conec1.nmesh; ++m) {

        // --- check fringe-point stencils ---
        int nh = conec1.nhole[m];

        if (nh > 0) {
            int jmaxm = grid2.mjmax[m];
            int kmaxm = grid2.mkmax[m];
            int lmaxm = grid2.mlmax[m];
            for (int n = 1; n <= nh; ++n) {
                int nserch = conec1.MHOLE(m, n);
                for (int nn = 1; nn <= nserch; ++nn) {
                    getgrd(m, grid1.x, grid1.y, grid1.z,
                           image.jimage, image.kimage, image.limage,
                           jmaxm, kmaxm, lmaxm);
                    int ioh = 0;
                    int ihh = 0;
                    int m1 = conec1.LHOLE(m, n, nn);
                    getint(m1, book2.ji, book2.ki, book2.li,
                           book2.jbpt, book2.kbpt, book2.lbpt,
                           book2.dxint, book2.dyint, book2.dzint);
                    int iset = MHBS(m, m1);
                    int is = IPNTR(m1, iset);
                    int ie = NPNTR(m1, iset);
                    int ic = ie - is + 1;
                    if (ic > 0) {
                        int jmax = grid2.mjmax[m1];
                        int kmax = grid2.mkmax[m1];
                        int lmax = grid2.mlmax[m1];
                        int jmaxx = jmax - 1;
                        int kmaxx = kmax - 1;
                        int lmaxx = lmax - 1;
                        if (jmax < 2) jmaxx = 1;
                        if (kmax < 2) kmaxx = 1;
                        if (lmax < 2) lmaxx = 1;
                        for (int i = is; i <= ie; ++i) {
                            if (book2.dxint[i] == 0. && book2.dyint[i] == 0. &&
                                book2.dzint[i] == 0.) {
                                jmaxx = jmax;
                                kmaxx = kmax;
                                lmaxx = lmax;
                            }
                            if (book2.ji[i] < 1 || book2.ji[i] > jmaxx ||
                                book2.ki[i] < 1 || book2.ki[i] > kmaxx ||
                                book2.li[i] < 1 || book2.li[i] > lmaxx) {
                                std::fprintf(funit(6),
                                    "  stopping...illegal target cell for fringe point %12d%12d%12d of mesh %12d\n",
                                    book2.jbpt[i], book2.kbpt[i], book2.lbpt[i], m);
                                std::fprintf(funit(6),
                                    "    target cell is %12d%12d%12d in mesh %12d\n",
                                    book2.ji[i], book2.ki[i], book2.li[i], m1);
                                std::fprintf(funit(6),
                                    "    legal range in mesh %12d is:  1<j<%12d  1<k<%12d  1<l<%12d\n",
                                    m1, jmaxx, kmaxx, lmaxx);
                                std::fprintf(funit(6),
                                    "    possible miscue in search routine\n");
                                std::exit(1);   // stop
                            }
                        } // 1600

                        for (int i = 1; i <= ic; ++i) {
                            int ii = i + is - 1;
                            double x0, y0, z0;
                            bound(book2.jbpt[ii], book2.kbpt[ii], book2.lbpt[ii],
                                  x0, y0, z0, jmaxm, kmaxm, lmaxm,
                                  grid1.x, grid1.y, grid1.z, ife);
                            surf.xb[i] = x0;
                            surf.yb[i] = y0;
                            surf.zb[i] = z0;
                        } // 131
                        getibl(m1, iblank, jmax, kmax, lmax);
                        getgrd(m1, grid1.x, grid1.y, grid1.z,
                               image.jimage, image.kimage, image.limage,
                               jmax, kmax, lmax);
                        ibchk(iworkc.ibck, iblank, ic, is,
                              book2.ji, book2.ki, book2.li, m1, jmax, kmax, lmax,
                              surf.xb, surf.yb, surf.zb,
                              jinew, kinew, linew,
                              grid1.x, grid1.y, grid1.z);
                        for (int i = 1; i <= ic; ++i) {
                            if (iworkc.ibck[i] > 0) {
                                // stencil contains interpolated boundary points
                                ioh    = ioh + 1;
                                iflag1 = iflag1 + 1;
                                int ii = i + is - 1;
                                pltpt.nill = pltpt.nill + 1;
                                pltpt.xill[pltpt.nill] = surf.xb[i];
                                pltpt.yill[pltpt.nill] = surf.yb[i];
                                pltpt.zill[pltpt.nill] = surf.zb[i];

                                // 702 format
                                std::fprintf(funit(8),
                                    " the fringe point jb,kb,lb %4d%4d%4d of mesh %3d\n"
                                    "    has a stencil from point ji,ki,li %4d%4d%4d which lies on an interp bndry of mesh %3d\n",
                                    book2.jbpt[ii], book2.kbpt[ii], book2.lbpt[ii], m,
                                    book2.ji[ii], book2.ki[ii], book2.li[ii], m1);

                                if (iholcor) {
                                    if (iworkc.ibck[i] > 1) {
                                        std::fprintf(funit(8),
                                            "       this illegal transfer was circumvented using zeroth order data transfer\n");
                                        std::fprintf(funit(8),
                                            "         old ji,ki,li: %12d%12d%12d\n",
                                            book2.ji[ii], book2.ki[ii], book2.li[ii]);
                                        std::fprintf(funit(8),
                                            "         new ji,ki,li: %12d%12d%12d\n",
                                            jinew[ii], kinew[ii], linew[ii]);
                                        book2.dxint[ii] = 0.;
                                        book2.dyint[ii] = 0.;
                                        book2.dzint[ii] = 0.;
                                        book2.ji[ii] = jinew[ii];
                                        book2.ki[ii] = kinew[ii];
                                        book2.li[ii] = linew[ii];
                                        iflag2 = iflag2 + 1;
                                    }
                                }
                            } else if (iworkc.ibck[i] < 0) {
                                // stencil contains hole and/or fringe points
                                ihh    = ihh + 1;
                                iflag1 = iflag1 + 1;
                                int ii = i + is - 1;
                                pltpt.nill = pltpt.nill + 1;
                                pltpt.xill[pltpt.nill] = surf.xb[i];
                                pltpt.yill[pltpt.nill] = surf.yb[i];
                                pltpt.zill[pltpt.nill] = surf.zb[i];

                                // 802 format
                                std::fprintf(funit(8),
                                    " the fringe point jb,kb,lb %4d%4d%4d of mesh %3d\n"
                                    "    has a stencil point ji,ki,li %4d%4d%4d which lies in a hole of mesh %3d\n",
                                    book2.jbpt[ii], book2.kbpt[ii], book2.lbpt[ii], m,
                                    book2.ji[ii], book2.ki[ii], book2.li[ii], m1);

                                if (iholcor) {
                                    if (iworkc.ibck[i] < -1) {
                                        std::fprintf(funit(8),
                                            "       this illegal transfer was circumvented using zeroth order data transfer\n");
                                        std::fprintf(funit(8),
                                            "         old ji,ki,li: %12d%12d%12d\n",
                                            book2.ji[ii], book2.ki[ii], book2.li[ii]);
                                        std::fprintf(funit(8),
                                            "         new ji,ki,li: %12d%12d%12d\n",
                                            jinew[ii], kinew[ii], linew[ii]);
                                        book2.dxint[ii] = 0.;
                                        book2.dyint[ii] = 0.;
                                        book2.dzint[ii] = 0.;
                                        book2.ji[ii] = jinew[ii];
                                        book2.ki[ii] = kinew[ii];
                                        book2.li[ii] = linew[ii];
                                        iflag2 = iflag2 + 1;
                                    }
                                }
                            }
                        } // 10
                        putint(chkst.nwr, m1, book2.ji, book2.ki, book2.li,
                               book2.jbpt, book2.kbpt, book2.lbpt,
                               book2.dxint, book2.dyint, book2.dzint);
                    }
                    MHHTRN(m, m1) = ihh + MHHTRN(m, m1);
                    MOHTRN(m, m1) = ioh + MOHTRN(m, m1);
                } // 141
            } // 14
        }

        // --- check boundary-point stencils ---
        int no = conec1.noutr[m];

        if (no > 0) {
            int jmaxm = grid2.mjmax[m];
            int kmaxm = grid2.mkmax[m];
            int lmaxm = grid2.mlmax[m];
            for (int n = 1; n <= no; ++n) {
                int nnn = 1;
                int nserch = conec1.MOUTR(m, nnn);
                for (int nn = 1; nn <= nserch; ++nn) {
                    getgrd(m, grid1.x, grid1.y, grid1.z,
                           image.jimage, image.kimage, image.limage,
                           jmaxm, kmaxm, lmaxm);
                    int ioo = 0;
                    int iho = 0;
                    nnn = 1;
                    int m1 = conec1.LOUTR(m, nnn, nn);
                    getint(m1, book2.ji, book2.ki, book2.li,
                           book2.jbpt, book2.kbpt, book2.lbpt,
                           book2.dxint, book2.dyint, book2.dzint);
                    int iset = MOBS(m, m1);
                    int is = IPNTR(m1, iset);
                    int ie = NPNTR(m1, iset);
                    int ic = ie - is + 1;
                    if (ic > 0) {
                        int jmax = grid2.mjmax[m1];
                        int kmax = grid2.mkmax[m1];
                        int lmax = grid2.mlmax[m1];
                        int jmaxx = jmax - 1;
                        int kmaxx = kmax - 1;
                        int lmaxx = lmax - 1;
                        if (jmax < 2) jmaxx = 1;
                        if (kmax < 2) kmaxx = 1;
                        if (lmax < 2) lmaxx = 1;
                        for (int i = is; i <= ie; ++i) {
                            if (book2.dxint[i] == 0. && book2.dyint[i] == 0. &&
                                book2.dzint[i] == 0.) {
                                jmaxx = jmax;
                                kmaxx = kmax;
                                lmaxx = lmax;
                            }
                            if (book2.ji[i] < 1 || book2.ji[i] > jmaxx ||
                                book2.ki[i] < 1 || book2.ki[i] > kmaxx ||
                                book2.li[i] < 1 || book2.li[i] > lmaxx) {
                                std::fprintf(funit(6),
                                    "  stopping...illegal target cell for boundary point %12d%12d%12d of mesh %12d\n",
                                    book2.jbpt[i], book2.kbpt[i], book2.lbpt[i], m);
                                std::fprintf(funit(6),
                                    "    target cell is %12d%12d%12d in mesh %12d\n",
                                    book2.ji[i], book2.ki[i], book2.li[i], m1);
                                std::fprintf(funit(6),
                                    "    legal range in mesh %12d is:  1<j<%12d  1<k<%12d  1<l<%12d\n",
                                    m1, jmaxx, kmaxx, lmaxx);
                                std::fprintf(funit(6),
                                    "    possible miscue in search routine\n");
                                std::exit(1);   // stop
                            }
                        } // 2600

                        for (int i = 1; i <= ic; ++i) {
                            int ii = i + is - 1;
                            double x0, y0, z0;
                            bound(book2.jbpt[ii], book2.kbpt[ii], book2.lbpt[ii],
                                  x0, y0, z0, jmaxm, kmaxm, lmaxm,
                                  grid1.x, grid1.y, grid1.z, ife);
                            surf.xb[i] = x0;
                            surf.yb[i] = y0;
                            surf.zb[i] = z0;
                        } // 231
                        getibl(m1, iblank, jmax, kmax, lmax);
                        getgrd(m1, grid1.x, grid1.y, grid1.z,
                               image.jimage, image.kimage, image.limage,
                               jmax, kmax, lmax);
                        ibchk(iworkc.ibck, iblank, ic, is,
                              book2.ji, book2.ki, book2.li, m1, jmax, kmax, lmax,
                              surf.xb, surf.yb, surf.zb,
                              jinew, kinew, linew,
                              grid1.x, grid1.y, grid1.z);
                        for (int i = 1; i <= ic; ++i) {
                            if (iworkc.ibck[i] > 0) {
                                // stencil contains interpolated boundary points
                                ioo    = ioo + 1;
                                iflag1 = iflag1 + 1;
                                int ii = i + is - 1;
                                pltpt.nill = pltpt.nill + 1;
                                pltpt.xill[pltpt.nill] = surf.xb[i];
                                pltpt.yill[pltpt.nill] = surf.yb[i];
                                pltpt.zill[pltpt.nill] = surf.zb[i];

                                // 502 format
                                std::fprintf(funit(8),
                                    " the boundary point jb,kb,lb %4d%4d%4d of mesh %3d\n"
                                    "    has a stencil from point ji,ki,li %4d%4d%4d which lies on an interp bndry of mesh %3d\n",
                                    book2.jbpt[ii], book2.kbpt[ii], book2.lbpt[ii], m,
                                    book2.ji[ii], book2.ki[ii], book2.li[ii], m1);

                                if (iholcor) {
                                    if (iworkc.ibck[i] > 1) {
                                        std::fprintf(funit(8),
                                            "       this illegal transfer was circumvented using zeroth order data transfer\n");
                                        std::fprintf(funit(8),
                                            "         old ji,ki,li: %12d%12d%12d\n",
                                            book2.ji[ii], book2.ki[ii], book2.li[ii]);
                                        std::fprintf(funit(8),
                                            "         new ji,ki,li: %12d%12d%12d\n",
                                            jinew[ii], kinew[ii], linew[ii]);
                                        book2.dxint[ii] = 0.;
                                        book2.dyint[ii] = 0.;
                                        book2.dzint[ii] = 0.;
                                        book2.ji[ii] = jinew[ii];
                                        book2.ki[ii] = kinew[ii];
                                        book2.li[ii] = linew[ii];
                                        iflag2 = iflag2 + 1;
                                    }
                                }
                            } else if (iworkc.ibck[i] < 0) {
                                // stencil contains hole and/or fringe points
                                iho    = iho + 1;
                                iflag1 = iflag1 + 1;
                                int ii = i + is - 1;
                                pltpt.nill = pltpt.nill + 1;
                                pltpt.xill[pltpt.nill] = surf.xb[i];
                                pltpt.yill[pltpt.nill] = surf.yb[i];
                                pltpt.zill[pltpt.nill] = surf.zb[i];

                                // 602 format
                                std::fprintf(funit(8),
                                    " the boundary point jb,kb,lb %4d%4d%4d of mesh %3d\n"
                                    "    has a stencil point ji,ki,li %4d%4d%4d which lies in a hole of mesh %3d\n",
                                    book2.jbpt[ii], book2.kbpt[ii], book2.lbpt[ii], m,
                                    book2.ji[ii], book2.ki[ii], book2.li[ii], m1);

                                if (iholcor) {
                                    if (iworkc.ibck[i] < -1) {
                                        std::fprintf(funit(8),
                                            "       this illegal transfer was circumvented using zeroth order data transfer\n");
                                        std::fprintf(funit(8),
                                            "         old ji,ki,li: %12d%12d%12d\n",
                                            book2.ji[ii], book2.ki[ii], book2.li[ii]);
                                        std::fprintf(funit(8),
                                            "         new ji,ki,li: %12d%12d%12d\n",
                                            jinew[ii], kinew[ii], linew[ii]);
                                        book2.dxint[ii] = 0.;
                                        book2.dyint[ii] = 0.;
                                        book2.dzint[ii] = 0.;
                                        book2.ji[ii] = jinew[ii];
                                        book2.ki[ii] = kinew[ii];
                                        book2.li[ii] = linew[ii];
                                        iflag2 = iflag2 + 1;
                                    }
                                }
                            }
                        } // 20
                        putint(chkst.nwr, m1, book2.ji, book2.ki, book2.li,
                               book2.jbpt, book2.kbpt, book2.lbpt,
                               book2.dxint, book2.dyint, book2.dzint);
                    }
                    MHOTRN(m, m1) = iho + MHOTRN(m, m1);
                    MOOTRN(m, m1) = ioo + MOOTRN(m, m1);
                } // 241
            } // 24
        }
    } // 31

    // --- output summary of stencil checks ---
    std::fprintf(funit(6), " \n");
    std::fprintf(funit(6),
        "     checking for four types of illegal transfer of information:\n");
    std::fprintf(funit(6),
        "     1. fringe points receive data from hole/fringe points\n");
    std::fprintf(funit(6),
        "     2. fringe points receive data from boundary points\n");
    std::fprintf(funit(6),
        "     3. boundary points receive data from hole/fringe points\n");
    std::fprintf(funit(6),
        "     4. boundary pts receive data from boundary points\n");

    for (int m = 1; m <= conec1.nmesh; ++m) {
        std::fprintf(funit(6), " \n");
        // write(6,*) '    checking stencils for mesh ',m
        std::fprintf(funit(6), "     checking stencils for mesh %12d\n", m);

        std::fprintf(funit(6), "       checking stencils of fringe points\n");
        int nh = conec1.nhole[m];
        for (int m1 = 1; m1 <= conec1.nmesh; ++m1) list[m1] = 0;
        if (nh <= 0) {
            std::fprintf(funit(6),
                "         there are no fringe pts in mesh %12d to check\n", m);
        } else {
            for (int n = 1; n <= nh; ++n) {
                int nserch = conec1.MHOLE(m, n);
                for (int nn = 1; nn <= nserch; ++nn) {
                    int m1 = conec1.LHOLE(m, n, nn);
                    if (list[m1] == 0) list[m1] = 1;
                }
            } // 612
            for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
                if (list[m1] > 0) {
                    if (conec1.nhole[m1] > 0) {
                        std::fprintf(funit(6),
                            "         1. %12d fringe pts in mesh %12d receive data from hole/fringe pts in mesh %12d\n",
                            MHHTRN(m, m1), m, m1);
                    } else {
                        std::fprintf(funit(6),
                            "         1. not applicable, no hole/fringe pts in mesh %12d\n", m1);
                    }
                }
            } // 501
            for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
                if (list[m1] > 0) {
                    if (conec1.noutr[m1] > 0) {
                        std::fprintf(funit(6),
                            "         2. %12d fringe pts in mesh %12d receive data from boundary pts in mesh %12d\n",
                            MOHTRN(m, m1), m, m1);
                    } else {
                        std::fprintf(funit(6),
                            "         2. not applicable, no boundary pts in mesh %12d\n", m1);
                    }
                }
            } // 701
        }

        std::fprintf(funit(6), "       checking stencils of boundary points\n");
        int no = conec1.noutr[m];
        for (int m1 = 1; m1 <= conec1.nmesh; ++m1) list[m1] = 0;
        if (no <= 0) {
            std::fprintf(funit(6),
                "         there are no boundary pts in mesh %12d to check\n", m);
        } else {
            for (int n = 1; n <= no; ++n) {
                int nnn = 1;
                int nserch = conec1.MOUTR(m, nnn);
                for (int nn = 1; nn <= nserch; ++nn) {
                    nnn = 1;
                    int m1 = conec1.LOUTR(m, nnn, nn);
                    if (list[m1] == 0) list[m1] = 1;
                }
            } // 812
            for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
                if (list[m1] > 0) {
                    if (conec1.nhole[m1] > 0) {
                        std::fprintf(funit(6),
                            "         3. %12d bndry pts in mesh %12d receive data from hole/fringe pts in mesh %12d\n",
                            MHOTRN(m, m1), m, m1);
                    } else {
                        std::fprintf(funit(6),
                            "         3. not applicable, no hole/fringe pts in mesh %12d\n", m1);
                    }
                }
            } // 801
            for (int m1 = 1; m1 <= conec1.nmesh; ++m1) {
                if (list[m1] > 0) {
                    if (conec1.noutr[m1] > 0) {
                        std::fprintf(funit(6),
                            "         4. %12d bndry pts in mesh %12d receive data from bndry pts in mesh %12d\n",
                            MOOTRN(m, m1), m, m1);
                    } else {
                        std::fprintf(funit(6),
                            "         4. not applicable, no boundary pts in mesh %12d\n", m1);
                    }
                }
            } // 901
        }
    } // 401

    if (iflag1 > 0) {
        std::fprintf(funit(6), " \n");
        std::fprintf(funit(6),
            "     ***CAUTION***CAUTION***CAUTION***CAUTION***CAUTION***CAUTION***CAUTION***\n");
        std::fprintf(funit(6),
            "     of the %12d illegal stencils detected above, %12d were circumvented\n",
            iflag1, iflag2);
        std::fprintf(funit(6),
            "     by  using zeroth order data transfer. illegal transfer of information\n");
        std::fprintf(funit(6),
            "     will no longer occur with these stencils, but accuracy and convergence\n");
        std::fprintf(funit(6),
            "     may suffer. a large geometric mismatch (see below, and also file 8)\n");
        std::fprintf(funit(6),
            "     indicates loss of accuracy. if possible, adjust maggie input parameters\n");
        std::fprintf(funit(6),
            "     and/or component grids to reduce the number of stencils needing correction\n");
        std::fprintf(funit(6),
            "     ***CAUTION***CAUTION***CAUTION***CAUTION***CAUTION***CAUTION***CAUTION***\n");
    }

    if (iflag2 != iflag1) {
        std::fprintf(funit(6), " \n");
        std::fprintf(funit(6),
            "     ***WARNINIG***WARNING***WARNING***WARNING***WARNING***WARNING***WARNING***\n");
        std::fprintf(funit(6),
            "                            DO NOT RUN FLOW SOLVER!!!\n");
        std::fprintf(funit(6),
            "     %12d illegal stencils could not be circumvented...must either adjust\n",
            iflag1 - iflag2);
        std::fprintf(funit(6),
            "     maggie input parameters and/or component grids to rectify this problem.\n");
        std::fprintf(funit(6),
            "                            DO NOT RUN FLOW SOLVER!!!\n");
        std::fprintf(funit(6),
            "     ***WARNINIG***WARNING***WARNING***WARNING***WARNING***WARNING***WARNING***\n");
    }
}

// ===========================================================================
//  subroutine ibchk(...)
// ===========================================================================
void ibchk(int* ib, int* iblank, int ic, int is, int* ji, int* ki, int* li,
           int m1, int jd, int kd, int ld,
           double* xb, double* yb, double* zb,
           int* jinew, int* kinew, int* linew,
           double* x, double* y, double* z)
{
    (void)m1;
    View3<int>    IBLANK(iblank, jd, kd, ld);
    View3<double> X(x, jd, kd, ld);
    View3<double> Y(y, jd, kd, ld);
    View3<double> Z(z, jd, kd, ld);

    double dmin0 = 1.e30;

    int ii = is - 1;
    for (int i = 1; i <= ic; ++i) {
        double dmin = dmin0;
        ib[i] = 0;
        ii = ii + 1;
        int j = ji[ii];
        int k = ki[ii];
        int l = li[ii];
        int jp1 = std::min(j + 1, jd);
        int kp1 = std::min(k + 1, kd);
        int lp1 = std::min(l + 1, ld);
        int i1 = IBLANK(j, k, l);
        int i2 = IBLANK(jp1, k, l);
        int i3 = IBLANK(jp1, kp1, l);
        int i4 = IBLANK(j, kp1, l);
        int i5 = IBLANK(j, k, lp1);
        int i6 = IBLANK(jp1, k, lp1);
        int i7 = IBLANK(jp1, kp1, lp1);
        int i8 = IBLANK(j, kp1, lp1);
        int ib1 = std::min({i1, i2, i3, i4, i5, i6, i7, i8});
        int ib2 = std::max({i1, i2, i3, i4, i5, i6, i7, i8});

        if (ib1 * ib2 != 1) {
            int jmin = 0, kmin = 0, lmin = 0;

            if (i1 == 1) {
                double d1 = (X(j, k, l) - xb[ii]) * (X(j, k, l) - xb[ii])
                          + (Y(j, k, l) - yb[ii]) * (Y(j, k, l) - yb[ii])
                          + (Z(j, k, l) - zb[ii]) * (Z(j, k, l) - zb[ii]);
                if (d1 < dmin) { dmin = d1; jmin = j; kmin = k; lmin = l; }
            }
            if (i2 == 1) {
                double d1 = (X(jp1, k, l) - xb[ii]) * (X(jp1, k, l) - xb[ii])
                          + (Y(jp1, k, l) - yb[ii]) * (Y(jp1, k, l) - yb[ii])
                          + (Z(jp1, k, l) - zb[ii]) * (Z(jp1, k, l) - zb[ii]);
                if (d1 < dmin) { dmin = d1; jmin = jp1; kmin = k; lmin = l; }
            }
            if (i3 == 1) {
                double d1 = (X(jp1, kp1, l) - xb[ii]) * (X(jp1, kp1, l) - xb[ii])
                          + (Y(jp1, kp1, l) - yb[ii]) * (Y(jp1, kp1, l) - yb[ii])
                          + (Z(jp1, kp1, l) - zb[ii]) * (Z(jp1, kp1, l) - zb[ii]);
                if (d1 < dmin) { dmin = d1; jmin = jp1; kmin = kp1; lmin = l; }
            }
            if (i4 == 1) {
                double d1 = (X(j, kp1, l) - xb[ii]) * (X(j, kp1, l) - xb[ii])
                          + (Y(j, kp1, l) - yb[ii]) * (Y(j, kp1, l) - yb[ii])
                          + (Z(j, kp1, l) - zb[ii]) * (Z(j, kp1, l) - zb[ii]);
                if (d1 < dmin) { dmin = d1; jmin = j; kmin = kp1; lmin = l; }
            }
            if (i5 == 1) {
                double d1 = (X(j, k, lp1) - xb[ii]) * (X(j, k, lp1) - xb[ii])
                          + (Y(j, k, lp1) - yb[ii]) * (Y(j, k, lp1) - yb[ii])
                          + (Z(j, k, lp1) - zb[ii]) * (Z(j, k, lp1) - zb[ii]);
                if (d1 < dmin) { dmin = d1; jmin = j; kmin = k; lmin = lp1; }
            }
            if (i6 == 1) {
                double d1 = (X(jp1, k, lp1) - xb[ii]) * (X(jp1, k, lp1) - xb[ii])
                          + (Y(jp1, k, lp1) - yb[ii]) * (Y(jp1, k, lp1) - yb[ii])
                          + (Z(jp1, k, lp1) - zb[ii]) * (Z(jp1, k, lp1) - zb[ii]);
                if (d1 < dmin) { dmin = d1; jmin = jp1; kmin = k; lmin = lp1; }
            }
            if (i7 == 1) {
                double d1 = (X(jp1, kp1, lp1) - xb[ii]) * (X(jp1, kp1, lp1) - xb[ii])
                          + (Y(jp1, kp1, lp1) - yb[ii]) * (Y(jp1, kp1, lp1) - yb[ii])
                          + (Z(jp1, kp1, lp1) - zb[ii]) * (Z(jp1, kp1, lp1) - zb[ii]);
                if (d1 < dmin) { dmin = d1; jmin = jp1; kmin = kp1; lmin = lp1; }
            }
            if (i8 == 1) {
                double d1 = (X(j, kp1, lp1) - xb[ii]) * (X(j, kp1, lp1) - xb[ii])
                          + (Y(j, kp1, lp1) - yb[ii]) * (Y(j, kp1, lp1) - yb[ii])
                          + (Z(j, kp1, lp1) - zb[ii]) * (Z(j, kp1, lp1) - zb[ii]);
                if (d1 < dmin) { dmin = d1; jmin = j; kmin = kp1; lmin = lp1; }
            }

            if (ib1 <= 0) {
                // stencil contains hole/fringe points
                if (dmin < dmin0) {
                    ib[i] = -2;
                    jinew[ii] = jmin;
                    kinew[ii] = kmin;
                    linew[ii] = lmin;
                } else {
                    ib[i] = -1;
                }
            }

            if (ib2 > 1) {
                // stencil contains interpolated boundary points
                if (dmin < dmin0) {
                    ib[i] = 2;
                    jinew[ii] = jmin;
                    kinew[ii] = kmin;
                    linew[ii] = lmin;
                } else {
                    ib[i] = 1;
                }
            }
        }
    } // 11
}
