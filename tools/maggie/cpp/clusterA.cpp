// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// clusterA.cpp — maggie "Cluster A: driver + setup + grid I/O"
// Faithful byte-exact C++ port of the Fortran-77 maggie subroutines:
//   maggie (main program), blckdat, setup, initia, rgrid, rp3d,
//   cellcen, wcelgrd, getgrd, getgrd2, branch
//
// Contract: mag_common.h pulls in mag_protos.h + tools_common.h + all COMMON
// globals. COMMON blocks are global structs; arrays are 1-based, column-major.
// Loop nesting / operation order / parenthesization kept identical to Fortran
// for FP byte-exactness (built with -ffp-contract=off).
//
// SIGNATURE: int main(int argc, char** argv)                              // program maggie
// SIGNATURE: void blckdat()
// SIGNATURE: void setup()
// SIGNATURE: void initia()
// SIGNATURE: void rgrid(int m,double* x,double* y,double* z,int jd,int kd,int ld)
// SIGNATURE: void rp3d(int m,double* x,double* y,double* z,int jd,int kd,int ld,int& ialph,int nmesh)
// SIGNATURE: void cellcen(double* x,double* y,double* z,int jd,int kd,int ld)
// SIGNATURE: void wcelgrd(int m,double* x,double* y,double* z,int* jimage,int* kimage,int* limage,int jd,int kd,int ld)
// SIGNATURE: void getgrd(int m,double* x,double* y,double* z,int* jimage,int* kimage,int* limage,int jd,int kd,int ld)
// SIGNATURE: void getgrd2(int m,double* x,double* y,double* z,int jd,int kd,int ld)
// SIGNATURE: void branch(int m1,int jd,int kd,int ld,int* jimage,int* kimage,int* limage,double* x1,double* y1,double* z1)
//
// mag_protos.h signature-vs-Fortran mismatches: NONE for this cluster.
//   - rp3d: mag_protos declares `int& ialph` (out/inout scalar); Fortran passes
//     ialph by reference and the routine assigns ialph=1 → matches (inout ok).
//   - branch: mag_protos declares jimage/kimage/limage as int*, x/y/z as double*
//     named x,y,z; Fortran names them x1,y1,z1. Same types/order → matches.
//   - wcelgrd/getgrd/getgrd2/rgrid: all match mag_protos.h exactly.

#include "mag_common.h"
#include <cstdlib>

// ---------------------------------------------------------------------------
// helper: build the scratch filename "temp_XXX.<m>" exactly as the Fortran
//   write(titl,'("temp_grd.",i1/i2/i3)') m does, returning length used.
//   Fortran pads titl(len+1:20) with blanks and opens titl(1:len).
// ---------------------------------------------------------------------------
static int build_scratch_name(char titl[21], const char* prefix, int m) {
    int len;
    if (m > 99) {
        len = 12;
        std::snprintf(titl, 21, "%s%3d", prefix, m);   // i3
    } else if (m > 9) {
        len = 11;
        std::snprintf(titl, 21, "%s%2d", prefix, m);   // i2
    } else {
        len = 10;
        std::snprintf(titl, 21, "%s%1d", prefix, m);   // i1
    }
    for (int i = len; i < 20; ++i) titl[i] = ' ';      // titl(len+1:20)=' '
    titl[20] = '\0';
    // filename actually opened is titl(1:len)
    return len;
}

// ===========================================================================
//  subroutine blckdat  — variables/scalars given initial or default values
// ===========================================================================
void blckdat() {
    // do 10 m=1,mdim : nsets(m)=0
    for (int m = 1; m <= mdim; ++m) book1.nsets[m] = 0;

    // do 20 n=1,mhldim ; do 20 m=1,mdim : ipntr(m,n)=0 ; npntr(m,n)=0
    for (int n = 1; n <= mhldim; ++n)
        for (int m = 1; m <= mdim; ++m) {
            book1.IPNTR(m, n) = 0;
            book1.NPNTR(m, n) = 0;
        }

    // do 30 m=1,mdim : ibpts(m)=0
    for (int m = 1; m <= mdim; ++m) book2.ibpts[m] = 0;

    // do 40 i=1,idim
    for (int i = 1; i <= idim; ++i) {
        book2.jbpt[i]  = 0;
        book2.kbpt[i]  = 0;
        book2.lbpt[i]  = 0;
        book2.ji[i]    = 0;
        book2.ki[i]    = 0;
        book2.li[i]    = 0;
        book2.dxint[i] = 0.0;
        book2.dyint[i] = 0.0;
        book2.dzint[i] = 0.0;
    }

    // do 50 i=1,mlen : iblank(i)=1
    for (int i = 1; i <= mlen; ++i) book3.iblank[i] = 1;

    // do 60 m=1,mdim : nhole(m)=0 ; noutr(m)=0
    for (int m = 1; m <= mdim; ++m) {
        conec1.nhole[m] = 0;
        conec1.noutr[m] = 0;
    }

    // do 70 m=1,mdim ; do 70 n=1,mhldim : moutr(m,n)=0 ; mhole(m,n)=0
    for (int m = 1; m <= mdim; ++m)
        for (int n = 1; n <= mhldim; ++n) {
            conec1.MOUTR(m, n) = 0;
            conec1.MHOLE(m, n) = 0;
        }

    // do 80 m=1,mdim ; do 80 n=1,mhldim ; do 80 nn=1,mdim
    for (int m = 1; m <= mdim; ++m)
        for (int n = 1; n <= mhldim; ++n)
            for (int nn = 1; nn <= mdim; ++nn) {
                conec1.LOUTR(m, n, nn) = 0;
                conec1.LHOLE(m, n, nn) = 0;
            }

    // do 90 m=1,mdim
    for (int m = 1; m <= mdim; ++m) {
        grdlim.jsrs[m] = 0;
        grdlim.jsre[m] = 0;
        grdlim.ksrs[m] = 0;
        grdlim.ksre[m] = 0;
        grdlim.lsrs[m] = 0;
        grdlim.lsre[m] = 0;
    }

    // do 100 m=1,mlen : x(m)=0 ; y(m)=0 ; z(m)=0
    for (int m = 1; m <= mlen; ++m) {
        grid1.x[m] = 0.0;
        grid1.y[m] = 0.0;
        grid1.z[m] = 0.0;
    }

    // do 110 m=1,mdim : mjmax=mkmax=mlmax=0
    for (int m = 1; m <= mdim; ++m) {
        grid2.mjmax[m] = 0;
        grid2.mkmax[m] = 0;
        grid2.mlmax[m] = 0;
    }

    // do 120 n=1,mdim ; do 120 m=1,mdim : ihbtyp(m,n)=0
    for (int n = 1; n <= mdim; ++n)
        for (int m = 1; m <= mdim; ++m) bound1.IHBTYP(m, n) = 0;

    // do 130 ip=1,ipmax
    for (int ip = 1; ip <= ipmax; ++ip) {
        bound1.jh1[ip] = 0;
        bound1.jh2[ip] = 0;
        bound1.kh1[ip] = 0;
        bound1.kh2[ip] = 0;
        bound1.lh1[ip] = 0;
        bound1.lh2[ip] = 0;
        bound2.jo1[ip] = 0;
        bound2.jo2[ip] = 0;
        bound2.ko1[ip] = 0;
        bound2.ko2[ip] = 0;
        bound2.lo1[ip] = 0;
        bound2.lo2[ip] = 0;
    }

    // do 140 m=1,mdim : iobtyp(m)=0 ; nobtyp(m)=0
    for (int m = 1; m <= mdim; ++m) {
        bound2.iobtyp[m] = 0;
        bound2.nobtyp[m] = 0;
    }

    return;
}

// ===========================================================================
//  subroutine setup  — set up for overlapped grid interpolation
// ===========================================================================
void setup() {
    double* x = grid1.x;   // common /grid1/ x(mlen),y(mlen),z(mlen)
    double* y = grid1.y;
    double* z = grid1.z;
    int* jimage = image.jimage;   // common /image/
    int* kimage = image.kimage;
    int* limage = image.limage;

    // read input data and initialize
    initia();

    // do 11 mesh=1,nmesh
    for (int mesh = 1; mesh <= conec1.nmesh; ++mesh) {

        int jdg = grid2.mjmax[mesh] + 1;
        int kdg = grid2.mkmax[mesh] + 1;
        int ldg = grid2.mlmax[mesh] + 1;

        // read in grid
        if (igrdtyp.ip3dgrd == 0) {
            // (cfl3d format)
            rgrid(mesh, x, y, z, jdg, kdg, ldg);
        } else {
            // (plot3d/tlns3d format)
            rp3d(mesh, x, y, z, jdg, kdg, ldg, igrdtyp.ialph, conec1.nmesh);
        }

        int npnts = jdg * kdg * ldg;

        if (npnts > mlen) {
            // format 601
            std::fprintf(funit(6),
                "%c%.*s\n\n%10s%s %s %s\n%30s%s%3d%5s%s%3d\n\n%.*s\n",
                '0', 130, "*******************************************************"
                          "*******************************************************"
                          "********************",
                "", "failure in setup due to number", "of total grid points exceeding allocate",
                "storage",
                "", "npnts =", npnts, "", "mlen =", mlen,
                130, "*******************************************************"
                     "*******************************************************"
                     "********************");
            std::fprintf(stderr, "setup\n");
            std::exit(1);
        }

        // check for branch cuts in grid
        if (mesh == 1) {
            std::fprintf(funit(6), "\n");                              // write(6,*)
            // write(6,*) '  *** beginning checks to identify branch',' cuts ***'
            std::fprintf(funit(6), "   *** beginning checks to identify branch cuts ***\n");
        }
        std::fprintf(funit(6), "\n");                                  // write(6,*)
        branch(mesh, jdg, kdg, ldg, jimage, kimage, limage, x, y, z);

        // transform grid to cell centered grid.
        cellcen(x, y, z, jdg, kdg, ldg);

        // write cell centered grid
        wcelgrd(mesh, x, y, z, jimage, kimage, limage, jdg, kdg, ldg);

    } // 11 continue

    frewind(1);   // rewind 1

    // recompute npnts (unused result, kept for fidelity)
    int npnts = 0;
    for (int m = 1; m <= conec1.nmesh; ++m) {
        npnts = (grid2.mjmax[m] + 1) * (grid2.mkmax[m] + 1) * (grid2.mlmax[m] + 1) + npnts;
    }
    (void)npnts;

    return;
}

// ===========================================================================
//  subroutine initia  — read connection and grid parameters
// ===========================================================================
void initia() {
    // 20a4 title -> 80 characters
    // Unit 5 is Fortran stdin; main() already consumed the first 5 header lines
    // from the same stdin stream, so continue reading from stdin here.
    tc::TextReader rd(stdin);      // unit 5 = stdin

    // read(5,10) (title(i),i=1,20)  ; write(6,11)
    std::string titleline = rd.line();
    // format 11: (/,1x,20a4)  -> blank line, then " " + first 80 chars of title
    {
        char t80[81];
        for (int i = 0; i < 80; ++i) t80[i] = ' ';
        t80[80] = '\0';
        for (int i = 0; i < 80 && i < (int)titleline.size(); ++i) t80[i] = titleline[i];
        std::fprintf(funit(6), "\n %s\n", t80);
    }

    // read(5,10) (skip header) ; write(6,21)
    rd.line();
    // format 21: (5x,5hnmesh,4x,6hiplt3d)
    std::fprintf(funit(6), "     nmesh    iplt3d\n");

    // read(5,*) nmesh,iplt3d
    {
        int tmp[2];
        rd.ints(tmp, 2);
        conec1.nmesh   = tmp[0];
        diagno.iplt3d  = tmp[1];
    }
    // write(6,20) nmesh,iplt3d      format 20: (8i10)
    std::fprintf(funit(6), "%10d%10d\n", conec1.nmesh, diagno.iplt3d);

    igrdtyp.ip3dgrd = 0;
    if (conec1.nmesh < 0) {
        igrdtyp.ip3dgrd = 1;
        conec1.nmesh = std::abs(conec1.nmesh);
    }

    if (conec1.nmesh > mdim) {
        // format 5: (5x,5hmdim=,i2,2x,6hnmesh=,i2,16hincrease mdim...,25hmust have mdim .ge. nmesh)
        std::fprintf(funit(6), "     mdim=%2d  nmesh=%2dincrease mdim...must have mdim .ge. nmesh\n",
                     mdim, conec1.nmesh);
        std::exit(0);   // stop
    }

    // read(5,10) (skip) ; write(6,101)
    rd.line();
    // format 101: (6x,4hjmax,6x,4hkmax,6x,4hlmax,5x,5hnoutr,5x,5hnhole,4x,6hnobtyp)
    std::fprintf(funit(6), "      jmax      kmax      lmax     noutr     nhole    nobtyp\n");

    int nout    = 0;
    int nob     = 0;
    int ihotest = 0;
    // do 100 n=1,nmesh
    for (int n = 1; n <= conec1.nmesh; ++n) {
        // read(5,91) mjmax,mkmax,mlmax,noutr,nhole,nobtyp   format 91: (7i10,l10)
        int v[6];
        {
            std::string l = rd.line();
            std::stringstream ss(l);
            for (int i = 0; i < 6; ++i) { std::string t; ss >> t; v[i] = std::atoi(t.c_str()); }
        }
        grid2.mjmax[n]  = v[0];
        grid2.mkmax[n]  = v[1];
        grid2.mlmax[n]  = v[2];
        conec1.noutr[n] = v[3];
        conec1.nhole[n] = v[4];
        bound2.nobtyp[n] = v[5];
        // write(6,91) mjmax,mkmax,mlmax,noutr,nhole,nobtyp
        // 7i10,l10 : only 6 ints written, then a logical field is unused (no value in list) -> not emitted
        std::fprintf(funit(6), "%10d%10d%10d%10d%10d%10d\n",
                     grid2.mjmax[n], grid2.mkmax[n], grid2.mlmax[n],
                     conec1.noutr[n], conec1.nhole[n], bound2.nobtyp[n]);
        grid2.mjmax[n] = grid2.mjmax[n] - 1;
        grid2.mkmax[n] = grid2.mkmax[n] - 1;
        grid2.mlmax[n] = grid2.mlmax[n] - 1;
        nout = nout + std::min(conec1.noutr[n], 1);
        nob  = nob + bound2.nobtyp[n];
        if (conec1.nhole[n] > 0) ihotest = ihotest + conec1.nhole[n];
        if (conec1.nhole[n] > mhldim || conec1.noutr[n] > mhldim) {
            // format 92
            std::fprintf(funit(6),
                "  stopping...number of holes/outer boundaries exceeds parameter mhldim\n");
            std::exit(0);
        }
    } // 100

    if (nob > ipmax) {
        // format 4051
        std::fprintf(funit(6),
            "  stopping...number of planes defining outer boundary points to be interpolated exceeds ipmax\n");
        std::exit(0);
    }

    // read(5,10) (skip) ; write(6,151)
    rd.line();
    // format 151: (6x,4hjsrs,...,4hlsre)
    std::fprintf(funit(6), "      jsrs      jsre      ksrs      ksre      lsrs      lsre\n");
    // do 150 n=1,nmesh
    for (int n = 1; n <= conec1.nmesh; ++n) {
        int s[6];
        rd.ints(s, 6);   // read(5,*) jsrs,jsre,ksrs,ksre,lsrs,lsre
        grdlim.jsrs[n] = s[0];
        grdlim.jsre[n] = s[1];
        grdlim.ksrs[n] = s[2];
        grdlim.ksre[n] = s[3];
        grdlim.lsrs[n] = s[4];
        grdlim.lsre[n] = s[5];

        if (grdlim.jsrs[n] <= 0)                                            grdlim.jsrs[n] = 1;
        if (grdlim.jsre[n] <= 0 || grdlim.jsre[n] > grid2.mjmax[n])         grdlim.jsre[n] = grid2.mjmax[n];
        if (grdlim.ksrs[n] <= 0)                                            grdlim.ksrs[n] = 1;
        if (grdlim.ksre[n] <= 0 || grdlim.ksre[n] > grid2.mkmax[n])         grdlim.ksre[n] = grid2.mkmax[n];
        if (grdlim.lsrs[n] <= 0)                                            grdlim.lsrs[n] = 1;
        if (grdlim.lsre[n] <= 0 || grdlim.lsre[n] > grid2.mlmax[n])         grdlim.lsre[n] = grid2.mlmax[n];

        // write(6,20) jsrs,jsre,ksrs,ksre,lsrs,lsre    format 20: (8i10)
        std::fprintf(funit(6), "%10d%10d%10d%10d%10d%10d\n",
                     grdlim.jsrs[n], grdlim.jsre[n], grdlim.ksrs[n],
                     grdlim.ksre[n], grdlim.lsrs[n], grdlim.lsre[n]);
    } // 150

    // do 200 n=1,nmesh ; do 200 m=1,nmesh : ihbtyp(m,n)=0
    for (int n = 1; n <= conec1.nmesh; ++n)
        for (int m = 1; m <= conec1.nmesh; ++m) bound1.IHBTYP(m, n) = 0;

    if (ihotest > 0) {
        int nt    = 0;
        int ntest = 0;
        int jp2   = 0;
        int ihtot = 0;
        // do 600 i=1,ihotest
        for (int i = 1; i <= ihotest; ++i) {

            // do 350 n=1,nmesh : icount(n)=0
            for (int n = 1; n <= conec1.nmesh; ++n) iseqc.icount[n] = 0;

            rd.line();                                     // read(5,10)
            // write(6,399) format 399: (1x,16hhole definition:)
            std::fprintf(funit(6), " hole definition:\n");

            rd.line();                                     // read(5,10)
            // write(6,401)
            std::fprintf(funit(6),
                "      mesh  hole no. ihbtyp(1) ihbtyp(2) ihbtyp(3) ihbtyp(4) ihbtyp(5)\n");

            // read(5,*) nn,nhle,(ihbtyp(nn,m),m=1,nmesh)
            int nn, nhle;
            {
                int need = 2 + conec1.nmesh;
                std::vector<int> vv(need);
                rd.ints(vv.data(), need);
                nn   = vv[0];
                nhle = vv[1];
                for (int m = 1; m <= conec1.nmesh; ++m) bound1.IHBTYP(nn, m) = vv[2 + (m - 1)];
            }

            if (nn < ntest) {
                // format 340
                std::fprintf(funit(6),
                    "\n\n stopping ... incorrect order for IHBTYP lines;\n"
                    " lines must be in ascending order according to mesh number\n"
                    " current mesh number being read =%3d ; previous mesh number =%3d\n",
                    nn, ntest);
                std::exit(0);
            } else {
                ntest = nn;
                bound1.ihplt[i] = 0;
                if (nhle < 0) {
                    bound1.ihplt[i] = 1;
                    nhle = std::abs(nhle);
                }
            }

            // write(6,20) nn,nhle,(ihbtyp(nn,m),m=1,nmesh)
            {
                std::string out;
                char b[16];
                std::snprintf(b, sizeof(b), "%10d", nn);   out += b;
                std::snprintf(b, sizeof(b), "%10d", nhle); out += b;
                for (int m = 1; m <= conec1.nmesh; ++m) {
                    std::snprintf(b, sizeof(b), "%10d", bound1.IHBTYP(nn, m)); out += b;
                }
                std::fprintf(funit(6), "%s\n", out.c_str());
            }

            int ih = 0;
            // do 360 m=1,nmesh : ih=ih+ihbtyp(nn,m)
            for (int m = 1; m <= conec1.nmesh; ++m) ih = ih + bound1.IHBTYP(nn, m);

            int jp1 = jp2 + 1;
            jp2 = jp2 + ih;
            bound1.ip1[i] = jp1;
            bound1.ip2[i] = jp2;

            ihtot = ihtot + ih;
            if (ihtot > ipmax) {
                // format 6051
                std::fprintf(funit(6),
                    "  stopping...number of planes defining holes exceeds ipmax\n");
                std::exit(0);
            }

            iseqc.iseq[1] = nt;
            // do 400 m=2,nmesh : iseq(m)=ihbtyp(nn,m-1)+iseq(m-1)
            for (int m = 2; m <= conec1.nmesh; ++m)
                iseqc.iseq[m] = bound1.IHBTYP(nn, m - 1) + iseqc.iseq[m - 1];

            rd.line();                                     // read(5,10)
            // write(6,605)
            std::fprintf(funit(6),
                "     imesh       jh1       jh2       kh1       kh2       lh1       lh2\n");

            // do 500 n=1,ih
            for (int nloop = 1; nloop <= ih; ++nloop) {
                int mm, j1, j2, k1, k2, l1, l2;
                {
                    int t7[7];
                    rd.ints(t7, 7);   // read(5,*) mm,j1,j2,k1,k2,l1,l2
                    mm = t7[0]; j1 = t7[1]; j2 = t7[2]; k1 = t7[3];
                    k2 = t7[4]; l1 = t7[5]; l2 = t7[6];
                }

                if (iseqc.icount[mm] > bound1.IHBTYP(nn, mm)) {
                    // format 501
                    std::fprintf(funit(6),
                        "\n\n stopping ... too many holes defined by mesh%3d for mesh =%3d\n",
                        mm, nn);
                    std::exit(0);
                }

                if (bound1.IHBTYP(nn, mm) == 0) {
                    // format 222 then 502
                    std::fprintf(funit(6),
                        "\n\n stopping ... inconsistent input for mesh m =%3d\n", nn);
                    std::fprintf(funit(6),
                        "\n attempted entry of JH1, JH2... values for mesh n =%3d but IHBTYP(m,n) = 0\n",
                        mm);
                    std::exit(0);
                } else {
                    iseqc.icount[mm] = iseqc.icount[mm] + 1;
                }

                if (j1 == 0) j1 = 1;
                if (j2 == 0) j2 = grid2.mjmax[mm] + 1;
                if (k1 == 0) k1 = 1;
                if (k2 == 0) k2 = grid2.mkmax[mm] + 1;
                if (l1 == 0) l1 = 1;
                if (l2 == 0) l2 = grid2.mlmax[mm] + 1;

                // check for inappropriate indices on components of hole surface
                if (std::abs(j1) > grid2.mjmax[mm] + 1) {
                    std::fprintf(funit(6), "   stopping: j1 = %4d jmax = %4d\n", j1, grid2.mjmax[mm] + 1);
                    std::exit(0);
                }
                if (std::abs(j2) > grid2.mjmax[mm] + 1) {
                    std::fprintf(funit(6), "   stopping: j2 = %4d jmax = %4d\n", j2, grid2.mjmax[mm] + 1);
                    std::exit(0);
                }
                if (std::abs(k1) > grid2.mkmax[mm] + 1) {
                    std::fprintf(funit(6), "   stopping: k1 = %4d kmax = %4d\n", k1, grid2.mkmax[mm] + 1);
                    std::exit(0);
                }
                if (std::abs(k2) > grid2.mkmax[mm] + 1) {
                    std::fprintf(funit(6), "   stopping: k2 = %4d kmax = %4d\n", k2, grid2.mkmax[mm] + 1);
                    std::exit(0);
                }
                if (std::abs(l1) > grid2.mlmax[mm] + 1) {
                    std::fprintf(funit(6), "   stopping: l1 = %4d lmax = %4d\n", l1, grid2.mlmax[mm] + 1);
                    std::exit(0);
                }
                if (std::abs(l2) > grid2.mlmax[mm] + 1) {
                    std::fprintf(funit(6), "   stopping: l2 = %4d lmax = %4d\n", l2, grid2.mlmax[mm] + 1);
                    std::exit(0);
                }

                // write(6,20) mm,j1,j2,k1,k2,l1,l2
                std::fprintf(funit(6), "%10d%10d%10d%10d%10d%10d%10d\n",
                             mm, j1, j2, k1, k2, l1, l2);
                iseqc.iseq[mm] = iseqc.iseq[mm] + 1;

                bound1.mh [iseqc.iseq[mm]] = mm;
                bound1.jh1[iseqc.iseq[mm]] = j1;
                bound1.jh2[iseqc.iseq[mm]] = j2;
                bound1.kh1[iseqc.iseq[mm]] = k1;
                bound1.kh2[iseqc.iseq[mm]] = k2;
                bound1.lh1[iseqc.iseq[mm]] = l1;
                bound1.lh2[iseqc.iseq[mm]] = l2;
            } // 500

            rd.line();                                     // read(5,10)
            rd.line();                                     // read(5,10)
            // write(6,129)
            std::fprintf(funit(6), "      search list for associated fringe points:\n");
            // write(6,131)
            std::fprintf(funit(6),
                "     nlist    lst(1)    lst(2)    lst(3)    lst(4)    ...  nlist\n");

            // read(5,*)mhole(nn,nhle),(lhole(nn,nhle,ll),ll=1,mhole(nn,nhle))
            {
                std::string l = rd.line();
                std::stringstream ss(l);
                int mval; std::string t;
                ss >> t; mval = std::atoi(t.c_str());
                conec1.MHOLE(nn, nhle) = mval;
                for (int ll = 1; ll <= mval; ++ll) {
                    ss >> t; conec1.LHOLE(nn, nhle, ll) = std::atoi(t.c_str());
                }
                // write(6,20) mhole,(lhole...)
                std::string out; char b[16];
                std::snprintf(b, sizeof(b), "%10d", mval); out += b;
                for (int ll = 1; ll <= mval; ++ll) {
                    std::snprintf(b, sizeof(b), "%10d", conec1.LHOLE(nn, nhle, ll)); out += b;
                }
                std::fprintf(funit(6), "%s\n", out.c_str());
            }
            // do 1502 ll=1,mhole(nn,nhle)
            for (int ll = 1; ll <= conec1.MHOLE(nn, nhle); ++ll) {
                if (conec1.LHOLE(nn, nhle, ll) == nn) {
                    // format 123
                    std::fprintf(funit(6),
                        "    stopping...mesh %3d not allowed in search list for mesh %3d fringe points\n",
                        nn, nn);
                    std::exit(0);
                }
            } // 1502

            rd.line();                                     // read(5,10)
            rd.line();                                     // read(5,10)
            // write(6,329)
            std::fprintf(funit(6), "      problem point options:\n");
            // write(6,331)
            std::fprintf(funit(6), "     iorph     ihole\n");
            // read(5,*) iorphh(nn,nhle),iholeh(nn,nhle)
            {
                int t2[2];
                rd.ints(t2, 2);
                snafu.IORPHH(nn, nhle) = t2[0];
                snafu.IHOLEH(nn, nhle) = t2[1];
            }
            // write(6,20) iorphh,iholeh
            std::fprintf(funit(6), "%10d%10d\n",
                         snafu.IORPHH(nn, nhle), snafu.IHOLEH(nn, nhle));

            nt = nt + ih;
            // do 580 n=1,nmesh
            for (int n = 1; n <= conec1.nmesh; ++n) {
                if (iseqc.icount[n] != bound1.IHBTYP(nn, n)) {
                    // format 222 then 581
                    std::fprintf(funit(6),
                        "\n\n stopping ... inconsistent input for mesh m =%3d\n", nn);
                    std::fprintf(funit(6),
                        "\n number of entries of JH1, JH2... values for mesh n =%3d is %3d but IHBTYP(m,n) =%3d\n",
                        n, iseqc.icount[n], bound1.IHBTYP(nn, n));
                    std::exit(0);
                }
            } // 580
        } // 600
    } else {
        rd.line();  // read(5,10)
        std::fprintf(funit(6), " hole definition:\n");                      // 399
        rd.line();
        std::fprintf(funit(6),
            "      mesh  hole no. ihbtyp(1) ihbtyp(2) ihbtyp(3) ihbtyp(4) ihbtyp(5)\n");  // 401
        rd.line();
        std::fprintf(funit(6),
            "     imesh       jh1       jh2       kh1       kh2       lh1       lh2\n");   // 605
        rd.line();
        std::fprintf(funit(6), "      search list for associated fringe points:\n");      // 129
        rd.line();
        std::fprintf(funit(6),
            "     nlist    lst(1)    lst(2)    lst(3)    lst(4)    ...  nlist\n");         // 131
        rd.line();
        std::fprintf(funit(6), "      problem point options:\n");                         // 329
        rd.line();
        std::fprintf(funit(6), "     iorph     ihole\n");                                 // 331
    }

    if (nob > 0) {
        iseqc.iseq[1] = 0;
        // do 650 m=2,nmesh : iseq(m)=nobtyp(m-1)+iseq(m-1)
        for (int m = 2; m <= conec1.nmesh; ++m)
            iseqc.iseq[m] = bound2.nobtyp[m - 1] + iseqc.iseq[m - 1];

        // do 700 nn=1,nmesh
        for (int nn = 1; nn <= conec1.nmesh; ++nn) {
            if (bound2.nobtyp[nn] > 0) {
                rd.line();                                 // read(5,10)
                rd.line();                                 // read(5,10)
                // write(6,699)
                std::fprintf(funit(6),
                    " definition of outer boundary points (cell center)to be interpolated:\n");
                // write(6,701)
                std::fprintf(funit(6),
                    "      mesh       jo1       jo2       ko1       ko2       lo1       lo2\n");

                // do 7001 n=1,nobtyp(nn)
                for (int n = 1; n <= bound2.nobtyp[nn]; ++n) {
                    int mm, j1, j2, k1, k2, l1, l2;
                    {
                        int t7[7];
                        rd.ints(t7, 7);   // read(5,*) mm,j1,j2,k1,k2,l1,l2
                        mm = t7[0]; j1 = t7[1]; j2 = t7[2]; k1 = t7[3];
                        k2 = t7[4]; l1 = t7[5]; l2 = t7[6];
                    }
                    if (bound2.nobtyp[mm] == 0) {
                        std::fprintf(funit(6),
                            "\n\n stopping ... inconsistent input for mesh m =%3d\n", mm);
                        std::fprintf(funit(6),
                            "\n attempted entry of JO1, JO2... values but NOBTYP(m) = 0\n");
                        std::exit(0);
                    }

                    if (j1 == 0) j1 = 1;
                    if (j2 == 0) j2 = grid2.mjmax[mm];
                    if (k1 == 0) k1 = 1;
                    if (k2 == 0) k2 = grid2.mkmax[mm];
                    if (l1 == 0) l1 = 1;
                    if (l2 == 0) l2 = grid2.mlmax[mm];

                    if (j1 < 0) j2 = 0;
                    if (k1 < 0) k2 = 0;
                    if (l1 < 0) l2 = 0;

                    // write(6,20) mm,j1,j2,k1,k2,l1,l2
                    std::fprintf(funit(6), "%10d%10d%10d%10d%10d%10d%10d\n",
                                 mm, j1, j2, k1, k2, l1, l2);
                    iseqc.iseq[mm] = iseqc.iseq[mm] + 1;

                    bound2.jo1[iseqc.iseq[mm]] = j1;
                    bound2.jo2[iseqc.iseq[mm]] = j2;
                    bound2.ko1[iseqc.iseq[mm]] = k1;
                    bound2.ko2[iseqc.iseq[mm]] = k2;
                    bound2.lo1[iseqc.iseq[mm]] = l1;
                    bound2.lo2[iseqc.iseq[mm]] = l2;
                } // 7001

                rd.line();                                 // read(5,10)
                rd.line();                                 // read(5,10)
                // write(6,229)
                std::fprintf(funit(6), "    search list for outer boundary points:\n");
                // write(6,231)
                std::fprintf(funit(6),
                    "     nlist    lst(1)    lst(2)    lst(3)    lst(4)    ...  nlist\n");

                int n = 1;

                // read(5,*)moutr(nn,n),(loutr(nn,n,ll),ll=1,moutr(nn,n))
                {
                    std::string l = rd.line();
                    std::stringstream ss(l);
                    int mval; std::string t;
                    ss >> t; mval = std::atoi(t.c_str());
                    conec1.MOUTR(nn, n) = mval;
                    for (int ll = 1; ll <= mval; ++ll) {
                        ss >> t; conec1.LOUTR(nn, n, ll) = std::atoi(t.c_str());
                    }
                    std::string out; char b[16];
                    std::snprintf(b, sizeof(b), "%10d", mval); out += b;
                    for (int ll = 1; ll <= mval; ++ll) {
                        std::snprintf(b, sizeof(b), "%10d", conec1.LOUTR(nn, n, ll)); out += b;
                    }
                    std::fprintf(funit(6), "%s\n", out.c_str());
                }
                // do 2502 ll=1,moutr(nn,n)
                for (int ll = 1; ll <= conec1.MOUTR(nn, n); ++ll) {
                    if (conec1.LOUTR(nn, n, ll) == nn) {
                        // format 223
                        std::fprintf(funit(6),
                            "    stopping...mesh %3d not allowed in search list for mesh %3d boundary points\n",
                            nn, nn);
                        std::exit(0);
                    }
                } // 2502

                rd.line();                                 // read(5,10)
                rd.line();                                 // read(5,10)
                // write(6,329) ; write(6,331)
                std::fprintf(funit(6), "      problem point options:\n");
                std::fprintf(funit(6), "     iorph     ihole\n");
                // read(5,*) iorpho(nn,n),iholeo(nn,n)
                {
                    int t2[2];
                    rd.ints(t2, 2);
                    snafu.IORPHO(nn, n) = t2[0];
                    snafu.IHOLEO(nn, n) = t2[1];
                }
                // write(6,20) iorpho,iholeo
                std::fprintf(funit(6), "%10d%10d\n",
                             snafu.IORPHO(nn, n), snafu.IHOLEO(nn, n));
            }
        } // 700
    } // else: nob<=0 (commented-out reads in Fortran, nothing done)

    int mb = 0;
    // do 750 m=1,nmesh
    for (int m = 1; m <= conec1.nmesh; ++m) {
        bound2.iobtyp[m] = 0;
        if (bound2.nobtyp[m] != 0) {
            bound2.iobtyp[m] = mb + 1;
            // do 725 nb=1,nobtyp(m)
            for (int nb = 1; nb <= bound2.nobtyp[m]; ++nb) {
                mb = mb + 1;
                if (bound2.jo1[mb] > grid2.mjmax[m] + 2 || bound2.jo2[mb] > grid2.mjmax[m] + 2 ||
                    bound2.ko1[mb] > grid2.mkmax[m] + 2 || bound2.ko2[mb] > grid2.mkmax[m] + 2 ||
                    bound2.lo1[mb] > grid2.mlmax[m] + 2 || bound2.lo2[mb] > grid2.mlmax[m] + 2 ||
                    bound2.jo1[mb] < -1                 || bound2.jo2[mb] < -1                 ||
                    bound2.ko1[mb] < -1                 || bound2.ko2[mb] < -1                 ||
                    bound2.lo1[mb] < -1                 || bound2.lo2[mb] < -1) {
                    // format 755
                    std::fprintf(funit(6),
                        "\n\n     value of outer boundary index exceeds constrained values on mesh %3d\n", m);
                    std::exit(0);
                }
            } // 725
        }
    } // 750

    // initialize iblank to 1 for each grid
    // do 1000 m=1,nmesh
    for (int m = 1; m <= conec1.nmesh; ++m) {
        int jd = grid2.mjmax[m];
        int kd = grid2.mkmax[m];
        int ld = grid2.mlmax[m];
        putibl(m, book3.iblank, jd, kd, ld);
    } // 1000

    // initialize the interpolation arrays for each grid
    int nwr = 0;
    // do 1100 m=1,nmesh
    for (int m = 1; m <= conec1.nmesh; ++m) {
        putint(nwr, m, book2.ji, book2.ki, book2.li,
               book2.jbpt, book2.kbpt, book2.lbpt,
               book2.dxint, book2.dyint, book2.dzint);
    } // 1100

    // consistency check
    int istop = 0;
    // do 1200 m=1,nmesh
    for (int m = 1; m <= conec1.nmesh; ++m) {
        int ntotal = (grid2.mjmax[m] + 1) * (grid2.mkmax[m] + 1) * (grid2.mlmax[m] + 1);
        if (ntotal > mlen) {
            // format 1205
            std::fprintf(funit(6),
                "0          fatal error in mesh = %3d\n"
                "          the total number of points in mesh is %6d max dimension is %6d\n",
                m, ntotal, mlen);
            istop = 1;
        }
    } // 1200

    // check each dimension
    if (conec1.nmesh > mdim) {
        // format 1305
        std::fprintf(funit(6),
            "0%.*s\n\n     failure in initia,number of grids input,nmesh=%3d   is greater than dimension, mdim =%3d\n\n%.*s\n",
            130, "*******************************************************"
                 "*******************************************************"
                 "********************",
            conec1.nmesh, mdim,
            130, "*******************************************************"
                 "*******************************************************"
                 "********************");
        istop = 1;
    } else {
        // do 1300 m=1,nmesh
        for (int m = 1; m <= conec1.nmesh; ++m) {
            if (grid2.mjmax[m] + 1 > jdim || grid2.mkmax[m] + 1 > kdim ||
                grid2.mlmax[m] + 1 > ldim) {
                // format 1310
                std::fprintf(funit(6),
                    "0%.*s\n\n     input failure in initia exceeds maximum dimension\n"
                    "                         mesh =%5d   jmax,jdim%5d%5d     kmax,kdim%5d%5d     lmax,ldim%5d%5d\n\n%.*s\n",
                    130, "*******************************************************"
                         "*******************************************************"
                         "********************",
                    m, grid2.mjmax[m] + 1, jdim, grid2.mkmax[m] + 1, kdim,
                    grid2.mlmax[m] + 1, ldim,
                    130, "*******************************************************"
                         "*******************************************************"
                         "********************");
                istop = 1;
            }
        } // 1300
    }

    if (istop == 1) { std::fprintf(stderr, "initia\n"); std::exit(1); }
    return;
}

// ===========================================================================
//  subroutine rgrid  — read a grid from unit 10 (CFL3D format) -> temp_grd.m
// ===========================================================================
void rgrid(int m, double* x, double* y, double* z, int jd, int kd, int ld) {
    char titl[21];

    // read(10) jmax,kmax,lmax
    int hdr[3];
    uread(10, hdr, 3 * sizeof(int));
    int jmax = hdr[0], kmax = hdr[1], lmax = hdr[2];

    if (jd != jmax || kd != kmax || ld != lmax) {
        // format 606
        std::fprintf(funit(6),
            "0%.*s\n\n     grid input failed due to inconsistent input dimensions:\n"
            "                         mesh number %5d\n"
            "                         jmax, tape jmax%5d%5d\n"
            "                         kmax, tape kmax%5d%5d\n"
            "                         lmax, tape lmax%5d%5d\n\n%.*s\n",
            130, "*******************************************************"
                 "*******************************************************"
                 "********************",
            m, jd, jmax, kd, kmax, ld, lmax,
            130, "*******************************************************"
                 "*******************************************************"
                 "********************");
        std::fprintf(stderr, "rgrid\n");
        std::exit(1);
    }

    // read(10) x,y,z   -- one record, 3*jd*kd*ld doubles, contiguous into x,y,z
    long n = (long)jd * kd * ld;
    {
        // The record holds x block, then y block, then z block, contiguously.
        // x,y,z are separate 1-based buffers, so read into a temp then split,
        // OR read each block directly since the Fortran record is one unit.
        // gfortran urec_read reads the whole record's leading bytes into buf.
        // We need a single contiguous buffer of 3n doubles.
        std::vector<double> tmp((size_t)(3 * n));
        uread(10, tmp.data(), (size_t)(3 * n) * sizeof(double));
        for (long i = 0; i < n; ++i) {
            x[1 + i] = tmp[(size_t)i];
            y[1 + i] = tmp[(size_t)(n + i)];
            z[1 + i] = tmp[(size_t)(2 * n + i)];
        }
    }

    int len = build_scratch_name(titl, "temp_grd.", m);
    std::string fname(titl, len);
    fopen_unit(30, fname.c_str(), "wb+");   // open unit 30 unformatted

    // write(iunit) x,y,z  -- one record with the three blocks contiguous
    {
        tc::RecordWriter rw;
        rw.put_bytes(&x[1], (size_t)n * sizeof(double));
        rw.put_bytes(&y[1], (size_t)n * sizeof(double));
        rw.put_bytes(&z[1], (size_t)n * sizeof(double));
        rw.flush(funit(30));
    }
    frewind(30);

    return;
}

// ===========================================================================
//  subroutine rp3d  — read a plot3d/tlns3d grid from unit 10 -> temp_grd.m
// ===========================================================================
void rp3d(int m, double* x, double* y, double* z, int jd, int kd, int ld,
          int& ialph, int nmesh) {
    char titl[21];

    // ialph - flag to interpret input geometry
    ialph = 1;

    if (m == 1) {
        // read(10) ngrid
        int ngrid;
        uread(10, &ngrid, sizeof(int));
        if (ngrid != nmesh) {
            // format 605
            std::fprintf(funit(6),
                " stopping, inconsistent input: nmesh, ngrid = %3d%3d\n", nmesh, ngrid);
            std::fprintf(stderr, "rp3d\n");
            std::exit(1);
        }
        // read(10) (ltest(mm),jtest(mm),ktest(mm),mm=1,nmesh)
        {
            std::vector<int> buf((size_t)(3 * nmesh));
            uread(10, buf.data(), (size_t)(3 * nmesh) * sizeof(int));
            for (int mm = 1; mm <= nmesh; ++mm) {
                dimenc.ltest[mm] = buf[(size_t)(3 * (mm - 1) + 0)];
                dimenc.jtest[mm] = buf[(size_t)(3 * (mm - 1) + 1)];
                dimenc.ktest[mm] = buf[(size_t)(3 * (mm - 1) + 2)];
            }
        }
    }

    if (jd != dimenc.jtest[m] || kd != dimenc.ktest[m] || ld != dimenc.ltest[m]) {
        // format 606
        std::fprintf(funit(6),
            "0%.*s\n\n     grid input failed due to inconsistent input dimensions:\n"
            "                         mesh number %5d\n"
            "                         jmax, tape jmax%5d%5d\n"
            "                         kmax, tape kmax%5d%5d\n"
            "                         lmax, tape lmax%5d%5d\n\n%.*s\n",
            130, "*******************************************************"
                 "*******************************************************"
                 "********************",
            m, jd, dimenc.jtest[m], kd, dimenc.ktest[m], ld, dimenc.ltest[m],
            130, "*******************************************************"
                 "*******************************************************"
                 "********************");
        std::fprintf(stderr, "rp3d\n");
        std::exit(1);
    }

    View3<double> X(x, jd, kd, ld), Y(y, jd, kd, ld), Z(z, jd, kd, ld);
    long n = (long)jd * kd * ld;

    if (ialph == 0) {
        // read(10) (((x(j,k,l),l),j),k), (((y...)), (((z...))  — l fastest
        std::vector<double> tmp((size_t)(3 * n));
        uread(10, tmp.data(), (size_t)(3 * n) * sizeof(double));
        long p = 0;
        for (int k = 1; k <= kd; ++k) for (int j = 1; j <= jd; ++j) for (int l = 1; l <= ld; ++l)
            X(j, k, l) = tmp[(size_t)(p++)];
        for (int k = 1; k <= kd; ++k) for (int j = 1; j <= jd; ++j) for (int l = 1; l <= ld; ++l)
            Y(j, k, l) = tmp[(size_t)(p++)];
        for (int k = 1; k <= kd; ++k) for (int j = 1; j <= jd; ++j) for (int l = 1; l <= ld; ++l)
            Z(j, k, l) = tmp[(size_t)(p++)];
    } else {
        // read(10) x-block, z-block, y-block (l fastest), then negate y
        std::vector<double> tmp((size_t)(3 * n));
        uread(10, tmp.data(), (size_t)(3 * n) * sizeof(double));
        long p = 0;
        for (int k = 1; k <= kd; ++k) for (int j = 1; j <= jd; ++j) for (int l = 1; l <= ld; ++l)
            X(j, k, l) = tmp[(size_t)(p++)];
        for (int k = 1; k <= kd; ++k) for (int j = 1; j <= jd; ++j) for (int l = 1; l <= ld; ++l)
            Z(j, k, l) = tmp[(size_t)(p++)];
        for (int k = 1; k <= kd; ++k) for (int j = 1; j <= jd; ++j) for (int l = 1; l <= ld; ++l)
            Y(j, k, l) = tmp[(size_t)(p++)];
        // do 10 l ; j ; k : y(j,k,l) = -y(j,k,l)
        for (int l = 1; l <= ld; ++l)
            for (int j = 1; j <= jd; ++j)
                for (int k = 1; k <= kd; ++k)
                    Y(j, k, l) = -Y(j, k, l);
    }

    int len = build_scratch_name(titl, "temp_grd.", m);
    std::string fname(titl, len);
    fopen_unit(30, fname.c_str(), "wb+");

    // write(iunit) x,y,z  -- whole arrays (column-major), one record
    {
        tc::RecordWriter rw;
        rw.put_bytes(&x[1], (size_t)n * sizeof(double));
        rw.put_bytes(&y[1], (size_t)n * sizeof(double));
        rw.put_bytes(&z[1], (size_t)n * sizeof(double));
        rw.flush(funit(30));
    }
    frewind(30);

    return;
}

// ===========================================================================
//  subroutine cellcen  — find the cell centers of the grid
// ===========================================================================
void cellcen(double* x, double* y, double* z, int jd, int kd, int ld) {
    View3<double> X(x, jd, kd, ld), Y(y, jd, kd, ld), Z(z, jd, kd, ld);

    int jmax = jd;
    int kmax = kd;
    int lmax = ld;

    // do 10 l=1,lmax-1 ; do 20 k=1,kmax-1 ; do 30 j=1,jmax-1
    for (int l = 1; l <= lmax - 1; ++l)
        for (int k = 1; k <= kmax - 1; ++k)
            for (int j = 1; j <= jmax - 1; ++j) {
                temp1c.X1(j, k, l) = (X(j, k, l) + X(j + 1, k, l) + X(j, k + 1, l) + X(j + 1, k + 1, l) +
                                      X(j, k, l + 1) + X(j, k + 1, l + 1) + X(j + 1, k, l + 1) +
                                      X(j + 1, k + 1, l + 1)) / 8.0;
                temp1c.Y1(j, k, l) = (Y(j, k, l) + Y(j + 1, k, l) + Y(j, k + 1, l) + Y(j + 1, k + 1, l) +
                                      Y(j, k, l + 1) + Y(j, k + 1, l + 1) + Y(j + 1, k, l + 1) +
                                      Y(j + 1, k + 1, l + 1)) / 8.0;
                temp1c.Z1(j, k, l) = (Z(j, k, l) + Z(j + 1, k, l) + Z(j, k + 1, l) + Z(j + 1, k + 1, l) +
                                      Z(j, k, l + 1) + Z(j, k + 1, l + 1) + Z(j + 1, k, l + 1) +
                                      Z(j + 1, k + 1, l + 1)) / 8.0;
            }

    // do 140 l=1,lmax-1 ; do 150 k=1,kmax-1 ; do 160 j=1,jmax-1
    for (int l = 1; l <= lmax - 1; ++l)
        for (int k = 1; k <= kmax - 1; ++k)
            for (int j = 1; j <= jmax - 1; ++j) {
                X(j, k, l) = temp1c.X1(j, k, l);
                Y(j, k, l) = temp1c.Y1(j, k, l);
                Z(j, k, l) = temp1c.Z1(j, k, l);
            }

    return;
}

// ===========================================================================
//  subroutine wcelgrd  — write cell-center coords for mesh m into temp_cen.m
// ===========================================================================
void wcelgrd(int m, double* x, double* y, double* z,
             int* jimage, int* kimage, int* limage, int jd, int kd, int ld) {
    char titl[21];

    // x,y,z dimensioned (jd,kd,ld); jimage/kimage/limage (jd-1,kd-1,ld-1)
    View3<double> X(x, jd, kd, ld), Y(y, jd, kd, ld), Z(z, jd, kd, ld);
    View3<int> JI(jimage, jd - 1, kd - 1, ld - 1);
    View3<int> KI(kimage, jd - 1, kd - 1, ld - 1);
    View3<int> LI(limage, jd - 1, kd - 1, ld - 1);

    int len = build_scratch_name(titl, "temp_cen.", m);
    std::string fname(titl, len);
    fopen_unit(30, fname.c_str(), "wb+");

    // write(iunit) (((x(j,k,l),j=1,jd-1),k=1,kd-1),l=1,ld-1), then y-block, z-block
    //   j fastest, k, then l  -- ONE record.
    {
        tc::RecordWriter rw;
        for (int l = 1; l <= ld - 1; ++l)
            for (int k = 1; k <= kd - 1; ++k)
                for (int j = 1; j <= jd - 1; ++j) rw.put_d(X(j, k, l));
        for (int l = 1; l <= ld - 1; ++l)
            for (int k = 1; k <= kd - 1; ++k)
                for (int j = 1; j <= jd - 1; ++j) rw.put_d(Y(j, k, l));
        for (int l = 1; l <= ld - 1; ++l)
            for (int k = 1; k <= kd - 1; ++k)
                for (int j = 1; j <= jd - 1; ++j) rw.put_d(Z(j, k, l));
        rw.flush(funit(30));
    }

    // write(iunit) jimage-block, kimage-block, limage-block  -- ONE record
    {
        tc::RecordWriter rw;
        for (int l = 1; l <= ld - 1; ++l)
            for (int k = 1; k <= kd - 1; ++k)
                for (int j = 1; j <= jd - 1; ++j) rw.put_i(JI(j, k, l));
        for (int l = 1; l <= ld - 1; ++l)
            for (int k = 1; k <= kd - 1; ++k)
                for (int j = 1; j <= jd - 1; ++j) rw.put_i(KI(j, k, l));
        for (int l = 1; l <= ld - 1; ++l)
            for (int k = 1; k <= kd - 1; ++k)
                for (int j = 1; j <= jd - 1; ++j) rw.put_i(LI(j, k, l));
        rw.flush(funit(30));
    }

    frewind(30);
    return;
}

// ===========================================================================
//  subroutine getgrd  — get cell-center grid for mesh m from temp_cen.m
// ===========================================================================
void getgrd(int m, double* x, double* y, double* z,
            int* jimage, int* kimage, int* limage, int jd, int kd, int ld) {
    char titl[21];

    int len = build_scratch_name(titl, "temp_cen.", m);
    std::string fname(titl, len);
    fopen_unit(30, fname.c_str(), "rb");

    long n = (long)jd * kd * ld;

    // read(iunit) x,y,z   -- one record; here dims are (jd,kd,ld)
    {
        std::vector<double> tmp((size_t)(3 * n));
        uread(30, tmp.data(), (size_t)(3 * n) * sizeof(double));
        for (long i = 0; i < n; ++i) {
            x[1 + i] = tmp[(size_t)i];
            y[1 + i] = tmp[(size_t)(n + i)];
            z[1 + i] = tmp[(size_t)(2 * n + i)];
        }
    }
    // read(iunit) jimage,kimage,limage  -- one record
    {
        std::vector<int> tmp((size_t)(3 * n));
        uread(30, tmp.data(), (size_t)(3 * n) * sizeof(int));
        for (long i = 0; i < n; ++i) {
            jimage[1 + i] = tmp[(size_t)i];
            kimage[1 + i] = tmp[(size_t)(n + i)];
            limage[1 + i] = tmp[(size_t)(2 * n + i)];
        }
    }
    frewind(30);
    fclose_unit(30);
    return;
}

// ===========================================================================
//  subroutine getgrd2  — get node-point grid for mesh m from temp_grd.m
// ===========================================================================
void getgrd2(int m, double* x, double* y, double* z, int jd, int kd, int ld) {
    char titl[21];

    int len = build_scratch_name(titl, "temp_grd.", m);
    std::string fname(titl, len);
    fopen_unit(30, fname.c_str(), "rb");

    long n = (long)jd * kd * ld;
    // read(iunit) x,y,z   -- one record
    {
        std::vector<double> tmp((size_t)(3 * n));
        uread(30, tmp.data(), (size_t)(3 * n) * sizeof(double));
        for (long i = 0; i < n; ++i) {
            x[1 + i] = tmp[(size_t)i];
            y[1 + i] = tmp[(size_t)(n + i)];
            z[1 + i] = tmp[(size_t)(2 * n + i)];
        }
    }
    frewind(30);
    fclose_unit(30);

    return;
}

// ===========================================================================
//  subroutine branch  — detect branch cuts, build jimage/kimage/limage maps
// ===========================================================================
void branch(int m1, int jd, int kd, int ld,
            int* jimage, int* kimage, int* limage,
            double* x1, double* y1, double* z1) {
    // x1,y1,z1 dimensioned (jd,kd,ld); image arrays (jd-1,kd-1,ld-1)
    View3<double> X1(x1, jd, kd, ld), Y1(y1, jd, kd, ld), Z1(z1, jd, kd, ld);
    View3<int> JIM(jimage, jd - 1, kd - 1, ld - 1);
    View3<int> KIM(kimage, jd - 1, kd - 1, ld - 1);
    View3<int> LIM(limage, jd - 1, kd - 1, ld - 1);

    int ibrnch[6 + 1];

    double eps = 1.e-7;

    // face names are cosmetic; not used in output. (face(1..6) in Fortran)

    // list-directed: write(6,*)'    checking mesh ',m1,' for branch cuts on',' boundaries'
    //   -> ' '(lead) + '    checking mesh ' + I12 + '  for branch cuts on' + ' boundaries'
    std::fprintf(funit(6), "     checking mesh %12d  for branch cuts on boundaries\n", m1);

    // default to no branch cuts
    for (int iface = 1; iface <= 6; ++iface) ibrnch[iface] = 0;
    // do 5 j=1,jd-1 ; k=1,kd-1 ; l=1,ld-1
    for (int j = 1; j <= jd - 1; ++j)
        for (int k = 1; k <= kd - 1; ++k)
            for (int l = 1; l <= ld - 1; ++l) {
                JIM(j, k, l) = j;
                KIM(j, k, l) = k;
                LIM(j, k, l) = l;
            }

    // ---- C-grid type branch cut on k=1 and k=kd boundaries ----
    if (kd > 2) {

        // branch cut in j-direction
        if (jd > 2) {
            // do 100 k=1,kd,kd-1
            for (int k = 1; k <= kd; k += (kd - 1)) {
                int iface = 3;
                if (k == kd) iface = 4;
                if (ibrnch[iface] > 0) continue;   // go to 100
                int kcheck = 0;
                int js0 = 100000, jq0 = 0, ls0 = 100000, lq0 = 0;
                // do 200 l=1,ld-1 ; js=1,jd-2 ; do 300 jq=js+1,jd
                for (int l = 1; l <= ld - 1; ++l)
                    for (int js = 1; js <= jd - 2; ++js)
                        for (int jq = js + 1; jq <= jd; ++jq) {
                            double delxyz = std::abs(X1(js, k, l) - X1(jq, k, l))
                                          + std::abs(Y1(js, k, l) - Y1(jq, k, l))
                                          + std::abs(Z1(js, k, l) - Z1(jq, k, l));
                            if (delxyz < eps) {
                                double del2 = std::abs(X1(js + 1, k, l) - X1(jq - 1, k, l))
                                            + std::abs(Y1(js + 1, k, l) - Y1(jq - 1, k, l))
                                            + std::abs(Z1(js + 1, k, l) - Z1(jq - 1, k, l));
                                if (del2 < eps) {
                                    if (k == 1) {
                                        if (JIM(js, 1, l) == js) {
                                            kcheck = kcheck + 1;
                                            JIM(js, 1, l)     = jq - 2;
                                            JIM(jq - 2, 1, l) = js;
                                        }
                                    }
                                    if (k == kd) {
                                        if (JIM(js, kd - 1, l) == js) {
                                            kcheck = kcheck + 1;
                                            JIM(js, kd - 2, l)     = jq - 2;
                                            JIM(jq - 2, kd - 2, l) = js;
                                        }
                                    }
                                    if (js < js0) js0 = js;
                                    if (jq > jq0) jq0 = jq;
                                    if (l < ls0)  ls0 = l;
                                    if (l > lq0)  lq0 = l;
                                }
                            }
                        } // 300 / 200

                int kk = 0;
                if (k == 1)  kk = 1;
                if (k == kd) kk = kd - 2;
                // do 201 l=1,ld-1 ; do 301 j=2,jd-2
                for (int l = 1; l <= ld - 1; ++l)
                    for (int j = 2; j <= jd - 2; ++j) {
                        if (JIM(j, kk, l) == j) {
                            if (JIM(j + 1, kk, l) != j + 1 && JIM(j - 1, kk, l) != j - 1) {
                                JIM(j, kk, l) = j - 1;
                            }
                        }
                    }

                if (kcheck > 0) {
                    ibrnch[iface] = 10;
                    // write(6,*)'      C-type branch cut detected on k = ',k,' boundary'
                    std::fprintf(funit(6), "       C-type branch cut detected on k = %12d boundary\n", k);
                    // write(6,*)'        for j = ',js0,' to j = ',jq0,' and for l = ',ls0,' to ',' l = ',lq0
                    std::fprintf(funit(6),
                        "         for j = %12d to j = %12d and for l = %12d to  l = %12d\n",
                        js0, jq0, ls0, lq0);
                }
            } // 100
        }

        // branch cut in l-direction
        if (ld > 2) {
            // do 110 k=1,kd,kd-1
            for (int k = 1; k <= kd; k += (kd - 1)) {
                int iface = 3;
                if (k == kd) iface = 4;
                if (ibrnch[iface] > 0) continue;   // go to 110
                int kcheck = 0;
                int js0 = 100000, jq0 = 0, ls0 = 100000, lq0 = 0;
                // do 210 j=1,jd-1 ; ls=1,ld-2 ; do 310 lq=ls+1,ld
                for (int j = 1; j <= jd - 1; ++j)
                    for (int ls = 1; ls <= ld - 2; ++ls)
                        for (int lq = ls + 1; lq <= ld; ++lq) {
                            double delxyz = std::abs(X1(j, k, ls) - X1(j, k, lq))
                                          + std::abs(Y1(j, k, ls) - Y1(j, k, lq))
                                          + std::abs(Z1(j, k, ls) - Z1(j, k, lq));
                            if (delxyz < eps) {
                                double del2 = std::abs(X1(j, k, ls + 1) - X1(j, k, lq - 1))
                                            + std::abs(Y1(j, k, ls + 1) - Y1(j, k, lq - 1))
                                            + std::abs(Z1(j, k, ls + 1) - Z1(j, k, lq - 1));
                                if (del2 < eps) {
                                    if (k == 1) {
                                        if (LIM(j, 1, ls) == ls) {
                                            kcheck = kcheck + 1;
                                            LIM(j, 1, ls)     = lq - 2;
                                            LIM(j, 1, lq - 2) = ls;
                                        }
                                    }
                                    if (k == kd) {
                                        if (LIM(j, kd - 1, ls) == ls) {
                                            kcheck = kcheck + 1;
                                            LIM(j, kd - 2, ls)     = lq - 2;
                                            LIM(j, kd - 2, lq - 2) = ls;
                                        }
                                    }
                                    if (ls < ls0) ls0 = ls;
                                    if (lq > lq0) lq0 = lq;
                                    if (j < js0)  js0 = j;
                                    if (j > jq0)  jq0 = j;
                                }
                            }
                        } // 310 / 210

                int kk = 0;
                if (k == 1)  kk = 1;
                if (k == kd) kk = kd - 2;
                // do 211 j=1,jd-1 ; do 311 l=2,ld-2
                for (int j = 1; j <= jd - 1; ++j)
                    for (int l = 2; l <= ld - 2; ++l) {
                        if (LIM(j, kk, l) == l) {
                            if (LIM(j, kk, l + 1) != l + 1 && LIM(j, kk, l - 1) != l - 1) {
                                LIM(j, kk, l) = l - 1;
                            }
                        }
                    }

                if (kcheck > 0) {
                    ibrnch[iface] = 30;
                    std::fprintf(funit(6), "       C-type branch cut detected on k = %12d boundary\n", k);
                    std::fprintf(funit(6),
                        "         for j = %12d to j = %12d and for l = %12d to  l = %12d\n",
                        js0, jq0, ls0, lq0);
                }
            } // 110
        }

    } // kd>2

    // ---- C-grid type branch cut on j=1 and j=jd boundaries ----
    if (jd > 2) {

        // branch cut in k-direction
        if (kd > 2) {
            // do 120 j=1,jd,jd-1
            for (int j = 1; j <= jd; j += (jd - 1)) {
                int iface = 1;
                if (j == jd) iface = 2;
                int ls0 = 100000, lq0 = 0, ks0 = 100000, kq0 = 0;
                int jcheck = 0;
                // do 220 l=1,ld-1 ; ks=1,kd-2 ; do 320 kq=ks+1,kd
                for (int l = 1; l <= ld - 1; ++l)
                    for (int ks = 1; ks <= kd - 2; ++ks)
                        for (int kq = ks + 1; kq <= kd; ++kq) {
                            double delxyz = std::abs(X1(j, ks, l) - X1(j, kq, l))
                                          + std::abs(Y1(j, ks, l) - Y1(j, kq, l))
                                          + std::abs(Z1(j, ks, l) - Z1(j, kq, l));
                            if (delxyz < eps) {
                                double del2 = std::abs(X1(j, ks + 1, l) - X1(j, kq - 1, l))
                                            + std::abs(Y1(j, ks + 1, l) - Y1(j, kq - 1, l))
                                            + std::abs(Z1(j, ks + 1, l) - Z1(j, kq - 1, l));
                                if (del2 < eps) {
                                    if (j == 1) {
                                        if (KIM(1, ks, l) == ks) {
                                            jcheck = jcheck + 1;
                                            KIM(1, ks, l)     = kq - 2;
                                            KIM(1, kq - 2, l) = ks;
                                        }
                                    }
                                    if (j == jd) {
                                        if (KIM(jd - 1, ks, l) == ks) {
                                            jcheck = jcheck + 1;
                                            KIM(jd - 2, ks, l)     = kq - 2;
                                            KIM(jd - 2, kq - 2, l) = ks;
                                        }
                                    }
                                    if (ks < ks0) ks0 = ks;
                                    if (kq > kq0) kq0 = kq;
                                    if (l < ls0)  ls0 = l;
                                    if (l > lq0)  lq0 = l;
                                }
                            }
                        } // 320 / 220

                int jj = 0;
                if (j == 1)  jj = 1;
                if (j == jd) jj = jd - 2;
                // do 221 l=1,ld-1 ; do 321 k=2,kd-2
                for (int l = 1; l <= ld - 1; ++l)
                    for (int k = 2; k <= kd - 2; ++k) {
                        if (KIM(jj, k, l) == k) {
                            if (KIM(jj, k + 1, l) != k + 1 && KIM(jj, k - 1, l) != k - 1) {
                                KIM(jj, k, l) = k - 1;
                            }
                        }
                    }

                if (jcheck > 0) {
                    ibrnch[iface] = 20;
                    std::fprintf(funit(6), "       C-type branch cut detected on j = %12d boundary\n", j);
                    std::fprintf(funit(6),
                        "         for k = %12d to k = %12d and for l = %12d to  l = %12d\n",
                        ks0, kq0, ls0, lq0);
                }
            } // 120
        }

        // branch cut in l-direction
        if (ld > 2) {
            // do 130 j=1,jd,jd-1
            for (int j = 1; j <= jd; j += (jd - 1)) {
                int iface = 1;
                if (j == jd) iface = 2;
                if (ibrnch[iface] > 0) continue;   // go to 130
                int ls0 = 100000, lq0 = 0, ks0 = 100000, kq0 = 0;
                int jcheck = 0;
                // do 230 k=1,kd-1 ; ls=1,ld-2 ; do 330 lq=ls+1,ld
                for (int k = 1; k <= kd - 1; ++k)
                    for (int ls = 1; ls <= ld - 2; ++ls)
                        for (int lq = ls + 1; lq <= ld; ++lq) {
                            double delxyz = std::abs(X1(j, k, ls) - X1(j, k, lq))
                                          + std::abs(Y1(j, k, ls) - Y1(j, k, lq))
                                          + std::abs(Z1(j, k, ls) - Z1(j, k, lq));
                            if (delxyz < eps) {
                                double del2 = std::abs(X1(j, k, ls + 1) - X1(j, k, lq - 1))
                                            + std::abs(Y1(j, k, ls + 1) - Y1(j, k, lq - 1))
                                            + std::abs(Z1(j, k, ls + 1) - Z1(j, k, lq - 1));
                                if (del2 < eps) {
                                    if (j == 1) {
                                        if (LIM(1, k, ls) == ls) {
                                            jcheck = jcheck + 1;
                                            LIM(1, k, ls)     = lq - 2;
                                            LIM(1, k, lq - 2) = ls;
                                        }
                                    }
                                    if (j == jd) {
                                        if (LIM(jd - 1, k, ls) == ls) {
                                            jcheck = jcheck + 1;
                                            LIM(jd - 2, k, ls)     = lq - 2;
                                            LIM(jd - 2, k, lq - 2) = ls;
                                        }
                                    }
                                    if (ls < ls0) ls0 = ls;
                                    if (lq > lq0) lq0 = lq;
                                    if (k < ks0)  ks0 = k;
                                    if (k > kq0)  kq0 = k;
                                }
                            }
                        } // 330 / 230

                int jj = 0;
                if (j == 1)  jj = 1;
                if (j == jd) jj = jd - 2;
                // do 231 k=1,kd-1 ; do 331 l=2,ld-2
                for (int k = 1; k <= kd - 1; ++k)
                    for (int l = 2; l <= ld - 2; ++l) {
                        if (LIM(jj, k, l) == l) {
                            if (LIM(jj, k, l + 1) != l + 1 && LIM(jj, k, l - 1) != l - 1) {
                                LIM(jj, k, l) = l - 1;
                            }
                        }
                    }

                if (jcheck > 0) {
                    ibrnch[iface] = 30;
                    std::fprintf(funit(6), "       C-type branch cut detected on j = %12d boundary\n", j);
                    std::fprintf(funit(6),
                        "         for k = %12d to k = %12d and for l = %12d to  l = %12d\n",
                        ks0, kq0, ls0, lq0);
                }
            } // 130
        }

    } // jd>2

    // ---- C-grid type branch cut on l=1 and l=ld boundaries ----
    if (ld > 2) {

        // branch cut in j-direction
        if (jd > 2) {
            // do 140 l=1,ld,ld-1
            for (int l = 1; l <= ld; l += (ld - 1)) {
                int iface = 5;
                if (l == ld) iface = 6;
                int js0 = 100000, jq0 = 0, ks0 = 100000, kq0 = 0;
                int lcheck = 0;
                // do 240 k=1,kd-1 ; js=1,jd-2 ; do 340 jq=js+1,jd
                for (int k = 1; k <= kd - 1; ++k)
                    for (int js = 1; js <= jd - 2; ++js)
                        for (int jq = js + 1; jq <= jd; ++jq) {
                            double delxyz = std::abs(X1(js, k, l) - X1(jq, k, l))
                                          + std::abs(Y1(js, k, l) - Y1(jq, k, l))
                                          + std::abs(Z1(js, k, l) - Z1(jq, k, l));
                            if (delxyz < eps) {
                                double del2 = std::abs(X1(js + 1, k, l) - X1(jq - 1, k, l))
                                            + std::abs(Y1(js + 1, k, l) - Y1(jq - 1, k, l))
                                            + std::abs(Z1(js + 1, k, l) - Z1(jq - 1, k, l));
                                if (del2 < eps) {
                                    if (l == 1) {
                                        if (JIM(js, k, 1) == js) {
                                            lcheck = lcheck + 1;
                                            JIM(js, k, 1)     = jq - 2;
                                            JIM(jq - 2, k, 1) = js;
                                        }
                                    }
                                    if (l == ld) {
                                        if (JIM(js, k, ld - 1) == js) {
                                            lcheck = lcheck + 1;
                                            JIM(js, k, ld - 2)     = jq - 2;
                                            JIM(jq - 2, k, ld - 2) = js;
                                        }
                                    }
                                    if (js < js0) js0 = js;
                                    if (jq > jq0) jq0 = jq;
                                    if (k < ks0)  ks0 = k;
                                    if (k > kq0)  kq0 = k;
                                }
                            }
                        } // 340 / 240

                int ll = 0;
                if (l == 1)  ll = 1;
                if (l == ld) ll = ld - 2;
                // do 241 k=1,kd-1 ; do 341 j=2,jd-2
                for (int k = 1; k <= kd - 1; ++k)
                    for (int j = 2; j <= jd - 2; ++j) {
                        if (JIM(j, k, ll) == j) {
                            if (JIM(j + 1, k, ll) != j + 1 && JIM(j - 1, k, ll) != j - 1) {
                                JIM(j, k, ll) = j - 1;
                            }
                        }
                    }

                if (lcheck > 0) {
                    ibrnch[iface] = 10;
                    std::fprintf(funit(6), "       C-type branch cut detected on l = %12d boundary\n", l);
                    std::fprintf(funit(6),
                        "         for j = %12d to j = %12d and for k = %12d to  k = %12d\n",
                        js0, jq0, ks0, kq0);
                }
            } // 140
        }

        // branch cut in k-direction
        if (kd > 2) {
            // do 150 l=1,ld,ld-1
            for (int l = 1; l <= ld; l += (ld - 1)) {
                int iface = 5;
                if (l == ld) iface = 6;
                if (ibrnch[iface] > 0) continue;   // go to 150
                int js0 = 100000, jq0 = 0, ks0 = 100000, kq0 = 0;
                int lcheck = 0;
                // do 250 j=1,jd-1 ; ks=1,kd-2 ; do 350 kq=ks+1,kd
                for (int j = 1; j <= jd - 1; ++j)
                    for (int ks = 1; ks <= kd - 2; ++ks)
                        for (int kq = ks + 1; kq <= kd; ++kq) {
                            double delxyz = std::abs(X1(j, ks, l) - X1(j, kq, l))
                                          + std::abs(Y1(j, ks, l) - Y1(j, kq, l))
                                          + std::abs(Z1(j, ks, l) - Z1(j, kq, l));
                            if (delxyz < eps) {
                                double del2 = std::abs(X1(j, ks + 1, l) - X1(j, kq - 1, l))
                                            + std::abs(Y1(j, ks + 1, l) - Y1(j, kq - 1, l))
                                            + std::abs(Z1(j, ks + 1, l) - Z1(j, kq - 1, l));
                                if (del2 < eps) {
                                    if (l == 1) {
                                        if (KIM(j, ks, 1) == ks) {
                                            lcheck = lcheck + 1;
                                            KIM(j, ks, 1)     = kq - 2;
                                            KIM(j, kq - 2, 1) = ks;
                                        }
                                    }
                                    if (l == ld) {
                                        if (KIM(j, ks, ld - 1) == ks) {
                                            lcheck = lcheck + 1;
                                            KIM(j, ks, ld - 2)     = kq - 2;
                                            KIM(j, kq - 2, ld - 2) = ks;
                                        }
                                    }
                                    if (ks < ks0) ks0 = ks;
                                    if (kq > kq0) kq0 = kq;
                                    if (j < js0)  js0 = j;
                                    if (j > jq0)  jq0 = j;
                                }
                            }
                        } // 350 / 250

                int ll = 0;
                if (l == 1)  ll = 1;
                if (l == ld) ll = ld - 2;
                // do 251 j=1,jd-1 ; do 351 k=2,kd-2
                for (int j = 1; j <= jd - 1; ++j)
                    for (int k = 2; k <= kd - 2; ++k) {
                        if (KIM(j, k, ll) == k) {
                            if (KIM(j, k + 1, ll) != k + 1 && KIM(j, k - 1, ll) != k - 1) {
                                KIM(j, k, ll) = k - 1;
                            }
                        }
                    }

                if (lcheck > 0) {
                    ibrnch[iface] = 20;
                    std::fprintf(funit(6), "       C-type branch cut detected on l = %12d boundary\n", l);
                    std::fprintf(funit(6),
                        "         for j = %12d to j = %12d and for k = %12d to  k = %12d\n",
                        js0, jq0, ks0, kq0);
                }
            } // 150
        }

    } // ld>2

    int iflag = 0;

    // ---- O-type branch cut between k=1 and k=kd boundaries ----
    if (kd > 2) {
        if (ibrnch[3] == 0 && ibrnch[4] == 0) {
            int jcount = 0;
            // do 400 l=1,ld ; j=1,jd
            for (int l = 1; l <= ld; ++l)
                for (int j = 1; j <= jd; ++j) {
                    double delxyz = std::abs(X1(j, 1, l) - X1(j, kd, l))
                                  + std::abs(Y1(j, 1, l) - Y1(j, kd, l))
                                  + std::abs(Z1(j, 1, l) - Z1(j, kd, l));
                    if (delxyz < eps) jcount = jcount + 1;
                }

            if (jcount == jd * ld) {
                ibrnch[3] = 50;
                ibrnch[4] = 50;
                // do 500 j=1,jd-1 ; l=1,ld-1
                for (int j = 1; j <= jd - 1; ++j)
                    for (int l = 1; l <= ld - 1; ++l) {
                        KIM(j, 1, l)      = kd - 2;
                        KIM(j, kd - 2, l) = 1;
                    }
                std::fprintf(funit(6),
                    "       O-type branch cut detected between k = 1 and k = %12d boundaries\n", kd);
                std::fprintf(funit(6),
                    "         for l = 1 to l = %12d and for j = 1 to j = %12d\n", ld, jd);
            } else if (jcount > 0 && jcount < jd * ld) {
                if (iflag > 0) {
                    std::fprintf(funit(6),
                        "       WARNING... %12d O-type branch cut points detected between k = 1 and k = %12d\n",
                        jcount, kd);
                    std::fprintf(funit(6),
                        "       however, this is not the entire j and l range as required\n");
                    std::fprintf(funit(6),
                        "       O-type image points not set...verify as correct\n");
                }
            }
        }
    }

    // ---- O-type branch cut between j=1 and j=jd boundaries ----
    if (jd > 2) {
        if (ibrnch[1] == 0 && ibrnch[2] == 0) {
            int kcount = 0;
            // do 410 l=1,ld ; k=1,kd
            for (int l = 1; l <= ld; ++l)
                for (int k = 1; k <= kd; ++k) {
                    double delxyz = std::abs(X1(1, k, l) - X1(jd, k, l))
                                  + std::abs(Y1(1, k, l) - Y1(jd, k, l))
                                  + std::abs(Z1(1, k, l) - Z1(jd, k, l));
                    if (delxyz < eps) kcount = kcount + 1;
                }

            if (kcount == kd * ld) {
                ibrnch[1] = 40;
                ibrnch[2] = 40;
                // do 510 k=1,kd-1 ; l=1,ld-1
                for (int k = 1; k <= kd - 1; ++k)
                    for (int l = 1; l <= ld - 1; ++l) {
                        JIM(1, k, l)      = jd - 2;
                        JIM(jd - 2, k, l) = 1;
                    }
                std::fprintf(funit(6),
                    "       O-type branch cut detected between j = 1 and j = %12d boundaries\n", jd);
                std::fprintf(funit(6),
                    "         for l = 1 to l = %12d and for k = 1 to k = %12d\n", ld, kd);
            } else if (kcount > 0 && kcount < kd * ld) {
                if (iflag > 0) {
                    std::fprintf(funit(6),
                        "       WARNING... %12d O-type branch cut points detected between j = 1 and j = %12d\n",
                        kcount, jd);
                    std::fprintf(funit(6),
                        "       however, this is not the entire k and l range as required\n");
                    std::fprintf(funit(6),
                        "       O-type image points not set...verify as correct\n");
                }
            }
        }
    }

    // ---- O-type branch cut between l=1 and l=ld boundaries ----
    if (ld > 2) {
        if (ibrnch[5] == 0 && ibrnch[6] == 0) {
            int lcount = 0;
            // do 420 j=1,jd ; k=1,kd
            for (int j = 1; j <= jd; ++j)
                for (int k = 1; k <= kd; ++k) {
                    double delxyz = std::abs(X1(j, k, 1) - X1(j, k, ld))
                                  + std::abs(Y1(j, k, 1) - Y1(j, k, ld))
                                  + std::abs(Z1(j, k, 1) - Z1(j, k, ld));
                    if (delxyz < eps) lcount = lcount + 1;
                }

            if (lcount == jd * kd) {
                ibrnch[5] = 60;
                ibrnch[6] = 60;
                // do 520 j=1,jd-1 ; k=1,kd-1
                for (int j = 1; j <= jd - 1; ++j)
                    for (int k = 1; k <= kd - 1; ++k) {
                        LIM(j, k, 1)      = ld - 2;
                        LIM(j, k, ld - 2) = 1;
                    }
                std::fprintf(funit(6),
                    "       O-type branch cut detected between l = 1 and l = %12d boundaries\n", ld);
                std::fprintf(funit(6),
                    "         for k = 1 to k = %12d and for j = 1 to j = %12d\n", kd, jd);
            } else if (lcount > 0 && lcount < jd * kd) {
                if (iflag > 0) {
                    std::fprintf(funit(6),
                        "       WARNING... %12d O-type branch cut points detected between l = 1 and l = %12d\n",
                        lcount, ld);
                    std::fprintf(funit(6),
                        "       however, this is not the entire j and k range as required\n");
                    std::fprintf(funit(6),
                        "       O-type image points not set...verify as correct\n");
                }
            }
        }
    }

    // message for no branch cuts on block boundaries
    int icount = 0;
    for (int iface = 1; iface <= 6; ++iface) {
        if (ibrnch[iface] != 0) icount = icount + 1;
    }
    if (icount == 0) {
        std::fprintf(funit(6), "       no branch cuts detected in this mesh\n");
        return;
    }
    return;
}

// ===========================================================================
//  program maggie  ->  int main
// ===========================================================================
int main(int argc, char** argv) {
    (void)argc; (void)argv;

    int itime = 0;
    cputim(itime);          // call cputim(0)

    pltpt.nbnd  = 0;
    pltpt.norph = 0;
    pltpt.nill  = 0;

    // epsc: tolerance on interpolation coefficients
    tol.epsc = .25e-3;

    // ---- open files (unit 5 = stdin) ----
    // read(5,*)  -- skip a line
    {
        char skip[8192];
        if (!std::fgets(skip, sizeof(skip), stdin)) { /* EOF: leave */ }
    }

    // read(5,'(a60)') grid  — 60-char filename, then trim trailing blanks to open
    auto read_a60 = [](char dst[81]) {
        char buf[8192];
        if (!std::fgets(buf, sizeof(buf), stdin)) { dst[0] = '\0'; return; }
        // strip newline
        size_t n = std::strlen(buf);
        while (n > 0 && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) buf[--n] = '\0';
        // read exactly a60: take up to 60 chars, blank-pad to 60
        for (int i = 0; i < 60; ++i) dst[i] = (i < (int)n ? buf[i] : ' ');
        dst[60] = '\0';
    };
    auto trimmed = [](const char* s) -> std::string {
        std::string t(s);
        size_t e = t.find_last_not_of(' ');
        if (e == std::string::npos) return std::string();
        return t.substr(0, e + 1);
    };

    read_a60(files.grid);
    fopen_unit(10, trimmed(files.grid).c_str(), "rb");   // status='old', unformatted

    read_a60(files.outpt);
    fopen_unit(6, trimmed(files.outpt).c_str(), "w");    // formatted

    read_a60(files.ovrlp);
    fopen_unit(2, trimmed(files.ovrlp).c_str(), "wb");   // unformatted

    read_a60(files.plt3d);
    fopen_unit(9, trimmed(files.plt3d).c_str(), "wb");   // unformatted

    // ---- output banner (unit 6) ----
    // format 83: 35(2h *)  -> " *" x35 (70 chars)
    const char* line83 = " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *";
    // format 87: 2(2h *),62x,2(2h *)
    const char* line87 = " * *                                                               * *";

    std::fprintf(funit(6), "%s\n", line83);                    // write(6,83)
    std::fprintf(funit(6), "%s\n", line83);                    // write(6,83)
    std::fprintf(funit(6), "%s\n", line87);                    // write(6,87)
    // 9900
    std::fprintf(funit(6),
        " * *           MAGGIE - CFL3D OVERSET-GRID PREPROCESSOR            * *\n");
    std::fprintf(funit(6), "%s\n", line87);                    // write(6,87)
    // 9990 (3 lines)
    std::fprintf(funit(6),
        " * *   VERSION 6.7 :  Computational Fluids Lab, Mail Stop 128,     * *\n");
    std::fprintf(funit(6),
        " * *                  NASA Langley Research Center, Hampton, VA    * *\n");
    std::fprintf(funit(6),
        " * *                  Release Date:  February  1, 2017.            * *\n");
    std::fprintf(funit(6), "%s\n", line87);                    // write(6,87)
    std::fprintf(funit(6), "%s\n", line83);                    // write(6,83)
    std::fprintf(funit(6), "%s\n", line83);                    // write(6,83)

    // format 88: (/19hinput/output files:)  -> blank line then text
    std::fprintf(funit(6), "\ninput/output files:\n");

    // write(6,'(''  '',a60)') grid/outpt/ovrlp/plt3d
    //   -> two spaces + the 60-char (blank-padded) filename field.
    // files.grid/outpt/ovrlp/plt3d already hold the exact a60 blank-padded text.
    std::fprintf(funit(6), "  %s\n", files.grid);
    std::fprintf(funit(6), "  %s\n", files.outpt);
    std::fprintf(funit(6), "  %s\n", files.ovrlp);
    std::fprintf(funit(6), "  %s\n", files.plt3d);

    // ---- set up and initialize ----
    blckdat();
    setup();

    // ---- compute interpolation data ----
    hole();
    outer();

    // ---- diagnostic checks ----
    diagnos();

    // ---- output connection file for flow solver ----
    output();

    itime = -1;
    cputim(itime);          // call cputim(-1)

    // stop -> return from main
    return 0;
}
