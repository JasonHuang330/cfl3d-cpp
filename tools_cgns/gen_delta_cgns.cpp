// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Jason Huang and cfl3d-cpp contributors.
// Part of cfl3d-cpp (https://github.com/JasonHuang330/cfl3d-cpp).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
// gen_delta_cgns.cpp — generate delta.cgns from delta.fmt (plot3d formatted grid)
// + the CFL3D deck's BC section, using the CGNS C API. Mirrors plot3dg_to_cgns.F's
// CFL3D-bctype -> CGNS-BCType mapping so the file round-trips through cgnstools.
//
// Self-consistent-verification path (no Fortran CGNS build): also stores the EXACT
// CFL3D numeric bctype in each BC node's name suffix "_bcNNNN" so the reader can
// recover it losslessly.
//
// build: g++ -O2 -std=c++17 gen_delta_cgns.cpp -I/opt/homebrew/include \
//        -L/opt/homebrew/lib -lcgns -o gen_delta_cgns
// run:   ./gen_delta_cgns <delta.fmt> <delta_cgns.inp> <delta.cgns>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cgnslib.h>

// --- read whitespace/comma-separated tokens, expanding Fortran "n*val" repeats ---
struct Tok {
    std::ifstream f;
    std::vector<std::string> pending;
    size_t pi = 0;
    explicit Tok(const char* fn) : f(fn) {}
    bool next(std::string& out) {
        while (pi >= pending.size()) {
            std::string line;
            if (!std::getline(f, line)) return false;
            pending.clear(); pi = 0;
            // split on space, tab, comma
            std::string cur;
            for (char c : line) {
                if (c==' '||c=='\t'||c==','||c=='\r') { if(!cur.empty()){pending.push_back(cur);cur.clear();} }
                else cur.push_back(c);
            }
            if(!cur.empty()) pending.push_back(cur);
        }
        std::string t = pending[pi++];
        // expand "n*val"
        auto star = t.find('*');
        if (star != std::string::npos) {
            int rep = std::atoi(t.substr(0,star).c_str());
            std::string val = t.substr(star+1);
            for (int r=1;r<rep;r++) pending.insert(pending.begin()+pi, val);
            t = val;
        }
        out = t; return true;
    }
    double nextd() { std::string s; if(!next(s)){fprintf(stderr,"EOF reading grid\n");exit(1);} return std::atof(s.c_str()); }
    int nexti()   { std::string s; if(!next(s)){fprintf(stderr,"EOF reading grid\n");exit(1);} return std::atoi(s.c_str()); }
};

static CGNS_ENUMT(BCType_t) map_bc(int ib) {
    switch (ib) {
        case 1000: return CGNS_ENUMV(BCInflowSupersonic);
        case 1001: return CGNS_ENUMV(BCSymmetryPlane);
        case 1002: return CGNS_ENUMV(BCExtrapolate);
        case 1003: return CGNS_ENUMV(BCFarfield);
        case 1005: case 1006: return CGNS_ENUMV(BCWallInviscid);
        case 1011: return CGNS_ENUMV(BCSymmetryPolar);
        case 1013: return CGNS_ENUMV(BCDegenerateLine);
        default:
            if (ib >= 2000 && ib < 3000) return CGNS_ENUMV(BCWallViscousHeatFlux);
            return CGNS_ENUMV(BCGeneral);
    }
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr,"usage: %s delta.fmt deck.inp out.cgns\n",argv[0]); return 1; }
    const char* gridfn = argv[1];
    const char* deckfn = argv[2];
    const char* outfn  = argv[3];

    // ---- read plot3d formatted grid (single block) ----
    Tok g(gridfn);
    int nblk = g.nexti();
    if (nblk != 1) { fprintf(stderr,"only single-block delta supported (nblk=%d)\n",nblk); return 1; }
    int id = g.nexti(), jd = g.nexti(), kd = g.nexti();
    long n = (long)id*jd*kd;
    printf("grid: %d x %d x %d  (%ld pts)\n", id, jd, kd, n);
    std::vector<double> x(n), y(n), z(n);
    for (long i=0;i<n;i++) x[i]=g.nextd();
    for (long i=0;i<n;i++) y[i]=g.nextd();
    for (long i=0;i<n;i++) z[i]=g.nextd();
    printf("grid read OK; x[0]=%.6g y[0]=%.6g z[0]=%.6g  x[last]=%.6g\n", x[0],y[0],z[0],x[n-1]);

    // ---- write CGNS ----
    remove(outfn);
    int fn, ier;
    ier = cg_open(outfn, CG_MODE_WRITE, &fn); if(ier){cg_error_print();return 1;}
    int B; cg_base_write(fn, "Base", 3, 3, &B);
    cgsize_t sz[9] = { id,jd,kd, id-1,jd-1,kd-1, 0,0,0 };
    int Z; ier = cg_zone_write(fn, B, "Zone1", sz, CGNS_ENUMV(Structured), &Z); if(ier){cg_error_print();return 1;}
    int C;
    cg_coord_write(fn,B,Z,CGNS_ENUMV(RealDouble),"CoordinateX",x.data(),&C);
    cg_coord_write(fn,B,Z,CGNS_ENUMV(RealDouble),"CoordinateY",y.data(),&C);
    cg_coord_write(fn,B,Z,CGNS_ENUMV(RealDouble),"CoordinateZ",z.data(),&C);
    printf("wrote grid coords\n");

    // ---- parse deck BC section and write BC_t nodes ----
    // The 6 face headers appear as "I0:", "IDIM:", "J0:", "JDIM:", "K0:", "KDIM:".
    // Each is followed by nseg lines: grid seg bctype a b c d ndata [+ ndata aux lines].
    std::ifstream dk(deckfn);
    std::string line;
    // face -> (fixed index plane, which two vary)  ; ipnts corner convention from plot3dg_to_cgns
    struct Face { const char* tag; const char* prefix; };
    Face faces[6] = {
        {"I0:","Ilo_Seg"},{"IDIM:","Ihi_Seg"},{"J0:","Jlo_Seg"},
        {"JDIM:","Jhi_Seg"},{"K0:","Klo_Seg"},{"KDIM:","Khi_Seg"}
    };
    int which = -1;          // current face (-1 = none)
    bool inBCzone = false;    // between the GRID/NBCI0 header and MSEQ
    int nbc = 0;
    while (std::getline(dk, line)) {
        if (line.find("MSEQ") != std::string::npos) break;
        // detect a face header line ("I0:  GRID SEGMENT ..." etc.)
        int hdr = -1;
        for (int m=0;m<6;m++) if (line.find(faces[m].tag) != std::string::npos) hdr=m;
        if (hdr >= 0) { which = hdr; inBCzone = true; continue; }
        if (!inBCzone || which < 0) continue;
        // parse a segment data line: gg seg ib a b c d nd
        std::istringstream ss(line);
        int gg,seg,ib,a,b,c,d,nd;
        if (!(ss>>gg>>seg>>ib>>a>>b>>c>>d>>nd)) continue;
        // skip ndata aux lines (a header line + a values line)
        if (nd>0)      { std::getline(dk,line); std::getline(dk,line); }
        else if (nd<0) { std::getline(dk,line); std::getline(dk,line); }
        if (ib==0) continue;   // 1-1 interface placeholder (handled by connectivity)
        {
            // build PointRange corners on this face plane (1-based CGNS vertex indices)
            if (a==0) a=1;
            cgsize_t p[6];
            // face plane fixed index: I0->i=1, IDIM->i=id, J0->j=1, JDIM->j=jd, K0->k=1, KDIM->k=kd
            int planeI= (which==0)?1:(which==1)?id:0;
            int planeJ= (which==2)?1:(which==3)?jd:0;
            int planeK= (which==4)?1:(which==5)?kd:0;
            // a..d are the two in-plane extents (order per deck: for I faces = j then k; J faces = i then k; K faces = i then j)
            int lo1=a, hi1=b, lo2=c, hi2=d;
            if (which<=1) { // I face: vary j(lo1..hi1), k(lo2..hi2)
                if(hi1==0)hi1=jd; if(hi2==0)hi2=kd;
                p[0]=planeI; p[1]=lo1; p[2]=lo2; p[3]=planeI; p[4]=hi1; p[5]=hi2;
            } else if (which<=3) { // J face: vary i, k
                if(hi1==0)hi1=id; if(hi2==0)hi2=kd;
                p[0]=lo1; p[1]=planeJ; p[2]=lo2; p[3]=hi1; p[4]=planeJ; p[5]=hi2;
            } else { // K face: vary i, j
                if(hi1==0)hi1=id; if(hi2==0)hi2=jd;
                p[0]=lo1; p[1]=lo2; p[2]=planeK; p[3]=hi1; p[4]=hi2; p[5]=planeK;
            }
            char name[96];
            snprintf(name,sizeof(name),"%s%04d_bc%d", faces[which].prefix, seg, ib);
            int bcidx;
            cgsize_t pr[2*3] = {p[0],p[1],p[2], p[3],p[4],p[5]};
            ier = cg_boco_write(fn,B,Z,name, map_bc(ib), CGNS_ENUMV(PointRange), 2, pr, &bcidx);
            if(ier){cg_error_print();return 1;}
            cg_boco_gridlocation_write(fn,B,Z,bcidx,CGNS_ENUMV(Vertex));
            nbc++;
        }
    }
    printf("wrote %d BC_t nodes\n", nbc);
    cg_close(fn);
    printf("delta.cgns written OK -> %s\n", outfn);
    return 0;
}
