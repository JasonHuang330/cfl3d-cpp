// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// mag_runtime.cpp — COMMON-block global instances, heap-backed big structs,
// Fortran unit-table I/O, wkstn helpers, and deterministic cputim.
#include "mag_common.h"

// ---- COMMON instances ----
Book1 book1;   Book2 book2;   Book3 book3;
Bound1 bound1; Bound2 bound2; Grdlim grdlim;
Grid1 grid1;   Grid2 grid2;   Conec1 conec1;
Image image;   Igrdtyp igrdtyp; Iseq iseqc;
Diagno diagno; Snafu snafu;   Pltpt pltpt;
Files files;   Tol tol;       Intrp1 intrp1;
Intrp2 intrp2; Where wherec;  Surf surf;
Surf2 surf2;   Norm normc;    Trace1 trace1;
Chkst chkst;   Chkpt chkpt;   Dimen dimenc;
Work workc;    Iwork iworkc;

Temp::Temp()  { p = new double[(size_t)ndim*ndim*6 + 1](); }
Temp::~Temp() { delete[] p; }
Temp tempc;

Temp1::Temp1() {
    x1p = new double[(size_t)jdim*kdim*ldim + 1]();
    y1p = new double[(size_t)jdim*kdim*ldim + 1]();
    z1p = new double[(size_t)jdim*kdim*ldim + 1]();
}
Temp1::~Temp1() { delete[] x1p; delete[] y1p; delete[] z1p; }
Temp1 temp1c;

// ---- Fortran unit table ----
static const int MAXUNIT = 100;
static FILE* g_units[MAXUNIT] = {0};

FILE* funit(int u) { return (u>=0 && u<MAXUNIT) ? g_units[u] : nullptr; }

void fopen_unit(int u, const char* name, const char* mode) {
    if (g_units[u]) std::fclose(g_units[u]);
    g_units[u] = std::fopen(name, mode);
    if (!g_units[u]) { std::fprintf(stderr,"cannot open %s (unit %d)\n",name,u); std::exit(1); }
}
void frewind(int u)     { if (g_units[u]) std::rewind(g_units[u]); }
void fclose_unit(int u) { if (g_units[u]) { std::fclose(g_units[u]); g_units[u]=nullptr; } }

void uwrite(int u, const void* p, size_t n) { tc::urec_write(g_units[u], p, n); }
void uread (int u, void* p, size_t n)       { tc::urec_read (g_units[u], p, n); }

// ---- wkstn.F: ismin / ismax (1-based, returns 1-based index) ----
// Fortran: do 100 nn=2,n,incx.  x is 1-based here (x[1..]).
int ismin(int n,double* x,int incx) {
    int idx = 1; double xmin = x[1];
    for (int nn=2; nn<=n; nn+=incx) {
        if (x[nn] < xmin) { idx = nn; xmin = x[nn]; }
    }
    return idx;
}
int ismax(int n,double* x,int incx) {
    int idx = 1; double xmax = x[1];
    for (int nn=2; nn<=n; nn+=incx) {
        if (x[nn] > xmax) { idx = nn; xmax = x[nn]; }
    }
    return idx;
}

// ---- deterministic cputim (mirrors patched Fortran cputim.F) ----
// Reports all timings as zero so unit-6 output is reproducible.
static float tim_[3+1][3+1]; // 1-based [1..3][1..3]
void cputim(int icall) {
    if (icall == 0) {
        for (int j=1;j<=3;++j) for (int i=1;i<=3;++i) tim_[i][j]=0.f;
    }
    for (int n=1;n<=3;++n) { tim_[n][1]=0.f; tim_[n][2]=0.f; tim_[n][3]=0.f; }
    if (icall < 0) {
        FILE* f = funit(6);
        std::fprintf(f, "\n\n"
          "                  time in seconds\n\n"
          "   node      user    system     total    wall clock\n");
        // format(" ",i6,3f10.2,f12.2)
        std::fprintf(f, " %6d%10.2f%10.2f%10.2f%12.2f\n",
                     0, (double)tim_[1][1], (double)tim_[2][1],
                     (double)(tim_[1][1]+tim_[2][1]), (double)tim_[3][1]);
        std::fprintf(f, " ------------------------------------\n");
        std::fprintf(f, " total:%10.2f%10.2f%10.2f\n",
                     (double)tim_[1][1], (double)tim_[2][1],
                     (double)(tim_[1][1]+tim_[2][1]));
        int nhrs  = (int)(tim_[3][1] / 3600.f);
        int nsecs = (int)(tim_[3][1] - nhrs*3600.f);
        int nmins = nsecs / 60;
        nsecs = nsecs - nmins*60;
        std::fprintf(f, "\n total run (wall) time = %4d hours %4d minutes %4d seconds\n",
                     nhrs, nmins, nsecs);
    }
}
