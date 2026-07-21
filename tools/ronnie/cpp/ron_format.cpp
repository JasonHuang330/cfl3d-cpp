// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// ron_format.cpp — byte-exact emulation of Fortran FORMAT output for ronnie.
#include "ron_common.h"

// Fortran Iw: right-justified integer in a field of width w (blanks if too wide
// gfortran prints '*' * w on overflow, but our values always fit).
std::string f_iw(long v,int w){
    char tmp[64];
    std::snprintf(tmp,sizeof(tmp),"%ld",v);
    std::string s(tmp);
    if ((int)s.size()>w) return std::string(w,'*');
    if ((int)s.size()<w) s=std::string(w-s.size(),' ')+s;
    return s;
}

// Fortran list-directed integer output (write(u,*) i): gfortran emits a
// leading blank followed by the minimal-width decimal representation.
std::string f_ild(long v){
    char tmp[64];
    std::snprintf(tmp,sizeof(tmp)," %ld",v);
    return std::string(tmp);
}

// Fortran list-directed real output (write(u,*) x) approximation. gfortran uses
// a leading blank + a general (G) field. Used only in isklt1-gated diagnostic
// buffer lines (unit 25, secondary output); not part of patch.bin.
std::string f_glist(double v){
    char tmp[64];
    std::snprintf(tmp,sizeof(tmp)," %.7g",v);
    return std::string(tmp);
}

// Fortran Fw.d: fixed-point, right-justified in width w.
std::string f_fw(double v,int w,int d){
    char fmt[16]; std::snprintf(fmt,sizeof(fmt),"%%.%df",d);
    char tmp[64]; std::snprintf(tmp,sizeof(tmp),fmt,v);
    std::string s(tmp);
    if ((int)s.size()<w) s=std::string(w-s.size(),' ')+s;
    return s;
}

// Banner shared by preronnie (unit 66) and ronnie (unit 25).
// isron=0 -> "PRERONNIE - RONNIE PREPROCESSOR"; 1 -> "RONNIE - CFL3D PATCHED-GRID
// PREPROCESSOR". Layout replicated from sizer.F/ronnie.F FORMAT 83/87/9900/9990.
//
// FORMAT 83 = 35(2h *) -> 35 copies of " *" = " * * ... *" (70 chars).
// A leading blank (Fortran carriage-control-free list starts col 1) — gfortran
// with default does NOT add a control char for these list writes; the "2h *"
// literal itself begins with a space so line = " * * ..." (35 times).
void ron_banner(int iunit,int isron,double memuse){
    FILE* f=FUNIT(iunit);
    // line of 35 " *"
    std::string star35; for(int i=0;i<35;++i) star35+=" *";
    // FORMAT 87: 2(2h *),62x,2(2h *)  ->  " * *" + 62 spaces + " * *"
    std::string blankrow=" * *"+std::string(62,' ')+" * *";
    std::fprintf(f,"%s\n",star35.c_str());
    std::fprintf(f,"%s\n",star35.c_str());
    std::fprintf(f,"%s\n",blankrow.c_str());
    // 9900 line: 2(2h *) + text + trailing + 4x + 2(2h *)
    if(isron){
        // ronnie: 38h"           RONNIE - CFL3D PATCHED-GRID" + 13h" PREPROCESSOR"
        //         + 11x + 2(2h *)
        std::string t=" * *";
        t+="           RONNIE - CFL3D PATCHED-GRID"; // 38 chars
        t+=" PREPROCESSOR";                            // 13 chars
        t+=std::string(11,' ');
        t+=" * *";
        std::fprintf(f,"%s\n",t.c_str());
    } else {
        // preronnie: 47h"                PRERONNIE - RONNIE PREPROCESSOR"
        //            + 11h"           " + 4x + 2(2h *)
        std::string t=" * *";
        t+="                PRERONNIE - RONNIE PREPROCESSOR"; // 47 chars
        t+="           ";                                      // 11 chars
        t+=std::string(4,' ');
        t+=" * *";
        std::fprintf(f,"%s\n",t.c_str());
    }
    std::fprintf(f,"%s\n",blankrow.c_str());
    // 9990 (three lines)
    {
        // line1: 2(2h *) + 43h"   VERSION 6.7 :  Computational Fluids Lab,"
        //        + 15h" Mail Stop 128," + 4x + 2(2h *)
        std::string t=" * *";
        t+="   VERSION 6.7 :  Computational Fluids Lab,"; // 43
        t+=" Mail Stop 128,";                              // 15
        t+=std::string(4,' ');
        t+=" * *";
        std::fprintf(f,"%s\n",t.c_str());
        // line2: 2(2h *) + 18x + 41h"NASA Langley Research Center, Hampton, VA"
        //        + 3x + 2(2h *)
        std::string t2=" * *";
        t2+=std::string(18,' ');
        t2+="NASA Langley Research Center, Hampton, VA"; // 41
        t2+=std::string(3,' ');
        t2+=" * *";
        std::fprintf(f,"%s\n",t2.c_str());
        // line3: 2(2h *) + 18x + 33h"Release Date:  February  1, 2017."
        //        + 11x + 2(2h *)
        std::string t3=" * *";
        t3+=std::string(18,' ');
        t3+="Release Date:  February  1, 2017."; // 33
        t3+=std::string(11,' ');
        t3+=" * *";
        std::fprintf(f,"%s\n",t3.c_str());
    }
    std::fprintf(f,"%s\n",blankrow.c_str());
    std::fprintf(f,"%s\n",star35.c_str());
    std::fprintf(f,"%s\n",star35.c_str());
    // memory allocation line (FORMAT 13, single precision):
    //   (/,' memory allocation: ',f12.6,' Mbytes, single precision')
    std::fprintf(f,"\n memory allocation: %s Mbytes, single precision\n",
                 f_fw(memuse,12,6).c_str());
}

// input/output files block (FORMAT 88 + the three '(''  '',a60)' writes)
void ron_iofiles(int iunit,const FStr& grid,const FStr& output,const FStr& patch){
    FILE* f=FUNIT(iunit);
    std::fprintf(f,"\ninput/output files:\n");
    // '(''  '',a60)' -> two spaces then 60-char field
    auto a60=[&](const FStr& s){ std::string v=s.s; v.resize(60,' ');
        std::fprintf(f,"  %s\n",v.c_str()); };
    a60(grid); a60(output); a60(patch);
}
