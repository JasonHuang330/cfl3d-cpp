// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// tools_common.h — shared helpers for translating CFL3D tools/ to C++.
// Covers Fortran unformatted record I/O, mixed-type record writing, text
// parsing (line reads + backspace + list-directed), and Fortran Ew.d output.
//
// Build tools with: g++ -std=c++17 -O2 -ffp-contract=off  (to match gfortran).

#pragma once
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>

namespace tc {

// ---- Fortran unformatted sequential record I/O -----------------------------
// gfortran wraps each record with a 4-byte little-endian length marker before
// AND after the payload. read copies up to nbytes and skips the remainder
// (Fortran READ consumes the whole record even when the list needs fewer values).
inline void urec_read(FILE* f, void* buf, size_t nbytes) {
    int32_t lead=0, trail=0;
    if (std::fread(&lead,4,1,f)!=1){ std::fprintf(stderr,"urec_read: marker\n"); std::exit(1);}
    size_t reclen=(size_t)lead, take=nbytes<reclen?nbytes:reclen;
    if (take && std::fread(buf,1,take,f)!=take){ std::fprintf(stderr,"urec_read: data\n"); std::exit(1);}
    if (reclen>take) std::fseek(f,(long)(reclen-take),SEEK_CUR);
    if (std::fread(&trail,4,1,f)!=1){ std::fprintf(stderr,"urec_read: marker\n"); std::exit(1);}
    (void)trail;
}
// returns false at EOF instead of exiting (for read-until-EOF loops)
inline bool urec_read_try(FILE* f, void* buf, size_t nbytes) {
    int32_t lead=0, trail=0;
    if (std::fread(&lead,4,1,f)!=1) return false;
    size_t reclen=(size_t)lead, take=nbytes<reclen?nbytes:reclen;
    if (take && std::fread(buf,1,take,f)!=take) return false;
    if (reclen>take) std::fseek(f,(long)(reclen-take),SEEK_CUR);
    if (std::fread(&trail,4,1,f)!=1) return false;
    (void)trail; return true;
}
inline void urec_write(FILE* f, const void* buf, size_t nbytes) {
    int32_t len=(int32_t)nbytes;
    std::fwrite(&len,4,1,f); std::fwrite(buf,1,nbytes,f); std::fwrite(&len,4,1,f);
}

// buffer mixed int/double fields into a single Fortran record
struct RecordWriter {
    std::vector<uint8_t> buf;
    void put_d(double v)  { const uint8_t* p=reinterpret_cast<uint8_t*>(&v); buf.insert(buf.end(),p,p+8); }
    void put_f(float v)   { const uint8_t* p=reinterpret_cast<uint8_t*>(&v); buf.insert(buf.end(),p,p+4); }
    void put_i(int32_t v) { const uint8_t* p=reinterpret_cast<uint8_t*>(&v); buf.insert(buf.end(),p,p+4); }
    void put_bytes(const void* p, size_t n){ const uint8_t* b=(const uint8_t*)p; buf.insert(buf.end(),b,b+n); }
    void flush(FILE* f) {
        int32_t len=(int32_t)buf.size();
        std::fwrite(&len,4,1,f); std::fwrite(buf.data(),1,buf.size(),f); std::fwrite(&len,4,1,f);
        buf.clear();
    }
};

// ---- text parsing (for formatted input files) ------------------------------
struct TextReader {
    FILE* f = nullptr;
    long prevpos = 0;
    explicit TextReader(FILE* fp): f(fp) {}
    std::string line() {                       // read(u,'(aN)') name
        prevpos = std::ftell(f);
        char buf[8192];
        if (!std::fgets(buf, sizeof(buf), f)) return std::string();
        std::string s(buf);
        while (!s.empty() && (s.back()=='\n'||s.back()=='\r')) s.pop_back();
        return s;
    }
    void backspace() { std::fseek(f, prevpos, SEEK_SET); }   // backspace(u)
    void ints(int* out, int n) {                             // read(u,*) list of ints
        std::string l = line(); std::stringstream ss(l);
        for (int i=0;i<n;++i){ std::string t; if(!(ss>>t)){out[i]=0;continue;}
            for(char&c:t) if(c==',')c=' '; out[i]=std::atoi(t.c_str()); }
    }
    void dbls(double* out, int n) {                          // read(u,*) list of reals
        std::string l = line(); std::stringstream ss(l);
        for (int i=0;i<n;++i){ std::string t; if(!(ss>>t)){out[i]=0;continue;}
            for(char&c:t){ if(c==',')c=' '; else if(c=='d'||c=='D')c='e'; } out[i]=std::atof(t.c_str()); }
    }
};

// ---- Fortran Ew.d formatted float output -----------------------------------
// Fortran normalizes the mantissa to [0.1,1.0) → "0.2000E+00"; C's %E uses
// [1,10) → "2.0000E-01". This emulates Fortran, right-justified in width w.
inline std::string fmtE(double v, int w, int d) {
    bool neg=false; if (v<0){ neg=true; v=-v; }
    int exp10=0; long im=0;
    if (v==0.0){ exp10=0; im=0; }
    else {
        exp10=(int)std::floor(std::log10(v))+1;
        double mant=v/std::pow(10.0,exp10), scale=std::pow(10.0,d);
        im=std::llround(mant*scale); long cap=std::llround(scale);
        if (im>=cap){ im=cap/10; ++exp10; }
        if (im!=0 && im<cap/10){ im*=10; --exp10; }
    }
    char digits[40]; std::snprintf(digits,sizeof(digits),"%0*ld",d,im);
    char es=exp10<0?'-':'+'; int ea=std::abs(exp10);
    std::string core="0."; core+=digits; core+="E"; core+=es;
    char eb[8]; std::snprintf(eb,sizeof(eb),"%02d",ea); core+=eb;
    if (neg) core="-"+core;
    // gfortran drops the optional leading '0' (0.ddd -> .ddd) when the full
    // form would overflow the field width.
    if ((int)core.size()>w){ size_t p=(size_t)(neg?1:0); if(p+1<core.size()&&core[p]=='0'&&core[p+1]=='.') core.erase(p,1); }
    if ((int)core.size()<w) core=std::string(w-core.size(),' ')+core;
    return core;
}

// ---- stdin prompt helpers --------------------------------------------------
inline void stdin_line(char* dst, size_t cap=81) {
    if (!std::fgets(dst,(int)cap,stdin)){ dst[0]='\0'; return; }
    size_t n=std::strlen(dst);
    while (n>0 && (dst[n-1]=='\n'||dst[n-1]=='\r'||dst[n-1]==' '||dst[n-1]=='\t')) dst[--n]='\0';
}
inline int stdin_int(){ int v; if(std::scanf("%d",&v)!=1) v=0; return v; }
inline double stdin_dbl(){ double v; if(std::scanf("%lf",&v)!=1) v=0; return v; }
inline void stdin_flush_line(){ int c; while((c=std::getchar())!='\n'&&c!=EOF){} }

} // namespace tc
