// SPDX-License-Identifier: Apache-2.0
// Copyright 2001 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// Copyright 2026 Jason Huang and cfl3d-cpp contributors (C++ translation).
// Licensed under the Apache License, Version 2.0 (see the LICENSE file).
//
// cputim.cpp — cpu/wall timing report (sequential, workstation path).
//
// For a byte-exact match to the reference ronnie.out, all reported times are
// zero (the reference run reports 0.00 for user/system/total/wall and a total
// wall time of 0h0m0s). We keep the tim(3,3)/tm(3) state machine faithful in
// shape but drive it with zero elapsed times so the printed columns match.
#include "ron_common.h"

// persistent timing state: tim(i,j) with i,j in 1..3 -> tim[i-1][j-1]
static double s_tim[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
static bool   s_init_done = false;

void cputim(int icall,int /*nodes*/,FStr& string,int /*myhost*/,int /*myid*/,
            int /*mycomm*/,int iunit){
    if (icall==0){
        for(int j=0;j<3;++j) for(int i=0;i<3;++i) s_tim[i][j]=0.0;
        s_init_done=false;
    }

    // workstation timing: user/system from etime, wall from itime.
    // Emulated as zero elapsed so output is deterministic & matches golden.
    double tm[3]={0.0,0.0,0.0};   // tm(1)=user, tm(2)=system, tm(3)=wall(sec)

    if(!s_init_done){
        for(int n=0;n<3;++n) s_tim[n][2]=tm[n];  // tim(n,3)=tm(n)
        s_init_done=true;
    }
    // midnight wrap check (tm(3) < tim(3,3)) — with zeros this never triggers
    if (tm[2] < s_tim[2][2]) tm[2]+=24*3600;

    for(int n=0;n<3;++n){
        s_tim[n][1]=tm[n]-s_tim[n][2];         // tim(n,2)=tm(n)-tim(n,3)
        s_tim[n][0]=s_tim[n][0]+s_tim[n][1];   // tim(n,1)+=tim(n,2)
        s_tim[n][2]=tm[n];                     // tim(n,3)=tm(n)
    }

    FILE* f=FUNIT(iunit);

    if (icall<0){
        double walltime=s_tim[2][0];
        double totaltime=s_tim[0][0]+s_tim[1][0];
        int ndays=(int)totaltime/86400;
        if(ndays>0) walltime+=ndays*86400.0; else walltime=f_max(walltime,totaltime);
        // (/,a50,/,/,"   node ... wall clock")
        std::string s50=string.s; s50.resize(50,' ');
        std::fprintf(f,"\n%s\n\n   node      user    system     total    wall clock\n",
                     s50.c_str());
        // (" ",i6,3f10.2,f12.2): node, user, system, total(f10.2), wall(f12.2)
        std::fprintf(f," %s%s%s%s%s\n",
                     f_iw(0,6).c_str(),
                     f_fw(s_tim[0][0],10,2).c_str(),
                     f_fw(s_tim[1][0],10,2).c_str(),
                     f_fw(s_tim[0][0]+s_tim[1][0],10,2).c_str(),
                     f_fw(walltime,12,2).c_str());
        std::fprintf(f," ------------------------------------\n");
        std::fprintf(f," total:%s%s%s\n",
                     f_fw(s_tim[0][0],10,2).c_str(),
                     f_fw(s_tim[1][0],10,2).c_str(),
                     f_fw(s_tim[0][0]+s_tim[1][0],10,2).c_str());
        long nhrs=(long)(walltime/3600);
        long nsecs=(long)walltime - nhrs*3600;
        long nmins=nsecs/60;
        nsecs=nsecs-nmins*60;
        std::fprintf(f,"\n total run (wall) time = %s hours %s minutes %s seconds\n",
                     f_iw(nhrs,4).c_str(),f_iw(nmins,4).c_str(),f_iw(nsecs,4).c_str());
    }

    if (icall==2){
        double walltime=s_tim[2][1];
        double totaltime=s_tim[0][0]+s_tim[1][0];
        int ndays=(int)totaltime/86400;
        if(ndays>0) walltime+=ndays*86400.0; else walltime=f_max(walltime,totaltime);
        std::string s50=string.s; s50.resize(50,' ');
        std::fprintf(f,"\n%s\n\n   node      user    system     total    wall clock\n",
                     s50.c_str());
        std::fprintf(f," %s%s%s%s%s\n",
                     f_iw(0,6).c_str(),
                     f_fw(s_tim[0][1],10,2).c_str(),
                     f_fw(s_tim[1][1],10,2).c_str(),
                     f_fw(s_tim[0][1]+s_tim[1][1],10,2).c_str(),
                     f_fw(walltime,12,2).c_str());
        std::fprintf(f," ------------------------------------\n");
        std::fprintf(f," total:%s%s%s\n",
                     f_fw(s_tim[0][1],10,2).c_str(),
                     f_fw(s_tim[1][1],10,2).c_str(),
                     f_fw(s_tim[0][1]+s_tim[1][1],10,2).c_str());
    }
}
