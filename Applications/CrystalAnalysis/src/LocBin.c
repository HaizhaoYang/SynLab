#include "mex.h" /* Always include this */
#include <math.h>
#include "matrix.h"
int findmax(double *A, int L) {
    int pos = 0;
    int i;
    for(i=0;i<L;i++) {
        if (A[pos]<A[i])
            pos = i;
    }
    return pos;
}

int findmax2(double *A, int st, int L) {
    int pos = st;
    int i;
    for(i=st;i<L+st;i++) {
        if (A[pos]<A[i])
            pos = i;
    }
    return pos;
}


/*This function computes the agl in num_wave directions*/
void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
int nrhs, const mxArray *prhs[]) /* Input variables */
{
    #define agl(ai,bi,j) agl[ai+dims[0]*(bi+dims[1]*j)]
    #define R(ai,bi,j) R[ai+dims[0]*(bi+dims[1]*j)]
    #define TTEng_1st(ai,bi,j) TTEng_1st[ai+dims[0]*(bi+dims[1]*j)]
    #define TTEng_2nd(ai,bi,j) TTEng_2nd[ai+dims[0]*(bi+dims[1]*j)]
    #define W_sec(ai,bi,j) W_sec[ai+dims[0]*(bi+dims[1]*j)]
    #define ss_energy(ci,di,ai,bi) ss_energy[ci+Nss[0]*(di+Nss[1]*(ai+Nss[2]*bi))]
    
    size_t ai, bi, ci, di, k, j, cnt;
    int num_wave, pos_temp;
    double *ss_energy, *agl, *R, *TTEng_1st, *TTEng_2nd, *W_sec;
    ss_energy = mxGetPr(prhs[0]);
    const mwSize *Nss = mxGetDimensions(prhs[0]);
    /*SPg = mxGetPr(prhs[1]);*/
    num_wave = mxGetScalar(prhs[1]);
    nrhs = 2;
    
    nlhs = 5;
    int ndim = 3, dims[3] = {Nss[2],Nss[3],num_wave}, numm = (int)Nss[1]/num_wave/8, numm_R = (int)floor(Nss[0]/4);
    plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[3] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[4] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    agl = mxGetPr(plhs[0]);
    R = mxGetPr(plhs[1]);
    TTEng_1st = mxGetPr(plhs[2]);
    TTEng_2nd = mxGetPr(plhs[3]);
    W_sec = mxGetPr(plhs[4]);
    
    int L = round(Nss[1]/num_wave), maxpos, st, ed, stR, edR, moveStep, pos, sgn, temp_di;
    double *temp, energy_sum, *temp_R, *temp2, *tempss, tpR;
    temp = (double *)mxMalloc(sizeof(double)*L);
    temp2 = (double *)mxMalloc(sizeof(double)*Nss[1]);
    tempss = (double *)mxMalloc(sizeof(double)*Nss[1]);
    temp_R = (double *)mxMalloc(sizeof(double)*Nss[0]);
    
    for (ai=0;ai<Nss[2];ai++) {
        for (bi=0;bi<Nss[3];bi++) {
            /*compute tempss*/
            for (di=0;di<Nss[1];di++) {
                for (ci=0;ci<Nss[0];ci++) {
                    if (ci==0)
                        temp2[di] = ss_energy(ci,di,ai,bi);
                    else
                        temp2[di] = temp2[di] + ss_energy(ci,di,ai,bi);
                }
            }
            moveStep = floor(L/4);
            /*maxpos = fmod(findmax(temp2,Nss[1]),L);*/
            maxpos = findmax(temp2,Nss[1]);
            maxpos = maxpos - floor(maxpos/L)*L;
            moveStep = floor(maxpos-L/2);
            /*if (maxpos>3*L/4) {*/
            if (moveStep>=0) {
                sgn = 1;
                for (cnt=0;cnt<Nss[1]-moveStep;cnt++) {
                    tempss[cnt] = temp2[cnt+moveStep];
                }
                for (cnt=Nss[1]-moveStep;cnt<Nss[1];cnt++) {
                    tempss[cnt] = temp2[cnt -(Nss[1]-moveStep)];
                }            }
            /*if (maxpos<L/4) {*/
            if (moveStep<0) {
                moveStep = -moveStep;
                sgn = -1;
                for (cnt = 0; cnt < Nss[1]-moveStep; cnt++) {
                    tempss[cnt+moveStep] = temp2[cnt];
                }
                for (cnt = Nss[1]-moveStep; cnt < Nss[1]; cnt++) {
                    tempss[cnt-(Nss[1]-moveStep)] = temp2[cnt];
                }            }/*
            if (maxpos>=L/4 & maxpos <=3*L/4) {
                sgn = 0;
                for (cnt=0;cnt<Nss[1];cnt++) {
                    tempss[cnt] = temp2[cnt];
                }
            }*/
            for (j=0;j<num_wave;j++) {
                for (cnt=0;cnt<L;cnt++) {
                    W_sec(ai,bi,j) = W_sec(ai,bi,j) + tempss[cnt+j*L];
                }
            }
            /*tempss[:] computed*/
            for (j=0;j<num_wave;j++) {
                /*compute temp[:]*/
            /*    for (k=0;k<L;k++) {
                    di = j*L + k;
                    temp[k] = tempss[di];
                }*/
                /*temp[:] computed*/
                maxpos = findmax2(tempss,j*L,L);
                ed = maxpos;
                cnt = 1;
                while (ed<Nss[1]-1 & cnt==1 & tempss[ed+1]>0 & ed-maxpos < numm) {
                    if (tempss[ed]>=tempss[ed+1])
                        ed++;
                    else
                        cnt = 0;
                }
                st = maxpos;
                cnt = 1;
                while (st>0 & cnt == 1 & tempss[st-1]>0 & maxpos-st < numm) {
                    if (tempss[st] >= tempss[st-1])
                        st--;
                    else
                        cnt = 0;
                }
                /*compute agl(ai,bi,j)*/
                energy_sum = 0;
                for (cnt = st; cnt <= ed; cnt++)
                    energy_sum = energy_sum + tempss[cnt];
                /*compute TTEng_1st(ai,bi,j)*/
                TTEng_1st(ai,bi,j) = energy_sum;
                for (cnt = st; cnt <= ed; cnt++) {
                    agl(ai,bi,j) = agl(ai,bi,j) + tempss[cnt]*(cnt+0.5)/energy_sum;
                    tempss[cnt] = 0;
                }
                /*agl(ai,bi,j) = fmod(agl(ai,bi,j) + sgn*moveStep,Nss[1]);*/
                agl(ai,bi,j) = agl(ai,bi,j) + sgn*moveStep;
             /*   if (agl(ai,bi,j)>=Nss[1]) {
                    agl(ai,bi,j) = agl(ai,bi,j) - Nss[1];
                }
                else {
                    if (agl(ai,bi,j)<0){
                        agl(ai,bi,j) = agl(ai,bi,j) + Nss[1];
                    }
                }*/
                /*compute R(ai,bi,j)*/
                st = st + sgn*moveStep;
                if (st>=Nss[1]) {
                    st = st - Nss[1];
                }
                else {
                    if (st<0){
                        st = st + Nss[1];
                    }
                }
                ed = ed + sgn*moveStep;
                if (ed>=Nss[1]) {
                    ed = ed - Nss[1];
                }
                else {
                    if (ed<0){
                        ed = ed + Nss[1];
                    }
                }
                /*st = (int)fmod(st + sgn*moveStep,Nss[1]);
                ed = (int)fmod(ed + sgn*moveStep,Nss[1]);
                mexPrintf("%d %d\n",st,ed);*/
                
                
                /*compute R*/
                /*compute temp_R*/
                if (1) {
                    if (st<=ed) {
                        for (di=st;di<=ed;di++) {
                            for (ci = 0;ci<Nss[0];ci++) {
                                if (di==st) {
                                    temp_R[ci] = ss_energy(ci,di,ai,bi);
                                }
                                else
                                    temp_R[ci] = temp_R[ci] + ss_energy(ci,di,ai,bi);
                            }
                        }
                        /* mexPrintf("%d %d %d\n",ai,bi,j);*/
                        maxpos = findmax(temp_R,Nss[0]);
                        /* mexPrintf("%d\n",maxpos);*/
                        ed = maxpos;
                        cnt = 1;
                        while (ed<Nss[0]-1 & cnt==1 & temp_R[ed+1]>0 & ed-maxpos < numm_R) {
                            if (temp_R[ed]>=temp_R[ed+1]) {
                                ed++;
                            }
                            else
                                cnt = 0;
                        }
                        st = maxpos;
                        cnt = 1;
                        while (st>0 & cnt == 1 & temp_R[st-1]>0 & maxpos-st < numm_R) {
                            if (temp_R[st] >= temp_R[st-1]) {
                                st--;
                            }
                            else
                                cnt = 0;
                        }
                        
                        /*mexPrintf("%d %d\n",st,ed);*/
                        energy_sum = 0;
                        for (cnt = st; cnt <= ed; cnt++) {
                            energy_sum = energy_sum + temp_R[cnt];
                        }
                        if (energy_sum>0) {
                            for (cnt = st; cnt <= ed; cnt++) {
                                R(ai,bi,j) = R(ai,bi,j) + temp_R[cnt]*(cnt)/energy_sum;
                            }
                        }
                        
                    }
                    else {
                        ed = ed + Nss[1];
                        for (di=st;di<=ed;di++) {
                            for (ci = 0;ci<Nss[0];ci++) {
                                if (di==st) {
                                    temp_R[ci] = ss_energy(ci,di,ai,bi);
                                }
                                else {
                                    if (di>=Nss[1]) {
                                        temp_R[ci] = temp_R[ci] + ss_energy(ci,di-Nss[1],ai,bi);
                                    }
                                    else
                                        temp_R[ci] = temp_R[ci] + ss_energy(ci,di,ai,bi);
                                }
                            }
                        }
                        maxpos = findmax(temp_R,Nss[0]);
                        ed = maxpos;
                        cnt = 1;
                        while (ed<Nss[0]-1 & cnt==1 & temp_R[ed+1]>0 & ed-maxpos < numm_R) {
                            if (temp_R[ed]>=temp_R[ed+1]) {
                                ed++;
                            }
                            else
                                cnt = 0;
                        }
                        st = maxpos;
                        cnt = 1;
                        while (st>0 & cnt == 1 & temp_R[st-1]>0 & maxpos-st < numm_R) {
                            if (temp_R[st] >= temp_R[st-1]) {
                                st--;
                            }
                            else
                                cnt = 0;
                        }
                        energy_sum = 0;
                        for (cnt = st; cnt <= ed; cnt++) {
                            energy_sum = energy_sum + temp_R[cnt];
                        }
                        if (energy_sum>0) {
                            for (cnt = st; cnt <= ed; cnt++) {
                                R(ai,bi,j) = R(ai,bi,j) + temp_R[cnt]*(cnt)/energy_sum;
                            }
                        }
                    }
                }
                
                
                
                /* mexPrintf("%d %d %d\n",ai,bi,j);*/
                
                /*compute TTEng_2nd(ai,bi,j)*/
                maxpos = findmax2(tempss,j*L,L);
                if (tempss[maxpos]>0) {
                    ed = maxpos;
                    cnt = 1;
                    while (ed<Nss[1]-1 & cnt==1 & tempss[ed+1]>0 & ed-maxpos < numm) {
                        if (tempss[ed]>=tempss[ed+1])
                            ed++;
                        else
                            cnt = 0;
                    }
                    st = maxpos;
                    cnt = 1;
                    while (st>0 & cnt == 1 & tempss[st-1]>0 & maxpos-st < numm) {
                        if (tempss[st] >= tempss[st-1])
                            st--;
                        else
                            cnt = 0;
                    }
                    energy_sum = 0;
                    for (cnt = st; cnt <= ed; cnt++)
                        energy_sum = energy_sum + tempss[cnt];
                    if (energy_sum<=TTEng_1st(ai,bi,j))
                        TTEng_2nd(ai,bi,j) = energy_sum;
                    else {
                        agl(ai,bi,j) = 0;
                        for (cnt = st; cnt <= ed; cnt++)
                            agl(ai,bi,j) = agl(ai,bi,j) + tempss[cnt]*(cnt+0.5)/energy_sum;
                        /*agl(ai,bi,j) = fmod(agl(ai,bi,j) + sgn*moveStep,Nss[1]);*/
                        agl(ai,bi,j) = agl(ai,bi,j) + sgn*moveStep;
                     /*   if (agl(ai,bi,j)>=Nss[1]) {
                            agl(ai,bi,j) = agl(ai,bi,j) - Nss[1];
                        }
                        else {
                            if (agl(ai,bi,j)<0){
                                agl(ai,bi,j) = agl(ai,bi,j) + Nss[1];
                            }
                        }*/
                        /*compute total energy*/
                        TTEng_2nd(ai,bi,j) = TTEng_1st(ai,bi,j);
                        TTEng_1st(ai,bi,j) = energy_sum;
                        
                        R(ai,bi,j) = 0;
                        /*compute R(ai,bi,j)*/
                        /*st = (int)fmod(st + sgn*moveStep,Nss[1]);*/
                        /*ed = (int)fmod(ed + sgn*moveStep,Nss[1]);*/
                        st = st + sgn*moveStep;
                        if (st>=Nss[1]) {
                            st = st - Nss[1];
                        }
                        else {
                            if (st<0){
                                st = st + Nss[1];
                            }
                        }
                        ed = ed + sgn*moveStep;
                        if (ed>=Nss[1]) {
                            ed = ed - Nss[1];
                        }
                        else {
                            if (ed<0){
                                ed = ed + Nss[1];
                            }
                        }
                        /*compute R*/
                        /*compute temp_R*/
                        if (1) {
                            if (st<=ed) {
                                for (di=st;di<=ed;di++) {
                                    for (ci = 0;ci<Nss[0];ci++) {
                                        if (di==st) {
                                            temp_R[ci] = ss_energy(ci,di,ai,bi);
                                        }
                                        else
                                            temp_R[ci] = temp_R[ci] + ss_energy(ci,di,ai,bi);
                                    }
                                }
                                maxpos = findmax(temp_R,Nss[0]);
                                ed = maxpos;
                                cnt = 1;
                                while (ed<Nss[0]-1 & cnt==1 & temp_R[ed+1]>0 & ed-maxpos < numm_R) {
                                    if (temp_R[ed]>=temp_R[ed+1]) {
                                        ed++;
                                    }
                                    else
                                        cnt = 0;
                                }
                                st = maxpos;
                                cnt = 1;
                                while (st>0 & cnt == 1 & temp_R[st-1]>0 & maxpos-st < numm_R) {
                                    if (temp_R[st] >= temp_R[st-1]) {
                                        st--;
                                    }
                                    else
                                        cnt = 0;
                                }
                                energy_sum = 0;
                                for (cnt = st; cnt <= ed; cnt++) {
                                    energy_sum = energy_sum + temp_R[cnt];
                                }
                                if (energy_sum>0) {
                                    for (cnt = st; cnt <= ed; cnt++) {
                                        R(ai,bi,j) = R(ai,bi,j) + temp_R[cnt]*(cnt)/energy_sum;
                                    }
                                }
                            }
                            else {
                                ed = ed + Nss[1];
                                for (di=st;di<=ed;di++) {
                                    for (ci = 0;ci<Nss[0];ci++) {
                                        if (di==st) {
                                            temp_R[ci] = ss_energy(ci,di,ai,bi);
                                        }
                                        else {
                                            if (di>=Nss[1]) {
                                                temp_R[ci] = temp_R[ci] + ss_energy(ci,di-Nss[1],ai,bi);
                                            }
                                            else
                                                temp_R[ci] = temp_R[ci] + ss_energy(ci,di,ai,bi);
                                        }
                                    }
                                }
                                maxpos = findmax(temp_R,Nss[0]);
                                ed = maxpos;
                                cnt = 1;
                                while (ed<Nss[0]-1 & cnt==1 & temp_R[ed+1]>0 & ed-maxpos < numm_R) {
                                    if (temp_R[ed]>=temp_R[ed+1]) {
                                        ed++;
                                    }
                                    else
                                        cnt = 0;
                                }
                                st = maxpos;
                                cnt = 1;
                                while (st>0 & cnt == 1 & temp_R[st-1]>0 & maxpos-st < numm_R) {
                                    if (temp_R[st] >= temp_R[st-1]) {
                                        st--;
                                    }
                                    else
                                        cnt = 0;
                                }
                                energy_sum = 0;
                                for (cnt = st; cnt <= ed; cnt++) {
                                    energy_sum = energy_sum + temp_R[cnt];
                                }
                                if (energy_sum>0) {
                                    for (cnt = st; cnt <= ed; cnt++) {
                                        R(ai,bi,j) = R(ai,bi,j) + temp_R[cnt]*(cnt)/energy_sum;
                                    }
                                }

                            }
                        }
                       
                    
                    }
                }
                else
                    TTEng_2nd(ai,bi,j) = 0;
                
                /* mexPrintf("%d %d %d\n",ai,bi,j);*/
            }
            for (j=0;j<num_wave;j++) {
                W_sec(ai,bi,j) = TTEng_1st(ai,bi,j)/W_sec(ai,bi,j);
            }
            /*end the loop for j=0,1,2*/
          /*  if (agl(ai,bi,0)<0) {
                tpR = agl(ai,bi,0) + Nss[1];
                agl(ai,bi,0) = agl(ai,bi,1);
                agl(ai,bi,1) = agl(ai,bi,2);
                agl(ai,bi,2) = tpR;
                tpR = R(ai,bi,0);
                R(ai,bi,0) = R(ai,bi,1);
                R(ai,bi,1) = R(ai,bi,2);
                R(ai,bi,2) = tpR;
                tpR = TTEng_1st(ai,bi,0);
                TTEng_1st(ai,bi,0) = TTEng_1st(ai,bi,1);
                TTEng_1st(ai,bi,1) = TTEng_1st(ai,bi,2);
                TTEng_1st(ai,bi,2) = tpR;
                tpR = TTEng_2nd(ai,bi,0);
                TTEng_2nd(ai,bi,0) = TTEng_2nd(ai,bi,1);
                TTEng_2nd(ai,bi,1) = TTEng_2nd(ai,bi,2);
                TTEng_2nd(ai,bi,2) = tpR;
                tpR = W_sec(ai,bi,0);
                W_sec(ai,bi,0) = W_sec(ai,bi,1);
                W_sec(ai,bi,1) = W_sec(ai,bi,2);
                W_sec(ai,bi,2) = tpR;

            }
            else if (agl(ai,bi,2)>=Nss[1]) {
                tpR = agl(ai,bi,2) - Nss[1];
                agl(ai,bi,2) = agl(ai,bi,1);
                agl(ai,bi,1) = agl(ai,bi,0);
                agl(ai,bi,0) = tpR;
                tpR = R(ai,bi,2);
                R(ai,bi,2) = R(ai,bi,1);
                R(ai,bi,1) = R(ai,bi,0);
                R(ai,bi,0) = tpR;
                tpR = TTEng_1st(ai,bi,2);
                TTEng_1st(ai,bi,2) = TTEng_1st(ai,bi,1);
                TTEng_1st(ai,bi,1) = TTEng_1st(ai,bi,0);
                TTEng_1st(ai,bi,0) = tpR;
                tpR = TTEng_2nd(ai,bi,2);
                TTEng_2nd(ai,bi,2) = TTEng_2nd(ai,bi,1);
                TTEng_2nd(ai,bi,1) = TTEng_2nd(ai,bi,0);
                TTEng_2nd(ai,bi,0) = tpR;
                tpR = W_sec(ai,bi,2);
                W_sec(ai,bi,2) = W_sec(ai,bi,1);
                W_sec(ai,bi,1) = W_sec(ai,bi,0);
                W_sec(ai,bi,0) = tpR;
            }*/
            /*mexPrintf("%f %f %f\n",agl(ai,bi,0),agl(ai,bi,1),agl(ai,bi,2));*/
        }
    }
    
    mxFree(temp);
    mxFree(temp2);
    mxFree(tempss);
return;
}
