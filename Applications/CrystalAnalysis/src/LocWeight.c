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


/*This function computes the local wave vectors via weighted average.
 
 By Haizhao Yang and Jianfeng Lu
 */
void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
#define agl(ai,bi,j) agl[ai+dims[0]*(bi+dims[1]*j)]
#define R(ai,bi,j) R[ai+dims[0]*(bi+dims[1]*j)]
#define TTEng_1st(ai,bi,j) TTEng_1st[ai+dims[0]*(bi+dims[1]*j)]
#define TTEng_2nd(ai,bi,j) TTEng_2nd[ai+dims[0]*(bi+dims[1]*j)]
#define W_sec(ai,bi,j) W_sec[ai+dims[0]*(bi+dims[1]*j)]
#define ss_energy(ci,di,ai,bi) ss_energy[ci+Nss[0]*(di+Nss[1]*(ai+Nss[2]*bi))]
#define ss_avgdx(ci,di,ai,bi) ss_avgdx[ci+Nss[0]*(di+Nss[1]*(ai+Nss[2]*bi))]
#define ss_avgdy(ci,di,ai,bi) ss_avgdy[ci+Nss[0]*(di+Nss[1]*(ai+Nss[2]*bi))]
#define pi 3.141592653589793
    
    size_t ai, bi, ci, di, k, j, cnt, cnt2;
    int num_wave, pos_temp;
    double *ss_energy, *ss_avgdx, *ss_avgdy, *agl, *R, *TTEng_1st, *TTEng_2nd, *W_sec; /* R2, agl2 */
    double temp_avgdx, temp_avgdy, sumEng;
    ss_energy = mxGetPr(prhs[0]);
    const mwSize *Nss = mxGetDimensions(prhs[0]);
    ss_avgdx = mxGetPr(prhs[1]);
    ss_avgdy = mxGetPr(prhs[2]);
    num_wave = mxGetScalar(prhs[3]);
    nrhs = 4;
    
    nlhs = 5;
    int ndim = 3, dims[3] = {Nss[2],Nss[3],num_wave}, numm = (int)Nss[1]/num_wave/8, numm_R = (int)floor(Nss[0]/4);
    plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[3] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[4] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    /* plhs[5] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL); */
    /* plhs[6] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL); */
    agl = mxGetPr(plhs[0]);
    R = mxGetPr(plhs[1]);
    TTEng_1st = mxGetPr(plhs[2]);
    TTEng_2nd = mxGetPr(plhs[3]);
    W_sec = mxGetPr(plhs[4]);
    /* agl2 = mxGetPr(plhs[5]); */
    /* R2 = mxGetPr(plhs[6]); */
    
    int L = round(Nss[1]/num_wave), maxpos, st, ed, stR, edR, moveStep, pos, sgn, temp_di, markPos;
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
            markPos = fmod(floor(maxpos/L),num_wave);
            maxpos = maxpos - markPos*L;
            moveStep = floor(maxpos-L/2);
            if (moveStep>=0) {
                sgn = 1;
                for (cnt=0;cnt<Nss[1]-moveStep;cnt++) {
                    tempss[cnt] = temp2[cnt+moveStep];
                }
                for (cnt=Nss[1]-moveStep;cnt<Nss[1];cnt++) {
                    tempss[cnt] = temp2[cnt -(Nss[1]-moveStep)];
                }            }
            if (moveStep<0) {
                moveStep = -moveStep;
                sgn = -1;
                for (cnt = 0; cnt < Nss[1]-moveStep; cnt++) {
                    tempss[cnt+moveStep] = temp2[cnt];
                }
                for (cnt = Nss[1]-moveStep; cnt < Nss[1]; cnt++) {
                    tempss[cnt-(Nss[1]-moveStep)] = temp2[cnt];
                }            }
            for (j=0;j<num_wave;j++) {
                for (cnt=0;cnt<L;cnt++) {
                    W_sec(ai,bi,j) = W_sec(ai,bi,j) + tempss[cnt+j*L];
                }
            }
            /*tempss[:] computed*/
            for (j=0;j<num_wave;j++) {
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
                energy_sum = 0;
                for (cnt = st; cnt <= ed; cnt++)
                    energy_sum = energy_sum + tempss[cnt];
                /*compute TTEng_1st(ai,bi,j)*/
                TTEng_1st(ai,bi,j) = energy_sum;
                
                /* /\*compute agl(ai,bi,j)*\/ */
                 for (cnt = st; cnt <= ed; cnt++) {
                /*     agl(ai,bi,j) = agl(ai,bi,j) + tempss[cnt]*cnt/energy_sum; */
                     tempss[cnt] = 0;
                 } 
                /* agl(ai,bi,j) = agl(ai,bi,j) + sgn*moveStep; */
                /*compute R(ai,bi,j)*/
                /* change (st,ed) back to original coordinate */
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
                int st_agl = st;
                int ed_agl = ed;
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
                        
                        /*
                         energy_sum = 0;
                         for (cnt = st; cnt <= ed; cnt++) {
                         energy_sum = energy_sum + temp_R[cnt];
                         }
                         if (energy_sum>0) {
                         for (cnt = st; cnt <= ed; cnt++) {
                         R(ai,bi,j) = R(ai,bi,j) + temp_R[cnt]*(cnt)/energy_sum;
                         }
                         } */
                        sumEng = 0;
                        temp_avgdx = 0; temp_avgdy = 0;
                        for (di=st_agl;di<=ed_agl;di++) {
                            for (ci=st;ci<=ed;ci++) {
                                temp_avgdx += ss_avgdx(ci,di,ai,bi);
                                temp_avgdy += ss_avgdy(ci,di,ai,bi);
                                sumEng += ss_energy(ci,di,ai,bi);
                            }
                        }
                        if (sumEng==0) {
                            agl(ai,bi,j) = 0;
                            R(ai,bi,j) = 0;
                            /*mexPrintf("%f %f, ",sumEng,TTEng_1st(ai,bi,j));*/
                        }else {
                            
                            temp_avgdx /= sumEng;
                            temp_avgdy /= sumEng;
                            R(ai,bi,j) = sqrt(temp_avgdx*temp_avgdx + temp_avgdy*temp_avgdy);
                            if (temp_avgdy>=0) {
                                agl(ai,bi,j) = acos(temp_avgdx/R(ai,bi,j));
                            }
                            else
                                agl(ai,bi,j) = pi - acos(temp_avgdx/R(ai,bi,j));
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
                        /*
                         energy_sum = 0;
                         for (cnt = st; cnt <= ed; cnt++) {
                         energy_sum = energy_sum + temp_R[cnt];
                         }
                         if (energy_sum>0) {
                         for (cnt = st; cnt <= ed; cnt++) {
                         R(ai,bi,j) = R(ai,bi,j) + temp_R[cnt]*(cnt)/energy_sum;
                         }
                         } */
                        ed_agl = ed_agl + Nss[1];
                        sumEng = 0;
                        temp_avgdx = 0; temp_avgdy = 0;
                        
                        for (di=st_agl;di<=ed_agl;di++) {
                            if (di>=Nss[1]) {
                                for (ci=st;ci<=ed;ci++) {
                                    temp_avgdx -= ss_avgdx(ci,di-Nss[1],ai,bi); /* "-" due to shifting to negative half plane */
                                    temp_avgdy -= ss_avgdy(ci,di-Nss[1],ai,bi);
                                    sumEng += ss_energy(ci,di-Nss[1],ai,bi);
                                }
                            }
                            else {
                                for (ci=st;ci<=ed;ci++) {
                                    temp_avgdx += ss_avgdx(ci,di,ai,bi);
                                    temp_avgdy += ss_avgdy(ci,di,ai,bi);
                                    sumEng += ss_energy(ci,di,ai,bi);
                                }
                            }
                        }
                        if (sumEng==0) {
                            agl(ai,bi,j) = 0;
                            R(ai,bi,j) = 0;
                            /*mexPrintf("%f %f, ",sumEng,TTEng_1st(ai,bi,j));*/
                        }else {
                            
                            temp_avgdx /= sumEng;
                            temp_avgdy /= sumEng;
                            R(ai,bi,j) = sqrt(temp_avgdx*temp_avgdx + temp_avgdy*temp_avgdy);
                            if (temp_avgdy>=0) {
                                agl(ai,bi,j) = acos(temp_avgdx/R(ai,bi,j));
                            }
                            else
                                agl(ai,bi,j) = pi - acos(temp_avgdx/R(ai,bi,j));
                        }
                    }
                }
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
                        /* agl2(ai,bi,j) = 0;
                         for (cnt = st; cnt <= ed; cnt++)
                         agl(ai,bi,j) = agl(ai,bi,j) + tempss[cnt]*cnt/energy_sum;
                         agl(ai,bi,j) = agl(ai,bi,j) + sgn*moveStep; */
                        /*compute total energy*/
                        TTEng_2nd(ai,bi,j) = TTEng_1st(ai,bi,j);
                        TTEng_1st(ai,bi,j) = energy_sum;
                        
                        R(ai,bi,j) = 0;
                        /* R2(ai,bi,j) = 0; */
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
                        int st_agl = st;
                        int ed_agl = ed;
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
                                /*
                                 energy_sum = 0;
                                 for (cnt = st; cnt <= ed; cnt++) {
                                 energy_sum = energy_sum + temp_R[cnt];
                                 }
                                 if (energy_sum>0) {
                                 for (cnt = st; cnt <= ed; cnt++) {
                                 R(ai,bi,j) = R(ai,bi,j) + temp_R[cnt]*(cnt)/energy_sum;
                                 }
                                 } */
                                sumEng = 0;
                                temp_avgdx = 0; temp_avgdy = 0;
                                
                                for (di=st_agl;di<=ed_agl;di++) {
                                    for (ci=st;ci<=ed;ci++) {
                                        temp_avgdx += ss_avgdx(ci,di,ai,bi);
                                        temp_avgdy += ss_avgdy(ci,di,ai,bi);
                                        sumEng += ss_energy(ci,di,ai,bi);
                                    }
                                }
                                temp_avgdx /= sumEng;
                                temp_avgdy /= sumEng;
                                
                                R(ai,bi,j) = sqrt(temp_avgdx*temp_avgdx + temp_avgdy*temp_avgdy);
                                if (temp_avgdy>=0) {
                                    agl(ai,bi,j) = acos(temp_avgdx/R(ai,bi,j));
                                }
                                else
                                    agl(ai,bi,j) = pi - acos(temp_avgdx/R(ai,bi,j));
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
                                /*
                                 energy_sum = 0;
                                 for (cnt = st; cnt <= ed; cnt++) {
                                 energy_sum = energy_sum + temp_R[cnt];
                                 }
                                 if (energy_sum>0) {
                                 for (cnt = st; cnt <= ed; cnt++) {
                                 R(ai,bi,j) = R(ai,bi,j) + temp_R[cnt]*(cnt)/energy_sum;
                                 }
                                 } */
                                ed_agl = ed_agl + Nss[1];
                                sumEng = 0;
                                temp_avgdx = 0; temp_avgdy = 0;
                                
                                for (di=st_agl;di<=ed_agl;di++) {
                                    if (di>=Nss[1]) {
                                        for (ci=st;ci<=ed;ci++) {
                                            temp_avgdx -= ss_avgdx(ci,di-Nss[1],ai,bi);
                                            
                                            temp_avgdy -= ss_avgdy(ci,di-Nss[1],ai,bi);
                                            sumEng += ss_energy(ci,di-Nss[1],ai,bi);
                                        }
                                    }
                                    else {
                                        for (ci=st;ci<=ed;ci++) {
                                            temp_avgdx += ss_avgdx(ci,di,ai,bi);
                                            temp_avgdy += ss_avgdy(ci,di,ai,bi);
                                            sumEng += ss_energy(ci,di,ai,bi);
                                        }
                                    }
                                }
                                temp_avgdx /= sumEng;
                                temp_avgdy /= sumEng;
                                
                                R(ai,bi,j) = sqrt(temp_avgdx*temp_avgdx + temp_avgdy*temp_avgdy);
                                if (temp_avgdy>=0) {
                                    agl(ai,bi,j) = acos(temp_avgdx/R(ai,bi,j));
                                }
                                else
                                    agl(ai,bi,j) = pi - acos(temp_avgdx/R(ai,bi,j));
                                
                                
                                
                                
                            }
                        }
                        
                        
                    }
                }
                else
                    TTEng_2nd(ai,bi,j) = 0;
                
            }
            for (j=0;j<num_wave;j++) {
                W_sec(ai,bi,j) = TTEng_1st(ai,bi,j)/W_sec(ai,bi,j);
                if (R(ai,bi,j)==0) {
                    R(ai,bi,j) = R(ai,bi,markPos);
                    agl(ai,bi,j) = agl(ai,bi,markPos) + (j-markPos)*pi/3;
                }
            }
        }
    }
    
    mxFree(temp);
    mxFree(temp2);
    mxFree(tempss);
    mxFree(temp_R);
    return;
}
