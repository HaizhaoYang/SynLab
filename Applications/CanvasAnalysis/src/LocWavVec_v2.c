/* This code estimates local wave vectors in a polar coordinate
 By Haizhao Yang */

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


/* This function computes two kinds of  angles and radii */
void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
#define agl(ai,bi) agl[ai+dims[0]*bi]
#define R(ai,bi) R[ai+dims[0]*bi]
#define ss_energy(ci,di,ai,bi) ss_energy[ci+Nss[0]*(di+Nss[1]*(ai+Nss[2]*bi))]
#define ss_avgdx(ci,di,ai,bi) ss_avgdx[ci+Nss[0]*(di+Nss[1]*(ai+Nss[2]*bi))]
#define ss_avgdy(ci,di,ai,bi) ss_avgdy[ci+Nss[0]*(di+Nss[1]*(ai+Nss[2]*bi))]
#define agl2(ai,bi) agl2[ai+dims[0]*bi]
#define R2(ai,bi) R2[ai+dims[0]*bi]
#define pi 3.141592653589793
    
    size_t ai, bi, ci, di, k, cnt, cnt2;
    int pos_st, pos_ed, pos_temp;
    double *ss_energy, *ss_avgdx, *ss_avgdy, *agl, *agl2, *R, *R2;
    double temp_avgdx, temp_avgdy, sumEng;
    ss_energy = mxGetPr(prhs[0]);
    const mwSize *Nss = mxGetDimensions(prhs[0]);
    ss_avgdx = mxGetPr(prhs[1]);
    ss_avgdy = mxGetPr(prhs[2]);
    pos_st = mxGetScalar(prhs[3])-1;
    pos_ed = mxGetScalar(prhs[4])-1;
    nrhs = 5;
    
    nlhs = 2;
    int ndim = 2, dims[2] = {Nss[2],Nss[3]}, numm_R = (int)floor(Nss[0]/4);
    plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
   /* plhs[2] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[3] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);*/
    /*agl = mxGetPr(plhs[0]);
    R = mxGetPr(plhs[1]);*/
    agl2 = mxGetPr(plhs[0]);
    R2 = mxGetPr(plhs[1]);
    
    int maxpos, st, ed, stR, edR, moveStep, pos, sgn, temp_di;
    double *temp, energy_sum, *temp_R, tpR;
    
    temp_R = (double *)mxMalloc(sizeof(double)*Nss[0]);
    if (pos_st<pos_ed) {
        int L_temp = pos_ed-pos_st+1;
        temp = (double *)mxMalloc(sizeof(double)*L_temp);
        for (ai=0;ai<Nss[2];ai++) {
            for (bi=0;bi<Nss[3];bi++) {
                /*compute temp*/
                for (di=pos_st;di<=pos_ed;di++) {
                    for (ci=0;ci<Nss[0];ci++) {
                        if (ci==0)
                            temp[di-pos_st] = ss_energy(ci,di,ai,bi);
                        else
                            temp[di-pos_st] = temp[di-pos_st] + ss_energy(ci,di,ai,bi);
                    }
                }
                /*temp[:] computed*/
                maxpos = findmax2(temp,0,L_temp);
                ed = maxpos;
                cnt = 1;
                while (ed<L_temp-1 & cnt==1 & temp[ed+1]>0) {
                    if (temp[ed]>=temp[ed+1])
                        ed++;
                    else
                        cnt = 0;
                }
                st = maxpos;
                cnt = 1;
                while (st>0 & cnt == 1 & temp[st-1]>0) {
                    if (temp[st] >= temp[st-1])
                        st--;
                    else
                        cnt = 0;
                }
                /*compute agl(ai,bi)*/
              /*  energy_sum = 0;
                for (cnt = st; cnt <= ed; cnt++)*/
                    energy_sum = energy_sum + temp[cnt];
                for (cnt = st; cnt <= ed; cnt++) {
                   /* agl(ai,bi) = agl(ai,bi) + temp[cnt]*cnt/energy_sum;*/
                    temp[cnt] = 0;
                }
               /* agl(ai,bi) += pos_st;*/
                st += pos_st;
                ed += pos_st;
                int st_agl = st;
                int ed_agl = ed;
                /*compute R(ai,bi)*/
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
               /* energy_sum = 0;
                for (cnt = st; cnt <= ed; cnt++) {
                    energy_sum = energy_sum + temp_R[cnt];
                }
                if (energy_sum>0) {
                    for (cnt = st; cnt <= ed; cnt++) {
                        R(ai,bi) = R(ai,bi) + temp_R[cnt]*cnt/energy_sum;
                    }
                }*/
                sumEng = 0;
                temp_avgdx = 0; temp_avgdy = 0;
                for (di=st_agl;di<=ed_agl;di++) {
                    for (ci=st;ci<=ed;ci++) {
                        temp_avgdx += ss_avgdx(ci,di,ai,bi);/**ss_energy(ci,di,ai,bi);*/
                        temp_avgdy += ss_avgdy(ci,di,ai,bi);/**ss_energy(ci,di,ai,bi);   //hz: we don't use ./ss_energy(ci,di,ai,bi) in SS_polar.c, so we don't need this step */
                        sumEng += ss_energy(ci,di,ai,bi);
                    }
                }
                
                temp_avgdx /= sumEng;
                temp_avgdy /= sumEng;
                R2(ai,bi) = sqrt(temp_avgdx*temp_avgdx + temp_avgdy*temp_avgdy);
                if (temp_avgdy>=0) {
                    agl2(ai,bi) = acos(temp_avgdx/R2(ai,bi));
                }
                else
                    agl2(ai,bi) = pi - acos(temp_avgdx/R2(ai,bi));
            }
            /*mexPrintf("%f %f %f\n",agl(ai,bi,0),agl(ai,bi,1),agl(ai,bi,2));*/
        }
    } else {
        int L_temp = pos_ed-pos_st+1+Nss[1];
        temp = (double *)mxMalloc(sizeof(double)*L_temp);
        for (ai=0;ai<Nss[2];ai++) {
            for (bi=0;bi<Nss[3];bi++) {
                /*compute temp*/
                for (di=pos_st;di<=pos_ed+Nss[1];di++) {
                    for (ci=0;ci<Nss[0];ci++) {
                        if (di<Nss[1]) {
                            if (ci==0)
                                temp[di-pos_st] = ss_energy(ci,di,ai,bi);
                            else
                                temp[di-pos_st] = temp[di-pos_st] + ss_energy(ci,di,ai,bi);
                        } else {
                            if (ci==0)
                                temp[di-pos_st] = ss_energy(ci,di-Nss[1],ai,bi);
                            else
                                temp[di-pos_st] = temp[di-pos_st] + ss_energy(ci,di-Nss[1],ai,bi);
                        }
                    }
                }
                /*temp[:] computed*/
                maxpos = findmax2(temp,0,L_temp);
                ed = maxpos;
                cnt = 1;
                while (ed<L_temp-1 & cnt==1 & temp[ed+1]>0) {
                    if (temp[ed]>=temp[ed+1])
                        ed++;
                    else
                        cnt = 0;
                }
                st = maxpos;
                cnt = 1;
                while (st>0 & cnt == 1 & temp[st-1]>0) {
                    if (temp[st] >= temp[st-1])
                        st--;
                    else
                        cnt = 0;
                }
                /*compute agl(ai,bi)*/
              /*  energy_sum = 0;
                for (cnt = st; cnt <= ed; cnt++)
                    energy_sum = energy_sum + temp[cnt];*/
                for (cnt = st; cnt <= ed; cnt++) {
                 /*   agl(ai,bi) = agl(ai,bi) + temp[cnt]*cnt/energy_sum;*/
                    temp[cnt] = 0;
                }
               /* agl(ai,bi) += pos_st;
                if (agl(ai,bi)>=Nss[1]) {
                    agl(ai,bi) -= Nss[1];
                }*/
                st += pos_st;
                ed += pos_st;
                int st_agl = st;
                int ed_agl = ed;
                /*compute R(ai,bi)*/
                for (di=st;di<=ed;di++) {
                    for (ci = 0;ci<Nss[0];ci++) {
                        if (di<Nss[1]) {
                            if (di==st) {
                                temp_R[ci] = ss_energy(ci,di,ai,bi);
                            }
                            else
                                temp_R[ci] = temp_R[ci] + ss_energy(ci,di,ai,bi);
                        }else {
                            if (di==st) {
                                temp_R[ci] = ss_energy(ci,di-Nss[1],ai,bi);
                            }
                            else
                                temp_R[ci] = temp_R[ci] + ss_energy(ci,di-Nss[1],ai,bi);
                        }
                    }
                }
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
             /*   energy_sum = 0;
                for (cnt = st; cnt <= ed; cnt++) {
                    energy_sum = energy_sum + temp_R[cnt];
                }
                if (energy_sum>0) {
                    for (cnt = st; cnt <= ed; cnt++) {
                        R(ai,bi) = R(ai,bi) + temp_R[cnt]*cnt/energy_sum;
                    }
                }*/
                sumEng = 0;
                temp_avgdx = 0; temp_avgdy = 0;
                for (di=st_agl;di<=ed_agl;di++) {
                    for (ci=st;ci<=ed;ci++) {
                        if (di<Nss[1]) {
                            temp_avgdx += ss_avgdx(ci,di,ai,bi);/**ss_energy(ci,di,ai,bi);*/
                            temp_avgdy += ss_avgdy(ci,di,ai,bi);/**ss_energy(ci,di,ai,bi);   //hz: we don't use ./ss_energy(ci,di,ai,bi) in SS_polar.c, so we don't need this step */
                            sumEng += ss_energy(ci,di,ai,bi);
                        }else {
                            temp_avgdx -= ss_avgdx(ci,di-Nss[1],ai,bi);/**ss_energy(ci,di,ai,bi);*/
                            temp_avgdy -= ss_avgdy(ci,di-Nss[1],ai,bi);/**ss_energy(ci,di,ai,bi);   //hz: we don't use ./ss_energy(ci,di,ai,bi) in SS_polar.c, so we don't need this step */
                            sumEng += ss_energy(ci,di-Nss[1],ai,bi);
                        }
                    }
                }
                
                temp_avgdx /= sumEng;
                temp_avgdy /= sumEng;
                R2(ai,bi) = sqrt(temp_avgdx*temp_avgdx + temp_avgdy*temp_avgdy);
                if (temp_avgdy>=0) {
                    agl2(ai,bi) = acos(temp_avgdx/R2(ai,bi));
                }
                else
                    agl2(ai,bi) = pi - acos(temp_avgdx/R2(ai,bi));
            }
            /*mexPrintf("%f %f %f\n",agl(ai,bi,0),agl(ai,bi,1),agl(ai,bi,2));*/
        }
    }
return;
}
