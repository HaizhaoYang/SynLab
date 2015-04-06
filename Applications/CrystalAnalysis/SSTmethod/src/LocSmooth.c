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
/*This function uses the symmetric property and computes the agl in one direction*/
void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
int nrhs, const mxArray *prhs[]) /* Input variables */
{
    #define agl(ai,bi) agl[ai+dims[0]*bi]
    #define TTEng_1st(ai,bi) TTEng_1st[ai+dims[0]*bi]
    #define TTEng_2nd(ai,bi) TTEng_2nd[ai+dims[0]*bi]
    #define ss_energy(ci,di,ai,bi) ss_energy[ci+Nss[0]*(di+Nss[1]*(ai+Nss[2]*bi))]
    
    size_t ai, bi, ci, di, k, j, cnt;
    int num_wave;
    double *ss_energy, *agl, *TTEng_1st, *TTEng_2nd;
    ss_energy = mxGetPr(prhs[0]);
    const mwSize *Nss = mxGetDimensions(prhs[0]);
    num_wave = mxGetScalar(prhs[1]);
    nrhs = 2;
    
    nlhs = 3;
    int ndim = 2, dims[2] = {Nss[2],Nss[3]}, numm = (int)(Nss[1]/num_wave/4);
    plhs[0] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    agl = mxGetPr(plhs[0]);
    TTEng_1st = mxGetPr(plhs[1]);
    TTEng_2nd = mxGetPr(plhs[2]);
    
    int L = Nss[1]/num_wave, maxpos, st, ed;
    double *temp, energy_sum, *tempss;
    temp = (double *)mxMalloc(sizeof(double)*L);
    tempss = (double *)mxMalloc(sizeof(double)*(2*numm+1));
    for (ai=0;ai<Nss[2];ai++) {
        for (bi=0;bi<Nss[3];bi++) {
            for (k=0;k<L;k++) {
                for (j=0;j<num_wave;j++) {
                    di = j*L + k;
                    for (ci=0;ci<Nss[0];ci++) {
                        if (ci==0 & j==0)
                            temp[k] = ss_energy(ci,di,ai,bi);
                        else
                            temp[k] = temp[k] + ss_energy(ci,di,ai,bi);
                    }
                }
            }
            /*temp[:] computed*/
            maxpos = findmax(temp,L);
            for (cnt=0;cnt<2*numm+1;cnt++) {
                tempss[cnt] = temp[(int)fmod(cnt+maxpos-numm,L)];
            }
            ed = numm;
            cnt = 1;
            while (ed<2*numm & cnt==1 & tempss[ed+1]>0) {
                if (tempss[ed]>=tempss[ed+1])
                    ed++;
                else
                    cnt = 0;
            }
            st = numm;
            cnt = 1;
            while (st>0 & cnt == 1 & tempss[st-1]>0) {
                if (tempss[st] >= tempss[st-1])
                    st--;
                else
                    cnt = 0;
            }
            /*compute agl(ai,bi)*/
            energy_sum = 0;
            for (cnt = st; cnt <= ed; cnt++)
                energy_sum = energy_sum + tempss[cnt];
            for (cnt = st; cnt <= ed; cnt++) {
                agl(ai,bi) = agl(ai,bi) + tempss[cnt]*(cnt+0.5)/energy_sum;
                temp[(int)fmod(cnt+maxpos-numm,L)] = 0;
            }
            agl(ai,bi) = fmod(agl(ai,bi)+maxpos-numm,L);/*any problem here??????*/
            
            /*compute TTEng_1st(ai,bi)*/
            TTEng_1st(ai,bi) = energy_sum;
            /*compute TTEng_2nd(ai,bi)*/
            maxpos = findmax(temp,L);
            if (temp[maxpos]>0) {
                for (cnt=0;cnt<2*numm+1;cnt++) {
                    tempss[cnt] = temp[(int)fmod(cnt+maxpos-numm,L)];
                }
                ed = numm;
                cnt = 1;
                while (ed<2*numm & cnt==1 & tempss[ed+1]>0) {
                    if (tempss[ed]>=tempss[ed+1])
                        ed++;
                    else
                        cnt = 0;
                }
                st = numm;
                cnt = 1;
                while (st>0 & cnt == 1 & tempss[st-1]>0) {
                    if (tempss[st] >= tempss[st-1])
                        st--;
                    else
                        cnt = 0;
                }
                /*compute agl(ai,bi)*/
                energy_sum = 0;
                for (cnt = st; cnt <= ed; cnt++)
                    energy_sum = energy_sum + tempss[cnt];
                if (energy_sum<=TTEng_1st(ai,bi))
                    TTEng_2nd(ai,bi) = energy_sum;
                else {
                    agl(ai,bi) = 0;
                    for (cnt = st; cnt <= ed; cnt++) {
                        agl(ai,bi) = agl(ai,bi) + tempss[cnt]*(cnt+0.5)/energy_sum;
                    }
                    agl(ai,bi) = fmod(agl(ai,bi)+maxpos-numm,L);
                    TTEng_2nd(ai,bi) = TTEng_1st(ai,bi);
                    TTEng_1st(ai,bi) = energy_sum;
                }
            }
            else
                TTEng_2nd(ai,bi) = 0;

        }
    }
    
    mxFree(temp);
    mxFree(tempss);
return;
}