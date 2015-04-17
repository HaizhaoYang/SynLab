/* This code computes the synchrosqueezed transform in a polar coordinate.
 by Haizhao Yang and Jianfeng Lu */

#include "mex.h" /* Always include this */
#include <math.h>
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
                 int nrhs, const mxArray *prhs[]) /* Input variables */
{
#define ccr(ai,bi) (ccr[ai+Nss[0]*bi])
#define cci(ai,bi) (cci[ai+Nss[0]*bi])
#define kk1(ai,bi) (kk1[ai+Nss[0]*bi])
#define kk2(ai,bi) (kk2[ai+Nss[0]*bi])
#define kb(loc1,loc2,ai,bi) kb[loc1+NB1*(loc2+NB2*(ai+Nss[0]*bi))]
#define avgdx(loc1,loc2,ai,bi) avgdx[loc1+NB1*(loc2+NB2*(ai+Nss[0]*bi))]
#define avgdy(loc1,loc2,ai,bi) avgdy[loc1+NB1*(loc2+NB2*(ai+Nss[0]*bi))]
    
    size_t ai, bi;
    int NB1, NB2, loc1, loc2, di = 0;
    double *kk1, *kk2, *ccr, *cci, *kb, *avgdx, *avgdy;
    double EXT, num_dir, da, dr, r, agl, R_low;
    double temp_energy;
    ccr = mxGetPr(prhs[0]);
    cci = mxGetPi(prhs[0]);
    kk1 = mxGetPr(prhs[1]);
    kk2 = mxGetPr(prhs[2]);
    EXT = mxGetScalar(prhs[3]);
    num_dir = mxGetScalar(prhs[4]);
    da = mxGetScalar(prhs[5]);
    dr = mxGetScalar(prhs[6]);
    NB1 = mxGetScalar(prhs[7]);
    NB2 = mxGetScalar(prhs[8]);
    R_low = mxGetScalar(prhs[9]);
    const mwSize *Nss = mxGetDimensions(prhs[0]);
    kb = mxGetPr(prhs[10]);
    avgdx = mxGetPr(prhs[11]);
    avgdy = mxGetPr(prhs[12]);
    nrhs = 13;
    
    nlhs = 0;
    
    for (ai=0;ai<Nss[0];ai++) {
        for (bi=0;bi<Nss[1];bi++) {
            if (kk1(ai,bi)<EXT) {
                r = sqrt(kk1(ai,bi)*kk1(ai,bi)+kk2(ai,bi)*kk2(ai,bi));
                if (kk1(ai,bi)>=0) {
                    agl = fmod(acos(kk2(ai,bi)/r),num_dir);
                }
                else
                    agl = fmod(3.1415926-acos(kk2(ai,bi)/r),num_dir);
                loc1 = round((r-R_low)/dr);
                loc2 = round(agl/da);
                
                temp_energy = ccr(ai,bi)*ccr(ai,bi) + cci(ai,bi)*cci(ai,bi);
                kb(loc1,loc2,ai,bi) = kb(loc1,loc2,ai,bi) + temp_energy;
                avgdx(loc1,loc2,ai,bi) = avgdx(loc1,loc2,ai,bi) + r * cos(agl) * temp_energy;
                avgdy(loc1,loc2,ai,bi) = avgdy(loc1,loc2,ai,bi) + r * sin(agl) * temp_energy;
                /* cannot use kk since those are symmetric */
            }
        }
    }
    return;
}
