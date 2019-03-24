#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    uint32_T *X;            /* input  */
    double *F;              /* output */
    int i, n, m;
    
    /* check number of variables */
    if(nrhs!=1) mexErrMsgTxt("There is one input variable.");
    if(nlhs >1) mexErrMsgTxt("There is one output variable.");
    
    /* initialize and verify input */
    if(!mxIsUint32(prhs[0])) mexErrMsgTxt("Input is not in uint32 format.");
    X = (uint32_T *)mxGetData(prhs[0]);
    n = mxGetNumberOfElements(prhs[0]);
    
    /* get maximum value of input */
    for(m=0, i=0; i<n; i++)
        if(X[i] > m)
            m = X[i];
    
    /* inialize output */
    plhs[0] = mxCreateDoubleMatrix(m+1, 1, mxREAL);
    F = mxGetPr(plhs[0]);
    
    /* get frequency distribution */
    for(i=0; i<n; i++)
        F[X[i]] += 1;
    
}
