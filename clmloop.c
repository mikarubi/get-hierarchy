#include "mex.h"        /* for mex   */
#include "time.h"       /* for qsort */

/*declare rand_mt functions*/
void init_genrand(unsigned long s);
double genrand_real2(void);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxArray *mxB, *mxHnm, *mxMb;
    double    *B,   *Hnm, dq, max_dq;
    uint32_T  *Mb, *U, ma, mb, mb_max, u, v,  i, j, n, flag;
    
    static int init_mt;
    if(!init_mt)
        init_mt=1, init_genrand(time(NULL));        /* initialize Mersenne Twister */
    
    if(nrhs || nlhs)
        mexErrMsgTxt("There are no inputs or outputs.");
    
    mxB   = mexGetVariable("caller", "B"  );
    mxHnm = mexGetVariable("caller", "Hnm");
    mxMb  = mexGetVariable("caller", "Mb" );
    
    if(!mxIsDouble(mxB))   mexErrMsgTxt("B is not in double format.");
    if(!mxIsDouble(mxHnm)) mexErrMsgTxt("Hnm is not in double format.");
    if(!mxIsUint32(mxMb))  mexErrMsgTxt("mxMb is not in uint32 format.");
    
    B     =  mxGetPr(mxB);
    Hnm   =  mxGetPr(mxHnm);
    Mb    =  (uint32_T *)mxGetData(mxMb);
    n     =  mxGetM(mxB);
    
    /* check community vector and convert to C indexing */
    for(i=0; i<n; i++){
        if(Mb[i]==0)
            mexErrMsgTxt("Mb must containt natural numbers");
        else
            Mb[i] -= 1;
    }
    
    U     = (uint32_T*)mxMalloc(n*sizeof(uint32_T));
    do{
        /* "inside-out" Fisher-Yates shuffle */
        for(i=0; i<n; i++){     /* numbers [0,i-1] are already shuffled at positions [0,i-1] */
            j = i*genrand_real2();      /* j is ma random position in [0,i-1] */
            if(j!=i)
                U[i] = U[j];            /* reassign number from position j to position i */
            U[j] = i;                   /* newly assign number i to position j */
        }
        
        flag = 0;
        for(i=0; i<n; i++){           	/* loop over all nodes in random order */
            u  = U[i];                 	/* node u */
            ma = Mb[u];                 /* current module ma of node u */
            
            max_dq = 0;
            for(mb=0; mb<n; mb++){     	/* loop over all modules */
                if(ma!=mb){
                    dq = Hnm[u + mb*n] - Hnm[u + ma*n] + B[u + u*n];
                    if(dq > max_dq + 1e-10){
                        max_dq = dq;
                        mb_max = mb;
                    }
                }
            }
            
            if(max_dq > 0){             /* if maximal increase is positive */
                flag = 1;
                mb = mb_max;
                Mb[u] = mb;             /* reassign module */
                
                for(v=0; v<n; v++){    	/* change node-to-module strengths */
                    Hnm[v + mb*n] += B[v + u*n];
                    Hnm[v + ma*n] -= B[v + u*n];
                }
            }
        }
    }while(flag);
    
    mexPutVariable("caller", "B",   mxB  );
    mexPutVariable("caller", "Hnm", mxHnm);
    mexPutVariable("caller", "Mb",  mxMb );
    
    mxFree(U);
}

