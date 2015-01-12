#include "mex.h"
#include <process.h>

extern "C" void gpuAdd(double* A, double* B, double *C);


void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
    int m,n;
    double *C, *A, *B;
    //if (mxGetNumberOfElements(prhs[0])!=5) mexErrMsgTxt("Wrong number of elements in A!");
    //if (mxGetNumberOfElements(prhs[1])!=5) mexErrMsgTxt("Wrong number of elements in B!");
    A = (double*)mxGetData(prhs[0]);
    B = (double*)mxGetData(prhs[1]);
    C=(double*)mxGetData(plhs[0]=mxCreateDoubleMatrix(1,5,mxREAL));
    for(m=0;m<1;m++){
        C[m] = A[m] + B[m];
        A[m]=0;
        B[m]=0;
    }
    mexPrintf("test text. %d, %d", nlhs, nrhs);
}