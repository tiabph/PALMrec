#include "mex.h"
#include "stdlib.h"
#include "math.h"

typedef long int32;
//result = LinkPoints(fitdata, linkdata, startpoints)
//format:[cx cy peak sx sy deltaR Nm b first-frame idx_in_fitdata]
void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
	double *pFitData;
	int32 *pLinkData, *pSPData;
	size_t fitnum, renum, fitdataSize;
	mxArray *Result;
	double *pResult;
	double tbuf[9];
	int32 m,n,p;
	
	pFitData = mxGetPr(prhs[0]);
	pLinkData = (int32 *)mxGetData(prhs[1]);
	pSPData = (int32 *)mxGetData(prhs[2]);
	
	fitnum = mxGetM(prhs[0]);
	fitdataSize = mxGetN(prhs[0]);
	renum = mxGetM(prhs[2]);
	
	if(fitdataSize!=9){
		mexPrintf("This function only accepts fitdata with size of X*9\n");
		return;
	}
	
	Result = mxCreateNumericMatrix(renum, fitdataSize+1, mxDOUBLE_CLASS, mxREAL);
	pResult = mxGetPr(Result);
	
	for(m=0;m<renum;m++){
		for(n=0;n<fitdataSize; n++){
			tbuf[n]=0.0;
		}
		p = pSPData[m]-1;
		while(p>=0){
			for(n=0;n<fitdataSize; n++){
				tbuf[n] += pFitData[p+n*fitnum];
			}
			p=pLinkData[p]-1;
		}
		for(n=0;n<fitdataSize; n++){
			tbuf[n] = tbuf[n] / pSPData[m + renum];//mean of raw data
		}
		tbuf[6] *= pSPData[m + renum];//sum of photon number
		tbuf[8] = pFitData[pSPData[m]-1 + 8*fitnum]; //first frame
		
		//copy data
		for(n=0;n<fitdataSize; n++){
			pResult[m + n*renum] = tbuf[n];
		}
		pResult[m + fitdataSize*renum] = pSPData[m];
	}
	
	plhs[0] = Result;
}