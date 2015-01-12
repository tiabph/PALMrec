/*threshold the wavelet result
*/
#include "mex.h"
#include "matrix.h"
#include "time.h"
#include "math.h"

#define _T float

_T CalMean(_T *pdata, int datalen);
_T CalStd(_T *pdata, int datalen, _T mean);

//result = det_DWT_thresh(img, level), threshold = mean + level*std
void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
	mxArray *Result;
	_T *pImg, *pResult;
	const int *psize;
	int imgsize, imglen, sizenum,m,n,base;
	_T mean,std,level,threshold;
	
	sizenum = mxGetNumberOfDimensions(prhs[0]);
	psize = mxGetDimensions(prhs[0]);
	imgsize = psize[0] * psize[1];
	Result = mxCreateNumericArray(sizenum,psize,mxSINGLE_CLASS,mxREAL);
	pImg = (_T *)mxGetData(prhs[0]);
	pResult = (_T *)mxGetData(Result);
	level = (_T)(mxGetPr(prhs[1])[0]);
	if(sizenum ==2){
		imglen = 1;
	} else{
		imglen = psize[2];
	}
	for(m=0; m<imglen; m++){
		base = m*imgsize;
		mean = CalMean(pImg + base, imgsize);
		std = CalStd(pImg + base, imgsize, mean);
		threshold = mean + std * level;
		for(n=0;n<imgsize;n++){
			pResult[base + n] = pImg[base + n] > threshold? pImg[base + n] : threshold;
		}
	}
	
	plhs[0] = Result;
}

_T CalMean(_T *pdata, int datalen){
	_T sum = 0.0f;
	int m;
	for(m=0; m<datalen; m++){
		sum += pdata[m];
	}
	sum = sum / datalen;
	return sum;
}

_T CalStd(_T *pdata, int datalen, _T mean){
	_T sum = 0.0f, temp;
	int m;
	for(m=0; m<datalen; m++){
		temp = pdata[m] - mean;
		sum += temp*temp;
	}
	sum = sum / (datalen-1);
	return sqrtf(sum);
}