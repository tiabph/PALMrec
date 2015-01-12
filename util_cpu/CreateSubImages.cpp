/*create sub-images from raw-image and detection info
*/
#include "mex.h"
#include "matrix.h"
#include "time.h"
#include "math.h"

#define _T float

//ROIResult = CreateSubImages(img, detectionbuf, fitl)
void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
	int img_height, img_width, imglen, detlen, resultSize[3];
	double *pDetBuf, cx, cy;
	_T *pResult, *pImg, temp, minInt;
	int fitl,rcx,rcy,rec_t,rec_b,rec_l,rec_r, frame;
	int m,n,mx,my,base;
	
	img_height = mxGetDimensions(prhs[0])[0];
	img_width = mxGetDimensions(prhs[0])[1];
	imglen = mxGetNumberOfDimensions(prhs[0])==2? 1: mxGetDimensions(prhs[0])[2];
	detlen = mxGetM(prhs[1]);
	fitl = mxGetPr(prhs[2])[0];
	pDetBuf = mxGetPr(prhs[1]);
	pImg = (_T *)mxGetData(prhs[0]);
	
	resultSize[0] = int(fitl*2+1);
	resultSize[1] = int(fitl*2+1);
	resultSize[2] = detlen;
	plhs[0] = mxCreateNumericArray(3, resultSize, mxSINGLE_CLASS, mxREAL);
	pResult = (_T *)mxGetData(plhs[0]);
	n=0;
	for(m=0; m<detlen; m++){
		frame = pDetBuf[m];
		cx = pDetBuf[m + detlen];
		cy = pDetBuf[m + detlen*2];
		rcx = int(cx+0.5);
		rcy = int(cy+0.5);
		rec_t = rcy - fitl -1;
		rec_b = rcy + fitl -1;
		rec_l = rcx - fitl -1;
		rec_r = rcx + fitl -1;
		
		if(rec_t >=0 && rec_l >=0 && rec_b<img_height && rec_r<img_width && frame<=imglen){
			//copy sub image
			base = (frame-1)*img_height*img_width;
			minInt = pImg[base + rec_t + rec_l*img_height];
			for(mx = rec_l; mx <=rec_r; mx++){
				for(my=rec_t; my<=rec_b; my++){
					temp = pImg[base + my + mx*img_height];
					pResult[n] = temp;
					minInt = minInt >temp? temp: minInt;
					n++;
				}
			}
			for(mx=1; mx<=resultSize[0]*resultSize[0]; mx++){
				pResult[n-mx] -= minInt;
			}
		} else{//fill with zeros
			for(mx=1; mx<=resultSize[0]*resultSize[0]; mx++){
				pResult[n] = 0.0f;
				n++;
			}
		}
	}
}