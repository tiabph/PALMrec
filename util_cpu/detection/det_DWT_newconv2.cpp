/*new version of conv2 function
using sparce matrix
*/
#include "mex.h"
#include "matrix.h"
#include "time.h"
#include "math.h"

#define _T float

void CalConv2_row(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing);
void CalConv2_col(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing);
void CalConv2(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing);
void CalWavelet(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int wllevel);
void CalWaveletW(_T *pimg, _T *presult_w2, _T *presult_w3, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen);
inline _T IUWT(_T A){
    return (((A+0.0848f)>0.0f)?(3.8247f*sqrtf(fabsf(A+0.0848f))):(-3.8247f*sqrtf(fabsf(A+0.0848f))));
}

void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
	_T *pImg, *pResult1, *pResult2, *pWT, *pbuf, *pbuf2;
	int m,n,imgsize[3] = {0,0,0};
	int imglen;
	const int *pImgSize;
	int imgcnt, imgoffset;
	pImg = (_T *)mxGetData(prhs[0]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);	
	imgsize[0] = m;
	imgsize[1] = n;
	pWT = (_T *)mxGetData(prhs[1]);
	if(mxGetNumberOfDimensions(prhs[0])==2){//single image
		imglen = 1;
		imgsize[2]=1;
		plhs[0] = mxCreateNumericArray(2,imgsize,mxSINGLE_CLASS,mxREAL);
		plhs[1] = mxCreateNumericArray(2,imgsize,mxSINGLE_CLASS,mxREAL);
	} else{//image stack
		pImgSize = mxGetDimensions(prhs[0]);
		m = pImgSize[0];
		n = pImgSize[1];	
		imgsize[0] = m;
		imgsize[1] = n;
		imglen = pImgSize[2];
		imgsize[2] = imglen;
		plhs[0] = mxCreateNumericArray(3,imgsize,mxSINGLE_CLASS,mxREAL);
		plhs[1] = mxCreateNumericArray(3,imgsize,mxSINGLE_CLASS,mxREAL);
	}

	pResult1 = (_T *)mxGetData(plhs[0]);
	pResult2 = (_T *)mxGetData(plhs[1]);
	pbuf = (_T *)malloc(m*n*5*sizeof(_T));
	
	for(imgcnt=0; imgcnt<imglen; imgcnt++){
		imgoffset = imgcnt*m*n;
		CalWaveletW(pImg + imgoffset, pResult1 + imgoffset, pResult2 + imgoffset, 
					pWT, pbuf, m, n, mxGetNumberOfElements(prhs[1]));
	}

	free(pbuf);
}

//calculate multi level wavelet into result, need a buffer 3X+2X size of image
void CalWaveletW(_T *pimg, _T *presult_w2, _T *presult_w3, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen){
	int m,n;
	int imgsize = img_w * img_h;
	_T *pr1, *pr2, *pw;
	CalWavelet(pimg, pbuf, pwt, pbuf+3*imgsize, img_w, img_h, wtlen, 3);
	
	//tranform A1 A2 A3
	for(m=0; m<3*imgsize; m++){
		pbuf[m] = IUWT(pbuf[m]);
	}
	//calculate W2 = A1-A2
	pr2 = pbuf + imgsize;
	pr1 = pbuf;
	pw = presult_w2;
	for(m=0; m<imgsize; m++){
		pw[m] = pr1[m] - pr2[m];
	}
	//calculate W3 = A2-A3
	pr2 = pbuf + imgsize*2;
	pr1 = pbuf + imgsize;
	pw = presult_w3;
	for(m=0; m<imgsize; m++){
		pw[m] = pr1[m] - pr2[m];
	}
}

//calculate multi level wavelet into result, need a buffer 2X size of image
void CalWavelet(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int wllevel){
	int m,n,levelcnt;
	_T *base_r, *base_w;
	for(levelcnt=1; levelcnt<=wllevel; levelcnt++){
		if(levelcnt==1){
			base_r = pimg;
		} else{
			base_r = presult + (levelcnt-2)*img_w*img_h;
		}
		base_w = presult + (levelcnt-1)*img_w*img_h;
		CalConv2(base_r, base_w, pwt, pbuf, img_w, img_h, wtlen, 1<<(levelcnt-1));
	}
}

//calculate the conv2 result, need a buffer 2X size of image
void CalConv2(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing){
	CalConv2_row(pimg, pbuf + img_w*img_h, pwt, pbuf, img_w, img_h, wtlen, spacing);
	CalConv2_col(pbuf + img_w*img_h, presult, pwt, pbuf, img_w, img_h, wtlen, spacing);
}

//row-direction conv2, need a buffer size of one line
void CalConv2_row(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing){
	_T crntWT;
	int crntOffset, crntBase;
	int hcnt,wcnt,m,n;
	
	for(hcnt=0; hcnt<img_h; hcnt++){
		crntBase = hcnt*img_w;
		//pos 0:
		crntWT = pwt[0];
		crntOffset = 0;
		for(m=0; m<img_w; m++){
			presult[m + crntBase] = pimg[m +crntBase]*crntWT;
		}
		
		//pose 1 to end:
		for(wcnt=1; wcnt<wtlen; wcnt++){
			crntWT = pwt[wcnt];
			crntOffset = wcnt * spacing;
			for(m=0; m<img_w; m++){
				pbuf[m] = pimg[m + crntBase]*crntWT;
			}
			//positive direction:
			for(m=0; m<(img_w - crntOffset); m++){
				presult[m + crntBase] += pbuf[m + crntOffset];
			}
			for(m=(img_w - crntOffset); m<img_w; m++){
				presult[m + crntBase] += pbuf[img_w*2 - (m + crntOffset) -1];
			}
			//negative direction:
			for(m=0; m<crntOffset; m++){
				presult[m + crntBase] += pbuf[crntOffset - m -1];
			}
			for(m=crntOffset; m<img_w; m++){
				presult[m + crntBase] += pbuf[m - crntOffset];
			}
		}
	}
}

//column-direction conv2, need a buffer size of one frame
void CalConv2_col(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing){
	_T crntWT;
	int crntOffset, crntBase1, crntBase2;
	int hcnt,wcnt,m,n,imgsize = img_w*img_h;
	
	//pos 0:
	crntWT = pwt[0];
	crntOffset = 0;
	for(m=0; m<imgsize; m++){
		presult[m] = pimg[m]*crntWT;
	}
	
	//pose 1 to end:
	for(wcnt=1; wcnt<wtlen; wcnt++){
		crntWT = pwt[wcnt];
		crntOffset = wcnt * spacing;
		for(m=0; m<imgsize; m++){
			pbuf[m] = pimg[m]*crntWT;
		}
		//positive direction:
		
		for(m=0; m<(img_h - crntOffset); m++){
			//presult[m] += pbuf[m + crntOffset];
			crntBase1 = m*img_w;
			crntBase2 = (m + crntOffset)*img_w;
			for(n=0;n<img_w;n++){
				presult[crntBase1 + n] += pbuf[crntBase2 + n];
			}
		}
		for(m=(img_h - crntOffset); m<img_h; m++){
			//presult[m] += pbuf[img_h - (m + crntOffset) -1];
			crntBase1 = m*img_w;
			crntBase2 = (img_h*2 - (m + crntOffset) -1)*img_w;
			for(n=0;n<img_w;n++){
				presult[crntBase1 + n] += pbuf[crntBase2 + n];
			}
		}
		
		//negative direction:
		for(m=0; m<crntOffset; m++){
			//presult[m] += pbuf[crntOffset - m -1];
			crntBase1 = m*img_w;
			crntBase2 = (crntOffset - m -1)*img_w;
			for(n=0;n<img_w;n++){
				presult[crntBase1 + n] += pbuf[crntBase2 + n];
			}
		}
		for(m=crntOffset; m<img_h; m++){
			//presult[m] += pbuf[m - crntOffset];
			crntBase1 = m*img_w;
			crntBase2 = (m - crntOffset)*img_w;
			for(n=0;n<img_w;n++){
				presult[crntBase1 + n] += pbuf[crntBase2 + n];
			}
		}
	}
}