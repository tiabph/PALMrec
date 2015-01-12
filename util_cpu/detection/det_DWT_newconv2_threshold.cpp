/*new version of conv2 function, with threshold to reduce memory use
using sparce matrix
*/
#include "mex.h"
#include "matrix.h"
#include "time.h"
#include "math.h"
#include "string.h"

#define _T float

_T CalStd(_T *pdata, int datalen, _T mean);
_T CalMean(_T *pdata, int datalen);
void ThresholdImage(_T *pdata, _T *presult, int datalen, _T level);
void CalConv2_row(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing);
void CalConv2_col(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing);
void CalConv2(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int spacing);
void CalWavelet(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen, int wllevel);
void CalWaveletW23(_T *pimg, _T *presult_w2, _T *presult_w3, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen);
void CalWaveletW2(_T *pimg, _T *presult_w2, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen);
void CalWavelet_Threshold(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen,
							_T level, int type);

inline _T IUWT(_T A){
    return (((A+0.0848f)>0.0f)?(3.8247f*sqrtf(fabsf(A+0.0848f))):(-3.8247f*sqrtf(fabsf(A+0.0848f))));
}

// result = det_DTW_newconv2_threshold(img, wtkernel, threshLevel, retType)
//retType(string): 'W2', 'W3', 'W2+W3', 'W2*W3'
void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
	_T *pImg, *pResult, *pWT, *pbuf, *pbuf2;
	int m,n,imgsize[3] = {0,0,0};
	int imglen;
	const int *pImgSize;
	int imgcnt, imgoffset, rettype;
	_T level;
	char pstr[20];
	
	level = ((_T *)mxGetData(prhs[2]))[0];
	rettype = -1;
	mxGetString(prhs[3], pstr, 20);
	if(strcmp(pstr, "W2") ==0){
		rettype = 1;
	}
	if(strcmp(pstr, "W3") ==0){
		rettype = 2;
	}
	if(strcmp(pstr, "W2+W3") ==0){
		rettype = 4;
	}
	if(strcmp(pstr, "W2*W3") ==0){
		rettype = 8;
	}
	if(rettype <0){
		mexPrintf("unknow ret type, need 'W2', 'W3', W2+W3', 'W2*W3'\n");
		return;
	}
	
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
		//plhs[1] = mxCreateNumericArray(2,imgsize,mxSINGLE_CLASS,mxREAL);
	} else{//image stack
		pImgSize = mxGetDimensions(prhs[0]);
		m = pImgSize[0];
		n = pImgSize[1];	
		imgsize[0] = m;
		imgsize[1] = n;
		imglen = pImgSize[2];
		imgsize[2] = imglen;
		plhs[0] = mxCreateNumericArray(3,imgsize,mxSINGLE_CLASS,mxREAL);
		//plhs[1] = mxCreateNumericArray(3,imgsize,mxSINGLE_CLASS,mxREAL);
	}

	pResult = (_T *)mxGetData(plhs[0]);
	//pResult2 = (_T *)mxGetData(plhs[1]);
	pbuf = (_T *)malloc(m*n*7*sizeof(_T));
	
	for(imgcnt=0; imgcnt<imglen; imgcnt++){
		imgoffset = imgcnt*m*n;
		//CalWaveletW(pImg + imgoffset, pResult1 + imgoffset, pResult2 + imgoffset, pWT, pbuf, m, n, mxGetNumberOfElements(prhs[1]));
		CalWavelet_Threshold(pImg + imgoffset, pResult + imgoffset, pWT,pbuf, m, n, 
							mxGetNumberOfElements(prhs[1]),level, rettype);
	}

	free(pbuf);
}

//calculate multi level wavelet into result, need a buffer 2X+3X+2X=7X size of image
//with threshold and wavelet operation
void CalWavelet_Threshold(_T *pimg, _T *presult, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen,
							_T level, int type){
	int imgsize = img_w*img_h;
	int m;
	
	switch(type){
		case 1://w2 only
			CalWaveletW2(pimg, pbuf, pwt, pbuf+2*imgsize, img_w, img_h, wtlen);
			ThresholdImage(pbuf, presult, imgsize, level);
		break;
		case 2://w3 only
			CalWaveletW23(pimg, pbuf, pbuf+imgsize, pwt, pbuf+2*imgsize, img_w, img_h, wtlen);
			ThresholdImage(pbuf+imgsize, presult, imgsize, level);
		break;
		case 4://w2+w3
			CalWaveletW23(pimg, pbuf, pbuf+imgsize, pwt, pbuf+2*imgsize, img_w, img_h, wtlen);
			ThresholdImage(pbuf, pbuf, imgsize, level);
			ThresholdImage(pbuf+imgsize, pbuf+imgsize, imgsize, level);
			for(m=0; m<imgsize; m++){
				presult[m] = pbuf[m] + pbuf[m+imgsize];
			}
		break;
		case 8://w2*w3
			CalWaveletW23(pimg, pbuf, pbuf+imgsize, pwt, pbuf+2*imgsize, img_w, img_h, wtlen);
			ThresholdImage(pbuf, pbuf, imgsize, level);
			ThresholdImage(pbuf+imgsize, pbuf+imgsize, imgsize, level);
			for(m=0; m<imgsize; m++){
				presult[m] = pbuf[m] * pbuf[m+imgsize];
			}
		break;
		default:
			mexPrintf("internal error in det_DTW_newconv2_threshold!\n");
	}
}

//calculate multi level wavelet into result, need a buffer 3X+2X size of image
void CalWaveletW2(_T *pimg, _T *presult_w2, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen){
	int m,n;
	int imgsize = img_w * img_h;
	_T *pr1, *pr2, *pw;
	CalWavelet(pimg, pbuf, pwt, pbuf+3*imgsize, img_w, img_h, wtlen, 2);
	
	//tranform A1 A2 A3
	for(m=0; m<2*imgsize; m++){
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
	/*pr2 = pbuf + imgsize*2;
	pr1 = pbuf + imgsize;
	pw = presult_w3;
	for(m=0; m<imgsize; m++){
		pw[m] = pr1[m] - pr2[m];
	}*/
}

//calculate multi level wavelet into result, need a buffer 3X+2X size of image
void CalWaveletW23(_T *pimg, _T *presult_w2, _T *presult_w3, _T *pwt, _T *pbuf, int img_w, int img_h, int wtlen){
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

void ThresholdImage(_T *pdata, _T *presult, int datalen, _T level){
	_T mean, std, threshold;
	int m;
	
	mean = CalMean(pdata, datalen);
	std = CalStd(pdata, datalen, mean);
	threshold = mean + std * level;
	for(m=0;m<datalen;m++){
		presult[m] = pdata[m] > threshold? pdata[m] : threshold;
	}
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