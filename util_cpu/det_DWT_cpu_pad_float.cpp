//cpu version of wavelet transform
// [W2 W3] = det_DWT_cpu_pad(img, wavelet1, wavelet2, wavelet3)
#include "mex.h"
#include "matrix.h"
#include "time.h"
#include "math.h"

#define _T float

_T IUWT(_T A){
    return (((A+0.0848f)>0.0f)?(3.8247f+sqrtf(fabsf(A+0.0848f))):0.0f);
}

void CalConv2(_T *pimg, _T *pwt, int img_h, int img_w, int imglen, int wtlen, _T *presult);
void ImgExpand(_T *pimg, int h, int w, int len, int pw, _T *presult);
void ImgPad(_T *pimg, int h, int w, int len, int pw);
void ImgDexpand(_T *pimg, int h, int w, int len, int pw, _T *presult);
void CalConv2_cr(_T *pimg, _T *pwt, int img_h, int img_w, int imglen, int wtlen, _T *presult);
void CalConv2_cr_nopad(_T *pimg, _T *pwt, int img_h, int img_w, int imglen, int wtlen, _T *presult);
void CalConv2_cr_nopad_frame(_T *pimg, _T *pwt, int img_h, int img_w, int wtlen, _T *presult, _T *pbuf);

void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
    const int *pImgSize = mxGetDimensions(prhs[0]);
    int Img_h = pImgSize[0], Img_w = pImgSize[1];
    int Img_len = (mxGetNumberOfDimensions(prhs[0])>2) ? pImgSize[2]:1;
    int plen1 = mxGetNumberOfElements(prhs[1]),plen2 = mxGetNumberOfElements(prhs[2]),plen3 = mxGetNumberOfElements(prhs[3]);
    mxArray *pW2;// = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),pImgSize,mxDOUBLE_CLASS,mxREAL);
    mxArray *pW3;// = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),pImgSize,mxDOUBLE_CLASS,mxREAL);
    _T *pImg;// = (_T *)mxGetPr(prhs[0]);
    _T *pW2d;// = (_T *)mxGetData(pW2);
    _T *pW3d;// = (_T *)mxGetData(pW3);
    _T *pWT1;// = (_T *)mxGetData(prhs[1]);
    _T *pWT2;// = (_T *)mxGetData(prhs[2]);
    _T *pWT3;// = (_T *)mxGetData(prhs[3]);
    _T *pA1;// = (_T *)malloc(Img_h*Img_w*8);
    _T *pA2;// = (_T *)malloc(Img_h*Img_w*8);
    _T *pBuf;// = (_T *)malloc(Img_h*Img_w*8);
    int m,n,i,imgcnt,offset;
    int clk_s,clk_e,time_expand=0,time_conv2=0;
    
    pW2 = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),pImgSize,mxSINGLE_CLASS,mxREAL);
    pW3 = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),pImgSize,mxSINGLE_CLASS,mxREAL);
    pImg = (_T *)mxGetPr(prhs[0]);
    pW2d = (_T *)mxGetData(pW2);
    pW3d = (_T *)mxGetData(pW3);
    pWT1 = (_T *)mxGetData(prhs[1]);
    pWT2 = (_T *)mxGetData(prhs[2]);
    pWT3 = (_T *)mxGetData(prhs[3]);
    pA1 = (_T *)malloc(Img_h*Img_w*sizeof(_T));
    pA2 = (_T *)malloc(Img_h*Img_w*sizeof(_T));
    pBuf = (_T *)malloc(Img_h*Img_w*sizeof(_T));
	for(imgcnt=0;imgcnt<Img_len;imgcnt++){
		offset = imgcnt*Img_h*Img_w;
		//cal W1->A1
		CalConv2_cr_nopad_frame(pImg+offset, pWT1, Img_h, Img_w, plen1, pA1, pBuf);
		
		//cal W2->A2
		CalConv2_cr_nopad_frame(pA1, pWT2, Img_h, Img_w, plen2, pA2, pBuf);
		for(n=0;n<Img_h*Img_w;n++){
			   pW2d[offset + n] = IUWT(pA1[n]) - IUWT(pA2[n]);
	             //pW2d[n] = pA1[n];
	             //pW3d[n] = pA2[n];
		}
		
		//cal W3->A1
		CalConv2_cr_nopad_frame(pA2, pWT3, Img_h, Img_w, plen3, pA1, pBuf);
		for(n=0;n<Img_h*Img_w;n++){
			   pW3d[offset + n] = IUWT(pA2[n]) - IUWT(pA1[n]);
		}
	}
    free(pA1);
    free(pA2);
    free(pBuf);
    plhs[0] = pW2;
    plhs[1] = pW3;
}


//calculation convolution of image
//consider padding and cutting the convolution result
void CalConv2(_T *pimg, _T *pwt, int img_h, int img_w, int imglen, int wtlen, _T *presult){
    int ch, cw, wh, ww, halflen = wtlen/2,offset;
    int imgh,imgw,cy;
    _T temp;
    int imgcnt;
    for(imgcnt=0;imgcnt<imglen;imgcnt++){
        offset = imgcnt * img_h * img_w;
        for(cw=halflen;cw<img_w-halflen;cw++){
            for(ch=halflen; ch<img_h-halflen; ch++){
                temp=0.0;
                for(ww=0;ww<wtlen;ww++){
                    for(wh=0;wh<wtlen;wh++){
                        imgh = ch + wh - halflen;
                        imgw = cw + ww - halflen;
                        temp=temp + pimg[offset + imgh + imgw*img_h] * pwt[(wtlen - wh -1) + (wtlen - ww -1)*wtlen];
                    }
                }
                presult[offset + ch + cw*img_h] = temp;
            }
        }
    }
}

void CalConv2_cr(_T *pimg, _T *pwt, int img_h, int img_w, int imglen, int wtlen, _T *presult){
    int ch, cw, wh, ww, halflen = wtlen/2,offset;
    int imgh,imgw,cy;
    _T temp;
    int imgcnt;
    for(imgcnt=0;imgcnt<imglen;imgcnt++){
        offset = imgcnt * img_h * img_w;
        for(cw=0;cw<img_w;cw++){
            for(ch=halflen; ch<img_h-halflen; ch++){
                temp=0.0;
                for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw;
                    temp=temp + pimg[offset + imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                presult[offset + ch + cw*img_h] = temp;
            }
        }
        for(cw=halflen;cw<(img_w-halflen);cw++){
            for(ch=0; ch<img_h; ch++){
                temp=0.0;
                for(ww=0;ww<wtlen;ww++){
                    imgh = ch;
                    imgw = cw + ww - halflen;
                    temp=temp + presult[offset + imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                pimg[offset + ch + cw*img_h] = temp;
            }
        }
    }
}

void CalConv2_cr_nopad(_T *pimg, _T *pwt, int img_h, int img_w, int imglen, int wtlen, _T *presult){
    int ch, cw, wh, ww, halflen = wtlen/2,offset;
    int imgh,imgw,cy;
    _T temp;
    int imgcnt;
    _T *pbuf = (_T *)malloc(img_h * img_w * sizeof(_T));
    for(imgcnt=0;imgcnt<imglen;imgcnt++){
        offset = imgcnt * img_h * img_w;
        //h-direction
        for(cw=0;cw<img_w;cw++){
            for(ch=0; ch< halflen; ch++){
                temp=0.0;
                for(wh=0;wh<halflen-ch;wh++){
                    imgh = halflen - ch - wh - 1;
                    imgw = cw;
                    temp=temp + pimg[offset + imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                for(wh=halflen-ch;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw;
                    temp=temp + pimg[offset + imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                pbuf[ch + cw*img_h] = temp;
            }
            for(ch=halflen; ch<img_h-halflen; ch++){
                temp=0.0;
                for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw;
                    temp=temp + pimg[offset + imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                pbuf[ch + cw*img_h] = temp;
            }
            for(ch=img_h-halflen; ch< img_h; ch++){
                temp=0.0;
                for(wh=0;wh<halflen+img_h-ch;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw;
                    temp=temp + pimg[offset + imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                for(wh=halflen+img_h-ch;wh<wtlen;wh++){
                    imgh = img_h + img_h - ch - wh + halflen -1;
                    imgw = cw;
                    temp=temp + pimg[offset + imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                pbuf[ch + cw*img_h] = temp;
            }
        }
        //w-direction
        for(ch=0;ch<img_h;ch++){
            for(cw=0; cw< halflen; cw++){
                temp=0.0;
                for(ww=0;ww<halflen-cw;ww++){
                    imgw = halflen - cw - ww - 1;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                for(ww=halflen-cw;ww<wtlen;ww++){
                    imgw = cw + ww - halflen;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                presult[offset + ch + cw*img_h] = temp;
            }
            
            for(cw=halflen; cw<img_w-halflen; cw++){
                temp=0.0;
                for(ww=0;ww<wtlen;ww++){
                    imgw = cw + ww - halflen;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                presult[offset + ch + cw*img_h] = temp;
            }
            
            for(cw=img_w-halflen; cw< img_w; cw++){
                temp=0.0;
                for(ww=0;ww<halflen+img_w-cw;ww++){
                    imgw = cw + ww - halflen;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                for(ww=halflen+img_w-cw;ww<wtlen;ww++){
                    imgw = img_w + img_w - cw - ww + halflen -1;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                presult[offset + ch + cw*img_h] = temp;
            }
        }
    }
    free(pbuf);
}

void CalConv2_cr_nopad_frame(_T *pimg, _T *pwt, int img_h, int img_w, int wtlen, _T *presult, _T *pbuf){
	int ch, cw, wh, ww, halflen = wtlen/2;
    int imgh,imgw,cy;
    _T temp;
    //int imgcnt;
    //_T *pbuf = (_T *)malloc(img_h * img_w * sizeof(_T));
    //for(imgcnt=0;imgcnt<imglen;imgcnt++){
        //offset = imgcnt * img_h * img_w;
        //h-direction
        for(cw=0;cw<img_w;cw++){
            for(ch=0; ch< halflen; ch++){
                temp=0.0;
                for(wh=0;wh<halflen-ch;wh++){
                    imgh = halflen - ch - wh - 1;
                    imgw = cw;
                    temp=temp + pimg[imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                for(wh=halflen-ch;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw;
                    temp=temp + pimg[imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                pbuf[ch + cw*img_h] = temp;
            }
            for(ch=halflen; ch<img_h-halflen; ch++){
                temp=0.0;
                for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw;
                    temp=temp + pimg[imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                pbuf[ch + cw*img_h] = temp;
            }
            for(ch=img_h-halflen; ch< img_h; ch++){
                temp=0.0;
                for(wh=0;wh<halflen+img_h-ch;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw;
                    temp=temp + pimg[imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                for(wh=halflen+img_h-ch;wh<wtlen;wh++){
                    imgh = img_h + img_h - ch - wh + halflen -1;
                    imgw = cw;
                    temp=temp + pimg[imgh + imgw*img_h] * pwt[(wtlen - wh -1)];
                }
                pbuf[ch + cw*img_h] = temp;
            }
        }
        //w-direction
        for(ch=0;ch<img_h;ch++){
            for(cw=0; cw< halflen; cw++){
                temp=0.0;
                for(ww=0;ww<halflen-cw;ww++){
                    imgw = halflen - cw - ww - 1;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                for(ww=halflen-cw;ww<wtlen;ww++){
                    imgw = cw + ww - halflen;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                presult[ch + cw*img_h] = temp;
            }
            
            for(cw=halflen; cw<img_w-halflen; cw++){
                temp=0.0;
                for(ww=0;ww<wtlen;ww++){
                    imgw = cw + ww - halflen;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                presult[ch + cw*img_h] = temp;
            }
            
            for(cw=img_w-halflen; cw< img_w; cw++){
                temp=0.0;
                for(ww=0;ww<halflen+img_w-cw;ww++){
                    imgw = cw + ww - halflen;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                for(ww=halflen+img_w-cw;ww<wtlen;ww++){
                    imgw = img_w + img_w - cw - ww + halflen -1;
                    imgh = ch;
                    temp=temp + pbuf[imgh + imgw*img_h] * pwt[(wtlen - ww -1)];
                }
                presult[ch + cw*img_h] = temp;
            }
        }
    //}
    //free(pbuf);
}

//padding functions
void ImgExpand(_T *pimg, int h, int w, int len, int pw, _T *presult){
	int m,n,imgcnt,offset,offset2;
    for(imgcnt=0;imgcnt<len;imgcnt++){
        offset = imgcnt*(h+pw+pw)*(w+pw+pw);
        offset2 = imgcnt*(h)*(w);
        for(m=0;m<w;m++){
            for(n=0;n<h;n++){
                presult[offset + n+pw + (m+pw)*(h+pw+pw)] = pimg[offset2 + n+m*h];
            }
        }
    }
}

void ImgPad(_T *pimg, int h, int w, int len, int pw){
	int ch,cw,chr,cwr;//current position and raw position
	int m,n,imgcnt,offset;
	ch = 0;
	cw = 0;
    for(imgcnt=0;imgcnt<len;imgcnt++){
        offset = imgcnt*(h+pw+pw)*(w+pw+pw);
        //corner
        for(m=0;m<pw;m++){
            for(n=0;n<pw;n++){
                pimg[offset + m*(h+pw+pw)+n] = pimg[offset + (pw+pw-m-1)*(h+pw+pw)+(pw+pw-n-1)];//UL:pimg[m][n] = pimg[pw*2-m-1][pw*2-n-1]
                pimg[offset + m*(h+pw+pw)+(h+pw+pw-n-1)] = pimg[offset + (pw+pw-m-1)*(h+pw+pw)+(h+n)];//DL:pimg[m][n] = pimg[pw*2-m-1][pw*2-n-1]
                pimg[offset + (w+pw+pw-m-1)*(h+pw+pw)+n] = pimg[offset + (w+m)*(h+pw+pw)+(pw+pw-n-1)];//UR:pimg[m][n] = pimg[pw*2-m-1][pw*2-n-1]
                pimg[offset + (w+pw+pw-m-1)*(h+pw+pw)+(h+pw+pw-n-1)] = pimg[offset + (w+m)*(h+pw+pw)+(h+n)];//DR:pimg[m][n] = pimg[pw*2-m-1][pw*2-n-1]
            }
        }
        //L R side
        for(m=0;m<pw;m++){
            for(n=pw;n<(h+pw);n++){
                pimg[offset + m*(h+pw+pw)+n] = pimg[offset + (pw+pw-m-1)*(h+pw+pw)+n];//L
                pimg[offset + (w+pw+pw-m-1)*(h+pw+pw)+n] = pimg[offset + (w+m)*(h+pw+pw)+n];//R
            }
        }
        //U D side
        for(m=pw;m<(w+pw);m++){
            for(n=0;n<pw;n++){
                pimg[offset + m*(h+pw+pw)+n] = pimg[offset + (m)*(h+pw+pw)+(pw+pw-n-1)];//U
                pimg[offset + m*(h+pw+pw)+(h+pw+pw-n-1)] = pimg[offset + (m)*(h+pw+pw)+(h+n)];//D
            }
        }
    }
}

void ImgDexpand(_T *pimg, int h, int w, int len, int pw, _T *presult){
	int m,n,imgcnt,offset,offset2;
    for(imgcnt=0;imgcnt<len;imgcnt++){
        offset = imgcnt*(h+pw+pw)*(w+pw+pw);
        offset2 = imgcnt*(h)*(w);
        for(m=0;m<w;m++){
            for(n=0;n<h;n++){
                presult[offset2 + n+m*h]=pimg[offset + n+pw + (m+pw)*(h+pw+pw)];
            }
        }
    }
}
