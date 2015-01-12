//cpu version of wavelet transform
// [W2 W3] = det_DWT_cpu(img, wavelet1, wavelet2, wavelet3)
#include "mex.h"
#include "matrix.h"

void CalConv2(double *pimg, double *pwt, int img_h, int img_w, int wtlen, double *presult);
void CalConv2_modified(double *pimg, double *pwt, int img_h, int img_w, int wtlen, double *presult);

void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
    const int *pImgSize = mxGetDimensions(prhs[0]);
    int Img_h = pImgSize[0], Img_w = pImgSize[1], Img_len = pImgSize[2];
    mxArray *pW2 = mxCreateNumericArray(3,pImgSize,mxDOUBLE_CLASS,mxREAL);
    mxArray *pW3 = mxCreateNumericArray(3,pImgSize,mxDOUBLE_CLASS,mxREAL);
    double *pImg = (double *)mxGetPr(prhs[0]);
    double *pW2d = (double *)mxGetData(pW2);
    double *pW3d = (double *)mxGetData(pW3);
    double *pWT1 = (double *)mxGetData(prhs[1]);
    double *pWT2 = (double *)mxGetData(prhs[2]);
    double *pWT3 = (double *)mxGetData(prhs[3]);
    double *pA1 = (double *)malloc(Img_h*Img_w*8);
    double *pA2 = (double *)malloc(Img_h*Img_w*8);
    int m,n;
    for(m=0;m<Img_len;m++){
        //copy data to A1
        for(n=0;n<Img_h*Img_w;n++){
            pA1[n] = pImg[n+m*Img_h*Img_w];
        }
        //cal W1
        CalConv2(pA1, pWT1, Img_h, Img_w, mxGetM(prhs[1]), pA2);
        //cal W2
        CalConv2(pA2, pWT2, Img_h, Img_w, mxGetM(prhs[2]), pA1);
        for(n=0;n<Img_h*Img_w;n++){
            pW2d[n+m*Img_h*Img_w] = pA2[n] - pA1[n];
        }
        //cal W3        
        CalConv2(pA1, pWT3, Img_h, Img_w, mxGetM(prhs[3]), pA2);
        for(n=0;n<Img_h*Img_w;n++){
            pW3d[n+m*Img_h*Img_w] = pA1[n] - pA2[n];
        }
    }
    free(pA1);
    free(pA2);
    plhs[0] = pW2;
    plhs[1] = pW3;
}


//calculation convolution of image
//consider padding and cutting the convolution result
void CalConv2(double *pimg, double *pwt, int img_h, int img_w, int wtlen, double *presult){
    int ch, cw, wh, ww, halflen = wtlen/2;
    int imgh,imgw,cy;
    double temp;
    for(cw=halflen;cw<img_w-halflen;cw++){
        for(ch=halflen; ch<img_h-halflen; ch++){
            temp=0.0;
            for(ww=0;ww<wtlen;ww++){
                for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    //if(imgh<0){
                    //    imgh = -imgh;
                    //}
                    //if(imgh>=img_h){
                    //    imgh = img_h+img_h-imgh-2;
                    //}
                    
                    //if(imgw<0){
                    //    imgw = -imgw;
                    //}
                    //if(imgw>=img_w){
                    //    imgw = img_w+img_w-imgw-2;
                    //}
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                }
            }
            presult[ch + cw*img_h] = temp;
        }
    }
}

void CalConv2_modified(double *pimg, double *pwt, int img_h, int img_w, int wtlen, double *presult){
    int ch, cw, wh, ww, halflen = wtlen/2, wlen = wtlen/2;
    int imgh,imgw,cy;
    double temp;
    
    for(cw=0;cw<wlen;cw++){
         for(ch=0;ch<wlen;ch++){//region 0
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    if(imgh<0){
                        imgh = -imgh;
                    }
                    if(imgw<0){
                        imgw = -imgw;
                    }
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }  
                }
                presult[ch + cw*img_h] = temp;
         }
         for(ch=wlen;ch<(img_h-wtlen);ch++){//region 1
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    if(imgw<0){
                        imgw = -imgw;
                    }
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }    
                }            
            presult[ch + cw*img_h] = temp;
         }
         for(ch=(img_h-wtlen);ch<img_h;ch++){//region 2
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    if(imgh>=img_h){
                        imgh = img_h+img_h-imgh-2;
                    }
                    if(imgw<0){
                        imgw = -imgw;
                    }
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }     
            }
            presult[ch + cw*img_h] = temp;
         }
    }
    for(cw=wlen;cw<(img_w-wtlen);cw++){
         for(ch=0;ch<wlen;ch++){//region 3
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    if(imgh<0){
                        imgh = -imgh;
                    }                    
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }          
            }
            presult[ch + cw*img_h] = temp;
         }
         for(ch=wlen;ch<(img_h-wtlen);ch++){//region 4
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }      
            }
            presult[ch + cw*img_h] = temp;
         }
         for(ch=(img_h-wtlen);ch<img_h;ch++){//region 5
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    if(imgh>=img_h){
                        imgh = img_h+img_h-imgh-2;
                    }
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }     
            }
            presult[ch + cw*img_h] = temp;
         }
    }
    for(cw=(img_w-wtlen);cw<img_w;cw++){
         for(ch=0;ch<wlen;ch++){//region 6
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    if(imgh<0){
                        imgh = -imgh;
                    }                 
                    if(imgw>=img_w){
                        imgw = img_w+img_w-imgw-2;
                    }   
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }   
            }
            presult[ch + cw*img_h] = temp;
         }
         for(ch=wlen;ch<(img_h-wtlen);ch++){//region 7
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    if(imgw>=img_w){
                        imgw = img_w+img_w-imgw-2;
                    }
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }    
            }
            presult[ch + cw*img_h] = temp;
         }
         for(ch=(img_h-wtlen);ch<img_h;ch++){//region 8
             temp=0.0;
             for(ww=0;ww<wtlen;ww++){
                 for(wh=0;wh<wtlen;wh++){
                    imgh = ch + wh - halflen;
                    imgw = cw + ww - halflen;
                    if(imgh>=img_h){
                        imgh = img_h+img_h-imgh-2;
                    }
                    if(imgw>=img_w){
                        imgw = img_w+img_w-imgw-2;
                    }
                    temp=temp + pimg[imgh + imgw*img_h]*pwt[wh + ww*wtlen];
                    }
            }
            presult[ch + cw*img_h] = temp;
         }
    }
}