#include "mex.h"
#include "matrix.h"
/*result = det_DWT_Findparticles(img, samplewidth)
 *find local maxima and calculation weight center
 *img: H*W*FRAME matrix
 *samplewidth: image sample window width, must be odd
 *result: FRAME*1 cell contains N*# matrix for [X Y H] formate
 *
 *LS Gu, 2014.9.6 tiabph@gmail.com
 */
int FindParticles(float *pimg, int img_h, int img_w, 
        float *buf_h, float *buf_w, float *buf_peak, int samplewidth);
//result = det_DWT_Findparticles(img, samplewidth)
void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
    const int *pImgSize = mxGetDimensions(prhs[0]);
    int swidth = (int)(*((double *)mxGetData(prhs[1])));
    int Img_h = pImgSize[0], Img_w = pImgSize[1];
    int Img_len = (mxGetNumberOfDimensions(prhs[0])>2) ? pImgSize[2]:1;
    float *pImg = (float *)mxGetPr(prhs[0]);
    mxArray* presult = mxCreateCellMatrix(Img_len, 1);
    mxArray* pmat = NULL;
    double *pmatd = NULL;
    
    int imgcnt,pnum,m;
    float *pmask = (float *)malloc(Img_h * Img_w * 8);
    float *buf_h = (float *)malloc(Img_h * Img_w * 8);
    float *buf_w = (float *)malloc(Img_h * Img_w * 8);
    float *buf_peak = (float *)malloc(Img_h * Img_w * 8);
    
    for(imgcnt=0;imgcnt<Img_len;imgcnt++){
        pnum = FindParticles(pImg+imgcnt*Img_h * Img_w, Img_h, Img_w, 
        buf_h, buf_w, buf_peak, swidth/2);
        pmat = mxCreateDoubleMatrix(pnum,3,mxREAL);
        pmatd = (double*)mxGetData(pmat);
        for(m=0;m<pnum;m++){
            pmatd[m] = buf_w[m];
        }
        for(m=0;m<pnum;m++){
            pmatd[m+pnum] = buf_h[m];
        }
        for(m=0;m<pnum;m++){
            pmatd[m+pnum+pnum] = buf_peak[m];
        }
        mxSetCell(presult,imgcnt,pmat);
    }
    
    free(pmask);
    free(buf_h);
    free(buf_w);
    free(buf_peak);
    plhs[0] = presult;
}


int FindParticles(float *pimg, int img_h, int img_w, 
        float *buf_h, float *buf_w, float *buf_peak, int samplewidth){
    int pnum = 0;
    int ch,cw,mh,mw,flag;
    double temp,sumx,sumy,sum;
    int m,n;
    for(cw=samplewidth; cw<(img_w-samplewidth); cw++){
        for(ch=samplewidth; ch<(img_h-samplewidth); ch++){
            flag=1;
            for(mw=(cw-samplewidth); mw<=(cw+samplewidth); mw++){
                for(mh=(ch-samplewidth); mh<=(ch+samplewidth); mh++){
                    if(pimg[ch + cw*img_h] <= pimg[mh + mw*img_h] && (mh!=ch || mw!=cw)){
                        flag=0;
                        break;
                    }
                }
                if(flag==0){
                    break;
                }
            }
            if(flag!=0){
                sumx=0;
                sumy=0;
                sum=0;
                for(m=-samplewidth;m<=samplewidth;m++){
                    for(n=-samplewidth;n<=samplewidth;n++){
                        temp = pimg[(ch+n) + (cw+m)*img_h];
                        sum=sum+temp;
                        sumx=sumx+temp*m;
                        sumy=sumy+temp*n;
                    }
                }
                sumx=sumx/sum;
                sumy=sumy/sum;
                buf_h[pnum] = sumy+ch+1;
                buf_w[pnum] = sumx+cw+1;
                buf_peak[pnum] = temp;
                
                pnum++;
            }
                    
        }
    }
    return pnum;
}