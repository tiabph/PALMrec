#include "mex.h"
#include "math.h"
/*
 PostProcessingFittingData(databuf, param), change the databuf
 */

double CalStd(float *pimg, int height, int width, 
        int top, int bottom, int left, int right, int fnum);

void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
    mxArray *dataBuf = (mxArray *)(int)(prhs[0]);
    mxArray *fittingParam, *fitInfo_raw, *detectionBuf, *img, *fitInfo = NULL, *fitCnt = NULL;
    double fitl, factor, gain, pixelSize;
    int pointCnt, imgHeight, imgWidth, imgLen;
    // local pointer
    double *pResult = NULL, *pDetectionBuf, *pFitInfo, *pFitInfo_raw,*pFitInfo_new;
    float *pImg;
    //other var
    double cx,cy,rcx,rcy,peak,sx,sy,photonNumber,I1,I2,s,b,deltaR;
    int fnumber, top, bottom, left, right, first_frame;
    int m,n,tbase1, tbase2,recnt;
    
    fittingParam = mxGetField(prhs[1], 0, "fitting");
    fitInfo_raw = mxGetField(dataBuf, 0, "fitInfo_raw");
    detectionBuf = mxGetField(dataBuf, 0, "detectionBuf");
    img = mxGetField(dataBuf, 0, "img");
//     fitInfo = mxCreateDoubleMatrix(1,1,mxREAL);
            
    fitl = mxGetPr(mxGetField(fittingParam, 0, "fitl"))[0];
    factor = mxGetPr(mxGetField(fittingParam, 0, "factor"))[0];
    gain = mxGetPr(mxGetField(fittingParam, 0, "gain"))[0];
    pixelSize = mxGetPr(mxGetField(fittingParam, 0, "pixelsize"))[0];
    pointCnt = int(mxGetPr(mxGetField(dataBuf, 0, "pointCnt"))[0] + 0.5);
    imgHeight = int(mxGetPr(mxGetField(dataBuf, 0, "height"))[0] + 0.5);
    imgWidth = int(mxGetPr(mxGetField(dataBuf, 0, "width"))[0] + 0.5);
    imgLen = int(mxGetPr(mxGetField(dataBuf, 0, "imglen"))[0] + 0.5);
    
    pDetectionBuf = mxGetPr(detectionBuf);
    pFitInfo_raw = mxGetPr(fitInfo_raw);
    pImg = (float *)(mxGetData(img));
       
    pFitInfo = (double *)malloc(pointCnt * 9 * sizeof(double));
    recnt = 0;
    for(m=0; m<pointCnt; m++){
        cx = pDetectionBuf[m + pointCnt*1]; //cx=round(databuf.detectionBuf(i,2));
        cy = pDetectionBuf[m + pointCnt*2]; // cy=round(databuf.detectionBuf(i,3));
        rcx = int(cx+0.5);                  // rcx=round(cx);
        rcy = int(cy+0.5);                  // rcy=round(cy);
        fnumber = int(pDetectionBuf[m]);    // n=databuf.detectionBuf(i,1);
        top = int(rcy - fitl +0.5);          // up=round(cy)-fitl; 
        bottom = int(rcy + fitl +0.5);      // bottom=round(cy)+fitl; 
        left = int(rcx - fitl +0.5);        // left=round(cx)-fitl;
        right = int(rcx + fitl +0.5);       // right=round(cx)+fitl;
        
//         if up>=1 && left>=1 && bottom<=databuf.height  && right<=databuf.width  && n<databuf.imglen ...
//             && fitInfo(i,1)~= fitl+1 && fitInfo(i,2)~= fitl+1 ...
//             && fitInfo(i,1)~= 0 && fitInfo(i,2)~= 0 ...
//             && fitInfo(i,1)~= fitl*2+1 && fitInfo(i,2)~= fitl*2+1
        if(top>=1 && left>=1 && bottom<=imgHeight && right <= imgWidth 
                && pFitInfo_raw[m] != fitl+1.0 && pFitInfo_raw[m + pointCnt] != fitl+1.0
                && pFitInfo_raw[m] != 0.0 && pFitInfo_raw[m + pointCnt] != 0.0
                && pFitInfo_raw[m] != fitl*2+1.0 && pFitInfo_raw[m + pointCnt] != fitl*2+1.0){
//     % ---------- position ----------
//             cx=rcx+fitInfo(i,1)-fitl-1;
//             cy=rcy+fitInfo(i,2)-fitl-1; 
            cx = rcx + pFitInfo_raw[m] - fitl - 1.0;
            cy = rcy + pFitInfo_raw[m + pointCnt] - fitl - 1.0;
//     % ---------- peak ----------
//             peak=fitInfo(i,3);
            peak = pFitInfo_raw[m + pointCnt *2];
//     % ---------- sigma x and y ----------
//             sx=fitInfo(i,5);
//             sy=fitInfo(i,5);
            sx=sy=pFitInfo_raw[m + pointCnt *4];
//     % ---------- photon number ----------
//             Nm=fitInfo(i,3);
            photonNumber = pFitInfo_raw[m + pointCnt *2];
//     % ---------- frame number ----------
//             first_frame = n;
            first_frame = int(fnumber+0.5);
//             if(first_frame >1)
//                 I1=(databuf.img(rcy-fitl:rcy+fitl,rcx-fitl:rcx+fitl,first_frame-1)); 
//             else
//                 I1 = NaN;
//             end
            if(first_frame > 1){
                I1 = CalStd(pImg, imgHeight, imgWidth, 
                        top-1, bottom-1, left-1, right-1, first_frame-1 -1);
            } else{
                I1=-1.0;
            }
//             if(first_frame < databuf.pointCnt)
//                 I2=(databuf.img(rcy-fitl:rcy+fitl,rcx-fitl:rcx+fitl,first_frame+1));
//             else
//                 I2 = NaN;
//             end
            if(first_frame < imgLen){
                I2 = CalStd(pImg, imgHeight, imgWidth, 
                        top-1, bottom-1, left-1, right-1, first_frame+1 -1);
            } else{
                I2=I1;
            }
            
            if(I1<0){
                I1=I2;
            } else if(I1>I2){
                I1=I2;
            }
//              
//             s=sx*a; 
            s=sx*pixelSize;
            b = I1 *factor/gain;
            deltaR=sqrt( (s*s+pixelSize*pixelSize/12.0)/photonNumber 
                    + 8*3.14159*s*s*s*s/pixelSize/pixelSize/photonNumber/photonNumber*b*b);
//             
//             
//             fitcnt = fitcnt+1;
//             bkgBuf(fitcnt,1,:) = I1(:);
//             bkgBuf(fitcnt,2,:) = I2(:);
//             tempbuf(fitcnt,1) = (s^2+a^2/12)/Nm;
//             tempbuf(fitcnt,2) = 8*pi*s^4/a^2/Nm^2;
//             databuf.fitInfo(fitcnt,:)=[cx cy peak sx sy 0 Nm 0 first_frame];
//         end
            pFitInfo[recnt + pointCnt * 0] = cx;
            pFitInfo[recnt + pointCnt * 1] = cy;
            pFitInfo[recnt + pointCnt * 2] = peak;
            pFitInfo[recnt + pointCnt * 3] = sx;
            pFitInfo[recnt + pointCnt * 4] = sy;
            pFitInfo[recnt + pointCnt * 5] = deltaR;
            pFitInfo[recnt + pointCnt * 6] = photonNumber;
            pFitInfo[recnt + pointCnt * 7] = b;
            pFitInfo[recnt + pointCnt * 8] = first_frame;
            recnt++;
        }
    }
    
    //set output 
    fitInfo = mxCreateDoubleMatrix(recnt,9,mxREAL);
    fitCnt = mxCreateDoubleMatrix(1,1,mxREAL);
    pFitInfo_new = mxGetPr(fitInfo);
    for(n=0;n<9;n++){
        tbase1 = recnt*n;
        tbase2 = pointCnt*n;
        for(m=0;m<recnt; m++){
            pFitInfo_new[m + tbase1] = pFitInfo[m + tbase2];
        }
    }
//     mexPrintf("parameter.fitting.fitl = %f\n", fitl);
    plhs[0] = dataBuf;
    if(mxGetField(dataBuf, 0, "fitInfo")==NULL){
        mxAddField(dataBuf, "fitInfo");
    }
    mxSetField(dataBuf, 0, "fitInfo", fitInfo);
    
    //fitCnt
    mxGetPr(fitCnt)[0] = recnt;
    if(mxGetField(dataBuf, 0, "fitCnt")==NULL){
        mxAddField(dataBuf, "fitCnt");
    }
    mxSetField(dataBuf, 0, "fitCnt", fitCnt);
    
    //free memory
    free(pFitInfo);
}

//calculate std of the image
double CalStd(float *pimg, int height, int width, 
        int top, int bottom, int left, int right, int fnum){
    int base = height*width*(fnum);
    double mean=0.0, pixcnt=0.0, sum=0.0;
    int m,n,cx,cy;
    //mexPrintf("%d,%d,%d,%d,%d,%d,%d\n", height, width, top, bottom, left, right, fnum);
    for(cx = left; cx<=right; cx++){
        for(cy=top; cy<=bottom;cy++){
            sum += pimg[base + cy + cx * height];
            pixcnt +=1.0;
        }
    }
    mean = sum / pixcnt;
    sum=0;
    for(cx = left; cx<=right; cx++){
        for(cy=top; cy<=bottom;cy++){
            sum += pow(pimg[base + cy + cx * height]-mean, 2);
        }
    }
    sum = sum / (pixcnt-1.0);
    
    return sqrt(sum);
}