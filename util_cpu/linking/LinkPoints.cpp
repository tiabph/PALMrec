#include "mex.h"
#include "stdlib.h"
#include "math.h"

typedef long int32;
//[result startpoints] = LinkPoints(fitData, [xidx yidx frameidx], [imgheight imgwidth imglen], [gap, dist])
void LinkPoints(double *px, double *py, double *pframe, 
				size_t pnum, double gap, double dist, 
				size_t imgWidth, size_t img_Height,
				int32 *presult);
int ListStartPoints(int32 *pLinkData, size_t pnum, int32 *presult);

void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
    double *pFitData, *pXdata, *pYdata, *pFrame;
    size_t xidx, yidx, frameidx, imgWidth, imgHeight, imgLen, pointCnt;
    double gap, dist;
    mxArray *Result, *Result_sp;
    int32 *pResult, *pspbuf, *pResult_sp;
    int spcnt,m,n;
    //check inputs
    if(nrhs<4){
		mexPrintf("need 3 parameters\nusage: result = LinkPoints(fitData, [xidx yidx frameidx], [imgheight imgwidth imglen], [gap, dist])\n");
		plhs[0] = NULL;
		return;
    }
    //get input data
    pFitData = mxGetPr(prhs[0]);
    pointCnt = mxGetM(prhs[0]);
    xidx = size_t(mxGetPr(prhs[1])[0]);
    yidx = size_t(mxGetPr(prhs[1])[1]);
    frameidx = size_t(mxGetPr(prhs[1])[2]);
    imgHeight = size_t(mxGetPr(prhs[2])[0]);
    imgWidth = size_t(mxGetPr(prhs[2])[1]);
    imgLen = size_t(mxGetPr(prhs[2])[2]);
    gap = mxGetPr(prhs[3])[0];
    dist = mxGetPr(prhs[3])[1];
    pXdata = pFitData + pointCnt*(xidx-1);
    pYdata = pFitData + pointCnt*(yidx-1);
    pFrame = pFitData + pointCnt*(frameidx-1);
    
    Result = mxCreateNumericMatrix(pointCnt, 1, mxINT32_CLASS, mxREAL);
    pResult = (int32 *)mxGetData(Result);
    
    LinkPoints(pXdata, pYdata, pFrame, 
				pointCnt, gap, dist, 
				imgWidth, imgHeight,
				pResult);
	plhs[0] = Result;
	//calculate start points
	if(nlhs>=2){
		pspbuf = (int32 *)mxMalloc(pointCnt*sizeof(int32)*2);
		spcnt = ListStartPoints(pResult, pointCnt, pspbuf);
		Result_sp = mxCreateNumericMatrix(spcnt, 2, mxINT32_CLASS, mxREAL);
		pResult_sp = (int32 *)mxGetData(Result_sp);
		for(m=0;m<spcnt; m++){
			pResult_sp[m] = pspbuf[m];
			pResult_sp[m + spcnt] = pspbuf[m + pointCnt];
		}
		plhs[1] = Result_sp;
		mxFree(pspbuf);
	}
}

void LinkPoints(double *px, double *py, double *pframe, 
				size_t pnum, double gap, double dist, 
				size_t imgWidth, size_t img_Height,
				int32 *presult){
	double *pBuf_x, *pBuf_y, *pBuf_frame;
	int32 *pBuf_index;
	int m,n,searchwidth=9;
	int tx,ty;
	int dx[] = {0, -1, 0, 0, 1, -1, -1, 1, 1}, dy[] = {0, 0, -1, 1, 0, -1, 1, -1, 1};
	int cx,cy,base, findflag, preidx;
	double tf, tdist;
	dist = dist*dist;
	//alloc mamory
	pBuf_x = (double *)mxMalloc(imgWidth*img_Height*sizeof(double));
	pBuf_y = (double *)mxMalloc(imgWidth*img_Height*sizeof(double));
	pBuf_frame = (double *)mxMalloc(imgWidth*img_Height*sizeof(double));
	pBuf_index = (int32 *)mxMalloc(imgWidth*img_Height*sizeof(int32));
	
	//init buffers
	for(m=0;m<img_Height; m++){
		for(n=0;n<imgWidth; n++){
			pBuf_index[n + m*imgWidth] = -1;
		}
	}
	for(m=0; m<pnum; m++){
		presult[m] = -1;
	}
	//search points
	for(m=0; m<pnum; m++){
		tx = int(px[m]);
		ty = int(py[m]);
		tf = pframe[m];
		findflag = 0;
		for(n=0;n<searchwidth; n++){
			cx = int(tx) + dx[n];
			cy = int(ty) + dy[n];
			base = cx + cy*imgWidth;
			if(cx>=0 && cy>=0 && cx<imgWidth && cy<img_Height && pBuf_index[base]>=0){
				tdist = (px[m] - pBuf_x[base])*(px[m] - pBuf_x[base]) 
						+ (px[m] - pBuf_x[base])*(px[m] - pBuf_x[base]);
				if(tdist < dist && abs(tf - pBuf_frame[base])<= gap && abs(tf - pBuf_frame[base])>0.5){
					findflag =1;
					preidx = pBuf_index[base]+1;
					pBuf_index[base] = -1;
					break;
				}
			}
		}
		if(findflag != 0){//find llink point
			presult[preidx] = m;
		}
		//reg point
		if(tx>=0 && ty>=0 && tx<imgWidth && ty<=img_Height){
			base = tx + ty*imgWidth;
			pBuf_x[base] = px[m];
			pBuf_y[base] = py[m];
			pBuf_frame[base] = pframe[m];
			pBuf_index[base] = m;
		}
	}
	//free memory
	mxFree(pBuf_x);
	mxFree(pBuf_y);
	mxFree(pBuf_frame);
	mxFree(pBuf_index);
}

int ListStartPoints(int32 *pLinkData, size_t pnum, int32 *presult){
	short *pFlagBuf = NULL;
	int m,n,cnt=0,deep;
	pFlagBuf = (short *)mxMalloc(pnum*sizeof(short));
	for(m=0;m<pnum;m++){
		pFlagBuf[m] = 0;
	}
	for(m=0;m<pnum;m++){
		if(pFlagBuf[m]==0){
			presult[cnt] = m+1;
			n=m;
			deep = 0;
			while(n>=0 && n<pnum){
				pFlagBuf[n] = 1;
				n = pLinkData[n]-1;
				deep++;
			}
			presult[cnt + pnum] = deep;
			cnt++;
		} 
	}
	mxFree(pFlagBuf);
	return cnt;
}