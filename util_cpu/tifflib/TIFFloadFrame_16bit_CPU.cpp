#include "mex.h"
#include "matrix.h"
#include ".\include\tiffio.h"
#include "windows.h"
/* imported function list:
 TIFFOpen
 TIFFClose
 TIFFReadEncodedStrip
 TIFFNumberOfStrips
 TIFFSetDirectory
 TIFFGetField
 TIFFNumberOfDirectories
 */
typedef TIFF* (* _TIFFOpen)(const char*, const char*);
typedef void (* _TIFFClose)(TIFF*);
typedef tsize_t (* _TIFFReadEncodedStrip)(TIFF*, tstrip_t, tdata_t, tsize_t);
typedef tstrip_t (* _TIFFNumberOfStrips)(TIFF*);
typedef int (* _TIFFSetDirectory)(TIFF*, tdir_t);
typedef int (* _TIFFGetField)(TIFF*, ttag_t, ...);
typedef tdir_t (* _TIFFNumberOfDirectories)(TIFF*);

HINSTANCE hDLL = NULL;
_TIFFOpen pTIFFOpen = NULL;
_TIFFClose pTIFFClose = NULL;
_TIFFReadEncodedStrip pTIFFReadEncodedStrip = NULL;
_TIFFNumberOfStrips pTIFFNumberOfStrips = NULL;
_TIFFSetDirectory pTIFFSetDirectory = NULL;
_TIFFGetField pTIFFGetField = NULL;
_TIFFNumberOfDirectories pTIFFNumberOfDirectories = NULL;


mxArray * GetImage(char * path, int imgidx);
int loadLibTiff();
void ReadFrame(TIFF *tif, int imgidx, uint16 *pbuf, 
        uint32 imageRowsPerStrip, uint32 imageWidth, uint32 imageLength);
mxArray * GetImageStack(char * path, int idxstart, int idxend);

void mexFunction(
     int nlhs, mxArray *plhs[],
     int nrhs, const mxArray *prhs[])
{
    int idx1=0,idx2=0;
    char path[1000];
    
    if(loadLibTiff()<0){
        return;
    }
    if(nrhs==0 || mxGetClassID(prhs[0])!= mxCHAR_CLASS){
        mexPrintf("Need file path");
        return;
    } else{
        mxGetNChars(prhs[0], path, 1000);
    }
    
    if(nrhs>=2){
        if(mxGetNumberOfElements(prhs[1]) ==1){
            idx1 = int(mxGetPr(prhs[1])[0])-1;
            plhs[0] = GetImage(path, idx1);
        } else{
            idx1 = int(mxGetPr(prhs[1])[0])-1;
            idx2 = int(mxGetPr(prhs[1])[1])-1;
            plhs[0] = GetImageStack(path, idx1, idx2);
        }
    } else{
        plhs[0] = GetImage(path, 0);
    }
    FreeLibrary(hDLL);
}

mxArray * GetImage(char * path, int imgidx){
    TIFF *tif=NULL;
	int stripsize, stripnum, imagesize;
	uint32 imageWidth, imageLength, TileWidth, TileLength, imageRowsPerStrip ;
	uint16 *pdata = NULL, *pResult;
    mxArray *Result = NULL;
	int m,n;
    
    tif = pTIFFOpen(path, "r");
    
	if(!tif){
		mexPrintf("open file failed!\n");
	} else{
        pTIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &imageRowsPerStrip);
        pTIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
        pTIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
        pTIFFSetDirectory(tif, imgidx);
        stripnum = pTIFFNumberOfStrips(tif);
        imagesize = imageWidth * imageLength;

        pdata = (uint16 *)malloc(imagesize * sizeof(uint16));


        n=0;//row cnt
        for(m=0;m<stripnum; m++){
            pTIFFReadEncodedStrip(tif, m, (void *)(pdata + m*imageRowsPerStrip*imageWidth), imageRowsPerStrip*imageWidth*sizeof(uint16));
        }
        
        pTIFFClose(tif);
        
        Result = mxCreateNumericMatrix(imageLength, imageWidth, mxUINT16_CLASS, mxREAL);
        pResult = (uint16 *)mxGetData(Result);
        for(m=0;m<imageLength; m++){
            for(n=0; n<imageWidth; n++){
                pResult[n*imageLength+m] = pdata[m*imageWidth +n];
            }
        }
        free(pdata);
    }
    return Result;
}

mxArray * GetImageStack(char * path, int idxstart, int idxend){
    TIFF *tif=NULL;
	int stripsize, stripnum, imagesize;
	uint32 imageWidth, imageLength, imageRowsPerStrip ;
	uint16 *pdata = NULL, *pResult;
    mxArray *Result = NULL;
	int m,n,base,idxcnt;
    mwSize reSize[3];
    
    
    tif = pTIFFOpen(path, "r");
    if(idxend<0){
        idxend = pTIFFNumberOfDirectories(tif)-1;
    }else if(idxend<idxstart){
        idxend = idxstart;
    }
    //mexPrintf("s:%d,e:%d,NOD:%d\n", idxstart, idxend, pTIFFNumberOfDirectories(tif));
    
	if(!tif){
		mexPrintf("open file failed!\n");
	} else{
        pTIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &imageRowsPerStrip);
        pTIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
        pTIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
        
        stripnum = pTIFFNumberOfStrips(tif);
        imagesize = imageWidth * imageLength;
        reSize[0] = imageLength;
        reSize[1] = imageWidth;
        reSize[2] = idxend-idxstart+1;
        
        Result = mxCreateNumericArray(3, reSize, mxUINT16_CLASS, mxREAL);
        pResult = (uint16 *)mxGetData(Result);
        pdata = (uint16 *)malloc(imagesize * sizeof(uint16)); 
        
        for(idxcnt = idxstart; idxcnt<=idxend; idxcnt++){
//             mexPrintf("%d\n",idxcnt);
            n=0;//row cnt
            ReadFrame(tif, idxcnt, pdata,imageRowsPerStrip, imageWidth, imageLength);
            base = idxcnt*imageLength*imageWidth;
            for(m=0;m<imageLength; m++){
                for(n=0; n<imageWidth; n++){
                    pResult[base + n*imageLength+m] = pdata[m*imageWidth +n];
                }
            }
        } 
        pTIFFClose(tif);
        free(pdata);
    }
    return Result;
}


//read one frame
void ReadFrame(TIFF *tif, int imgidx, uint16 *pbuf, 
        uint32 imageRowsPerStrip, uint32 imageWidth, uint32 imageLength){
	int stripsize, stripnum, imagesize;
	int m;

    pTIFFSetDirectory(tif, imgidx);
    stripnum = pTIFFNumberOfStrips(tif);
    imagesize = imageWidth * imageLength;

    for(m=0;m<stripnum; m++){
        pTIFFReadEncodedStrip(tif, m, (void *)(pbuf + m*imageRowsPerStrip*imageWidth), imageRowsPerStrip*imageWidth*sizeof(uint16));
    }
}

int loadLibTiff(){
    HINSTANCE hDLL;
    hDLL=LoadLibrary("libtiff.dll");
    if(hDLL){
//         mexPrintf("%d\n",
         pTIFFOpen = (_TIFFOpen)GetProcAddress(hDLL,"TIFFOpen");
//         mexPrintf("%d\n",
         pTIFFClose = (_TIFFClose)GetProcAddress(hDLL,"TIFFClose");
//         mexPrintf("%d\n",
         pTIFFReadEncodedStrip = (_TIFFReadEncodedStrip)GetProcAddress(hDLL,"TIFFReadEncodedStrip");
//         mexPrintf("%d\n",
         pTIFFNumberOfStrips = (_TIFFNumberOfStrips)GetProcAddress(hDLL,"TIFFNumberOfStrips");
//         mexPrintf("%d\n",
         pTIFFSetDirectory = (_TIFFSetDirectory)GetProcAddress(hDLL,"TIFFSetDirectory");
//         mexPrintf("%d\n",
         pTIFFGetField = (_TIFFGetField)GetProcAddress(hDLL,"TIFFGetField"); 
//         mexPrintf("%d\n",
         pTIFFNumberOfDirectories = (_TIFFNumberOfDirectories)GetProcAddress(hDLL,"TIFFNumberOfDirectories"); 
//          FreeLibrary(hDLL);
         return 0;
    } else{
        mexPrintf("Can not find libtiff.dll!\n");
        return -1;
    }
}
