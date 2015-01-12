%% PALM reconstruction script
%init path
addpath('.\detection');
addpath('.\PALM');
addpath('.\common');
%init parameters
filepath = '.\test\';
filename = 'a647_big.tif';
filename_out = 'a647_big_rec.tif';
outpath = '.\test\';

param.filepath = filepath;
param.filename = filename;

%% load image data
tic
databuf = palmRec_LoadImage([], param);
disp(['--load image data time: ' num2str(toc)])

%% detection
tic
param.detection.threshold = 2;
param.detection.type = 1;
param.detection.windowWidth = 5;

databuf = palmRec_FindParticles(databuf, param);
disp(['--detection time: ' num2str(toc)])


%% create sub_img (ROI), format [frame, y, x, height]
tic
param.fitting.fitl = 3;

databuf = palmRec_CreateSubImages(databuf, param);
disp(['--create sub_img time: ' num2str(toc)])

%% fit sub_img
temptime = tic;
blockSize = 10000;
param.fitting.factor = 1;
param.fitting.gain = 1;
param.fitting.sigma = 1.5;

databuf = palmRec_Fitting(databuf, param);
disp(['--fit sub_img time: ' num2str(toc(temptime))])

%% post process fit info
tic
param.fitting.pixelsize = 100;

databuf = palmRec_PostProcessFittingData(databuf, param);
disp(['--post process fit info time: ' num2str(toc)])

%% linking
tic
param.linking.gap = 1;
param.linking.neighbor = 8;
databuf = palmRec_Linking(databuf, param);
disp(['--linking time: ' num2str(toc)])

%% reconstruction
tic
param.reconstruction.amp = 16;
param.reconstruction.recType = 'DOT';%DOT or SPOT

databuf = palmRec_Reconstruction(databuf, param);
disp(['--reconstruction time: ' num2str(toc)])

%% write image file
tic
databuf = palmRec_Save(databuf, param);
disp(['--write image file time: ' num2str(toc)])