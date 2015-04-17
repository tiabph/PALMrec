function [databuf, timeResult, param] = PALMrec_cmd(inifile, imgfile, varargin)
%% PALM reconstruction script
%init path
addpath('.\detection');
addpath('.\PALM');
addpath('.\common');
addpath('.\util_cpu');
addpath('.\drift');
addpath('.\cmdline');
%init parameters
% inifile = '.\cmdline\default.ini';
param = ReadINI(inifile);
if(nargin>1)
    param = parseArgin(varargin, param);
end
if ~isempty(imgfile)
    param.fullpath = imgfile;
end
param = CheckParam(param);

timestamp = tic();
%% load image data
tic
% param.filepath = filepath;
% param.filename = filename;
databuf = palmRec_LoadImage_cpu([], param);
disp(['--load image data time: ' num2str(toc)])
timeResult.load = toc();

%% detection
tic
% param.detection.threshold = 2;
% param.detection.type = 1;
% param.detection.type_str = 'W2';
% param.detection.windowWidth = 5;

databuf = palmRec_FindParticles_cpu(databuf, param);
disp(['--detection time: ' num2str(toc)])
timeResult.detection = toc();

%% create sub_img (ROI), format [frame, y, x, height]
tic
% param.fitting.fitl = 3;%half width

databuf = palmRec_CreateSubImages_cpu(databuf, param);
disp(['--create sub_img time: ' num2str(toc)])
timeResult.createROI = toc();

%% fit sub_img
temptime = tic;
blockSize = 10000;
% param.fitting.factor = 1;
% param.fitting.gain = 1;
% param.fitting.sigma = 1.5;

databuf = palmRec_Fitting(databuf, param);
disp(['--fit sub_img time: ' num2str(toc(temptime))])
timeResult.fitting = toc(temptime);

%% post process fit info
tic
% param.fitting.pixelsize = 100;

databuf = palmRec_PostProcessFittingData_cpu(databuf, param);
disp(['--post process fit info time: ' num2str(toc)])
timeResult.postfit = toc();

%% drift correction
tic
% param.drift.type = 'file';
% param.drift.path = param.filepath;
% param.drift.file = beadsfile;

databuf = palmRec_DriftCorrection(databuf, param);

disp(['--drift correction time: ' num2str(toc)])
timeResult.drift = toc();

%% linking
tic
% param.linking.gap = 1;
% param.linking.neighbor = 8;
databuf = palmRec_Linking_cpu(databuf, param);
disp(['--linking time: ' num2str(toc)])
timeResult.linking = toc();
databuf = palmRec_Filter(databuf, param);

%% reconstruction
tic
% param.reconstruction.amp = 16;
% param.reconstruction.recType = 'DOT';%DOT or SPOT

databuf = palmRec_Reconstruction(databuf, param);
disp(['--reconstruction time: ' num2str(toc)])
timeResult.reconstruction = toc();

%% write image file
tic
databuf = palmRec_Save(databuf, param);
disp(['--write image file time: ' num2str(toc)])
timeResult.write = toc();
timeResult.total = toc(timestamp);

disp('--------------- TIME REPORT ---------------');
disp(timeResult);

end