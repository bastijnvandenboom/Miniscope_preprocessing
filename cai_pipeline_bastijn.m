%% Script to analyze calcium imaging data

% all scripts are based on NoRMCorre motion correction and
% CNMF-E for cell extraction

% first step is run NoRMCorre
% select the tiff file to analyze, function will create two directories
% one to put output NoRMCorre in (if true), other to put CNMF-E video and
% results in. 

tic;

normcorre_bastijn(true); 

%normcorre_bastijn(save_png)
% save_png: true to save figures

% second step is CNMF-E 
% it will choose the data from the previous step
% you will have to specify some variables

cnmf_e_bastijn(15,2,1,7,17,0.85,20,5,true);
%bulsara 20 pnr 10

%cnmf_e_bastijn(15,2,1,7,17,0.8,7.4,5,true);
%cnmf_e_bastijn(15,2,1,7,17,0.85,20,10,true) working on this for 14481

%Note: fix true (one should write explicit true)

%cnmf_e_bastijn(Fs,ssub,tsub,gSig,gSiz,min_corr,min_pnr,bd,saveim)
%Fs = 15;            % frame rate
%ssub = 2;           % spatial downsampling factor
%tsub = 1;           % temporal downsampling factor IMPORTANT, use 1 in the end
%gSig = 10;           % [mouse deep brain: 7, rat 10] width of the gaussian kernel, which can approximates the average neuron shape
%gSiz = 17;          % [mouse deep brain: 17, rat 25]maximum diameter of neurons in the image plane. larger values are preferred.
%min_corr = 0.8;     % minimum local correlation for a seeding pixel [str2double(response(1))]
%min_pnr = 7.4;       % minimum peak-to-noise ratio for a seeding pixel [str2double(response(2))]
%bd = 10;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
%saveim = true;      %save images


fprintf('Total time to run imaging pipeline:     %.2f seconds\n', toc);