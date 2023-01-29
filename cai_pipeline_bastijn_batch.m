 %% BATCH Script to analyze calcium imaging data

% NoRMCorre motion correction and CNMF-E for cell extraction
% run all animals in 1 folder, you only have to select main folder

% directory structure: animal folder - .tiff file
% function will create two directories in the animal folder
% one to put output NoRMCorre in (if true), other to put CNMF-E video and
% results in.

% adjustment in batch for CNMF-E is that it will use try-catch until it
% succesfully analized CNMF-E. The problem was that we had regularly errors
% due to min_pnr being too low. Before you realize this occured, you lost
% hours of analyses time. Therefore, the current script increases the pnr
% until no error occurs (or until it is out of values!).

% the vector min_pnr_vec will be used to run CNMF-E until no error occurs

%normcorre_bastijn(save_png)
% save_png: true to save figures

% second step is CNMF-E
% it will choose the data from the previous step
% you will have to specify some variables

%cnmf_e_bastijn_batch(Fs,ssub,tsub,gSig,gSiz,min_corr,min_pnr,bd,saveim)
%Fs = 15;            % frame rate
%ssub = 2;           % spatial downsampling factor
%tsub = 1;           % temporal downsampling factor IMPORTANT, use 1 in the end
%gSig = 10;           % [mouse deep brain: 7, rat 10] width of the gaussian kernel, which can approximates the average neuron shape
%gSiz = 17;          % [mouse deep brain: 17, rat 25]maximum diameter of neurons in the image plane. larger values are preferred.
%min_corr = 0.8;     % minimum local correlation for a seeding pixel [str2double(response(1))]
%min_pnr = 7.4;       % minimum peak-to-noise ratio for a seeding pixel [str2double(response(2))]
%bd = 10;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
%saveim = true;      %save images


%double check function
%why does cnmf-e show incorrect min_pnr?

function cai_pipeline_bastijn_batch

clear all

% get rid of warnings
warning('off')

% vector with PNR values to try
pnr_vec = [4 5 6 7.4 9 11 15 18 20 23 26 30 35 40]; 
corr_vec = [.8 .8 .8 .8 .8 .8 .8 .8 .8 .8 .8 .85 .85 .9]; 

tic;

% ask user to select main directory
dialog_title = 'Select directory with individual animal folders with inside .tiff files';
name_dir = uigetdir('',dialog_title);
cd([name_dir]);

%find animal folders
names = dir(name_dir);
names(ismember( {names.name}, {'.', '..'})) = []; % delete . and ..

%go to folder and run for loop
for i = 1:length(names)
    
    clearvars -except i pnr_vec names corr_vec
    close all
    
    dir_nm = [names(i).folder '\' names(i).name];
    cd(dir_nm);
    files_dir = dir('*.tif');
    file_nm = files_dir.name;
    name = [dir_nm '\' files_dir.name];
    
    
    % run NoRMCorre
    normcorre_bastijn_batch(dir_nm, name, file_nm, true);
    
    % go to folder with tiff file
    cd([dir_nm, '\CNMF-E']);    
    
    % run CNMF-E
    clearvars -except i pnr_vec names corr_vec
    
    j = 1;
	try
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %1
    catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %2
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %3
    catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %4
    catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %5
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %6
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %7
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %8
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %9
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %10
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %11
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %12
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %13
	catch; try
        j=j+1;
        cnmf_e_bastijn_batch(15,2,1,7,17,corr_vec(j),pnr_vec(j),5,true);    %14
    catch
        disp('Could not run CNMF-E, please run manually')
    end; end; end; end; end; end; end; end; end; end; end; end; end; end
end

toc

fprintf('Total time to run batch imaging pipeline:     %.2f seconds\n', toc);

