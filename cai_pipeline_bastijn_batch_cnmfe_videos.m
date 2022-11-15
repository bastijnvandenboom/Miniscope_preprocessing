 %% BATCH Script to analyze calcium imaging data

% CNMF-E for cell extraction and save demixed and denoised videos
% run all animals in 1 folder, you only have to select main folder

% directory structure: animal folder - CNMF-E - _NoRMCorre.tiff file
% Function will save CNMF-E video and results

% adjustment in batch for CNMF-E is that it will use try-catch until it
% succesfully analized CNMF-E. The problem was that we had regularly errors
% due to min_pnr being too low. Before you realize this occured, you lost
% hours of analyses time. Therefore, the current script increases the pnr
% until no error occurs (or until it is out of values!).

% the vector min_pnr_vec will be used to run CNMF-E until no error occurs

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

function cai_pipeline_bastijn_batch_cnmfe_videos

clear all

% vector with PNR values to try
pnr_vec = [5 6 7.4 11 15 18 20 23];    %[7.4 9 11 13 15 18 20 23]

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
    
    clearvars -except i pnr_vec names
    close all
    
    % go to folder with tiff file
    cd([names(i).folder '\' names(i).name '\CNMF-E']);    
    
    % run CNMF-E
    clearvars -except i pnr_vec names
    
    j = 1;
    try
        cnmf_e_bastijn_batch_videos(15,2,1,7,17,0.8,pnr_vec(j),5,true,true);    %1
    catch
        try
            j=j+1;
            cnmf_e_bastijn_batch_videos(15,2,1,7,17,0.8,pnr_vec(j),5,true,true);    %2
        catch
            try
                j=j+1;
                cnmf_e_bastijn_batch_videos(15,2,1,7,17,0.85,pnr_vec(j),5,true,true);    %3
            catch
                try
                    j=j+1;
                    cnmf_e_bastijn_batch_videos(15,2,1,7,17,0.85,pnr_vec(j),5,true,true);    %4
                catch
                    try
                        j=j+1;
                        cnmf_e_bastijn_batch_videos(15,2,1,7,17,0.85,pnr_vec(j),5,true,true);    %5
                    catch
                        try
                            j=j+1;
                            cnmf_e_bastijn_batch_videos(15,2,1,7,17,0.85,pnr_vec(j),5,true,true);    %6
                        catch
                            try
                                j=j+1;
                                cnmf_e_bastijn_batch_videos(15,2,1,7,17,0.85,pnr_vec(j),5,true,true);    %7
                            catch
                                j=j+1;
                                cnmf_e_bastijn_batch_videos(15,2,1,7,17,0.9,pnr_vec(j),5,true,true);    %8
                            end
                        end
                    end
                end
            end
        end
    end
end

toc

fprintf('Total time to run batch imaging pipeline:     %.2f seconds\n', toc);

