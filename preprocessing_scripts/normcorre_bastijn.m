function normcorre_bastijn(save_png)

% NoRMCorre
% function normcorre_bastijn(save_png)
% save_png: true to save figures
% rigid: true to save rigid video
% nonrigid: true to save nonrigid video
% Script is based on demo_1p.m
% applying the NoRMCorre motion correction algorithm on 
% 1-photon widefield imaging data
% adjusted for our miniscope recordings
% batch processing, runs faster, no plots inbetween
% adjusted by Bastijn van den Boom 2019-02-07

%% start parallel pool

close all;
gcp;

fprintf('Parallel pool started 1/16\n')

%% read data, convert to double, create directories

%name = './15407_DS_HF.avi';        %.avi fine as well
%addpath(genpath('../../NoRMCorre'));

normcorre_choose_data;

cd([dir_nm]);

mkdir NoRMCorre;
mkdir CNMF-E;

cd([dir_nm, '\NoRMCorre']);

%  info = imfinfo(name);

% imageStack=[];
% numberofImages = length(info);
% for k = 1:numberofImages
%     currentImage = imread(name,k,'Info',info);
%     Yf(:,:,k) = currentImage;
% end
 Yf = read_file(name);
Yf = single(Yf);
[d1,d2,T] = size(Yf);

fprintf('Data loaded and converted to matrix 2/16\n')
%% perform some sort of deblurring/high pass filtering

if (0)    
    hLarge = fspecial('average', 40);
    hSmall = fspecial('average', 2); 
    for t = 1:T
        Y(:,:,t) = filter2(hSmall,Yf(:,:,t)) - filter2(hLarge, Yf(:,:,t));
    end
    %Ypc = Yf - Y;
    bound = size(hLarge,1);
else
    gSig = 7; 
    gSiz = 17; 
    psf = fspecial('gaussian', round(2*gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    %Y = imfilter(Yf,psf,'same');
    %bound = 2*ceil(gSiz/2);
    Y = imfilter(Yf,psf,'symmetric');
    bound = 0;
end

fprintf('High pass filtering 3/16\n')
%% first try out rigid motion correction
    % exclude boundaries due to high pass filtering effects
options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);

fprintf('Rigid motion correction 4/16\n')
%% register using the high pass filtered data and apply shifts to original data
tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
    % exclude boundaries due to high pass filtering effects
tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset
    % apply shifts on the whole movie
    
fprintf('Register data and apply to original 5/16\n')
%% compute metrics 
[cY,mY,vY] = motion_metrics(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
[cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);

[cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
[cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);

fprintf('Compute rigid motion metrics 6/16\n')
%% plot rigid shifts and metrics

%save_png = true; %save png image of motion

shifts_r = squeeze(cat(3,shifts1(:).shifts));
figure;
    subplot(311); plot(shifts_r);
        title('Rigid shifts','fontsize',14,'fontweight','bold');
        legend('y-shifts','x-shifts');
    subplot(312); plot(1:T,cY,1:T,cM1);
        title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
        legend('raw','rigid');
    subplot(313); plot(1:T,cYf,1:T,cM1f);
        title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
        legend('raw','rigid');

if save_png
    file_name = file_nm;
    file_name(find(file_nm=='.',1,'last'):end) = [];    
    png_nm = [file_name, '_rigid_shift'];
    saveas(gcf, png_nm, 'bmp');
end

fprintf('Plot rigid movements in X, Y, Z direction 7/16\n')
%% now apply non-rigid motion correction
% non-rigid motion correction is likely to produce very similar results
% since there is no raster scanning effect in wide field imaging

options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
    'grid_size',[128,128]*2,'mot_uf',4,'correct_bidir',false, ...
    'overlap_pre',32,'overlap_post',32,'max_shift',20);

tic; [M2,shifts2,template2] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_nr,template1); toc % register filtered data
tic; Mpr = apply_shifts(Yf,shifts2,options_nr,bound/2,bound/2); toc % apply the shifts to the removed percentile

fprintf('Apply non-rigid correction 8/16\n')
%% compute metrics

[cM2,mM2,vM2] = motion_metrics(M2,options_nr.max_shift);
[cM2f,mM2f,vM2f] = motion_metrics(Mpr,options_nr.max_shift);

fprintf('Compute non-rigid motion corrections 9/16\n')
%% plot shifts        

%save_png = true; %save png image of motion

shifts_r = squeeze(cat(3,shifts1(:).shifts));
shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_nr(:,2,:))';
shifts_y = squeeze(shifts_nr(:,1,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients for filtered data','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[],'XLim',[0,T-3])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

if save_png
    file_name = file_nm;
    file_name(find(file_nm=='.',1,'last'):end) = [];    
    png_nm = [file_name, '_nonrigid_shift'];
    saveas(gcf, png_nm, 'bmp');
end

fprintf('Plot shift after correcting 10/16\n')
%% display downsampled data

tsub = 1;       %5

Y_ds = downsample_data(Y,'time',tsub);
Yf_ds = downsample_data(Yf,'time',tsub);
M1_ds = downsample_data(M1,'time',tsub);
M1f_ds = downsample_data(Mr,'time',tsub);
M2_ds = downsample_data(M2,'time',tsub);
M2f_ds = downsample_data(Mpr,'time',tsub);
nnY_ds = quantile(Y_ds(:),0.0005);
mmY_ds = quantile(Y_ds(:),0.9995);
nnYf_ds = quantile(Yf_ds(:),0.0005);
mmYf_ds = quantile(Yf_ds(:),0.99995);

fprintf('Display downsampled data if required 11/16\n')


%% save nonrigid corrected video in .tif for CNMF-E 
 
% Important! make sure directory is not open in windows explorer!
 
fprintf('If you get the error "Error using imwrite (line 454)", you probably have the directory open in windows explorer.\n')
fprintf('Close the directory and run this part again.\n')

cd([dir_nm, '\CNMF-E']);
vid_nm = [file_name, '_NoRMCorre', '.tif'];

%Rigid

% if rigid
%     M1f_ds_tif = M1f_ds;
%     vid_nm = [file_name, '_NormCorreRigid', '.tif'];
%     M1f_ds_tif = mat2gray(M1f_ds_tif);
%         for i = 1 : size(M1f_ds_tif,3);
%             imwrite(M1f_ds_tif(:,:,i),vid_nm,'Compression','none','WriteMode','append');
%         end
% end

%Non-rigid
% if nonrigid
%     M2f_ds_tif = M2f_ds;
%     vid_nm = [file_name, '_NormCorreNonrigid', '.tif'];
%     M2f_ds_tif = mat2gray(M2f_ds_tif);
%     for i = 1 : size(M2f_ds_tif,3);
%         imwrite(M2f_ds_tif(:,:,i),vid_nm,'Compression','none','WriteMode','append');
%     end
% end
 
% lets go for fast_tiff_write! aint nobody got time for waiting to save a tiff

% 3d matrix to 2d cell vector & double to uint8
M2f_ds = mat2gray(M2f_ds);

for i = 1:size(M2f_ds,3)
    M2f_ds_tif{i} = im2uint8(M2f_ds(:,:,i));
end

tic

fTIF = Fast_Tiff_Write(vid_nm);
for k = 1:length(M2f_ds_tif)
    fTIF.WriteIMG(M2f_ds_tif{k}');
end
fTIF.close

fprintf('Saved tiff file in %.2s seconds. \n', toc)
fprintf('Saved to .tif. Lets go to CNMF-E! 14/16\n')



%% total time to run script

clear; clc; close all;  


end


