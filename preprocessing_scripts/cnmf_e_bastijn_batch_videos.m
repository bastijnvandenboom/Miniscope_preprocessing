function cnmf_e_bastijn_batch_videos(Fs,ssub,tsub,gSig,gSiz,min_corr,min_pnr,bd,saveim,savevid)

% CNMF-E analysis for Miniscope recordings
% function
% cnmf_e_bastijn(Fs,ssub,tsub,gSig,gSiz,min_corr,min_pnr,bd,saveim)
% Fs = frame rate
% ssub = spatial downsampling factor
% tsub = temporal downsampling factor IMPORTANT, use 1 in the end
% gSig = width of the gaussian kernel, which can approximates the average neuron shape
% gSiz = maximum diameter of neurons in the image plane. larger values are preferred.
% saveim = save images
% min_corr = minimum local correlation for a seeding pixel [str2double(response(1))]
% min_pnr = minimum peak-to-noise ratio for a seeding pixel [str2double(response(2))]
% bd = number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
% saveim = save images
% savevid = save video of denoised and demixed signals
% script based on cnmfe_v2 from Tycho
% adjusted by Bastijn van den Boom 2019-02-07

% Fs = 15;            % frame rate
% ssub = 2;           % spatial downsampling factor
% tsub = 1;           % temporal downsampling factor IMPORTANT, use 1 in the end
% gSig = 10;           % [mouse deep brain: 7, rat 10] width of the gaussian kernel, which can approximates the average neuron shape
% gSiz = 17;          % [mouse deep brain: 17, rat 25]maximum diameter of neurons in the image plane. larger values are preferred.
% saveim = true;        %save images
% min_corr = 0.8;     % minimum local correlation for a seeding pixel [str2double(response(1))]
% min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel [str2double(response(2))]
% bd = 10;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
% savevid = true        % save video

tic;

%% clear workspace
close all;  
global  d1 d2 numFrame ssub tsub sframe num2read Fs neuron neuron_ds ...
    neuron_full Ybg_weights; %#ok<NUSED> % global variables, don't change them manually

fprintf('Continue! 1\n')

%% select data and map it to the RAM

global data Ysiz d1 d2 numFrame; 

%% select file 

dir_nm = pwd
find_file = dir(fullfile(dir_nm, '*NoRMCorre.tif'))
nam = ([dir_nm '\' find_file.name])

    % use pre-specified file 
    if exist(nam, 'file')
        [dir_nm, file_nm, file_type] = fileparts(nam);
    else
        dir_nm = 0; 
    end



%% convert the data to mat file
nam_mat = [dir_nm, filesep, file_nm, '.mat'];
if strcmpi(file_type, '.mat')
    fprintf('The selected file is *.mat file\n');
elseif  exist(nam_mat, 'file')
    % the selected file has been converted to *.mat file already
    fprintf('The selected file has been replaced with its *.mat version.\n');
elseif or(strcmpi(file_type, '.tif'), strcmpi(file_type, '.tiff'))
    % convert
    tic;
    fprintf('converting the selected file to *.mat version...\n');
    nam_mat = tif2mat(nam);
    fprintf('Time cost in converting data to *.mat file:     %.2f seconds\n', toc);
elseif or(strcmpi(file_type, '.h5'), strcmpi(file_type, '.hdf5'))
    fprintf('the selected file is hdf5 file\n'); 
    temp = h5info(nam);
    dataset_nam = ['/', temp.Datasets.Name];
    dataset_info = h5info(nam, dataset_nam);
    dims = dataset_info.Dataspace.Size;
    ndims = length(dims);
    d1 = dims(2); 
    d2 = dims(3); 
    numFrame = dims(end);
    Ysiz = [d1, d2, numFrame]; 
    fprintf('\nThe data has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1, d2, numFrame, prod(Ysiz)*8/(2^30));
    return; 
else
    fprintf('The selected file type was not supported yet! email me to get support (zhoupc1988@gmail.com)\n');
    return;
end

%% information of the data 
data = matfile(nam_mat);
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames

fprintf('\nThe data has been mapped to RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1, d2, numFrame, prod(Ysiz)*8/(2^30));

cd([dir_nm]);

fprintf('Continue! 2\n')

%% create Source2D class object for storing results and parameters
%Fs = 15;            % frame rate
%ssub = 2;           % spatial downsampling factor
%tsub = 1;           % temporal downsampling factor IMPORTANT, use 1 in the end
%gSig = 10;           % [mouse deep brain: 7, rat 10] width of the gaussian kernel, which can approximates the average neuron shape
%gSiz = 17;          % [mouse deep brain: 17, rat 25]maximum diameter of neurons in the image plane. larger values are preferred.
neuron_full = Sources2D('d1',d1,'d2',d2, ... % dimensions of datasets
    'ssub', ssub, 'tsub', tsub, ...  % downsampleing
    'gSig', gSig,...    % sigma of the 2D gaussian that approximates cell bodies
    'gSiz', gSiz);      % average neuron size (diameter)
neuron_full.Fs = Fs;         % frame rate
% neuron_full.Fs = Fs/tsub; 
tau_adjusted = (Fs); % to accomodate variable frame rate decay

% with dendrites or not 
with_dendrites = false;  %true necessary for Tycho A_C_display
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    neuron_full.options.search_method = 'dilate'; 
    neuron_full.options.bSiz = 80;
else
    % determine the search locations by selecting a round area
    neuron_full.options.search_method = 'ellipse';
    neuron_full.options.dist = 5;
end

merge_thr = [1e-1, 0.85, 0];     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
dmin = 1;

fprintf('Continue! 3\n')

%% options for running deconvolution 
neuron_full.options.deconv_flag = true; 
neuron_full.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'optimize_smin', true, ...
    'max_tau', tau_adjusted);    % maximum decay time (unit: frame);

fprintf('Continue! 4')

%% downsample data for fast and better initialization
sframe=1;						% user input: first frame to read (optional, default:1)
num2read= numFrame;             % user input: how many frames to read   (optional, default: until the end)

tic;
cnmfe_load_data;
fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);

Y = neuron.reshape(Y, 1);       % convert a 3D video into a 2D matrix

fprintf('Continue! 5\n')

%% compute correlation image and peak-to-noise ratio image.
cnmfe_show_corr_pnr;    % this step is not necessary, but it can give you some...
                        % hints on parameter selection, e.g., min_corr & min_pnr
% ask for min_corr and min_pnr values

% prompt={'Enter min_corr:','Enter min_pnr:'};
% name='Input min corr and min pnr';
% numlines=1;
% defaultanswer={'0.8','7.4'};
% response = inputdlg(prompt,name,numlines,defaultanswer);

fprintf('Continue! 6\n')

%% initialization of A, C

% play around with the variable bd

%saveim = true; %save images

% parameters
debug_on = false;   % visualize the initialization procedue. 
save_avi = false;   %save the initialization procedure as an avi movie. 
patch_par = [1,1]*1; %1;  % divide the optical field into m X n patches and do initialization patch by patch. It can be used when the data is too large 
K = []; % maximum number of neurons to search within each patch. you can use [] to search the number automatically

%min_corr = 0.8;     % minimum local correlation for a seeding pixel [str2double(response(1))]
%min_pnr = 10;       % minimum peak-to-noise ratio for a seeding pixel [str2double(response(2))]
min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
%bd = 10;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
neuron.updateParams('min_corr', min_corr, 'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, 'bd', bd);
neuron.options.nk = 1;  % number of knots for detrending 

% greedy method for initialization
tic;
[center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi);
fprintf('Time cost in initializing neurons:     %.2f seconds\n', toc);

% show results
figure;
imagesc(Cn, [0.1, 0.95]);
hold on; plot(center(:, 2), center(:, 1), 'or');
colormap; axis off tight equal;
title(sprintf('Correlation: %02.2f, PNR: %02.2f, gSig: %02.2f, & gSiz: %02.2f',min_corr,min_pnr,gSig,gSiz),'fontweight','bold','fontsize',14); 

if saveim
    Corr_nm = [file_nm, '_Correlation'];
    saveas(gcf, Corr_nm, 'bmp');
end

% sort neurons
[~, srt] = sort(max(neuron.C, [], 2), 'descend');
neuron.orderROIs(srt);
neuron_init = neuron.copy();

fprintf('Continue! 7\n')

%% iteratively update A, C and B
% parameters, merge neurons
display_merge = false;          % visually check the merged neurons
view_neurons = false;           % view all neurons

% parameters, estimate the background
spatial_ds_factor = 2;      % spatial downsampling factor. it's for faster estimation
thresh = 18;     % threshold for detecting frames with large cellular activity. (mean of neighbors' activity  + thresh*sn)

bg_neuron_ratio = 1.5;  % spatial range / diameter of neurons

% parameters, estimate the spatial components
update_spatial_method = 'hals';  % the method for updating spatial components {'hals', 'hals_thresh', 'nnls', 'lars'}
Nspatial = 50;       % this variable has different meanings: 
                    %1) udpate_spatial_method=='hals' or 'hals_thresh',
                    %then Nspatial is the maximum iteration 
                    %2) update_spatial_method== 'nnls', it is the maximum
                    %number of neurons overlapping at one pixel 
               
% parameters for running iteratiosn 
nC = size(neuron.C, 1);    % number of neurons 

maxIter = 1;        % maximum number of iterations 
miter = 1; 
while miter <= maxIter
    %% merge neurons
    cnmfe_quick_merge;              % run neuron merges
    cnmfe_merge_neighbors;          % merge neurons if two neurons' peak pixels are too close 
    
    %% udpate background 
    % estimate the background
    tic;
    cnmfe_update_BG;
    fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
    % neuron.playMovie(Ysignal); % play the video data after subtracting the background components.
    
    %% update spatial & temporal components
    tic;
    for m=1:2  
        %temporal
        neuron.updateTemporal_endoscope(Ysignal);
        cnmfe_quick_merge;              % run neuron merges
        %spatial
        neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
        
        % post process the spatial components (you can run either of these two operations, or both of them)
        neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        %neuron.compactSpatial();    % run this line if neuron shapes are circular 
        cnmfe_merge_neighbors; 
        
        % stop the iteration when neuron numbers are unchanged. 
        if isempty(merged_ROI)
            break;
        end
    end
    fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);
    
    %% pick neurons from the residual (cell 4).
    if miter==1
%         seed_method = 'auto'; 
%         [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par, seed_method, debug_on); % method can be either 'auto' or 'manual'
        neuron.options.seed_method = 'auto'; % methods for selecting seed pixels {'auto', 'manual'}
        [center_new, Cn_res, pnr_res] = neuron.pickNeurons(Ysignal - neuron.A*neuron.C, patch_par); % method can be either 'auto' or 'manual'
    
    end
    
    %% stop the iteration 
    temp = size(neuron.C, 1); 
    if or(nC==temp, miter==maxIter)
        break; 
    else
        miter = miter+1; 
        nC = temp; 
    end
end

fprintf('Continue! 8\n')

%% Here we will ask whether to apply results to full resolution


flag = 1;

% apply results to the full resolution
if(flag)
    if or(ssub>1, tsub>1)
        neuron_ds = neuron.copy();  % save the result
        neuron = neuron_full.copy();
        cnmfe_full;
        neuron_full = neuron.copy();
    end
else
end

% delete some neurons and run CNMF-E iteration 
% neuron.orderROIs('decay_time');  % you can also use {'snr', 'mean', 'decay_time'} 
% neuron.viewNeurons([], neuron.C_raw); 
tic;
cnmfe_update_BG;
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc);
%update spatial & temporal components
tic;

for m=1:2
    %temporal
    neuron.updateTemporal_endoscope(Ysignal);
    cnmfe_quick_merge;              % run neuron merges

    %spatial
    neuron.updateSpatial_endoscope(Ysignal, Nspatial, update_spatial_method);
    neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    neuron.compactSpatial(); 
    cnmfe_merge_neighbors; 
end
fprintf('Time cost in updating spatial & temporal components:     %.2f seconds\n', toc);

b0 = mean(Y,2)-neuron.A*mean(neuron.C,2); 
Ybg = bsxfun(@plus, Ybg, b0-mean(Ybg, 2)); 
neuron.orderROIs('snr'); 

fprintf('Continue! 9\n')

%% display neurons
dir_neurons = sprintf('%s%s%s_neurons%s', dir_nm, filesep, file_nm, filesep);
if exist('dir_neurons', 'dir')
    temp = cd();
    cd(dir_neurons);
    delete *;
    cd(temp);
else
    mkdir(dir_neurons);
end
%neuron.viewNeurons([], neuron.C_raw, dir_neurons);
%close(gcf); 

fprintf('Continue! 10\n')


%% save results 

Cn = imresize(Cn,[d1,d2]);
neuron.Cn = Cn;
results = neuron.obj2struct(); 
results_name = [dir_nm, filesep, file_nm, '_results.mat'];

results.options.bd = bd;
results.options.min_pnr = min_pnr;
results.options.min_corr = min_corr;
results.options.ssub = ssub;
results.options.tsub = tsub;
results.options.gSig = gSig;
results.options.gSiz = gSiz;
results.Fs = Fs;

save(results_name, 'results');

fprintf('Continue! 11\n')


%% save correlation map

%saveim = true;

figure;
imagesc(Cn, [0.1, 0.95]);
colormap; axis off tight equal;

if saveim
    Corr_nm = [file_nm, '_Correlation'];
    saveas(gcf, Corr_nm, 'bmp');
end


%% plot all traces

%saveim = true; %save images

%%% all identified traces (C_raw)%%%
figure
C_raw_norm = neuron.obj2struct.C_raw;
for i = 1: size(C_raw_norm, 1)
    C_raw_norm(i, :) = normalize(C_raw_norm(i, :));
end
plot((C_raw_norm + (1: size(C_raw_norm, 1))')')
axis tight
axis square
title('Identified traces (C raw)')

if saveim
    saveas(gcf, 'C_raw', 'bmp');
end


% plot all fitted traces (C)
figure
C_norm = neuron.obj2struct.C;
for i = 1: size(C_norm, 1)
    C_norm(i, :) = normalize(C_norm(i, :));
end
plot((C_norm + (1: size(C_norm, 1))')')
axis tight
axis square
title('Fitted traces (C)')

if saveim
    saveas(gcf, 'C', 'bmp');
end

fprintf('Continue! 13\n')

%% sort cells based on max spiking

%saveim = true; %save images

% cut some traces
C_raw_norm_cut = C_raw_norm;
% C_raw_norm_cut = C_raw_norm(1:10,:);
%C_raw_norm_cut([9],:) = [];

%  sort traces
[M,I] = max(C_raw_norm_cut, [], 2);
[~,idx] = sort(I, 'descend');
C_raw_norm_cut = C_raw_norm_cut(idx,:);

figure
plot((C_raw_norm_cut(1:size(C_raw_norm_cut),:) + (1:size(C_raw_norm_cut))')')
axis tight
axis square
title('Identified traces sorted based on max signal')

if saveim
    saveas(gcf, 'C_raw_sorted', 'bmp');
end

fprintf('Great job! 14\n')


%% make some videos

if ~exist('t_begin', 'var')
    t_begin = 1;
end
% if ~exist('t_end', 'var')
    t_end = size(neuron.C, 2);
% end
if ~exist('kt', 'var')
    kt = 1;
end

%% data preparation
Y = neuron.reshape(Y, 2); 
Yac = neuron.reshape(neuron.A*neuron.C, 2);
Ybg = neuron.reshape(Ybg, 2);
Ysignal = neuron.reshape(Ysignal, 2);
figure('position', [0,0, 600, 400]);

if ~exist('center_ac', 'var')
    center_ac = nanmedian(max(neuron.A,[],1)'.*max(neuron.C,[],2));
end
range_res = [-1,1]*center_ac;
if ~exist('range_ac', 'var')
    range_ac = center_ac*1.01+range_res;
end
if ~exist('range_Y', 'var')
    if ~exist('multi_factor', 'var')
        temp = quantile(Y(randi(numel(Y), 10000,1)), [0.01, 0.98]);
        multi_factor = ceil(diff(temp)/diff(range_ac));
%     else
%         temp = quantile(Y(randi(numel(Y), 10000,1)), 0.01);
    else
        temp = quantile(Y(randi(numel(Y), 10000,1)), 0.01); 
    end
    center_Y = temp(1) + multi_factor*center_ac;
    range_Y = center_Y + range_res*multi_factor;
end


%% create avi file for demixed video
if savevid == true
%     if ~exist('avi_filename', 'var')
        avi_filename =[dir_nm, filesep, file_nm, '_demixed'];
%     end
    avi_file = VideoWriter(avi_filename);
    if ~isnan(neuron.Fs)
        avi_file.FrameRate= neuron.Fs/kt;
    end
    avi_file.open();
end

%% add pseudo color to denoised signals
[K, T]=size(neuron.C);
% draw random color for each neuron
% tmp = mod((1:K)', 6)+1;
Y_mixed = zeros(neuron.options.d1*neuron.options.d2, T, 3);
temp = prism;
% temp = bsxfun(@times, temp, 1./sum(temp,2));
col = temp(randi(64, K,1), :);
for m=1:3
    Y_mixed(:, :, m) = neuron.A* (diag(col(:,m))*neuron.C);
end
Y_mixed = uint16(Y_mixed/(1*center_ac)*65536);
%% play and save
% ax_y =   axes('position', [0.015, 0.51, 0.3, 0.42]);
% ax_bg=   axes('position', [0.015, 0.01, 0.3, 0.42]);
% ax_signal=    axes('position', [0.345, 0.51, 0.3, 0.42]);
% ax_denoised =    axes('position', [0.345, 0.01, 0.3, 0.42]);
% ax_res =    axes('position', [0.675, 0.51, 0.3, 0.42]);
% ax_mix =     axes('position', [0.675, 0.01, 0.3, 0.42]);

for m=t_begin:kt:t_end
%     axes(ax_y); cla;
%     imagesc(Ybg(:, :,m)+Ysignal(:, :, m), range_Y);
%     %     set(gca, 'children', flipud(get(gca, 'children')));
%     title('Raw data');
%     axis equal off tight;
%     
%     axes(ax_bg); cla;
%     imagesc(Ybg(:, :, m),range_Y);
%     %     set(gca, 'children', flipud(get(gca, 'children')));
%     axis equal off tight;
%     title('Background');
%     
%     axes(ax_signal); cla;
%     imagesc(Ysignal(:, :, m), range_ac); hold on;
%     %     set(gca, 'children', flipud(get(gca, 'children')));
%     title(sprintf('(Raw-BG) X %d', multi_factor));
%     axis equal off tight;
% %     
%     axes(ax_denoised); cla;
%     imagesc(Yac(:, :, m), range_ac);
%     %     imagesc(Ybg(:, :, m), [-50, 50]);
%     title(sprintf('Denoised X %d', multi_factor));
%     axis equal off tight;
    
%     axes(ax_res); cla;
%     imagesc(Ysignal(:, :, m)-Yac(:, :, m), range_res);
%     %     set(gca, 'children', flipud(get(gca, 'children')));
%     title(sprintf('Residual X %d', multi_factor));
%     axis equal off tight;
%     %         subplot(4,6, [5,6,11,12]+12);
    
%     axes(ax_mix); cla;
    imagesc(neuron.reshape(Y_mixed(:, m,:),2));  hold on;
    title('Demixed');
%     text(1, 10, sprintf('Time: %.2f second', m/neuron.Fs), 'color', 'w', 'fontweight', 'bold');
    
    axis equal tight off;
    %     box on; set(gca, 'xtick', []);
    %     set(gca, 'ytick', []);
    
    drawnow(); 
    if savevid == true
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [400, 600]);
        avi_file.writeVideo(temp);
    end
end

if savevid == true
    avi_file.close();
end

fprintf('Great job! You saved demixed video. 14\n')


%% create avi file for denoised signal
if savevid == true
%     if ~exist('avi_filename', 'var')
        avi_filename =[dir_nm, filesep, file_nm, '_denoised'];
%     end
    avi_file = VideoWriter(avi_filename);
    if ~isnan(neuron.Fs)
        avi_file.FrameRate= neuron.Fs/kt;
    end
    avi_file.open();
end

%% play and save
% ax_y =   axes('position', [0.015, 0.51, 0.3, 0.42]);
% ax_bg=   axes('position', [0.015, 0.01, 0.3, 0.42]);
% ax_signal=    axes('position', [0.345, 0.51, 0.3, 0.42]);
% ax_denoised =    axes('position', [0.345, 0.01, 0.3, 0.42]);
% ax_res =    axes('position', [0.675, 0.51, 0.3, 0.42]);
% ax_mix =     axes('position', [0.675, 0.01, 0.3, 0.42]);

for m=t_begin:kt:t_end
%     axes(ax_y); cla;
%     imagesc(Ybg(:, :,m)+Ysignal(:, :, m), range_Y);
%     %     set(gca, 'children', flipud(get(gca, 'children')));
%     title('Raw data');
%     axis equal off tight;
%     
%     axes(ax_bg); cla;
%     imagesc(Ybg(:, :, m),range_Y);
%     %     set(gca, 'children', flipud(get(gca, 'children')));
%     axis equal off tight;
%     title('Background');
%     
%     axes(ax_signal); cla;
%     imagesc(Ysignal(:, :, m), range_ac); hold on;
%     %     set(gca, 'children', flipud(get(gca, 'children')));
%     title(sprintf('(Raw-BG) X %d', multi_factor));
%     axis equal off tight;
% %     
%     axes(ax_denoised); cla;
    imagesc(Yac(:, :, m), range_ac);
    %     imagesc(Ybg(:, :, m), [-50, 50]);
    title(sprintf('Denoised X %d', multi_factor));
    axis equal off tight;
    
%     axes(ax_res); cla;
%     imagesc(Ysignal(:, :, m)-Yac(:, :, m), range_res);
%     %     set(gca, 'children', flipud(get(gca, 'children')));
%     title(sprintf('Residual X %d', multi_factor));
%     axis equal off tight;
%     %         subplot(4,6, [5,6,11,12]+12);
    
%     axes(ax_mix); cla;
%     imagesc(neuron.reshape(Y_mixed(:, m,:),2));  hold on;
%     title('Demixed');
%      text(1, 10, sprintf('Time: %.2f second', m/neuron.Fs), 'color', 'w', 'fontweight', 'bold');
%     
%     axis equal tight off;
%     %     box on; set(gca, 'xtick', []);
%     %     set(gca, 'ytick', []);
    
    drawnow(); 
    if savevid == true
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [400, 600]);
        avi_file.writeVideo(temp);
    end
end

if savevid == true
    avi_file.close();
end

fprintf('Great job! You saved denoised video. 14\n')


%% total time to run script

fprintf('Total time to run CNMF-e:     %.2f seconds\n', toc);

end
