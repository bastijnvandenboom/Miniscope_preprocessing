%% Plot Cn (correlation plot)
Cn = results.Cn;

figure;
imagesc(Cn, [0.1, 0.95]);
colormap; axis off tight equal;



%% plot few cells
C_raw_norm_cut_max = results_cut.C_raw;
%  sort traces
[M,I] = max(C_raw_norm_cut_max, [], 2);
[~,idx] = sort(I, 'descend');
C_raw_norm_cut_max = C_raw_norm_cut_max(idx,:);

numCells = 30;
C_raw_norm_cut2 = C_raw_norm_cut_max(1:numCells, :);


figure
plot((C_raw_norm_cut2 + (1: size(C_raw_norm_cut2, 1))')')
axis tight
axis square
title('Sample identified traces of good cells (C raw)')

if saveim
    saveas(gcf, 'C_raw_goodCells_Max_sample30', 'bmp');
end

%% This part is only for heatmap

saveim = true; %save images

C_raw = results.C_raw;

%  sort traces
[M,I] = max(C_raw, [], 2);
[~,idx]  = sort(I, 'ascend');
C_raw = C_raw(idx,:);

C_raw_cut = C_raw;

%Cut out some shit cells
shit_cells = [18 34]; %79 80 95 96 %42
C_raw_cut(shit_cells,:) = [];

% Fill in
height = 71 ; 
start1 = 535;
length1 = 290;
start2 = 2147;
length2 = 290;

figure
imagesc(C_raw_cut),colormap('parula')  %jet or parula
rectangle('Position',[start1 0 length1 height])
rectangle('Position',[start2 0 length2 height])

if saveim
    saveas(hFig, 'C_raw_goodCells_heatmap', 'bmp');
end

%old stuff
% figure
% bwr = @(n)interp1([1 2 3], [0 0 1; 1 1 1; 1 0  0], linspace(1, 3, n), 'linear')
% hFig = imagesc(C_raw_cut);
% colormap(bwr(64));

% hHM = HeatMap(C_raw_cut, 'Colormap' ,jet)
% hFig = plot(hHM);

%% working on a new way to plot contours using variable A
% figure
% neuron.image(neuron.A(:, 1))    %plot 1 neurons
% neuron.viewNeurons(1);          %view spatial and temporal component (all of them neuron.viewNeurons([]))

A = results.A;
C = results.C; 
C_raw = results.C_raw;
S = results.S;

% % to fix
% figure1
% neuron.image(A(:,1))


%% normalize raw signals to look at the cells

%%% Normalize all identified traces (C_raw)%%%
C_raw_norm = C_raw;

for i = 1: size(C_raw_norm, 1)
    C_raw_norm(i, :) = normalize(C_raw_norm(i, :));
end

C_raw_norm_cut = C_raw_norm;

%% check all cells and delete if shit

%C_raw_norm_cut =  results_cut.C_raw;
i = 1;

%C_raw_norm_cut = C_raw_cut;


for i = 1:10:size(C_raw_norm_cut,1)
    C_raw_norm_part = C_raw_norm_cut(i:i+9,:);
    
    figure
    plot((C_raw_norm_part + (1: size(C_raw_norm_part, 1))')')
    axis tight
    
    axis square
    title(i)
end

% find these values! what are the shitty cells?
shit_cells = [1 2 3 6 10 11 17 19 21 22 24 27 29 30 31 38 39 40 45 47 48 49 51 52 54 55 61 63 66 67 68 72 73 75 80 81 83 84 87 88 89 90 92 93 94 95 100 101 103 104 105 107 108 111 114 116 117 118 120 121 122 123 124 125 127 128 129 130 131:136 138 139 141:144 147:163] ;

%Conny rec1 [82 91 101 104 105 107 108 109 112 116]
%Eevee rec1 [71 73 95 115 116]
%Nessie rec1 [44 151 161]
%Conny rec4 [37 57 58 59 61 62 64 65 66 67 68 69 71]
%Conny rec5 [72 73 74 79 80 90 92 97 103 108 109 110]
%Eevee rec2 [90 102 128 137 138 142 144 146 150 151 152 153 154 155 156 157 158 159 161 162 163] 
%Eevee rec3 [61 66 91 92 96 97 106 107 108 109 110] 

%15407 DBS shit_cells = [1 39 43] ;
%14481 DBS shit_cells = [120 121 122 123] ;
%20333 shit_cells = [61 62 63 64 65 66 67] ;



%% Save all good normalized cells in figure

saveim = true; %save images

C_raw_norm_cut(shit_cells,:) = [];

figure
plot((C_raw_norm_cut + (1: size(C_raw_norm_cut, 1))')')
axis tight
axis square
title('Sample identified traces of good cells (C raw)')

if saveim
    saveas(gcf, 'C_raw_goodCells', 'bmp');
end

%% fitted traces

% Now normalize fitted traces (C) and cut out the shit cells
C_norm = results.C;

for i = 1: size(C_norm, 1)
    C_norm(i, :) = normalize(C_norm(i, :));
end

C_norm_cut = C_norm;
C_norm_cut(isnan(C_norm_cut))=0;
C_norm_cut(shit_cells,:) = [];


figure
plot((C_norm_cut + (1: size(C_norm_cut, 1))')')
axis tight
axis square
title('Fitted traces of good cells (C)')

if saveim
    saveas(gcf, 'C_goodCells', 'bmp');
end

%% sort cells based on max spiking

saveim = true; %save images

C_raw_norm_cut_max = C_raw_norm_cut;

%  sort traces
[M,I] = max(C_raw_norm_cut_max, [], 2);
[~,idx] = sort(I, 'descend');
C_raw_norm_cut_max = C_raw_norm_cut_max(idx,:);

figure
plot((C_raw_norm_cut_max(1:size(C_raw_norm_cut_max),:) + (1:size(C_raw_norm_cut_max))')')
axis tight
axis square
title('Identified traces sorted based on max signal')

if saveim
    saveas(gcf, 'C_raw_goodCells_sorted', 'bmp');
end



%% do the same for spatial (A) footprint

% take out shit cells
A_cut = A;
A_cut(:, shit_cells) = []; 


% to fix

% figure
% %spy(A(:,1), [222 358]);
% %gplot(A(:,1),[222 358]);
% imagesc(sparse(A(:,100)))
% colorbar
% 
% imagesc(results.Cn)

%% do the same for spikes (S)

S_cut = S;
S_cut(shit_cells, :) = []; 
S_norm_cut = S_cut;

for i = 1: size(S_norm_cut, 1)
    S_norm_cut(i, :) = normalize(S_norm_cut(i, :));
end

figure
plot((S_norm_cut(1:size(S_norm_cut),:) + (1:size(S_norm_cut))')')
axis tight
axis square
title('Spiking of good cells (S)')

if saveim
    saveas(gcf, 'S_goodCells', 'bmp');
end


%% save into new results file

results_cut = results;
results_cut.A = A_cut;
results_cut.C = C_norm_cut;
results_cut.C_raw = C_raw_norm_cut;
results_cut.S = S_norm_cut;

save('msCamAll_small_NormCorre_CropGray_results_cut.mat', 'results_cut');
%now save results_cut



%% plot data with window over grooming bouts

%Conny rec1
% rectangle('Position',[412 0 24 107]) %[left_down_X Y pixels_to_right
% pixels_up]
% rectangle('Position',[2940 0 2082 107])
% rectangle('Position',[8062 0 76 107])

%Eevee rec1
% rectangle('Position',[0 0 238 113])
% rectangle('Position',[696 0 16 113])

%Nessie rec1
% rectangle('Position',[8024 0 44 160])
% rectangle('Position',[8204 0 696 160])

%Conny rec4
% rectangle('Position',[778 0 324 61])

%Conny rec5
% rectangle('Position',[980 0 174 99])
% rectangle('Position',[2557 0 115 99])

%Eevee rec2
% rectangle('Position',[960 0 69 144])
% rectangle('Position',[1655 0 433 144])

%Eevee rec3
%rectangle('Position',[831 0 551 100])
%rectangle('Position',[2244 0 70 100])

%Eevee rec6
%rectangle('Position',[0 0 1300 42])
%rectangle('Position',[1564 0 236 42])


%DBS SHORT 1M
% rectangle('Position',[900 0 900 43])
% rectangle('Position',[3600 0 900 43])
% 
% rectangle('Position',[900 0 900 120])
% rectangle('Position',[3600 0 900 120])
% 
% rectangle('Position',[900 0 900 61])
% rectangle('Position',[3600 0 900 61])


A = results_cut.A;
C = results_cut.C;
C_raw = results_cut.C_raw;
S = results_cut.S;

height = 102;
start1 = 900;
length1 = 900;
start2 = 3500;
length2 = 900;


saveim = true;

%plot figure C_raw
figure
plot((C_raw + (1: size(C_raw, 1))')')
rectangle('Position',[start1 0 length1 height])
rectangle('Position',[start2 0 length2 height])
axis tight
axis square
title('Identified traces of good cells (C raw)')


if saveim
    saveas(gcf, 'C_raw_goodCells_groom', 'bmp');
end



%plot figure c sorted based on max freq
C_raw_norm_cut_max = C_raw;

%  sort traces
[M,I] = max(C_raw_norm_cut_max, [], 2);
[~,idx] = sort(I, 'descend');
C_raw_norm_cut_max = C_raw_norm_cut_max(idx,:);

figure
plot((C_raw_norm_cut_max(1:size(C_raw_norm_cut_max),:) + (1:size(C_raw_norm_cut_max))')')
rectangle('Position',[start1 0 length1 height])
rectangle('Position',[start2 0 length2 height])
axis tight
axis square
title('Identified traces sorted based on max signal')

if saveim
    saveas(gcf, 'C_raw_goodCells_sorted_groom', 'bmp');
end


%plot figure C
figure
plot((C + (1: size(C, 1))')')
rectangle('Position',[start1 0 length1 height])
rectangle('Position',[start2 0 length2 height])
axis tight
axis square
title('Fitted traces of good cells (C)')

if saveim
    saveas(gcf, 'C_goodCells_groom', 'bmp');
end


%plot figure with spikes
figure
plot((S(1:size(S),:) + (1:size(S))')')
rectangle('Position',[start1 0 length1 height])
rectangle('Position',[start2 0 length2 height])
axis tight
axis square
title('Identified traces sorted based on max signal')

if saveim
    saveas(gcf, 'S_goodCells_groomNoNorm', 'bmp');
end


%% cut out those grooming events

% for now, just cut all grooming events and put next to each other
% cut 10sec pre grooming as baseline/non-grooming
% for convenience, start with spikes (and calculate herz)


%Conny rec1 groom 2940:5022, nongroom 857:2939
%Eevee rec1 groom 1:238, nongroom 239:476
%Nessie rec1 groom 8204:8900, RIGHT BEFORE IS ALSO GROOM, SO PICK BEGIN SESSION
%nongroom 1:697

%Conny rec4 groom 778:1077, nongroom 478:777

%conny rec5.1 groom 980:1129, nongroom 830:979
%conny rec5.2 groom 2557:2646, nongroom 2467:2556

%Eevee rec2.1 groom 960:1019, nongroom 900:959
%Eevee rec2.2 groom 1655:2088, nongroom 1222:1654

%Eevee rec3.1 groom 831:1370, nongroom 291:830
%Eevee rec3.2 groom 2244:2303, nongroom 2184:2243

%Eevee rec6.1 groom 1:1290, nongroom 2711:4000
%Eevee rec6.2 groom 1564:1773, nongroom 1354:1563

%!!!!! - ADJUST - !!!!!
%grooming range
gr = [1:450];
%nongrooming range
ngr = [451:719];

S = results_cut.S;
S_groom = full(S);

%!!!!! - ADJUST - !!!!!
%shitcells = (142:143);     %fill in cells with NaN values
%S_groom(shitcells, :) = [];       %cells have no value, get rid of them

S_groom_temp = S_groom(:, gr); %select grooming frames  1:238

S_groom_ave = mean(S_groom_temp, 'omitnan');     %sum grooming over time
S_groom_SEM = std(S_groom_temp, 'omitnan')/sqrt(size(S_groom_temp,1)); %SEM over time

%!!!!! - ADJUST - !!!!!
%bin in 30frames (1sec)
%to fix! now you have to manually add NaN's at the end to get to round
%number
%S_groom_ave(326:330) = NaN;      % divide by 30 using [] or NaN    
%S_groom_SEM(326:330) = NaN;

S_groom_ave_bin = arrayfun(@(X) ( nanmean(S_groom_ave([X:(X+29)]))), [1:30:size(S_groom_ave,2)-29] );
S_groom_SEM_bin = arrayfun(@(X) ( nanmean(S_groom_SEM([X:(X+29)]))), [1:30:size(S_groom_SEM,2)-29] );
S_groom_time = [S_groom_ave_bin; S_groom_SEM_bin];      %variable grooming per sec

S_groom_ave_hz = nansum(S_groom_ave_bin) / size(S_groom_ave_bin,2) %calculate spike Hz  



%not grooming
S_notGroom = full(S);

%!!!!! - ADJUST - !!!!!
%S_notGroom(shitcells, :) = [];       %cells have no value, get rid of them

S_notGroom_temp = S_notGroom(:, ngr); %select non-sgrooming (10s pre = 300frames)

S_notGroom_ave = mean(S_notGroom_temp, 'omitnan');     %sum non-grooming over time
S_notGroom_SEM = std(S_notGroom_temp, 'omitnan')/sqrt(size(S_notGroom_temp,1)); %SEM over time

%!!!!! - ADJUST - !!!!!
%bin in 30frames (1sec)
%to fix! now you have to manually add NaN's at the end to get to round
%number
%S_notGroom_ave(326:330) = NaN;
%S_notGroom_SEM(326:330) = NaN;

S_notGroom_ave_bin = arrayfun(@(X) ( nanmean(S_notGroom_ave([X:(X+29)]))), [1:30:size(S_notGroom_ave,2)-29] );
S_notGroom_SEM_bin = arrayfun(@(X) ( nanmean(S_notGroom_SEM([X:(X+29)]))), [1:30:size(S_notGroom_SEM,2)-29] );
S_notGroom_time = [S_notGroom_ave_bin; S_notGroom_SEM_bin];      %variable non-grooming per sec


S_notGroom_ave_hz = nansum(S_notGroom_ave_bin) / size(S_notGroom_ave_bin,2)    %calculate spike Hz


%final output
S_final_grooming = padcat(S_groom_ave_hz, S_notGroom_ave_hz, S_groom_time(1,:), S_groom_time(2,:), S_notGroom_time(1,:), S_notGroom_time(2,:));

save('S_final_grooming.mat', 'S_final_grooming')







%% cut out those grooming events and do the same thing with C_raw

% for now, just cut all grooming events and put next to each other
% cut 10sec pre grooming as baseline/non-grooming
% for convenience, start with spikes (and calculate herz)


%Conny

%grooming
C_raw = results_cut.C_raw;
C_raw_groom = full(C_raw);

%C_raw_groom(56:60, :) = [];       %cells have no value, get rid of them
C_raw_groom_temp = C_raw_groom(:, gr); %select grooming frames

C_raw_groom_ave = mean(C_raw_groom_temp, 'omitnan');     %sum grooming over time
C_raw_groom_SEM = std(C_raw_groom_temp, 'omitnan')/sqrt(size(C_raw_groom_temp,1)); %SEM over time

%!!!!! - ADJUST - !!!!!
%bin in 30frames (1sec)
%to fix! now you have to manually add NaN's at the end to get to round
%number
%C_raw_groom_ave(326:330) = NaN;
%C_raw_groom_SEM(326:330) = NaN;

C_raw_groom_ave_bin = arrayfun(@(X) ( nanmean(C_raw_groom_ave([X:(X+29)]))), [1:30:size(C_raw_groom_ave,2)-29] );
C_raw_groom_SEM_bin = arrayfun(@(X) ( nanmean(C_raw_groom_SEM([X:(X+29)]))), [1:30:size(C_raw_groom_SEM,2)-29] );
C_raw_groom_time = [C_raw_groom_ave_bin; C_raw_groom_SEM_bin];      %variable grooming per sec

C_raw_groom_ave_hz = nansum(C_raw_groom_ave_bin) / size(C_raw_groom_ave_bin,2) %calculate spike Hz  



%not grooming
C_raw_notGroom = full(C_raw);

%C_raw_notGroom(56:60, :) = [];       %cells have no value, get rid of them
C_raw_notGroom_temp = C_raw_notGroom(:, ngr); %select non-grooming (10s pre = 300frames)

C_raw_notGroom_ave = mean(C_raw_notGroom_temp, 'omitnan');     %sum non-grooming over time
C_raw_notGroom_SEM = std(C_raw_notGroom_temp, 'omitnan')/sqrt(size(C_raw_notGroom_temp,1)); %SEM over time

%!!!!! - ADJUST - !!!!!
%bin in 30frames (1sec)
%to fix! now you have to manually add NaN's at the end to get to round
%number
%C_raw_notGroom_ave(326:330) = NaN;
%C_raw_notGroom_SEM(326:330) = NaN;

C_raw_notGroom_ave_bin = arrayfun(@(X) ( nanmean(C_raw_notGroom_ave([X:(X+29)]))), [1:30:size(C_raw_notGroom_ave,2)-29] );
C_raw_notGroom_SEM_bin = arrayfun(@(X) ( nanmean(C_raw_notGroom_SEM([X:(X+29)]))), [1:30:size(C_raw_notGroom_SEM,2)-29] );
C_raw_notGroom_time = [C_raw_notGroom_ave_bin; C_raw_notGroom_SEM_bin];      %variable non-grooming per sec


C_raw_notGroom_ave_hz = nansum(C_raw_notGroom_ave_bin) / size(C_raw_notGroom_ave_bin,2)    %calculate spike Hz


%final output
C_raw_final_grooming = padcat(C_raw_groom_ave_hz, C_raw_notGroom_ave_hz, C_raw_groom_time(1,:), C_raw_groom_time(2,:), C_raw_notGroom_time(1,:), C_raw_notGroom_time(2,:));

save('C_raw_final_grooming.mat', 'C_raw_final_grooming')









