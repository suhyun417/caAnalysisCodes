% genFig_registerSessions.m
%
% 2022/10 SHP
% - visualization of longitudinal registration results
% - same neuron's responses across trials over days

clear all;


%% settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/NIFVAULT/PROJECTS/parksh';
        dirProcdata = '/Volumes/NIFVAULT/PROCDATA/parksh';
        dirRawdata = '/Volumes/rawdata/parksh';
    else % on virtual machine
        dirProjects = '/nifvault/projects/parksh';
        dirProcdata = '/nifvault/procdata/parksh';
        dirRawdata = '/nifvault/rawdata/parksh';
    end
end

addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));
% gcp; % for parallel processingls

dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');


%% Session info & optional parameters
setSubj = {'Tabla', 1; 'Max', 3};

iSubj = 2; % 1;

nameSubj = setSubj{iSubj,1}; %'Max'; % 'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
FOV_ID = setSubj{iSubj,2}; %3; %1; %3; %1;
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);


%% load saved files
% cells pooled across days
fname_stack = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellAcrossDay.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_stack, 'cellIDAcrossDay'); %, 'stackCellCenter')

% cell quality info 
fname_cellQC = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellQC.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_cellQC, 'infoCells')

% translational shift across days
fname_shifts = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_shifts.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID));  
load(fname_shifts, 'shifts')

% aligned cells TS and spatial info
fname_caTSFOV = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_DFLsorted.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
load(fname_caTSFOV, 'cellTS', 'cellPix')

%% Stacked cells across days: playing with examples
indCellValid = find(cat(1, cellTS.nTrial1_total)>8);

% tempATrial = cat(2, cellPix(indCellValid).repPix);
% tempATrial(~isnan(tempATrial)) = 10;
tempA = cat(2, cellPix.repPix);
tempA(~isnan(tempA)) = 1;
tempA(:, indCellValid) = tempA(:, indCellValid).*10;

imgCells = sum(tempA, 2, 'omitnan');
imgCells_2d = reshape(imgCells, size(infoCells(1).imgFOV));

figure;
set(gcf, 'Color', 'w')
imagesc(imgCells_2d)
colormap(turbo)

% A_temp = full(reshape(neuron.A(:,i),d1,d2));

figure;
set(gcf, 'Color', 'w')
cmap_cell = colormap(hsv(length(cellPix)));
for iCell = 1:length(cellPix)
    plot(cellPix(iCell).contourCell{1}(1,1:end), cellPix(iCell).contourCell{1}(2,1:end), ...
        'Color', cmap_cell(iCell, :), 'linewidth', 1); hold on;
    text(cellPix(iCell).contourCell{1}(1,end), cellPix(iCell).contourCell{1}(2,end), num2str(iCell), ...
        'color', 'k')
end
set(gca, 'YDir', 'reverse', 'XLim', [0-20 size(infoCells(1).imgFOV, 2)+20], 'YLim', [0-20 size(infoCells(1).imgFOV, 1)+20])



%%
isCell = ~isnan(cellIDAcrossDay);
temp = sum(isCell, 2);
indCell_long = find(temp>4);


%% plot
fig_supercellmovie = figure;
set(fig_supercellmovie, 'Position', [1500 1000 1200 400], 'Color', 'w')
sp(1) = subplot('Position', [0.05 0.25 0.3 0.6]);
sp(2) = subplot('Position', [0.4 0.55 0.55 0.35]);
sp(3) = subplot('Position', [0.4 0.05 0.55 0.35]);

for iC = 1:length(indCellValid)
    
    iCell = indCellValid(iC); %70;
    
    figure(fig_supercellmovie);
    axes(sp(1)); cla;
    A_temp = cellPix(iCell).repPix./max(cellPix(iCell).repPix);
    contourf(reshape(A_temp,d1,d2), [0,0]+1); hold on;
    text(cellPix(iCell).contourCell{1}(1,end), cellPix(iCell).contourCell{1}(2,end),...
        num2str(iCell), 'Color', 'k')
        
    
    axes(sp(2));
    imagesc(zscore(cellTS(iCell).matTS_movie1, 0, 2));
    %     set(gca, 'XTickLabel', 20:20:120, 'YTick', cellTS(iCell).nTrial1_set , 'YTickLabel', setDateSession(cellTS(iCell).idAcrossSession(:,1)),'TickDir', 'out', 'Box', 'off')
    title(sprintf('Cell %d/%d: Mov 1', iCell, length(cellTS)))
    %     colormap(hot)
    
    axes(sp(3));
    imagesc(zscore(cellTS(iCell).matTS_movie2, 0, 2));
    %     set(gca, 'XTickLabel', 20:20:120, 'YTick', cellTS(iCell).nTrial2_set, 'YTickLabel', setDateSession(cellTS(iCell).idAcrossSession(:,1)),'TickDir', 'out', 'Box', 'off')
    title('Mov 2')
    xlabel('Time (s)')
    %     colormap(hot)
    
    set(sp(1), 'YDir', 'reverse', 'XLim', [0-20 size(infoCells(1).imgFOV, 2)+20], 'YLim', [0-20 size(infoCells(1).imgFOV, 1)+20])
    colormap(sp(1), 'gray')
    set(sp(2:3), 'TickDir', 'out', 'XTickLabel', 20:20:120)
    colormap(sp(2), 'hot')
    colormap(sp(3), 'hot')


    input('')
end

% isCell = ~isnan(stackCellCenter);
% catCell = sum(isCell, 3);
% figure;
% imagesc(catCell)
% % catCell(8:14, 102:107)
% % stackCellCenter(8:14, 102:107, :)

% % prep the DFL ts from all the sessions
% for iS = 1:length(setDateSession)
%     dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
%     dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
%     load(fullfile(dirProcdata_session, 'DFL_ts_tML'));
%     
%     resultsDFL(iS).tS_session = tS_session;
% end
% clear tS_session 
% 
% fig_supercellmovie = figure;
% set(fig_supercellmovie, 'Position', [1500 1000 800 600], 'Color', 'w')
% 
% for iC = 1:length(indCell_long)
%     
%     cellID = cellIDAcrossDay(indCell_long(iC),:);
%     
%     
%     % cellID = [86 75 96 69 55 61 63 nan]; %[6 10 13 10 12 8 6 nan]; %[109 93 113 nan nan 76 74 nan]; %[91 81 104 75 nan 64 nan nan]; %[69 53 68 nan 35 nan nan nan]; %[21 12 21 18 nan 12 12 16]; %[48 33 45 34 nan 32 33 32]; %[34 20 32 nan nan 28 23 25]; %[118 nan 28 nan nan 17 17 nan]; %[83 68 92 61 48 55 58 50]; %from Max FOV3
%     % nameSubj = 'Max'; %'Tabla';
%     % FOV_ID = 3; %1;
%     % [infoSession, opts] = readInfoSession(nameSubj, FOV_ID);
%     %
%     % [c, ia, indRun] = unique(infoSession.(1), 'sorted');
%     % setDateSession = c(2:end); % 1st one is always empty
%     % nSession = length(setDateSession);
%     
%     tempMatTS1 = []; tempMatTS2 = []; nTrial = [];
%     for iS = 1:length(cellID)
%         
%         if isnan(cellID(iS))
%             continue;
%         end
%         
%         %     dateSession = setDateSession{iS}; %'20191113'; % '20191125'; %'20191113'; %'20191125';
%         %     % datestr(datenum(dateSession, 'yyyymmdd'), 'yymmdd') % for bhv files
%         %
%         %     dirProcdata_session = fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
%         %
%         %     load(fullfile(dirProcdata_session, 'DFL_ts_tML'));
%         
%         %     figure(200);
%         %     subplot(2,1,1)
%         %     title('Mov 1')
%         %     plot(squeeze(tS_session(1).matTS_norm(:, cellID(iS), :)))
%         %     hold on
%         
%         tempMatTS1 = cat(1, tempMatTS1, squeeze(resultsDFL(iS).tS_session(1).matTS_C_raw_zscore(:, cellID(iS), :))');
%         nTrial = cat(1, nTrial, size(tempMatTS1, 1));
%         
%         %     subplot(2,1,2)
%         %     title('Mov 2')
%         %     plot(squeeze(tS_session(2).matTS_norm(:, cellID(iS), :)))
%         %     hold on
%         
%         tempMatTS2 = cat(1, tempMatTS2, squeeze(resultsDFL(iS).tS_session(2).matTS_C_raw_zscore(:, cellID(iS), :))');
%     end
%     
%     tempValidS = ~isnan(cellID);
%     
%     figure(fig_supercellmovie);
%     set(gcf, 'color', 'w')
%     subplot(2,1,1)
%     imagesc(tempMatTS1);
%     set(gca, 'XTickLabel', 20:20:120, 'YTick', nTrial, 'YTickLabel', setDateSession(tempValidS),'TickDir', 'out', 'Box', 'off')
% %     title(sprintf('Cell %d/%d: Mov 1', iC, length(indCell_long)))
%     colormap(hot)
%     subplot(2,1,2)
%     imagesc(tempMatTS2);
%     set(gca, 'XTickLabel', 20:20:120, 'YTick', nTrial, 'YTickLabel', setDateSession(tempValidS),'TickDir', 'out', 'Box', 'off')
%     title('Mov 2')
%     xlabel('Time (s)')
%     colormap(hot)
%     
% %     figure(200);
% %     set(gcf, 'color', 'w')
% %     subplot(2,1,1)
% %     plot(tempMatTS1'); axis tight
% %     set(gca, 'XTickLabel', 20:20:120, 'TickDir', 'out', 'Box', 'off')
% %     title(sprintf('Cell %d/%d: Mov 1', iC, length(indCell_long)))
% % %     colormap(hot)
% %     subplot(2,1,2)
% %     plot(tempMatTS2'); axis tight
% %     set(gca, 'XTickLabel', 20:20:120, 'TickDir', 'out', 'Box', 'off')
% %     title('Mov 2')
% %     xlabel('Time (s)')
% % %     colormap(hot)
%     
%     input('')
%     
% end


% %% Visualization of cell contours across days
% cellColor = hsv(nSession);
% % cellColor = [166,206,227;...
% % 31,120,180;...
% % 178,223,138;...
% % 51,160,44;...
% % 251,154,153;...
% % 227,26,28;...
% % 253,191,111;...
% % 255,127,0;...
% % 202,178,214;...
% % 106,61,154;...
% % 255,255,153;...
% % 177,89,40]./255; 
% fig_acSession = figure;
% set(fig_acSession, 'Color', 'w')
% 
% for iSession = 1:nSession
% % iSession = 1; 
% dateSession = setDateSession{iSession}; %'20191113'; %setDateSession{iSession};
% 
% dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
% dirPreproc = fullfile(dirProcdata_session, '_preproc');
% 
% 
% %%
% addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
% cnmfe_setup;
% d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
% 
% load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
% 
% validIndCell = [];
% validIndCell(:,1) = 1:length(neuron.ids);
% if strcmpi(nameSubj, 'max')
%     load(fullfile(dirProcdata_session, 'validIndCell.mat'), 'indCell')
%     validIndCell = indCell.validCell;
% end
% 
% % % % sort the cells based on PNR
% % % [a, sortedIndCell] = sort(infoCells(iSession).pnrs, 'descend');
% % [sortedIndCell] = orderROIs(neuron, 'circularity');
% 
% % color cells to indicate sessions
% 
% thr = 0.3;
% widthContour = 1.5;
% Coor = neuron.get_contours(thr); 
% infoCells(iSession).coor_0p3 = Coor;
% 
% figure(fig_acSession);
% for i = 1:size(Coor, 1)
%     %         cont = medfilt1(Coor{i}')';
%     cont = Coor{i};
%     if size(cont,2) > 1 % "shifts" values are in the order of image matrix dimensions: first dimension is vertical ("y") and second dimension is horizontal ("x")
%         plot(cont(1,1:end)+shifts(iSession, 2), cont(2,1:end)+shifts(iSession, 1), 'Color', cellColor(iSession, :), 'linewidth', widthContour); hold on;
%     end
% end
% end
% axis tight
% 
% for thr = [0.2:0.1:0.5]
%     locbad=[]; nrows=[]; ncols=[];
% Coor = neuron.get_contours(thr); 
% [nrows, ncols] = cellfun(@size, Coor);
% locbad = find(ncols<2); % the ones that get_contours couldn't get a reasonable localized contour at this threshold
% length(locbad)
% end
% 
% figure;
% for iC = 1:length(locbad)    
%     i = locbad(iC);
%     A_temp = full(reshape(neuron.A(:,i),d1,d2));
%     A_temp = medfilt2(A_temp,[3,3]);
%     A_temp = A_temp(:);
%     [temp,ind] = sort(A_temp(:).^2,'ascend');
%     temp =  cumsum(temp);
%     ff = find(temp > (1-thr)*temp(end),1,'first');
%     if ~isempty(ff)
%         CC{i} = contour(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor', 'b'); hold on;
%     end
%     set(gca, 'YDir', 'reverse', 'XLim', [0 d2], 'YLim', [0 d1])
%     input('')
% end

%% order_ROIs
function [srt] = orderROIs(obj, srt)
%% order neurons
% srt: sorting order
nA = sqrt(sum(obj.A.^2));
nr = length(nA);
if nargin<2
    srt='srt';
end
K = size(obj.C, 1);

if ischar(srt)
    if strcmpi(srt, 'decay_time')
        % time constant
        if K<0
            disp('Are you kidding? You extracted 0 neurons!');
            return;
        else
            taud = zeros(K, 1);
            for m=1:K
                temp = ar2exp(obj.P.kernel_pars(m));
                taud(m) = temp(1);
            end
            [~, srt] = sort(taud);
        end
    elseif strcmp(srt, 'mean')
        if obj.options.deconv_flag
            temp = mean(obj.C,2)'.*sum(obj.A);
        else
            temp = mean(obj.C,2)'.*sum(obj.A)./obj.P.neuron.sn';
        end
        [~, srt] = sort(temp, 'descend');
    elseif strcmp(srt, 'sparsity_spatial')
        temp = sqrt(sum(obj.A.^2, 1))./sum(abs(obj.A), 1);
        [~, srt] = sort(temp);
    elseif strcmp(srt, 'sparsity_temporal')
        temp = sqrt(sum(obj.C_raw.^2, 2))./sum(abs(obj.C_raw), 2);
        [~, srt] = sort(temp, 'descend');
    elseif strcmp(srt, 'circularity')
        % order neurons based on its circularity
        tmp_circularity = zeros(K,1);
        for m=1:K
            [w, r] = nnmf(obj.reshape(obj.A(:, m),2), 1);
            ky = sum(w>max(w)*0.3);
            kx = sum(r>max(r)*0.3);
            tmp_circularity(m) = abs((kx-ky+0.5)/((kx+ky)^2));
        end
        [~, srt] = sort(tmp_circularity, 'ascend');
    elseif strcmpi(srt, 'pnr')
        pnrs = max(obj.C, [], 2)./std(obj.C_raw-obj.C, 0, 2);
        [~, srt] = sort(pnrs, 'descend');
    elseif strcmpi(srt, 'temporal_cluster')
        obj.orderROIs('pnr');
        dd = pdist(obj.C_raw, 'cosine');
        tree = linkage(dd, 'complete');
        srt = optimalleaforder(tree, dd);
    elseif strcmpi(srt, 'spatial_cluster')
        obj.orderROIs('pnr');
        A_ = bsxfun(@times, obj.A, 1./sqrt(sum(obj.A.^2, 1)));
        temp = 1-A_' * A_;
        dd = temp(tril(true(size(temp)), -1));
        dd = reshape(dd, 1, []);
        tree = linkage(dd, 'complete');
        srt = optimalleaforder(tree, dd);
    else %if strcmpi(srt, 'snr')
        snrs = var(obj.C, 0, 2)./var(obj.C_raw-obj.C, 0, 2);
        [~, srt] = sort(snrs, 'descend');
    end
end
obj.A = obj.A(:, srt);
obj.C = obj.C(srt, :);

%
% figure;


thr = 0.2;
CC = cell(size(neuron.A, 2),1);
CR = cell(size(neuron.A, 2),2);
% cmap_cell = cool(size(neuron.A, 2));
for iCC = 1:size(neuron.A ,2)
    i = iCC; %sortedIndCell(iCC);
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contour(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor',cellColor(iSession, :), 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off


% title(sprintf('%s: %s', nameSubj, dateSession))





% get the contours and image field of view
% neuron_b = neuron.batches{1}.neuron;

% Generate cell location map within FOV
thr = 0.5; % the lower the smaller (more centralized) the contour
cellColor = [1 1 1];
widthContour = 1;
[d1,d2] = size(neuron.Cn);

figure;
imagesc(zeros(d1, d2)); % background
colormap(gray);
caxis([0 0.1]);
hold on;

CC = cell(size(neuron.A, 2),1);
CR = cell(size(neuron.A, 2),2);
for i = 1:size(neuron.A ,2)
    A_temp = full(reshape(neuron.A(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        CC{i} = contourf(reshape(A_temp,d1,d2), [0,0]+A_temp(ind(ff)), 'LineColor',cellColor, 'linewidth', widthContour);
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{i,1} = [ii,jj]';
        CR{i,2} = A_temp(fp)';
    end
    hold on;
end
axis off
title(sprintf('%s: %s', nameSubj, dateSession))
% %save
% print(gcf, fullfile(dirFig, sprintf('SourceFOV_solidWhite_bkgdBlack_thr%s_%s_%s', strrep(num2str(thr),'.', 'p'), nameSubj, dateSession)), '-depsc');

% end
%

% thr = 0.3; % the lower the smaller (more centralized) the contour
Coor = neuron.get_contours(thr);
imgFOV = neuron.Cn.*neuron.PNR;

%         figure;
% neuron.show_contours([], [], imgFOV, 'true');

figure;
[center] = neuron.estCenter();
imagesc(imgFOV); colormap(gray);
hold on
plot(center(:,2), center(:, 1), 'r.');
text(center(:,2)+1, center(:,1), num2str([1:length(neuron.ids)]'), 'Color', 'w');




%% Possible code
% 1. Function performing longitudinal registration
% - Load reference image
% - For each following sessions
%   - load session image
%   - perform the rigid-body registration
%   - save the shifts (and registration parameters and other outcomes)
%       - for each animal and each FOV (e.g. Tabla_FOV1_shifts.mat)
% 2. Another script to extract cells across daily sessions using the shifts
% - Maybe it's better to pool in the pixel space. imagine a grid in n x n pixel
% resolution (e.g. 3 pixels? I can choose the criterion) and make a list of
% grid in columnar way. Then gather cell IDs from the REF and SES 
% 3. Some kind of cell sorting/exclusion need to be done to exclude
% non-valid sources (e.g. not circular, too low SNR/PNR). 
%   Example (from Sources2D.m):
%       nA = sqrt(sum(obj.A.^2));
%       nr = length(nA);
%       K = size(obj.C, 1);
%       % order neurons based on its circularity
%       tmp_circularity = zeros(K,1);
%       for m=1:K
%           [w, r] = nnmf(obj.reshape(obj.A(:, m),2), 1);
%           ky = sum(w>max(w)*0.3);
%           kx = sum(r>max(r)*0.3);
%           tmp_circularity(m) = abs((kx-ky+0.5)/((kx+ky)^2));
%       end
%       [~, srt] = sort(tmp_circularity, 'ascend');
% 0. The cell exclusion can be done first before pooling cells across days.
% 0. So from now on we will always load the valid cell index computed from
% the above code #3 and apply that to the time series etc. Dealing with
% indices will be tricky and requires extra attention.

% ## Other things
% %% trim spatial components
%         function [ind_small] = trimSpatial(obj, thr, sz)
%             % remove small nonzero pixels
%             if nargin<2;    thr = 0.01; end
%             if nargin<3;    sz = 5; end
%             
%             se = strel('square', sz);
%             ind_small = false(size(obj.A, 2), 1);
%             for m=1:size(obj.A,2)
%                 ai = obj.A(:,m);
%                 ai_open = imopen(obj.reshape(ai,2), se);
%                 
%                 temp = full(ai_open>max(ai)*thr);
%                 l = bwlabel(obj.reshape(temp,2), 4);   % remove disconnected components
%                 [~, ind_max] = max(ai_open(:));
%                 
%                 ai(l(:)~=l(ind_max)) = 0;
%                 if sum(ai(:)>0) < obj.options.min_pixel %the ROI is too small
%                     ind_small(m) = true;
%                 end
%                 obj.A(:, m) = ai(:);
%             end
%             %             ind_small = find(ind_small);
%             %             obj.delete(ind_small);
%         end
% %% keep spatial shapes compact
%         function compactSpatial(obj)
%             for m=1:size(obj.A, 2)
%                 ai = obj.reshape(obj.A(:, m), 2);
%                 ai = circular_constraints(ai);
%                 obj.A(:, m) = ai(:);
%             end
%         end




