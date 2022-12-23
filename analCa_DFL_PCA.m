% analCa_DFL_PCA.m
%
% 2022/12/23 SHP modified from genFig_dataMining_DFL.m
% perform PCA on dynamic face localizer (movie) responses and save the
% results of PCA

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

for iSubj = 1:2 %2; %1;
    
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
    load(fname_caTSFOV, 'cellTS', 'cellPix', 'resultsDFL')
    
    %%
    % block information
    paramPCA.condName = {'Scene Motion', 'Human Faces', 'Monkey Bodies', 'Monkey Faces', 'Object Motion', 'Optic Flow'}';
    paramPCA.condOrder = [1 2 6 3 5 4; 5 6 4 1 3 2]';
    [b, indCondReorder] = sort(paramPCA.condOrder, 1);
    
    indCellValid = find(cat(1, cellTS.nTrial1_total)>8); %
    % % indCellValid_trial = find(cat(1, cellTS.nTrial1_total)>8); % cells that have more than 8 trials for movie 1
    % %
    % for ii = 1:length(cellTS)
    %     minsnrmovie1(ii, 1) = min(cellTS(ii).snr_movie1);
    % end
    %
    % % indCellValid_snr = find(minsnrmovie1>0.1);
    %
    % % indCellValid = intersect(indCellValid_trial, indCellValid_snr);
    
    clear matAvg* resultsPCA_*
    for iCell = 1:length(indCellValid)
        
        matAvgTS1(:, iCell) = mean(cellTS(indCellValid(iCell)).matTS_movie1)'; % now it's scaled dF
        matAvgTS2(:, iCell) = mean(cellTS(indCellValid(iCell)).matTS_movie2)'; %
        
        steAvgTS1(:, iCell) = std((cellTS(indCellValid(iCell)).matTS_movie1)./sqrt(size(cellTS(indCellValid(iCell)).matTS_movie1, 1)-1))'; % now it's scaled dF
        steAvgTS2(:, iCell) = std((cellTS(indCellValid(iCell)).matTS_movie2)./sqrt(size(cellTS(indCellValid(iCell)).matTS_movie2, 1)-1))';
        
    end
    matAvgTS(1).matTS = matAvgTS1;
    matAvgTS(1).steTS = steAvgTS1;
    matAvgTS(2).matTS = matAvgTS2;
    matAvgTS(2).steTS = steAvgTS2;
    
    for iMovie = 1:2
        
        matAvgTS(iMovie).matAvgTS_block = reshape(matAvgTS(iMovie).matTS(1:1200, :)', size(matAvgTS(iMovie).matTS, 2), 200, 6);
        matAvgTS(iMovie).matAvgTS_block_reorder = matAvgTS(iMovie).matAvgTS_block(:,:,indCondReorder(:,iMovie));
        
        % matAvgTS1_block = reshape(matAvgTS1(1:1200, :)', size(matAvgTS1, 2), 200, 6);
        % matAvgTS2_block = reshape(matAvgTS2(1:1200, :)', size(matAvgTS2, 2), 200, 6);
        %
        % matAvgTS_block_reorder{1} = matAvgTS1_block(:,:,indCondReorder(:,1));
        % matAvgTS_block_reorder{2} = matAvgTS2_block(:,:,indCondReorder(:,2));
        
        
        [coeffRun, scoreRun, latentRun, tsquaredRun, explainedRun] = pca(zscore(matAvgTS(iMovie).matTS)');
        [sortedScore, indCellRun] = sort(scoreRun(:,1:10), 'descend');
        
        resultsPCA_run(iMovie).explained = explainedRun;
        resultsPCA_run(iMovie).coeff = coeffRun(:, 1:10);
        resultsPCA_run(iMovie).score = scoreRun(:, 1:10);
        resultsPCA_run(iMovie).indCellSorted = indCellRun;
        
        for iB = 1:6
            
            [coeff, score, latent, tsquared, explained] = pca(zscore(matAvgTS(iMovie).matAvgTS_block_reorder(:,:,iB)')'); %pca(zscore(matAvgTS1)');
            [sortedScore, indCell] = sort(score(:,1:10), 'descend');
            
            %     explained(1:10)
            
            %             figure;
            %             set(gcf, 'Color', 'w', 'Position', [700 700 885 415])
            %             subplot(1,3,1);
            %             plot(coeff(:,1:3))
            %             legend('PC1', 'PC2', 'PC3')
            %             set(gca, 'XTick', 0:100:200, 'XTickLabel', 0:10:20)
            %             xlabel('Time (s)')
            %             title(sprintf('%s', paramCond.condName{iB}))
            %             subplot(1,3,2)
            %             imagesc(zscore(matAvgTS_block_reorder{iMovie}(indCell(:,1),:,iB)')')
            %             colormap(hot)
            %             set(gca, 'XTick', 0:50:200, 'XTicklabel', 0:5:20)
            %             title('PC1 sorted')
            %             xlabel('Time (s)')
            %             ylabel('Cells (sorted)')
            %             subplot(1,3,3)
            %             imagesc(zscore(matAvgTS_block_reorder{iMovie}(indCell(:,2),:,iB)')')
            %             colormap(hot)
            %             set(gca, 'XTick', 0:50:200, 'XTicklabel', 0:5:20)
            %             title('PC2 sorted')
            %             xlabel('Time (s)')
            %             ylabel('Cells (sorted)')
            %
            %             print(gcf, fullfile(dirFig, sprintf('%s_FOV%d_DFLmovie%d_PCA_BlockID%d', nameSubj, FOV_ID, iMovie, iB)), '-depsc')
            
            resultsPCA_block(iB, iMovie).explained = explained;
            resultsPCA_block(iB, iMovie).coeff = coeff(:, 1:10);
            resultsPCA_block(iB, iMovie).score = score(:, 1:10);
            resultsPCA_block(iB, iMovie).indCellSorted = indCell;
            
        end
    end
    
    resultsPCA(iSubj).resultsPCA_run = resultsPCA_run;
    resultsPCA(iSubj).resultsPCA_block = resultsPCA_block;
    
    resultsPCA(iSubj).matAvgTS = matAvgTS;
    
end


save(fullfile(dirProcdata, '_marmoset/invivoCalciumImaging/', 'DFL_TS_PCA.mat'), 'resultsPCA', 'paramPCA')



