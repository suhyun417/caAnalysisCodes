% getValidIndCell_Max.m
%
% For Max's cells where the FOV included inner cannula and edge of the lens,
% this script saves the indices of cells that are valid for each session
% This needs to be done just one time.
%   - 1. Load the ROIs drawn manually in ImageJ (using ReadImageJROI function
%   from "https://github.com/DylanMuir/ReadImageJROI")
%   - 2. Select the cells based on the center position of the cells
%   - 3. Save the "validIndCells.mat" to the main procdata_session directory

clear all;

addpath('/projects/parksh/_toolbox/ReadImageJROI/');
addpath('/projects/parksh/_toolbox/NoRMCorre/'); % to save the figure, I use "loadtiff" included in this package

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS/parksh';
    dirProcdata = '/Volumes/PROCDATA/parksh';
    dirRawdata = '/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/projects/parksh';
    dirProcdata = '/procdata/parksh';
    dirRawdata = '/rawdata/parksh';
end


nameSubj = 'Max';
% get session info
[infoSession, opts] = readInfoSession(nameSubj);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

for iSession = 1:nSession
    dateSession = setDateSession{iSession};
    
    dirProcdata_session = fullfile(dirProcdata, '/_marmoset/invivoCalciumImaging/', nameSubj, 'Session', dateSession);
    dirPreproc = fullfile(dirProcdata_session, '_preproc');
    
    %% Read source data and compute center location
    addpath(fullfile(dirProjects, '/_toolbox/CNMF_E/'));
    cnmfe_setup;
    d_sources2D = dir(fullfile(dirProcdata_session, 'Sources2D_all*'));
    
    load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
    
    [center] = neuron.estCenter();
    
    %% Read ROI data
    d_roi = dir(fullfile(dirProcdata_session, '*.zip'));
    [excAreaROIs] = ReadImageJROI(fullfile(d_roi.folder, d_roi.name));
    
    %% Get the cell indices
    validIndCell = [];
    for iRoi = 1:length(excAreaROIs)
        areaExc = excAreaROIs{iRoi}.vnRectBounds; % [nTop nLeft nBottom nRight] area needs to be excluded
        
        indCellExc_x1 = center(:,2)>areaExc(2);
        indCellExc_x2 = center(:,2)<areaExc(4);
        indCellExc_y1 = center(:,1)>areaExc(1);
        indCellExc_y2 = center(:,1)<areaExc(3);
        validIndCell(:,iRoi) = sum(cat(2, indCellExc_x1, indCellExc_x2, indCellExc_y1, indCellExc_y2), 2)<4;
        
    end
    
    indCell = struct([]);
    indCell(1).validCell = find(min(validIndCell, [], 2)>0);
    indCell(1).excCell = find(min(validIndCell, [], 2)<1);
    
    save(fullfile(dirProcdata_session, 'validIndCell.mat'), 'indCell', 'excAreaROIs');
    
    % save the figure with locations of included & excluded cells
    d_file = dir(fullfile(dirPreproc, '*_mc.tif'));
    imgFOV = loadtiff(fullfile(d_file(1).folder, d_file(1).name), 10, 1); % load just one frame (10th frame)
    
    figCenter1 = figure;
    set(figCenter1, 'Color', 'w', 'Position', [100 100 430 620])
    imagesc(imgFOV); colormap(gray); hold on;
    plot(center(indCell.validCell, 2), center(indCell.validCell, 1), 'b.', 'MarkerSize', 10)
    hold on;
    plot(center(indCell.excCell, 2), center(indCell.excCell, 1), 'r.', 'MarkerSize', 10)
    setCell = cat(1, indCell.validCell, indCell.excCell);
    text(center(setCell, 2)+1, center(setCell, 1), num2str(setCell), 'Color', 'k')
    axis off
    title(sprintf('%s %s: valid cells in blue', nameSubj, dateSession))
    print(figCenter1, sprintf(fullfile(dirProcdata_session, '%s_%s_validCells_1'), nameSubj, dateSession), '-depsc')
    
    figCenter2 = figure;
    set(figCenter2, 'Color', 'w', 'Position', [100 100 430 620])
    plot(center(indCell.validCell, 2), center(indCell.validCell, 1), 'b.', 'MarkerSize', 10)
    hold on;
    plot(center(indCell.excCell, 2), center(indCell.excCell, 1), 'r.', 'MarkerSize', 10)
    set(gca, 'YDir', 'reverse');
    xlim([1 size(neuron.Cn, 2)])
    ylim([1 size(neuron.Cn, 1)])
    setCell = cat(1, indCell.validCell, indCell.excCell);
    text(center(setCell, 2)+1, center(setCell, 1), num2str(setCell), 'Color', 'k');
    title(sprintf('%s %s: valid cells in blue', nameSubj, dateSession));
    print(figCenter2, sprintf(fullfile(dirProcdata_session, '%s_%s_validCells_2'), nameSubj, dateSession), '-depsc')
    
end


