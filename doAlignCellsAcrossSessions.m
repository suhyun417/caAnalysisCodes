function [] = doAlignCellsAcrossSessions(nameSubj, FOV_ID, flagSaveFile)

% 2022/10/24 SHP
%   -Align cells across daily sessions for a given FOV
%   using the shifts saved in SUBJ_FOV#_shifts.mat
%   modified from doAssembleCellCenterInfoAcrossSessions.m
%   1. Load "neuron.A" from each session, spatial smoothing, apply threshold, 
%   apply shifts, then assign the cell ID to the selected (above threshold)
%    indices. Save this cell ID spatial component matrix for each
%   session.
%   2. Start from Cell 1 of Session 1, gather the cell ID from each session 
%   from the corresponding spatial location from the
%   matrix saved in Step #1. If there is at least one cell matched to the
%   location of this particular cell, record the cell ID(s) in a separate
%   matrix (cellIDAcrossDay).
%   -Results are saved as SUBJ_FOV#_cellAcrossDay.mat
% 


%% settings
flagBiowulf = 0; %1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/NIFVAULT/projects/parksh';
        dirProcdata = '/Volumes/NIFVAULT/procdata/parksh';
        dirRawdata = '/Volumes/NIFVAULT/rawdata/parksh';
    else % on virtual machine
        dirProjects = '/nifvault/projects/parksh';
        dirProcdata = '/nifvault/procdata/parksh';
        dirRawdata = '/nifvault/rawdata/parksh';
    end
end
% dirFig = fullfile(dirProjects, '0Marmoset/Ca/_labNote/_figs/');

addpath(fullfile(dirProjects, '_toolbox/TIFFstack'));
addpath(fullfile(dirProjects, '_toolbox/NoRMCorre/'));
addpath(fullfile(dirProjects, '_toolbox/Fast_Tiff_Write/'));
addpath(fullfile(dirProjects, '_toolbox/imagetools/'));
addpath(fullfile(dirProjects, '_toolbox/CNMF_E/'));
cnmfe_setup;
% gcp; % for parallel processingls


%% Session info & optional parameters
% nameSubj = 'Tabla';
% FOV_ID = 1;
% get session info
[infoSession, opts] = readInfoSession(nameSubj, FOV_ID);

[c, ia, indRun] = unique(infoSession.(1), 'sorted');
setDateSession = c(2:end); % 1st one is always empty
nSession = length(setDateSession);

%% load necessary files 
% registration info
fname_shifts = fullfile(dirProjects, sprintf('0Marmoset/Ca/tempData/%s_FOV%d_shifts.mat', nameSubj, FOV_ID));
load(fname_shifts, 'shifts')

% cell quality info 
fname_cellQC = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellQC.mat',...
    nameSubj, FOV_ID, nameSubj, FOV_ID)); 
load(fname_cellQC, 'infoCells')

%% apply shift to each cell's spatial component, then save the registered spatial info in matrix
% stackCellCenter = [];
cellAcrossDay = struct([]);
for iSession = 1:nSession
    
    fprintf(1, '\n Processing %s session %d/%d...\n', nameSubj, iSession, nSession)

    d_sources2D = dir(fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/Session/%s/Sources2D_all*',...
        nameSubj, setDateSession{iSession})));
    load(fullfile(d_sources2D(1).folder, d_sources2D(1).name));
    
    thr = 0.2; %0.3; % the lower the smaller (more centralized) the contour   
    [d1,d2] = size(neuron.Cn); %d1: FOV vertical, d2: FOV horizontal
%     Coor = neuron.get_contours(thr); 
    
%     figure;
%     for i = 1:size(Coor, 1)
%         %         cont = medfilt1(Coor{i}')';
%         cont = Coor{i};
%         if size(cont,2) > 1 % "shifts" values are in the order of image matrix dimensions: first dimension is vertical ("y") and second dimension is horizontal ("x")
%             plot(cont(1,1:end)+shifts(iSession, 2), cont(2,1:end)+shifts(iSession, 1), 'Color', cellColor(iSession, :), 'linewidth', widthContour); hold on;
%         end
%         set(gca, 'YDir', 'reverse', 'XLim', [0 d2], 'YLim', [0 d1])
%     end
    
%     curCanvas = NaN(size(neuron.Cn));
    stackCell = NaN(size(neuron.A));
    invalidCell = []; % any cell that may not pass the below procedure
    validCell = intersect(intersect(infoCells(iSession).indCellValid_spatialCluster_0p2, infoCells(iSession).indCellValid_fov), ...
        infoCells(iSession).indCellValid_snr);
    
    for i = 1:size(neuron.A ,2)
        if ~ismember(i, validCell)
            continue;
        end
        A_temp = full(reshape(neuron.A(:,i),d1,d2));
        A_temp = medfilt2(A_temp,[3,3]);
        A_temp = A_temp(:);
        [temp,ind] = sort(A_temp(:).^2,'ascend');
        temp =  cumsum(temp);
        ff = find(temp > (1-thr)*temp(end),1,'first'); % cumulative index
        
        if A_temp(ind(ff)) < 0 
            invalidCell = cat(1, invalidCell, i);
            continue;
        end
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        
        % apply shifts from across-session registration
        clear ii_shift jj_shift fp_shift
        ii_shift = round(ii+shifts(iSession,1)); 
        jj_shift = round(jj+shifts(iSession,2)); % 
        
        if min([ii_shift;jj_shift])<=0 || max(ii_shift)>d1 || max(jj_shift)>d2 % if out of range
            locvalidpix = sum(cat(2, ii_shift>0, ii_shift<d1, jj_shift>0, jj_shift<d2), 2)>3;
            ii_shift = ii_shift(locvalidpix);
            jj_shift = jj_shift(locvalidpix);
        end
        fp_shift = sub2ind([d1, d2], ii_shift, jj_shift);
                
%         curCanvas(fp_shift) = i;
        stackCell(fp_shift, i) = i;
    end
    
    cellAcrossDay(iSession).stackCell = stackCell;
    cellAcrossDay(iSession).invalidCell = invalidCell;
end
  
 
%% Make a session - cell id matrix
idMatrix = [];
for iS = 1:nSession
    for iCell = 1:size(cellAcrossDay(iS).stackCell, 2)
        cellCount = size(idMatrix, 1);
        idMatrix(cellCount+1, 1) = iS;
        idMatrix(cellCount+1, 2) = iCell;
    end
end

%% create a checkbox
flagDone = NaN(size(idMatrix, 1), 1);
for iCell = 1:size(idMatrix,1)
    idSession = idMatrix(iCell, 1);
    idCell = idMatrix(iCell, 2);
    if ismember(idCell, intersect(intersect(infoCells(idSession).indCellValid_spatialCluster_0p2, infoCells(idSession).indCellValid_fov), ...
        infoCells(idSession).indCellValid_snr))...
            && ~ismember(idCell, cellAcrossDay(idSession).invalidCell)
        flagDone(iCell) = 0;
    end
end

%% For each cell, check the corresponding cells from other sessions
bigStackCell = cat(2, cellAcrossDay(:).stackCell); % pixel by cells (from all the sessions)
bigStackCell(:, isnan(flagDone)) = NaN; % invalid cells
cellIDAcrossDay = [];
countValidCell = 0;
for iCell = 1:length(flagDone)
%     idSession = idMatrix(iCell, 1);
%     idCell = idMatrix(iCell, 2);
    
    if isnan(flagDone(iCell)) || flagDone(iCell) % if it's not valid or checked already, skip
        continue;
    end
    
    countValidCell = countValidCell + 1;    

    % current cell's spatial location
    curSite = ~isnan(bigStackCell(:,iCell));
    
    % look for other cells occupying here
    candidateCells = ~isnan(min(bigStackCell(curSite, :)));
      
    if sum(candidateCells, 2)>1      
        
        tempID = NaN(1, length(cellAcrossDay)); % one row for each
        
        locCandidateCells = find(candidateCells>0);
        
        % check whether it's more than one cell from each session
        tempSession = idMatrix(candidateCells, 1);
        [a b c] = unique(tempSession);
        if length(c) > length(a) % if there are more than one cell from one session captured
            
            tempSizeCell = sum(~isnan(bigStackCell(curSite, candidateCells)));
            
            indValid_candidateCells = [];
            for iSS = 1:length(a) % unique sessions
                if sum(tempSession==a(iSS))>1 % more than one cell from this session
                    setDup = find(tempSession==a(iSS));
                    [~, ind] = max(tempSizeCell(setDup)); % the one overlaps with the current cell most
                    indValid_candidateCells = cat(2, indValid_candidateCells, setDup(ind));
                else
                    indValid_candidateCells = cat(2, indValid_candidateCells, b(iSS));
                end
            end
                   
            locCandidateCells = locCandidateCells(indValid_candidateCells);            
%             fprintf(1, 'iCell %d: more than one cell found', iCell);
%             idMatrix(candidateCells, :)
%             input('')
        else
            locCandidateCells = find(candidateCells>0);
        %         tempIDCandidate = idMatrix(candidateCells, :);
        end   
        
        tempID(idMatrix(locCandidateCells, 1)) = idMatrix(locCandidateCells, 2);
        
        cellIDAcrossDay(countValidCell, :) = tempID;
        flagDone(locCandidateCells) = 1;
        bigStackCell(:, locCandidateCells) = NaN;
        
    else % if no overlapping cells
        if contains(lower(nameSubj), 'tab') && iCell == 211 % manual fixing for Tabla's cell
            locCandidateCells = [211 295 423 522 706 1014 1128];
            tempID = [NaN 87 67 79 83 NaN 49 NaN NaN 72 76 NaN];
            
            cellIDAcrossDay(countValidCell, :) = tempID;
            flagDone(locCandidateCells) = 1;
            bigStackCell(:, locCandidateCells) = NaN;
            
        elseif contains(lower(nameSubj), 'tab') && iCell == 212 % manual fixing for Tabla's cell
            locCandidateCells = [212 300 745];
            tempID = [NaN 88 72 NaN NaN NaN 88 NaN NaN NaN NaN NaN];
            
            cellIDAcrossDay(countValidCell, :) = tempID;
            flagDone(locCandidateCells) = 1;
            bigStackCell(:, locCandidateCells) = NaN;
            
        else        
            tempID = NaN(1, length(cellAcrossDay));
            tempID(idMatrix(iCell, 1)) = idMatrix(iCell, 2); % filling this cell
                        
            cellIDAcrossDay(countValidCell, :) = tempID;
            flagDone(iCell) = 1;
            bigStackCell(:, iCell) = NaN;
        end
    end
end

%
paramAlignCells.idMatrix = idMatrix;
paramAlignCells.flagDone = flagDone;



if flagSaveFile
%     fname_stack = fullfile(dirProjects, sprintf('0Marmoset/Ca/tempData/%s_FOV%d_stackedCenter.mat', nameSubj, FOV_ID));
    fname_stack = fullfile(dirProcdata, sprintf('_marmoset/invivoCalciumImaging/%s/FOV%d/%s_FOV%d_cellAcrossDay.mat', nameSubj, FOV_ID, nameSubj, FOV_ID));
    save(fname_stack, 'cellAcrossDay', 'paramAlignCells', 'cellIDAcrossDay')
    fprintf(1, '\n Saving files for %s FOV %d: in total %d/%d cells aligned \n', nameSubj, FOV_ID, countValidCell, length(flagDone))
end
  
    
    

end % end of function