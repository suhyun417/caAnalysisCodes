% computeDFF_meanF.m
%
% compute delta F over F using mean F over time (within run) as a
% denominator

addpath('/projects/parksh/_toolbox/TIFFstack');
addpath('/projects/parksh/_toolbox/Fast_Tiff_Write/');
addpath('/projects/parksh/_toolbox/imagetools/');

setSubj = {'Max'}; % {'Tabla', 'Max'}; % 'Max'; % 'Tabla';
dateSession = '20191223'; % '20191113'; %'20191125'; %'20191126'; %'20191114'; %'20191112';

for iSubj = 1:length(setSubj)
    nameSubj = setSubj{iSubj};
%     dirRawData_session = fullfile('/rawdata/parksh/calciumImaging/', [dateSession, '_', nameSubj]); %20180529_Hoppy';
    dirPreProcData_session = fullfile('/procdata/parksh/_marmoset/invivoCalciumImaging/', nameSubj, dateSession, '_preproc');   
%     % Optional params
%     flagSaveFigure = 1;
%     dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';
    
    % List of runs to process
    d_mat = dir(fullfile(dirPreProcData_session, '*RigidMC.mat'));
%     [listRun{1:length(d_mat)}] = deal(d_mat.name);
    
    nRun = length(d_mat);
    for iRun = 9:nRun
        [p, nameRun, ext] = fileparts(d_mat(iRun).name);
        saveFileName = fullfile(dirPreProcData_session, [nameRun '_dFF.mat']);
        saveFileName_tif = fullfile(dirPreProcData_session, [nameRun '_dFF.tif']); 
        
        Mr = loadtiff(fullfile(dirPreProcData_session, [nameRun, '.tif']));
%         load(fullfile(dirPreProcData_session, d_mat(iRun).name), 'Mr');
        
        Y_avg = mean(Mr, 3);
        Y = (Mr - repmat(Y_avg, 1, 1, size(Mr,3)))./repmat(Y_avg, 1, 1, size(Mr,3));
        Y_std = std(Y, [], 3);
        Ysiz = size(Y)'; % [d1, d2, T]'; % following CNMFe's mat file convention
        
        save(saveFileName, 'Y_avg', 'Y', 'Y_std', 'Ysiz');
        fprintf(1, ':: Run #%d/%d: dFF of run %s from session %s_%s is saved as .mat file\n', iRun, nRun, dateSession, nameSubj, d_mat(iRun).name)
        fastTiffStackWrite(saveFileName_tif, single(Y));
        fastTiffStackWrite(fullfile(dirPreProcData_session, [nameRun '_dFF_std.tif']), Y_std);
        fprintf(1, ':: Run #%d/%d: dFF of run %s from session %s_%s is saved as .tif file\n', iRun, nRun, dateSession, nameSubj, d_mat(iRun).name)
        
        clear Mr Y_avg Y_dFF Y_std
    end
end
    
        


% % Mr = Mr(:,:,1:600); % get rid of severe head motion if needed
% 
% avgMR = mean(Mr, 3);
% dFF = (Mr - repmat(avgMR, 1, 1, size(Mr,3)))./repmat(avgMR, 1, 1, size(Mr,3));
% aa = strsplit(filename, '.');
% filehead = aa{1};
% saveFileName = [filehead, '_dFF'];
% saveastiff(dFF, fullfile(dirProcData_session, saveFileName))