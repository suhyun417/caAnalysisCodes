function [directory] = setDir_shp()

% flagBiowulf = 0; %1; %0;
% 
% if flagBiowulf
%     directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
%     directory.dirFig = '/data/parks20/analysis/_figs';
% else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.dirProjects = '/Volumes/VNDLab_Data/projects/parksh'; %'/Volumes/NIFVAULT/PROJECTS/parksh';
        directory.dirProcdata = '/Volumes/VNDLab_Data/procdata/parksh'; %'/Volumes/NIFVAULT/PROCDATA/parksh';
        directory.dirRawdata = '/Volumes/VNDLab_Data/rawdata/parksh'; %'/Volumes/rawdata/parksh';
    else % on virtual machine
        directory.dirProjects = '/VNDLab_Data/projects/parksh'; %'/nifvault/projects/parksh';
        directory.dirProcdata = '/VNDLab_Data/procdata/parksh'; %'/nifvault/procdata/parksh';
        directory.dirRawdata = '/VNDLab_Data/rawdata/parksh'; %'/nifvault/rawdata/parksh';
    end
% end