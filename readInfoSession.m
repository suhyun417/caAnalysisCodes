function [infoT, opts] = readInfoSession(nameSubj, FOV_ID)

% Read a spreadsheet from a .xls file containing calcium imaging session information
% into a table, using "readtable" function
%   Usage: [infoT, opts] = readInfoSession(nameSubj);
%   Input
%       - nameSubj: full name of subject (string)
%       - FOV_ID: ID of FOV, in number (double)
%   Output
%       - infoT: table containing the full contents
%       - opts: options for the table
% 2020/01/22 SHP
% 2021/12/13 SHP: modify to retrieve the FOV info from the input


ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/VNDLab_Data/projects/parksh'; %'/Volumes/NIFVAULT/PROJECTS/parksh';
    dirProcdata = '/Volumes/VNDLab_Data/procdata/parksh'; %'/Volumes/NIFVAULT/PROCDATA/parksh';
    dirRawdata = '/Volumes/VNDLab_Data/rawdata/parksh'; %'/Volumes/rawdata/parksh';
else % on virtual machine
    dirProjects = '/VNDLab_Data/projects/parksh'; %'/nifvault/projects/parksh';
    dirProcdata = '/VNDLab_Data/procdata/parksh'; %'/nifvault/procdata/parksh';
    dirRawdata = '/VNDLab_Data/rawdata/parksh'; %'/nifvault/rawdata/parksh';
end


% get session info from spreadsheet
fname_infoSession = fullfile(dirProjects, sprintf('0Marmoset/Ca/infoSession_%s_FOV%d.xlsx', nameSubj, FOV_ID));

% switch lower(nameSubj)
%     case 'tabla'
%         fname_infoSession = fullfile(dirProjects, '0Marmoset/Ca/infoSession_Tabla_FOV1.xlsx');
% %         fname_infoSession = fullfile(dirProjects, '0Marmoset/Ca/infoSession_Tabla_FOV3.xlsx');
%     case 'max'
%         fname_infoSession = fullfile(dirProjects, '0Marmoset/Ca/infoSession_Max_FOV3.xlsx');
% %         fname_infoSession = fullfile(dirProjects, '0Marmoset/Ca/infoSession_Max_FOV2.xlsx');
% end

opts = detectImportOptions(fname_infoSession);
opts = setvartype(opts, {'Date', 'ImagingFilename', 'ExpName', 'MLFilename', 'stimulus', 'notes'}, 'char');
opts = setvartype(opts, 'flagPreproc', 'double');
infoT = readtable(fname_infoSession, opts);