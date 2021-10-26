function [] = doConcatRuns(listFileName, fname_cat)

% Concatenate files and save it as .mat and .tif
%   Usage:  [] = doConcatRuns(listFileName, fnameOut);
%   Input:
%       - listFileName: Cell array containing the list of filenames to be
%       concatenated. Each filename should contain the full path
%       - fnameOut: filename for the output (concatenated) file with
%       full path
%
% 2020/01/22 SHP

% clear all;
addpath('/projects/parksh/_toolbox/TIFFstack');
addpath('/projects/parksh/_toolbox/imagetools/')

fname_catmat = [fname_cat '.mat'];
fname_cattif = [fname_cat, '.tif'];
opts_tiff.append = true;
opts_tiff.big = true;
% d = dir(fname_catmat);
% 
% if isempty(d)
%     fprintf(1, 'Concatenated .mat file for all runs does not exist. Creating one now...\n')
    
    %     Y = [];
    %     data = matfile(fullfile(dirPreproc, fname_catmat), 'Writable', true);
    count = 0; % for concatenation
%     fprintf(1, 'The concatenated file is being saved as %s...\n', fname_catmat)
    
    nFile = length(listFileName);
    %      d_all = dir(fullfile(dirPreproc, '*_sDS_cat.mat'));
    
    for iFile = 1:length(listFileName)
        
        d_file = dir(listFileName{iFile});
        
        fprintf(1, '      Loading file #%d/%d: %s \n', iFile, length(listFileName), d_file.name)
        Yf_cat = loadtiff(fullfile(d_file.folder, d_file.name));
        Yf_cat = single(Yf_cat);
        
        % write the file into a single mat file
        [d1, d2, T] = size(Yf_cat);
        tic;
        if iFile == 1
            Y = Yf_cat;
            save(fname_catmat, 'Y', '-v7.3');
        else
            data = matfile(fname_catmat, 'Writable', true);
            data.Y(:, :, count+1:count+T) = Yf_cat;
        end
        toc
        
%         saveastiff(Yf_cat, fname_cattif, opts_tiff); 
        %         Y = cat(3, Y, Mr);
        clear Y Yf_cat
        count = count+T;
        
    end
    data = matfile(fname_catmat, 'Writable', true);
    Ysiz  = size(data, 'Y')'; % [d1, d2, T]'; % following CNMFe's mat file convention
    save(fname_catmat, 'Ysiz', '-append'); %, '-v7.3')
    fprintf(1, '...Done!\n')
%     fprintf(1, 'Saving data in TIFF now...')
    
    %     tic; saveastiff(data.Y, fullfile(dirPreproc, [fname, '.tif'])); toc;
%     tic; fastTiffStackWrite([fname_cat, '.tif'], data.Y); toc;
%     fprintf(1, '...Done!\n')
    %     Ysiz = size(Y)'; % [d1, d2, T]'; % following CNMFe's mat file convention
    
    %     fprintf(1, 'The concatenated file is being saved as %s...\n', fname_catmat)
    %     save(fullfile(dirPreproc, fname_catmat), 'Y', 'Ysiz', '-v7.3')
    
end


