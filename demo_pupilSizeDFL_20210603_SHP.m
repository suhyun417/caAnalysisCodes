% demo_pupilSizeDFL_20210603_SHP.m
%
% Read pupil size data from DFL session and start from there
%

setSubjName = {'Tabla', 'Max'};
iSubj = 1; %
% for iSubj = 1:length(setSubjName)
    clear Results*
    
    nameSubj =  setSubjName{iSubj}; %'Max'; %'Tabla'; %'Max'; %'Tabla'; %'Max'; %'Tabla';
    % dateSession = '20191113'; %'20191125';
    
    % get session info
    [infoSession, opts] = readInfoSession(nameSubj);
    S = table2struct(infoSession);
    
    % setExpName = {S.ExpName}';
    setMLFilename = {S.MLFilename}';
    
    indDFLRuns = contains(setMLFilename, 'DFL') & cat(1, S.flagPreproc) > 0 & contains({S.stimulus}', 'set1_1'); %% containing "DFL" in filename AND flagPreproc value of 1
    setFilename = setMLFilename(indDFLRuns);
    
    % Tabla
    % 191113: eye signal not good (another range of x/y/pupil present during the first half)
    % 191114 & 191118: eye signal good. Length of data points are 6-7ms shorter than 2min
    % except the last run of 191118
    % 191119: starting to lose x signal + y & pupil is bad in general
    % 191120: x is lost, y& pupil was good in the first run, but not in the
    % others
    % 191121: x is lost, length of data points are 6-7ms shorter than 2min, y
    % & pupil bad in some runs
   % 191125: x is lost, for the first file signals look okay & length is fine 
   
    
    for iFile = 1:length(setFilename)
        filename = strcat(setFilename{iFile}, '.bhv2');
        
        dateSession = filename(1:6);
        
        if str2num(dateSession) < 191121
            dirBHV = '/archive_rawdata1/parksh/behavior/MonkeyLogic_Ca/'; %
        else
            dirBHV = '/rawdata/parksh/behavior/MonkeyLogic_Ca/'; %
        end
        
        % filename = '191121_Tabla_Ca_BPM_123909.bhv2'; % change it to a file you have
        
        %% Read the file
        data = mlread(fullfile(dirBHV, filename)); % mlread(filename);
        
        %% eye data during stimulus on 
        % Event Code Numbers & Names : TASK_START = 10; FP_ON = 20;
        % WAIT_FOR_TR = 30; MOVIE_ON = 40; REWARD = 90; 
        % TRIG onset = 900; TRIG offset = 990; (TTL from ML to Inscopix DAQ On & Off)
        locStimOn = find(data.BehavioralCodes.CodeNumbers == 40);
        time_stimOn = floor(data.BehavioralCodes.CodeTimes(locStimOn))
        data.AnalogData
        
        figure;
        plot(data.AnalogData.Eye, '.')
        line([time_stimOn time_stimOn], get(gca, 'YLim'), 'Color', 'm')
        title(sprintf('%s: x & y gaze', filename))
        axis tight
        xlabel('Time')
        ylabel('Voltage')
        legend('X', 'Y')
        
        figure;
        plot(data.AnalogData.General.Gen1, 'g.')
        line([time_stimOn time_stimOn], get(gca, 'YLim'), 'Color', 'm')
        title(sprintf('%s: pupil size change', filename))
        axis tight
        xlabel('Time')
        ylabel('Voltage')
        
        input('')
    end
       
        durMovie_ms = 120000;
        xData = data.AnalogData.Eye(time_stimOn+1:time_stimOn+durMovie_ms, 1);
        yData = data.AnalogData.Eye(time_stimOn+1:time_stimOn+durMovie_ms, 2);
        pData = data.AnalogData.General.Gen1(time_stimOn+1:time_stimOn+durMovie_ms);
        
        
        %% Reshape the eye data into matrix of time x conditions
        x_condition = reshape(xData, 20000, 6);
        y_condition = reshape(yData, 20000, 6);
        pupil_condition = reshape(pData, 20000, 6);
        
        % quick check on eye data from different conditions
        for iCond = 1:6
            figure(101);
            subplot(1, 2, 1); cla;
            plot(x_condition(:, iCond), y_condition(:, iCond), '-')
            xlabel('x')
            ylabel('y')
%             title(sprintf('%s: condition ID %d', filename, iCond))
            
            subplot(1, 2, 2); cla;
            plot(x_condition(:, iCond), 'b.'); hold on;
            plot(y_condition(:, iCond), 'r.');
            plot(pupil_condition(:, iCond), 'g.');
            title(sprintf('%s: condition ID %d', filename, iCond))
            legend('x', 'y', 'pupil')
            xlabel('Time (ms)')
            
            input('')
        end
        
        % check x,y gaze & pupil change together
        for iCond = 1:6
            figure(200);
            clf;
            plot(x_condition(:, iCond), 'b.'); hold on;
            plot(y_condition(:, iCond), 'r.');
            plot(pupil_condition(:, iCond), 'g.');
            title(sprintf('%s: condition ID %d', filename, iCond))
            input('')
        end
        
        
        
%         % histogram to check criterion for selecting valid eye signal
%         figure;
%         histogram(x_condition);
        
        
        %% NOTE: information about video contents should be added during the analysis
        % (contents vary depending on the videofile (movie 1_1 vs. movie
        % 5_1 have different order) 
        % movie 1_1: scene, human face, optic flow, marmoset body, cars, marmoset faces
      
        
        %% NOTE: some noise can be filtered out (e.g. wdenoise)
        
                
        
    end
    
% end


















