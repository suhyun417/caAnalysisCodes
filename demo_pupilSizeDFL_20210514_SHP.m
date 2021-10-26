% demo_pupilSizeDFL_20210514_SHP.m
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
    
    indDFLRuns = contains(setMLFilename, 'DFL') & cat(1, S.flagPreproc)>0; %% containing "DFL" in filename AND flagPreproc value of 1
    setFilename = setMLFilename(indDFLRuns);
    
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
        
        %% pupil data
        xData = data.AnalogData.Eye(:,1);
        yData = data.AnalogData.Eye(:,2);
        pData = data.AnalogData.General.Gen1;
        
        % temporary plot to check the data
        selPeriod = 20001:21000; % time points that we want to look at
        for i = 1:length(xData)-1
            figure(101);
            plot(xData(i:i+1), yData(i:i+1), 'b-');
            xlim([8 10])
            ylim([2 7])
            hold on;
%             input('')
        end
            
        figure;
        plot(xData, yData, '-');
        plot(xData(selPeriod), yData(selPeriod), '-');
        
        
        
        %% NOTE: some noise possibly need to be filtered out
        %% Check the noise using FFT
        Fs = 1000; % sampling rate
        T = 1/Fs;
        L = length(xData);             % Length of signal
        t = (0:L-1)*T; % time vector
        
        Y = fft(xData);
        % Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
        
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        % Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
        
        f = Fs*(0:(L/2))/L;
        figure;
        plot(f,P1)
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        
        %% NOTE: information about video contents should be added during the analysis
        % (contents vary depending on the videofile (movie 1_1 vs. movie
        % 5_1 have different order) HF, MF, SCN, OBJ motion, OF, MB
        
        
    end
    
% end


















