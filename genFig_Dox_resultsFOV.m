% genFig_Dox_resultsFOV.m
%
% 2020/04/24 SHP

setNameSubj = {'Tabla', 'Max'};
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

orderROI{1} = [1 6 2 4 5 7 3];
orderROI{2} = [2 3 4 1];
cMap_ROI{1} = hsv(7).*0.7;
cMap_ROI{2} = hsv(4).*0.7;

for iSubj = 1:2
nameSubj = setNameSubj{iSubj};
dirDox = sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Dox', nameSubj);
fname = sprintf('doxResults_FOVROI_%s.mat', nameSubj);
load(fullfile(dirDox, fname))

setDate = cat(1, paramSession.dateSession);
setDate_vec = datenum(setDate, 'yyyymmdd');
xAxis_date = setDate_vec - setDate_vec(1);

mF_ROI = cat(1, resultsROI.meanF_ROI);
mF = cat(1, resultswholeFOV.meanF);

legendCell = strcat('ROI', string(num2cell(1:size(mF_ROI, 2))));
legendCell{size(mF_ROI, 2)+1} = 'Whole FOV';

figure;
set(gcf, 'Color', 'w', 'Position', [1300 200 440 500])
set(gca, 'ColorOrder', cMap_ROI{iSubj}); hold on;
plot(xAxis_date, mF_ROI(:, orderROI{iSubj}), 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on;
plot(xAxis_date, mF, 'ko-', 'LineWidth', 4, 'MarkerFaceColor', 'w', 'MarkerSize', 8);
plot(xAxis_date(12), max(mF_ROI(12,:))+15, 'rv', 'MarkerFaceColor', 'r')
text(xAxis_date(12)+1, max(mF_ROI(12,:))+10, '2nd Dox (2mg/kg)', 'Color', 'r', 'VerticalAlignment', 'bottom')
legend(legendCell, 'Location', 'Best')
title(sprintf('Average fluorescence change: %s', nameSubj))
xlabel('Days after 1st doxycycline administration')
ylabel('Fluorescence (a.u.)')
set(gca, 'Box', 'off', 'TickDir', 'out')
print(gcf, fullfile(dirFig, sprintf('%s_dox_resultsFOVROI_Reorder', nameSubj)), '-depsc');

end

setSession{1} = 1:12;
setSession{2} = 12:16;
for iSubj = 1:2
    nameSubj = setNameSubj{iSubj};
    dirDox = sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Dox', nameSubj);
    fname = sprintf('doxResults_FOVROI_%s.mat', nameSubj);
    load(fullfile(dirDox, fname))
    for iDose = 1:2
        setDate_vec=[];
        
        dox(iSubj, iDose).mF = cat(1, resultswholeFOV(setSession{iDose}).meanF);
        dox(iSubj, iDose).mF_norm = dox(iSubj, iDose).mF./dox(iSubj, iDose).mF(1);
        dox(iSubj, iDose).mF_ROI = cat(1, resultsROI(setSession{iDose}).meanF_ROI);
        dox(iSubj, iDose).mF_ROI_norm = dox(iSubj, iDose).mF_ROI./dox(iSubj, iDose).mF_ROI(1,:);
        dox(iSubj, iDose).setDate = cat(1, paramSession(setSession{iDose}).dateSession);
        setDate_vec = datenum(dox(iSubj, iDose).setDate, 'yyyymmdd');
        dox(iSubj, iDose).xAxis_date = setDate_vec - setDate_vec(1);
    end
end

% Dose dependency: Whole FOV
fig_dose = figure;
set(fig_dose, 'Color', 'w', 'Position', [1300 200 440 500]);
cMap_subj = [0.5 0 0; 0 0 0.8]; %jet(size(dox,1));
lineStyle_dose = {'-', ':'};
% markerfacecolor_dose = {'
for iSubj = 1:2
    for iDose = 1:2
        figure(fig_dose);
        plot(dox(iSubj, iDose).xAxis_date, dox(iSubj, iDose).mF_norm.*100, 'o-', 'Color', cMap_subj(iSubj,:), 'LineStyle', lineStyle_dose{iDose}, ...
            'LineWidth', 3, 'MarkerFaceColor', 'w', 'MarkerSize', 8);
        hold on;
    end
end
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2)
ylim([70 102])
xlabel('Days after doxycycline administration')
ylabel('Relative fluorescence (%)')
% L=legend('Tabla: high dose', 'Tabla: low dose', 'Max: high dose', 'Max: low dose');
print(fig_dose, fullfile(dirFig, 'dox_resultsWholeFOV_dose_noLegend'), '-depsc');

% Dose dependency: ROI
% setMarkers = {'o', '^', 'sq', 'x', '+'};
orderROI{1} = [1 6 2 4 5 7 3];
orderROI{2} = [2 3 4 1];
cMap_ROI{1} = hsv(7).*0.7;
cMap_ROI{2} = hsv(4).*0.7;
lineStyle_dose = {'-', ':'};

ylim_subj{1} = [60 100];
ylim_subj{2} = [70 102];

for iSubj = 1:2
    fig_dose_ROI = figure;
    set(fig_dose_ROI, 'Color', 'w', 'Position', [200 200 440 500]);
    set(gca, 'ColorOrder', cMap_ROI{iSubj});
    hold on;
    
    for iDose = 1:2
        figure(fig_dose_ROI);
        plot(dox(iSubj, iDose).xAxis_date, dox(iSubj, iDose).mF_ROI_norm(:, orderROI{iSubj}).*100, 'o-', 'LineStyle', lineStyle_dose{iDose}, ...
            'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8);
        hold on;
    end
    xlabel('Days after doxycycline administration')
    ylabel('Relative fluorescence (%)')
    title(sprintf('Average fluorescence change for each ROI: %s', setNameSubj{iSubj}))
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2)
    ylim(ylim_subj{iSubj});
    print(fig_dose_ROI, fullfile(dirFig, sprintf('dox_resultsFOVROI_dose_%s', setNameSubj{iSubj})), '-depsc');
end




