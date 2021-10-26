% genFig_Dox_resultsNeuron.m
%
% 2020/04/24 SHP

setNameSubj = {'Tabla', 'Max'};
dirFig = '/projects/parksh/0Marmoset/Ca/_labNote/_figs/';

setSession{1} = 1:12;
setSession{2} = 12:16;

for iSubj = 1:2
    nameSubj = setNameSubj{iSubj};
    dirDox = sprintf('/procdata/parksh/_marmoset/invivoCalciumImaging/%s/Dox', nameSubj);
    fname = sprintf('doxResults_neuron_%s.mat', nameSubj);
    load(fullfile(dirDox, fname))
    
     for iDose = 1:2
        setDate_vec=[];
        
        dox_neuron(iSubj, iDose).nCell = cat(1, resultsDox_neuron(setSession{iDose}).nCell);
        dox_neuron(iSubj, iDose).nCell_norm = dox_neuron(iSubj, iDose).nCell./dox_neuron(iSubj, iDose).nCell(1);
        
%         dox_neuron(iSubj, iDose).min_pnr = cat(1, resultsDox_neuron_paramCNMFE(setSession{iDose}).min_pnr);
%         dox_neuron(iSubj, iDose).min_pnr_norm = dox_neuron(iSubj, iDose).min_pnr./dox_neuron(iSubj, iDose).min_pnr(1);
%         dox_neuron(iSubj, iDose).min_corr = cat(1, resultsDox_neuron_paramCNMFE(setSession{iDose}).min_corr);
%         dox_neuron(iSubj, iDose).min_corr_norm = dox_neuron(iSubj, iDose).min_corr./dox_neuron(iSubj, iDose).min_corr(1);
        
        dox_neuron(iSubj, iDose).setDate = cat(1, resultsDox_neuron(setSession{iDose}).dateSession);
        setDate_vec = datenum(dox_neuron(iSubj, iDose).setDate, 'yyyymmdd');
        dox_neuron(iSubj, iDose).xAxis_date = setDate_vec - setDate_vec(1);
     end
    
end


% Dose dependency: number of cells
fig_dose = figure;
set(fig_dose, 'Color', 'w', 'Position', [1300 200 440 500]);
cMap_subj = [0.5 0 0; 0 0 0.8]; %jet(size(dox,1));
lineStyle_dose = {'-', ':'};
% markerfacecolor_dose = {'
for iSubj = 1:2
    for iDose = 1:2
        figure(fig_dose);
        plot(dox_neuron(iSubj, iDose).xAxis_date, dox_neuron(iSubj, iDose).nCell_norm.*100, 'o-', 'Color', cMap_subj(iSubj,:), 'LineStyle', lineStyle_dose{iDose}, ...
            'LineWidth', 3, 'MarkerFaceColor', 'w', 'MarkerSize', 8);
        hold on;
    end
end
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2)
ylim([70 102])
xlabel('Days after doxycycline administration')
ylabel('Number of neurons')
L=legend('Tabla: high dose', 'Tabla: low dose', 'Max: high dose', 'Max: low dose');


