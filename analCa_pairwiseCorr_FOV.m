% 예시 데이터 로딩
load('~/data/procdata/parksh/_marmoset/invivoCalciumImaging/XYPC_tabla.mat')
% load('~/data/procdata/parksh/_marmoset/invivoCalciumImaging/XYPC_max.mat')

scale = 3.125;
bin_edges = 0:5:300;

% 예: PC에 대해 시각화
[c1, sim1, counts1] = correlation_vs_distance(x, y, PC1, bin_edges, scale);
[c2, sim2, counts2] = correlation_vs_distance(x, y, PC2, bin_edges, scale);

% % optional
% disp(table(c1', sim1', counts1', ...
%     'VariableNames', {'Distance', 'Similarity', 'N_paris'}))

% Plot
figure;
plot(c1, sim1, 'bo-', 'LineWidth', 2); hold on
plot(c2, sim2, 'ro-', 'LineWidth', 2);

xlabel('Pairwise Distance (µm)');
ylabel('Mean Correlation');
title('Similarity vs Spatial Distance');
grid on;
legend('PC1', 'PC2')

% 
% function [bin_centers, sim_by_bin] = correlation_vs_distance(x, y, PC, bin_edges, scale)
% % pairwise correlation as a function of distance between two neurons
% % input
% % x, y
% % PC
% % bin_edges
% % scale
% % output
% % bin_centers: center of each distance bin (um)
% % sim_by_bin: mean correlation within bin
% 
% if nargin < 5
%     scale = 3.125;
% end
% 
% % conversion
% X = [x(:) * scale, y(:) * scale];
% n = size(X, 1);
% PC = PC(:);
% 
% % 거리 행렬 계산
% dist_mat = squareform(pdist(X));
% sim_mat = 1 - squareform(pdist(PC, 'correlation'));
% 
% % 자기 자신은 제외
% sim_mat(eye(n) == 1) = NaN;
% 
% % 거리 bin 계산
% bin_centers = 0.5 * (bin_edges(1:end-1) + bin_edges(2:end));
% sim_by_bin = nan(size(bin_centers));
% 
% for i = 1:length(bin_centers)
%     mask = dist_mat >= bin_edges(i) & dist_mat < bin_edges(i+1);
%     values = sim_mat(mask);
%     sim_by_bin(i) = mean(values, 'omitnan');
% end
% end

