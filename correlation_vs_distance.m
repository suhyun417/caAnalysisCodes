function [bin_centers, sim_by_bin, bin_counts] = correlation_vs_distance(x, y, PC, bin_edges, scale)
% pairwise correlation as a function of distance between two neurons
% input
% x, y
% PC
% bin_edges
% scale
% output
% bin_centers: center of each distance bin (um)
% sim_by_bin: mean correlation within bin

if nargin < 5
    scale = 3.125;
end

% conversion
X = [x(:) * scale, y(:) * scale];
n = size(X, 1);
PC = PC(:);

% 거리 행렬 계산
dist_mat = squareform(pdist(X));
sim_mat = zscore(PC)*zscore(PC)'; %1 - squareform(pdist(PC, 'correlation')); % similarity (1 - squareform of correlation distance)

% 자기 자신은 제외
sim_mat(eye(n) == 1) = NaN;

% 거리 bin 계산
bin_centers = 0.5 * (bin_edges(1:end-1) + bin_edges(2:end));
sim_by_bin = nan(size(bin_centers));
bin_counts = zeros(size(bin_centers));

for i = 1:length(bin_centers)
    mask = dist_mat >= bin_edges(i) & dist_mat < bin_edges(i+1);
    values = sim_mat(mask);
    if ~isempty(values) && any(~isnan(values))
        sim_by_bin(i) = mean(values, 'omitnan');
        bin_counts(i) = sum(~isnan(values));
    end
end

end