
%% Plot latency
clear ; close all ; clc

root = 'S:\Public\BC\Mario Data';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select files', root, 'MultiSelect', 'on');
FILE = string(FILE);
n_file = length(FILE);

DATA = cell(n_file,1);
for n = 1:n_file
    fpath = fullfile(PATH, FILE(n));
    matdata = load(fpath, 'latency');
    DATA{n} = matdata.latency;
end

ALL = cat(1, DATA{:});
ALL = sort(ALL, 'descend');
badI = ALL > 0.2;
disp(['Bad = ' num2str(sum(badI))])
ALL = ALL(~badI);
med = 1000*median(ALL);

%% Plot
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 2*[2 2 2 1.2])
movegui(fig, 'center')
clear ax h

cc = jet(n_file);

bin_size = 0.05;
bin_range = 20;
bins = 0:bin_size:bin_range;

ax = subplot(1,1,1); cla ; hold on
h = gobjects(n_file+1,1);
for n = 1:n_file
    h(n) = histogram(1000*DATA{n}, bins, 'FaceColor', cc(n,:));
end
h(n_file+1) = histogram(1000*ALL, bins, 'FaceColor', 'k');
xline(med, '--r');
title(['median = ' num2str(med) ' ms'])

set(h, 'Normalization', 'probability', 'EdgeColor', 'none')

delete(h(1:end-1))

set(ax, 'Color', 'none', 'LineWidth', 0.75,...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')

xlim([4 10])
xlabel('latency (ms)')
ylabel('probability')

ax.YLim(1) = -0.05*ax.YLim(2);



