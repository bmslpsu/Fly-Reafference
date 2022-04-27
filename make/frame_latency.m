function [] = frame_latency()
%% frame_latency.:
root = 'E:\DATA\Reafferent';
[FILE,PATH] = uigetfile({'*.mat'},'Select data', root, 'MultiSelect', 'on');
FILE = string(FILE);
n_file = length(FILE);

gamma = nan(n_file,1);
ALL = cell(n_file,1);
for n = 1:length(FILE)
    finfo = strsplit(FILE(n), '_');
    gamma(n) = str2double(finfo(3));
    ALL{n} = load(fullfile(PATH,FILE(n)),'FLY','U','N');
end
[gamma, gI] = sort(gamma);
ALL = ALL(gI);

%% Get latency
clearvars -except root FILE ALL gamma
clc

n_gamma = length(gamma);
Latency = cell(n_gamma,1);
for g = 1:n_gamma
   Latency{g} = 1000*ALL{g}.FLY.base.camera_times;
end
Latency_all = cat(2, Latency{:});
med = median(Latency_all(:));

%% Plot
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 2 2])
movegui(fig, 'center')
clear ax h

bin_size = 0.05;
bin_range = 20;
bins = 5:bin_size:bin_range;

ax = subplot(1,1,1); cla ; hold on ; title(['median = ' num2str(med) ' ms'])
h = histogram(Latency_all, bins, 'FaceColor', 'k');
xline(med, '--r');

set(h, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.4)
set(ax, 'Color', 'none', 'LineWidth', 0.75,...
    'XGrid', 'off', 'YGrid', 'off', 'Box', 'off')

% xlim([4 10])
xlabel('latency (ms)')
ylabel('probability')

ax.YLim(1) = -0.05*ax.YLim(2);


end