function [] = SS_time()
%% SS_time:
root = 'E:\DATA\Reafferent';
[FILE,PATH] = uigetfile({'*.mat'},'Select data', root, 'MultiSelect', 'on');
FILE = string(FILE);
n_file = length(FILE);

gamma = nan(n_file,1);
ALL = cell(n_file,1);
for n = 1:length(FILE)
    finfo = strsplit(FILE(n), '_');
    gamma(n) = str2double(finfo(3));
    ALL{n} = load(fullfile(PATH,FILE(n)),'FLY','DATA','U','N');
end

%% Body, error full time
clearvars -except ALL gamma
gI = 1;

clss = ["base", "learn", "relearn"];
TimePeriods = [30 905 305]; % times for [baseline, learn, relearn]
n_clss = length(clss);
time_space = 30;

bin_size = 1;
bin_range = 50;
bins = -(bin_range + bin_size/2):bin_size:(bin_range + bin_size/2);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1.5*[2 2 5 2])
clear ax h
cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [1 0.6 1];

% func = 37.5*sin(2*pi*0.5*tt);
ax(1,1) = subplot(2,4,1:3); cla ; hold on ; ylabel('Body (°)')
    %plot(tt, func, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
    yline(37.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    yline(-37.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    for n = 1:n_clss
        if n > 1
            window_shift = sum(TimePeriods(1:n-1));
        else
            window_shift = 0;
        end
        time_shift = ALL{gI}.FLY.(clss(n)).time(:,1) + time_space*(n-1) + window_shift;
        xx = [time_shift(1) time_shift(end) time_shift(end) time_shift(1)];
        yy = 100*[-1 -1 1 1];
        patch(xx, yy, cc.(clss(n)), 'FaceAlpha', 0.1, 'EdgeColor', 'none')
        
        h.body(1,n) = plot(time_shift, nanmean(ALL{gI}.FLY.(clss(n)).body,2), 'Color', cc.(clss(n)));
    end
    
ax(1,2) = subplot(2,4,4); cla ; hold on
    for n = 1:n_clss
        h.hist(1,n) = histogram(ALL{gI}.FLY.(clss(n)).body, bins, ...
            'FaceColor', cc.(clss(n)), 'EdgeColor', cc.(clss(n)));
    end
    
ax(2,1) = subplot(2,4,5:7); cla ; hold on ; ylabel('Error (°)')
    for n = 1:n_clss
        if n > 1
            window_shift = sum(TimePeriods(1:n-1));
        else
            window_shift = 0;
        end
        time_shift = ALL{gI}.FLY.(clss(n)).time(:,1) + time_space*(n-1) + window_shift;
        h.error(1,n) = plot(time_shift, nanmean(ALL{gI}.FLY.(clss(n)).error,2), 'Color', cc.(clss(n)));
    end
    xlabel('time (s)')
    
ax(2,2) = subplot(2,4,8); cla ; hold on
    for n = 1:n_clss
        h.hist(2,n) = histogram(ALL{gI}.FLY.(clss(n)).error, bins, ...
            'FaceColor', cc.(clss(n)), 'EdgeColor', cc.(clss(n)));
    end
    
set(ax, 'Color', 'none', 'LineWidth', 1)
set(h.hist, 'Normalization', 'probability', 'EdgeColor', 'none', ...
    'Orientation', 'horizontal', 'DisplayStyle', 'bar', 'FaceAlpha', 0.4)

linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
linkaxes(ax(:,1), 'xy')

set([h.body , h.error], 'LineWidth', 0.5)

set(ax(:,1), 'XLim', [-20 1300], 'XTick', 0:100:1300)
set(ax, 'YLim', 55*[-1 1], 'YTick', -50:10:50)
% ax(1,2).XLim(1) = -0.05*ax(1,2).XLim(2);
% ax(2,2).XLim(1) = -0.05*ax(2,2).XLim(2);

set(ax(:,2), 'YColor', 'none')
set(ax(1,1), 'XColor', 'none')

end