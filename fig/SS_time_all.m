function [] = SS_time_all()
%% SS_time_all.:
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

%% Body, error full time
clearvars -except ALL gamma
clc

gI = 1;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1.5*[2 2 5 5])
clear ax h
cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [0 0.6 0.1];

clss = ["base", "learn", "relearn"];
n_clss = length(clss);
time_space = 30;
TimePeriods = [110 810 110];
ax(1,1) = subplot(5,1,1); cla ; hold on ; ylabel('Body (°)')
    for n = 1:n_clss
        if n > 1
            window_shift = sum(TimePeriods(1:n-1));
        else
            window_shift = 0;
        end
        time_shift = ALL{gI}.FLY.(clss(n)).time(:,1) + time_space*(n-1) + window_shift;
        h.body(1,n) = plot(time_shift, nanmean(ALL{gI}.FLY.(clss(n)).body,2), 'Color', cc.(clss(n)));
        
        A = max(ALL{gI}.FLY.(clss(n)).function, [], 'all');
        plot(time_shift,  A*ones(size(time_shift)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        plot(time_shift, -A*ones(size(time_shift)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    end
    
ax(2,1) = subplot(5,1,2); cla ; hold on ; ylabel('Error (°)')
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
    
ax(3,1) = subplot(5,1,3); cla ; hold on ; ylabel('Gain')
    for n = 1:n_clss
        if n > 1
            window_shift = sum(TimePeriods(1:n-1));
        else
            window_shift = 0;
        end
        time_shift = ALL{gI}.FLY.(clss(n)).H(1).time + time_space*(n-1) + window_shift;
        Hgain_all = [ALL{gI}.FLY.(clss(n)).H(:).gain];
        h.gain(:,n) = plot(time_shift, Hgain_all, '.', ...
            'Color', cc.(clss(n)), 'MarkerSize', 8, 'markerFaceColor', 'none');
        
        for f = 1:ALL{gI}.N.fly
            time_shift = ALL{gI}.FLY.(clss(n)).H(f).fit_line.gain.x + time_space*(n-1) + window_shift;
            plot(time_shift, ALL{gI}.FLY.(clss(n)).H(f).fit_line.gain.fit, 'k', 'LineWidth', 1)
        end
    end

ax(4,1) = subplot(5,1,4); cla ; hold on ; ylabel('Phase (°)')
    for n = 1:n_clss
        if n > 1
            window_shift = sum(TimePeriods(1:n-1));
        else
            window_shift = 0;
        end
        time_shift = ALL{gI}.FLY.(clss(n)).H(1).time + time_space*(n-1) + window_shift;
        Hphase_all = rad2deg([ALL{gI}.FLY.(clss(n)).H(:).phase]);
        h.phase(:,n) = plot(time_shift, Hphase_all, '.', ...
            'Color', cc.(clss(n)), 'MarkerSize', 8, 'markerFaceColor', 'none');
        
        for f = 1:ALL{gI}.N.fly
            time_shift = ALL{gI}.FLY.(clss(n)).H(f).fit_line.phase.x + time_space*(n-1) + window_shift;
            plot(time_shift, rad2deg(ALL{gI}.FLY.(clss(n)).H(f).fit_line.phase.fit), 'k', 'LineWidth', 1)
        end
    end
    
ax(5,1) = subplot(5,1,5); cla ; hold on ; ylabel('R2')
    for n = 1:n_clss
        if n > 1
            window_shift = sum(TimePeriods(1:n-1));
        else
            window_shift = 0;
        end
        time_shift = ALL{gI}.FLY.(clss(n)).H(1).time + time_space*(n-1) + window_shift;
        HR2_all = [ALL{gI}.FLY.(clss(n)).H(:).R2];
        h.R2(:,n) = plot(time_shift, HR2_all, '.', ...
            'Color', cc.(clss(n)), 'MarkerSize', 8, 'markerFaceColor', 'none');
    end
    
set(ax, 'Color', 'none', 'LineWidth', 1)
linkaxes(ax, 'x')
linkaxes(ax(1:2,1), 'y')

set([h.body , h.error], 'LineWidth', 0.5)

set(ax, 'XLim', [-20 1100], 'XTick', 0:100:1300)
set(ax(1:2), 'YLim', 100*[-1 1], 'YTick', -100:25:100)
set(ax(3), 'YLim', [0 max([h.gain(:).YData])])
set(ax(4), 'YLim', [min([h.phase(:).YData]) max([h.phase(:).YData])])
set(ax(5), 'YLim', [0 1])

% set(ax(:,2), 'YColor', 'none')
set(ax(1,1), 'XColor', 'none')

end