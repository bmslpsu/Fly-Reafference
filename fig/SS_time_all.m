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

gI = 4;

cc.base = [0.5 0.5 0.5];
% cc.learn = [0 0.7 1];
cc.learn = cc.base;
% cc.relearn = [1 0.4 0.6];
cc.relearn = cc.base;
cc.fly = (hsv(5));

mkrsz = 7;

% Patch
xx = [0, 1100, 1100, 0];
n_std = 1;
facealpha = 0;

clss = ["base", "learn", "relearn"];
n_clss = length(clss);
time_space = 30;
TimePeriods = [110 810 110];
window_shift = nan(1,n_clss);
for n = 1:n_clss
    if n > 1
        window_shift(n) = sum(TimePeriods(1:n-1));
    else
        window_shift(n) = 0;
    end
end

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 3 6])
clear ax h

ax(1,1) = subplot(5,1,1); cla ; hold on ; ylabel('Body (°)') ; title(gamma(gI))
    for n = 1:n_clss
        time_shift = ALL{gI}.FLY.(clss(n)).time(:,1) + time_space*(n-1) + window_shift(n);
        h.body(1,n) = plot(time_shift, nanmean(ALL{gI}.FLY.(clss(n)).body,2), 'Color', cc.(clss(n)));
        
        A = max(ALL{gI}.FLY.(clss(n)).function, [], 'all');
        plot(time_shift,  A*ones(size(time_shift)), '--', 'Color', [1 0.5 0.5], 'LineWidth', 2)
        plot(time_shift, -A*ones(size(time_shift)), '--', 'Color', [1 0.5 0.5], 'LineWidth', 2)
        
        win_xx = [time_shift(1) time_shift(end) time_shift(end) time_shift(1)];
        win_yy = 100*[-1 -1 1 1];
        patch(win_xx, win_yy, cc.(clss(n)), 'FaceAlpha', 0.2, 'LineStyle', 'none')
    end
    
ax(2,1) = subplot(5,1,2); cla ; hold on ; ylabel('Error (°)')
    for n = 1:n_clss
        time_shift = ALL{gI}.FLY.(clss(n)).time(:,1) + time_space*(n-1) + window_shift(n);
        h.error(1,n) = plot(time_shift, nanmean(ALL{gI}.FLY.(clss(n)).error,2), 'Color', cc.(clss(n)));
    end
    xlabel('time (s)')
    
ax(3,1) = subplot(5,1,3); cla ; hold on ; ylabel('Gain')
    for n = 1:n_clss
        time_shift = ALL{gI}.FLY.(clss(n)).H(1).time + time_space*(n-1) + window_shift(n);
        Hgain_all = [ALL{gI}.FLY.(clss(n)).H(:).gain];
        Hgain_mean = nanmean(Hgain_all(:));
        Hgain_std = nanstd(Hgain_all(:));
        
        yy = [Hgain_mean-n_std*Hgain_std, Hgain_mean-n_std*Hgain_std, ...
            Hgain_mean+n_std*Hgain_std, Hgain_mean+n_std*Hgain_std];
        patch(xx, yy, cc.(clss(n)), 'FaceAlpha', facealpha, 'LineStyle', 'none')
        h.gain(:,n) = plot(time_shift, Hgain_all, '.', ...
            'Color', cc.(clss(n)), 'MarkerSize', mkrsz, 'markerFaceColor', 'none');
        
        for f = 1:ALL{gI}.N.fly
            time_shift = ALL{gI}.FLY.(clss(n)).H(f).fit_line.gain.x + time_space*(n-1) + window_shift(n);
            plot(time_shift, ALL{gI}.FLY.(clss(n)).H(f).fit_line.gain.fit, ...
                'Color', cc.fly(f,:), 'LineWidth', 1)
        end
    end

ax(4,1) = subplot(5,1,4); cla ; hold on ; ylabel('Phase (°)')
    for n = 1:n_clss
        time_shift = ALL{gI}.FLY.(clss(n)).H(1).time + time_space*(n-1) + window_shift(n);
        Hphase_all = rad2deg([ALL{gI}.FLY.(clss(n)).H(:).phase]);
     	Hphase_mean = nanmean(Hphase_all(:));
        Hphase_std = nanstd(Hphase_all(:));
        yy = [Hphase_mean-n_std*Hphase_std, Hphase_mean-n_std*Hphase_std, ...
            Hphase_mean+n_std*Hphase_std, Hphase_mean+n_std*Hphase_std];
        patch(xx, yy, cc.(clss(n)), 'FaceAlpha', facealpha, 'LineStyle', 'none')
        h.phase(:,n) = plot(time_shift, Hphase_all, '.', ...
            'Color', cc.(clss(n)), 'MarkerSize', mkrsz, 'markerFaceColor', 'none');
        
        for f = 1:ALL{gI}.N.fly
            time_shift = ALL{gI}.FLY.(clss(n)).H(f).fit_line.phase.x + time_space*(n-1) + window_shift(n);
            plot(time_shift, rad2deg(ALL{gI}.FLY.(clss(n)).H(f).fit_line.phase.fit), ...
                'Color', cc.fly(f,:), 'LineWidth', 1)
        end
    end
    
ax(5,1) = subplot(5,1,5); cla ; hold on ; ylabel('R2')
    for n = 1:n_clss
        time_shift = ALL{gI}.FLY.(clss(n)).H(1).time + time_space*(n-1) + window_shift(n);
        HR2_all = [ALL{gI}.FLY.(clss(n)).H(:).R2];
        h.R2(:,n) = plot(time_shift, HR2_all, '.', ...
            'Color', cc.(clss(n)), 'MarkerSize', mkrsz, 'markerFaceColor', 'none');
    end
    
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8)
linkaxes(ax, 'x')
linkaxes(ax(1:2,1), 'y')

set([h.body , h.error], 'LineWidth', 0.5)

set(ax, 'XLim', [-20 1200], 'XTick', 0:100:1200)
set(ax(1:2), 'YLim', 100*[-1 1], 'YTick', -100:50:100)
set(ax(3), 'YLim', [0 max([h.gain(:).YData])], 'YTick', 0:1:10)
set(ax(4), 'YLim', [min([h.phase(:).YData]) max([h.phase(:).YData])])
set(ax(5), 'YLim', [0 1])

set(ax(3), 'YLim', [0 6])
set(ax(4), 'YLim', [-150 30])

% set(ax(:,2), 'YColor', 'none')
set(ax(1:end-1,1), 'XColor', 'none')

end