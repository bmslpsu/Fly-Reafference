function [] = SS_gamma_stats()
%% SS_gamma_stats.:
root = 'E:\DATA\Reafferent\processed';
[FILE,PATH] = uigetfile({'*.mat'},'Select data', root, 'MultiSelect', 'off');
FILE = string(FILE);

load(fullfile(PATH,FILE), 'Data', 'gamma', 'N_fly');

%% Figure
clearvars -except Data gamma N_fly
clc

n_gamma = length(gamma);

period = string(categories(Data.period));
n_period = length(period);

metric = ["H_gain_slope", "H_phase_slope", "H_gain", "H_phase", "H_compensation_error"];
n_metric = length(metric);

showbox = false;

cc.base = [0.7 0 0];
cc.learn = [0 0.3 1];
cc.relearn = cc.base;
cc.fly = (hsv(5));
cc_all = [cc.base ; cc.learn ; cc.relearn];
cc.gamma = parula(n_gamma);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 1.5 n_metric*1.2])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_metric,1);
for m = 1:n_metric
    ax(m) = subplot(n_metric,1,m); cla ; hold on
    ylabel([char(metric(m)) ' ()'], 'Interpreter', 'none')

    % Boxplot
    if showbox
        bx = boxplot(Data.(metric(m)), te, ...
            'Width', 0.6, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 8);

        h = get(bx(5,:),{'XData','YData'});
        for c = 1:size(h,1)
           patch(h{c,1}, h{c,2}, cc_all(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end

        set(findobj(ax(m),'tag','Median'), 'Color', 'k','LineWidth', 1);
        set(findobj(ax(m),'tag','Box'), 'Color', 'none');
        set(findobj(ax(m),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
        set(findobj(ax(m),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
        ax(m).Children = ax(m).Children([end 1:end-1]);
    end
    
    % Plot data points & lines connecting periods
    for g = 1:n_gamma
       for f = 1:N_fly(g)
           fly_gamma_vector = (Data.gamma == gamma(g)) ...
               & (Data.fly == f);
           fly_table = Data(fly_gamma_vector,:);
           plot(1:n_period, fly_table.(metric(m)), '.-', 'Color', cc.gamma(g,:), 'LineWidth', 0.5, ...
               'MarkerFaceColor', 'none', 'MarkerSize', 10)
       end        
    end
    
    % Set Y-limits
    if strcmp(metric(m), 'H_gain')
        ax(m).YLim = 1.05*max(abs(ax(m).YLim))*[0 1];
        yticks(0:1:5)
    elseif strcmp(metric(m), 'H_phase')
        ax(m).YLim = 1.05*max(abs(ax(m).YLim))*[-1 0];
        yticks(-140:20:0)
    elseif strcmp(metric(m), 'H_compensation_error')
        ax(m).YLim = [0 1.1];
        yticks(0:0.2:1)
    else
        ax(m).YLim = 1.05*max(abs(ax(m).YLim))*[-1 1];
    end
    
    % Indicate periods
    win_sz = 0.3;
    for n = 1:n_period
        xx = [n-win_sz, n+win_sz, n+win_sz, n-win_sz];
        yy = [ax(m).YLim(1) ax(m).YLim(1) ax(m).YLim(2) ax(m).YLim(2)];
        hp = patch(xx, yy, cc.(period(n)), 'FaceAlpha', 0.1, 'LineStyle', 'none');
        uistack(hp, 'bottom')
    end
end
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'Box', 'off', 'XTickLabelRotation', 45)
set(ax, 'XLim', [0.5 n_period+0.5], 'XColor', 'none')
set(ax(1), 'YTick', -0.5:0.1:0.5)
set(ax(2), 'YTick', -20:10:20)

end