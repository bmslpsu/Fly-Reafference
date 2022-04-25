function [] = SS_learning_rate()
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

%% Make table
clearvars -except ALL gamma
clc

n_gamma = length(gamma);

period = ["base", "learn", "relearn"];
n_period = length(period);

metric = ["gain", "phase"];
n_metric = length(metric);

N_fly = cellfun(@(x) x.N.fly, ALL);
n_all = sum(n_gamma*n_period*mean(N_fly));
T = table(nan(n_all,1), categorical(repmat(period(1),[n_all,1])), nan(n_all,1), nan(n_all,1), ...
    'VariableNames',{'gamma','period', 'fly', 'slope'});

Slope = [];
for m = 1:n_metric
    Slope.(metric(m)) = T;
    if strcmp(metric(m), 'phase')
        scale = 180/pi;
    else
        scale = 1;
    end
    p = 1;
    for g = 1:n_gamma
        for n = 1:n_period
            for f = 1:N_fly(g)
            	fly_offset = sum(N_fly(1:g-1));
                Slope.(metric(m)).gamma(p) = gamma(g);
                Slope.(metric(m)).fly(p) = f + 0*fly_offset;
                Slope.(metric(m)).period(p) = categorical(period(n));
                Slope.(metric(m)).slope(p) = scale*60*ALL{g}.FLY.(period(n)).H(f).fit_line.(metric(m)).slope;
                p = p + 1;
            end
        end
    end
end

%% Figure
clc

cc.base = [0.7 0 0];
cc.learn = [0 0.3 1];
cc.relearn = cc.base;
cc.fly = (hsv(5));
cc_all = [cc.base ; cc.learn ; cc.relearn];
cc.gamma = parula(n_gamma);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 1.5 n_metric*1.5])
clear ax h
ax = gobjects(n_metric,1);
for m = 1:n_metric
    ax(m) = subplot(n_metric,1,m); cla ; hold on ; ylabel([char(metric(m)) ' learning rate ()'])

%     bx = boxplot(Slope.(metric(m)).slope, Slope.(metric(m)).period, ...
%         'Width', 0.6, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 8);
% 
%     h = get(bx(5,:),{'XData','YData'});
%     for c = 1:size(h,1)
%        patch(h{c,1}, h{c,2}, cc_all(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
%     end
% 
%     set(findobj(ax(m),'tag','Median'), 'Color', 'k','LineWidth', 1);
%     set(findobj(ax(m),'tag','Box'), 'Color', 'none');
%     set(findobj(ax(m),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
%     set(findobj(ax(m),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
%     ax(m).Children = ax(m).Children([end 1:end-1]);
    
%     % Plot data points
%     period_vector = findgroups(Slope.(metric(m)).period);
%     gamma_vector = findgroups(Slope.(metric(m)).gamma);
%     for p = 1:n_all
%        plot(period_vector(p), Slope.(metric(m)).slope(p), '.', ...
%            'Color', cc.gamma(gamma_vector(p),:),...
%            'MarkerFaceColor', 'none', 'MarkerSize', 8)
%     end
    
    % Plot lines connecting periods
    for g = 1:n_gamma
       for f = 1:N_fly(g)
           fly_gamma_vector = (Slope.(metric(m)).gamma == gamma(g)) ...
               & (Slope.(metric(m)).fly == f);
           fly_table = Slope.(metric(m))(fly_gamma_vector,:);
           plot(1:n_period, fly_table.slope, '.-', 'Color', cc.gamma(g,:), 'LineWidth', 0.5, ...
               'MarkerFaceColor', 'none', 'MarkerSize', 10)
       end        
    end
    
    ax(m).YLim = 1.1*max(abs(ax(m).YLim))*[-1 1];
    %plot(G_keep + jitter, time_const_keep, '.', 'Color', [0.5 0.5 0.5 0.2], 'MarkerSize', 2)
end
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'Box', 'off', 'XTickLabelRotation', 45)
set(ax, 'XLim', [0.5 n_period+0.5])
set(ax(1), 'YTick', -0.5:0.1:0.5)
set(ax(2), 'YTick', -20:10:20)

%% Stats
P = [];
stats = [];
for m = 1:n_metric
    G = {Slope.(metric(m)).period, Slope.(metric(m)).gamma};
    [P.(metric(m)),~,stats.(metric(m))] = anovan(Slope.(metric(m)).slope, G);
    %[P.(metric(m)),~,stats.(metric(m))] = kruskalwallis(Slope.(metric(m)).slope, G{1});
end

%%
[results,~,~,gnames] = multcompare(stats.phase, 'Alpha', 0.05, "Dimension",[1]);


end