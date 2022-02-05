function [] = SS_gain_phase()
%% SS_gain_phase:
% root = 'E:\DATA\Reafferent';
root = 'Q:\OneDrive - PSU\OneDrive - The Pennsylvania State University\Research\Data';
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

%% Get data
clc
clearvars -except ALL gamma
clss = ["base", "learn", "relearn"];
switch_clss = ["learn", "base", "learn"];
TimePeriods = [30 905 305]; % times for [baseline, learn, relearn]

n_clss = length(clss);
n_clm = n_clss + 1;
n_gamma = length(ALL);

metric = ["gain", "phase", "compensation_error"];
n_metric = length(metric);

for g = 1:n_gamma
    ALL{g}.stats_data = [];
    for m = 1:n_metric
        for f = 1:ALL{g}.N.fly
            ALL{g}.stats_data(f).(metric(m)) = [];
            ALL{g}.stats_data(f).G = [];

            % Main data from baseline, learning, realearning
            for n = 1:n_clss
                main_data = ALL{g}.FLY.(clss(n)).H(f).(metric(m));
                nanI = ~isnan(main_data);
                main_data = main_data(nanI);
                ALL{g}.stats_data(f).(metric(m)) = [ALL{g}.stats_data(f).(metric(m)) ; main_data];
                ALL{g}.stats_data(f).G = [ALL{g}.stats_data(f).G ; n*ones(size(main_data))];
            end

            % Prediction data
            for n = 1:n_clss
                pred_data = ALL{g}.FLY.(clss(n)).H_prediction(f).(metric(m));
                nanI = ~isnan(pred_data);
                pred_data = pred_data(nanI);
                ALL{g}.stats_data(f).(metric(m)) = [ALL{g}.stats_data(f).(metric(m)) ; pred_data];
                ALL{g}.stats_data(f).G = [ALL{g}.stats_data(f).G ; (n_clss + n)*ones(size(pred_data))];
            end

            % Stats
            [~,~,stats] = kruskalwallis(ALL{g}.stats_data(f).(metric(m)), ALL{g}.stats_data(f).G, 'off');
            P = multcompare(stats, 'CType', 'bonferroni', 'Display', 'off');
            P = P(:,[1:2,end]);
            ALL{g}.stats_data(f).P.(metric(m)) = P;
        end
    end
end

%% Plot
gI = 1; % gamma #
fI = 1; % fly #

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1.5*[2 2 9 4])
movegui(fig, 'center')
clear ax h
cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [1 0.6 1];

ax = gobjects(n_metric, n_clm);
for m = 1:n_metric
    % Time
    for n = 1:n_clss
        subI = n_clm*(m-1) + n;
        ax(m,n) = subplot(n_metric,n_clm,subI); cla ; hold on ; ylabel(metric(m), 'Interpreter', 'none')
            % Prediction
            plotdata = nanmean(ALL{gI}.FLY.(switch_clss(n)).H_prediction(fI).(metric(m))) * ones(1,2);
            errdata = nanstd(ALL{gI}.FLY.(switch_clss(n)).H_prediction(fI).(metric(m))) * ones(1,2);
            if strcmp(metric(m), 'phase')
                plotdata = rad2deg(plotdata);
                errdata = rad2deg(errdata);
            end
            h = shadedErrorBar([0 TimePeriods(n)], plotdata, errdata, ...
                'lineProps', {'Color', 0.7*cc.(switch_clss(n)), 'LineWidth', 0.1});
        
            % Raw data
            plotdata = ALL{gI}.FLY.(clss(n)).H(fI).(metric(m));
            if strcmp(metric(m), 'phase')
                plotdata = rad2deg(plotdata);
            end
            plot(ALL{gI}.FLY.(clss(n)).H(fI).time, plotdata, '.', ...
                'Color', cc.(clss(n)), 'MarkerSize', 9)
            
            % Fit line
            plotdata = ALL{gI}.FLY.(clss(n)).H(fI).fit_line.(metric(m)).fit;
            if strcmp(metric(m), 'phase')
                plotdata = rad2deg(plotdata);
            end
            plot(ALL{gI}.FLY.(clss(n)).H(fI).fit_line.(metric(m)).x, plotdata, '--k', 'LineWidth', 1)
            
            % Slope & p-value
            slope = num2str(ALL{gI}.FLY.(clss(n)).H(fI).fit_line.(metric(m)).slope);
            pval = num2str(ALL{gI}.FLY.(clss(n)).H(fI).fit_line.(metric(m)).P);
            tstring = ['slope = ' num2str(slope) ', p = ' num2str(pval)];
            title(tstring)
  
        ax(m,n).XLim(1) = -0.0*TimePeriods(n);
    end
    
    % Boxplot
    subI = subI + 1;
    ax(m,n_clm) = subplot(n_metric,n_clm,subI); cla ; hold on ; ylabel(metric(m), 'Interpreter', 'none')
        plotdata = ALL{g}.stats_data(fI).(metric(m));
        if strcmp(metric(m), 'phase')
            plotdata = rad2deg(plotdata);
        end
        
        % Plot
        bx = boxplot(plotdata, ALL{g}.stats_data(fI).G, ...
            'Width', 0.8, 'Symbol', '', 'OutlierSize', 0.5);
        h = get(bx(5,:),{'XData','YData'});
        for c = 1:3
           patch(h{c,1},h{c,2}, cc.(clss(c)), 'EdgeColor', 'none', 'FaceAlpha', 1);
           patch(h{c+3,1},h{c+3,2}, 0.7*cc.(switch_clss(c)), 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
        bax = ax(m,n_clm);
        set(findobj(bax,'tag','Median'), 'Color', 'w','LineWidth', 0.75);
        set(findobj(bax,'tag','Box'), 'Color', 'none');
        set(findobj(bax,'tag','Upper Whisker'), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        set(findobj(bax,'tag','Lower Whisker'), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
        bax.Children = bax.Children([end 1:end-1]);
end

set(ax, 'Color', 'none', 'LineWidth', 1, 'Box', 'off')
set(ax(:,1:end-1), 'XGrid', 'on', 'YGrid', 'on')

set(ax(1,:), 'YLim', [0 1.5])
set(ax(2,:), 'YLim', [-40 10])
set(ax(3,:), 'YLim', [0 1])
% set(ax(4,:), 'YLim', [0 1])

for m = 1:n_metric
   linkaxes(ax(m,:), 'y') 
end

set(ax(:,2:end-0), 'YColor', 'none')
set(ax(1:end-1,:), 'XColor', 'none')
set(ax(end,end), 'XColor', 'none')

end