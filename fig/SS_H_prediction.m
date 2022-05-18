function [] = SS_H_prediction()
%% SS_H_prediction.:
root = 'E:\DATA\Reafferent\processed';
[FILE,PATH] = uigetfile({'*.mat'},'Select data', root, 'MultiSelect', 'off');
FILE = string(FILE);

load(fullfile(PATH,FILE), 'Data', 'gamma');

%% Make table
clearvars -except Data gamma N_fly
clc

n_gamma = length(gamma);
T_pred = cell(n_gamma, 1);
for g = 1:n_gamma
    T_gamma = Data(Data.gamma==gamma(g),:);
    T_base = T_gamma(T_gamma.period == categorical("base"),:);
    T_learn = T_gamma(T_gamma.period == categorical("learn"),:);

    Hlearn_gain = T_learn.H_gain;
    Hlearn_phase = T_learn.H_phase;
    Hlearn_compensation_error = T_learn.H_compensation_error;
    
    Glearn_gain = T_learn.G_gain;
    Glearn_phase = T_learn.G_phase;

    T_comb = [T_base(:,[1,3:11]), table(Hlearn_gain, Hlearn_phase, Hlearn_compensation_error, ...
                                        Glearn_gain, Glearn_phase)];
    T_pred{g} = T_comb;
end

%% Figure
metric = ["gain", "phase", "compensation_error"];
n_metric = length(metric);

clss = ["Hpred", "Hlearn"];
n_clss = length(clss);

showbox = true;

cc.base = [0.7 0 0];
cc.learn = [0 0.3 1];
cc.relearn = [0.5 0 1];
cc.fly = (hsv(5));
% cc_all = [cc.base ; cc.learn ; cc.relearn];
cc_all = [cc.learn ; cc.relearn];
% cc.gamma = parula(n_gamma);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 n_gamma*1.2 n_metric*1.2])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_metric,n_gamma);

P = cell(n_gamma,1);
k = 1;
for m = 1:n_metric
    for g = 1:n_gamma
        ax(m,g) = subplot(n_metric,n_gamma,k); cla ; hold on        
        data = nan(size(T_pred{g},1),n_clss);
        G = nan(size(T_pred{g},1),n_clss);
        for n = 1:n_clss
           name = [char(clss(n)) '_' char(metric(m))];
           data(:,n) = T_pred{g}.(name);
           G(:,n) = n*ones(size(data(:,n)));
        end
        [~,P{g}.(metric(m))] = ttest(data(:,1), data(:,2));
        %[P{g}.(metric(m))] = ranksum(data(:,2), data(:,3));
        
      	title([num2str(gamma(g)) ', p = ' num2str(P{g}.(metric(m)))])
        
        if g == 1
            ylabel(metric(m), 'Interpreter', 'none')
        end
        
        % Boxplot
        if showbox
            bx = boxplot(data, ...
                'Width', 0.6, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 8);

            h = get(bx(5,:),{'XData','YData'});
            for c = 1:size(h,1)
               patch(h{c,1}, h{c,2}, cc_all(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            end

            set(findobj(ax(m,g),'tag','Median'), 'Color', 'w','LineWidth', 1);
            set(findobj(ax(m,g),'tag','Box'), 'Color', 'none');
            set(findobj(ax(m,g),'tag','Upper Whisker'), 'Color', 'k', 'LineStyle', '-');
            set(findobj(ax(m,g),'tag','Lower Whisker'), 'Color', 'k', 'LineStyle', '-');
            ax(m,g).Children = ax(m,g).Children([end 1:end-1]);
        end
        
        plot(G', data', '.-', 'Color', [0.6*[1 1 1] 1], 'LineWidth', 0.5, ...
               'MarkerFaceColor', 'none', 'MarkerSize', 7)
           
        % Set Y-limits
        if strcmp(metric(m), 'gain')
            ax(m,g).YLim(1) = 0;
            %yticks(0:0.5:5)
        elseif strcmp(metric(m), 'phase')
            ax(m,g).YLim(2) = 0;
            %ax(m,g).YLim(1) = -120;
            %yticks(-140:10:0)
        elseif strcmp(metric(m), 'compensation_error')
            ax(m,g).YLim = [0 1.1];            
        else
            ax(m).YLim = 1.05*max(abs(ax(m).YLim))*[-1 1];
        end
        
        k = k + 1;
    end
end
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8, 'Box', 'off', 'XTickLabelRotation', 45)
set(ax, 'XColor', 'none')

% set(ax(1,1), 'YLim', [0 1.5], 'YTick', 0:0.5:1.5)
% set(ax(1,2), 'YLim', [0 2], 'YTick', 0:0.5:2)
% set(ax(1,3), 'YLim', [0 4], 'YTick', 0:1:4)
% set(ax(1,4), 'YLim', [0 5], 'YTick', 0:1:5)

set(ax(1,1), 'YLim', [0 1], 'YTick', 0:0.5:1)
set(ax(1,2), 'YLim', [0 2], 'YTick', 0:0.5:2)
set(ax(1,3), 'YLim', [0 4], 'YTick', 0:1:4)
set(ax(1,4), 'YLim', [0 5], 'YTick', 0:1:5)

set(ax(2,1), 'YLim', [-30 0], 'YTick', -30:10:0)
set(ax(2,2), 'YLim', [-40 0], 'YTick', -40:10:0)
set(ax(2,3), 'YLim', [-60 0], 'YTick', -60:10:0)
set(ax(2,4), 'YLim', [-120 0], 'YTick', -120:20:0)

end