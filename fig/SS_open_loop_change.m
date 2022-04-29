function [] = SS_open_loop_change()
%% SS_open_loop_change.:
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
    
    L = (Glearn_gain.*exp(deg2rad(Glearn_phase)*1j)) ./ ...
        (T_base.G_gain.*exp(deg2rad(T_base.G_phase)*1j));
    L_gain = abs(L);
    L_phase = rad2deg(angle(L));

    T_comb = [T_base(:,[1,3:11]), table(Hlearn_gain, Hlearn_phase, Hlearn_compensation_error, ...
                                        Glearn_gain, Glearn_phase, L_gain, L_phase)];
    T_pred{g} = T_comb;
end

T_pred = cat(1, T_pred{:});

%% Figure
clc
metric = ["gain", "phase"];
n_metric = length(metric);

clss = ["L"];

showbox = false;

cc.gamma = parula(n_gamma);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1.2*[2 2 2 n_metric*1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_metric,1);

P = cell(n_gamma,1);
k = 1;
for m = 1:n_metric
    ax(m,1) = subplot(n_metric,1,m); cla ; hold on ; ylabel(metric(m), 'Interpreter', 'none')

    name = [char(clss) '_' char(metric(m))];
    data = T_pred.(name);
    G = T_pred.gamma;

    % Boxplot
    if showbox
        bx = boxplot(data, G, 'Positions', G, ...
            'Width', 0.2, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 6);

        h = get(bx(5,:),{'XData','YData'});
        for c = 1:size(h,1)
           patch(h{c,1}, h{c,2}, cc.gamma(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        end

        set(findobj(ax(m,1),'tag','Median'), 'Color', 'w','LineWidth', 1);
        set(findobj(ax(m,1),'tag','Box'), 'Color', 'none');
        set(findobj(ax(m,1),'tag','Upper Whisker'), 'Color', 'k', 'LineStyle', '-');
        set(findobj(ax(m,1),'tag','Lower Whisker'), 'Color', 'k', 'LineStyle', '-');
        ax(m,1).Children = ax(m,1).Children([end 1:end-1]);
    end
    
    for g = 1:n_gamma
        gI = G == gamma(g);
        med = mean(data(gI));
        h.fly(g) = plot(G(gI), data(gI), '.', 'Color', [0.6*[1 1 1] 1], 'LineWidth', 0.5, ...
               'MarkerFaceColor', 'none', 'MarkerSize', 8);
      	h.med(g) = plot(gamma(g), med, '.', 'Color', cc.gamma(g,:), ...
            'MarkerFaceColor', 'none', 'MarkerSize', 12);
    end
end
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8, 'Box', 'on', 'XTickLabelRotation', 45)
set(ax, 'XGrid', 'on', 'YGrid', 'on')
% set(ax, 'XColor', 'none')
% set(ax(1,:), 'YLim', [0 6], 'YTick', 0:1:5)
% set(ax(2,:), 'YLim', [-120 0], 'YTick', -120:20:0)

axes(ax(1))
gamma_num = (-2:0.01:4)';
Ls_num = 1 ./ (1 - gamma_num);
plot(gamma_num, Ls_num, 'k', 'LineWidth', 1)
xline(1, '--k');
yline(1, '--r');

xlim([min(gamma_num) max(gamma_num)])
xticks(-10:0.5:10)
ylim(3*[-1 1])
yticks(-10:1:10)

uistack([h.fly, h.med], 'top')
xlabel('\gamma')

end