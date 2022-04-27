function [] = SS_H_prediction_fly()
%% SS_H_prediction_fly.:
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
clearvars -except root FILE ALL gamma
clc

n_gamma = length(gamma);

metric = ["gain", "phase"];
n_metric = length(metric);

Data = cell(n_gamma,1);
for g = 1:n_gamma
    n_fly = ALL{g}.N.fly;
    %Data{g} = cell(n_fly,1);
    for f = 1:n_fly
        for m = 1:n_metric
            base = ALL{g}.FLY.base.H(f).(metric(m));
            pred = ALL{g}.FLY.base.H_prediction(f).(metric(m));
            learn = ALL{g}.FLY.learn.H(f).(metric(m));
            
            base = base(~isnan(base));
            pred = pred(~isnan(pred));
            learn = learn(~isnan(learn));
            
            comb = [base ; pred ; learn];

            name = ['H_' char(metric(m))];
            Data{g}(f).(name) = comb;
            
            name = ['P_' char(metric(m))];
            [~,P] = ttest2(pred, learn);
            Data{g}(f).(name) = P;
            
            if m == n_metric
                G = [1*ones(size(base)) ; 2*ones(size(pred)) ; 3*ones(size(learn))];
                Data{g}(f).G = categorical(G, [1 2 3], {'natural', 'predicted', 'experiemental'});
            end
        end
    end
end

%% Plot
close all

cc.base = [0.7 0 0];
cc.learn = [0 0.3 1];
cc.relearn = [0.5 0 1];
cc_all = [cc.base ; cc.learn ; cc.relearn];

fig = gobjects(n_gamma,1);
ax = cell(n_gamma,1);
for g = 1:n_gamma
    fig(g) = figure;
    n_fly = ALL{g}.N.fly;
    ax{g} = gobjects(n_metric,n_fly);
    p = 1;
    for m = 1:n_metric
        for f = 1:n_fly
            ax{g}(m,f) = subplot(n_metric,n_fly,p); hold on ; ylabel(metric(m))
            
            name = ['P_' char(metric(m))];
            P = 0.001*round(Data{g}(f).(name) / 0.001);
            title(P)
            
            name = ['H_' char(metric(m))];
            bx = boxplot(Data{g}(f).(name), Data{g}(f).G, ...
                'Width', 0.6, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 8);
            
            h = get(bx(5,:),{'XData','YData'});
            for c = 1:size(h,1)
               patch(h{c,1}, h{c,2}, cc_all(c,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            end
            
            set(findobj(ax{g}(m,f),'tag','Median'), 'Color', 'w','LineWidth', 1);
            set(findobj(ax{g}(m,f),'tag','Box'), 'Color', 'none');
            set(findobj(ax{g}(m,f),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
            set(findobj(ax{g}(m,f),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
            ax{g}(m,f).Children = ax{g}(m,f).Children([end 1:end-1]); 
            
            p = p + 1;
        end
    end
    set(fig(g), 'Color', 'w', 'Units', 'inches', 'Position', 1*[4 4 n_fly*1.6 n_metric*2], ...
        'Name', ['gamma = ' num2str(gamma(g))])
    set(ax{g}, 'Color', 'none', 'LineWidth', 0.75, 'Box', 'off', 'FontSize', 8, ...
        'XTickLabelRotation', 45)
    linkaxes(ax{g}, 'x')
    linkaxes(ax{g}(1,:), 'y')
    linkaxes(ax{g}(2,:), 'y')
    set(ax{g}(:,2:end), 'YColor', 'none')
    set(ax{g}(1:end-1,:), 'XColor', 'none')
    movegui(fig(g), 'center')
end

end