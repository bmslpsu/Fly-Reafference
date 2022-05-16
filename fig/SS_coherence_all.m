function [] = SS_coherence_all()
%% SS_coherence_all.:
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

%% Coherence stats
clearvars -except ALL gamma
clc

n_gamma = length(gamma);

clss = ["base", "learn"];
n_clss = length(clss);

fv = (0:0.02:50)';
f_crit = 0.5;
f_critI = fv == f_crit;

% Stats = cell(n_gamma,1);
% G = cell(n_gamma,1);
% for g = 1:n_gamma
%     for n = 1:n_clss
%         flydata = cat(2, ALL{g}.FLY.(clss(n)).all.body_raw_coherence{:});
%      	Stats{g} = [Stats{g} ; flydata(f_critI,:)'];
%         G{g} = [G{g} ; n*ones(size(flydata,2))'];
%     end
% end

Stats = [];
for g = 1:n_gamma
    for n = 1:n_clss
        flydata = cat(2, ALL{g}.FLY.(clss(n)).all.body_raw_coherence{:});
     	Stats(g).(clss(n)) = flydata(f_critI,:);
    end
end

for g = 1:n_gamma 
    [~,Stats(g).P] = ttest(Stats(g).base, Stats(g).learn);
end

%% Coherence
cc.base = repmat([0.5 0.5 0.5], [n_gamma, 1]);
% cc.learn = [1 0 0];
% cc.learn = [parula(n_gamma-1) ; [1 0 0]];
cc.learn = repmat([55 156 235]./255, [n_gamma, 1]);

% cc.learn = repmat([1 0 0], [n_gamma, 1]);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 1.8 n_gamma*1.6 1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(1,n_gamma);

f_win = 0.015;
xx = [f_crit-f_win, f_crit+f_win, f_crit+f_win, f_crit-f_win];
yy = [-0.1 -0.1 1 1];
for g = 1:n_gamma
    ax(g) = subplot(1,n_gamma,g); cla ; hold on
    title(['\gamma = ' num2str(gamma(g)) ', p = ' num2str(Stats(g).P,2)])
        patch(xx, yy, [0.7 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none')
        for n = 1:n_clss
            flydata = cat(2, ALL{g}.FLY.(clss(n)).all.body_raw_coherence{:});
            plot(fv, flydata, 'Color', [cc.(clss(n))(g,:) 1], 'LineWidth', 0.75)
        end
    if g == 1
        ylabel('Coheherence')
        xlabel('Frequency (hz)')
    end
end
    
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8)
set(ax, 'XScale', 'log', 'XLim', [0.3 1], 'XTick', [0.1:0.1:0.9 1:3], 'YLim', [-0.04 1])
set(ax(2:end), 'YColor', 'none', 'XTickLabels', [])
% set(ax, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on')
linkaxes(ax, 'xy')

end