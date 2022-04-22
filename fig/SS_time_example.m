function [] = SS_time_example()
%% SS_time_example.:
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

clss = ["base", "learn", "relearn"];
t_length = 6; % [s] time length
flyI = [1 2 3 1]; % fly number index
t_start = [30 14  10  10;...  % [s] start times
           30 20 104  6;
           8  0  18  0];

n_clss = length(clss);
n_gamma = length(gamma);

cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [0 0.6 0.1];
cc.error = [112 96 167]./255;

lw = 1;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 n_gamma*2.5 n_clss*1.5])
clear ax h
subI = 1;
ax = gobjects(n_clss, n_gamma);
A = nan(n_clss, n_gamma);
for n = 1:n_clss
    for g = 1:n_gamma
        ax(n,g) = subplot(n_clss, n_gamma, subI); cla ; hold on
        if (n == 1)
            title(gamma(g))
        end
        if (g == 1)
           ylabel(clss(n)) 
        end
        
        A(n,g) = max(ALL{g}.FLY.(clss(n)).function(:,1));
        
        time = ALL{g}.FLY.(clss(n)).time(:,1);
        Fs = 1 ./ mean(diff(time));
        winI = round(Fs*t_start(n,g) + (round(0*Fs):round(t_length*Fs)) + 1);
        win_time = time(winI);
        win_time = win_time - win_time(1);
        
        plot(win_time, ALL{g}.FLY.(clss(n)).function(winI,1), ...
            'Color', [0.5 0.5 0.5], 'LineWidth', lw)
        plot(win_time, ALL{g}.FLY.(clss(n)).body(winI,flyI(g)), ...
            'Color', cc.(clss(n)),  'LineWidth', lw)
        plot(win_time, ALL{g}.FLY.(clss(n)).error(winI,flyI(g)), '--', ...
            'Color', cc.error,  'LineWidth', lw)        
        
        ylim(5*A(n,g)*[-1 1])
        yticks([0 A(n,g)])
        ax(n,g).XLim(1) = -0.2;
        subI = subI + 1;
    end
end
linkaxes(ax, 'x')
set(ax, 'Color', 'none', 'LineWidth', 1)
set(ax(1:end-1,:), 'XColor', 'none')

end