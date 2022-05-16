function [] = SS_time_raw()
%% SS_time_raw.:
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

%%
clearvars -except ALL gamma
clc

n_gamma = length(gamma);

clss = ["base", "learn"];
n_clss = length(clss);

time_end = 80;

Fs = 100;
Fc = 0.1;
[b,a] = butter(5, Fc / (Fs/2), 'low');
Data = cell(n_gamma,1);
for g = 1:n_gamma
    for n = 1:n_clss
        time = ALL{g}.FLY.(clss(n)).time(:,1);
        endI = find(time >= time_end, 1, 'first');
        time = time(1:endI,:);
        pos = cat(2, ALL{g}.FLY.(clss(n)).all.body_raw{:});
        pos = pos(1:endI,:);
        
        vel = pos;
        vel_filt = pos;
        for c = 1:size(vel_filt,2)
            vel(:,c) = central_diff(pos(:,c), 1/Fs);
            vel_filt(:,c) = filtfilt(b, a, vel(:,c));
        end
        vel_sign = sign((mean(vel_filt,1)));
        pos_norm = pos .* vel_sign;
        %vel_norm = vel .* vel_sign;
        vel_filt_norm = vel_filt .* vel_sign;
        
        for c = 1:size(vel_filt,2)
            [out] = fit_exp(time, pos_norm(:,c), false);
            Data{g}.(clss(n)).slope(1,c) = out.slope;
            Data{g}.(clss(n)).R2(1,c) = out.R2;
        end
        
        Data{g}.(clss(n)).time = time;
        Data{g}.(clss(n)).pos = pos_norm;
        Data{g}.(clss(n)).vel = vel_filt_norm;
        Data{g}.(clss(n)).mean_vel = mean((vel_filt_norm),1);
    end
end

%% Raw time
cc.base = repmat([0.5 0.5 0.5], [n_gamma, 1]);
% cc.learn = [parula(n_gamma-1) ; [1 0 0]];
cc.learn = repmat([55 156 235]./255, [n_gamma, 1]);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 n_gamma*1.6 1.5])
movegui(fig, 'center')
clear ax h
ax = gobjects(1,n_gamma);

for g = 1:n_gamma
    ax(g) = subplot(1,n_gamma,g); cla ; hold on
    title(['\gamma = ' num2str(gamma(g)) ])
        for n = 1:n_clss
            plot(Data{g}.(clss(n)).time, Data{g}.(clss(n)).pos, ...
                'Color', [cc.(clss(n))(g,:) 1], 'LineWidth', 1)
        end
    if g == 1
        ylabel('Body (°)')
        xlabel('Time (s)')
    end
end
    
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8)
set(ax, 'XLim',[0 time_end], 'XTick', 0:20:1000)

all_pos = cellfun(@(x) x.learn.pos(:), Data, 'UniformOutput', false);
all_pos = cat(1, all_pos{:});
set(ax, 'YLim', max(all_pos)*[-0.1 1], 'YTick', 10e4*[-2:0.1:3])
set(ax(2:end), 'YColor', 'none', 'XTickLabels', [])
linkaxes(ax, 'xy')

%% Mean Velocity
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 1.5 3])
movegui(fig, 'center')
clear ax h
ax = gobjects(n_clss,1);
cc_all = [0.5 0.5 0.5];
P = [];
for n = 1:n_clss
    plotdata = [];
    G = [];
    for g = 1:n_gamma
        data = Data{g}.(clss(n)).mean_vel';
        plotdata = [plotdata ;  data];
        if g < 5
            gg = 1;
        else
            gg = 2;
        end
        G = [G ; g*ones(size(data))];
    end
    [~,P.(clss(n))] = ttest2(plotdata(G==1), plotdata(G==2));
    
    ax(n) = subplot(n_clss,1,n) ; cla ; hold on ; title(clss(n))
    ylabel('Velocity ()')
        bx = boxplot(plotdata, G, ...
            'Width', 0.6, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 8);
        
        h = get(bx(5,:),{'XData','YData'});
        for c = 1:size(h,1)
           patch(h{c,1}, h{c,2}, cc_all, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end

        set(findobj(ax(n),'tag','Median'), 'Color', 'w','LineWidth', 1);
        set(findobj(ax(n),'tag','Box'), 'Color', 'none');
        set(findobj(ax(n),'tag','Upper Whisker'), 'Color', 'k', 'LineStyle', '-');
        set(findobj(ax(n),'tag','Lower Whisker'), 'Color', 'k', 'LineStyle', '-');
        ax(n).Children = ax(n).Children([end 1:end-1]);
        
        plot(findgroups(G), plotdata, '.', ...
            'Color', cc.learn(1,:), 'LineWidth', 0.5, ...
               'MarkerFaceColor', 'none', 'MarkerSize', 7)
        
end
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'FontSize', 8)
set(ax, 'Box', 'off')
set(ax(1), 'XColor', 'none')

P

end