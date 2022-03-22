function [] = make_data_sine()
%% make_data_sine:
%

gamma_folder = 0.5;
fpath = 'E:\EXPERIMENTS\MAGNO\Experiment_reafferent_sine';
root = fullfile(fpath, ['gamma=' num2str(gamma_folder)]);
filename = ['SS_gamma_' num2str(gamma_folder)];

% Select files
[D,I,N,U,~,~,~,basename] = GetFileData(root,'*.mat', false, 'fly', 'trial', 'gamma');

%% Get Data
warning('off', 'signal:findpeaks:largeMinPeakHeight')
close all
clc

% Saccade detection parameters
sacd.showplot = false;
sacd.Fc_detect = [10 nan];
sacd.Fc_ss = [nan nan];
sacd.amp_cut = 25;
sacd.dur_cut = inf;
sacd.thresh = [0, 1, 4, 300];
sacd.true_thresh = 275;
sacd.sacd_length = 2;
sacd.pks = [];
sacd.min_pkdist = 0.1;
sacd.min_pkwidth = 0.03;
sacd.min_pkprom = 75;
sacd.min_pkthresh = 0;
sacd.boundThresh = [0.2 60];
sacd.direction = 0;

Fs = 90;
Fc = 3;
[b_low,a_low] = butter(5, Fc/(Fs/2), 'low');
[b_high,a_high] = butter(5, 0.3/(Fs/2), 'high');

win_size = 10;
overlap = 5;
R2_cut = 0.5;

DATA = [D(:,1:3) , splitvars(table(num2cell(zeros(N.file,15))))];
DATA.Properties.VariableNames(4:end) = {'time', 'function', 'body', 'body_raw', 'body_intrp', ...
    'error', 'display', 'function_display', 'error_display', 'head', 'body_saccade', 'camera_times', ...
    'H', 'G', 'H_prediction'};
TimePeriods = [125 905 905]; % times for [baseline, learn, relearn]
% TimePeriods = [30 905 305]; % times for [baseline, learn, relearn]
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
        
    % Get system properties
    gamma = D.gamma(n); % gamma (coupling gain between fly body and display)
    %gain = -1 + gamma; % total gain (natural + coupling)
    
    % Load DAQ, body, head, & wing data
	all_data = load(fullfile(root, [basename{n} '.mat']),...
        'data','t_p', 'bAngles', 't_v', 'ampl', 'freq'); % load data
    
    % Make time vector
    tintrp = (0:1/Fs:TimePeriods(D.trial(n)))';
    
    % Recreate the input sine function
    freq = all_data.freq;
    A = all_data.ampl;
    Phi = 0;
    func = A*sin(2*pi*freq*tintrp + Phi);
    func = func - mean(func);

    % Sync video and display data
    daq_time = all_data.t_p;
    
    display = all_data.data(:,1);
    display = 3.75*(round(96*(display ./ 10)) - 1);
    display = rad2deg(unwrap(deg2rad(display)));
    
    dx_pat = diff(display);
    syncI = find(abs(dx_pat) > 5, 1, 'first') + 1;
    sync_time = daq_time(syncI);
    
    daq_time_sync = daq_time - sync_time;
    display_intrp = interp1(daq_time_sync, display, tintrp, 'nearest');
    
    % Remove saccades from body
    vid_time = all_data.t_v;
    DATA.camera_times{n} = diff(vid_time);
    body_intrp = interp1(vid_time, all_data.bAngles, tintrp, 'pchip');
    body_scd = saccade_all(body_intrp, tintrp, sacd.thresh, sacd.true_thresh, sacd.Fc_detect, ...
                                sacd.Fc_ss, sacd.amp_cut, sacd.dur_cut, sacd.direction, sacd.pks, sacd.sacd_length, ...
                                sacd.min_pkdist, sacd.min_pkwidth, sacd.min_pkprom, ...
                                sacd.min_pkthresh, sacd.boundThresh, false);
    
%     figure (1)
%     pause
%     clf

   	% Body
    body_filt = filtfilt(b_low, a_low, body_scd.position);
    body_filt = filtfilt(b_high, a_high, body_filt);
    body_filt = body_filt - mean(body_filt);
    
    % Error
    error = func + body_filt*(-1 + gamma);
    error_filt = filtfilt(b_low, a_low, error);
    
    % Display
    display_filt = filtfilt(b_low, a_low, display_intrp);
    display_filt = filtfilt(b_high, a_high, display_filt);
    display_filt = display_filt - mean(display_filt);
    
    % Reconstruct the input sine function
    pat_func = display_filt - body_filt*(gamma);
    
    % Remove window around saccades from data
    body_rmv = body_filt;
    error_rmv = error_filt;
    display_rmv = display_filt;
    if ~isempty(body_scd.saccades_all)
        body_rmv(body_scd.saccades_all.Index) = nan;
        error_rmv(body_scd.saccades_all.Index) = nan;
        display_rmv(body_scd.saccades_all.Index) = nan;
    end
    
    % Store time domain data
    DATA.time{n} = tintrp;
    DATA.function{n} = func;
    DATA.body{n} = body_rmv;
    DATA.body_intrp{n} = body_filt;
    DATA.body_raw{n} = body_intrp - mean(body_intrp);
    DATA.body_saccade{n} = body_scd;
    DATA.error{n} = error_rmv;
    DATA.display{n} = display_rmv;
    DATA.error_display{n} = display_rmv - body_filt;
    DATA.function_display{n} = pat_func;
    
    % Sliding LSSA
    body_lssa = lssq_sliding(DATA.time{n}, DATA.body{n}, freq, win_size, overlap, R2_cut, false, false);
    input_lssa = lssq_sliding(DATA.time{n}, DATA.function{n}, freq, win_size, overlap, R2_cut, false, false);
    %error_lssa = lssq_sliding(DATA.time{n}, DATA.error{n}, freq, win_size, overlap, R2_cut, false, true);
    
    % Get the gain, phase, & complex response of the closed-loop response
    H.time = body_lssa.time;
    H.R2 = body_lssa.R2;
    H.gain = body_lssa.magnitude ./ input_lssa.magnitude;
    H.phase = body_lssa.phase - input_lssa.phase;
    H.complex = H.gain .* exp(1i*H.phase);
    H.compensation_error = abs( (1 + 0*1i) + H.complex*(-1 + gamma));
    
    % Calculate the open-loop response
    G.complex = H.complex ./ (1 + H.complex*(-1 + gamma));
    G.gain = abs(G.complex);
    G.phase = angle(G.complex);
    
    % Calculate the predicted closed-loop response
    if D.trial(n)==2 % learning phase where gamma does not equal 0
        gamma_pred = 0; % make prediciton for gamma=0
    else % baseline or relearn phase where gamma=0
        gamma_pred = gamma_folder; % make prediciton for gamma
    end
    H_prediction.complex = G.complex ./ (1 - G.complex*(-1 + gamma_pred));
    H_prediction.gain = abs(H_prediction.complex);
    H_prediction.phase = angle(H_prediction.complex);
    H_prediction.compensation_error = abs( (1 + 0*1i) + H_prediction.complex*(-1 + gamma_pred));
    
    % Fit linear trends to the transforms
    fit_line_fields = ["gain", "phase", "compensation_error" ];
    for k = 1:length(fit_line_fields)
        fd = fit_line_fields(k);
        H.fit_line.(fd) = fit_line(H.time, H.(fd), false);
        if ~strcmp(fd, fit_line_fields(3)) % no compensation error for G
            G.fit_line.(fd) = fit_line(H.time, G.(fd), false);
        end
    end
    
    % Store system ID data
    DATA.H{n} = H;
    DATA.G{n} = G;
    DATA.H_prediction{n} = H_prediction;
     
    % Plot
    %[out] = fit_sine(DATA.time{n}, DATA.body{n}, 0.5, true);
    figure (2) ; clf
    subplot(1,1,1) ; cla ; hold on
        plot(DATA.time{n}, DATA.function{n}, 'k')
        %plot(DATA.time{n}, DATA.body_raw{n}, 'Color', [0.5 0.5 0.5 0.7])
        plot(DATA.time{n}, DATA.function_display{n}, '-', 'Color', 'b')
        plot(DATA.time{n}, DATA.body{n}, 'r')
        plot(DATA.time{n}, DATA.error{n}, 'g')
        plot(DATA.time{n},DATA.error_display{n}, 'c')
	pause
end

%% Sort data based on condition
FLY = [];
names = string(DATA.Properties.VariableNames(4:end));
clss = ["base", "learn", "relearn"];
for n = 1:N.trial
    FLY.(clss(n)).all = DATA(DATA.trial == n, :);
    for d = 1:length(names)
        FLY.(clss(n)).(names(d)) = cat(2, FLY.(clss(n)).all.(names(d)){:});
    end
end

%% Camera times
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2 1.5])
ax = subplot(1,1,1); cla ; hold on
bin_size = 0.05;
bin_range = 50;
bins = 0:bin_size:(bin_range + bin_size/2);
camera_times_all = 1000*(cat(1,DATA.camera_times{:}));
h = histogram(camera_times_all, bins, 'Normalization', 'probability', ...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
xline(median(camera_times_all), 'k', 'LineWidth', 1);
title(median(camera_times_all))

ylabel('probability')
xlabel('processing time (ms)')

set(ax, 'Color', 'none', 'LineWidth', 1, 'XLim', [10 15], 'XTick', 10:1:15)
set(ax, 'YLim', [-0.005 0.08])

%% Body, error time
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 5 4])
clear ax h
lw = 1;
cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [0 0.6 0.1];

winI = Fs*0 + (round(0*Fs):round(10*Fs)) + 1;
win_time = FLY.base.time(winI);
win_time = win_time - win_time(1);
cc.error = [112 96 167]./255;
ax(1,1) = subplot(3,1,1); cla ; hold on ; ylabel('(°)')
    plot(win_time, FLY.base.function(winI,1), 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
    plot(win_time, nanmean(FLY.base.body(winI,1),2), 'Color', cc.base,  'LineWidth', lw)
    plot(win_time, nanmean(FLY.base.error(winI,:),2), 'Color', cc.error,  'LineWidth', lw)
    
ax(2,1) = subplot(3,1,2); cla ; hold on ; ylabel('(°)')
    plot(win_time, FLY.base.function(winI,1), 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
    plot(win_time, nanmean(FLY.learn.body(winI,1),2), 'Color', cc.learn,  'LineWidth', lw)
    plot(win_time, nanmean(FLY.learn.error(winI,1),2), 'Color', cc.error,  'LineWidth', lw)
    
ax(3,1) = subplot(3,1,3); cla ; hold on ; ylabel('(°)') ; xlabel('time (s)')
    plot(win_time, FLY.base.function(winI,1), 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
    plot(win_time, nanmean(FLY.relearn.body(winI,1),2), 'Color', cc.relearn,  'LineWidth', lw)
    plot(win_time, nanmean(FLY.relearn.error(winI,1),2), 'Color', cc.error,  'LineWidth', lw)
    
set(ax, 'Color', 'none', 'LineWidth', 1)
set(ax(:,1), 'XLim', [-0.3 round(range(win_time))])
set(ax, 'YLim', 55*[-1 1], 'YTick', -50:25:50)
set(ax(1:end-1), 'XColor', 'none')

%% Body, error full time
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1.5*[2 2 5 2])
clear ax h
cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [0 0.6 0.1];

clss = ["base", "learn", "relearn"];
n_clss = length(clss);
time_space = 30;

bin_size = 1;
bin_range = 200;
bins = -(bin_range + bin_size/2):bin_size:(bin_range + bin_size/2);

% func = 37.5*sin(2*pi*0.5*tt);
ax(1,1) = subplot(2,4,1:3); cla ; hold on ; ylabel('Body (°)')
    %plot(tt, func, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
    yline(37.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    yline(-37.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    for n = 1:n_clss
        if n > 1
            window_shift = sum(TimePeriods(1:n-1));
        else
            window_shift = 0;
        end
        time_shift = FLY.(clss(n)).time(:,1) + time_space*(n-1) + window_shift;
        h.body(1,n) = plot(time_shift, nanmean(FLY.(clss(n)).body,2), 'Color', cc.(clss(n)));
    end
    
ax(1,2) = subplot(2,4,4); cla ; hold on
    for n = 1:n_clss
        h.hist(1,n) = histogram(FLY.(clss(n)).body, bins, 'FaceColor', cc.(clss(n)));
    end
    
ax(2,1) = subplot(2,4,5:7); cla ; hold on ; ylabel('Error (°)')
    for n = 1:n_clss
        if n > 1
            window_shift = sum(TimePeriods(1:n-1));
        else
            window_shift = 0;
        end
        time_shift = FLY.(clss(n)).time(:,1) + time_space*(n-1) + window_shift;
        h.error(1,n) = plot(time_shift, nanmean(FLY.(clss(n)).error,2), 'Color', cc.(clss(n)));
    end
    xlabel('time (s)')
    
ax(2,2) = subplot(2,4,8); cla ; hold on
    for n = 1:n_clss
        h.hist(2,n) = histogram(FLY.(clss(n)).error, bins, 'FaceColor', cc.(clss(n)));
    end
    
set(ax, 'Color', 'none', 'LineWidth', 1)
set(h.hist, 'Normalization', 'probability', 'EdgeColor', 'none', 'Orientation', 'horizontal')

linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
linkaxes(ax(:,1), 'xy')

set([h.body , h.error], 'LineWidth', 0.5)

set(ax(:,1), 'XLim', [-20 1300], 'XTick', 0:100:1300)
set(ax, 'YLim', 125*[-1 1], 'YTick', -100:25:100)
% ax(1,2).XLim(1) = -0.05*ax(1,2).XLim(2);
% ax(2,2).XLim(1) = -0.05*ax(2,2).XLim(2);

set(ax(:,2), 'YColor', 'none')
set(ax(1,1), 'XColor', 'none')

%% Gain, phase & compensation error in time
clss = ["base", "learn", "relearn"];
metric = ["gain", "phase", "compensation_error", "R2"];
n_clss = length(clss);
n_metric = length(metric);
data = [];
for k = 1:n_metric
    for f = 1:N.fly
        data(f).(metric(k)) = [];
        data(f).G = [];
       	for n = 1:n_clss
            temp = FLY.(clss(n)).H(f).(metric(k));
            nanI = ~isnan(temp);
            temp = temp(nanI);
            data(f).(metric(k)) = [data(f).(metric(k)) ; temp];
            data(f).G = [data(f).G ; n*ones(size(temp))];
        end
    end
end


fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1.5*[2 2 3 4])
movegui(fig, 'center')
clear ax h

cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [0 0.6 0.1];
cc.R2 = [0 0 0];

rng(1)
spread = 0.3;
p = 1;
for k = 1:n_metric
    for f = 1:N.fly
        ax(k,f) = subplot(n_metric,N.fly,p); cla ; hold on
        if f == 1
            ylabel(metric(k), 'Interpreter', 'none')
        end
        
        if k == 1
           title(['fly ' num2str(f)]) 
        end
        
        plotdata = data(f).(metric(k));
        if strcmp(metric(k), 'phase')
            plotdata = rad2deg(plotdata);
        end
        
        % Plot
        bx = boxplot(plotdata, data(f).G, ...
            'Width', 0.5, 'Symbol', '', 'OutlierSize', 0.5);
        h = get(bx(5,:),{'XData','YData'});
        for c = 1:size(h,1)
           patch(h{c,1},h{c,2}, cc.(clss(c)), 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
        set(findobj(ax(k,f),'tag','Median'), 'Color', 'w','LineWidth', 1);
        set(findobj(ax(k,f),'tag','Box'), 'Color', 'none');
        set(findobj(ax(k,f),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
        set(findobj(ax(k,f),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
        ax(k,f).Children = ax(k,f).Children([end 1:end-1]);
        
        %jitter = rand(size(data(f).G)) * spread - (spread/2);
        %plot(data(f).G + jitter, plotdata, '.', 'Color', [0.5 0.5 0.5 0.1], 'MarkerSize', 2)
        
        p =  p + 1;
    end
end

set(ax, 'Color', 'none', 'LineWidth', 1, 'Box', 'off', 'XColor', 'none')
set(ax(1,:), 'YLim', [0 1.3])
set(ax(2,:), 'YLim', [-30 5])
set(ax(3,:), 'YLim', [0 1])
set(ax(4,:), 'YLim', [0 1])
% set(ax(:,1), 'XLim', [-0.3 round(range(win_time))])
% set(ax, 'YLim', 55*[-1 1], 'YTick', -50:25:50)
% set(ax(1:end-1), 'XColor', 'none')

%% SAVE
disp('Saving...')
savedir = 'E:\DATA\Reafferent';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'DATA', 'FLY', 'D', 'I', 'U', 'N', '-v7.3')
disp('SAVING DONE')
end