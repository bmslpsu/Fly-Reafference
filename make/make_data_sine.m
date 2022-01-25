function [] = make_data_sine()
%% make_data_sine:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%

gamma_folder = -1;
fpath = 'E:\EXPERIMENTS\MAGNO\Experiment_reafferent_sine';
root = fullfile(fpath, ['gamma=' num2str(gamma_folder)]);
filename = ['data_gamma=' num2str(gamma_folder)];

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root,'*.mat', false, 'fly', 'trial', 'gamma');

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
sacd.sacd_length = 3;
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
R2_cut = 0.7;

DATA = [D(:,1:3) , splitvars(table(num2cell(zeros(N.file,15))))];
DATA.Properties.VariableNames(4:end) = {'time', 'function', 'body', 'body_raw', 'body_intrp', ...
    'error', 'display', 'function_display', 'error_display', 'head', 'body_saccade', 'camera_times', ...
    'H', 'G', 'H_prediction'};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
        
    % Get system properties
    gamma = D.gamma(n); % gamma (coupling gain between fly body and display)
    gain = -1 + gamma; % total gain (natural + coupling)
    
    % Make time vector
    switch D.trial(n)
        case 1
            T = 30;
        case 2
            T = 905;
        case 3
            T = 305;
        otherwise
            error('file name error')
    end
   	tintrp = (0:1/Fs:T)';
    
    % Recreate the input sine function
    freq = 0.5;
    A = 37.5;
    Phi = 0;
    func = A*sin(2*pi*freq*tintrp + Phi);
    func = func - mean(func);
    
    % Load DAQ, body, head, & wing data
	all_data = load(fullfile(root, [basename{n} '.mat']),'data','t_p', 'bAngles', 't_v'); % load data

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
    H.gain = body_lssa.magnitude ./ input_lssa.magnitude;
    H.phase = body_lssa.phase - input_lssa.phase;
    H.complex = H.gain .* exp(1i*H.phase);
    H.compensation_error = abs( (1 + 0*1i) + gain*H.complex );
    
    % Calculate the open-loop response
    G.complex = H.complex ./ (1 - H.complex*(1 - gamma));
    G.gain = abs(G.complex);
    G.phase = angle(G.complex);
    
    % Calculate the predicted closed-loop response for gamma=0
    H_prediction.complex = G.complex ./ (1 + G.complex);
    H_prediction.gain = abs(H_prediction.complex);
    H_prediction.phase = angle(H_prediction.complex);
    H_prediction.compensation_error = abs( (1 + 0*1i) - H_prediction.complex);
    
    % Store system ID data
    DATA.H{n} = H;
    DATA.G{n} = G;
    DATA.H_prediction{n} = H_prediction;
     
%     % Plot
%     %[out] = fit_sine(DATA.time{n}, DATA.body{n}, 0.5, true);
%     figure (2) ; clf
%     subplot(1,1,1) ; cla ; hold on
%         plot(DATA.time{n}, DATA.function{n}, 'k')
%         %plot(DATA.time{n}, DATA.body_raw{n}, 'Color', [0.5 0.5 0.5 0.7])
%         %ssplot(DATA.time{n}, DATA.function_display{n}, '-', 'Color', 'b')
%         plot(DATA.time{n}, DATA.body{n}, 'r')
%         plot(DATA.time{n}, DATA.error{n}, 'g')
%         %plot(DATA.time{n},DATA.error_display{n}, 'c')
% 	%pause
end

%% Sort data based on condition
FLY = [];
names = string(DATA.Properties.VariableNames(4:end));
clss = ["base", "learn", "relearn"];
for n = 1:N.trial
    FLY.(clss(n)).all = DATA(DATA.trial == n, :);
    for d = 1:length(names)
        FLY.(clss(n)).(names(d)) = cat(2, FLY.(clss(n)).all.(names(d)){:});
        %FLY.(clss(n)).(names(d)) = FLY.(clss(n)).all.(names(d));
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
xline(median(camera_times_all), 'k', 'LineWidth', 1)
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

winI = 1 + round(0*Fs):round(10*Fs);
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
lw = 0.5;
cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [0 0.6 0.1];

bin_size = 1;
bin_range = 50;
bins = -(bin_range + bin_size/2):bin_size:(bin_range + bin_size/2);

tt = (0:0.001:40+900+20+300)';
func = 37.5*sin(2*pi*0.5*tt);
ax(1,1) = subplot(2,4,1:3); cla ; hold on ; ylabel('Body (°)')
    %plot(tt, func, 'Color', [0.5 0.5 0.5], 'LineWidth', lw)
    yline(37.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    yline(-37.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    plot(FLY.base.time, nanmean(FLY.base.body,2), 'Color', cc.base,  'LineWidth', lw)
    plot(FLY.learn.time + 40, nanmean(FLY.learn.body,2), 'Color', cc.learn, 'LineWidth', lw)
    plot(FLY.relearn.time + 40 + 900 + 20, nanmean(FLY.relearn.body,2), 'Color', cc.relearn, 'LineWidth', lw)
    
ax(1,2) = subplot(2,4,4); cla ; hold on
    h.hist(1,1) = histogram(FLY.base.body, bins, 'FaceColor', cc.base);
    h.hist(1,2) = histogram(FLY.learn.body, bins, 'FaceColor', cc.learn);
    h.hist(1,3) = histogram(FLY.relearn.body, bins, 'FaceColor', cc.relearn);
    
ax(2,1) = subplot(2,4,5:7); cla ; hold on ; ylabel('Error (°)')
    plot(FLY.base.time, nanmean(FLY.base.error,2), 'Color', cc.base,  'LineWidth', lw)
    plot(FLY.learn.time + 40, nanmean(FLY.learn.error,2), 'Color', cc.learn, 'LineWidth', lw)
    plot(FLY.relearn.time + 40 + 900 + 20, nanmean(FLY.relearn.error,2), 'Color', cc.relearn, 'LineWidth', lw)
    xlabel('time (s)')
    
ax(2,2) = subplot(2,4,8); cla ; hold on
    h.hist(2,1) = histogram(FLY.base.error, bins, 'FaceColor', cc.base);
    h.hist(2,2) = histogram(FLY.learn.error, bins, 'FaceColor', cc.learn);
    h.hist(2,3) = histogram(FLY.relearn.error, bins, 'FaceColor', cc.relearn);
    
    
set(ax, 'Color', 'none', 'LineWidth', 1)
set(h.hist, 'Normalization', 'probability', 'EdgeColor', 'none', 'Orientation', 'horizontal')
linkaxes(ax(:,1), 'xy')
linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')

set(ax(:,1), 'XLim', [-20 40*2+900 + 300])
set(ax, 'YLim', 55*[-1 1], 'YTick', -50:10:50)
% ax(1,2).XLim(1) = -0.05*ax(1,2).XLim(2);
% ax(2,2).XLim(1) = -0.05*ax(2,2).XLim(2);

set(ax(:,2), 'YColor', 'none')
set(ax(1,1), 'XColor', 'none')

%% Gain, phase & compensation error in time
clss = ["base", "learn", "relearn"];
metric = ["gain", "phase", "compensation_error"];
n_metric = length(metric);
data = [];
for k = 1:n_metric
    data(f).G = [];
    data(f).(metric(k)) = [];
    for n = 1:length(clss)
       	for f = 1:N.fly
            temp = FLY.(clss(n)).H(f).(metric(k));
            nanI = ~isnan(temp);
            temp = temp(nanI);
            data(f).(metric(k)) = [data(f).(metric(k)) ; temp];
            data(f).G = [data(f).G ; n*ones(size(temp))];
        end
    end
end


fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 5 4])
clear ax h
lw = 1;
cc.base = [0.9 0 0];
cc.learn = [0 0.7 1];
cc.relearn = [0 0.6 0.1];

p = 1;
for k = 1:n_metric
    for f = 1:N.fly
        ax(k,f) = subplot(n_metric,N.fly,p); cla ; hold on
        if f == 1
            ylabel(metric(k))
        end
        
        if k == 1
           title(['fly ' num2str(f)]) 
        end
        
        boxplot(data(f).(metric(k)), data(f).G)
        
        
        p =  p + 1;
    end
end

set(ax, 'Color', 'none', 'LineWidth', 1)
% set(ax(:,1), 'XLim', [-0.3 round(range(win_time))])
% set(ax, 'YLim', 55*[-1 1], 'YTick', -50:25:50)
% set(ax(1:end-1), 'XColor', 'none')


%% SAVE
% disp('Saving...')
% savedir = 'E:\Reafferent_Gain_experiment\processed';
% save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
%     'DATA', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
% disp('SAVING DONE')
end