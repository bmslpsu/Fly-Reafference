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
fpath = 'Q:\OneDrive - PSU\OneDrive - The Pennsylvania State University\Research\Data\Experiment_reafferent_sine';
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

DATA = [D(:,1:3) , splitvars(table(num2cell(zeros(N.file,10))))];
DATA.Properties.VariableNames(4:end) = {'time', 'function', 'body', 'body_raw', ...
    'error', 'display', 'function_display', 'error_display', 'head', 'body_saccade'};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    gamma = D.gamma(n); % gamma (coupling gain between fly body and display)
    
    % Make time vector
    switch D.trial(n)
        case 1
            T = 30;
        case 2
            T = 600;
        case 3
            T = 300;
        otherwise
            error('file name error')
    end
   	tintrp = (0:1/Fs:T)';
    
    % Recreate the input sine function
    freq = 0.5;
    func = (37.5 + 0*3.75)*sin(2*pi*freq*tintrp);
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
    
    % Store data
    DATA.time{n} = tintrp;
    DATA.function{n} = func;
    DATA.body{n} = body_rmv;
    DATA.body_raw{n} = body_intrp - mean(body_intrp);
    DATA.body_saccade{n} = body_scd;
    DATA.error{n} = error_rmv;
    DATA.display{n} = display_rmv;
    DATA.error_display{n} = display_rmv - body_filt;
    DATA.function_display{n} = pat_func;
    
    % Sliding LSSA
    out = lssq_sliding(DATA.time{n}, DATA.body{n}, freq, win_size, overlap, true);
    
%     % Plot
%     %[out] = fit_sine(DATA.time{n}, DATA.body{n}, 0.5, true);
%     figure (2) ; clf
%     subplot(1,1,1) ; cla ; hold on
%         plot(DATA.time{n}, DATA.function{n}, 'k')
%         plot(DATA.time{n}, DATA.body_raw{n}, 'Color', [0.5 0.5 0.5 0.7])
%         plot(DATA.time{n}, DATA.function_display{n}, '-', 'Color', 'b')
%         plot(DATA.time{n}, DATA.body{n}, 'r')
%         plot(DATA.time{n}, DATA.error{n}, 'g')
%         plot(DATA.time{n},DATA.error_display{n}, 'c')
% 	%pause
end

%% Sort data based on condition
FLY = [];
names = string(DATA.Properties.VariableNames(4:end-1));
clss = ["base", "learn", "relearn"];
for n = 1:N.trial
    FLY.(clss(n)).all = DATA(DATA.trial == n, :);
    for d = 1:length(names)
        FLY.(clss(n)).(names(d)) = cat(2, FLY.(clss(n)).all.(names(d)){:});
        %FLY.(clss(n)).(names(d)) = FLY.(clss(n)).all.(names(d));
    end
end

%% SAVE
% disp('Saving...')
% savedir = 'E:\Reafferent_Gain_experiment\processed';
% save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
%     'DATA', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
% disp('SAVING DONE')
end