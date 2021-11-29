function [] = make_data(rootdir)
%% make_data:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%

gain = 1;
root = ['E:\Reafferent_Gain_experiment\gain=' num2str(gain)];
filename = ['data_gain=' num2str(gain)];

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root,'*.mat',false);

%% Get Data
close all
clc

% Saccade detection parameters
sacd.showplot = false;
sacd.Fc_detect = [10 nan];
sacd.Fc_ss = [nan nan];
sacd.amp_cut = 4;
sacd.dur_cut = inf;
sacd.thresh = [70 , 2, 1, 0];
sacd.true_thresh = 250;
sacd.sacd_length = nan;
sacd.pks = [];
sacd.min_pkdist = 0.1;
sacd.min_pkwidth = 0.03;
sacd.min_pkprom = 75;
sacd.min_pkthresh = 0;
sacd.boundThresh = [0.2 60];
sacd.direction = 0;


Fs = 93;
Fc = 10;
[b,a] = butter(5, Fc/(Fs/2), 'low');

DATA = [D , splitvars(table(num2cell(zeros(N.file,5))))];
DATA.Properties.VariableNames(5:end) = {'time', 'pattern','body','head','error'};
for n = 1:N.file
    %disp(kk)
    disp(basename{n})
    
    % Load DAQ, body, head, & wing data
	all_data = load(fullfile(root, [basename{n} '.mat']),'data','t_p', 'feedback', 'bAngles', 't_v'); % load data
    
    % Get synced frame times and pattern data
    trig = all_data.data(:,3);
    trig = 0.02*round(trig ./ 0.02);
    dxtrig = diff(trig);
    dxtrig = [dxtrig(1) ; dxtrig];
    [~,pks] = findpeaks(dxtrig, 'MinPeakHeight', 0.02, 'NPeaks', 1);
    pks_time = all_data.t_p(pks);
    time_sync = all_data.t_p - pks_time(1);

    vid_time = all_data.t_v;
    pat = all_data.data(:,1);
    pat = 3.75*round(96*(pat ./ 10));
    pat = rad2deg(unwrap(deg2rad(pat)));
    pat_pos = pat  - pat(1);
    body = all_data.bAngles - mean(all_data.bAngles) + mean(pat_pos);
    
    switch D.trial(n)
        case 1
            T = 60;
        case 2
            T = 600;
        case 3
            T = 600;
        otherwise
            error('file name error')
    end
   	tintrp = (0.5:1/Fs:T+0.5)';
    
    body_norm = interp1(vid_time, body, tintrp);
    pat_norm = interp1(time_sync, pat_pos, tintrp);
    body_filt = filtfilt(b, a, body_norm);
    pat_filt = filtfilt(b, a, pat_norm);
    %body_filt = body_filt - body_filt(1);
    pat_filt = pat_filt - pat_filt(1);
    
    %body_filt = rad2deg(wrapToPi(deg2rad(body_filt)));
    
    DATA.time{n} = tintrp - tintrp(1);
    DATA.pattern{n} = pat_filt;
    DATA.body{n} = body_filt;
    DATA.error{n} = body_filt - pat_filt;
    
%     body_scd = saccade_all(DATA.body{n}, DATA.time{n}, sacd.thresh, sacd.true_thresh, sacd.Fc_detect, ...
%                                 sacd.Fc_ss, sacd.amp_cut, sacd.dur_cut, sacd.direction, sacd.pks, sacd.sacd_length, ...
%                                 sacd.min_pkdist, sacd.min_pkwidth, sacd.min_pkprom, ...
%                                 sacd.min_pkthresh, sacd.boundThresh, sacd.showplot);
end

%% Sort data based on condition
FLY = [];
data = ["time", "body", "pattern"];
class = ["base", "learn", "relearn"];
for n = 1:N.trial
    FLY.(class(n)).all = DATA(DATA.trial == n, :);
    for d = 1:length(data)
        FLY.(class(n)).(data(d)) = cat(2, FLY.(class(n)).all.(data(d)){:});
    end
end

%% SAVE
disp('Saving...')
savedir = 'E:\Reafferent_Gain_experiment\processed';
save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
    'DATA', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end