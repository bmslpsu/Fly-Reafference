function [] = make_data(rootdir)
%% make_data:
%
%   INPUTS:
%       rootdir    	:   root directory
%
%   OUTPUTS:
%       -
%

% load(fullfile(path, file), 'bAngles', 't_v', 'data', 't_p', 'feedback')

root = 'E:\EXPERIMENTS\MAGNO\Experiment_body_reafferent\gain=0';

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(root,'*.mat',false);

%% Get Data
close all
clc

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

R = 180;
sz = 10;
hbins = (-R - sz/2):sz:(R + sz/2); 
cc = hsv(N.fly);

fig = figure (1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 6])
clear ax h
ax = gobjects(N.trial,2);
for n = 1:N.trial
    subI1 = (3*n - 2):(3*n - 1);
    ax(n,1) = subplot(N.trial,3,subI1); cla ; hold on
    title(class(n))
        h.trial{n} = plot(FLY.(class(n)).time, FLY.(class(n)).body, '-');
        set(h.trial{n}, {'Color'}, num2cell(cc,2))
    ax(n,1).XLim(1) = -0.02*FLY.(class(n)).time(end);
    
    ax(n,2) = subplot(N.trial,3,3*n); cla ; hold on
        for f = 1:N.fly
            h.hist(n,f) = histogram(FLY.(class(n)).body(:,f), hbins, 'Normalization', 'probability', ...
                'EdgeColor', 'none', 'FaceColor', cc(f,:));
        end
        ax(n,2).YLim(1) = -0.02;
end
cellfun(@(x) set(x, 'LineWidth', 0.5, 'MarkerFaceColor', 'none', 'MarkerSize', 2), h.trial)
set(ax, 'Color', 'none', 'LineWidth', 1)
% set(ax(1:end-1,1), 'XColor', 'none')
% set(ax(:,2), 'XLim', 200*[-1 1], 'XTick', -180:90:180)
% set(ax(:,1), 'YLim', 200*[-1 1], 'YTick', -180:90:180)

linkaxes(ax(:,2), 'xy')

XLabelHC = get(ax(N.trial,1), 'XLabel');
set(XLabelHC, 'String', 'time (s)')

YLabelHC = get(ax(2,1), 'YLabel');
set(YLabelHC, 'String', 'body (°)')

YLabelHC = get(ax(2,2), 'YLabel');
set(YLabelHC, 'String', 'probability')

leg = legend(h.trial{1}, string(1:N.fly));
leg.Title.String = 'Fly #';
leg.Position = [0.1666    0.4353    0.0997    0.1345];
% set(ax, 'XLim', [0 100])

%% SAVE
disp('Saving...')
% savedir = 'E:\DATA\Magno_Data\Multibody';
% save(fullfile(savedir, [filename '_' datestr(now,'mm-dd-yyyy') '.mat']), ...
%     'FUNC', 'DATA', 'GRAND', 'FLY', 'D', 'I', 'U', 'N', 'T', '-v7.3')
disp('SAVING DONE')
end