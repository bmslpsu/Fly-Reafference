function [] = load_data()
%% load_dataload_data: loads the videuo, tracking & DAQ data for one trial
%
%   INPUT:
%       root     : root directory
%
%   OUTPUT:
%       FLY         : structure containing fly video & kinematics
%

clear ; close all ; clc

export = true;

root = 'Q:\OneDrive - PSU\OneDrive - The Pennsylvania State University\Research\Manuscripts\Reafferent\data\test';
rootpat = 'Q:\OneDrive - PSU\OneDrive - The Pennsylvania State University\Git\Arena\Patterns';

% Select pattern file
[FILE.pat, PATH.pat] = uigetfile({'*.mat'}, ...
    'Select pattern file', rootpat, 'MultiSelect','off');
pattern_data = load(fullfile(PATH.pat,FILE.pat),'pattern'); % load pattern

% Select .mat file
[FILE.mat, PATH] = uigetfile({'*.mat'}, ...
    'Select mat file', root, 'MultiSelect','off');

[~,basename,~] = fileparts(FILE.mat);
FILE.vid = [basename '.mp4'];

% Load .mat data
data = load(fullfile(PATH, FILE.mat));

% Video reader
Vin =  VideoReader(fullfile(PATH, FILE.vid));

%%
% Calculate frame rate and set time of video
Fs = 1 ./ mean(diff(data.t_v)); % frame rate [hz]
T = 10; % time of output video
n_frame = ceil(10*Fs); % # of frames of output video

% Select start frame
figure ; clf
plot(data.t_v, data.bAngles)
title('Select start point')
grid on
[x,~] = ginput(1); % select start point
close
start_time = 2*round(x./2); % start time to nearest 2 seconds
start_frame = round(start_time*Fs); % start frame
all_frames = start_frame:(start_frame + n_frame - 1); % all frames
all_time = data.t_v(all_frames); % all times
norm_time = all_time - all_time(1); % normalized time

FILE.movie = [basename '_movie_frame_' num2str(start_frame) '.mp4'];

% Load the video data
vid = uint8(zeros(Vin.Height, Vin.Width, n_frame));
for n = 1:n_frame
    frame = read(Vin, all_frames(n));
    vid(:,:,n) = frame(:,:,1);
end

%% Retrack the body
[B, imgstats] = bodytracker_vid(vid, 0, true, false);
centroid = cat(1,imgstats.Centroid);
centroid = medfilt1(centroid, 5, [], 1);
centroid_med = median(centroid,1);
body_radius = 0.2*median(cat(1,imgstats.MajorAxisLength));
offset = -2;
heading = centroid + body_radius*[cosd(B + offset) -sind(B + offset)];
close

% Low-pass filter body angles
[b, a] = butter(3, 10 / (Fs/2), 'low');
B = filtfilt(b, a, B);

% Recreate the input sine wave
R = data.ampl*sin(2*pi*data.freq*all_time);

% Compute the visual motion
gamma = data.gain;
V = R + gamma*B;

% Sync video and display data
D = data.data(:,1);
D = 3.75*(round(96*(D ./ 10)) - 1);
D = rad2deg(unwrap(deg2rad(D)));
dx = diff(D);
syncI = find(abs(dx) > 5, 1, 'first') + 1;
sync_time = data.t_p(syncI);
daq_time_sync = data.t_p - sync_time;
D_int = interp1(daq_time_sync, D, data.t_v, 'nearest');

% Get actual visual motion
V_exp = D_int(all_frames);

% Check
figure ; clf ; hold on ; title('Click if blue & green curves match')
plot(R - mean(R), 'k', 'LineWidth', 1)
plot(B - mean(B), 'r', 'LineWidth', 1)
plot(V_exp - mean(V_exp), 'g', 'LineWidth', 1)
plot(V - mean(V), 'b', 'LineWidth', 1)
axis tight
pause
close

% Get pattern image
pat_image = 255*pattern_data.pattern.Pats(1,:,1,5);
dim = size(vid);
pat_radius = floor(max(dim(1:2))/1.9);
pat_thickness = 6;

%% Make output video
if export
    movie_path = fullfile(PATH, 'movie');
    mkdir(movie_path)
    Vout = VideoWriter(fullfile(movie_path, FILE.movie),'MPEG-4');
    Vout.FrameRate = Fs;
    Vout.Quality = 100;
    open(Vout)
end

cc.B = [0.9 0 0];
cc.R = [0.8 0.8 0.8];
cc.V = [0 0.2 0.8];
lw = 3;
flip_vid = flipud(vid);

fig = figure (1); clf % main figure window for display & export
set(fig, 'Color', 'k', 'Renderer', 'OpenGL');
fig.Position(3:4) = [1320, 500];
movegui(fig, 'center')
clear ax ax_pat h

ax(1,1) = subplot(1,16,1:6); cla ; hold on; axis image
    title(['\gamma = ' num2str(gamma)], 'FontSize', 20, 'Color', 'w')
    freezeColors
    colormap(ax(1,1), gray)
    ax_off = 30;
    ylim([-ax_off , dim(1) + ax_off])
    xlim([-ax_off , dim(2) + ax_off])
    H.frame = imshow(flip_vid(:,:,1));
    H.heading = [];
    
    ax_pat(1) = axes; axis image
    set(ax_pat(1), 'Position', ax(1).Position, ...
        'XLim', ax(1).XLim, 'YLim', ax(1).XLim, 'Color', 'none')
    colormap(ax_pat(1), [0 0 0; cc.V])
    
    ax_pat(2) = axes; axis image
    set(ax_pat(2), 'Position', ax(1).Position, ...
        'XLim', ax(1).XLim, 'YLim', ax(1).XLim, 'Color', 'none')
    colormap(ax_pat(2), [0 0 0; cc.R])
        
ax(1,2) = subplot(1,16,8:16); cla ; hold on
    scale = 0.7;
    shift_up = 0.8*(1 - scale)*ax(1,2).Position(4);
    ax(1,2).Position(4) = scale*ax(1,2).Position(4);
    ax(1,2).Position(2) = shift_up + ax(1,2).Position(2);
    
    H.V = animatedline('Color', cc.V, 'LineWidth', lw);
    H.B = animatedline('Color', cc.B, 'LineWidth', lw);
    H.R = animatedline('Color', cc.R, 'LineWidth', lw);
    
    leg = legend([H.B H.R H.V], 'Body', 'Stimulus', 'Visual motion = Stimulus + \gamma Body');
    leg.Position = [0.4200    0.0427    0.2780    0.1600];
    set(leg', 'TextColor', 'w', 'Box', 'off')
    
    xlabel('Time (s)')
    ylabel('Angle (°)')
    
set(ax(2:end), 'FontSize', 16, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', ...
    'FontWeight', 'bold','LineWidth', 2, 'XLim', [-0.3 round(norm_time(end))])
set(ax(2:end), 'XTick', 0:10)
set(ax(2:end), 'YLim', 10*ceil(max(abs(B-mean(B)))./10)*[-1 1])
set(ax(2:end), 'YTick', -200:20:200)
set(ax_pat, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')

ax(2).YLabel.FontSize = 20;
ax(2).XLabel.FontSize = 20;

fig_sz = fig.Position(3:4);
for n = 1:n_frame
    set(fig, 'CurrentAxes', ax(1))
        frame = flip_vid(:,:,n);
        set(H.frame, 'CData', frame)
        delete(H.heading)
        H.heading(1) = plot([centroid(n,2) heading(n,2)], [centroid(n,1) heading(n,1)], '-', ...
            'Color', cc.B, 'LineWidth', 4, 'MarkerSize', 7);
     	H.heading(2) = plot(centroid(n,2), centroid(n,1), '.', ...
            'Color', cc.B, 'MarkerFaceColor', 'none', 'LineWidth', 4, 'MarkerSize', 25);
        
        
        set(fig, 'CurrentAxes', ax_pat(1))
        cla(ax_pat(1))
        theta = -deg2rad(V(n)) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
        x = 0.9*pat_radius * cos(theta) + centroid_med(1);
        y = 0.9*pat_radius * sin(theta) + centroid_med(2);
        z = zeros(1,length(x));
        hs = surface(ax_pat(1),[x;x],[y;y],[z;z],[pat_image;pat_image], ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', pat_thickness);
        
        cla(ax_pat(2))
        theta = -deg2rad(R(n)) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
        x = pat_radius * cos(theta) + centroid_med(1);
        y = pat_radius * sin(theta) + centroid_med(2);
        z = zeros(1,length(x));
        hs = surface(ax_pat(2),[x;x],[y;y],[z;z],[pat_image;pat_image], ...
            'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', pat_thickness);
        
    set(fig, 'CurrentAxes', ax(2)); hold on
        addpoints(H.B, norm_time(n), B(n) - mean(B))
        addpoints(H.R, norm_time(n), R(n))
        addpoints(H.V, norm_time(n), V(n) - mean(V))
        
    drawnow
    
    if export
        fig_frame = getframe(fig);
        fig_frame.cdata = fig_frame.cdata(...
            round(0.02*fig_sz(2)):end-round(0.02*fig_sz(2)), ...
            round(0.11*fig_sz(1)):end-round(0.06*fig_sz(1)), :);
        writeVideo(Vout, fig_frame);
    end
end
if export
    close(Vout)
end

end