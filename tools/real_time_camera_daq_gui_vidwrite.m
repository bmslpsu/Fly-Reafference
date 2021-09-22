function [vidData, t_v, bAngles, data, t_p, feedback] = ...
    real_time_camera_daq_gui_vidwrite(vid, clmode, gain, n_frame, showplot, savpath, vidwrite)
%% real_time_camera_daq_gui
% clc ; close all
daqreset
if nargin < 7
    vidwrite = true;
    if nargin < 6
        savpath = 'C:\Users\boc5244\Documents\MATLAB\test';
        if nargin < 5
            showplot = true;
        end
    end
end
pathparts = strsplit(savpath,filesep);
savedir = fullfile(pathparts{1:end-1});
% fname = pathparts{end};
vidpath = [savpath '.mp4'];
Vwriter = VideoWriter(vidpath,'MPEG-4');
Vwriter.FrameRate = 60;
open(Vwriter)

Panel_com('stop')
% patID = 2;
yPos = 5;
% Panel_com('set_pattern_id', patID);	% set pattern
% Panel_com('set_position', [40, yPos])

%% Set camera properties
% [vid,src] = Basler_acA640_750um();
stop(vid)
triggerconfig(vid, 'manual');
% ROI.x = 304;
% ROI.y = 300;
% ROI.xoff = (round(672 - ROI.x)/2);
% ROI.yoff = (round(512 - ROI.y)/2);
% vid.ROIPosition = [ROI.xoff ROI.yoff ROI.x ROI.y];
% roi = vid.ROIPosition;

%% Configure DAQ's
[s,~,~] = daq_realtime(10000, [1 2 3]);
logpath = fullfile(savedir, 'log.bin');
delete(logpath)
logid = fopen(logpath, 'w');
lh = addlistener(s.Read, 'DataAvailable', @(src, event, fid) storeData(src, event, logid));

%% Create variables to store data
% vidData = uint8(zeros(roi(4), roi(3), n_frame));
vidData = [];
t_v = zeros(n_frame,1);
bAngles = zeros(n_frame,1);
arena_gain = nan(n_frame,1);
arena_pos = nan(n_frame,1);

%% Loop to set initial heading angle
disp('Set heading')
[prev_ang] = set_heading(vid, false);

Panel_com('stop')    
pause(3)

%% Track in real time
n = 1;
% prev_ang = nan;
disp('start')
switch clmode % position or velocity feedback
    case 'pos'
        Panel_com('set_mode', [0,0]) % 0=open,1=closed,2=fgen,3=vmode,4=pmode
    case 'vel'
        Panel_com('set_mode', [0,0]) % 0=open,1=closed,2=fgen,3=vmode,4=pmode
        Panel_com('start')
    otherwise
       error('mode must be "pos" or "vel"') 
end

startBackground(s.Read);
start(vid) % start the video aquisition
while n <= n_frame
    if ~mod(n,100) || (n == 1) || (n == n_frame)
        disp(n)
    end
    [frame, ~] = getsnapshot(vid); % get snapshot
    [norm_ang, imgprop] = getbodyang(frame, prev_ang); % get body angle in frame
    
    arena_ang = 180 - norm_ang;
    bAngles(n) = arena_ang;
    if n == 1
        tic
        t_v(n) = 0;
    else
        t_v(n) = toc; 
    end
        
    switch clmode % position or velocity feedback
        case 'pos'
            panel_ang = gain * arena_ang;
            wrap_ang = rad2deg(wrapTo2Pi(deg2rad(panel_ang)));
            panel_pos = round(95*(wrap_ang / 360)) + 1;
            Panel_com('set_position', [panel_pos, yPos])
            arena_pos(n) = panel_pos;
        case 'vel' 
            if n > 2
                vel = mean( diff(bAngles(n-2:n)) ) ./ mean( diff(t_v(n-2:n)) );
            else
             	vel = 0;
            end
            panel_gain = gain * (vel ./ 3.75); % normalized panel gain
            if panel_gain > 100
                warning('panel gain saturated')
                panel_gain = 100;
            elseif panel_gain < -100
                warning('panel gain saturated')
                panel_gain = -100;
            end
            panel_gain = round(panel_gain);
            %disp(panel_gain)
            Panel_com('stop')
            Panel_com('send_gain_bias', [panel_gain 0 0 0])
            Panel_com('start')
            arena_gain(n) = panel_gain;
    end
    
    if showplot
        if n > 1
            delete([h.cent h.heading])
            set(H, 'cdata', frame)
        else
            H = imshow(frame); hold on
        end
        cent = imgprop.Centroid; % centroid
        L = imgprop.MajorAxisLength / 3; % long axis
        head = cent + L*[sind(norm_ang), -cosd(norm_ang)];
        heading = [cent ; head];
        h.heading = plot(heading(:,1), heading(:,2), '-r', 'LineWidth', 1.5);
        h.cent = plot(cent(1), cent(2), '.r', 'MarkerSize', 15);
        drawnow
    end
    
    if vidwrite
        writeVideo(Vwriter, frame);
    end
    %vidData(:,:,n) = frame;
    prev_ang = norm_ang;
    n = n + 1;
end
stop(vid) % stop video aquisition
flushdata(vid) % flush all the image data stored in the memory buffer
stop(s.Read)
fclose(logid);
close(Vwriter)
Panel_com('stop')

logid_read = fopen(logpath, 'r');
[data,~] = fread(logid_read, [4,inf], 'double');
data = data';
t_p = data(:,1);
data = data(:,2:end);
feedback.arena_gain = arena_gain;
feedback.arena_pos = arena_pos;

fclose(logid_read);
delete(logpath)

fs = 1 / mean(diff(t_v));
disp(['Frame rate = ' num2str(fs)])

close(gcf)

end

%% DAQ listener
function storeData(src, event, fid)
    %[data, timestamps, ~] = read(src, src.ScansAvailableFcnCount, "OutputFormat", "Matrix");
    data = [event.TimeStamps, event.Data]';
    fwrite(fid,data,'double');
end