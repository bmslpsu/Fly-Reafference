function [FLY, DATA, FILE, PATH] = load_trial_data(rootdir, vidFPS, showplot)
%% load_trial_data: loads the videuo, tracking & DAQ data for one trial
%
%   INPUT:
%       rootdir     : root directory
%
%   OUTPUT:
%       FLY         : structure containing fly video & kinematics
%       DATA    	: structure containing raw data
%

if nargin < 3
    showplot = true;
end

pat_ypos = 5;

if ~isfolder(rootdir)
    dirflag = false;
    [rootdir,mainfile,mainext] = fileparts(rootdir);
    FILE.raw = [mainfile , mainext];
else
    dirflag = true;
end


% Select pattern file
rootpat = 'C:\Users\boc5244\Documents\GitHub\Arena\Patterns';
[FILE.pat, PATH.pat] = uigetfile({'*.mat'}, ...
    'Select pattern file', rootpat, 'MultiSelect','off');

% Create data paths
if dirflag 
    [FILE.raw, PATH.raw] = uigetfile({'*.mat'}, ... % select fly file manually
        'Select file (only from given root directory)', ...
        rootdir, 'MultiSelect','off');
end

% Make movie path
PATH.movie = fullfile(PATH.raw, 'movie');
mkdir(PATH.movie);

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
FILE.mat = [FILE.basename '.mat'];
FILE.vid = [FILE.basename '.mp4'];

% Get gamma value
finfo = split(FILE.basename, '_');
gamma = str2double(finfo(6));

% Load data
disp('Loading Data ...')
DATA.main = load(fullfile(PATH.raw,FILE.mat),'data','t_p','vidData','t_v', 'bAngles');
DATA.pattern = load(fullfile(PATH.pat,FILE.pat),'pattern');

% Set up video reader
FLY.vidIn = VideoReader(fullfile(PATH.raw, FILE.vid));

% Set up video writer
FLY.vidOut = VideoWriter(fullfile(PATH.movie, FILE.vid),'MPEG-4');
FLY.vidOut.Quality = 100;
FLY.vidOut.FrameRate = vidFPS;

%% Sync video & display data & filter
% Times
daq_time = DATA.main.t_p;
fly_time = DATA.main.t_v;
stop_time = 5*floor( max(fly_time) ./ 5);
[~,stopI] = min(abs(fly_time - stop_time));
main_time = fly_time(1:stopI);

% Display & body
body = DATA.main.bAngles;
display = DATA.main.data(:,1);
display = 3.75*(round(96*(display ./ 10)));
display = rad2deg(unwrap(deg2rad(display)));
display = hampel(display);

% Find start time & sync DAQ to body
dx_pat = diff(display);
syncI = find(abs(dx_pat) > 5, 1, 'first') + 1;
sync_time = daq_time(syncI);

daq_time = daq_time(syncI:end);
daq_time_sync = daq_time - sync_time;

display = display(syncI:end);

% Filter
Fs = 1 / mean(diff(daq_time));
Fc_low = 3;
[b_low,a_low] = butter(5, Fc_low/(Fs/2), 'low');
display_filt = filtfilt(b_low, a_low, display);

Fs = 1 / mean(diff(fly_time));
Fc_low = 20;
[b_low,a_low] = butter(5, Fc_low/(Fs/2), 'low');
body_filt = filtfilt(b_low, a_low, body);

% Interpolate
display_intrp = interp1(daq_time_sync, display_filt, main_time, 'nearest');
body_intrp = interp1(fly_time, body_filt, main_time, 'pchip');

% Reconstruct input function & display
func_real = display_intrp - gamma*body_intrp;

% Error
error = display_intrp - body_intrp;

%% Receate input function 7 display
func = 37.5*sin(2*pi*0.5*main_time);
display_theoretical = func + gamma*body_intrp;

%% Store data
FLY.time = main_time;
FLY.func = func;
FLY.body = body_intrp;
FLY.display = display_intrp;
FLY.display_theoretical = display_theoretical;
FLY.error = error;
FLY.Fs = Fs;
FLY.func_real = func_real;
FLY.ypos = pat_ypos;

%% Figure
if showplot
    clf ; hold on
    plot(FLY.time, FLY.func, 'k')
    plot(FLY.time, FLY.func_real, 'b')
    plot(FLY.time, FLY.body, 'r')
    plot(FLY.time, FLY.display, 'g')
    plot(FLY.time, FLY.display_theoretical, 'm')
end

%% Get kinematic locations
% switch type % load registered video
%     case 'magno'
%         FLY.body_centroid   = cat(1,DATA.body.imgstats(TRIG.range).Centroid); % centroid of image
%         FLY.body_length     = cat(1,DATA.body.imgstats(TRIG.range).MajorAxisLength); % long axis of image
%     case 'rigid'
%         FLY.body_centroid   = repmat([FLY.raw_yP/2, FLY.raw_xP/2], [FLY.nframe , 1]);
%         FLY.body_length     = (FLY.raw_xP/2)*ones(FLY.nframe, 1);
% end
% 
% FLY.head_length = scaleHead * mean(FLY.body_length) / 6;
% FLY.body_glob = DATA.head.head_mask.global;
% FLY.head_hinge = fliplr(fix(padsize)) + DATA.head.head_mask.move_points.rot + [shiftR 0];
% FLY.head_tip   = FLY.head_hinge + FLY.head_length*[sind(FLY.int_head + FLY.body_glob) , ...
%                     -cosd(FLY.int_head + FLY.body_glob)];
% 
% FLY.body_length = scaleBody*FLY.body_length;
% FLY.body_tip = FLY.body_centroid + repmat((FLY.body_length/2),1,2).* ...
%                                     [sind(FLY.int_body), -cosd(FLY.int_body)];
% FLY.reverse = FLY.body_centroid + repmat((FLY.body_length/2),1,2).* ...
%                                      -[sind(FLY.int_body), -cosd(FLY.int_body)];
% FLY.heading = [FLY.body_centroid , FLY.body_tip];
% FLY.abdomen = [FLY.body_centroid , FLY.reverse];

end