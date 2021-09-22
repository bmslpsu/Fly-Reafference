function [] = Experiment_body_reafferent_control(Fn)
%% Experiment_body_reafferent_control: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
beep on
% Fn = 0;

%% Set directories & experimental parameters
gain = 1;
root = ['C:\BC\Experiment_body_reafferent\gain = ' num2str(gain)];
mkdir(root)

%% EXPERIMENTAL PARAMETERS
patID = 1;
yPos = 5;
mode = 'pos';
fps = 100;
gain_all = [0 gain 0];
n_frame = [1*60, 10*60, 10*60] .* fps;
n_rep = length(gain_all); % # of repetitions

%% Camera Setup
[vid,src] = Basler_acA640_750um();

%% EXPERIMENT LOOP
disp('Start Experiment:')
% Panel_com('quiet_mode_off')
for n = 1:n_rep
    fprintf('Trial: %i \n', n)
    fname = ['fly_' num2str(Fn) '_trial_' num2str(n) '_gain_' num2str(gain_all(n)) '_' mode '.mat'];
    %preview(vid) % open video preview window
    
    Panel_com('stop')    
	pause(1) % pause between buffer & experiment
    
    % EXPERIMENT SETUP
    disp('Play Stimulus:')
    Panel_com('set_pattern_id', patID);	% set pattern
    Panel_com('set_position', [49, yPos]); % set starting position (xpos,ypos)
    %Panel_com('set_mode', [0,0]); % 0=open,1=closed,2=fgen,3=vmode,4=pmode
	
    % START EXPERIMENT & DATA COLLECTION
    vidpath = fullfile(root, fname(1:end-4));
    [vidData, t_v, bAngles, data, t_p, feedback] = ...
        real_time_camera_daq_gui_vidwrite(vid, mode, gain_all(n), n_frame(n), false, vidpath, true);
    
    for b = 1:5
        beep
        pause(1)
    end
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v', 'bAngles', 'feedback');
end

delete(vid)
disp('Done');
daqreset
imaqreset
end