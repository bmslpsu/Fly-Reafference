function [] = montage_realtime_sine()
%% montage_realtime: makes movie for fly in magnetic tether
%
% 	Includes fly video, registered video, body tracking, head tracking, 
%   wing tracking,& pattern position
%
%   INPUT:
%       rootdir     : directory containing all files
%       vidFs       : video display FPS
%       export      : boolean (1=export video to images)
%
%   OUTPUT:
%       MOV         : structure containing movie 
%

clear ; clc ; close all
root = 'E:\EXPERIMENTS\MAGNO\Experiment_reafferent_sine\gamma=-1';
export = false;
vidFPS = 88;
pat_ypos = 5;

[FLY, DATA, ~, ~] = load_trial_data(root, vidFPS, true);
n_frame = FLY.vidIn.NumFrames;

%% Set parameters
pat_win = -40*[1 1];
pat_radius = 105;
pat_thickness = 5;

%% Load first frame
init_frame = read(FLY.vidIn,1);
dim = size(init_frame);

%% Get centroid
SE = strel('disk', 10, 8);
centroid = nan(n_frame, 2);
tic
parfor n = 1:n_frame
    disp(n)
    frame = read(FLY.vidIn, n);
    frame = fliplr(medfilt2(frame(:,:,1), [10 10]));
    bw = imbinarize(frame);
    bw = imerode(bw, SE);
    bw = imdilate(bw, SE);
    [y,x] = find(bw);
    centroid(n,:) = [mean(x) mean(y)];
%     if n == 1
%         H = imshow(bw);
%     else
%         set(H, 'CData', bw)
%     end
%     drawnow
end
toc

%% Compute heading
[~,stopI] = min(abs(DATA.main.t_v - 30));
centroid = centroid(1:stopI,:);
heading = centroid + 50*[-sind(FLY.body) , cosd(FLY.body)];

%% Make Movie
FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL','Position', 0.8*[100, 100, 2000, 600]);
% set(FIG, 'Visible','off');
linewidth = 1.25; % displayed line width
fontsize = 14;
clear h ax ax_pat h_gaze
ax(1) = subplot(1,10,1:3) ; cla; hold on; axis square % for raw fly & pattern vid
        %title('Arena Frame','Color','w','FontSize', 1.2*fontsize)
        ax(1).XLim = [-pat_win(1) pat_win(1)+dim(1)];
        ax(1).YLim = [-pat_win(2) pat_win(2)+dim(1)];
        ax_pat(1) = axes; axis image
        set(ax_pat(1), 'Position', ax(1).Position, 'XLim', ax(1).XLim, 'YLim', ax(1).XLim)
        %ax(1).Position(3:4) = 0.9*ax(1).Position(3:4);
ax(2) = subplot(1,10,5:10)  ; cla ; hold on
        ylabel('(°)','Color','w','FontSize',fontsize)
        xlabel('Time (s)','Color','w','FontSize',fontsize)
     	h.func = animatedline('Color', [0.7 0.7 0.7],'LineWidth',linewidth); % for pattern angle
        h.display = animatedline('Color','g','LineWidth',linewidth, 'LineStyle', '--'); % for pattern angle
        h.body = animatedline('Color','r','LineWidth',linewidth); % for body angle
        leg = legend([h.body h.func h.display], 'body', 'input', 'display', 'Box', 'off', ...
            'TextColor', 'w', 'FontWeight', 'bold', 'FontSize', 14);
        ax(2).Position(4) = 0.7*ax(2).Position(4);
        ax(2).Position(2) = 2.2*ax(2).Position(2);
        %ax(2).XLim(1) = -5;
        
set(ax_pat, 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
set(ax(2:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5, 'XLim', [-2 round(FLY.time(end))])
% set(ax(end),'XTick', 0:2:round(FLY.int_time(end)))

gs = 2;
cmap = [zeros(gs,1), linspace(0,1,gs)', zeros(gs,1)];
colormap(cmap)
pat_image = 255*DATA.pattern.pattern.Pats(1,:,1,pat_ypos);

fig_sz = FIG.Position(3:4);
if export
    open(FLY.vidOut)
end

disp('Exporting Video...')
tic
for n = 1:n_frame
    if ~mod(n,100) || (n == 1)
        disp(n)
    end
    
    % Get frames
    frame = read(FLY.vidIn, n);
    frame = fliplr(frame);
    frame = 1.0*imadjust(medfilt2(frame(:,:,1), [3 3])); % raw frame
    
    % Display video
    set(FIG, 'CurrentAxes', ax(1)); hold on
        if n == 1
            H.frame = imshow(frame); axis image
        else
            set(H.frame, 'CData', frame);
            delete(H.heading)
            delete(H.centroid)
        end
        H.heading = plot([centroid(n,1) heading(n,1)], ...
            [centroid(n,2) heading(n,2)], 'r', 'LineWidth', 2);
        H.centroid = plot(centroid(n,1), centroid(n,2), 'r.', 'MarkerSize', 25);
        %plot(mean(FLY.body_tip(win,1)), mean(FLY.body_tip(win,2)), 'c.', 'MarkerSize', 15)

        if pat_ypos ~= 1
            pat_pos = FLY.display_theoretical(n);
            % Make pattern ring
            set(FIG, 'CurrentAxes', ax_pat(1))
            cla(ax_pat(1))
            theta = -deg2rad(pat_pos) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
            x = pat_radius * cos(theta) + dim(2)/2 - 2;
            y = pat_radius * sin(theta) + dim(1)/2;
            z = zeros(1,length(x));
            hs = surface(ax_pat(1),[x;x],[y;y],[z;z],[pat_image;pat_image], ...
                'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', pat_thickness);
        end
    
    % Body & pattern plot
 	set(FIG, 'CurrentAxes', ax(2))
        addpoints(h.func, FLY.time(n), FLY.func(n))
        addpoints(h.body, FLY.time(n), FLY.body(n) - mean(FLY.body(1)))
        addpoints(h.display, FLY.time(n), FLY.display_theoretical(n) - mean(FLY.display_theoretical))
        
    drawnow
    
    if export
        fig_frame = getframe(FIG);
        fig_frame.cdata = fig_frame.cdata(...
            round(0.1*fig_sz(2)):end-round(0.1*fig_sz(2)), ...
            round(0.12*fig_sz(1)):end-round(0.07*fig_sz(1)), :);
        writeVideo(FLY.vidOut,fig_frame);
    end
end
toc

if export
 	disp('Saving...')
    pause(1)
    close(FLY.vidOut)
end
disp('DONE')
end