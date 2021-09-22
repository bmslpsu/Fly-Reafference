function [norm_ang] = set_heading(vid, auto)
%% set_heading: moniter fly heading in real time and click to flip initial heading estimate by 180 deg
%
%   INPUT:
%       vid : video object
%
%   OUTPUT:
%       prev_ang : last angle [deg]
%

if nargin < 2
    auto = false;
end

global flip loop
flip = false;
loop = true;

yPos = 5;

prev_ang = nan;
[frame, ~] = getsnapshot(vid); % get snapshot
H = imshow(frame); hold on
h.heading = [];
h.cent = [];
% fig = figure ; clf
uicontrol('Style','pushbutton','String','Flip','Callback',@flipcheck,'Position', [10 150 60 30])
uicontrol('Style','pushbutton','String','Skip','Callback',@skipcheck,'Position', [10 190 60 30])
start(vid) % start the video aquisition
while loop
    [frame, ~] = getsnapshot(vid); % get snapshot
    [norm_ang, imgprop] = getbodyang(frame, prev_ang); % get body angle in frame
    
    if flip
        norm_ang = norm_ang + 180;
        flip = false;
    end
    
    arena_ang = 180 - norm_ang;
    wrap_ang = rad2deg(wrapTo2Pi(deg2rad(arena_ang)));
    V = round(95*(wrap_ang / 360)) + 1;
    %Panel_com('set_position', [V, yPos])
    
    delete([h.cent h.heading])
    set(H, 'cdata', frame)
    cent = imgprop.Centroid; % centroid
    L = imgprop.MajorAxisLength / 3; % long axis
    head = cent + L*[sind(norm_ang), -cosd(norm_ang)];
    heading = [cent ; head];
    h.heading = plot(heading(:,1), heading(:,2), '-r', 'LineWidth', 1.5);
    h.cent = plot(cent(1), cent(2), '.r', 'MarkerSize', 15);
    drawnow
    
    prev_ang = norm_ang;
    if auto
        loop = false;
    end
end
stop(vid) % stop video aquisition
flushdata(vid) % flush all the image data stored in the memory buffer
close(gcf)

end

%% Callback function for flip and skip
function flipcheck(~,~)
    global flip
    flip = true;
end

function skipcheck(~,~)
    global loop
    loop = false;
    close(gcf)
end