function [norm_ang, imgprop] = getbodyang(frame,  prev_ang)
%% getbodyang: gets the body angle of a fly in a frame
%
%   INPUTS:
%       frame   	: 2D frame (uint8)
%       prev_ang   	: previous angle [°]
%
%   OUTPUTS:
%       bAngle   	: body angle [°]
%       imgprop  	: image properties
%

offset  = 90;  % shift the reference frame so 0° is the top vertical axis in video [°]
dthresh = 160; % threshold for detecting >180° flips in ellipse orientation, or angles > 360°
shift   = 0;   % shift to keep angle wrapped [°] (dynamic)

% SE_erode = strel('disk',5,4); % erosion mask
% SE_dilate = strel('disk',12,4); % dilation mask

bw = imbinarize(frame); % binarize
bw = medfilt2(bw,[9 9]); % median filter
%bw = imerode(bw, SE_erode); % erode

% Get image properties
imgprop = regionprops(bw,'Centroid','Area','BoundingBox','Orientation', ...
    'MajorAxisLength','MinorAxisLength'); % image reigon properties

% Pick out the body (object with largest area)
[~,sort_area] = sort([imgprop.Area],'descend');
bodyI = sort_area(1); 
imgprop = imgprop(bodyI);

% Get the angle
norm_ang = -(imgprop.Orientation - offset + shift); % normalized, unwrapped angle [°]

% Check for changes in angle > 180°. Correct for the ellipse fit and unwrap angles.
if ~isnan(prev_ang)
    dbody = norm_ang - prev_ang; % change in body angle between frames [°]
    magd = abs(dbody); % magnitude of change [°]
    signd = sign(dbody); % direction of change
    if magd > dthresh % 180° or more flip
        flipn = round(magd/180); % how many 180° shifts we need
        shift = -signd*flipn*180; % shift amount [°]
        norm_ang = norm_ang + shift; % normalized, unwrapped angle [°]
    end
else
    % Make start angle positive
    if norm_ang < 0
        norm_ang = norm_ang + 360;
    end
end

end