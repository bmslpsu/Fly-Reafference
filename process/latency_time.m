function [latency] = latency_time(root, n_use_frame)
%% getbodyang: gets the body angle of a fly in a frame
%
%   INPUTS:
%       root      	: root folder
%       n_use_frame	: how many frames to use. If empty, use all frames
%
%   OUTPUTS:
%       latency   	: vector of latency times in seconds
%

if nargin < 2
    n_use_frame = [];
end

[FILE,PATH] = uigetfile({'*.mp4'}, 'Select files', root, 'MultiSelect', 'on');
FILE = string(FILE);
n_file = length(FILE);

if isempty(n_use_frame)
   n_use_frame = V.NumFrames;
end

savedir = fullfile(PATH, 'latency');
mkdir(savedir)
for n = 1:n_file
    fpath = fullfile(PATH, FILE(n));
    V = VideoReader(fpath);
    
    latency = nan(n_use_frame,1);
    ang = nan;
    for f = 1:n_use_frame
        if ~mod(f,10)
            fprintf('%i \n', f)
        end
        frame = read(V, f);
        tic
        ang = getbodyang(frame(:,:,1), ang);
        latency(f) = toc;
    end
    fdata = split(FILE(n), '.');
    basename = char(fdata(1));
    savepath = fullfile(savedir, [basename '.mat']);
    save(savepath, '-v7.3', 'latency')
end

end