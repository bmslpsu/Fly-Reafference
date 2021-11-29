function [W_norm, T_norm, W, T] = pull_windows(data, sort_vector, win_length, fs, norm, showplot)
%% pull_windows: pulls out part of signal based on logical values
%
% 	Find where the window signal is high & pull out data
%
%   INPUT:
%       data        : data signal to pull out windows from
%       sort_vector	: window vector, true values are windows, false values are between windows
%       win_length  : window length in # of indicies, leave empty or nan for automatic detection
%                     or NaN for automatic median length
%       fs          : sampling frequency to create time vectors, leave empty if not needed
%       norm        : (boolean) normalize windows to pre-window
%       showplot  	: BOOLEAN show plot
%
%   OUTPUT:
%       W           : cell array with extracted windows
%       T           : cell array with time vector corresponding to each window
%       T_norm   	: cell array with time vector corresponding to each window normalized to start time
%

if nargin < 6
    showplot = false;
    if nargin < 5
        norm = [false false];
        if nargin < 4
            fs = 1;
            if nargin < 3
                win_length = [];
            end
        end
    end
end

% Check inputs
dim = size(data);
n_point = length(sort_vector);
assert(dim(1) == n_point, 'signal and windows vector must have the same number of rows')

if isempty(fs) % no sampling frequency specified
   fs = 1;
end

if isempty(win_length) % no window length specified
    win_length = [0 0];
elseif length(win_length) == 1 % Just window lenght, no pre-window
    win_length = [0 win_length];
end

if length(norm) == 1
    norm = repmat(norm, [1 2]);
end

assert( (isa(win_length(2), 'double') && (win_length(2) >= 1)) || ...
                (win_length(2) == 0) || isnan(win_length(2)), ...
                 ['"win_length" must be a positive scalar double', ...
                 ' > 1 or empty for automatic window length detection', ...
                 ' or NaN for automatic median window length'])
             
% Time vector
time = (1:n_point)' ./ fs;

% Find window start points
sort_vector = logical(round(sort_vector));
p = [1, -1]; % extract positive and negative windows
W = cell(2,1);
W_norm = cell(2,1);
T = cell(2,1);
T_norm = cell(2,1);
startI = cell(2,1);
for k = 1:2 % extract positive and negative windows
    sort_diff = diff(p(k)*double(sort_vector));
    sort_diff = [sort_diff ; sort_diff(end)];
    [~,sI] = findpeaks(sort_diff); % start points
    n_win = length(sI); % # of windows

    % Find window end points & add pre-window if specified
    if win_length(2) == 0 % find window lengths automatically
        [~,eI] = findpeaks(-sort_diff); % end points
        if length(eI) < n_win % if last window does not end, just use the last point of array
           eI = [eI ; n_point];
        end
    elseif isnan(win_length(2)) % find window lengths automatically, but use median windows size for all
        [~,eI] = findpeaks(-sort_diff); % end points
        if length(eI) < n_win % if last window does not end, just use the last point of array
           eI = [eI ; n_point];
        end
        med_win_sz = median(eI - sI); % median window size
        eI = sI + med_win_sz;
    else % set windows lengths based on constant size
        eI = sI + win_length(2);
        if any(eI > n_point)
           warning('some windows end points are outside signal length') 
        end
        eI(eI > n_point) = n_point;
    end
    winI = eI - sI;
    
    % Pre-window size
    if win_length(1) > 0
       preI = sI - win_length(1);
    else
       preI = sI;
    end

    % Get windows
    w = cell(1,n_win);
    w_norm = cell(1,n_win);
    t = cell(1,n_win);
    t_norm = cell(1,n_win);
    for n = 1:n_win
        win_span = (sI(n):eI(n))';
        pre_span = (preI(n):sI(n)-1)';
        all_span = [pre_span ; win_span];
        
        % Only window
        w{n} = data(win_span,:,:);
        t{n} = win_span ./ fs;
        
        % Normalized window
        w_norm{n} = data(all_span,:,:);
        switch norm(k)
            case true
                w_base = mean(data(pre_span,:,:));
            case false
                w_base = 0;
            otherwise
                error('"norm" must be true or false')
        end
        w_norm{n} = w_norm{n} - w_base;        
        t_norm{n} = [(-length(pre_span):-1)' ; (win_span - win_span(1))] ./ fs;
    end

    if length(unique(winI)) == 1 % if all windows are the same length, concatenate cells
        w = cat(2, w{:});
        w_norm = cat(2, w_norm{:});
        t = cat(2, t{:});
        t_norm = cat(2, t_norm{:});
    end
    
    W{k,1} = w;
 	W_norm{k,1} = w_norm;
    T{k,1} = t;
    T_norm{k,1} = t_norm;
end

if showplot
    fig = figure (304); clf
    set(fig,'Color','w','Units','inches')
    clear ax

    ax(1) = subplot(3,1,1) ; hold on ; axis tight
    yyaxis right
    ax(1) = gca; set(ax(1), 'YColor', 'g')
        xx = [time ; time(end) ; time(1)];
        yy = double([sort_vector ; sort_vector(1); sort_vector(1)]);
        patch(xx, yy, 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.2)
        %plot(time, sort_vector, 'g')
    
    yyaxis left
    ax(2) = gca; set(ax(2), 'YColor', 'k')
        plot(time, data, 'k-')
        plot(T{1}, W{1}, 'r-')
     	plot(T{2}, W{2}, 'b-')
  
    ax(3) = subplot(3,1,2) ; hold on
        yline(0, 'k--')
        plot(T_norm{1}, W_norm{1}, 'r')
        plot(mean(T_norm{1},2), mean(W_norm{1},2), 'k', 'LineWidth', 1)

    ax(4) = subplot(3,1,3) ; hold on
        yline(0, 'k--')
        plot(T_norm{2}, W_norm{2}, 'b')
        plot(mean(T_norm{2},2), mean(W_norm{2},2), 'k', 'LineWidth', 1)

    set(ax, 'Color', 'none', 'LineWidth', 1)
    linkaxes(ax(2:end),'y')
    linkaxes(ax(3:end),'x')
end
end
