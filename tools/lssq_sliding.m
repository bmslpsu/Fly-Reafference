function [out] = lssq_sliding(t, x, f, win_size, overlap, showplot)
%% lssq_sliding: perform a sliding window least-squares-spectral-analysis (LSSA)
%
%   INPUTS:
%       t           : time [s], must be evenly spaced and monotonically increasing
%       x           : signal, can have Nan's
%       f           : [nx1] frequencies to fit [Hz]
%       win_size   	: window size for LSSA
%       overlap   	: window overlap
%       showplot  	: (boolean) show plot if true
%
%   OUTPUTS:
%       out         : analysis output
%

if nargin < 6
    showplot = true;
end

n_point = length(t); % # of data points
T = t(end) - t(1); % total time
ts = mean(diff(t)); % sampling time
fs = 1 / ts; % sampling rate

win_index = round(fs * win_size); % window size in indicies
overlap_index = round(fs * overlap); % overlap size in indicies
shift = win_size - overlap; % how much to shift in time for each window
shift_index = shift*fs; % shift window index size
n_win = floor(T / shift); % total # of windows for given time

out.time = nan(n_win,1);
out.magnitude = nan(n_win,1);
out.phase = nan(n_win,1);
out.fit = cell(n_win,1);

% Sliding LSSQ
for k = 1:n_win
    %disp(k)
    startI = (1 + shift_index*(k-1)); % shift start point by overlap size
    endI = startI + win_index; % end point extends from start to the window size
    winI = startI:endI; % window
    
    if any(winI > n_point)
        warning('Window too small')
        %winI = winI(winI <= n_point);
    else
        t_win = t(winI);
        x_win = x(winI);
        fit_out = fit_sine(t_win, x_win, f, false);
        out.time(k) = mean(t_win);
        out.magnitude(k) = fit_out.mag;
        out.phase(k) = fit_out.phase;
        out.fit{k} = fit_out.fit;
    end
end

if showplot
    ax(1) = subplot(2,1,1); cla ; hold on
        plot(out.time, out.magnitude, 'ob-', 'LineWidth', 1)
        ylim([0 5*ceil(max(out.magnitude)/5)])
    ax(2) = subplot(2,1,2); cla ; hold on
        plot(out.time, out.phase, 'or-', 'LineWidth', 1)
        ylim([-pi pi])
        
	linkaxes(ax, 'x')
    set(ax, 'XLim', [0 t(end)])
end

end