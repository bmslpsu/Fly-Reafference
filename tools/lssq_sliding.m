function [out] = lssq_sliding(t, x, f, win_size, overlap, n_detrend, R2_cut, showfits, showplot)
%% lssq_sliding: perform a sliding window least-squares-spectral-analysis (LSSA)
%
%   INPUTS:
%       t           : time [s], must be evenly spaced and monotonically increasing
%       x           : signal, can have Nan's
%       f           : [nx1] frequencies to fit [Hz]
%       win_size   	: window size for LSSA
%       overlap   	: window overlap
%       n_detrend   : order of polynomial to detrend with
%       showfits  	: (boolean) show fit plots
%       showplot  	: (boolean) show final magnitude & phase plot if true
%
%   OUTPUTS:
%       out         : analysis output
%

if nargin < 9
    showplot = true;
    if nargin < 8
        showfits = false;
        if nargin < 7
           n_detrend = [];
            if nargin < 6
               R2_cut = 0.7; 
            end
        end
    end
end

n_point = length(t); % # of data points
T = t(end) - t(1); % total time
ts = mean(diff(t)); % sampling time
fs = 1 / ts; % sampling rate

win_index = round(fs * win_size); % window size in indicies
%overlap_index = round(fs * overlap); % overlap size in indicies
shift = win_size - overlap; % how much to shift in time for each window
shift_index = shift*fs; % shift window index size
n_win = floor(T / shift) - 1; % total # of windows for given time

out.time = nan(n_win,1);
out.magnitude = nan(n_win,1);
out.phase = nan(n_win,1);
out.dc = nan(n_win,1);
out.R2 = nan(n_win,1);
out.fit = cell(n_win,1);

% Sliding LSSQ
for k = 1:n_win
    %disp(k)
    startI = round((1 + shift_index*(k-1))); % shift start point by overlap size
    endI = round(startI + win_index); % end point extends from start to the window size
    winI = startI:endI; % window
    
    if any(winI > n_point)
        %warning('Window too small')
        %winI = winI(winI <= n_point);
    else
        % Get data in window to fit
        t_win = t(winI);
        x_win = x(winI);
        
        % Detrend window
        if ~isempty(n_detrend)
            x_win = detrend(x_win, n_detrend);
        end
        
        % Fit data
        fit_out = fit_sine(t_win, x_win, f, showfits);
        
        % Store fit data
        out.time(k) = mean(t_win);
        out.fit{k} = fit_out.fit; 

        % Make sure fit is good enough
        if fit_out.R2 < R2_cut
            warning(['R^2 under ' num2str(R2_cut) ': skipping '])
        else
            out.magnitude(k) = fit_out.mag;
            out.phase(k) = fit_out.phase;
            out.dc(k) = fit_out.dc;
            out.R2(k) = fit_out.R2;
        end
        
        if showfits
            pause
            clf
        end
    end
end

if showplot
    ax(1) = subplot(4,1,1); cla ; hold on
        plot(out.time, out.magnitude, 'ob-', 'LineWidth', 1)
        ylim([0 5*ceil(max(out.magnitude)/5)])
        ylabel('magnitude')
    ax(2) = subplot(4,1,2); cla ; hold on
        plot(out.time, rad2deg(out.phase), 'or-', 'LineWidth', 1)
        ylim(180*[-1 1])
        ylabel('phase (rad)')
    ax(3) = subplot(4,1,3); cla ; hold on
        plot(out.time, out.dc, 'og-', 'LineWidth', 1)
        %ylim([-pi pi])
        ylabel('DC')
    ax(4) = subplot(4,1,4); cla ; hold on
        plot(out.time, out.R2, 'ok-', 'LineWidth', 1)
        ylim([0 1])
        ylabel('R^2')
        xlabel('time (s)')
        
	linkaxes(ax, 'x')
    set(ax, 'LineWidth', 1, 'XLim', [0 t(end)])
end

end