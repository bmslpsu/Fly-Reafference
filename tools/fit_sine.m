function [out] = fit_sine(t, x, f, showplot)
%% fit_sine: 
%
%   INPUTS:
%       x           : fit data
%       t           : time data
%       f           : frequency
%       showplot  	: (boolean) show debug plot if true
%
%   OUTPUTS:
%       out         : fit structure
%

if nargin < 4
    showplot = true;
end

nanI = isnan(x);
t = t(~nanI);
x = x(~nanI);

yu = max(x);
yl = min(x);
yr = (yu - yl);
% ya = sqrt(yr);
% yz = x-yu+(yr/2);
% ym = mean(x);
fit = @(b,x)  b(1).*(sin(2*pi*f*x + b(2)));
% fit = @(b,x)  b(1).*cos(2*pi*f*x) + b(2).*cos(2*pi*f*x) + b(3);
fcn = @(b) sum((fit(b,t) - x).^2);
s = fminsearch(fcn, [yr; 0]);

ts = mean(diff(t));
xp = min(t):(ts / 4):max(t);

% if s(1) >=0 
%     out.mag = s(1)';
%     out.phase = s(2);
% else
%     out.mag = -s(1)';
%     out.phase = s(2) + pi;
% end

% get magnitude and phase
out.mag = s(1);
out.phase = wrapToPi(s(2));

% out.dc = s(3);
out.x = xp;
out.fit = fit(s,xp);

% out.gain = sqrt(s(1)^(2) + s(2)^(2));
% out.phase = rad2deg( atan2(s(2), s(1)) );
% out.dc = s(3);
if showplot
    figure
    hold on
    plot(t, x, '.b')
    plot(out.x, out.fit, 'r', 'LineWidth', 1)
    grid
end

end