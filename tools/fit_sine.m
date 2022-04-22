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

% Remove Nan's
nanI = isnan(x);
t = t(~nanI);
x = x(~nanI);

% Define function to fit
fit = @(b,x)  b(1).*(sin(2*pi*f*x + b(2))) + b(3);

% Get approximant values for starting points
yu = max(x);
yl = min(x);
yr = (yu - yl) / 2;
ym = mean(x);
x0 = [yr ; 0 ; ym];

% Define cost function
fcn = @(b) sum((fit(b,t) - x).^2);

% Minimize cost fucntion and solve for coefficents
opt = optimset('MaxIter', 1500);
s = fminsearch(fcn, x0, opt);

% Get magnitude, phase, & DC
out.mag = s(1);
out.phase = wrapToPi(s(2));
out.dc = s(3);

% Correct for negative gain fits
if out.mag < 0
    out.mag = -out.mag;
    out.phase = out.phase - pi;
end

% Calculate fit accuracy with R^2
x_fit = fit(s, t);
R = corrcoef(x, x_fit);
R2 = R(1,2).^(2);
out.R2 = R2;

% Create fit curve
ts = mean(diff(t));
xp = min(t):(ts / 4):max(t);
out.x = xp;
out.fit = fit(s,xp);

% Plot
if showplot
    hold on
    title(['R^2 = ' num2str(R2)])
    plot(t, x, '.k')
    plot(out.x, out.fit, 'r', 'LineWidth', 1)
    grid
end

end