function [out] = fit_line(x, y, showplot)
%% fit_line: 
%
%   INPUTS:
%       x           : x-data
%       y           : y-data
%       showplot  	: (boolean) show debug plot if true
%
%   OUTPUTS:
%       out         : fit structure
%

if nargin < 3
    showplot = true;
end

% Remove Nan's
nanI = isnan(y);
x = x(~nanI);
y = y(~nanI);

% Define function to fit
fit = @(b,x)  b(2) + b(1).*x; % linear 1st order polynomial

% Define cost function
fcn = @(b) sum((fit(b,x) - y).^2);

% Minimize cost fucntion and solve for coefficents
opt = optimset('MaxIter', 1500);
x0 = [0 0];
s = fminsearch(fcn, x0, opt);

% Get slope and y-intercept
out.slope = s(1);
out.intercept = s(2);

% Calculate fit accuracy with R^2
x_fit = fit(s, x);
[R,P] = corrcoef(y, x_fit);
R2 = R(1,2).^(2);
out.R2 = R2;
out.P = P(1,2);

% Create fit curve
% ts = mean(diff(x));
% xp = min(x):max(x);
out.x = x;
out.fit = fit(s,out.x);

% Plot
if showplot
    hold on
    title(['R^2 = ' num2str(R2)])
    plot(x, y, '.k')
    plot(out.x, out.fit, 'r', 'LineWidth', 1)
    grid
end

end