function [] = SS_example()
%% SS_example:


f = 0.5;
A = 37.5;
phi = 0;
fs = 100;
ts = 1 / fs;
T = 1/f;
tt = 3*T;

t = (0:ts:tt)';
r = A*sin(2*pi*f*t + phi);

gain = 0.9;
phase = -15;

y = gain*A*sin(2*pi*f*t + phi + deg2rad(phase));
e = r - y;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 4 2])
movegui(fig, 'center')

ax = subplot(1,1,1) ; cla ; hold on
    plot(t, r, 'k', 'LineWidth', 1)
    plot(t, y, 'r', 'LineWidth', 1)
    plot(t, e, 'b', 'LineWidth', 1)
    
    set(ax, 'Color', 'none', 'XLim', [-0.05*tt 1.05*tt])

end