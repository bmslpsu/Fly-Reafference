function [cmap] = purples(n)
%% purples: make purple colormap
% n: # levels in colormap
t = linspace(0,1,n-1)';
w = [235 220 255] / 255;
p = [100 9 200] / 255;
cmap = (1-t)*w + t*p;
cmap = [cmap ; [40 0 100]/255];
% colormap(cmap)
% pcolor(repmat(0:100,[20,1]))
% shading interp

end