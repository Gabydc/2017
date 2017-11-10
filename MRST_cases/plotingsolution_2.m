function [h1]=plotingsolution_2(G,W,name,x,n)
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% clf;
% Plot the grid
subplot(1,2,n)
h1=plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
plotCellData(G,x)

% Plot the wells
% plotWell(G, W);
view(30,50)

axis equal tight; colormap(jet(128));
title(name);
