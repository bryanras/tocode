
% This is a file for generating a figure and adding certain properties.


% Generate the picture.
plot(huge(:,1),huge(:,2),'+-',huge(:,ah*2-1),huge(:,ah*2),'+-');

a1=gca;
set(a1, 'FontSize', [12]);
set(a1, 'XTick', [-2 0 2], 'YTick', [-2 0 2] );
set(a1, 'Xlim', [-3 3], 'Ylim', [-4 4] );
set(a1, 'XTickLabel', [-2 0 2], 'YTickLabel', [-2 0 2]);

title(sprintf('\\lambda=%g',lams(ah)));
