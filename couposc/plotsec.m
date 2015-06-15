function plotsec(huge,pts,isec)

%
% Just a little plotter.  Feel free to delete.
%
% Usage: plotsec(huse,pts,isec), where isec is any phi_1 = constant section.
%

% Section, start.
ss=(isec-1)*pts(2)+1;

seccy=[ss:isec*pts(2),ss];

plot(huge(seccy,3),huge(seccy,4))
hold on

plot(huge(seccy,11),huge(seccy,12))
plot(huge(seccy,19),huge(seccy,20))
plot(huge(seccy,43),huge(seccy,44))
plot(huge(seccy,63),huge(seccy,64))
plot(huge(seccy,79),huge(seccy,80))
plot(huge(seccy,83),huge(seccy,84))

hold off

