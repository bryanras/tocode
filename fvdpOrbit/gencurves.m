function [] = gencurves(huge,lbl)

%
% Usage: gencurves(huge);
%
% Yet another figure generator.
%
% Set lbl = 0 to suppress x,y,z, labels
%

if nargin < 2 || isempty(lbl)
	lbl = 1;
end

% Tack on a last row.
huge(end+1,:) = huge(1,:);

% Define colors.
colors = 'kcrgmb';
nc = length(colors);

% Create figure and add the first line.
figure;
view([322.5 30]);

% Now add lines.
for ii=1:size(huge,2)/3

	% Get the color.
	clr = mod(ii-1,nc)+1;
	line(huge(:,3*ii-2),huge(:,3*ii-1),huge(:,3*ii), ...
			'Color',colors(clr))

end

if lbl
	xlabel('x_1','FontSize',12);
	ylabel('x_2','FontSize',12);
	zlabel('x_3','FontSize',12,'Rotation',0);
end
