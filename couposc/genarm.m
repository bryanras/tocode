function genarm(xx,pts,crd,vert)

% This is a function for plotting trunk plots of a two-torus in R^4.
% 
% Usage:
%
%	genarm(xx,pts,crd,vert)
%
% 
% xx = points themselves.
% crd = [coordinate1 coordinate2]
% vert = vertical component -- either 1 or 2, for phi1 or phi2
%
% For example, to plot an i-section, starting at 1, with increasing
% phi2 on the vertical axis, write 
% genarm(xx,pts,[1,2],2)
%
% A zero in the vertical component creates a planar plot.
%
% This function does no error checking, so your start point had better be right.
% If you want to make the mesh see-though, type "hidden off".
% If you want to remove the grid lines, type "grid off".

% I'll make this as simple as possible.
if vert==2
	scn=[1:pts(2):((pts(1)-1)*pts(2)+1),1];
elseif vert==1
	scn=[1:pts(2),1];
else
	disp(' Man, you''re stupid.');
	return;
end
		
% More conditionals, but I don't feel like being slick today.
if vert==1
	bumpit=pts(2);
else 
	bumpit=1;
end


% Q: How many sections are there going to be? A: pts(vert)+1
for ii=1:pts(vert)
	
	% Set the other sections.
	ZZ(:,ii)=ii*ones(length(scn),1);
	XX(:,ii)=xx(scn,crd(1)); YY(:,ii)=xx(scn,crd(2));

	% Update.
	scn=scn+bumpit;
end

% Tack on the final section at the top.
ZZ(:,pts(vert)+1)=ZZ(:,pts(vert))+1;
XX(:,pts(vert)+1)=XX(:,1); 	YY(:,pts(vert)+1)=YY(:,1);

% Plot it.
a2=mesh(XX,YY,ZZ);

% Set some defaults.
a1=gca;
colormap([0,0,0]);
set(a1, 'FontSize', [12]);
set(a1, 'XTick', [-1 0 1], 'YTick', [-1 0 1], 'ZTick', []);
%set(a1, 'XTickLabel', [-1 0 1], 'YTickLabel', [-1 0 1], 'ZTickLabel', []);
set(a1, 'Projection', 'perspective');
set(a1, 'XLim', [-1,1], 'YLim', [-1,1]);
pbaspect([1 1 1.8])

% Add the labels.
xlab=sprintf('x_{%d}',mod(crd(1)-1,4)+1);
ylab=sprintf('x_{%d}',mod(crd(2)-1,4)+1);
philab=sprintf('\\phi_%d',vert);
h1=xlabel(xlab); h2=ylabel(ylab); h3=zlabel(philab);

% Set the axis label defaults.
set([h1,h2,h3],'FontSize',[12],'FontWeight','bold');
set(h3,'Rotation',[0]);

% This is for all the plane-intersection shit.  See 8/14/03.
if 0 == 0
	set(a2, 'FaceAlpha', [0.5], 'EdgeColor', [0.7,0.7,0.7] );

	% Get the intersection lines.
	[y1,y2] = harint(xx,pts);
	y1=sortrows(y1,vert); y2=sortrows(y2,vert);

	hold on
	a3 = plot3( y1(:,crd(1)+2), y1(:,crd(2)+2), y1(:,vert));
	a4 = plot3( y2(:,crd(1)+2), y2(:,crd(2)+2), y2(:,vert));
	set(a3, 'LineWidth', [3.0]); 
	set(a4, 'LineWidth', [3.0], 'LineStyle', '-.', 'Color', [0 0.5 0]);
	hold off

end
