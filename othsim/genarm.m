function genarm(pts,crd,com,fx,xx)

% This is a function for plotting trunk plots of a three-torus.
% 
% Usage:
%
%	genarm(pts,crd,com,fx,xx)
%
% 
% xx = points themselves.
% crd = [coordinate1 coordinate2]
% com = [horizontal component, vertical component]
% fx  = point at which to fix the third angular component
%
% For example, to plot an i-section with increasing
% k on the vertical axis, the value of j must be fixed, say at 1:
%
%        genarm(pts,[1,2],[1,3],1,xx)
%
%
% Be careful -- this function does no error checking.
%
% If you want to make the mesh see-though, type "hidden off".
% If you want to remove the grid lines, type "grid off".

% Get the fixed coordinate.
com(3) = 1+2+3 - sum(com);
pts(length(pts)+1)=1;

% Here is the increment along each seciton, plus 
% the bump value between sections, and the increment used when starting.
inc 	= prod( pts(com(1)+1:length(pts)) );
bumpit	= prod( pts(com(2)+1:length(pts)) );
sinc 	= prod( pts(com(3)+1:length(pts)) );

% Here is the start value.
start = (fx-1)*sinc+1;

% Here is the first section with closed circle.
scn = start:inc:start+inc*(pts(com(1))-1);
scn(length(scn)+1) = scn(1);
		
% Q: How many sections are there going to be? A: pts(com(2)) (+1)
for ii=1:pts(com(2))
	
	% Set the other sections.
	ZZ(:,ii)=ii*ones(length(scn),1);
	XX(:,ii)=xx(scn,crd(1)); YY(:,ii)=xx(scn,crd(2));

	% Update.
	scn=scn+bumpit;

end

% Tack on the final section at the top.
ZZ(:,pts(com(2))+1)=ZZ(:,pts(com(2)))+1;
XX(:,pts(com(2))+1)=XX(:,1); 	YY(:,pts(com(2))+1)=YY(:,1);

% Plot it.
mesh(XX,YY,ZZ);

% Set some defaults.
colormap([0,0,0]);
a1=gca;
set(a1, 'FontSize', [12]);
set(a1, 'XTick', [-1 0 1], 'YTick', [-1 0 1], 'ZTick', []);
set(a1, 'XTickLabel', [-1 0 1], 'YTickLabel', [-1 0 1], 'ZTickLabel', []);
set(a1, 'Projection', 'perspective');
set(a1, 'XLim', [-1.5,1.5], 'YLim', [-1.5,1.5]);
pbaspect([1 1 1.8])

% Add the labels.
xlab=sprintf('x_{%d}',mod(crd(1)-1,8)+1);
ylab=sprintf('x_{%d}',mod(crd(2)-1,8)+1);
philab=sprintf('\\phi_%d',com(2));
h1=xlabel(xlab); h2=ylabel(ylab); h3=zlabel(philab);

% Set the axis label defaults.
set([h1,h2,h3],'FontSize',[12],'FontWeight','bold')
set(h3,'Rotation',[0]);


