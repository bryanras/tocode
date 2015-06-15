function han = torit(xx,ptsj,op)

%
% Usage: han = torit(xx,ptsj,op)
%
% This function plots a pretty picture of a torus.
%
% ptsj is the number of points along the meridian.
%
% If you want to produce a "cut-away" torus, use op=1. Default is op=0.
%
% han is the handle to the object.
%

% Get defaults.
if nargin < 3 || isempty(op)
	op = false;
end

% Get the number of meridia.
ptsi = round(size(xx,1)/ptsj);

% Get the slice to plot. If the user wants a cut-away, plot a piece of it.
if op
	slc=[ round(22*ptsi/45):ptsi, 1:round(12*ptsi/45)]; 
else
	slc = [1:ptsi,1];
end

% Generate a mesh.
XX=reshape(xx(:,1),ptsj,ptsi);
YY=reshape(xx(:,2),ptsj,ptsi);	
ZZ=reshape(xx(:,3),ptsj,ptsi);	


% Close the meridia. Possibly remove a slice.
XX(ptsj+1,:)=XX(1,:);
YY(ptsj+1,:)=YY(1,:);
ZZ(ptsj+1,:)=ZZ(1,:);
XX = XX(:,slc);
YY = YY(:,slc);
ZZ = ZZ(:,slc);

% if ~op
% 	XX(:,ptsi+1)=XX(:,1);
% 	YY(:,ptsi+1)=YY(:,1);
% 	ZZ(:,ptsi+1)=ZZ(:,1);
% end


% Here's the plot.
han = mesh(XX,YY,ZZ);
colormap([0,0,0]);

% Set some axis defaults.
axis([-6, 6, -6, 6, -2.5, 2.5]);

% Get the axis handle.
a1=gca;
set(a1, 'FontSize', 12);
set(a1,'XLim', [-6, 6], 'YLim', [-6, 6], 'ZLim', [-2.5, 2.5]);
set(a1, 'XTick', [-5 0 5], 'YTick', [-5 0 5], 'ZTick', [-2 0 2]);
set(a1, 'XTickLabel', {'','0','5'}, ...
	'YTickLabel', {'','0','5'}, 'ZTickLabel', [-2 0 2]);
set(a1, 'Projection', 'perspective');
set(a1,'Position',[0.07,0.05,0.9,0.94]);


% Add labels.
add3dx(gcf);

if nargout < 1
	clear han;
end
