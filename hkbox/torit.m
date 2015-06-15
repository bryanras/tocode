function torit(xx,ptsi,ptsj,typ,open)

% This function plots a pretty picture of a torus.
%
% The usage is torit(xx,ptsi,ptsj,typ,open), where "typ" determines if the 
% torus is a mesh (typ=1) or a surface (typ=2).  See "colormap".
%
% If you want to make the mesh see-though, type "hidden off".
% If you want to remove the grid lines, type "grid off".
%
% If you want to leave the torus open at the ends, use open=1.
%
% Default is typ = 1, open = 0.
%

% Set some defaults.
if ~exist('open','var')
    open=1;
end
if ~exist('typ','var')
    typ = 1;
end


% Generate a mesh.
XX=reshape(xx(:,1),ptsj,ptsi);
YY=reshape(xx(:,2),ptsj,ptsi);	
ZZ=reshape(xx(:,3),ptsj,ptsi);	

% Add on the last bit.
XX(ptsj+1,:)=XX(1,:);
YY(ptsj+1,:)=YY(1,:);
ZZ(ptsj+1,:)=ZZ(1,:);
if open ~= 1
	XX(:,ptsi+1)=XX(:,1);
	YY(:,ptsi+1)=YY(:,1);
	ZZ(:,ptsi+1)=ZZ(:,1);
end

if typ==2 
	colormap('default');
	surf(XX,YY,ZZ);
else
	mesh(XX,YY,ZZ);
	colormap([0,0,0]);
end

% Set some axis defaults.
%axis([-6, 6, -6, 6, -2.5, 2.5]);
axis([-1.2, 1.2, -1.2, 1.2, 0, 1.5]);  % For a blow-up.

% Get the goddamn axis handle.
a1=gca;
set(a1, 'FontSize', [12]);
set(a1, 'XTick', [-4 0 4], 'YTick', [-4 0 4], 'ZTick', [-1 0 1]);
set(a1, 'XTickLabel', [-4 0 4], 'YTickLabel', [-4 0 4], 'ZTickLabel', [-1 0 1]);
set(a1, 'Projection', 'perspective');

% Add the labels.
h1=xlabel('x_1'); h2=ylabel('x_2'); h3=zlabel('x_3');

% Set the axis label defaults.
set([h1,h2,h3],'FontSize',[12],'FontWeight','bold')
set(h3,'Rotation',[0]);

