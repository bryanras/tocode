function []=genFigs

%
% Usage: genFigs
%
%
% Script for generating some figures that I use in a paper.
% Just a place to store some commands.
%


% The following are the lambdas at which we want to form plots.
%lam1: Periodic orbit curves, w = sqrt(0.84), all on same plot
%lam2: Tori, w = sqrt(0.84), one plot per lambda.
%lam3: Periodic orbit curves, w = sqrt(0.78), all on same plot
%lam4: Tori, w = sqrt(0.78), one plot per lambda.
lam1 = [0,0.1,0.2,0.3,0.34,0.3410,0.3415]';
lam2 = [0.3,0.3415]';
lam3 = [0,0.1,0.2,0.3,0.35,0.38,0.39,0.3945]';
lam4 = [0.38,0.388]';

% This just gets rid of a stupid mLint error.
dvec = [1,1]; huge = [1,1];

% Add plotting function. There has to be a way to do this cleanly.
ccd = pwd;
cd('../fvdpOrbit');
fcnh=@gencurves;
cd(ccd);

% Use this to keep all the figure handles.
fh = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% w = sqrt(0.84);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load orbits.
load ../fvdpOrbit/55ptSqrt84;

% Rename and extract.
hugeo = huge;
dveco = dvec(dvec(:,3)==-1,1);

% Go ahead and plot the orbits while we're in the correct directory.
% Get the lambda values from the diary vector.
idx = interp1(dveco,(1:length(dveco))',lam1,'nearest');

% Check to see that those lambda values were actually in
% the vector.
if max(abs(lam1-dveco(idx)))>10*eps
	disp(['Error: Couldn''t find all lambdas in periodic', ...
			' orbit data, w=sqrt(0.84)']);
end

% Plot them regardless.
cls = reshape([3*idx-2,3*idx-1,3*idx]',length(idx)*3,1);
feval(fcnh,hugeo(:,cls),0);

%%%%%%%%%%%% Prettify the plot.
a1= gca;
grid on;
set(a1, 'FontSize', 12);
set(a1,'XLim', [-6, 6], 'YLim', [-6, 6], 'ZLim', [-2.5, 2.5]);
set(a1, 'XTick', [-5 0 5], 'YTick', [-5 0 5], 'ZTick', [-2 0 2]);
set(a1, 'XTickLabel', {'','0','5'}, ...
	'YTickLabel', {'','0','5'}, 'ZTickLabel', [-2 0 2]);
set(a1, 'Projection', 'perspective');
set(a1,'Position',[0.07,0.05,0.9,0.94]);

fh(end+1) = gcf;
add3dx(fh(end));

% Create textarrow
annotation(...
  fh(end),'textarrow',...
  [0.2518 0.3143],[0.6667 0.5952],...
  'String',{'\lambda=0'},...
  'FontName','helvetica',...
  'FontSize',12);
 
% Create textarrow
annotation(...
  fh(end),'textarrow',...
  [0.4464286, 0.307453],[0.228609523809524, 0.326661618045035],...
  'String',{'\lambda=0.3415'},...
  'FontName','helvetica',...
  'FontSize',12);
 

%%%%%% End Prettification.


% Go back to original directory and load tori.
load 45x45Sqrt84;

% Clean up.
dvec = dvec(dvec(:,end)==-1,1);
hugeo = [hugeo;hugeo(1,:)];

% Get the indices corresponding to the requested lambdas.
idx = interp1(dvec,(1:length(dvec))',lam2,'nearest');
idxo = interp1(dveco,(1:length(dveco))',lam2,'nearest');

% Check to make sure they are there.
if any(abs(lam2-dveco(idxo))>10*eps)
	disp(['Error: Couldn''t find all lambdas in periodic', ...
		' orbit data, w=sqrt(0.84)']);
end
if any(abs(lam2-dvec(idx))>10*eps)
	disp(['Error: Couldn''t find all lambdas in torus', ...
		' data, w=sqrt(0.84)']);
end
	
	
% Different plot for each lambda in lam2.
for ii=1:length(idx)

	fh(end+1)=figure;

	% Get the section.
	scn = (3*idx(ii)-2:3*idx(ii));
	scno = (3*idxo(ii)-2:3*idxo(ii));

	% Plot torus.
	han = torit(huge(:,scn),45,1);

	% Plot line.
	lh = line(hugeo(:,scno(1)),hugeo(:,scno(2)),hugeo(:,scno(3)));

	% Prettify.
	colormap([0 0 0]);
	set(han,'EdgeColor',[0.3,0.3,0.3],'FaceAlpha',0.5);
	set(lh,'LineStyle','-.','LineWidth',3,'Color','b');

	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% w = sqrt(0.87);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load orbits.
load ../fvdpOrbit/55ptSqrt78;

% Rename and extract.
hugeo = huge;
dveco = dvec(dvec(:,3)==-1,1);

% Go ahead and plot the orbits while we're in the correct directory.
% Get the lambda values from the diary vector.
idx = interp1(dveco,(1:length(dveco))',lam3,'nearest');

% Check to see that those lambda values were actually in
% the vector.
if max(abs(lam3-dveco(idx)))>10*eps
	disp(['Error: Couldn''t find all lambdas in periodic', ...
			' orbit data, w=sqrt(0.78)']);
end

% Plot them regardless.
cls = reshape([3*idx-2,3*idx-1,3*idx]',length(idx)*3,1);
feval(fcnh,hugeo(:,cls),0);

%%%%%%%%%%%%%%%%%% Prettify the plot.
a1= gca;
grid on;
set(a1, 'FontSize', 12);
set(a1,'XLim', [-6, 6], 'YLim', [-6, 6], 'ZLim', [-2.5, 2.5]);
set(a1, 'XTick', [-5 0 5], 'YTick', [-5 0 5], 'ZTick', [-2 0 2]);
set(a1, 'XTickLabel', {'','0','5'}, ...
	'YTickLabel', {'','0','5'}, 'ZTickLabel', [-2 0 2]);
set(a1, 'Projection', 'perspective');
set(a1,'Position',[0.07,0.05,0.9,0.94]);

fh(end+1) = gcf;
add3dx(fh(end));

% Create textarrow
annotation(...
  fh(end),'textarrow',...
  [0.2465 0.309],[0.6691 0.5976],...
  'String',{'\lambda=0'},...
  'FontName','helvetica',...
  'FontSize',12);
 
% Create textarrow
annotation(...
  fh(end),'textarrow',...
  [0.5375 0.3857],[0.181 0.2881],...
  'String',{'\lambda=0.3945'},...
  'FontName','helvetica',...
  'FontSize',12);
 

%%%%%%%%%%%%%%%%%%% End Prettification.


% Go back to original directory and load tori.
load 105x45Sqrt78;

% Clean up.
dvec = dvec(dvec(:,end)==-1,1);
hugeo = [hugeo;hugeo(1,:)];

% Get the indices corresponding to the requested lambdas.
idx = interp1(dvec,(1:length(dvec))',lam4,'nearest');
idxo = interp1(dveco,(1:length(dveco))',lam4,'nearest');

% Check to make sure they are there.
if any(abs(lam4-dveco(idxo))>10*eps)
	disp(['Error: Couldn''t find all lambdas in periodic', ...
		' orbit data, w=sqrt(0.78)']);
end
if any(abs(lam4-dvec(idx))>10*eps)
	disp(['Error: Couldn''t find all lambdas in torus', ...
		' data, w=sqrt(0.78)']);
end
	
	
% Different plot for each lambda in lam4.
for ii=1:length(idx)

	fh(end+1)=figure;

	% Get the section.
	scn = (3*idx(ii)-2:3*idx(ii));
	scno = (3*idxo(ii)-2:3*idxo(ii));

	% Plot torus.
	han = torit(huge(:,scn),45,1);

	% Plot line.
	lh = line(hugeo(:,scno(1)),hugeo(:,scno(2)),hugeo(:,scno(3)));

	% Prettify.
	colormap([0 0 0]);
	set(han,'EdgeColor',[0.3,0.3,0.3],'FaceAlpha',0.5);
	set(lh,'LineStyle','-.','LineWidth',3,'Color','b');
	
end

% These lines save everything. Be careful.
%saveas(fh(1),'fvdpOrbits84','epsc');
%saveas(fh(2),'fvdp843','epsc');
%saveas(fh(3),'fvdp843415','epsc');
%saveas(fh(4),'fvdpOrbits78','epsc');
%saveas(fh(5),'fvdp78380','epsc');
%saveas(fh(6),'fvdp78388','epsc');

