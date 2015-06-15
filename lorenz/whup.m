function [xx,nb]=whup(pts)

%
% Usage: [xx,nb] = whup(pts)
%
% This is just an initialization function.
% It works on an interpolation premise, and it requires the 
% approximation in the file XXpt.
%
% The normals come from a process as follows:
% 
% 1. Get the tangents through center-difference approxmation.
% 2. Project sdvec=[1,0,0] onto each normal and normalize. Thats nb(ii,1:3);
% 3. Cross the tangents for nb(ii,4:6)
% 

% Load up the initial approxmation.
xx=load('200pt');
pty = size(xx,1);			% This should be XX of XXpt fame.
xx(pty+1,:)=xx(1,:);

% Get the arc length between points.
al(1) = 0;
for ii=2:pty+1 al(ii)=al(ii-1)+norm(xx(ii,:)-xx(ii-1,:)); end
al = al/al(pty+1);

% Now we interpolate.
xx=interp1(al,xx,0:1/pts:(1-1/pts));


% OK, now here's the hard part.  We need to get the normal vectors.
% I see no way to avoid a loop.  Oh, well.  This isn't the most intensive
% part of the code anyway.
nb=zeros(pts,6);
sdvec=[0,0,1];
for ii=1:pts

	ip = mod(ii,pts)+1; im = mod(ii-2,pts)+1;      % Forward and backward.

	% Get the tangents.
	tn = xx(ip,:)-xx(im,:);
	tn=tn/norm(tn);

	% Project. I could write this more efficiently, but I want to preserve
	% the philosophy in the code.
	nb(ii,1:3) = sdvec-tn*(tn*sdvec');

	if norm(nb(ii,1:3))<0.05 disp('damn'); end 		% Just a warning.

	nb(ii,1:3) = nb(ii,1:3)/norm(nb(ii,1:3));

	% Get the second normal vector.
	nb(ii,4:6) = cross(nb(ii,1:3),tn);
	nb(ii,4:6) = nb(ii,4:6)/norm(nb(ii,4:6));
	
end
