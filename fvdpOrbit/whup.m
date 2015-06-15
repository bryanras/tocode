function [xx,nb]=whup(pts)

%
% Usage: [xx,nb] = whup(pts)
%
% This is just an initialization function.
%
% The initial orbit and vectors are really simple.
% 

% Radius.
rad=5;

% Angles.
phi = (0:pts-1)'*2*pi/pts;

% Points.
xx=rad*[cos(phi),sin(phi),zeros(pts,1)];

% Normal vectors. One will always be straight up.
nb=[zeros(pts,2),ones(pts,1),xx/rad]; 
 
