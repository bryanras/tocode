function N=locn(xll,xlr,xul,xur)

% This function calculates a normal vector to a rectangle of points in R^3.

% Get the tangent vectors.
xtani = xur-xll;
xtanj = xul-xlr;

% Assign a normal vector.
N = cross(xtani,xtanj);
N=N/norm(N);

