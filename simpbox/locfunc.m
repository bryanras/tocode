function G=locfunc(xxl,lam)

% This is a local function using the box scheme.

% The rows of xxl are xx+ and xx respectively.

xtan=xxl(1,:)-xxl(2,:);

nn=[xtan(2), -xtan(1)];
nn=nn/norm(nn);

G= nn*funcy(0.5*(xxl(1,:)+xxl(2,:)),lam);


