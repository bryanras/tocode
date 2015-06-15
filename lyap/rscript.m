f='vdp'; Df='Dvdp'; v0=[1;0];
tol=1.e-8; maxit=10;
N=pts; h=0.01;

% Set up the original guess.
load sguess;

% Do some spline interpolation.
ss(size(ss,1)+1,:) = ss(1,:);
newref = 1+(size(ss,1)-1)*[0:(N-1)]/N;

s = interp1(ss,newref,'cubic');
s=s';


% Initial guess for splines.
tau = 6.66;


[s,tau,YY,flag]=newtms(s,tau,N,h,v0,f,Df,tol,maxit);
M=full(monod(YY));
eig(M)
