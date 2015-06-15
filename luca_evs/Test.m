n=14; p=7; Imax=40; tol=1.e-14; 
% call to a test problem: 1<= prob <= 4
%prob=3; A=TestProb(n,p,prob); 
% for prob=5 need n=3; p=301.
n=3; p=301; Imax=40; tol=1.e-14; 
prob=5; A=TestProb(n,p,prob);
for alg=1:4,
   disp(sprintf('algorithm used %g',alg))   
if alg==1, [T,U]=ProdQRIter(A,tol,Imax); % Timo's algorithm 
elseif alg==2,   [T,U,PT]=AlgL1(A); % algorithm Luca 1
elseif alg==3 [T,U]=AlgL2(A,tol,Imax); % algorithm Luca 2
else  % algorithm Luca 3;  
   [T,U,PT]=AlgL3(A);
end

PA=eye(n); for j=1:p, PA=PA*A(:,:,j); end
if alg == 1, PT=eye(n); for j=1:p, PT=PT*T(:,:,j); end, end
if alg == 3, PT=eye(n); for j=1:2, PT=PT*T(:,:,j); end, end
norm(PA-U*PT*U')
% careful when looking at results of error criteria below
% because sorting etc. is confusing for complex entries
% besides, sort is not reliable in Matlab

% use Schur to get e.values of PA and of PT if alg=2 or 4
[UA,PA]=schur(PA); for i=1:n-2, PA(i,i+2:n)=0; 
   if PA(i+1,i) == 0, PA(i,i+1)=0; end; end;
if alg==2, for i=1:n-2, PT(i,i+2:n)=0; 
      if PT(i+1,i) == 0, PT(i,i+1)=0; end; end; 
elseif alg==4, for i=1:n-2, PT(i,i+2:n)=0; 
      if PT(i+1,i) == 0, PT(i,i+1)=0; end; end;
end;
ea=sort(eig(PA)); 
if alg==1, ed=sort(diag(PT)); else ed=sort(eig(PT)); end
%i=1; while i < n
%if ~isreal(ea(i)) == 1 
%   if imag(ea(i)) > 0, temp=ea(i); ea(i)=ea(i+1); ea(i+1)=temp; end
%   i=i+2; else, i=i+1; end, end
%i=1; while i < n
%if ~isreal(ed(i)) == 1 & (abs(real(ed(i))-real(ed(i+1))) <= tol ...
%   | abs(real(ed(i))-real(ed(i+1)))/abs(real(ed(i))) <= tol)
%   if imag(ed(i)) > 0, temp=ed(i); ed(i)=ed(i+1); ed(i+1)=temp; end
%   i=i+2; else, i=i+1; end, end
[ea ed]
abserr=norm(ea-ed); relerr=abserr/norm(ed);
disp(sprintf('evs-abserr=%g   evs-relerr=%g',abserr, relerr)) 
end