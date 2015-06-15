function A=TestProb(n,p,prob)

if prob == 1
   % test problem #1
   A=randn(n,n,p);
elseif prob == 2
   % test problem #2
   A=zeros(n,n,p);
   b=randn(n,n); [q,r]=qr(b);
   for i=1:p, 
   if mod(i,2), 
      for j=1:n, A(j,j,i)=j; end, 
   else, 
%      for j=1:n, A(j,j,i)=1/j; end; 
      for j=1:n, A(j,j,i)=10^(-j); end; 
   end
   A(:,:,i)=q'*A(:,:,i)*q;   end
elseif prob == 3
   % test problem #3
   A=randn(n,n,p); 
   %for i=1:p, A(:,:,i)=triu(A(:,:,i)); end
   %b=randn(n,n); [q,r]=qr(b);
   for i=1:p, [q,r]=qr(A(:,:,i)); A(:,:,i)=r; end;
   for i=1:p, 
   if mod(i,2), 
      for j=1:n, A(j,j,i)=j; end, 
   else, 
      for j=1:n, A(j,j,i)=1/j; end; 
   end
   A(:,:,i)=q*A(:,:,i)*q';   end
elseif prob == 4
   % test problem #4
   A=randn(n,n,p); B=randn(n,n,p); A=A+sqrt(-1)*B;
else 
   % test problem #5, n=3, p=301
   as=[1 1 1; 0 -0.9 1; 0 0 1.1]; 
   q=1/3*[2 -2 1; 2 1 -2; 1 2 2]; b=q*as*q';
   for i=2:p-1, A(:,:,i)=b; end;
   A(:,:,p)=q*as; A(:,:,1)=as*q';
end