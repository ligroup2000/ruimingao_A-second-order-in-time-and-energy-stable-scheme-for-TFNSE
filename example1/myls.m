function [x,res]=myls(A,b,eps)
% solve the rank deficient least squares problem by SVD
% x is the LS solution, res is the residue

[m,n]=size(A);
[U,S,V]=svd(A,0);

if nargin < 3
    eps=1e-12;
end
s=diag(S);
%figure
%semilogy(s,'r.')
r=sum(s>eps);


x=zeros(n,1);
for i=1:r
    x=x+(U(:,i)'*b)/s(i)*V(:,i);
end

if (nargout>1)
%     res=0;
%     for i=r+1:n
%         res=res+(U(:,i)'*b)^2;
%     end
%     res=sqrt(res);
    res = norm(A*x-b)/norm(b);
end


end