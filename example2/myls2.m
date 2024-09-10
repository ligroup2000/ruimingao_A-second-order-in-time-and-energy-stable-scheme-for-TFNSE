function x=myls2(A,b,eps)
% solve the rank deficient least squares problem by SVD
% x is the LS solution, res is the residue

[m,n]=size(A);
[Q,R]=qr(A,0);

if nargin < 3
    eps=1e-13;
end
s=diag(R);

%figure
%semilogy(s,'r.')
r=sum(abs(s)>eps);

Q = Q(:, 1:r);
R = R(1:r,1:r);
b1 = b(r+1:m+r);

x= R\(Q.'*b1);

%res = norm(A(:,1:r)*x-b1)/norm(b1)
end