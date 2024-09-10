function x = myls_qr(A,b,eps)
%myls_qr 此处提供此函数的摘要
%   此处提供详细说明
[m,~] = size(A);
[Q,R] = qr(A,0);
if nargin<3
    eps = 1e-13;
end
s = diag(R);
r = sum(abs(s)>eps);
Q = Q(:,1:r);
R = R(1:r,1:r);
b1 = b(r+1:m+r);
x =R\(Q.'*b1);
end
