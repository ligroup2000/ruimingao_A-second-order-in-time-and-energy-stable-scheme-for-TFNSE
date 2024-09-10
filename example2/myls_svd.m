function [x,res] = myls_svd(A,b,eps)
%myls_svd 此处提供此函数的摘要
%   此处提供详细说明
[~,n] = size(A);
[U,S,V] = svd(A,0);
if nargin<3
    eps = 1e-12;
end
s = diag(S);
r = sum(s>eps);
x = zeros(n,1);
for i=1:r
    x = x + (U(:,i)'*b)/s(i)*V(:,i);
end
if (nargout>1)
    res = norm(A*x - b)/norm(b);
end
end
