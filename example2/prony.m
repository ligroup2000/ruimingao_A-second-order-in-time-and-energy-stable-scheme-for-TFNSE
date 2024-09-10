function [wsnew,xsnew] = prony(xs,ws)
%prony 此处提供此函数的摘要
%   此处提供详细说明
M = length(xs);
errbnd = 1d-12;
h = zeros(2*M,1);
for j=1:2*M
    h(j) = xs.^(j-1)*ws';
end
C = h(1:M);
R = h(M:2*M-1);
H = hankel(C,R);
b = -h;
q = myls_qr(H,b,errbnd);
r = length(q);
A = zeros(2*M,r);
Coef = [1;flipud(q)];
xsnew = roots(Coef);
for j=1:2*M
    A(j,:) = xsnew.^(j-1);
end
wsnew = myls_svd(A,h,errbnd);
ind =find(real(xsnew)>=0);
p = length(ind);
assert(sum(abs(wsnew(ind))<1d-15)==p)
ind = find(real(xsnew)<0);
xsnew = xsnew(ind);
wsnew = wsnew(ind);
end
