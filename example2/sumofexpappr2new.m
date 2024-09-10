function [xs,ws,nexp] = sumofexpappr2new(alpha,reps,dt,Tfinal)
%sumofexpappr2new Caputo分数阶导数的核函数指数和逼近
% Caputo 分数阶导数核函数t^{-alpha}的指数和逼近的 MATLAB 程序代码
% Caputo 分数阶导数的核函数t^{-alpha}可以用指数和逼近.
% 应用这个逼近公式，可以对 Caputo 分数阶导数给出快速的 L1 逼近等.
% 测试: [xs,ws,nexp] =sumofexpappr2new(1.2,10^(-3),0.01,1);

delta = dt/Tfinal;
h  = 2*pi/(log(3) + alpha*log(1/cos(1)) + log(1/reps));
tlower = 1/alpha*log(reps*gamma(1 + alpha));
if alpha>=1
    tupper = log(1/delta) + log(log(1/reps)) + log(alpha) + 1/2;
else
    tupper = log(1/delta) + log(log(1/reps));
end
M = floor(tlower/h);
N = ceil(tupper/h);
n1 = M:-1;
xs1 = -exp(h*n1);
ws1 = h/gamma(alpha)*exp(alpha*h*n1);
[ws1new,xs1new] = prony(xs1,ws1);
n2 = 0:N;
xs2 = -exp(h*n2);
ws2 = h/gamma(alpha)*exp(alpha*h*n2);
xs = [-real(xs1new);-real(xs2.')];
ws = [real(ws1new);real(ws2.')];
xs = xs/Tfinal;
ws = ws/Tfinal^alpha;
nexp = length(ws);
return
end
