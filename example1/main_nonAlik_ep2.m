clear;clc;

%convergence rate calculation
Ly=2*pi; Lx=2*pi; N=8; 
nu=1; alpha=[0.4,0.6,0.8]; 
T=1;

M=[8;16;32;64];
EU1=zeros(length(M),length(alpha));  EU2=EU1; EP=EU1;

for m=1:length(alpha)

    delta=2/alpha(m);
    
    for k=1:length(M)
    
        [error_u1,error_u2,error_p,Mass]=nonF21_sigma_NS_2D(Lx,Ly,T,M(k),N,alpha(m),delta,nu);
        EU1(k,m)=error_u1(end); EU2(k,m)=error_u2(end); EP(k,m)=error_p(end);
        
    end
    
end

for m=1:length(alpha)
    order_u1(:,m)=log(EU1(1:end-1,m)./EU1(2:end,m))./(log(M(2:end)./M(1:end-1)));
    order_u2(:,m)=log(EU2(1:end-1,m)./EU2(2:end,m))./(log(M(2:end)./M(1:end-1)));
    order_p(:,m)=log(EP(1:end-1,m)./EP(2:end,m))./(log(M(2:end)./M(1:end-1)));
end

format shortE
disp(EU1)
disp(order_u1)
disp(EU2)
disp(order_u2)
disp(EP)
disp(order_p)


