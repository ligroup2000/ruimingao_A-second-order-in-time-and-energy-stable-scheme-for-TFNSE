clear;clc;

%convergence rate calculation
Ly=2*pi; Lx=2*pi; Lz=2*pi;
N=8; 
nu=1; alpha=[0.4;0.6;0.8]; 
T=1;

M=[8;16;32;64];
EU=zeros(length(M),length(alpha));  EP=EU;

for m=1:length(alpha) 

    for k=1:length(M)
    
        [error_u1,error_u2,error_u3,error_p,Mass]=nonF21_sigma_NS_2D(Lx,Ly,Lz,T,M(k),N,alpha(m),nu);
        EU(k,m)=max(max([error_u1,error_u2,error_u3])); 
        EP(k,m)=max(error_p);
        
    end
    
end

for m=1:length(alpha)
    order_u(:,m)=log(EU(1:end-1,m)./EU(2:end,m))./(log(M(2:end)./M(1:end-1)));
    order_p(:,m)=log(EP(1:end-1,m)./EP(2:end,m))./(log(M(2:end)./M(1:end-1)));
end

format shortE
disp(EU)
disp(order_u)
disp(EP)
disp(order_p)
