clear;clc;

%convergence rate calculation
Ly=2*pi; Lx=2*pi; N=8; 
nu=1; 
alpha=[0.4,0.6,0.8]  ; 
T=1;

M=[4000;8000;16000;32000;64000];
T1=zeros(length(M),length(alpha)); 
T2=T1;

for m=1:length(alpha)

    delta=2;
    
    for k=1:length(M)    
        [time_L21]=nonF21_sigma_NS_2D(Lx,Ly,T,M(k),N,alpha(m),delta,nu);
        [time_fast]=nonF21_NS(Lx,Ly,T,M(k),N,delta,alpha(m),nu);
        T1(k,m)=time_L21; 
        T2(k,m)=time_fast;
    end
    
end


disp(T1)
disp(T2)
save T1