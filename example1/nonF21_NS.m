function [time]=nonF21_NS(Lx,Ly,T,M,N,delta,alpha,nu)
% function [error_u1,error_u2,error_p,Mass]=nonF21_NS(Lx,Ly,T,M,N,delta,alpha,nu)
% function nonF21_NS
% clear;clc;
% 快速Aliknov格式计算Navier-Stokes方程
%算例2：奇异性真解
%% 给定初值
% Ly=2*pi; Lx=2*pi; N=8; 
% nu=1; alpha=0.4; 
% T=1;  M=20;
% delta=2;

sigma=1-alpha/2;
up=Ly; bottom=0; right=Lx; left=0;
Nx=N; Ny=N;
h1=abs(up-bottom)/Ny; h2=abs(right-left)/Nx;

Ymesh=bottom:h1:up-h1;
Xmesh=left:h2:right-h2;

coe_fun=@(s1,s2,t,eta,weight) (s1-eta)/weight*exp(-weight*(t-s1))-exp(-weight*(t-s1))/(weight^2)...
-(s2-eta)/weight*exp(-weight*(t-s2))+exp(-weight*(t-s2))/(weight^2);  
%% 时间网格划分
Tmesh=zeros(1,M+1);
for k=1:M+1
    Tmesh(k)=T*power((k-1)/M,delta);
end
tau=[0,Tmesh(2:end)-Tmesh(1:end-1)];
tau_min=min(tau(2:end));

if tau_min<1e-14
    tau_min=1e-14;
end
ep=1e-8;
[S,omega,N_exp] = sumofexpappr2new(alpha,ep,tau_min,T);
% [omega, zez] = updateweights(S,tau_min);

%% 相关矩阵计算
Qx=1i*2*pi/Lx*([0:1:Nx/2-1,-Nx/2:1:-1]); 
Qy=1i*2*pi/Ly*([0:1:Ny/2-1,-Ny/2:1:-1]);
Kx=kron(ones(Ny,1),Qx); Ky=kron(ones(1,Nx),Qy(:));
Kxxyy=Kx.^2+Ky.^2;
Khat=1./Kxxyy;
Khat(1,1)=0;

%% 计算真解
U1_exact=fun_u1(Xmesh,Ymesh,Tmesh,alpha);
U2_exact=fun_u2(Xmesh,Ymesh,Tmesh,alpha);
P_exact=fun_p(Xmesh,Ymesh,Tmesh,alpha);

%% 给定初值
U1=zeros(Nx,Ny,M+2); U2=zeros(Nx,Ny,M+2); P=zeros(Nx,Ny,M+2);
U1(:,:,1)=U1_exact(:,:,1); U10=U1(:,:,1);
U2(:,:,1)=U2_exact(:,:,1); U20=U2(:,:,1);
P(:,:,1)=P_exact(:,:,1);

temp=Kx.*fft2(U10)+Ky.*fft2(U20); temp=real(ifft2(temp));
Mass=sum(sum(temp));
    
%% 计算第一步
t=sigma*Tmesh(2)+(1-sigma)*Tmesh(1);
a0=power(t,1-alpha)/(tau(2)*gamma(2-alpha));
Kmatrix=a0*ones(Ny,Nx)-nu*sigma*Kxxyy;
U10=U1(:,:,1); U20=U2(:,:,1); 
error=1; number=1;

while error>ep

    %中间点的非线性项
    U1_sigma=sigma*U10+(1-sigma)*U1(:,:,1);
    U2_sigma=sigma*U20+(1-sigma)*U2(:,:,1);
    F1=U1_sigma.*real(ifft2(Kx.*fft2(U1_sigma)))...
        +U2_sigma.*real(ifft2(Ky.*fft2(U1_sigma))); F1(1,1)=0;
    F2=U1_sigma.*real(ifft2(Kx.*fft2(U2_sigma)))...
        +U2_sigma.*real(ifft2(Ky.*fft2(U2_sigma))); F2(1,1)=0;
    
    force_term1=f1(sigma*Tmesh(2)+(1-sigma)*Tmesh(1),Xmesh,Ymesh,alpha,nu); 
    ff1=fft2(force_term1);
    force_term2=f2(sigma*Tmesh(2)+(1-sigma)*Tmesh(1),Xmesh,Ymesh,alpha,nu); 
    ff2=fft2(force_term2);
    
    P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2)-Kx.*ff1-Ky.*ff2);
    
    %计算新的迭代值
    U1_new=((a0*ones(Ny,Nx)+nu*(1-sigma)*Kxxyy).*fft2(U1(:,:,1))...
        -fft2(F1)-Kx.*P_f+ff1)./Kmatrix;
    U2_new=((a0*ones(Ny,Nx)+nu*(1-sigma)*Kxxyy).*fft2(U2(:,:,1))...
        -fft2(F2)-Ky.*P_f+ff2)./Kmatrix;
    
    U1_new=real(ifft2(U1_new)); U2_new=real(ifft2(U2_new));

    error=max(max(abs([U10-U1_new;U20-U2_new])));
    U10=U1_new; U20=U2_new;

    number=number+1;
    if number>100
        disp(2);
        break;
    end
end

U1(:,:,2)=U10;   U2(:,:,2)=U20;
F1=U10.*real(ifft2(Kx.*fft2(U10)))+U20.*real(ifft2(Ky.*fft2(U10)));
F2=U10.*real(ifft2(Kx.*fft2(U20)))+U20.*real(ifft2(Ky.*fft2(U20)));
P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2)-Kx.*ff1-Ky.*ff2);

P(:,:,2)=real(ifft2(P_f));

temp=Kx.*fft2(U10)+Ky.*fft2(U20); temp=real(ifft2(temp));
Mass=[Mass;sum(sum(temp))];
%% 循环
ep=1e-10;  
U1_before=zeros(Ny,Nx,N_exp);  
U2_before=zeros(Ny,Nx,N_exp);  

t1=clock;
for i=3:M+1
    
    U10=U1(:,:,i-1); U20=U2(:,:,i-1); 
    error=1; number=1;
    temp1=zeros(Nx,Ny); temp2=temp1;
    t=sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1);
    tau_sigma=sigma*tau(i)+(1-sigma)*tau(i-1);  
    Bn=zeros(N_exp,3);
    
    for k=1:N_exp
        
        Bn(k,1)=(coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i-1),S(k))+coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i-2),S(k)))/(tau(i)*(tau(i)+tau(i-1)));
        Bn(k,2)=-(coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i),S(k))+coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i-2),S(k)))/(tau(i)*tau(i-1));
        Bn(k,3)=(coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i),S(k))+coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i-1),S(k)))/(tau(i-1)*(tau(i)+tau(i-1)));

        temp1=temp1+omega(k)*(exp(-S(k)*tau_sigma)*U1_before(:,:,k)+Bn(k,2)*U1(:,:,i-1)+Bn(k,3)*U1(:,:,i-2));
        temp2=temp2+omega(k)*(exp(-S(k)*tau_sigma)*U2_before(:,:,k)+Bn(k,2)*U2(:,:,i-1)+Bn(k,3)*U2(:,:,i-2));
    end
    
    mid=omega.*Bn(:,1);
    a0=power(t-Tmesh(i-1),1-alpha)/(tau(i)*gamma(2-alpha));
    Kmatrix=(a0+sum(mid)/gamma(1-alpha))*ones(Ny,Nx)-nu*sigma*Kxxyy;
    
    while error>ep
        
         %中间点的非线性项
        U1_sigma=sigma*U10+(1-sigma)*U1(:,:,i-1);
        U2_sigma=sigma*U20+(1-sigma)*U2(:,:,i-1);
        F1=U1_sigma.*real(ifft2(Kx.*fft2(U1_sigma)))...
            +U2_sigma.*real(ifft2(Ky.*fft2(U1_sigma))); F1(1,1)=0;
        F2=U1_sigma.*real(ifft2(Kx.*fft2(U2_sigma)))...
            +U2_sigma.*real(ifft2(Ky.*fft2(U2_sigma))); F2(1,1)=0;
        
        force_term1=f1(sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1),Xmesh,Ymesh,alpha,nu); 
        ff1=fft2(force_term1);
        force_term2=f2(sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1),Xmesh,Ymesh,alpha,nu); 
        ff2=fft2(force_term2);

        P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2)-Kx.*ff1-Ky.*ff2);
        
        %计算新的迭代值
        U1_new=(-fft2(temp1/gamma(1-alpha))+(a0*ones(Ny,Nx)+nu*(1-sigma)*Kxxyy).*fft2(U1(:,:,i-1))...
            -fft2(F1)-Kx.*P_f+ff1)./Kmatrix;
        U2_new=(-fft2(temp2/gamma(1-alpha))+(a0*ones(Ny,Nx)+nu*(1-sigma)*Kxxyy).*fft2(U2(:,:,i-1))...
            -fft2(F2)-Ky.*P_f+ff2)./Kmatrix;
    
        U1_new=real(ifft2(U1_new)); U2_new=real(ifft2(U2_new));
        
        error=max(max(abs([U10-U1_new;U20-U2_new])));
        U10=U1_new; U20=U2_new;

        number=number+1;
        if number>100
            disp(i);
            break;
        end
    end
        
    U1(:,:,i)=U10; U2(:,:,i)=U20;
    F1=U10.*real(ifft2(Kx.*fft2(U10)))+U20.*real(ifft2(Ky.*fft2(U10)));
    F2=U10.*real(ifft2(Kx.*fft2(U20)))+U20.*real(ifft2(Ky.*fft2(U20)));
    force_term1=f1(Tmesh(i),Xmesh,Ymesh,alpha,nu); 
    ff1=fft2(force_term1);
    force_term2=f2(Tmesh(i),Xmesh,Ymesh,alpha,nu); 
    ff2=fft2(force_term2);
    P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2)-Kx.*ff1-Ky.*ff2);
    
    P(:,:,i)=real(ifft2(P_f));
    
    temp=Kx.*fft2(U10)+Ky.*fft2(U20);  temp=real(ifft2(temp));
    
    Mass=[Mass;sum(sum(temp))];
    
     for k=1:N_exp
        U1_before(:,:,k)=exp(-S(k)*tau_sigma)*U1_before(:,:,k)+Bn(k,1)*U1(:,:,i)+Bn(k,2)*U1(:,:,i-1)+Bn(k,3)*U1(:,:,i-2);
        U2_before(:,:,k)=exp(-S(k)*tau_sigma)*U2_before(:,:,k)+Bn(k,1)*U2(:,:,i)+Bn(k,2)*U2(:,:,i-1)+Bn(k,3)*U2(:,:,i-2);
    end   
    
end

t2=clock;
time=etime(t2,t1);
U1(:,:,end)=[]; U2(:,:,end)=[]; P(:,:,end)=[];


% %% 计算误差
% error1=abs(U1-U1_exact);
% error2=abs(U2-U2_exact);
% errorp=abs(P-P_exact);
% error_u1=zeros(length(Tmesh),1); error_u2=error_u1; error_p=error_u1;
% 
% for j=1:length(Tmesh)
%     error_u1(j)=sum(sum(error1(:,:,j).^2))+...
%         sum(sum(real(ifft2(Kx.*fft2(error1(:,:,j)))).^2))+...
%         sum(sum(real(ifft2(Ky.*fft2(error1(:,:,j)))).^2)); error_u1(j)=sqrt(h1*h2*error_u1(j));
%    error_u2(j)=sum(sum(error2(:,:,j).^2))+...
%         sum(sum(real(ifft2(Kx.*fft2(error2(:,:,j)))).^2))+...
%         sum(sum(real(ifft2(Ky.*fft2(error2(:,:,j)))).^2)); error_u2(j)=sqrt(h1*h2*error_u2(j));
%     error_p(j)=sum(sum(sum(errorp(:,:,j).^2))); error_p(j)=sqrt(h1*h2*error_p(j));
% end

% disp(max(error_u1)); disp(max(error_u2)); disp(max(error_p))
end
