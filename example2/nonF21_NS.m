function [time]=nonF21_NS(Lx,Ly,Lz,T,M,N,alpha,delta,nu)
% function [error_u1,error_u2,error_p,Mass]=nonF21_NS(Lx,Ly,Lz,T,M,N,alpha,nu)
% function nonF21_sigma_NS_2D
% clear;clc; M=100，200，400，800
%%%%%for 3D，有真解三维方程
%% 给定初值
% Ly=2*pi; Lx=2*pi; Lz=2*pi;
% N=8; 
% nu=100; alpha=0.5; 
% T=1; M=640;

up=Ly; bottom=0; 
right=Lx; left=0;
height=Lz; down=0;

Nx=N; Ny=N; Nz=N;
h1=abs(up-bottom)/Ny; 
h2=abs(right-left)/Nx;
h3=abs(height-down)/Nz;

sigma=1-alpha/2;
% delta=2/alpha;

Ymesh=bottom:h1:up-h1;
Xmesh=left:h2:right-h2;
Zmesh=down:h3:height-h3;

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

%% 相关矩阵计算
Qx=1i*2*pi/Lx*([0:1:Nx/2-1,-Nx/2:1:-1]); 
Qy=1i*2*pi/Ly*([0:1:Ny/2-1,-Ny/2:1:-1]);
Qz=1i*2*pi/Lz*([0:1:Nx/2-1,-Nx/2:1:-1]); 
Kx=kron(ones(Ny,1),Qx); Kx=repmat(Kx,1,1,Nz);
Ky=kron(ones(1,Nx),Qy(:)); Ky=repmat(Ky,1,1,Nz);
Kz=reshape(Qz,1,1,Nz); Kz=repmat(Kz,Ny,Nx,1);
Kxxyy=Kx.^2+Ky.^2+Kz.^2;
Khat=1./Kxxyy;
Khat(1,1,1)=0;

%% 计算真解
U1_exact=fun_u1(Xmesh,Ymesh,Zmesh,Tmesh,alpha);
U2_exact=fun_u2(Xmesh,Ymesh,Zmesh,Tmesh,alpha);
U3_exact=fun_u3(Xmesh,Ymesh,Zmesh,Tmesh,alpha);
P_exact=fun_p(Xmesh,Ymesh,Zmesh,Tmesh,alpha);

%% 给定初值
U1=zeros(Nx,Ny,Nz,M+2); U2=zeros(Nx,Ny,Nz,M+2); U3=zeros(Nx,Ny,Nz,M+2);
P=zeros(Nx,Ny,Nz,M+2);

U1(:,:,:,1)=U1_exact(:,:,:,1); U10=U1(:,:,:,1);
U2(:,:,:,1)=U2_exact(:,:,:,1); U20=U2(:,:,:,1);
U3(:,:,:,1)=U3_exact(:,:,:,1); U30=U3(:,:,:,1);
P(:,:,:,1)=P_exact(:,:,:,1);

temp=Kx.*fftn(U10)+Ky.*fftn(U20)+Kz.*fftn(U30); 
temp=real(ifftn(temp));
Mass=sum(sum(sum(temp)));
    
%% 计算第一步
t=sigma*Tmesh(2)+(1-sigma)*Tmesh(1);
An=zeros(2-1,1);
An(1)=1/(tau(2)*gamma(2-alpha))*power(t,1-alpha);
a0=An(1);
Kmatrix=a0*ones(Ny,Nx,Nz)-nu*sigma*Kxxyy;
error=1; number=1; ep=1e-8;

while error>ep

    %中间点的非线性项
    U1_sigma=sigma*U10+(1-sigma)*U1(:,:,:,1);
    U2_sigma=sigma*U20+(1-sigma)*U2(:,:,:,1);
    U3_sigma=sigma*U30+(1-sigma)*U3(:,:,:,1);
    
    F1=U1_sigma.*real(ifftn(Kx.*fftn(U1_sigma)))...
        +U2_sigma.*real(ifftn(Ky.*fftn(U1_sigma)))...
        +U3_sigma.*real(ifftn(Kz.*fftn(U1_sigma))); F1(1,1,1)=0;
    F2=U1_sigma.*real(ifftn(Kx.*fftn(U2_sigma)))...
        +U2_sigma.*real(ifftn(Ky.*fftn(U2_sigma)))...
        +U3_sigma.*real(ifftn(Kz.*fftn(U2_sigma))); F2(1,1,1)=0;
    F3=U1_sigma.*real(ifftn(Kx.*fftn(U3_sigma)))...
        +U2_sigma.*real(ifftn(Ky.*fftn(U3_sigma)))...
        +U3_sigma.*real(ifftn(Kz.*fftn(U3_sigma))); F3(1,1,1)=0;
    
    force_term1=f1(sigma*Tmesh(2)+(1-sigma)*Tmesh(1),Xmesh,Ymesh,Zmesh,alpha,nu); 
    ff1=fftn(force_term1);
    force_term2=f2(sigma*Tmesh(2)+(1-sigma)*Tmesh(1),Xmesh,Ymesh,Zmesh,alpha,nu); 
    ff2=fftn(force_term2);
    force_term3=f3(sigma*Tmesh(2)+(1-sigma)*Tmesh(1),Xmesh,Ymesh,Zmesh,alpha,nu); 
    ff3=fftn(force_term3);
    
    P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3)-Kx.*ff1-Ky.*ff2-Kz.*ff3);
    
    %计算新的迭代值
    U1_new=((a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U1(:,:,:,1))...
        -fftn(F1)-Kx.*P_f+ff1)./Kmatrix;
    U2_new=((a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U2(:,:,:,1))...
        -fftn(F2)-Ky.*P_f+ff2)./Kmatrix;
    U3_new=((a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U3(:,:,:,1))...
        -fftn(F3)-Kz.*P_f+ff3)./Kmatrix;
    
    U1_new=real(ifftn(U1_new)); 
    U2_new=real(ifftn(U2_new));
    U3_new=real(ifftn(U3_new)); 

    error=max(max(max(abs([U10-U1_new;U20-U2_new;U30-U3_new]))));
    U10=U1_new; U20=U2_new; U30=U3_new;

    number=number+1;
    if number>100
        disp(2);
        break;
    end
end

U1(:,:,:,2)=U10;   U2(:,:,:,2)=U20;  U3(:,:,:,2)=U30;
F1=U10.*real(ifftn(Kx.*fftn(U10)))+U20.*real(ifftn(Ky.*fftn(U10)))+U30.*real(ifftn(Kz.*fftn(U10)));
F2=U10.*real(ifftn(Kx.*fftn(U20)))+U20.*real(ifftn(Ky.*fftn(U20)))+U30.*real(ifftn(Kz.*fftn(U20)));
F3=U10.*real(ifftn(Kx.*fftn(U30)))+U20.*real(ifftn(Ky.*fftn(U30)))+U30.*real(ifftn(Kz.*fftn(U30)));
P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3)-Kx.*ff1-Ky.*ff2-Kz.*ff3);

P(:,:,:,2)=real(ifftn(P_f));

temp=Kx.*fftn(U10)+Ky.*fftn(U20)+Kz.*fftn(U30); temp=real(ifftn(temp));
Mass=[Mass;sum(sum(sum(temp)))];

%% 循环
ep=1e-8; 

U1_before=zeros(Ny,Nx,Nz,N_exp);  
U2_before=zeros(Ny,Nx,Nz,N_exp); 
U3_before=zeros(Ny,Nx,Nz,N_exp);

t1=clock;
for i=3:M+1
    
    U10=U1(:,:,:,i-1); U20=U2(:,:,:,i-1); U30=U3(:,:,:,i-1); 
    error=1; number=1;
    t=sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1);
    temp1=zeros(Ny,Nx,Nz); temp2=zeros(Ny,Nx,Nz); temp3=zeros(Ny,Nx,Nz);
    tau_sigma=sigma*tau(i)+(1-sigma)*tau(i-1);
    Bn=zeros(N_exp,3);
    
    for k=1:N_exp
        
        Bn(k,1)=(coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i-1),S(k))+coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i-2),S(k)))/(tau(i)*(tau(i)+tau(i-1)));
        Bn(k,2)=-(coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i),S(k))+coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i-2),S(k)))/(tau(i)*tau(i-1));
        Bn(k,3)=(coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i),S(k))+coe_fun(Tmesh(i-1),Tmesh(i-2),t,Tmesh(i-1),S(k)))/(tau(i-1)*(tau(i)+tau(i-1)));

        temp1=temp1+omega(k)*(exp(-S(k)*tau_sigma)*U1_before(:,:,:,k)+Bn(k,2)*U1(:,:,:,i-1)+Bn(k,3)*U1(:,:,:,i-2));
        temp2=temp2+omega(k)*(exp(-S(k)*tau_sigma)*U2_before(:,:,:,k)+Bn(k,2)*U2(:,:,:,i-1)+Bn(k,3)*U2(:,:,:,i-2));
        temp3=temp3+omega(k)*(exp(-S(k)*tau_sigma)*U3_before(:,:,:,k)+Bn(k,2)*U3(:,:,:,i-1)+Bn(k,3)*U3(:,:,:,i-2));
    end
    
    mid=omega.*Bn(:,1);
    a0=power(t-Tmesh(i-1),1-alpha)/(tau(i)*gamma(2-alpha));
    Kmatrix=(a0+sum(mid)/gamma(1-alpha))*ones(Ny,Nx,Nz)-nu*sigma*Kxxyy;
       
    while error>ep  

        %中间点的非线性项
        U1_sigma=sigma*U10+(1-sigma)*U1(:,:,:,i-1);
        U2_sigma=sigma*U20+(1-sigma)*U2(:,:,:,i-1);
        U3_sigma=sigma*U30+(1-sigma)*U3(:,:,:,i-1);
        
        F1=U1_sigma.*real(ifftn(Kx.*fftn(U1_sigma)))...
            +U2_sigma.*real(ifftn(Ky.*fftn(U1_sigma)))...
            +U3_sigma.*real(ifftn(Kz.*fftn(U1_sigma))); F1(1,1,1)=0;
        F2=U1_sigma.*real(ifftn(Kx.*fftn(U2_sigma)))...
            +U2_sigma.*real(ifftn(Ky.*fftn(U2_sigma)))...
            +U3_sigma.*real(ifftn(Kz.*fftn(U2_sigma))); F2(1,1,1)=0;
        F3=U1_sigma.*real(ifftn(Kx.*fftn(U3_sigma)))...
            +U2_sigma.*real(ifftn(Ky.*fftn(U3_sigma)))...
            +U3_sigma.*real(ifftn(Kz.*fftn(U3_sigma))); F3(1,1,1)=0;
        
        force_term1=f1(sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1),Xmesh,Ymesh,Zmesh,alpha,nu); 
        ff1=fftn(force_term1);
        force_term2=f2(sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1),Xmesh,Ymesh,Zmesh,alpha,nu); 
        ff2=fftn(force_term2);
        force_term3=f3(sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1),Xmesh,Ymesh,Zmesh,alpha,nu); 
        ff3=fftn(force_term3);

        P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3)-Kx.*ff1-Ky.*ff2-Kz.*ff3);
        
        %计算新的迭代值
        U1_new=(-fftn(temp1/gamma(1-alpha))+(a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U1(:,:,:,i-1))...
            -fftn(F1)-Kx.*P_f+ff1)./Kmatrix;
        U2_new=(-fftn(temp2/gamma(1-alpha))+(a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U2(:,:,:,i-1))...
            -fftn(F2)-Ky.*P_f+ff2)./Kmatrix;
        U3_new=(-fftn(temp3/gamma(1-alpha))+(a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U3(:,:,:,i-1))...
            -fftn(F3)-Kz.*P_f+ff3)./Kmatrix;
    
        U1_new=real(ifftn(U1_new)); 
        U2_new=real(ifftn(U2_new));
        U3_new=real(ifftn(U3_new));
        
        error=max(max(max(abs([U10-U1_new;U20-U2_new;U30-U3_new]))));
        U10=U1_new; U20=U2_new; U30=U3_new;

        number=number+1;
        if number>100
            disp(i);
            break;
        end
    end
        
    U1(:,:,:,i)=U10; U2(:,:,:,i)=U20; U3(:,:,:,i)=U30;
    F1=U10.*real(ifftn(Kx.*fftn(U10)))+U20.*real(ifftn(Ky.*fftn(U10)))+U30.*real(ifftn(Kz.*fftn(U10)));
    F2=U10.*real(ifftn(Kx.*fftn(U20)))+U20.*real(ifftn(Ky.*fftn(U20)))+U30.*real(ifftn(Kz.*fftn(U20)));    
    F3=U10.*real(ifftn(Kx.*fftn(U30)))+U20.*real(ifftn(Ky.*fftn(U30)))+U30.*real(ifftn(Kz.*fftn(U30))); 
    
    force_term1=f1(Tmesh(i),Xmesh,Ymesh,Zmesh,alpha,nu); 
    ff1=fftn(force_term1);
    force_term2=f2(Tmesh(i),Xmesh,Ymesh,Zmesh,alpha,nu); 
    ff2=fftn(force_term2);
    force_term3=f3(Tmesh(i),Xmesh,Ymesh,Zmesh,alpha,nu); 
    ff3=fftn(force_term3);
    
    P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3)-Kx.*ff1-Ky.*ff2-Kz.*ff3);
    
    P(:,:,:,i)=real(ifftn(P_f));
    
    temp=Kx.*fftn(U10)+Ky.*fftn(U20)+Kz.*fftn(U30);  
    temp=real(ifftn(temp));
    
    Mass=[Mass;sum(sum(sum(temp)))];
    
    for k=1:N_exp
        U1_before(:,:,:,k)=exp(-S(k)*tau_sigma)*U1_before(:,:,:,k)+Bn(k,1)*U1(:,:,:,i)+Bn(k,2)*U1(:,:,:,i-1)+Bn(k,3)*U1(:,:,:,i-2);
        U2_before(:,:,:,k)=exp(-S(k)*tau_sigma)*U2_before(:,:,:,k)+Bn(k,1)*U2(:,:,:,i)+Bn(k,2)*U2(:,:,:,i-1)+Bn(k,3)*U2(:,:,:,i-2);
        U3_before(:,:,:,k)=exp(-S(k)*tau_sigma)*U3_before(:,:,:,k)+Bn(k,1)*U3(:,:,:,i)+Bn(k,2)*U3(:,:,:,i-1)+Bn(k,3)*U3(:,:,:,i-2);
    end
    
end

t2=clock;
time=etime(t2,t1);
U1(:,:,:,end)=[]; U2(:,:,:,end)=[]; U3(:,:,:,end)=[]; P(:,:,:,end)=[];

end
