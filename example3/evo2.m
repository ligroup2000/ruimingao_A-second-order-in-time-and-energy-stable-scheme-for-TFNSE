% function [Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,N,nu,eta,alpha,delta)
% function fracNS_Cal_energy
clear;clc;

%% 给定初值
Ly=1; Lx=1; N=128; 
nu=1e-4; alpha=0.8; 
T=1.2; 
eta=10;

up=Ly; bottom=0; right=Lx; left=0;
Nx=N; Ny=N;
h1=abs(up-bottom)/Ny; h2=abs(right-left)/Nx;
sigma=1-alpha/2;
delta=2/alpha;

Ymesh=bottom:h1:up-h1;
Xmesh=left:h2:right-h2;

%% 时间网格划分
% T=1;
T0=0.01; M=20; 
Tmesh=zeros(M+1,1);
for k=1:M+1
    Tmesh(k)=T0*power((k-1)/M,delta);
end
tau=[0;Tmesh(2:end)-Tmesh(1:end-1)];

%% 相关矩阵计算
Qx=1i*2*pi/Lx*([0:1:Nx/2-1,-Nx/2:1:-1]); 
Qy=1i*2*pi/Ly*([0:1:Ny/2-1,-Ny/2:1:-1]);
Kx=kron(ones(Ny,1),Qx); Ky=kron(ones(1,Nx),Qy(:));
Kxxyy=Kx.^2+Ky.^2;
Khat=1./Kxxyy;
Khat(1,1)=0;

%% 给定初值
U1=zeros(Nx,Ny,M+2); U2=zeros(Nx,Ny,M+2); P=zeros(Nx,Ny,M+2);
Evolu=zeros(Nx,Ny,M+2);
U1(:,:,1)=initial_u1(Xmesh,Ymesh,Tmesh,alpha,100);
U2(:,:,1)=initial_u2(Xmesh,Ymesh,Tmesh,alpha);
U10=U1(:,:,1); U20=U2(:,:,1);
F1=U10.*real(ifft2(Kx.*fft2(U10)))+U20.*real(ifft2(Ky.*fft2(U10)));
F2=U10.*real(ifft2(Kx.*fft2(U20)))+U20.*real(ifft2(Ky.*fft2(U20)));

P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2)); 
P(:,:,1)=real(ifft2(P_f));

temp=Kx.*fft2(U10)+Ky.*fft2(U20); temp=real(ifft2(temp));
Mass=sum(sum(temp));
Evolu(:,:,1)=Kx.*fft2(U20)-Ky.*fft2(U10);
bas_energy=(Ly*Lx/Nx/Ny)*(sum(sum(U1(:,:,1).^2))+sum(sum(U2(:,:,1).^2))); 
Energy=bas_energy/2;
    
%% 计算第一步
t=sigma*Tmesh(2)+(1-sigma)*Tmesh(1);
An=zeros(2-1,1);
Star=zeros(2-1,1);
An(1)=1/(tau(2)*gamma(2-alpha))*power(t,1-alpha);
a0=An(1);
Star(1,:)=An;
Kmatrix=a0*ones(Ny,Nx)-nu*sigma*Kxxyy;
U10=U1(:,:,1); U20=U2(:,:,1); 
error=1; number=1; ep=1e-8;

while error>ep

    %中间点的非线性项
    U1_sigma=sigma*U10+(1-sigma)*U1(:,:,1);
    U2_sigma=sigma*U20+(1-sigma)*U2(:,:,1);
    F1=U1_sigma.*real(ifft2(Kx.*fft2(U1_sigma)))...
        +U2_sigma.*real(ifft2(Ky.*fft2(U1_sigma))); F1(1,1)=0;
    F2=U1_sigma.*real(ifft2(Kx.*fft2(U2_sigma)))...
        +U2_sigma.*real(ifft2(Ky.*fft2(U2_sigma))); F2(1,1)=0;
    
    P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));
    
    %计算新的迭代值
    U1_new=((a0*ones(Ny,Nx)+nu*(1-sigma)*Kxxyy).*fft2(U1(:,:,1))...
        -fft2(F1)-Kx.*P_f)./Kmatrix;
    U2_new=((a0*ones(Ny,Nx)+nu*(1-sigma)*Kxxyy).*fft2(U2(:,:,1))...
        -fft2(F2)-Ky.*P_f)./Kmatrix;
    
    U1_new=real(ifft2(U1_new)); U2_new=real(ifft2(U2_new));

    error=max(max(abs([U10-U1_new;U20-U2_new])));
    U10=U1_new; U20=U2_new;

    number=number+1;
    if number>100
        disp(2);
        break;
    end
end

Ph=zeros(2-1,1);
Ph(1)=1/An(1);    

bas_energy=(sum(sum(U10.^2))+sum(sum(U20.^2)))/2;
for k=1:2-1
    bas_energy=bas_energy+nu*Ph(2-k)*...
        (sum(sum(real(ifft2(Kx.*fft2(sigma*U1(:,:,k+1)+(1-sigma)*U1(:,:,k)))).^2))+...
    sum(sum(real(ifft2(Ky.*fft2(sigma*U2(:,:,k+1)+(1-sigma)*U2(:,:,k)))).^2)));
end

Energy=[Energy;(Ly*Lx/Nx/Ny)*bas_energy];

U1(:,:,2)=U10;  U2(:,:,2)=U20;
Evolu(:,:,2)=Kx.*fft2(U20)-Ky.*fft2(U10);
F1=U10.*real(ifft2(Kx.*fft2(U10)))+U20.*real(ifft2(Ky.*fft2(U10)));
F2=U10.*real(ifft2(Kx.*fft2(U20)))+U20.*real(ifft2(Ky.*fft2(U20)));

P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));

P(:,:,2)=real(ifft2(P_f));

temp=Kx.*fft2(U10)+Ky.*fft2(U20);  temp=real(ifft2(temp));

Mass=[Mass;sum(sum(temp))];
%% 循环
ep=1e-12; 
for i=3:M+1
    
    U10=U1(:,:,i-1); U20=U2(:,:,i-1); 
    error=1; number=1;
    t=sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1);
    
    An=zeros(i-1,1); Bn=zeros(i-1,1);
    for k=1:i-1
        An(i-k)=1/(tau(k+1)*gamma(2-alpha))*(power(t-Tmesh(k),1-alpha)-power(t-min(t,Tmesh(k+1)),1-alpha));
        Bn(i-k)=2/(tau(k+1)^2*gamma(2-alpha))*((Tmesh(k)-(Tmesh(k)+Tmesh(k+1))/2)*power(t-Tmesh(k),1-alpha)-(Tmesh(k+1)-(Tmesh(k)+Tmesh(k+1))/2)*power(t-Tmesh(k+1),1-alpha))...
            +2/(tau(k+1)^2*gamma(3-alpha))*(power(t-Tmesh(k),2-alpha)-power(t-Tmesh(k+1),2-alpha));
    end
    Bn(1)=0;
    
    Cn=zeros(i-1,1);
    if i==3
        Cn(1)=An(1)+1/(tau(3)/tau(2)*(1+tau(3)/tau(2)))*Bn(2);
        Cn(2)=An(2)-1/(1+tau(3)/tau(2))*Bn(2);
    else 
        Cn(1)=An(1)+1/(tau(i)/tau(i-1)*(1+tau(i)/tau(i-1)))*Bn(2);
        for g=2:i-2
            Cn(g)=An(g)+1/(tau(i-g)/tau(i-g-1)*(1+tau(i-g)/tau(i-g-1)))*Bn(g+1)-1/(1+tau(i-g+1)/tau(i-g))*Bn(g);
        end
        Cn(i-1)=An(i-1)-1/(1+tau(3)/tau(2))*Bn(i-1);
    end
    
    a0=Cn(1);
    temp1=zeros(Ny,Nx); temp2=zeros(Ny,Nx);
    Star=[Star,zeros(i-2,1)];
    Star=[Star;Cn'];
    
    for s=2:i
        temp1=temp1+Cn(i-s+1)*(U1(:,:,s)-U1(:,:,s-1));
        temp2=temp2+Cn(i-s+1)*(U2(:,:,s)-U2(:,:,s-1));
    end
    
    Kmatrix=a0*ones(Ny,Nx)-nu*sigma*Kxxyy;
       
    while error>ep  

        %中间点的非线性项
        U1_sigma=sigma*U10+(1-sigma)*U1(:,:,i-1);
        U2_sigma=sigma*U20+(1-sigma)*U2(:,:,i-1);
        F1=U1_sigma.*real(ifft2(Kx.*fft2(U1_sigma)))...
            +U2_sigma.*real(ifft2(Ky.*fft2(U1_sigma))); F1(1,1)=0;
        F2=U1_sigma.*real(ifft2(Kx.*fft2(U2_sigma)))...
            +U2_sigma.*real(ifft2(Ky.*fft2(U2_sigma))); F2(1,1)=0;
        
        P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));
        
        %计算新的迭代值
        U1_new=(-fft2(temp1)+nu*(1-sigma)*Kxxyy.*fft2(U1(:,:,i-1))...
            -fft2(F1)-Kx.*P_f)./Kmatrix;
        U2_new=(-fft2(temp2)+nu*(1-sigma)*Kxxyy.*fft2(U2(:,:,i-1))...
            -fft2(F2)-Ky.*P_f)./Kmatrix;
    
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
    P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));
    
    P(:,:,i)=real(ifft2(P_f));
    
    temp=Kx.*fft2(U10)+Ky.*fft2(U20);  temp=real(ifft2(temp));
    
    Mass=[Mass;sum(sum(temp))];
    Evolu(:,:,i)=Kx.*fft2(U20)-Ky.*fft2(U10);
    
    %计算能量
    
    Ph=zeros(i-1,1); Ph(1)=1/Star(i-1,1);
    for j=2:i-1
        ss=0;
        for k=1:j-1
            ss=ss+(Star(i-k,j-k)-Star(i-k,j-k+1))*Ph(k);
        end
        Ph(j)=ss/Star(i-j,1);
    end
    
    bas_energy=(sum(sum(U10.^2))+sum(sum(U20.^2)))/2;
    for k=1:i-1
        bas_energy=bas_energy+...
            nu*Ph(i-k)*(sum(sum(real(ifft2(Kx.*fft2(sigma*U1(:,:,k+1)+(1-sigma)*U1(:,:,k)))).^2))+...
            sum(sum(real(ifft2(Ky.*fft2(sigma*U2(:,:,k+1)+(1-sigma)*U2(:,:,k)))).^2)));
    end

    Energy=[Energy;(Ly*Lx/Nx/Ny)*bas_energy];
    
end

% U1(:,:,end)=[]; U2(:,:,end)=[]; P(:,:,end)=[];
%% 变步长
tau_min=1e-4; tau_max=1e-2; i=M+1;
while (Tmesh(i)<T)
    
    i=i+1; U1(:,:,i)=zeros(Nx,Ny); U2(:,:,i)=zeros(Nx,Ny);
    U10=U1(:,:,i-1); U20=U2(:,:,i-1); 
    error=1; number=1;
    sub=abs(Energy(i-1)-Energy(i-2))/(tau(i-1)^2);
    tau=[tau;max(tau_min,tau_max/(sqrt(1+eta*sub)))];
    Tmesh=[Tmesh;Tmesh(end)+tau(end)];
    t=sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1);
    
    if Tmesh(end)+tau(end)>T
        break;
    end
    
    An=zeros(i-1,1); Bn=zeros(i-1,1);
    for k=1:i-1
        An(i-k)=1/(tau(k+1)*gamma(2-alpha))*(power(t-Tmesh(k),1-alpha)-power(t-min(t,Tmesh(k+1)),1-alpha));
        Bn(i-k)=2/(tau(k+1)^2*gamma(2-alpha))*((Tmesh(k)-(Tmesh(k)+Tmesh(k+1))/2)*power(t-Tmesh(k),1-alpha)-(Tmesh(k+1)-(Tmesh(k)+Tmesh(k+1))/2)*power(t-Tmesh(k+1),1-alpha))...
            +2/(tau(k+1)^2*gamma(3-alpha))*(power(t-Tmesh(k),2-alpha)-power(t-Tmesh(k+1),2-alpha));
    end
    Bn(1)=0; Bn=real(Bn);
    
    Cn=zeros(i-1,1);
    if i==3
        Cn(1)=An(1)+1/(tau(3)/tau(2)*(1+tau(3)/tau(2)))*Bn(2);
        Cn(2)=An(2)-1/(1+tau(3)/tau(2))*Bn(2);
    else 
        Cn(1)=An(1)+1/(tau(i)/tau(i-1)*(1+tau(i)/tau(i-1)))*Bn(2);
        for g=2:i-2
            Cn(g)=An(g)+1/(tau(i-g)/tau(i-g-1)*(1+tau(i-g)/tau(i-g-1)))*Bn(g+1)-1/(1+tau(i-g+1)/tau(i-g))*Bn(g);
        end
        Cn(i-1)=An(i-1)-1/(1+tau(3)/tau(2))*Bn(i-1);
    end
    
    a0=Cn(1);
    temp1=zeros(Ny,Nx); temp2=zeros(Ny,Nx);
    Star=[Star,zeros(i-2,1)];
    Star=[Star;Cn'];
    
    for s=2:i
        temp1=temp1+Cn(i-s+1)*(U1(:,:,s)-U1(:,:,s-1));
        temp2=temp2+Cn(i-s+1)*(U2(:,:,s)-U2(:,:,s-1));
    end
    
    Kmatrix=a0*ones(Ny,Nx)-nu*sigma*Kxxyy;
       
    while error>ep  

        %中间点的非线性项
        U1_sigma=sigma*U10+(1-sigma)*U1(:,:,i-1);
        U2_sigma=sigma*U20+(1-sigma)*U2(:,:,i-1);
        F1=U1_sigma.*real(ifft2(Kx.*fft2(U1_sigma)))...
            +U2_sigma.*real(ifft2(Ky.*fft2(U1_sigma))); F1(1,1)=0;
        F2=U1_sigma.*real(ifft2(Kx.*fft2(U2_sigma)))...
            +U2_sigma.*real(ifft2(Ky.*fft2(U2_sigma))); F2(1,1)=0;

        P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));
        
        %计算新的迭代值
        U1_new=(-fft2(temp1)+nu*(1-sigma)*Kxxyy.*fft2(U1(:,:,i-1))...
            -fft2(F1)-Kx.*P_f)./Kmatrix;
        U2_new=(-fft2(temp2)+nu*(1-sigma)*Kxxyy.*fft2(U2(:,:,i-1))...
            -fft2(F2)-Ky.*P_f)./Kmatrix;
    
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
    temp=Kx.*fft2(U10)+Ky.*fft2(U20);  temp=real(ifft2(temp));
    
    Mass=[Mass;sum(sum(temp))];
    Evolu(:,:,i)=Kx.*fft2(U20)-Ky.*fft2(U10);
    
    %计算能量
    Ph=zeros(i-1,1); Ph(1)=1/Star(i-1,1);
    for j=2:i-1
        ss=0;
        for k=1:j-1
            ss=ss+(Star(i-k,j-k)-Star(i-k,j-k+1))*Ph(k);
        end
        Ph(j)=ss/Star(i-j,1);
    end
    
    bas_energy=(sum(sum(U10.^2))+sum(sum(U20.^2)))/2;
    for k=1:i-1
        bas_energy=bas_energy+nu*Ph(i-k)*...
            (sum(sum(real(ifft2(Kx.*fft2(sigma*U1(:,:,k+1)+(1-sigma)*U1(:,:,k)))).^2))+...
            sum(sum(real(ifft2(Ky.*fft2(sigma*U2(:,:,k+1)+(1-sigma)*U2(:,:,k)))).^2)));
    end
    
    Energy=[Energy;(Ly*Lx/Nx/Ny)*bas_energy];
    
    Tmesh(i)=Tmesh(i-1)+tau(i);
    
end

%% 计算最后一步
i=length(Tmesh); tau(end)=T-Tmesh(i);  
Tmesh(i)=T;
t=sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1);

U10=U1(:,:,i-1); U20=U2(:,:,i-1); U1(:,:,i)=zeros(Nx,Ny); U2(:,:,i)=zeros(Nx,Ny);
error=1; number=1;

An=zeros(i-1,1); Bn=zeros(i-1,1);
for k=1:i-1
    An(i-k)=1/(tau(k+1)*gamma(2-alpha))*(power(t-Tmesh(k),1-alpha)-power(t-min(t,Tmesh(k+1)),1-alpha));
    Bn(i-k)=2/(tau(k+1)^2*gamma(2-alpha))*((Tmesh(k)-(Tmesh(k)+Tmesh(k+1))/2)*power(t-Tmesh(k),1-alpha)-(Tmesh(k+1)-(Tmesh(k)+Tmesh(k+1))/2)*power(t-Tmesh(k+1),1-alpha))...
        +2/(tau(k+1)^2*gamma(3-alpha))*(power(t-Tmesh(k),2-alpha)-power(t-Tmesh(k+1),2-alpha));
end
Bn(1)=0; Bn=real(Bn);

Cn=zeros(i-1,1);
if i==3
    Cn(1)=An(1)+1/(tau(3)/tau(2)*(1+tau(3)/tau(2)))*Bn(2);
    Cn(2)=An(2)-1/(1+tau(3)/tau(2))*Bn(2);
else 
    Cn(1)=An(1)+1/(tau(i)/tau(i-1)*(1+tau(i)/tau(i-1)))*Bn(2);
    for g=2:i-2
        Cn(g)=An(g)+1/(tau(i-g)/tau(i-g-1)*(1+tau(i-g)/tau(i-g-1)))*Bn(g+1)-1/(1+tau(i-g+1)/tau(i-g))*Bn(g);
    end
    Cn(i-1)=An(i-1)-1/(1+tau(3)/tau(2))*Bn(i-1);
end

a0=Cn(1);
temp1=zeros(Ny,Nx); temp2=zeros(Ny,Nx);
Star=[Star,zeros(i-2,1)];
Star=[Star;Cn'];

for s=2:i
        temp1=temp1+Cn(i-s+1)*(U1(:,:,s)-U1(:,:,s-1));
        temp2=temp2+Cn(i-s+1)*(U2(:,:,s)-U2(:,:,s-1));
end

Kmatrix=a0*ones(Ny,Nx)-nu*sigma*Kxxyy;

while error>ep  

    %中间点的非线性项
    U1_sigma=sigma*U10+(1-sigma)*U1(:,:,i-1);
    U2_sigma=sigma*U20+(1-sigma)*U2(:,:,i-1);
    F1=U1_sigma.*real(ifft2(Kx.*fft2(U1_sigma)))...
        +U2_sigma.*real(ifft2(Ky.*fft2(U1_sigma))); F1(1,1)=0;
    F2=U1_sigma.*real(ifft2(Kx.*fft2(U2_sigma)))...
        +U2_sigma.*real(ifft2(Ky.*fft2(U2_sigma))); F2(1,1)=0;

    P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));

    %计算新的迭代值
    U1_new=(-fft2(temp1)+nu*(1-sigma)*Kxxyy.*fft2(U1(:,:,i-1))...
        -fft2(F1)-Kx.*P_f)./Kmatrix;
    U2_new=(-fft2(temp2)+nu*(1-sigma)*Kxxyy.*fft2(U2(:,:,i-1))...
        -fft2(F2)-Ky.*P_f)./Kmatrix;

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
Evolu(:,:,i)=Kx.*fft2(U20)-Ky.*fft2(U10);
temp=Kx.*fft2(U10)+Ky.*fft2(U20);  temp=real(ifft2(temp));
Mass=[Mass;sum(sum(temp))];

%计算能量

Ph=zeros(i-1,1); Ph(1)=1/Star(i-1,1);
for j=2:i-1
    ss=0;
    for k=1:j-1
        ss=ss+(Star(i-k,j-k)-Star(i-k,j-k+1))*Ph(k);
    end
    Ph(j)=ss/Star(i-j,1);
end

bas_energy=(sum(sum(U10.^2))+sum(sum(U20.^2)))/2;
for k=1:i-1
    bas_energy=bas_energy+nu*Ph(i-k)*...
        (sum(sum(real(ifft2(Kx.*fft2(sigma*U1(:,:,k+1)+(1-sigma)*U1(:,:,k)))).^2))+...
        sum(sum(real(ifft2(Ky.*fft2(sigma*U2(:,:,k+1)+(1-sigma)*U2(:,:,k)))).^2)));
end

Energy=[Energy;(Ly*Lx/Nx/Ny)*bas_energy];

Tmesh(i)=Tmesh(i-1)+tau(i);

save("Evolution_for_ep1.mat","Evolu");
save("Tmesh_ep1.mat",Tmesh);
% end
