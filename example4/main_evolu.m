% function [Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,N,nu,eta,alpha,delta)
% function fracNS_Cal_energy
%% 给定初值
Ly=pi; Lx=pi; Lz=pi;
N=32; 
nu=5e-3; alpha=0.8; 
T=20;
eta=10;

up=Ly; bottom=-pi; 
right=Lx; left=-pi;
height=Lz; down=-pi;

Nx=N; Ny=N; Nz=N;
h1=abs(up-bottom)/Ny; 
h2=abs(right-left)/Nx;
h3=abs(height-down)/Nz;

sigma=1-alpha/2;
delta=2/alpha;

Ymesh=bottom:h1:up-h1;
Xmesh=left:h2:right-h2;
Zmesh=down:h3:height-h3;
sq=Ly*Lx*Lz/Nx/Ny/Nz;

%% 时间网格划分
% T=1;
T0=0.1; M=20; 
Tmesh=zeros(M+1,1);
for k=1:M+1
    Tmesh(k)=T0*power((k-1)/M,delta);
end
tau=[0;Tmesh(2:end)-Tmesh(1:end-1)];

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

%% 给定初值
U1=zeros(Nx,Ny,Nz,M+2); U2=zeros(Nx,Ny,Nz,M+2); U3=zeros(Nx,Ny,Nz,M+2);
P=zeros(Nx,Ny,Nz,M+2);

U1(:,:,:,1)=initial_u3(Xmesh,Ymesh,Zmesh,0,alpha);
U2(:,:,:,1)=initial_u2(Xmesh,Ymesh,Zmesh,0,alpha);
U3(:,:,:,1)=initial_u3(Xmesh,Ymesh,Zmesh,0,alpha);

U10=U1(:,:,:,1); U20=U2(:,:,:,1); U30=U3(:,:,:,1);
F1=U10.*real(ifftn(Kx.*fftn(U10)))+U20.*real(ifftn(Ky.*fftn(U10)))+U30.*real(ifftn(Kz.*fftn(U10)));
F2=U10.*real(ifftn(Kx.*fftn(U20)))+U20.*real(ifftn(Ky.*fftn(U20)))+U30.*real(ifftn(Kz.*fftn(U20)));
F3=U10.*real(ifftn(Kx.*fftn(U30)))+U20.*real(ifftn(Ky.*fftn(U30)))+U30.*real(ifftn(Kz.*fftn(U30)));

P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3));

P(:,:,:,1)=real(ifftn(P_f));temp=Kx.*fftn(U10)+Ky.*fftn(U20)+Kz.*fftn(U30); 
temp=real(ifftn(temp));
Mass=sum(sum(sum(temp)));

bas_energy=sq*(sum(sum(sum(U1(:,:,:,1).^2)))+sum(sum(sum(U2(:,:,:,1).^2)))+sum(sum(sum(U3(:,:,:,1).^2)))); 
Energy=bas_energy/2;

W1(:,:,:,1)=Ky.*fftn(U30)-Kz.*fftn(U20);
TT=0;
    
%% 计算第一步
t=sigma*Tmesh(2)+(1-sigma)*Tmesh(1);
An=zeros(2-1,1);
Star=zeros(2-1,1);
An(1)=1/(tau(2)*gamma(2-alpha))*power(t,1-alpha);
a0=An(1);
Star(1,:)=An;
Kmatrix=a0*ones(Ny,Nx,Nz)-nu*sigma*Kxxyy;
U10=U1(:,:,:,1); U20=U2(:,:,:,1); U30=U3(:,:,:,1);
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
    
    P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3));
    
    %计算新的迭代值
    U1_new=((a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U1(:,:,:,1))...
        -fftn(F1)-Kx.*P_f)./Kmatrix;
    U2_new=((a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U2(:,:,:,1))...
        -fftn(F2)-Ky.*P_f)./Kmatrix;
    U3_new=((a0*ones(Ny,Nx,Nz)+nu*(1-sigma)*Kxxyy).*fftn(U3(:,:,:,1))...
        -fftn(F3)-Kz.*P_f)./Kmatrix;
    
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
W1(:,:,:,2)=Ky.*fftn(U30)-Kz.*fftn(U20);
TT=[TT;Tmesh(2)];
Ph=zeros(2-1,1);
Ph(1)=1/An(1);    

bas_energy=sum(sum(sum(U1(:,:,:,2).^2)))+sum(sum(sum(U2(:,:,:,2).^2)))+sum(sum(sum(U3(:,:,:,2).^2)));
bas_energy=bas_energy/2;
for k=1:2-1
    bas_energy=bas_energy+nu*Ph(2-k)*...
        (sum(sum(sum(real(ifftn(Kx.*fftn(sigma*U1(:,:,:,k+1)+(1-sigma)*U1(:,:,:,k)))).^2)))+...
    sum(sum(sum(real(ifftn(Ky.*fftn(sigma*U2(:,:,:,k+1)+(1-sigma)*U2(:,:,:,k)))).^2)))+...
    sum(sum(sum(real(ifftn(Kz.*fftn(sigma*U3(:,:,:,k+1)+(1-sigma)*U3(:,:,:,k)))).^2))));
end

Energy=[Energy;sq*bas_energy];


temp=Kx.*fftn(U10)+Ky.*fftn(U20)+Kz.*fftn(U30); temp=real(ifftn(temp));
Mass=[Mass;sum(sum(sum(temp)))];

%% 循环
ep=1e-8; 
for i=3:M+1
    
    U10=U1(:,:,:,i-1); U20=U2(:,:,:,i-1); U30=U3(:,:,:,i-1);
    error=1; number=1;
    t=sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1);
    
    An=zeros(i-1,1); Bn=zeros(i-1,1);
    for k=1:i-1
        An(i-k)=1/(tau(k+1)*gamma(2-alpha))*(power(t-Tmesh(k),1-alpha)-power(t-min(t,Tmesh(k+1)),1-alpha));
        Bn(i-k)=2/(tau(k+1)^2*gamma(2-alpha))*((Tmesh(k)-(Tmesh(k)+Tmesh(k+1))/2)*power(t-Tmesh(k),1-alpha)-(Tmesh(k+1)-(Tmesh(k)+Tmesh(k+1))/2)*power(t-Tmesh(k+1),1-alpha))...
            +2/(tau(k+1)^2*gamma(3-alpha))*(power(t-Tmesh(k),2-alpha)-power(t-Tmesh(k+1),2-alpha));
    end
    Bn(1)=0;Bn=real(Bn);
    
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
    temp1=zeros(Ny,Nx,Nz); temp2=zeros(Ny,Nx,Nz); temp3=zeros(Ny,Nx,Nz);
    Star=[Star,zeros(i-2,1)];
    Star=[Star;Cn'];
    
    for s=2:i
        temp1=temp1+Cn(i-s+1)*(U1(:,:,:,s)-U1(:,:,:,s-1));
        temp2=temp2+Cn(i-s+1)*(U2(:,:,:,s)-U2(:,:,:,s-1));
        temp3=temp3+Cn(i-s+1)*(U3(:,:,:,s)-U3(:,:,:,s-1));
    end
    
    Kmatrix=a0*ones(Ny,Nx,Nz)-nu*sigma*Kxxyy;
       
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


        P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3));
        
        %计算新的迭代值
        U1_new=(-fftn(temp1)+nu*(1-sigma)*Kxxyy.*fftn(U1(:,:,:,i-1))...
            -fftn(F1)-Kx.*P_f)./Kmatrix;
        U2_new=(-fftn(temp2)+nu*(1-sigma)*Kxxyy.*fftn(U2(:,:,:,i-1))...
            -fftn(F2)-Ky.*P_f)./Kmatrix;
        U3_new=(-fftn(temp3)+nu*(1-sigma)*Kxxyy.*fftn(U3(:,:,:,i-1))...
            -fftn(F3)-Kz.*P_f)./Kmatrix;
    
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
    
    temp=Kx.*fftn(U10)+Ky.*fftn(U20)+Kz.*fftn(U30);  
    temp=real(ifftn(temp));
    
    Mass=[Mass;sum(sum(sum(temp)))];
    
    %计算能量
    
    Ph=zeros(i-1,1); Ph(1)=1/Star(i-1,1);
    for j=2:i-1
        ss=0;
        for k=1:j-1
            ss=ss+(Star(i-k,j-k)-Star(i-k,j-k+1))*Ph(k);
        end
        Ph(j)=ss/Star(i-j,1);
    end
    
    bas_energy=sum(sum(sum(U1(:,:,:,i).^2)))+sum(sum(sum(U2(:,:,:,i).^2)))+sum(sum(sum(U3(:,:,:,i).^2)));
    bas_energy=bas_energy/2;
    for k=1:i-1
        bas_energy=bas_energy+nu*Ph(i-k)*...
            (sum(sum(sum(real(ifftn(Kx.*fftn(sigma*U1(:,:,:,k+1)+(1-sigma)*U1(:,:,:,k)))).^2)))+...
        sum(sum(sum(real(ifftn(Ky.*fftn(sigma*U2(:,:,:,k+1)+(1-sigma)*U2(:,:,:,k)))).^2)))+...
        sum(sum(sum(real(ifftn(Kz.*fftn(sigma*U3(:,:,:,k+1)+(1-sigma)*U3(:,:,:,k)))).^2))));
    end
    Energy=[Energy;sq*bas_energy];
    
end

W1(:,:,:,3)=Ky.*fftn(U30)-Kz.*fftn(U20);
TT=[TT;Tmesh(M+1)];
%% 变步长
tau_min=1e-3; tau_max=1e-1; i=M+1; n_ev=3;
while (Tmesh(i)<T)
    
    i=i+1; 
    U1(:,:,:,i)=zeros(Nx,Ny,Nz); 
    U2(:,:,:,i)=zeros(Nx,Ny,Nz); 
    U3(:,:,:,i)=zeros(Nx,Ny,Nz);
    U10=U1(:,:,:,i-1); U20=U2(:,:,:,i-1); U30=U3(:,:,:,i-1);
    
    error=1; number=1;
    sub=abs(Energy(i-1)-Energy(i-2))^2/(tau(i-1)^2);
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
    temp1=zeros(Ny,Nx,Nz); temp2=zeros(Ny,Nx,Nz); temp3=zeros(Ny,Nx,Nz);
    Star=[Star,zeros(i-2,1)];
    Star=[Star;Cn'];
    
    for s=2:i
        temp1=temp1+Cn(i-s+1)*(U1(:,:,:,s)-U1(:,:,:,s-1));
        temp2=temp2+Cn(i-s+1)*(U2(:,:,:,s)-U2(:,:,:,s-1));
        temp3=temp3+Cn(i-s+1)*(U3(:,:,:,s)-U3(:,:,:,s-1));
    end
    
    Kmatrix=a0*ones(Ny,Nx,Nz)-nu*sigma*Kxxyy;
       
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


        P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3));
        
        %计算新的迭代值
        U1_new=(-fftn(temp1)+nu*(1-sigma)*Kxxyy.*fftn(U1(:,:,:,i-1))...
            -fftn(F1)-Kx.*P_f)./Kmatrix;
        U2_new=(-fftn(temp2)+nu*(1-sigma)*Kxxyy.*fftn(U2(:,:,:,i-1))...
            -fftn(F2)-Ky.*P_f)./Kmatrix;
        U3_new=(-fftn(temp3)+nu*(1-sigma)*Kxxyy.*fftn(U3(:,:,:,i-1))...
            -fftn(F3)-Kz.*P_f)./Kmatrix;
    
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
    
    temp=Kx.*fftn(U10)+Ky.*fftn(U20)+Kz.*fftn(U30);  
    temp=real(ifftn(temp));
    
    Mass=[Mass;sum(sum(sum(temp)))];    
    
    %计算能量
    Ph=zeros(i-1,1); Ph(1)=1/Star(i-1,1);
    for j=2:i-1
        ss=0;
        for k=1:j-1
            ss=ss+(Star(i-k,j-k)-Star(i-k,j-k+1))*Ph(k);
        end
        Ph(j)=ss/Star(i-j,1);
    end
    
    bas_energy=sum(sum(sum(U1(:,:,:,i).^2)))+sum(sum(sum(U2(:,:,:,i).^2)))+sum(sum(sum(U3(:,:,:,i).^2)));
    bas_energy=bas_energy/2;
    UU1=sigma*U1(:,:,:,2:end)+(1-sigma)*U1(:,:,:,1:end-1);
    UU2=sigma*U2(:,:,:,2:end)+(1-sigma)*U2(:,:,:,1:end-1);
    UU3=sigma*U3(:,:,:,2:end)+(1-sigma)*U3(:,:,:,1:end-1);
    for k=1:i-1
        bas_energy=bas_energy+nu*Ph(i-k)*...
            (sum(sum(sum(real(ifftn(Kx.*fftn(UU1(:,:,:,k)))).^2)))+...
        sum(sum(sum(real(ifftn(Ky.*fftn(UU2(:,:,:,k)))).^2)))+...
        sum(sum(sum(real(ifftn(Kz.*fftn(UU3(:,:,:,k)))).^2))));
    end
    Energy=[Energy;sq*bas_energy];
    
    Tmesh(i)=Tmesh(i-1)+tau(i);
    if mod(i,10)==1
        n_ev=n_ev+1;
        W1(:,:,:,n_ev)=Ky.*fftn(U30)-Kz.*fftn(U20);
        TT=[TT;Tmesh(i)];
    end
    
end

%% 计算最后一步
i=length(Tmesh); tau(end)=T-Tmesh(i);  
Tmesh(i)=T;
t=sigma*Tmesh(i)+(1-sigma)*Tmesh(i-1);

U10=U1(:,:,:,i-1); U20=U2(:,:,:,i-1); U30=U3(:,:,:,i-1);
U1(:,:,:,i)=zeros(Nx,Ny,Nz); 
U2(:,:,:,i)=zeros(Nx,Ny,Nz); 
U3(:,:,:,i)=zeros(Nx,Ny,Nz);
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
temp1=zeros(Ny,Nx,Nz); temp2=zeros(Ny,Nx,Nz); temp3=zeros(Ny,Nx,Nz);
Star=[Star,zeros(i-2,1)];
Star=[Star;Cn'];

for s=2:i
    temp1=temp1+Cn(i-s+1)*(U1(:,:,:,s)-U1(:,:,:,s-1));
    temp2=temp2+Cn(i-s+1)*(U2(:,:,:,s)-U2(:,:,:,s-1));
    temp3=temp3+Cn(i-s+1)*(U3(:,:,:,s)-U3(:,:,:,s-1));
end

Kmatrix=a0*ones(Ny,Nx,Nz)-nu*sigma*Kxxyy;

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


    P_f=-Khat.*(Kx.*fftn(F1)+Ky.*fftn(F2)+Kz.*fftn(F3));

    %计算新的迭代值
    U1_new=(-fftn(temp1)+nu*(1-sigma)*Kxxyy.*fftn(U1(:,:,:,i-1))...
        -fftn(F1)-Kx.*P_f)./Kmatrix;
    U2_new=(-fftn(temp2)+nu*(1-sigma)*Kxxyy.*fftn(U2(:,:,:,i-1))...
        -fftn(F2)-Ky.*P_f)./Kmatrix;
    U3_new=(-fftn(temp3)+nu*(1-sigma)*Kxxyy.*fftn(U3(:,:,:,i-1))...
        -fftn(F3)-Kz.*P_f)./Kmatrix;

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
    
temp=Kx.*fftn(U10)+Ky.*fftn(U20)+Kz.*fftn(U30);  
temp=real(ifftn(temp));

Mass=[Mass;sum(sum(sum(temp)))];   

%计算能量

Ph=zeros(i-1,1); Ph(1)=1/Star(i-1,1);
for j=2:i-1
    ss=0;
    for k=1:j-1
        ss=ss+(Star(i-k,j-k)-Star(i-k,j-k+1))*Ph(k);
    end
    Ph(j)=ss/Star(i-j,1);
end

bas_energy=sum(sum(sum(U1(:,:,:,i).^2)))+sum(sum(sum(U2(:,:,:,i).^2)))+sum(sum(sum(U3(:,:,:,i).^2)));
bas_energy=bas_energy/2;
for k=1:i-1
    bas_energy=bas_energy+nu*Ph(i-k)*...
        (sum(sum(sum(real(ifftn(Kx.*fftn(sigma*U1(:,:,:,k+1)+(1-sigma)*U1(:,:,:,k)))).^2)))+...
    sum(sum(sum(real(ifftn(Ky.*fftn(sigma*U2(:,:,:,k+1)+(1-sigma)*U2(:,:,:,k)))).^2)))+...
    sum(sum(sum(real(ifftn(Kz.*fftn(sigma*U3(:,:,:,k+1)+(1-sigma)*U3(:,:,:,k)))).^2))));
end
Energy=[Energy;sq*bas_energy];

Tmesh(i)=Tmesh(i-1)+tau(i);

n_ev=n_ev+1;
W1(:,:,:,n_ev)=Ky.*fftn(U30)-Kz.*fftn(U20);
TT=[TT;Tmesh(i)];

save("W1.mat","W1");
save("TT.mat","TT");
% end
