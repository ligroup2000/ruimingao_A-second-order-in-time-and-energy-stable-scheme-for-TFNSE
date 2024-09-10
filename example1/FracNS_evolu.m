function [U1,U2,Tmesh]=FracNS_evolu(T,Ly,Lx,N,nu,eta,alpha)

%% 给定初值
Ly=2; Lx=2; N=8; delta=1;
nu=0.001; alpha=0.5;
delta=(2-alpha)/alpha;
up=Ly; bottom=0; right=Lx; left=0;
Nx=N; Ny=N;
h1=abs(up-bottom)/Ny; h2=abs(right-left)/Nx;

Ymesh=bottom:h1:up-h1;
Xmesh=left:h2:right-h2;

%% 时间网格划分
% T=1;
T0=0.01; M=3e2; 
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
U1(:,:,1)=initial_u1(Xmesh,Ymesh,Tmesh,alpha);
U2(:,:,1)=initial_u2(Xmesh,Ymesh,Tmesh,alpha);
    
%% 循环
ep=1e-8; 
for i=2:M+1
    
    U10=U1(:,:,i-1); U20=U2(:,:,i-1); 
    error=1; number=1;
    An=zeros(i-1,1);
    
    for k=1:i-1
      An(i-k)=(power(Tmesh(i)-Tmesh(k),1-alpha)-power(Tmesh(i)-Tmesh(k+1),1-alpha))/(gamma(2-alpha)*tau(k+1));
    end
        
    temp1=zeros(Nx,Ny); temp2=temp1;
    for k=2:i
        temp1=temp1+An(i-k+1)*(fft2(U1(:,:,k))-fft2(U1(:,:,k-1)));
        temp2=temp2+An(i-k+1)*(fft2(U2(:,:,k))-fft2(U2(:,:,k-1)));
    end
    
    Kmatrix=An(1)*ones(Nx,Ny)-nu*Kxxyy;
    
    while error>ep
        
        F1=U10.*real(ifft2(Kx.*fft2(U10)))+U20.*real(ifft2(Ky.*fft2(U10))); F1(1,1)=0;
        F2=U10.*real(ifft2(Kx.*fft2(U20)))+U20.*real(ifft2(Ky.*fft2(U20))); F2(1,1)=0;
        P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));
        
        U1_new=(-temp1-fft2(F1)-Kx.*P_f)./Kmatrix;
        U2_new=(-temp2-fft2(F2)-Ky.*P_f)./Kmatrix;        
        error=real(max(max(abs([fft2(U10)-U1_new;fft2(U20)-U2_new]))));        
        
        U1_new=real(ifft2(U1_new)); U2_new=real(ifft2(U2_new));
        U10=U1_new; U20=U2_new;
        
        number=number+1;
        if number>100
%             error
            disp(i);
            break;
        end
    end
        
    U1(:,:,i)=U10; U2(:,:,i)=U20;   
    disp(Tmesh(end))
    pcolor(U10);
    shading interp
    pause(0.001)
end

U1(:,:,end)=[]; U2(:,:,end)=[]; P(:,:,end)=[];
%% 变步长
tau_min=1e-4; tau_max=1e-3; i=M+1;
while (Tmesh(i)<T)
    
    i=i+1; U1(:,:,i)=zeros(Nx,Ny); U2(:,:,i)=zeros(Nx,Ny);
    U10=U1(:,:,i-1); U20=U2(:,:,i-1); 
    error=1; number=1;
    An=zeros(i-1,1);
    sub=sum(sum((U1(:,:,i-1)-U1(:,:,i-2)).^2))+sum(sum((U2(:,:,i-1)-U2(:,:,i-2)).^2));
    sub=sub/(tau(i-1)^2);
    tau=[tau;max(tau_min,tau_max/(sqrt(1+eta*sub)))];
    Tmesh=[Tmesh;Tmesh(end)+tau(end)];
    
    if Tmesh(end)+tau(end)>T
        break;
    end
    
    for k=1:i-1
      An(i-k)=(power(Tmesh(i)-Tmesh(k),1-alpha)-power(Tmesh(i)-Tmesh(k+1),1-alpha))/(gamma(2-alpha)*tau(k+1));
    end
        
    temp1=zeros(Nx,Ny); temp2=temp1;
    for k=2:i
        temp1=temp1+An(i-k+1)*(fft2(U1(:,:,k))-fft2(U1(:,:,k-1)));
        temp2=temp2+An(i-k+1)*(fft2(U2(:,:,k))-fft2(U2(:,:,k-1)));
    end
    
    Kmatrix=An(1)*ones(Nx,Ny)-nu*Kxxyy;
    
    while error>ep
        
        F1=U10.*real(ifft2(Kx.*fft2(U10)))+U20.*real(ifft2(Ky.*fft2(U10))); F1(1,1)=0;
        F2=U10.*real(ifft2(Kx.*fft2(U20)))+U20.*real(ifft2(Ky.*fft2(U20))); F2(1,1)=0;
        P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));
        
        U1_new=(-temp1-fft2(F1)-Kx.*P_f)./Kmatrix;
        U2_new=(-temp2-fft2(F2)-Ky.*P_f)./Kmatrix;
        U1_new=real(ifft2(U1_new)); U2_new=real(ifft2(U2_new));
        
        error=max(max(abs([U10-U1_new;U20-U2_new])));
        U10=U1_new; U20=U2_new;
        
        number=number+1;
        if number>100
            error
            disp(i);
            break;
        end
    end
        
    U1(:,:,i)=U10; U2(:,:,i)=U20; 
    
    disp(Tmesh(end))
    pcolor(U10);
    shading interp
    pause(0.001)
end

%% 计算最后一步

i=length(Tmesh); tau(end)=T-Tmesh(i);  
Tmesh(i)=T;

U10=U1(:,:,i-1); U20=U2(:,:,i-1); U1(:,:,i)=zeros(Nx,Ny); U2(:,:,i)=zeros(Nx,Ny);
error=1; number=1;
An=zeros(i-1,1);

for k=1:i-1
  An(i-k)=(power(Tmesh(i)-Tmesh(k),1-alpha)-power(Tmesh(i)-Tmesh(k+1),1-alpha))/(gamma(2-alpha)*tau(k+1));
end

temp1=zeros(Nx,Ny); temp2=temp1;
for k=2:i
    temp1=temp1+An(i-k+1)*(fft2(U1(:,:,k))-fft2(U1(:,:,k-1)));
    temp2=temp2+An(i-k+1)*(fft2(U2(:,:,k))-fft2(U2(:,:,k-1)));
end

Kmatrix=An(1)*ones(Nx,Ny)-nu*Kxxyy;

while error>ep

    F1=U10.*real(ifft2(Kx.*fft2(U10)))+U20.*real(ifft2(Ky.*fft2(U10))); F1(1,1)=0;
    F2=U10.*real(ifft2(Kx.*fft2(U20)))+U20.*real(ifft2(Ky.*fft2(U20))); F2(1,1)=0;
    P_f=-Khat.*(Kx.*fft2(F1)+Ky.*fft2(F2));

    U1_new=(-temp1-fft2(F1)-Kx.*P_f)./Kmatrix;
    U2_new=(-temp2-fft2(F2)-Ky.*P_f)./Kmatrix;
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
disp(Tmesh(end))
pcolor(U10);
shading interp
pause(0.001)

end

