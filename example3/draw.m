% clear;clc;
% 
% load('Evolution_for_ep1.mat')
load('Tmesh_ep1.mat')

Ly=1; Lx=1; N=128; 

up=Ly; bottom=0; right=Lx; left=0;
Nx=N; Ny=N;
h1=abs(up-bottom)/Ny; h2=abs(right-left)/Nx;

Ymesh=bottom:h1:up-h1;
Xmesh=left:h2:right-h2;

K=[1,500,1000,1550,2550,4000];
subplot(2,3,1)
k=K(1);
u11=real(ifft2(Evolu(:,:,k)));
disp(Tmesh(k))

pcolor(Xmesh,Ymesh,u11);
shading interp;
title('T=0');

subplot(2,3,2)
k=K(2);
u11=real(ifft2(Evolu(:,:,k)));
disp(Tmesh(k))

pcolor(Xmesh,Ymesh,u11);
shading interp;
title('T=0.27');

subplot(2,3,3)
k=K(3);
u11=real(ifft2(Evolu(:,:,k)));
disp(Tmesh(k))

pcolor(Xmesh,Ymesh,u11);
shading interp;
title('T=0.37');

subplot(2,3,4)
k=K(4);
u11=real(ifft2(Evolu(:,:,k)));
disp(Tmesh(k))

pcolor(Xmesh,Ymesh,u11);
shading interp;
title('T=0.44');

subplot(2,3,5)
k=K(5);
u11=real(ifft2(Evolu(:,:,k)));
disp(Tmesh(k))

pcolor(Xmesh,Ymesh,u11);
shading interp;
title('T=0.55');

subplot(2,3,6)
k=K(6);
u11=real(ifft2(Evolu(:,:,k)));
disp(Tmesh(k))

pcolor(Xmesh,Ymesh,u11);
shading interp;
title('T=0.70');
