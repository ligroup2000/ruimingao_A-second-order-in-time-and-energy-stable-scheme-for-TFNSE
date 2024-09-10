%绘制能量变化
clear;clc;

Ly=pi; Lx=pi; Lz=pi; N=8; 
nu=0.1; alpha=0.8;
T=50; eta=1e1;
% 
[Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,Lz,N,nu,eta,alpha);
save energy81

alpha=0.5;

[Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,Lz,N,nu,eta,alpha);
save energy51
% 


load('energy51.mat')
E5=Energy; Tmesh5=Tmesh;
T5=Tmesh(2:end)-Tmesh(1:end-1); T5(end)=[];
M5=Mass;
E5(1)=[]; Tmesh5(1)=[]; M5(1)=[];

load('energy81.mat')
E8=Energy; Tmesh8=Tmesh;
T8=Tmesh(2:end)-Tmesh(1:end-1); T8(end)=[];
M8=Mass;
E8(1)=[]; Tmesh8(1)=[]; M8(1)=[];
% 
figure(1)
plot(Tmesh5,E5,'b-','Linewidth',2)
hold on
plot(Tmesh8,E8,'g-','Linewidth',2)
hold off
xlabel('Time')
ylabel('Energy')
set(gca,'FontSize',14)
legend('\alpha=0.5','\alpha=0.8')


figure(3)
subplot(1,2,1)
semilogy(Tmesh5,abs(M5))
set(gca,'YLim',[10^(-20) 10^(-10)])
xlabel('Time')
ylabel('$D(\textbf{u}_{N}^{n})$','Interpreter','latex')
set(gca,'FontSize',14)

subplot(1,2,2)
semilogy(Tmesh8,abs(M8))
set(gca,'YLim',[10^(-20) 10^(-10)])
xlabel('Time')
ylabel('$D(\textbf{u}_{N}^{n})$','Interpreter','latex')
set(gca,'FontSize',14)


figure(5)
plot(Tmesh5(2:end),T5,'b-','Linewidth',2)
hold on
plot(Tmesh8(2:end),T8,'g-','Linewidth',2)
hold off
legend('\alpha=0.5','\alpha=0.8')
xlabel('Time')
ylabel('Step-size')
set(gca,'FontSize',14)
