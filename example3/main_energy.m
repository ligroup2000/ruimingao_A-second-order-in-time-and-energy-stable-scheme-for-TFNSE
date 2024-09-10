%绘制能量变化
clear;clc;

Ly=1; Lx=1; N=8; 
nu=0.1; 
T=50; eta=1e1;

alpha=0.8;
[Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,N,nu,eta,alpha);
save energy81

alpha=0.6;

[Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,N,nu,eta,alpha);
save energy61

alpha=0.4;

[Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,N,nu,eta,alpha);
save energy41

load('energy41.mat');
E4=Energy; Tmesh4=Tmesh;
T4=Tmesh(2:end)-Tmesh(1:end-1); T4(end)=[];
M4=Mass;
E4(1)=[]; Tmesh4(1)=[]; M4(1)=[];

load('energy61.mat'); 
E6=Energy; Tmesh6=Tmesh;
T6=Tmesh(2:end)-Tmesh(1:end-1); T6(end)=[];
M6=Mass;
E6(1)=[]; Tmesh6(1)=[]; M6(1)=[];

load('energy81.mat');
E8=Energy; Tmesh8=Tmesh;
T8=Tmesh(2:end)-Tmesh(1:end-1); T8(end)=[];
M8=Mass;
E8(1)=[]; Tmesh8(1)=[]; M8(1)=[];

figure(1)
plot(Tmesh4,E4,'r-','Linewidth',2)
hold on
plot(Tmesh6,E6,'b-','Linewidth',2)
hold on
plot(Tmesh8,E8,'g-','Linewidth',2)
hold off
xlabel('Time')
ylabel('Energy')
set(gca,'FontSize',14)
legend('\alpha=0.4','\alpha=0.6','\alpha=0.8')
% % 
axes('Position',[0.20,0.40,0.2,0.4]);
plot(Tmesh4(1:77),E4(1:77),'r-','Linewidth',2)
hold on
plot(Tmesh6(1:100),E6(1:100),'b-','Linewidth',2)
hold on
plot(Tmesh8(1:110),E8(1:110),'g-','Linewidth',2)
%set(gca,'YTick',[0:1:2]);
legend('\alpha=0.4','\alpha=0.6','\alpha=0.8')
set(gca,'FontSize',8)


figure(2)
subplot(1,3,1)
semilogy(Tmesh4(1:end),abs(M4(1:end)))
set(gca,'YLim',[10^(-25) 10^(-16)])
xlabel('Time')
ylabel('$D(\textbf{u}_{N}^{n})$','Interpreter','latex')
set(gca,'FontSize',14)

subplot(1,3,2)
semilogy(Tmesh6,abs(M6))
set(gca,'YLim',[10^(-25) 10^(-16)])
xlabel('Time')
ylabel('$D(\textbf{u}_{N}^{n})$','Interpreter','latex')
set(gca,'FontSize',14)

subplot(1,3,3)
semilogy(Tmesh8,abs(M8))
set(gca,'YLim',[10^(-25) 10^(-16)])
xlabel('Time')
ylabel('$D(\textbf{u}_{N}^{n})$','Interpreter','latex')
set(gca,'FontSize',14)

% 
% 
figure(5)
plot(Tmesh4(2:end),T4,'r-','Linewidth',1)
hold on
plot(Tmesh6(2:end),T6,'b-','Linewidth',2)
hold on
plot(Tmesh8(2:end),T8,'g-','Linewidth',2)
hold off
% set(gca,'YLim',[10^(-4) 10^(0)])
legend('\alpha=0.4','\alpha=0.6','\alpha=0.8')
xlabel('Time')
ylabel('Step-size')
set(gca,'FontSize',14)