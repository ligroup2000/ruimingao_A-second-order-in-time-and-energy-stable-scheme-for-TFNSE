clear;clc;

Ly=2*pi; Lx=2*pi; Lz=2*pi;
N=8; 
nu=0.1; 
T=50;
eta=1;

alpha=0.8;
[Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,Lz,N,nu,eta,alpha);
save("E8.mat","Energy");
save("M8.mat","Mass");
save("T8.mat","Tmesh");

alpha=0.6;

[Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,Lz,N,nu,eta,alpha);
save("E6.mat","Energy");
save("M6.mat","Mass");
save("T6.mat","Tmesh");

alpha=0.4;

[Energy,Mass,Tmesh]=fracNS_Cal_energy(T,Ly,Lx,Lz,N,nu,eta,alpha);
save("E4.mat","Energy");
save("M4.mat","Mass");
save("T4.mat","Tmesh");
% 

load('E4.mat'); load('M4.mat'); load('T4.mat');
E4=Energy; Tmesh4=Tmesh;
T4=Tmesh(2:end)-Tmesh(1:end-1); T4(end)=[];
M4=Mass;
E4(1)=[]; Tmesh4(1)=[]; M4(1)=[];

load('E6.mat'); load('M6.mat'); load('T6.mat');
E6=Energy; Tmesh6=Tmesh;
T6=Tmesh(2:end)-Tmesh(1:end-1); T6(end)=[];
M6=Mass;
E6(1)=[]; Tmesh6(1)=[]; M6(1)=[];

load('E8.mat'); load('M8.mat'); load('T8.mat');
E8=Energy; Tmesh8=Tmesh;
T8=Tmesh(2:end)-Tmesh(1:end-1); T8(end)=[];
M8=Mass;
E8(1)=[]; Tmesh8(1)=[]; M8(1)=[];
% 
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
% 
axes('Position',[0.20,0.40,0.2,0.4]);
plot(Tmesh3(1:208),E3(1:208),'r-','Linewidth',2)
hold on
plot(Tmesh5(1:372),E5(1:372),'b-','Linewidth',2)
hold on
plot(Tmesh8(1:607),E8(1:607),'g-','Linewidth',2)
%set(gca,'YTick',[0:1:2]);
legend('\alpha=0.3','\alpha=0.5','\alpha=0.8')
set(gca,'FontSize',8)


% 
figure(2)
subplot(1,3,1)
semilogy(Tmesh4,abs(M4))
set(gca,'YLim',[10^(-25) 10^(-10)])
xlabel('Time')
ylabel('$D(\textbf{u}_{N}^{n})$','Interpreter','latex')
set(gca,'FontSize',14)

subplot(1,3,2)
semilogy(Tmesh6,abs(M6))
set(gca,'YLim',[10^(-25) 10^(-10)])
xlabel('Time')
ylabel('$D(\textbf{u}_{N}^{n})$','Interpreter','latex')
set(gca,'FontSize',14)


subplot(1,3,3)
semilogy(Tmesh8,abs(M8))
set(gca,'YLim',[10^(-25) 10^(-10)])
xlabel('Time')
ylabel('$D(\textbf{u}_{N}^{n})$','Interpreter','latex')
set(gca,'FontSize',14)



figure(5)
plot(Tmesh4(2:end),T4,'r-','Linewidth',1)
hold on
plot(Tmesh6(2:end),T6,'b-','Linewidth',1)
hold on
plot(Tmesh8(2:end),T8,'g-','Linewidth',1)
hold off
xlabel('time')
ylabel('Step-size')
set(gca,'FontSize',14)
set(gca,'YLim',[10^(-4) 10^(0)])
axes('Position',[0.25,0.45,0.2,0.4]);
plot(Tmesh4(1:25),T4(1:25),'r-','Linewidth',2)
hold on
plot(Tmesh6(1:25),T6(1:25),'b-','Linewidth',2)
hold on
plot(Tmesh8(1:25),T8(1:25),'g-','Linewidth',2)
%set(gca,'YTick',[0:1:2]);
legend('\alpha=0.4','\alpha=0.6','\alpha=0.8')
set(gca,'FontSize',8)
legend('\alpha=0.4','\alpha=0.6','\alpha=0.8')

