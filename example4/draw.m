load('TT.mat')
load('W1.mat')

%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 1;     % ��ͼ����
TightPlot.RowNumber = 1;    % ��ͼ����
TightPlot.GapW = 0.1;  % ��ͼ֮������Ҽ��
TightPlot.GapH = 0.1;   % ��ͼ֮������¼��
TightPlot.MarginsLower = 0.1;   % ��ͼ��ͼƬ�·��ļ��
TightPlot.MarginsUpper = 0.1;  % ��ͼ��ͼƬ�Ϸ��ļ��
TightPlot.MarginsLeft = 0.2;   % ��ͼ��ͼƬ�󷽵ļ��
TightPlot.MarginsRight = 0.2;  % ��ͼ��ͼƬ�ҷ��ļ��

%% drawing
Ly=pi; Lx=pi; Lz=pi;
N=32;

up=Ly; bottom=-pi; 
right=Lx; left=-pi;
height=Lz; down=-pi;

Nx=N; Ny=N; Nz=N;
h1=abs(up-bottom)/Ny; 
h2=abs(right-left)/Nx;
h3=abs(height-down)/Nz;

Ymesh=bottom:h1:up;
Xmesh=left:h2:right;
Zmesh=down:h3:height;

[Xmesh,Ymesh,Zmesh]=meshgrid(Xmesh,Ymesh,Zmesh);
K=[1,12,14,15,16,18];


k=K(1);
W=real(ifftn(W1(:,:,:,k)));
disp(TT(k))


WW=zeros(N+1,N+1,N+1);
WW(1:N,1:N,1:N)=W;  
WW(1:N,1:N,N+1)=W(1:N,1:N,1);  
WW(1:N,1+N,1:N+1)=WW(1:N,1,1:N+1);  
WW(1+N,1:N+1,1:N+1)=WW(1,1:N+1,1:N+1);  
WW(N+1,N+1,1:N+1)=WW(1,1,1:N+1);
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]); 
subplot(3,2,1)
isosurface(Xmesh,Ymesh,Zmesh,WW,1);
set(gca,'XLim',[-pi,pi])
set(gca,'YLim',[-pi,pi])
set(gca,'ZLim',[-pi,pi])
set(gca,'xtick',[],'ytick',[],'ztick',[])
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$y$','interpreter','latex','FontSize',15)
zlabel('$z$','interpreter','latex','FontSize',15)
box on;  
view(-60,30);
brighten(0.5);       %����
camlight right;      %��Դλ��
lighting phong;      %����ģʽ
title('T=0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=K(2);
W=real(ifftn(W1(:,:,:,k)));
disp(TT(k))

WW=zeros(N+1,N+1,N+1);
WW(1:N,1:N,1:N)=W;  
WW(1:N,1:N,N+1)=W(1:N,1:N,1);  
WW(1:N,1+N,1:N+1)=WW(1:N,1,1:N+1);  
WW(1+N,1:N+1,1:N+1)=WW(1,1:N+1,1:N+1);  
WW(N+1,N+1,1:N+1)=WW(1,1,1:N+1);
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]); 
subplot(3,2,2)
isosurface(Xmesh,Ymesh,Zmesh,WW,1);
set(gca,'XLim',[-pi,pi])
set(gca,'YLim',[-pi,pi])
set(gca,'ZLim',[-pi,pi])
set(gca,'xtick',[],'ytick',[],'ztick',[])
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$y$','interpreter','latex','FontSize',15)
zlabel('$z$','interpreter','latex','FontSize',15)
box on;  
view(-60,30);
% brighten(0.5);       %����
camlight right;      %��Դλ��
lighting phong;      %����ģʽ
title('T=1.55')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=K(3);
W=real(ifftn(W1(:,:,:,k)));
disp(TT(k))

WW=zeros(N+1,N+1,N+1);
WW(1:N,1:N,1:N)=W;  
WW(1:N,1:N,N+1)=W(1:N,1:N,1);  
WW(1:N,1+N,1:N+1)=WW(1:N,1,1:N+1);  
WW(1+N,1:N+1,1:N+1)=WW(1,1:N+1,1:N+1);  
WW(N+1,N+1,1:N+1)=WW(1,1,1:N+1);
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]); 
subplot(3,2,3)
isosurface(Xmesh,Ymesh,Zmesh,WW,1);
set(gca,'XLim',[-pi,pi])
set(gca,'YLim',[-pi,pi])
set(gca,'ZLim',[-pi,pi])
set(gca,'xtick',[],'ytick',[],'ztick',[])
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$y$','interpreter','latex','FontSize',15)
zlabel('$z$','interpreter','latex','FontSize',15)
box on;  
view(-60,30);
% brighten(0.5);       %����
camlight right;      %��Դλ��
lighting phong;      %����ģʽ
title('T=3.14')

%%%%%%%%%%%%%%%%%%%%%%%%%%
k=K(4);
W=real(ifftn(W1(:,:,:,k)));
disp(TT(k))

WW=zeros(N+1,N+1,N+1);
WW(1:N,1:N,1:N)=W;  
WW(1:N,1:N,N+1)=W(1:N,1:N,1);  
WW(1:N,1+N,1:N+1)=WW(1:N,1,1:N+1);  
WW(1+N,1:N+1,1:N+1)=WW(1,1:N+1,1:N+1);  
WW(N+1,N+1,1:N+1)=WW(1,1,1:N+1);
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]); 
subplot(3,2,4)
isosurface(Xmesh,Ymesh,Zmesh,WW,1);
set(gca,'XLim',[-pi,pi])
set(gca,'YLim',[-pi,pi])
set(gca,'ZLim',[-pi,pi])
set(gca,'xtick',[],'ytick',[],'ztick',[])
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$y$','interpreter','latex','FontSize',15)
zlabel('$z$','interpreter','latex','FontSize',15)
box on;  
view(-60,30);
% brighten(0.5);       %����
camlight right;      %��Դλ��
lighting phong;      %����ģʽ
title('T=4.07')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=K(5);
W=real(ifftn(W1(:,:,:,k)));
disp(TT(k))

WW=zeros(N+1,N+1,N+1);
WW(1:N,1:N,1:N)=W;  
WW(1:N,1:N,N+1)=W(1:N,1:N,1);  
WW(1:N,1+N,1:N+1)=WW(1:N,1,1:N+1);  
WW(1+N,1:N+1,1:N+1)=WW(1,1:N+1,1:N+1);  
WW(N+1,N+1,1:N+1)=WW(1,1,1:N+1);
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]); 
subplot(3,2,5)
isosurface(Xmesh,Ymesh,Zmesh,WW,1);
set(gca,'XLim',[-pi,pi])
set(gca,'YLim',[-pi,pi])
set(gca,'ZLim',[-pi,pi])
set(gca,'xtick',[],'ytick',[],'ztick',[])
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$y$','interpreter','latex','FontSize',15)
zlabel('$z$','interpreter','latex','FontSize',15)
box on;  
view(-60,30);
% brighten(0.5);       %����
camlight right;      %��Դλ��
lighting phong;      %����ģʽ
title('T=5.03')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=K(6);
W=real(ifftn(W1(:,:,:,k)));
disp(TT(k))

WW=zeros(N+1,N+1,N+1);
WW(1:N,1:N,1:N)=W;  
WW(1:N,1:N,N+1)=W(1:N,1:N,1);  
WW(1:N,1+N,1:N+1)=WW(1:N,1,1:N+1);  
WW(1+N,1:N+1,1:N+1)=WW(1,1:N+1,1:N+1);  
WW(N+1,N+1,1:N+1)=WW(1,1,1:N+1);
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]); 
subplot(3,2,6)
isosurface(Xmesh,Ymesh,Zmesh,WW,1);
set(gca,'XLim',[-pi,pi])
set(gca,'YLim',[-pi,pi])
set(gca,'ZLim',[-pi,pi])
set(gca,'xtick',[],'ytick',[],'ztick',[])
xlabel('$x$','interpreter','latex','FontSize',15)
ylabel('$y$','interpreter','latex','FontSize',15)
zlabel('$z$','interpreter','latex','FontSize',15)
box on;  
view(-60,30);
% brighten(0.5);       %����
camlight right;      %��Դλ��
lighting phong;      %����ģʽ
title('T=7.00')