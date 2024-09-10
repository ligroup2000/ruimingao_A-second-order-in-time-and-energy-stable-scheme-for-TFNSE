% load('TT.mat')
% load('W1.mat')

%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 1;     % ��ͼ����
TightPlot.RowNumber = 1;    % ��ͼ����
TightPlot.GapW = 0.05;  % ��ͼ֮������Ҽ��
TightPlot.GapH = 0.05;   % ��ͼ֮������¼��
TightPlot.MarginsLower = 0.11;   % ��ͼ��ͼƬ�·��ļ��
TightPlot.MarginsUpper = 0.11;  % ��ͼ��ͼƬ�Ϸ��ļ��
TightPlot.MarginsLeft = 0.24;   % ��ͼ��ͼƬ�󷽵ļ��
TightPlot.MarginsRight = 0.41;  % ��ͼ��ͼƬ�ҷ��ļ��

%% drawing
Ly=2*pi; Lx=2*pi; Lz=2*pi;
N=64;

up=Ly; bottom=0; 
right=Lx; left=0;
height=Lz; down=0;

Nx=N; Ny=N; Nz=N;
h1=abs(up-bottom)/Ny; 
h2=abs(right-left)/Nx;
h3=abs(height-down)/Nz;

Ymesh=bottom:h1:up;
Xmesh=left:h2:right;
Zmesh=down:h3:height;

[Xmesh,Ymesh,Zmesh]=meshgrid(Xmesh,Ymesh,Zmesh);

k=9;
W=real(ifftn(W1(:,:,:,k)));
disp(TT(k))

figure(1)
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
isosurface(Xmesh,Ymesh,Zmesh,WW,1);
set(gca,'XLim',[0,2*pi])
set(gca,'YLim',[0,2*pi])
set(gca,'ZLim',[0,2*pi])
set(gca,'xtick',[],'ytick',[],'ztick',[])
xlabel('$x$','interpreter','latex','FontSize',20)
ylabel('$y$','interpreter','latex','FontSize',20)
zlabel('$z$','interpreter','latex','FontSize',20)
box on;  
view(-60,30);
brighten(0.5);       %����
camlight right;      %��Դλ��
lighting phong;      %����ģʽ