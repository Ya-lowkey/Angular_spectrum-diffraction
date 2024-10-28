%Script file: Angular Spectrum diffraction
% 
%Purpose :
%    By numerical simulation, the diffraction pattern of Fresnel approximation 
%    is calculated by Angular Spectrum method at different receiving plane distances
%
%Record of revisions:
%                Date                 Programmer            Description of change
%              =========             =============          =====================
%              01/29/23               梦想是优秀社畜              Original code
%
%% 
clc;clear;set(0,'defaultfigurecolor','w'); %将画布默认底色改为白色
milli = 1e-3;micron = 1e-6;nano = 1e-9; % 设置长度单位（毫米 微米 纳米）
%% 衍射屏形状
N = 200; % 采样点数
imag = imread('cloud.bmp'); %读取不规则形状的光孔图片
m = 22; % 孔径大小
imag = imresize(imag,[m m]); %设置孔径大小设为22*22（mm）
imag = padarray(imag,[(N-m)/2,(N-m)/2],255); %在矩阵周围补值（取纯白），使其符合衍射屏大小
imag = im2bw(imag,0.5); %设定阈值使扩展补值后矩阵通孔处值为0，其余为255（二值图像）
imag = ~imag;
imag = double(imag);

%% 参数设置
pix = 4.65*micron; %像素宽(mm)
lambda = 532*nano; % 入射波长
k = 2*pi/lambda; % 波数
z1 = 10*milli; %观察屏距离1(mm)
z2 = 15*milli; %观察屏距离2(mm)

%% 坐标系建立
originpoint = ceil((N+1)/2); 
totalsize = pix*N;
% 空间坐标系
xaxis=(1-originpoint:N-originpoint)*pix; 
x=repmat(xaxis,[length(xaxis),1]);
yaxis=(1-originpoint:N-originpoint)*pix;
y=repmat(yaxis',[1,length(xaxis)]);
% 角谱衍射坐标系
dfx=1/totalsize;
dx=(1-originpoint:N-originpoint)*dfx;
dfy=1/totalsize;
dy=(1-originpoint:N-originpoint)*dfy;
[dx,dy] = meshgrid(dx,dy); % 生成网格矩阵

%% 角谱衍射
beam1 = imag.*exp(-1j*k*(x.^2+y.^2)/(2*z1)); 
Beam1 = fftshift(fft2(fftshift(beam1))); 
h1 = exp(1j*k*z1*sqrt(1-(lambda*dx).^2-(lambda*dy).^2)); %空间传递函数
H1 = Beam1.* h1; 
A1 = fftshift(ifft2(fftshift(H1)));
B1 = (A1/max(max(A1))); %振幅归一化
E1 = abs(A1); %光强

%% 绘图
figure(1); %绘制观察屏距离10mm的图像
subplot(2,2,1),imshow(imag,[]),axis on ;title('孔径图像'); %绘制孔径图像
subplot(2,2,2),mesh(dx,dy,E1),title('光强三维分布'); %光强三维分布图
subplot(2,2,3),plot(dx(1,:),E1(N/2,:)),title('光强二维分布'); %光强二维分布图
pos=axis; %取得当前坐标轴的范围，即[xmin xmax ymin ymax]
xlabel('长度/mm','position',[1.2*pos(2) -0.01*pos(4)]); %设置x轴单位标签在图的右下方
ylabel('光强/cd','position',[1.09*pos(1) 0.9*pos(4)]); %设置y轴单位标签在图的左上方
subplot(2,2,4),imshow(abs(B1));axis image;colormap gray;title('衍射图像'); %绘制衍射图像，经测试后发现gray的线性映射成像效果最好

%% 不同距离对比
beam2 = imag.*exp(-1j*k*(x.^2+y.^2)/(2*z2));
Beam2 = fftshift(fft2(fftshift(beam2)));
H2 = Beam2.*exp(1j*k*z2*sqrt(1-(lambda*dx).^2-(lambda*dy).^2));
A2 = fftshift(ifft2(fftshift(H2)));
B2 = (A2/max(max(A2))); %振幅归一化
E2 = abs(A2); %光强
figure(2); %绘制观察屏距离15mm的图像
subplot(2,2,1),imshow(imag,[]),axis on ; %绘制孔径图像
subplot(2,2,2),mesh(x,y,E2),title('光强三维分布'); %光强三维分布图
subplot(2,2,3),plot(x(1,:),E2(N/2,:)),title('光强二维分布'); %光强二维分布图
pos=axis; %取得当前坐标轴的范围，即[xmin xmax ymin ymax]
xlabel('长度/mm','position',[1.2*pos(2) -0.01*pos(4)]); %设置x轴单位标签在图的右下方
ylabel('光强/cd','position',[1.09*pos(1) 0.93*pos(4)]); %设置y轴单位标签在图的左上方
subplot(2,2,4),imshow(abs(B2));axis image; colormap gray;title('衍射图像'); %绘制衍射图像，经测试后发现gray的线性映射成像效果最好

figure(3); % 作图对不同观察屏距离衍射的光强二维分布进行定量分析
subplot(1,2,1),plot(dx(1,:),E1(N/2,:)),title(['观察屏距离为',num2str(z1),'mm的光强二维分布']);pos=axis; %取得当前坐标轴的范围，即[xmin xmax ymin ymax]
xlabel('长度/mm','position',[1.2*pos(2) -0.01*pos(4)]); %设置x轴单位标签在图的右下方
subplot(1,2,2),plot(dx(1,:),E2(N/2,:)),title(['观察屏距离为',num2str(z2),'mm的光强二维分布']);pos=axis; %取得当前坐标轴的范围，即[xmin xmax ymin ymax]
xlabel('长度/mm','position',[1.2*pos(2) -0.01*pos(4)]); %设置x轴单位标签在图的右下方
