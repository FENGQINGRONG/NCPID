close all;
clc
clear variables

path(path,'operater1')
path(path,'frameletlib');


%% To do data gray

x = imread('.\images\satellite.jpg');
    if ndims( x ) == 3
    lam = im2double( rgb2gray( x ) );
     else
    lam = im2double( x );
    end
    
img = lam;
[row,col] = size(img);
 scale = 20; 
orig = scale*img;
snr=@(x,I)10*log10(norm(I(:))^2 /norm(I(:) - x(:))^2 );
%% ����˻�ͼ��ʱ����ģ��
siz = 15;  sigma = 1.8;
psf = fspecial('gaussian',siz,sigma);
z = conv2c(orig,psf);
%% Add poisson noise
%% fixes pseudo-random noise
randn('seed',1); rand('seed', 1); 
randn('state',1); rand('state', 1);

peak=100;
g = poissrnd(peak*z)/peak;
g(isnan(g)) = 0;
%% ����ʼģ���
siz0 = 15;  sigma0 =1.8;
 psf0 = fspecial('gaussian',siz0,sigma0);
params.alpha        = 1;
params.lambda       = 400;
params.gamma        =1;
params.tau          = 1.05;
params.mu           = [0.2,0.2,0.2];
params.Nit          = 50;
params.tol          = 1.0e-5;
params.frame = 1;
params.Level = 4;
params.wLevel = 1/2;

upro=zeros(size(g));
[u_L1L2,psfpro,out] = poisson_L1L2_blind(g,orig,psf0,params);
SNR_L1L2=snr(u_L1L2, orig);




