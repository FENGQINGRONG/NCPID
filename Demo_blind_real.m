close all;
clc
clear variables
path(path,'operater1')
path(path,'frameletlib');

[ filename, pathname ] = uigetfile('./images/*.*', 'load image');
x = imread( fullfile( pathname, filename ) );

    if ndims( x ) == 3
    x = im2double( rgb2gray( x ) );
     else
    x = im2double( x );
    end
     g=x;

%% ����ʼģ���
siz0 = 15;  sigma0 = 1.8;
psf0 = fspecial('gaussian',siz0,sigma0);
[Mh,Nh]=size(psf0);

params.lambda = 0.001;
params.gamma = 10/(params.lambda);
params.gamma1 = 0.0008;
params.maxit = 180;

%The following parameters can be fixed
params.mu = 0.06;
params.frame = 1;
params.Level = 4;
params.wLevel = 1/2;
params.tol = 1e-4;

[upro,w3pro,psfpro,out] = PIDSB_RLBlindDeconPro(g,psf0,params);
%%
params.alpha        = 1;
params.lambda       = 200;%200
params.gamma        =1;%0.006*params.lambda;
params.tau          = 1.05;
a=0.00002;
params.mu           = [a,a,a];
params.Nit          = 50;
params.tol          = 1.0e-3;
params.frame = 1;
params.Level = 4;
params.wLevel = 1/2;

[u_L1L2,psfpro] = poisson_L1L2_blind(upro,g,psf0,params);



