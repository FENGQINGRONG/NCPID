
close all;
clc
clear variables
addpath('.\data')
addpath('.\operater')

maxValuelist    = 500;
psf=  fspecial('average',7); 

%% To do data gray   
  
    Img = imread('.\data\baby.tif');

      
for mm = 1:length(maxValuelist)
    max_dim  = length(maxValuelist);
    MaxValue = maxValuelist(mm);
    
    %%
    params.psf          = psf; 
    K = BlurMatrix(psf, size(Img));
    params.K = K;
    params.H = K;
    Img = double(Img);
    Img = MaxValue * Img / max(Img(:));
    snr=@(x,Img)10*log10(norm(Img(:))^2 /norm(Img(:) - x(:))^2 );
    
    params.Img = Img;
    stream = RandStream('mt19937ar', 'Seed', 88);
    RandStream.setGlobalStream(stream);
    
    Blr = K * Img; 
    Blr = max(0, Blr);
    Bn = poissrnd(Blr);
%     figure2=imshow(Bn,[]);
 
    % get PSNR and SSIM values for degraded image
    j=1;
    psnr_noisy = psnr(Bn, Img, MaxValue);
    snr_noisy = snr(Bn, Img);
    
    ssim_noisy = ssim(Bn, Img, 'DynamicRange', MaxValue);
    ReErr_noisy =norm(Bn(:)-Img(:))/norm(Img(:));
    
%% NCPID('NON-AC')
params.alpha        =1;
params.lambda       = 1000;
params.mu           = [0.005,0.005,0.02];
params.MaxValue     = MaxValue;
params.Nit          = 50;
params.tol          = 1.0e-3;
params.frame       =1;
params.Level       =4;
    tic 
    out = poisson_L1L2(Bn,Bn,params); 
    toc

    u_L1L2=out.sol;
    psnr_L1L2 = psnr(out.sol, Img, MaxValue); 
    snr_L1L2 = snr(out.sol, Img);
    ssim_L1L2 = ssim(out.sol, Img, 'DynamicRange', MaxValue);
    ReErr_L1L2 =norm(out.sol(:)-Img(:))/norm(Img(:));

%% NCPID('AC')
params.alpha        = 1;
params.lambda       = 600;
params.tau          = 1.05;
params.mu           = [10,0.01,0.5];
params.MaxValue = MaxValue;

params.frame       =1;
params.Level       =4;
    tic 
    out = poisson_L1L2_fast(Bn,Bn,params); 
    toc


 u_L1L2_fast=out.sol;
    psnr_L1L2_fast = psnr(out.sol, Img, MaxValue); 
    snr_L1L2_fast = snr(out.sol, Img); 
    ssim_L1L2_fast = ssim(out.sol, Img, 'DynamicRange', MaxValue);
    ReErr_L1L2_fast =norm(out.sol(:)-Img(:))/norm(Img(:));

end
