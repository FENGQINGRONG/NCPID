
close all;
clc
clear variables
cd('C:\Users\Administrator\Desktop\3.7ÊµÑé_nonblind')
addpath('.\operater')
%addpath('.\compare_method')
maxValuelist    = [500];%[350 300 200 100];
% alphaValueList    =0.1];%[25 20 15 10];
% H   =  fspecial('gaussian',11, 3); 
% H   =  fspecial('average',11); %11x11µÄ¾ØÕó´óÐ¡Îª1/121=0.0083
psf             = fspecial('gaussian', [25,25],1.6); 
% figure1=imshow(psf,[]);
cur_file = ['baby.tif', 'boat.pgm', 'bridge.pgm', 'buty(256).png', 'einstein.pgm', 'house.pgm', 'lena.pgm', 'leopard.pgm', 'parrot.png', 'peppers.pgm', 'tulips.pgm'];

 datasave_time=zeros(10,2);
%% To do data gray
folder = 'C:\Users\Administrator\Desktop\3.7ÊµÑé_nonblind\Test image_nonblind\';
for picnum = 1:3
        if picnum==1    
        filepaths = dir(fullfile(folder,'*.png'));
        end
        if picnum==2    
        filepaths = dir(fullfile(folder,'*.tif'));
        end
        if picnum==3    
        filepaths = dir(fullfile(folder,'*.bmp'));
        end
%         if picnum==4    
%         filepaths = dir(fullfile(folder,'*.pgm'));
%         end
 
for ii = 1 : length(filepaths)

    filename = filepaths(ii).name;
  
    Img = imread(fullfile(folder, filename));
% display(sprintf('denoise processing of %s...', cur_file));
% for ii = 1:10  
        if picnum==2    
        ii=ii+7;
        end
        if picnum==3    
        ii=ii+8;
        end
%         if picnum==4    
%         ii=ii+10;
%         end
for mm = 1:length(maxValuelist)
    max_dim  = length(maxValuelist);
    MaxValue = maxValuelist(mm);
    
    %%
%     params = ParamSet(MaxValue);
    params.psf          = psf; % For denoising. (K is the identity operator)

%  params.alpha = alphaValueList(mm);  % regularization parameters
    
%     Img = imread(strcat(ima_dir, filesep, cur_file(ii))); %gray-scale image-
%     figure, imshow(Img, []), title('observed image')
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




%% L1L2
params.alpha        = 1;
params.lambda       = 600;
params.mu           = [1e-3,1e-3,0.05];
params.MaxValue     = MaxValue;
params.Nit          = 50;
params.tol          = 1.0e-3;
params.frame       =1;
params.Level       =4;


    t0 = cputime;
    out = poisson_L1L2(Bn,Bn,params); 
    t =  cputime - t0;
    datasave_time(ii,mm)  = t;
    
    mm=mm+1;
    
    
params.alpha        = 1;
params.lambda       = 2000;
params.tau           =1.07;
params.mu           = [10,0.05,0.1];    
    t1 = cputime;
    out = poisson_L1L2_fast(Bn,Bn,params); 
    tt =  cputime - t1;
    datasave_time(ii,mm)  = tt; 
    

end
end
end
 save('.\Result\nonblind_gauss_imgvscpu.mat')
