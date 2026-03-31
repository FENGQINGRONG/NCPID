function [u,w3,psf,out] = PIDSB_RLBlindDeconPro(g,ker,params)
%% If you find the code is useful to you, please cite our paper: 
% (1) Luxin Yan, Houzhang Fang, and Sheng Zhong,
% Blind image deconvolution with spatially adaptive total variation regularization
% Optics Letters, Vol. 37, No. 14, 2012.
% (2) Houzhang Fang, Luxin Yan, Hai Liu, and Yi Chang,
% Blind Poissonian images deconvolution with framelet regularization
% Optics Letters, Vol. 38, No. 4, 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Houzhang Fang, houzhangfang@xidian.edu.cn
% Copyright (C) 2015 Houzhang Fang. All rights reserved.
%----------------------------------------------------------------------

tol = params.tol;
mu = params.mu;
lambda = params.lambda;
gamma1 = params.gamma1;
frame = params.frame;
Level = params.Level;
wLevel = params.wLevel;
maxiter = params.maxit;

[m,n] = size(g);
Fh = psf2otf(ker,[m,n]);

%% Decomposition (D) and reconstruction (R) framelet operator
[D,R] = GenerateFrameletFilter(frame);
nD = length(D);
params.nD = nD;

%% Compute the weighted thresholding parameters.
muLevel = getwThresh(mu,wLevel,Level,D);
params.muLevel = muLevel;

%% Initialization
w1 = real(ifft2(fft2(g).*Fh));
w2 = FraDecMultiLevel(g,D,Level); % compute Dg
w3 = zeros(m,n);
b1 = zeros(m,n);
b2 = w2 ;
b3 = zeros(m,n);
u = g;

%% Calcaulate the energy functional value
out.NSDE = [];
out.pt = [];
out.YY = [];

if isvar('oimg')
    out.ppsnr = [];
    out.nmse = [];
end

loop = 0;
StopFlag = 0;

%% main loop
while StopFlag == 0
    loop = loop + 1;
    %% Save old values
    uold = u;
    
    %% Update PSF
    scale = conv2c(u,ker);
    
    %% Avoid division by zero
    index = find(~scale);
    scale(index) = eps;
    hh = g./scale;
    hh(index) = 0;
    hhh = real(otf2psf(conj(fft2(u)) .* fft2(hh), size(ker)));
    psf = ker.*hhh;
    
    %% Compute the divergence
    divergpsf = divergence(psf);
    divipsf = 1 - gamma1 * divergpsf;
    index = find(~divipsf);
    divipsf(index) = eps;
    
    psf1 = psf./divipsf;
    psf1(index) = 0;
    psf = hhh.*psf1;
    
    psf = max(psf,0);
    psf = psf / sum(psf(:));  % check kernel sums to 1
       
    Fh = psf2otf(psf,[m,n]);
    
    %% Update the image 
    [u,w1,w2,w3,b1,b2,b3] = oneloopiteration_f_frame(Fh,g,w1,w2,w3,b1,b2,b3,D,R,params);
    
    %% Calculate the relative error NSDE
    NSDE = norm(w3(:) - uold(:))/norm(w3(:));
    out.NSDE = [out.NSDE NSDE];
    
     if isvar('oimg')
%         ppsnr = psnr(w3,oimg);
%         out.ppsnr = [out.ppsnr ppsnr];
%          nmse_n = nmse(oimg,w3);
%          out.nmse = [out.nmse; nmse_n];
    end
    YY = w3;
    
    %% Iteration stop criterion
    disp(['loop: ' num2str(loop+1) ', NSDE: ' num2str(NSDE)]);
    if NSDE < tol
        StopFlag = 1;
        fprintf('Iteration convergence at loop = %d\n',loop);
    elseif loop >= maxiter
        StopFlag = 2;
        fprintf('Maximal iteration numbers is reached!!!\n');
    end
end

%% ------------------SUBFUNCTION----------------------------- 
function [u,w1,w2,w3,b1,b2,b3] = oneloopiteration_f_frame(Fh,g,w1,w2,w3,b1,b2,b3,D,R,params)
% Author: Houzhang Fang, houzhangfang@xidian.edu.cn
% Copyright (C) 2015 Houzhang Fang. All rights reserved.
%-------------------------------------------------------------------------
 % Calculate the image 
lambda = params.lambda;
gamma = params.gamma;
Level = params.Level;
muLevel = params.muLevel;
nD = params.nD;

%%% Solve u.
Ft1 = conj(Fh).*fft2(w1 - b1);

for ki=1:Level
    for ji=1:nD-1
        for jj=1:nD-1
            C{ki}{ji,jj} = w2{ki}{ji,jj} - b2{ki}{ji,jj};
        end
    end
end
Ft2 = fft2(FraRecMultiLevel(C,R,Level));
Ft3 = fft2(w3 - b3);
nume = Ft1 + Ft2 + Ft3;
deno = 2 + abs(Fh).^2;
Fu = nume./deno;
u = real(ifft2(Fu));

%%% Solve w1
Hu = real(ifft2(Fu.*Fh));
w1 = 0.5*( b1 + Hu - gamma + sqrt( ( b1 + Hu - gamma).^2 + 4*gamma*g ) );

%%% Solve w3
C = FraDecMultiLevel(u,D,Level); % compute Du^{k+1}
for ki=1:Level
    for ji=1:nD-1
        for jj=1:nD-1
            w2{ki}{ji,jj} = wthresh(C{ki}{ji,jj} + b2{ki}{ji,jj},'s',muLevel{ki}{ji,jj}*lambda*gamma);
        end
    end
end

%%% Solve w3
w3 = max(b3 + u, 0);

%%% Update b1
b1 = b1 + Hu - w1;

%%% Update b2
for ki=1:Level
    for ji=1:nD-1
        for jj=1:nD-1
            if ((ji~=1)||(jj~=1))||(ki==Level)
                deltab = C{ki}{ji,jj} - w2{ki}{ji,jj};
                b2{ki}{ji,jj} = b2{ki}{ji,jj} + deltab;
            end
        end
    end
end

%%% Update b3
b3 = b3 + u - w3;
return

%% ------------------SUBFUNCTION----------------------------- 
function diverg = divergence(x,hs)
% Author: Houzhang Fang, houzhangfang@xidian.edu.cn
% Copyright (C) 2015 Houzhang Fang. All rights reserved.
%----------------------------------------------------------------------
if ~exist('hs','var')
    hs = 1;
end
dxm = 1/hs * conv2(x,[-1 1 0],'same');
dxp = 1/hs * conv2(x,[0 -1 1],'same');
dym = 1/hs * conv2(x,[-1 1 0]','same');
dyp = 1/hs * conv2(x,[0 -1 1]','same');
den1 = sqrt(dxp.^2 + m(dym,dyp).^2);
index = find(~den1);
den1(index) = eps;
div1 = dxp ./ den1;
div1(index) = 0;

den2 = sqrt(dyp.^2 + m(dxm,dxp).^2);
index = find(~den2);
den2(index) = eps;
div2 = dyp ./ den2;
div2(index) = 0;

term1 = conv2(div1,[-1 1 0],'same');
term2 = conv2(div2,[-1 1 0]','same');
diverg = 1/hs .* term1 + 1/hs .* term2;

function c = m(a,b)
c = (sign(a) + sign(b)) ./ 2 .* min(abs(a),abs(b));