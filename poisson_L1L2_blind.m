function [u,psf,out] = poisson_L1L2_blind(f,orig,ker,opts)

%%%.
% Overlapping Group Sparse and Nonconvex Second-order Total Variation Priors for Restoring Poisson Noisy Images
%

Nit             = opts.Nit;
tol             = opts.tol;
alpha           = opts.alpha;
lambda          = opts.lambda; % The regularization parameter.
gamma           = opts. gamma;
mu              = opts.mu;
tau             = opts.tau;
frame           = opts.frame;
Level           = opts.Level;
wLevel          = opts.wLevel;
relError        = zeros(Nit,1); % Compute error relative to the previous iteration.


[D,R]=GenerateFrameletFilter(frame);
W  = @(x) FraDecMultiLevel(x,D,Level); % Frame decomposition
WT = @(x) FraRecMultiLevel(x,R,Level); % Frame reconstruction
%% initizing image u and kernerl K

[row, col] = size(f);
Fk = psf2otf(ker,[row, col]);
%% initialization of the auxiliary and dual variables
% Lagrange multipliers
c=zeros(row, col);
b=zeros(row, col);
e=W(f);

% variable initialization
w=real(ifft2(fft2(f).*Fk));
z=zeros(row, col);
d=W(f);
u=zeros(size(f));
%% initialization of energy and break

%%

%% Deconvolution

iter=1;
tstart = tic;
for i=1:20
    beta=-u./(norm(u,2)+1e-10);
for  it = 1:Nit
    
    u_old   = u;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% Update PSF    
    
    C=CoeffOper1('-',d,e);%d-b
    lhs     = mu(1) * abs(Fk).^2 +mu(2)+mu(3);
    rhs     = mu(1) * conj(Fk).*fft2(w - c) +...
              mu(2) * fft2((z - b))+mu(3) * fft2(WT(C))-alpha*fft2(beta);% conj(otf6)-e
     Fu       = rhs./lhs;
     u       = real(ifft2(Fu));
    % w sub-problem
%     Ku = K * u;
     Ku = real(ifft2(Fu.*Fk));
    temp = mu(1) * (Ku + c) -lambda;
    w = (temp + sqrt(temp.^2 + 4 * mu(1) * lambda .* f))/ (2 * mu(1));
     
     %%************** d sub-problem***********
    muLevel = getwThresh(1/mu(3),wLevel,Level,D);
    Cpe=CoeffOper1('+',W(u),e);%d = prox_l1(W(u)+e, alpha/mu(3))
    d=CoeffOper1('vs',Cpe,muLevel);%     d = prox_l1(Cpe, alpha/mu(3));
    
     
         %***** z sub-problem **
    z = max(u+b,0);
%     z = min(z,MaxValue);
   

    %******* c b e Lagrange multiplier update ******************
    c = c + Ku - w;
    b = b + u - z;
    deltab=CoeffOper1('-',W(u),d);
    e=CoeffOper1('+',e,deltab);
    
   
  
    %***** Some statistics ***
    relError(it)    = norm(u - u_old, 'fro')/norm(u_old+1e-10, 'fro');
    relerr_1=relError(it);%
    relerr_2 =norm(u-orig,'fro')/norm(orig,'fro');
    relerr_3 = sqrt(sum((u-u_old).^2))/max(sqrt(sum(u_old.^2)),1);
    relerr_4 = log10(norm(u-u_old, 2)/(norm(u_old, 2)+1e-10));%效果好
    itererr(iter,1:4)=[relerr_1,relerr_2,relerr_3,relerr_4];
    iter=iter+1;
    
    if relError(it) < tol
        break;
    end
      mu=min(tau*mu,1e15);
end
   %% update psf   
       scale = conv2c(u,ker);
    
    %% Avoid division by zero
    index = find(~scale);
    scale(index) = eps;
    kk = f./scale;
    kk(index) = 0;
    kkk = real(otf2psf(conj(fft2(u)) .* fft2(kk), size(ker)));
    psf = ker.*kkk;
    
    %% Compute the divergence
    divergpsf = divergence(psf);
    divipsf = 1 - gamma/lambda * divergpsf;
    index = find(~divipsf);
    divipsf(index) = eps;
    
    psf1 = psf./divipsf;
    psf1(index) = 0;
    psf = kkk.*psf1;
    
    psf = max(psf,0);
    psf = psf / sum(psf(:));  % check kernel sums to 1
       
    Fk = psf2otf(psf,[row, col]);


    outrelerr_1= norm(u - u_old, 'fro')/norm(u_old+1e-10, 'fro');
    outrelerr_2 =norm(u-orig,'fro')/norm(orig,'fro');
    outrelerr_3 = sqrt(sum((u-u_old).^2))/max(sqrt(sum(u.^2)),1);
    outrelerr_4 = log10(norm(u-u_old, 2)/(norm(u_old, 2)+1e-10));
    outitererr(i,1:4)=[outrelerr_1,outrelerr_2,outrelerr_3,outrelerr_4];
end
out.sol                 = u;
out.relativeError       = itererr;
out.outrelativeError    = outitererr;
out.OverallItration     = iter-1; %No of itr to converge
end

%%
function [D,R]=GenerateFrameletFilter(frame)
if frame==0          %Haar Wavelet
    D{1}=[0 1 1]/2;
    D{2}=[0 1 -1]/2;
    D{3}='cc';
    R{1}=[1 1 0]/2;
    R{2}=[-1 1 0]/2;
    R{3}='cc';
elseif frame==1      %Piecewise Linear Framelet
    D{1}=[1 2 1]/4;
    D{2}=[1 0 -1]/4*sqrt(2);
    D{3}=[-1 2 -1]/4;
    D{4}='ccc';
    R{1}=[1 2 1]/4;
    R{2}=[-1 0 1]/4*sqrt(2);
    R{3}=[-1 2 -1]/4;
    R{4}='ccc';
elseif frame==3      %Piecewise Cubic Framelet
    D{1}=[1 4 6 4 1]/16;
    D{2}=[1 2 0 -2 -1]/8;
    D{3}=[-1 0 2 0 -1]/16*sqrt(6);
    D{4}=[-1 2 0 -2 1]/8;
    D{5}=[1 -4 6 -4 1]/16;
    D{6}='ccccc';
    R{1}=[1 4 6 4 1]/16;
    R{2}=[-1 -2 0 2 1]/8;
    R{3}=[-1 0 2 0 -1]/16*sqrt(6);
    R{4}=[1 -2 0 2 -1]/8;
    R{5}=[1 -4 6 -4 1]/16;
    R{6}='ccccc';
end
end
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
end

function c = m(a,b)
c = (sign(a) + sign(b)) ./ 2 .* min(abs(a),abs(b));
end
