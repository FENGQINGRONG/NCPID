function out = poisson_L1L2(U0,Img, opts)

%%%.
% Overlapping Group Sparse and Nonconvex Second-order Total Variation Priors for Restoring Poisson Noisy Images
%

 u               = U0;
[row, col] 		= size(Img);


Nit             = opts.Nit;
tol             = opts.tol;
alpha           = opts.alpha;
lambda          = opts.lambda; % The regularization parameter.
mu              = opts.mu;
% tau             = opts.tau;
MaxValue        = opts.MaxValue;
frame           = opts.frame;
Level           = opts.Level;
relError        = zeros(Nit,1); % Compute error relative to the previous iteration.

% Level=4;
% frame=1;
wLevel = 1/2;
[D,R]=GenerateFrameletFilter(frame);
W  = @(x) FraDecMultiLevel(x,D,Level); % Frame decomposition
WT = @(x) FraRecMultiLevel(x,R,Level); % Frame reconstruction
% Lagrange multipliers
c=zeros(row, col);
b=c;
e=W(u);

% variable initialization

z=zeros(row, col);
d=W(u);

psf = opts.psf;
K = opts.K;
w=K*u;
% w=zeros(row, col);
%**************Initialize Lagrange Multipliers***********************

eigK            = psf2otf(psf, [row col]); %In the fourier domain
eigKtK          = abs(eigK).^2;
II              = abs(psf2otf([1],[row col]).^2);




muLevel = getwThresh(1/mu(3),wLevel,Level,D);
%%
for i=1:20
  beta=-u./norm(u,2);
for k = 1:Nit
    
    %********* u sub-problem (use FFT's)********   
    u_old   = u;
    C=CoeffOper('-',d,e);%d-b
    lhs     = mu(1) * eigKtK +(mu(2)+mu(3))*II;
    rhs     = mu(1) * conj(eigK).*fft2(w - c) +...
              fft2(mu(2) * (z - b)+mu(3) * WT(C))-alpha*fft2(beta);% conj(otf6)-e
     u       = rhs./lhs;
     u       = real(ifft2(u));
    % w sub-problem
    Ku = K * u;
    temp = mu(1) * (Ku + c) -lambda;
    w = (temp + sqrt(temp.^2 + 4 * mu(1) * lambda .* Img))/ (2 * mu(1));
     
   
    
     
     %%************** d sub-problem***********
    Cpe=CoeffOper('+',W(u),e);%d = prox_l1(W(u)+e, alpha/mu(3))
    d=CoeffOper('vs',Cpe,muLevel);%     d = prox_l1(Cpe, alpha/mu(3));
    
     
         %***** z sub-problem **
    z = max(u+b,0);
    z = min(z,MaxValue); 

    %******* c b e Lagrange multiplier update ******************
    c = c + K * u - w;
    b = b + u - z;
    deltab=CoeffOper('-',W(u),d);
    e=CoeffOper('+',e,deltab);
%     e = e + W(u) - d;    
    %***** Some statistics ***
    relError(k)    = norm(u - u_old, 'fro')/norm(u_old, 'fro');
    
    if relError(k) < tol
        break;
    end
%   mu=min(tau*mu,1e15);

end
out.sol                 = u;
out.relativeError       = relError(1:k);
out.OverallItration     = size(out.relativeError,2); %No of itr to converge
end
end


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

