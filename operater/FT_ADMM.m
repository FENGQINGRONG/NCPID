function u=FT_ADMM(f,ker)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   lambda*||Wu||_0+1/2||Au-f||_2^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nit             = opts.Nit;
tol             = opts.tol;
gamma           = opts.gamma;
mu              = opts.mu;
eta             = opts.eta;
lambda          = opts.lambda; % The regularization parameter.
rho             = opts.rho;
maxit       = opts.MaxValue;
frame           = opts.frame;
Level           = opts.Level;
relError        = zeros(Nit,1); % Compute error relative to the previous iteration.

[D,R]=GenerateFrameletFilter(frame);nD=length(D);
W  = @(x) FraDecMultiLevel(x,D,Level); % Frame decomposition
WT = @(x) FraRecMultiLevel(x,R,Level); % Frame reconstruction
[m,n]=size(f);
cker=rot90(ker,2);

blur=@(f,k)imfilter(f,k,'circular');
A = @(x) blur(x,ker); AT = @(x) blur(x,cker); % Convolution operator

eigenP=eigenofP(ker,(eta+rho)/gamma,m,n);

FT = @(x)fft2(x);
IFT = @(x)ifft2(x);

%initialization
w=zeros(m,n);
u=w;
z=w;
% b=W(f);
% d=b;

p1=zeros(m,n);
p2=zerso(m,n);

for nstep=1:maxit
    % solve w
    Au = A(u);
    temp = gamma * Au +p1 -mu;
    w = (temp + sqrt(temp.^2 + 4 * gamma * mu .* Img))/ (2 * gamma);
    
    ksi1=w-p1/gamma;
    ksi2=z-p2/eta;
 % solve u
     % solve half 
     
    %     d=CoeffOper('l0',Cpb,muLevel);
        C=DALHT(W(u),d,lambda/gamma,rho/gamma,0);  
        u=IFT(FT(AT(ksi1)+eta/gamma*ksi2+rho/gamma*WT(C))./eigenP);
    % solve z
    z=max(u+p2/eta,0);

   
    
    % update p1 p2
    p1 = p1 + gamma*(A(u) - w);
    p2 = p2 + eta*(u - z);
    
     relError(k)    = norm(u - u_old)/norm(u_old);
    
    if relError(k) < tol
        break;
    end

end


function eigenP=eigenofP(ker,mu,m,n)

[nker,mker]=size(ker);
tmp=zeros(m,n);tmp(1:nker,1:mker)=ker;
tmp=circshift(tmp,[-floor(nker/2),-floor(mker/2)]);
eigenP=abs(fft2(tmp)).^2+mu;

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

function muLevel=getwThresh(mu,wLevel,Level,D)
nfilter=1;
nD=length(D);
if wLevel<=0
    for ki=1:Level
        for ji=1:nD-1
            for jj=1:nD-1
                muLevel{ki}{ji,jj}=mu*nfilter*norm(D{ji})*norm(D{jj});
            end
        end
        nfilter=nfilter*norm(D{1});
    end
else
    for ki=1:Level
        for ji=1:nD-1
            for jj=1:nD-1
                if ji==1 && jj==1
                    muLevel{ki}{ji,jj}=0;
                else
                    muLevel{ki}{ji,jj}=mu*nfilter;
                end
            end
        end
        %         muLevel{ki}{1,2}=0;muLevel{ki}{2,1}=0;
        %         muLevel{ki}{2,2}=0;
        %muLevel{ki}{2,3}=0;muLevel{ki}{3,2}=0;
        %muLevel{ki}{3,3}=0;
        nfilter=nfilter*wLevel;
    end
end