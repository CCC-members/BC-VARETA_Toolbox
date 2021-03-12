function [eigenX] = higgs_eigendecomposition(X,param)
use_gpu    = param.use_gpu;
eigreg     = param.eigreg;

%% perform singular value decomposition
X          = (X + X')/2;
if(use_gpu) 
    gpuX   = gpuArray(X);
    [U,D]  = eig(gpuX);
    d      = diag(D);
    clear gpuX X D
    [U,d]  = gather(U,d);
else
    [U,D]  = eig(X);
    d      = diag(D);
    clear X D
end

%% regularize singular values
dmax       = max(d);
dmin       = min(d);
if (dmin == 0 || dmin/dmax < eigreg)
    d      = (d + eigreg*dmax)/(1 + eigreg);
end

%% perform singular value composition
eigenX.U   = U;
eigenX.d   = d;
eigenX.X   = U*spdiags(d,0,length(d),length(d))*U';
end