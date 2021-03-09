function [sigma2xi,Sigmajj,Thetajj] = higgs_maximization(Psixixi,Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param)
p       = param.p;
q       = param.q;
axi     = param.axi;
aj      = param.aj;
penalty = param.penalty;
%% Source precision/covariance matrix estimator
if penalty == 0 % naive
    if(~param.use_gpu)
        [U,D]             = eig(Psijj);
        V                 = Iq/U;
        d                 = diag(D);
        if min(d) < 0
            d             = d + abs(min(d));
        end
        dfix              = d + max(d)*eigreg;
        [U,dfix,V]        = gather(U,dfix,V);
        clear D d
    else
        gpuPsijj          = gpuArray(Psijj);
        [U,D]             = eig(gpuPsijj);
        V                 = Iq/U;
        d                 = diag(D);
        if min(d) < 0
            d             = d + abs(min(d));
        end
        dfix              = d + max(d)*eigreg;
        [U,dfix,V]        = gather(U,dfix,V);
        clear D d
        clear gpuPsijj
    end
    Thetajj               = V*spdiags(1./dfix,0,q,q)*U;
    Thetajj               = (Thetajj + Thetajj')/2;
    Sigmajj               = U*spdiags(dfix,0,q,q)*V;
    Sigmajj               = (Sigmajj + Sigmajj')/2;
elseif penalty == 1 % lasso
    [Thetajj,Sigmajj]     = higgs_lasso_caller(Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
elseif penalty == 2 % ridge
    if(~param.use_gpu)
        [U,D]             = eig(Psijj);
        V                 = Iq/U;
        d                 = diag(D);
        if min(d) < 0
            d             = d + abs(min(d));
        end
        dfix              = d + max(d)*eigreg;
        [U,dfix,V]        = gather(U,dfix,V);
        clear D d
    else
        gpuPsijj          = gpuArray(Psijj);
        [U,D]             = eig(gpuPsijj);
        V                 = Iq/U;
        d                 = diag(D);
        if min(d) < 0
            d             = d + abs(min(d));
        end
        dfix              = d + max(d)*eigreg;
        [U,dfix,V]        = gather(U,dfix,V);
        clear D d
        clear gpuPsijj
    end
    d_frob                = (sqrt(dfix.^2 + 4*aj^2) - dfix)/(2*aj^2);
    Thetajj               = U*spdiags(d_frob,0,q,q)*V;
    Thetajj               = (Thetajj + Thetajj')/2;
    Sigmajj               = V*spdiags(1./d_frob,0,q,q)*U;
    Sigmajj               = (Sigmajj + Sigmajj')/2;
end
%%
%% Residual precision/variance estimator
% sigma2xi          = (sum(abs(diag(Psixixi)))/p + axi)*ones(p,1);
sigma2xi          = abs(diag(Psixixi)) + axi;
end