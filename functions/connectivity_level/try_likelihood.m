function [llh] = try_likelihood(mask,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param)
%% Perform ntry iterations
ntry          = param.ntry;
param.penalty = 2;
llh           = zeros(1,ntry);
for try_iter = 1:ntry
    %% Expectation
    [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_post] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,param);
    %%
    %% Maximization
    [sigma2xi,Sigmajj,Thetajj]            = higgs_maximization(Psixixi,Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
    Thetajj(mask)                         = 0;
    if(~param.use_gpu)
        [U,D]                 = eig(Thetajj);
        V                     = Iq/U;
        d                     = diag(D);
        if min(d) < 0
            d                 = d + abs(min(d));
        end
        dfix                  = d + max(d)*eigreg;
        [U,dfix,V]            = gather(U,dfix,V);
        clear D d
    else
        gpuThetajj            = gpuArray(Thetajj);
        [U,D]                 = eig(gpuThetajj);
        V                     = Iq/U;
        d                     = diag(D);
        if min(d) < 0
            d                 = d + abs(min(d));
        end
        dfix                  = d + max(d)*eigreg;
        [U,dfix,V]            = gather(U,dfix,V);
        clear D d
        clear gpuThetajj
    end
    Thetajj         = U*spdiags(dfix,0,q,q)*V;
    Thetajj         = (Thetajj + Thetajj')/2;
    Sigmajj         = V*spdiags(1./dfix,0,q,q)*U;
    Sigmajj         = (Sigmajj + Sigmajj')/2;
    %% Compute likelihood
    [llh(try_iter)] = higgs_likelihood(Sxixi,sigma2xi,Sjj,Thetajj,Sigmajj_post,param);
    llh0            = llh(try_iter);
    sigma2xi0       = sigma2xi;
    Sigmajj0        = Sigmajj;
end
end