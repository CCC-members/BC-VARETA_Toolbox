function [llh] = higgs_try_likelihood(mask,Svv,Lvj,sigma2xi0,Sigmajj0,llh0,param)
%% Perform ntry iterations
q             = param.q;
ntry          = param.ntry;
param.penalty = 2;
llh           = zeros(1,ntry);
for try_iter = 1:ntry
    %% Expectation
    [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_post] = higgs_expectation(Svv,Lvj,sigma2xi0,Sigmajj0,param);
    %%
    %% Maximization
    [sigma2xi,theta2xi,Sigmajj,Thetajj]  = higgs_maximization(Psixixi,Psijj,Svv,Lvj,sigma2xi0,Sigmajj0,llh0,param);
    Thetajj.X(mask)             = 0;
    [Thetajj]                   = higgs_eigendecomposition(Thetajj.X,param);
    Sigmajj.U                   = Thetajj.U;
    Sigmajj.d                   = 1./Thetajj.d;
    Sigmajj.X                   = Sigmajj.U*spdiags(Sigmajj.d,0,q,q)*Sigmajj.U';
    %% Compute likelihood
    [llh(try_iter)]             = higgs_likelihood(Sxixi,theta2xi,Sjj,Thetajj,Sigmajj_post,param);
    llh0                        = llh(try_iter);
    sigma2xi0                   = sigma2xi;
    Sigmajj0                    = Sigmajj;
end
end