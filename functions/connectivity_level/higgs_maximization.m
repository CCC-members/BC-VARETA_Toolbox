function [sigma2xi,theta2xi,Sigmajj,Thetajj] = higgs_maximization(Psixixi,Psijj,Svv,Lvj,sigma2xi0,Sigmajj0,llh0,param)
q       = param.q;
axi     = param.axi;
aj      = param.aj;
penalty = param.penalty;
%% Source precision/covariance matrix estimator
if penalty == 0 % naive
    Thetajj.U             = Psijj.U;
    Thetajj.d             = 1./Psijj.d;
    Thetajj.X             = Thetajj.U*spdiags(Thetajj.d,0,q,q)*Thetajj.U';
    Sigmajj               = Psijj;
elseif penalty == 1 % lasso
    [Thetajj,Sigmajj]     = higgs_lasso_caller(Psijj,Svv,Lvj,sigma2xi0,Sigmajj0,llh0,param);
elseif penalty == 2 % ridge
    Thetajj.U             = Psijj.U;
    Thetajj.d             = (sqrt(Psijj.d.^2 + 4*aj^2) - Psijj.d)/(2*aj^2);
    Thetajj.X             = Thetajj.U*spdiags(Thetajj.d,0,q,q)*Thetajj.U';
    Sigmajj.U             = Thetajj.U;
    Sigmajj.d             = 1./Thetajj.d;
    Sigmajj.X             = Sigmajj.U*spdiags(Sigmajj.d,0,q,q)*Sigmajj.U';
end
%%
%% Residual precision/variance estimator
% sigma2xi          = (sum(abs(diag(Psixixi)))/p + axi)*ones(p,1);
sigma2xi          = abs(diag(Psixixi)) + axi;
theta2xi          = 1./sigma2xi;
end