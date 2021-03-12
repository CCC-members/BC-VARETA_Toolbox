function [sigma2xi0,Sigmajj0,llh0] = higgs_warm_start(Svv,Lvj,Thetajj,Sigmajj,param)
%% Warm-start by cross-validating diagonal initial values
% Grid space for Cross-validation
p               = param.p;
grid_sigma2xi   = (-5:0.25:5);
grid_sigma2j    = (-5:0.25:5);
%% Cross-validation cycle
llh = higgs_crossvalidation(grid_sigma2xi,grid_sigma2j,Sigmajj,Thetajj,Svv,Lvj,param); 
%% optimal value
[val]           = max(llh(:));
[id_xi,id_j]    = find(llh == val);
sigma2xi        = 10^(grid_sigma2xi(id_xi));
sigma2xi        = sigma2xi*ones(p,1);
theta2xi        = (1/sigma2xi)*ones(p,1);
sigma2j         = 10^(grid_sigma2j(id_j));
Sigmajj.X       = sigma2j*Sigmajj.X;
Sigmajj.d       = sigma2j*Sigmajj.d;
Thetajj.X       = (1/sigma2j)*Thetajj.X;
Thetajj.d       = (1/sigma2j)*Thetajj.d;
%% Recompute likelihood
% Expectation
[Sxixi,Psixixi,Sjj,Psijj,Sigmajj_pst] = higgs_expectation(Svv,Lvj,sigma2xi,Sigmajj,param);
% Compute likelihood
llh0            = higgs_likelihood(Sxixi,theta2xi,Sjj,Thetajj,Sigmajj_pst,param);
sigma2xi0       = sigma2xi;
Sigmajj0        = Sigmajj;
end