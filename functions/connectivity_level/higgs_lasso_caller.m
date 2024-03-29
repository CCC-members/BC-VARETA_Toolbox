function [Thetajj,Sigmajj,llh_inner] = higgs_lasso_caller(Psijj,Svv,Lvj,sigma2xi0,Sigmajj0,llh0,param)
%% Compute hg_lasso_lqa
%  it first corrects the eigenvalues of the effective empirical covariance Psijj and normalize it by the infimun
%  norm Psijjfix, the scaling guarantees more stable computations of the hg-lasso-lqa, continues with performing
%  hg-lasso-lqa Thetajj_lasso and the unbiased statistics Thetajj_unb and Thetajj_var
m               = param.m;
q               = param.q;
aj              = param.aj;
Ajj             = param.Ajj;
Psijj_fix       = Psijj;
Psijj_fix.X     = Psijj.X/sqrt(sum(Psijj.d.^2)/q);
Psijj_fix.d     = Psijj.d/sqrt(sum(Psijj.d.^2)/q);
param.a         = aj;
param.A         = Ajj;
param.maxiter   = param.maxiter_inner;
[Thetajj_lasso] = hg_lasso_lqa(Psijj_fix,param);
Thetajj_unb     = 2*Thetajj_lasso.X - Thetajj_lasso.X*Psijj_fix.X*Thetajj_lasso.X;
Thetajj_unb     = abs(Thetajj_unb);
Thetajj_unb     = (Thetajj_unb + Thetajj_unb')/2;
Thetajj_var     = sqrt(abs(diag(Thetajj_lasso.X))*abs(diag(Thetajj_lasso.X))' + abs(Thetajj_lasso.X).^2);
Thetajj_var     = Thetajj_var - diag(diag(Thetajj_var));
Thetajj_var     = (Thetajj_var + Thetajj_var')/2;
clear Thetajj_lasso Psijj_fix
%% Precompute hg_lasso_lqa equivalent to Frobenious solution
%  this is the closest solution to the hermitian graphical lasso with parameter aj and equivalent to the first
%  iteration of the lqa with gamma == 1, if for a given iteration the lasso model biases the computations and
%  the global hyperparameters posterior likelihood does not increases in ntry iterations higgs-lasso still goes
%  on with higgs-ridge prewarming, this approximation is a result of the empirical study of higgs-ridge with
%  parameter aj^2 (aj same as lasso aj = sqrt(log(q)/m)) in terms of stability and similarity to the unbiased
%  statistics given by Jankova conditions,
Thetajj.U       = Psijj.U;
Thetajj.d       = (sqrt(Psijj.d.^2 + 4*aj^2) - Psijj.d)/(2*aj^2);
Thetajj.X       = Thetajj.U*spdiags(Thetajj.d,0,q,q)*Thetajj.U';
%% try likelihood evolution for the mask given by all thresholds
% evaluate performance of Rayleigh thresholds in terms of the global hyperparameters posterior likelihood, the 
% overlapping of the distributions for null hypothesis values (Rayleigh) and the alternative introduces ambiguity 
% in the decision for a given Rayleigh threshold we move rth between the value corresponding to the peak of the 
% Rayleigh distribution (rth1)
[rth]           = higgs_lasso_stabilizer(Psijj,Thetajj,Thetajj_unb,Thetajj_var,Svv,Lvj,sigma2xi0,Sigmajj0,llh0,param);
Thetajj.X(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var)))) = 0;
[Thetajj]       = higgs_eigendecomposition(Thetajj.X,param);
Sigmajj.U       = Thetajj.U;
Sigmajj.d       = 1./Thetajj.d;
Sigmajj.X       = Sigmajj.U*spdiags(Sigmajj.d,0,q,q)*Sigmajj.U';
end
