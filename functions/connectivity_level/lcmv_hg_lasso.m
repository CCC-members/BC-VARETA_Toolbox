function [Thetajj,Sjj,Tjv] = lcmv_hg_lasso(Svv,Lvj,param)
p             = size(Lvj,1);
scaleLvj      = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
Lvj           = Lvj/scaleLvj;
scaleV        = (sum(abs(diag(Svv)))/p);
Svv           = Svv/scaleV;
gamma         = param.gamma;
%%
[Tjv,T1jv,Wout]   = mkfilt_lcmv(Lvj,Svv,gamma);
Tjv               = Tjv';
%%
Sjj               = higgs_eigendecomposition(Tjv*Svv*Tjv',param);
[Thetajj,Sigmajj] = twostep_lasso_caller(Sjj,param);
Thetajj           = Thetajj.X;
Sjj               = Sjj.X;
end