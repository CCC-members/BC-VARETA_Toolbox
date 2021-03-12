function [Thetajj,Sjj,gamma_grid,gamma,gcv,Tjv] = eloreta_hg_lasso(Svv,Lvj,param)
p             = size(Lvj,1);
Ip            = eye(p);
scaleLvj      = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
Lvj           = Lvj/scaleLvj;
scaleV        = (sum(abs(diag(Svv)))/p);
Svv           = Svv/scaleV;
gamma1        = param.gamma1;
gamma2        = param.gamma2;
delta_gamma   = param.delta_gamma;
gamma_grid    = gamma1:delta_gamma:gamma2;
gcv           = zeros(length(gamma_grid),1);
count         = 1;
%%

for gamma = gamma_grid
   [Tjv,Wout]       = mkfilt_eloreta(Lvj,gamma);
    Tjv              = Tjv';
    Txiv             = Ip - Lvj*Tjv;
    gcv(count)       = (1/p)*sum(abs(diag(Txiv*Svv*Txiv')))/((1/p)*sum(abs(diag(Txiv))))^2;
    count             = count + 1;
end
%%
[gcv_opt,idx_gamma]       = min(gcv);
gamma                     = gamma_grid(idx_gamma);
[Tjv,Wout]                = mkfilt_eloreta(Lvj,gamma);
Tjv                       = Tjv';
Sjj                       = higgs_eigendecomposition(Tjv*Svv*Tjv',param);
[Thetajj,Sigmajj]         = twostep_lasso_caller(Sjj,param);
Thetajj                   = Thetajj.X;
Sjj                       = Sjj.X;
end