function [rth] = higgs_lasso_stabilizer(Thetajj_unb,Thetajj_var,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param)
rth1      = param.rth1;
rth2      = param.rth2;
m         = param.m;
ntry      = param.ntry;
%% try likelihood evolution for the mask given by all thresholds
rth_grid  = rth1:0.1:rth2;
llh_grid  = zeros(length(rth_grid),ntry);
for rth_count = 1:length(rth_grid)
    %% Zeros mask
    rth                = rth_grid(rth_count);
    mask               = find(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var))));
    %% Perform ntry iterations
    [llh_grid(rth_count,:)] = try_likelihood(mask,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
end
[rth]     = check_likelihood_trend(llh0,llh_grid,rth_grid,ntry);
end