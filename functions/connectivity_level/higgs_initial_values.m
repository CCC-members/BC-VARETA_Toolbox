function [Svv,Lvj,scale,scaleLvj,sigma2xi0,Sigmajj0,llh0,param] = higgs_initial_values(Svv,Lvj,param)
%% Parameters
prew            = param.prew;
m               = param.m;
p               = param.p;
Op              = param.Op;
q               = param.q;
Iq              = param.Iq;
aj              = param.aj;
penalty         = param.penalty;
%%
%% Lead Field and Data scaling by frobenious norm equivalent to the largest singular values or the
% average variance of the observations, helps to treat the problen in a standard scale for every case,
% a projected identity matrix Svv = Lvj*Iq*Ljv under this conditions will generate a signal with
% (sum(abs(diag(Svv)))/p) = 1 if we fix the value of param.axi = 1E-2 the inferior admisible noise
% will be 1% of the average signal amplitude
scaleLvj        = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
Lvj             = Lvj/scaleLvj;
scaleV          = (sum(abs(diag(Svv)))/p);
Svv             = Svv/scaleV;
%%
%% Decide not-prewarming (prew = 0) or prewarming (prew = 1), prew = 1 is recomended for severe
%  ill-conditioning, whereas the number of observed variables is of the same order of magnitude
%  that the number of hidden variables use prew = 0, in both cases the scale of Svv is readjusted
%  by and initial computation of Sjj making the maximum singular values of Sjj and Sigmajj_pst
%  identique,this will avoid biasing the EM initial step
if prew == 0 % not-prewarming
    %  - if prew = 0 Sjj is precomputed by the first EM Expectation step with unitary initialization
    %  Sigmajj0 = Iq sigma2xi0 = 1, it rescales Svv and axi by scaleJ the structure of Thetajj and Sigmajj
    %  is the identity matrix    
    Sigmajj0.X      = Iq;
    sigma2xi0       = Op;
    [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_post] = higgs_expectation(Svv,Lvj,sigma2xi0,Sigmajj0,param);
    scaleJ          = mean(abs(diag(Psijj.X)))/mean(abs(diag(Sigmajj_post.X)));
    Svv             = Svv/scaleJ;
    param.axi       = param.axi/scaleJ;
    scale           = scaleLvj^2/(scaleV*scaleJ);
    Thetajj.d       = diag(Iq);
    Thetajj.X       = Iq;
    Sigmajj.d       = diag(Iq);
    Sigmajj.X       = Iq;
elseif prew == 1 % prewarming
    %  - if prew = 1 we use an initial subspace two-step solution based on the crossspectral enet-ssbl
    %  (equivalent to an univariate higgs) and hermitian graphical model, this provides a multivariate
    %  initialization to the EM or prewarming, crossspectral enet-ssbl also adjusts the scale by scaleJ
    %  for every penalty a two-step graphical computation is performed to produce the initial structure
    %  of Thetajj and Sigmajj    
    disp("-->> Running prewarming, source connectivity. This operation may take a long time."); 
    parcellation        = [];
    counter             = 1;
    for ii = 1:1:q
        parcellation{counter} = ii;
        counter               = counter + 1;
    end
    param.Nsamp         = m;
    param.parcellation  = parcellation;
    param.W             = eye(size(Lvj,2));
    param.flag          = "-->> Running prewarming, source activity";
    [s2j,sigma2j_post,Tjv,Svv,scaleJ,~] = higgs_sSSBLpp(Svv,Lvj,param);
%     [Psijj]               = higgs_eigendecomposition(Tjv*Svv*Tjv' + W*sigma2jW - sigma2j_post0*LvjWsigma2jW,param);
    [Psijj]              = higgs_eigendecomposition(Tjv*Svv*Tjv',param);
    if penalty == 0 % naive       
        Thetajj.U       = Psijj.U;
        Thetajj.d       = 1./Psijj.d;
        Thetajj.X       = Thetajj.U*spdiags(Thetajj.d,0,q,q)*Thetajj.U';
        Sigmajj         = Psijj;
        clear Sjj
    elseif penalty == 1 % lasso
        [Thetajj,Sigmajj] = twostep_lasso_caller(Psijj,param);
        clear Sjj
    elseif penalty == 2 % ridge
        Thetajj.U       = Psijj.U;
        Thetajj.d       = (sqrt(Psijj.d.^2 + 4*aj^2) - Psijj.d)/(2*aj^2);
        Thetajj.X       = Thetajj.U*spdiags(Thetajj.d,0,q,q)*Thetajj.U';
        Sigmajj.U       = Thetajj.U;
        Sigmajj.d       = 1./Thetajj.d;
        Sigmajj.X       = Sigmajj.U*spdiags(Sigmajj.d,0,q,q)*Sigmajj.U';
        clear Sjj
    end
    param.axi           = param.axi/scaleJ;
    scale               = scaleLvj^2/(scaleV*scaleJ);    
end
%% General warm-start by cross-validating diagonal initial values
%  search for the optimal initial values sigma2xi0 Sigmajj0 by cross-validating sigma2xi and sigma2j*Sigmajj,
%  Sigmajj is defined above in any of the cases prew = 0 (Sigmajj = Iq) or prew = 1 (Sigmajj is determined by
%  the twostep_lasso_caller), the optimal criteria is due to the hyperparameters posterior distribution

[sigma2xi0,Sigmajj0,llh0] = higgs_warm_start(Svv,Lvj,Thetajj,Sigmajj,param);
end
