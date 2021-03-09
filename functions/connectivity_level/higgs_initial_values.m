function [Svv,Lvj,Ljv,scale,scaleLvj,sigma2xi0,Sigmajj0,llh0,param] = higgs_initial_values(Svv,Lvj,param)
%% Parameters
prew            = param.prew;
m               = param.m;
p               = param.p;
Op              = param.Op;
q               = param.q;
Iq              = param.Iq;
aj              = param.aj;
penalty         = param.penalty;
eigreg          = param.eigreg;
%%
%% Lead Field and Data scaling by frobenious norm equivalent to the largest singular values or the
% average variance of the observations, helps to treat the problen in a standard scale for every case,
% a projected identity matrix Svv = Lvj*Iq*Ljv under this conditions will generate a signal with
% (sum(abs(diag(Svv)))/p) = 1 if we fix the value of param.axi = 1E-2 the inferior admisible noise
% will be 1% of the average signal amplitude

scaleLvj        = sqrt(sum(abs(diag(Lvj*Lvj')))/p);
Lvj             = Lvj/scaleLvj;
Ljv             = Lvj';
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
    
    Sigmajj0        = Iq;
    sigma2xi0       = Op;
    [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_post] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,param);
    scaleJ          = (sum(abs(diag(Sjj)))/q)/(sum(abs(diag(Sigmajj_post)))/q);
    Svv             = Svv/scaleJ;
    param.axi       = param.axi/scaleJ;
    scale           = scaleLvj^2/(scaleV*scaleJ);
    Thetajj         = Iq;
    Sigmajj         = Iq;
    clearvars Sjj;
elseif prew == 1 % prewarming
    %  - if prew = 1 we use an initial subspace two-step solution based on the crossspectral enet-ssbl
    %  (equivalent to an univariate higgs) and hermitian graphical model, this provides a multivariate
    %  initialization to the EM or prewarming, crossspectral enet-ssbl also adjusts the scale by scaleJ
    %  for every penalty a two-step graphical computation is performed to produce the initial structure
    %  of Thetajj and Sigmajj
    
    disp('-->> prewarming');
    parcellation        = [];
    counter             = 1;
    for ii = 1:1:q
        parcellation{counter} = ii;
        counter               = counter + 1;
    end
    param.Nsamp         = m;
    param.parcellation  = parcellation;
    param.W             = eye(length(Lvj));
    param.flag          = "-->> Running prewarming, source activity";
    [sigma2j,sigma2j_post,Tjv,Svv,~,scaleJ,~] = sSSBLpp(Svv,Lvj,param);
    Tvj                 = Tjv';
    Sjj                 = Tjv*Svv*Tvj;
    clear Tvj Tjv;
    Sjj                 = (Sjj + Sjj')/2;
    disp("-->> Running prewarming, source connectivity. This operation may take a long time.");
    if penalty == 0 % naive
        if(~param.use_gpu)
            [U,D]             = eig(Sjj);
            V                 = Iq/U;
            d                 = diag(D);
            if min(d) < 0
                d             = d + abs(min(d));
            end
            dfix              = d + max(d)*eigreg;
            [U,dfix,V]        = gather(U,dfix,V);
            clear D d
            clear Sjj
        else
            gpuSjj            = gpuArray(Sjj);
            [U,D]             = eig(gpuSjj);
            V                 = Iq/U;
            d                 = diag(D);
            if min(d) < 0
                d             = d + abs(min(d));
            end
            dfix              = d + max(d)*eigreg;
            [U,dfix,V]        = gather(U,dfix,V);
            clear D d
            clear gpuSjj Sjj
        end
        Thetajj               = V*spdiags(1./dfix,0,q,q)*U;
        Thetajj               = (Thetajj + Thetajj')/2;
        Sigmajj               = U*spdiags(dfix,0,q,q)*V;
        Sigmajj               = (Sigmajj + Sigmajj')/2;
    elseif penalty == 1 % lasso
        [Thetajj,Sigmajj]     = twostep_lasso_caller(Sjj,param);
    elseif penalty == 2 % ridge
        if(~param.use_gpu)
            [U,D]             = svd(Sjj);
            V                 = Iq/U;
            d                 = diag(D);
            if min(d) < 0
                d             = d + abs(min(d));
            end
            dfix              = d + max(d)*eigreg;
            [U,dfix,V]        = gather(U,dfix,V);
            clear D d
            clear Sjj
        else
            gpuSjj            = gpuArray(Sjj);
            [U,D]             = eig(gpuSjj);
            V                 = Iq/U;
            d                 = diag(D);
            if min(d) < 0
                d             = d + abs(min(d));
            end
            dfix              = d + max(d)*eigreg;
            [U,dfix,V]        = gather(U,dfix,V);
            clear D d
            clear gpuSjj
            clear Sjj
        end
        d_frob                = (sqrt(dfix.^2 + 4*aj^2) - dfix)/(2*aj^2);
        Thetajj               = U*spdiags(d_frob,0,q,q)*V;
        Thetajj               = (Thetajj + Thetajj')/2;
        Sigmajj               = V*spdiags(1./d_frob,0,q,q)*U;
        Sigmajj               = (Sigmajj + Sigmajj')/2;
    end
    param.axi           = param.axi/scaleJ;
    scale               = scaleLvj^2/(scaleV*scaleJ);    
end
%% General warm-start by cross-validating diagonal initial values
%  search for the optimal initial values sigma2xi0 Sigmajj0 by cross-validating sigma2xi and sigma2j*Sigmajj,
%  Sigmajj is defined above in any of the cases prew = 0 (Sigmajj = Iq) or prew = 1 (Sigmajj is determined by
%  the twostep_lasso_caller), the optimal criteria is due to the hyperparameters posterior distribution

[sigma2xi0,Sigmajj0,llh0] = higgs_warm_start(Svv,Lvj,Ljv,Thetajj,Sigmajj,param);
end
