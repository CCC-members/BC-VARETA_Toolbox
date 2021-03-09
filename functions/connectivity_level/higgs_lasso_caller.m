function [Thetajj,Sigmajj,llh_inner] = higgs_lasso_caller(Psijj,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param)
%% Effectutates the higgs-lasso maximization step
%  attempts to solve two ambigities in the higgs-lasso maximization 1st- the hg-lasso sparsity biases the results
%  in local and global computations of Thetajj,
q              = param.q;
m              = param.m;
aj             = param.aj;
Ajj            = param.Ajj;
nu             = param.nu;
rth1           = param.rth1;
rth2           = param.rth2;
maxiter_inner  = param.maxiter_inner;
ntry           = param.ntry;
eigreg         = param.eigreg;
%% Compute hg_lasso_lqa
%  it first corrects the eigenvalues of the effective empirical covariance Psijj and normalize it by the infimun
%  norm Psijjfix, the scaling guarantees more stable computations of the hg-lasso-lqa, continues with performing
%  hg-lasso-lqa Thetajj_lasso and the unbiased statistics Thetajj_unb and Thetajj_var
if(~param.use_gpu)
    [U,D]          = eig(Psijj);
    V              = Iq/U;
    d              = diag(D);
    if min(d) < 0
        d          = d + abs(min(d));
    end
    dfix           = d + max(d)*eigreg;
    [U,dfix,V]     = gather(U,dfix,V);
    clear D d
else
    gpuPsijj       = gpuArray(Psijj);
    [U,D]          = eig(gpuPsijj);
    V              = Iq/U;
    d              = diag(D);
    if min(d) < 0
        d          = d + abs(min(d));
    end
    dfix           = d + max(d)*eigreg;
    [U,dfix,V]     = gather(U,dfix,V);
    clear gpuPsijj
end
Psijj_fix      = U*spdiags(dfix,0,q,q)*V;
Psijj_fix      = Psijj_fix/sqrt(sum(abs(diag(Psijj_fix*Psijj_fix')))/length(Psijj_fix));
[Thetajj_lasso,llh_inner] = hg_lasso_lqa1(Psijj_fix,m,Ajj,aj,nu,maxiter_inner,param);
Thetajj_unb    = 2*Thetajj_lasso - Thetajj_lasso*Psijj_fix*Thetajj_lasso;
Thetajj_unb    = (Thetajj_unb + Thetajj_unb')/2;
Thetajj_var    = sqrt(abs(diag(Thetajj_lasso))*abs(diag(Thetajj_lasso))' + abs(Thetajj_lasso).^2);
Thetajj_var    = (Thetajj_var + Thetajj_var')/2;
%% Precompute hg_lasso_lqa equivalent to ridge solution
%  this is the closest solution to the hermitian graphical lasso with parameter aj and equivalent to the first
%  iteration of the lqa with gamma == 1, if for a given iteration the lasso model biases the computations and
%  the global hyperparameters posterior likelihood does not increases in ntry iterations higgs-lasso still goes
%  on with higgs-ridge prewarming, this approximation is a result of the empirical study of higgs-ridge with
%  parameter aj^2 (aj same as lasso aj = sqrt(log(q)/m)) in terms of stability and similarity to the unbiased
%  statistics given by Jankova conditions,

d_frob         = (sqrt(dfix.^2 + 4*aj^2) - dfix)/(2*aj^2);
Thetajj        = U*spdiags(d_frob,0,q,q)*V';
Thetajj        = (Thetajj + Thetajj')/2;
%% try likelihood evolution for the mask given by all thresholds
% evaluate performance of Rayleigh thresholds in terms of the global hyperparameters posterior likelihood, the 
% overlapping of the distributions for null hypothesis values (Rayleigh) and the alternative introduces ambiguity 
% in the decision for a given Rayleigh threshold we move rth between the value corresponding to the peak of the 
% Rayleigh distribution (rth1)

if ntry == 0
    %% check only partial likelihood
    rth_grid         = rth1:0.1:rth2;
    llhjj_grid       = zeros(length(rth_grid),1);
    for rth_count = 1:length(rth_grid)
        %% Zeros mask
        rth                   = rth_grid(rth_count);
        mask                  = find(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var))));
        Thetajj_mask          = Thetajj;
        Thetajj_mask(mask)    = 0;
        if(~param.use_gpu)
            [U,D]             = eig(Thetajj_mask);
            V                 = Iq/U;
            d                 = diag(D);
            if min(d) < 0
                d             = d + abs(min(d));
            end
            dfix              = d + max(d)*eigreg;
            [U,dfix,V]        = gather(U,dfix,V);
            clear D d
        else
            gpuThetajj_mask   = gpuArray(Thetajj_mask);
            [U,D]             = eig(gpuThetajj_mask);
            V                 = Iq/U;
            d                 = diag(D);
            if min(d) < 0
                d             = d + abs(min(d));
            end
            dfix              = d + max(d)*eigreg;
            [U,dfix,V]        = gather(U,dfix,V);
            clear D d
            clear gpuThetajj_mask
        end
        Thetajj_mask          = U*spdiags(dfix,0,q,q)*V;
        Thetajj_mask          = (Thetajj_mask + Thetajj_mask')/2;
        llhjj_grid(rth_count) = sum(log(abs(d_fix))) - sum(abs(sum(Thetajj_mask.*transpose(Psijj_fix),2))) - aj*sum(abs(Thetajj_mask(:)));
    end
    [val,id_rth]     = max(llhjj_grid);
    id_rth           = min(id_rth);
    rth              = rth_grid(id_rth);
else
    %% check global likelihood
    [rth]            = higgs_lasso_stabilizer(Thetajj_unb,Thetajj_var,Svv,Lvj,Ljv,sigma2xi0,Sigmajj0,llh0,param);
end
Thetajj(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var)))) = 0;
if(~param.use_gpu)
    [U,D]             = eig(Thetajj);
    V                 = Iq/U;
    d                 = diag(D);
    if min(d) < 0
        d             = d + abs(min(d));
    end
    dfix              = d + max(d)*eigreg;
    [U,dfix,V]        = gather(U,dfix,V);
    clear D d
else
    gpuThetajj        = gpuArray(Thetajj);
    [U,D]             = eig(gpuThetajj);
    V                 = Iq/U;
    d                 = diag(D);
    if min(d) < 0
        d             = d + abs(min(d));
    end
    dfix              = d + max(d)*eigreg;
    [U,dfix,V]        = gather(U,dfix,V);
    clear D d
end
Thetajj        = U*spdiags(dfix,0,q,q)*V;
Thetajj        = (Thetajj + Thetajj')/2;
Sigmajj        = V*spdiags(1./dfix,0,q,q)*U;
Sigmajj        = (Sigmajj + Sigmajj')/2;
end
