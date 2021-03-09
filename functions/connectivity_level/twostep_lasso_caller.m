function [Thetajj,Sigmajj] = twostep_lasso_caller(Psijj,param)
m              = param.m;
q              = param.q;
aj             = param.aj;
Ajj            = param.Ajj;
nu             = param.nu;
rth1           = param.rth1;
rth2           = param.rth2;
eigreg         = param.eigreg;
maxiter_inner  = param.maxiter_inner;
%% Compute hg_lasso_lqa
disp("Compute hg_lasso_lqa. This operation may take a long time.");
if(~param.use_gpu)
    [U,D]             = eig(Psijj);
    V                 = Iq/U;
    d                 = diag(D);
    if min(d) < 0
        d             = d + abs(min(d));
    end
    dfix              = d + max(d)*eigreg;
    [U,dfix,V]        = gather(U,dfix,V);
    clear D d
else
    gpuPsijj          = gpuArray(Psijj);
    [U,D]             = eig(gpuPsijj);
    V                 = Iq/U;
    d                 = diag(D);
    if min(d) < 0
        d             = d + abs(min(d));
    end
    dfix              = d + max(d)*eigreg;
    [U,dfix,V]        = gather(U,dfix,V);
    clear D d
end
Psijj_fix      = U*spdiags(dfix,0,q,q)*V;
Psijj_fix      = (Psijj_fix + Psijj_fix')/2;
Psijj_fix      = Psijj_fix/sqrt(sum(abs(diag(Psijj_fix*Psijj_fix')))/length(Psijj_fix));
[Thetajj_lasso,~] = hg_lasso_lqa1(Psijj_fix,m,Ajj,aj,nu,maxiter_inner,param);
Thetajj_unb    = 2*Thetajj_lasso - Thetajj_lasso*Psijj_fix*Thetajj_lasso;
Thetajj_unb    = (Thetajj_unb + Thetajj_unb')/2;
Thetajj_var    = sqrt(abs(diag(Thetajj_lasso))*abs(diag(Thetajj_lasso))' + abs(Thetajj_lasso).^2);
clearvars Thetajj_lasso Psijjfix;
Thetajj_var    = (Thetajj_var + Thetajj_var')/2;
%% Precompute hg_lasso_lqa equivalent to Frobenious solution
d_frob         = (sqrt(dfix.^2 + 4*aj^2) - dfix)/(2*aj^2);
Thetajj        = U*spdiags(d_frob,0,q,q)*V;
clearvars U V;
Thetajj        = (Thetajj + Thetajj')/2;
%% try likelihood evolution for the mask given by all thresholds

run_bash_mode       = param.run_bash_mode;
if(~run_bash_mode)
    process_waitbar = waitbar(0,'Please wait...');
end
%% check only partial likelihood
rth_grid         = rth1:0.1:rth2;
llhjj_grid       = zeros(length(rth_grid),1);
fprintf(1,'-->> Computing optimal Ryleigh threshold: %3d%%\n',0);
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
    end
    Thetajj_mask          = U*spdiags(dfix,0,q,q)*V;
    Thetajj_mask          = (Thetajj_mask + Thetajj_mask')/2;
    llhjj_grid(rth_count) = sum(log(abs(dfix))) - sum(abs(diag(Thetajj_mask*Psijj_fix))) - aj*sum(abs(Thetajj_mask(:)));
    fprintf(1,'\b\b\b\b%3.0f%%',(rth_count/length(rth_grid))*100 - 1);
    if(~run_bash_mode)
           waitbar(rth_count/length(rth_grid),process_waitbar,strcat("Computing optimal Ryleigh threshold: ",num2str(fix((rth_count/length(rth_grid))*100)-1),"%"));
    end
end
[val,id_rth]     = max(llhjj_grid);
clearvars llhjj_grid;
id_rth           = min(id_rth);
rth              = rth_grid(id_rth);
Thetajj(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var)))) = 0;
clearvars Thetajj_var Thetajj_unb;
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
Thetajj          = U*spdiags(dfix,0,q,q)*V;
Sigmajj          = V*spdiags(1./dfix,0,q,q)*U;
clearvars U D V;
Thetajj          = (Thetajj + Thetajj')/2;
Sigmajj          = (Sigmajj + Sigmajj')/2;
fprintf(1,'\b\b\b\b%3.0f%%',100);
fprintf(1,'\n');
if(~run_bash_mode)
    waitbar(1,process_waitbar,strcat("Computing optimal Ryleigh threshold: ",num2str(100),"%"));
    delete(process_waitbar)
end
end
