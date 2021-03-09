function [Theta,llh] = hg_lasso_lqa1(Psi,m,A,a,nu,maxiter,param)
%% Hermitian Graphical LASSO by Local Quadratic Approximation
%  computations based on Singula Value Decomposition
%% hg-lasso initialization 
q                 = length(Psi);
Iq                = eye(q);
m2                = m^2;
a2                = a^2;
A2                = A.^2;
Psi_inv           = Iq/Psi; 
Psi_inv           = (Psi_inv + Psi_inv')/2; 
use_gpu           = param.use_gpu;
run_bash_mode     = param.run_bash_mode;
if(use_gpu)
    gpuPsi_st         = gpuArray(zeros(q));
    gpuTheta_tmp      = gpuArray(zeros(q));
end
idx               = (A > 0);
idx0              = (A == 0);
gamma2            = zeros(q);
llh               = zeros(maxiter,1);
Theta             = Psi_inv;
%% Main cycle
if(~run_bash_mode)
    process_waitbar = waitbar(0,'Please wait...');
end
fprintf(1,'-->> Running Hermitian Graphical LASSO: %3d%%\n',0);
for k_inner = 1:maxiter
    %% Estimation of variances Gamma of Gaussian Mixtures prior    
    DET               = 1 + 4*m2*a2*A2(idx).*abs(Theta(idx)).^2;
    gamma2(idx)       = (sqrt(DET) - 1)./(2*m*a2*A2(idx));
    gamma2(idx0)      = m*abs(Theta(idx0)).^2;
    %%
    %% Standarization of the Empirical Covariance Matrix
    ninf              = max(gamma2(:));
    gamma2ninf        = gamma2/ninf;
    st_factor1        = nu*ninf^(-1/2)*gamma2ninf.^(1/2).*Psi_inv; % lqa corrected factor 1
    st_factor2        = (1 + gamma2ninf*nu); % lqa corrected factor 2
    %%
    Psi_st_inv        = st_factor1./st_factor2; % Inverted Standard Empirical Covariance Matrix
    Psi_st            = Iq/Psi_st_inv; % Standard Empirical Covariance Matrix
    Psi_st            = (Psi_st + Psi_st')/2; 
    %%
    %% Standard Precision Matrix estimator
    if(~use_gpu)
        [U,D,V]           = svd(Psi_st);
    else
        gpuPsi_st(:,:)    = Psi_st;
        [U,D,V]           = svd(gpuPsi_st);
        [U,D,V]           = gather(U,D,V);
    end
    d                 = abs(diag(D));
    Theta_st          = (1/2)*U*spdiags(sqrt(d.^2 + 4) - d,0,q,q)*V';
    Theta_st          = (Theta_st + Theta_st')/2;
    %%
    %% Unstandarized Precision Matrix estimator
    Theta_tmp         = gamma2.^(1/2).*Theta_st;
    if(~use_gpu)
        [~,D,~]           = svd(Theta_tmp);
    else
        gpuTheta_tmp(:,:) = Theta_tmp;
        [U,D,V]           = svd(gpuTheta_tmp);
        [~,D,~]           = gather(U,D,V);
    end
    d                 = abs(diag(D));
    llh_tmp           = sum(log(d)) - sum(abs(diag(Theta_tmp*Psi))) -a*sum(abs(A(:).*Theta_tmp(:)));
    if k_inner == 1
        Theta         = Theta_tmp;
        llh(k_inner)  = llh_tmp;
    elseif k_inner > 1 && llh_tmp >= llh(k_inner-1)
        Theta         = Theta_tmp;
        llh(k_inner)  = llh_tmp;
    else
        llh(k_inner:end) = llh(k_inner-1);
        if(~run_bash_mode)
           waitbar((maxiter-1)/(maxiter),process_waitbar,strcat("Running Hermitian Graphical LASSO: ",num2str(fix(((maxiter-1)/maxiter)*100)),"%"));
        end
        fprintf(1,'\b\b\b\b%3.0f%%',(k_inner/maxiter)*100);
        break
    end
    fprintf(1,'\b\b\b\b%3.0f%%',(k_inner/maxiter)*100);
end

fprintf(1,'\b\b\b\b%3.0f%%',100);
fprintf(1,'\n');
if(~run_bash_mode)
    waitbar(1,process_waitbar,strcat("Running Hermitian Graphical LASSO: ",num2str(100),"%"));
    delete(process_waitbar)
end

end