function [Theta,llh] = hg_lasso_lqa(Psi,param)
%% Hermitian Graphical LASSO by Local Quadratic Approximation
%  computations based on Singula Value Decomposition
%% hg-lasso initialization 
a                 = param.a;
A                 = param.A;
maxiter           = param.maxiter;
nu                = param.nu;
m                 = param.m;
q                 = param.q;
run_bash_mode     = param.run_bash_mode;
m2                = m^2;
a2                = a^2;
A2                = A.^2;
Psi_inv.U         = Psi.U;
Psi_inv.d         = 1./Psi.d;
Psi_inv.X         = Psi_inv.U*spdiags(Psi_inv.d,0,q,q)*Psi_inv.U';
Theta             = Psi_inv;
idx               = (A > 0);
idx0              = (A == 0);
gamma2            = zeros(q);
llh               = zeros(maxiter,1);
%% Main cycle
if(~run_bash_mode)
    process_waitbar = waitbar(0,'Please wait...');
end
fprintf(1,'-->> Running Hermitian Graphical LASSO: %3d%%\n',0);
for k_inner = 1:maxiter
    %% Estimation of variances Gamma of Gaussian Mixtures prior    
    det               = 1 + 4*m2*a2*A2(idx).*abs(Theta.X(idx)).^2;
    gamma2(idx)       = (sqrt(det) - 1)./(2*m*a2*A2(idx));
    gamma2(idx0)      = m*abs(Theta.X(idx0)).^2;
    %%
    %% Standarization of the Empirical Covariance Matrix
    ninf              = max(gamma2(:));
    gamma2ninf        = gamma2/ninf;
    st_factor1        = nu*ninf^(-1/2)*gamma2ninf.^(1/2).*Psi_inv.X; % lqa corrected factor 1
    st_factor2        = (1 + gamma2ninf*nu); % lqa corrected factor 2
    [Psi_st_inv]      = higgs_eigendecomposition(st_factor1./st_factor2,param); % eigendecomposition inverted standard empirical covariance matrix
    Psi_st.U          = Psi_st_inv.U; % eigenvectors of standard empirical covariance matrix
    Psi_st.d          = 1./Psi_st_inv.d; % eigenvalues of standard empirical covariance matrix
    %%
    %% Standard Precision Matrix estimator
    Theta_st.U        = Psi_st.U; % eigenvectors of standard precision matrix Riccati estimator
    Theta_st.d        = (sqrt(Psi_st.d.^2 + 4) - Psi_st.d)/2; % eigenvalues of standard precision matrix Riccati estimator
    Theta_st.X        = Theta_st.U*spdiags(Theta_st.d,0,q,q)*Theta_st.U'; % standard precision matrix Riccati estimator
    %%
    %% Unstandarized Precision Matrix estimator
    [Theta_tmp]       = higgs_eigendecomposition(gamma2.^(1/2).*Theta_st.X,param); % precision matrix eigendecomposition
    llh_tmp           = sum(log(Theta_tmp.d)) - abs(sum(sum(Theta_tmp.X.*transpose(Psi.X),2))) - a*sum(abs(A(:).*Theta_tmp.X(:))); % hglasso likelihood
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