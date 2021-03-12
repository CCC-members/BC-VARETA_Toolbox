function [Theta,llh] = hg_lasso_lqa2(Psi,m,A,a,nu,maxiter)
%% Hermitian Graphical LASSO by Local Quadratic Approximation
%  computations based on Eigenvalue and Eigenvector Decomposition 
%% hg-lasso initialization  
q                 = length(Psi);
Iq                = eye(q);
m2                = m^2;
a2                = a^2;
A2                = A.^2;
Psi_inv           = Iq/Psi; 
Psi_inv           = (Psi_inv + Psi_inv')/2; 
idx               = (A > 0);
idx0              = (A == 0);
gamma2            = zeros(q);
llh               = zeros(maxiter,1);
Theta             = Psi_inv;
%% Main cycle
for k_inner = 1:maxiter
    %% Estimation of variances Gamma of Gaussian Mixtures prior
    DET              = 1 + 4*m2*a2*A2(idx).*abs(Theta(idx)).^2;
    gamma2(idx)      = (sqrt(DET) - 1)./(2*m*a2*A2(idx)); 
    gamma2(idx0)     = m*abs(Theta(idx0)).^2;
    %%
    %% Standarization of the Empirical Covariance Matrix
    ninf             = max(gamma2(:));
    gamma2ninf       = gamma2/ninf;
    st_factor1       = nu*ninf^(-1/2)*gamma2ninf.^(1/2).*Psi_inv; % lqa corrected factor 1
    st_factor2       = (1 + gamma2ninf*nu); % lqa corrected factor 2
    %%
    Psi_st_inv       = st_factor1./st_factor2; % Inverted Standard Empirical Covariance Matrix
    Psi_st           = Iq/Psi_st_inv; % Standard Empirical Covariance Matrix
    Psi_st           = (Psi_st + Psi_st')/2;    
    %%
    %% Standard Precision Matrix estimator
    [U,D]            = eig(Psi_st);
    d                = abs(diag(D));
    Theta_st         = (1/2)*U*diag(sqrt(d.^2 + 4) - d)*U';
    Theta_st         = (Theta_st + Theta_st')/2;
    %%
    %% Unstandarized Precision Matrix estimator
    Theta_tmp        = gamma2.^(1/2).*Theta_st;
    %%
    [U,D,V]          = eig(Theta_tmp);
    d                = abs(diag(D));
    llh_tmp          = sum(log(d)) - sum(abs(diag(Theta_tmp*Psi))) -a*sum(abs(A(:).*Theta_tmp(:)));
    if k_inner == 1
        Theta        = Theta_tmp;
        llh(k_inner) = llh_tmp;
    elseif k_inner > 1 && llh_tmp >= llh(k_inner-1)
        Theta        = Theta_tmp;
        llh(k_inner) = llh_tmp;
    else
        llh(k_inner:end) = llh(k_inner-1);
        break
    end
end
end