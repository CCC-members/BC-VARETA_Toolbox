function [s2j,sigma2j_post,Tjv,Svv,Lvj,sigma2j,sigma2j_post0,Lvjsigma2j,scaleJ,scaleLvj] = sSSBLpp_ultralarge(Svv,Lvj,parcellation)
% Elastic Net_Sparse Bayesian Learning
%% Static parameters
[p,q]         = size(Lvj);
Ip            = spdiags(ones(p,1),0,p,p);
a             = length(parcellation);
maxiter_outer = 60;
maxiter_inner = 1;
s_alpha       = q;
r_alpha       = q;
s_rho         = 1;
r_rho         = 1;

%% Declaring dinamic variables
sigma2j_post  = zeros(q,1);
s2j           = zeros(q,1);
etha          = zeros(q,1);

%% Initial values
sigma2j       = 1E0*ones(q,1);
sigma2x       = 1E0;
alpha         = 1E0;
rho           = 1E0;

%% Scaling Lead Field
Ljv           = Lvj';
scaleLvj      = sqrt(trace(Lvj*Ljv)/p);
Lvj           = Lvj/scaleLvj;
Ljv           = Ljv/scaleLvj;

%% Scaling data using a first pass solution
sigma2jLjv    = spdiags(2*sigma2j,0,q,q)*Ljv;
sigma2j_post0 = sigma2jLjv/(Lvj*sigma2jLjv+sigma2x*Ip);
Lvjsigma2j    = sigma2jLjv';
% Compute diagonals of the posterior covariance
for count_gen = 1:q
    sigma2j_post(count_gen) = 2*sigma2j(count_gen) - sigma2j_post0(count_gen,:)*Lvjsigma2j(:,count_gen);
end
% Iterative Transfer Operator
Tjv           = (1/sigma2x).*(sigma2jLjv-sigma2j_post0*(Lvjsigma2j*Ljv));
% Compute empirical covariance
SvvTvj        = Svv*Tjv';
for count_gen = 1:q
    s2j(count_gen) = abs(Tjv(count_gen,:)*SvvTvj(:,count_gen));
end
scaleJ        = mean(abs(s2j + sigma2j_post))/mean(abs(sigma2j_post));
Svv           = Svv/scaleJ;

%% Outer cycle
process_waitbar = waitbar(0,'Please wait...');
for cont1 = 1:maxiter_outer
    waitbar((cont1)/(maxiter_outer),process_waitbar,strcat('sSSBLpp loop # ',num2str(cont1)));
    disp(['sSSBLpp loop # ',num2str(cont1)])
    %% Inner cycle
    for cont11 = 1:maxiter_inner
        sigma2jLjv                  = spdiags(2*sigma2j,0,q,q)*Ljv;
        sigma2j_post0               = sigma2jLjv/(Lvj*sigma2jLjv+sigma2x*Ip);
        Lvjsigma2j                  = sigma2jLjv';
        % Compute diagonals of the posterior covariance
        for count_gen = 1:q
            sigma2j_post(count_gen) = 2*sigma2j(count_gen)-sigma2j_post0(count_gen,:)*Lvjsigma2j(:,count_gen);
        end
        % Iterative Transfer Operator
        Tjv                         = (1/sigma2x).*(sigma2jLjv-sigma2j_post0*(Lvjsigma2j*Ljv));
        % Compute empirical covariance
        SvvTvj                      = Svv*Tjv';
        for count_gen = 1:q
            s2j(count_gen)          = abs(Tjv(count_gen,:)*SvvTvj(:,count_gen));
        end
        % Update Gammas
        for area = 1:a
            idx_area                = parcellation{area};
            etha(idx_area)          = sqrt((1./4).^2+(alpha*rho).*sum(s2j(idx_area) + sigma2j_post(idx_area)))-1./4;
        end
        idx_etha                    = find((s2j + sigma2j_post)<0);
        etha(idx_etha)              = 0;
        gamma                       = rho + etha;
        sigma2j_bar                 = etha./gamma;
        sigma2j                     = (1/(2*alpha))*sigma2j_bar;
    end
    %% Update alpha
    idx_alpha  = find(sigma2j_bar > 0);
    alpha      = (length(idx_alpha)/2 + s_alpha)/(sum((s2j(idx_alpha) + sigma2j_post(idx_alpha))./(sigma2j_bar(idx_alpha))) + r_alpha);
    %% Update rho
    f_aux      = @(k_aux) r_rho + sum(ones(q,1)./(1-sigma2j_bar))/q - (s_rho - 1/2)/k_aux - trascendent_term(k_aux);
    rho        = fzero(f_aux,[0.000001 70000]);
    sigma2j    = (1/(2*alpha))*sigma2j_bar;
end
delete(process_waitbar);
end