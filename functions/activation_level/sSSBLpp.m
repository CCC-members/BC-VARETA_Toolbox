function [s2j,sigma2j,Tjv,Svv,scaleJ,scaleLvj] = sSSBLpp(Svv,Lvj,param)
W             = param.W;
parcellation  = param.parcellation;
run_bash_mode = param.run_bash_mode;
flag          = param.flag;

% Elastic Net_Sparse Bayesian Learning
%% Standardization
LvjW          = Lvj*W;
WLjv          = LvjW';

%% Static parameters
[p,q]         = size(LvjW);
Ip            = spdiags(ones(p,1),0,p,p);
a             = length(parcellation);
maxiter_outer = 60;
maxiter_inner = 1;
s_alpha       = q;
r_alpha       = q;
s_rho         = 1;
r_rho         = 1;

%% Scaling Lead Field
scaleLvj      = sqrt(trace(LvjW*WLjv)/p);
LvjW          = LvjW/scaleLvj;
WLjv          = WLjv/scaleLvj;

%% Initial values
sigma2j_post  = zeros(q,1);
s2j           = zeros(q,1);
etha          = zeros(q,1);
sigma2j       = 1E0*ones(q,1);
sigma2x       = 1E0;
alpha         = 1E0;
rho           = 1E0;

%% Scaling data using a first pass solution
sigma2jWLjv   = spdiags(2*sigma2j,0,q,q)*WLjv;
sigma2j_post0 = sigma2jWLjv/(LvjW*sigma2jWLjv+sigma2x*Ip);
LvjWsigma2j   = sigma2jWLjv';
% Compute diagonals of the posterior covariance
for count_gen = 1:q
    sigma2j_post(count_gen) = 2*sigma2j(count_gen) - sigma2j_post0(count_gen,:)*LvjWsigma2j(:,count_gen);
end
% Iterative Transfer Operator
Tjv           = (1/sigma2x).*(sigma2jWLjv-sigma2j_post0*(LvjWsigma2j*WLjv));
% Compute empirical covariance
SvvTvj        = Svv*Tjv';
for count_gen = 1:q
    s2j(count_gen) = abs(Tjv(count_gen,:)*SvvTvj(:,count_gen));
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
% Update alpha
idx_alpha     = find(sigma2j_bar > 0);
alpha         = (length(idx_alpha)/2 + s_alpha)/(sum((s2j(idx_alpha) + sigma2j_post(idx_alpha))./(sigma2j_bar(idx_alpha))) + r_alpha);
    
% scaleJ        = mean(abs(s2j + sigma2j_post))/mean(abs(sigma2j_post));
% scaleJ        = mean(abs(s2j))/max(abs(sigma2j_post));
scaleJ        = mean(sigma2j);
Svv           = Svv/scaleJ;

%% Initial values
sigma2j_post  = zeros(q,1);
s2j           = zeros(q,1);
etha          = zeros(q,1);
sigma2j       = 1E0*ones(q,1);
sigma2x       = 1E0;
alpha         = 1E0;
rho           = 1E0;

%% Outer cycle
if(~run_bash_mode)
    process_waitbar = waitbar(0,'Please wait...');
end
disp(flag);
fprintf(1,strcat('-->> sSSBL++ (',param.str_band,') process: %3d%%\n'),0);
for cont1 = 1:maxiter_outer
    %% Inner cycle
    for cont11 = 1:maxiter_inner
        sigma2jWLjv                 = spdiags(2*sigma2j,0,q,q)*WLjv;
        sigma2j_post0               = sigma2jWLjv/(LvjW*sigma2jWLjv+sigma2x*Ip);
        LvjWsigma2j                 = sigma2jWLjv';
        % Compute diagonals of the posterior covariance
        for count_gen = 1:q
            sigma2j_post(count_gen) = 2*sigma2j(count_gen) - sigma2j_post0(count_gen,:)*LvjWsigma2j(:,count_gen);
        end
        % Iterative Transfer Operator
        Tjv                         = (1/sigma2x).*(sigma2jWLjv-sigma2j_post0*(LvjWsigma2j*WLjv));
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
    
    if(~run_bash_mode)
        text = replace(param.str_band,'_','-');
        waitbar((cont1)/(maxiter_outer),process_waitbar,strcat("sSSBL++ (",text,") process: ",num2str(fix((cont1/maxiter_outer)*100)-1),"%"));
    end
    fprintf(1,'\b\b\b\b%3.0f%%',(cont1)/(maxiter_outer)*100);
end
fprintf(1,'\n');
%% Destandardization using a final pass solution
disp("-->> Running sSSBL++ destandardization.");
sigmajW                     = spdiags(sqrt(2*sigma2j),0,q,q)*W';
Wsigma2jW                   = sum(sigmajW.^2,1)';
sigma2jWLjv                 = spdiags(2*sigma2j,0,q,q)*WLjv;
sigma2j_post0               = (W*sigma2jWLjv)/(LvjW*sigma2jWLjv+sigma2x*Ip);
LvjWsigma2j                 = sigma2jWLjv';
LvjWsigma2jW                = LvjWsigma2j*W';

% Only save the diagonals of the Posterior Covariance
for count_gen = 1:size(Lvj,2)
    sigma2j_post(count_gen) = Wsigma2jW(count_gen) - sigma2j_post0(count_gen,:)*LvjWsigma2jW(:,count_gen);
end

% Iterative Transference Operator
Tjv                        = (1/sigma2x).*(W*sigma2jWLjv-sigma2j_post0*(LvjWsigma2j*WLjv));
if(~run_bash_mode)
    text = replace(param.str_band,'_','-');
    waitbar(1,process_waitbar,strcat("sSSBL++ (",text,") process: ",num2str(100),"%"));
end

% Compute 'miu' for all slices of 'V'
SvvTvj                     = Svv*Tjv';
for count_gen = 1:size(Lvj,2)
    s2j(count_gen)         = abs(Tjv(count_gen,:)*SvvTvj(:,count_gen));
end
fprintf(1,'\n');
if(~run_bash_mode)
    delete(process_waitbar);
end
end