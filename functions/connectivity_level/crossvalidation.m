function llh = crossvalidation(grid_sigma2xi,grid_sigma2j,Sigmajj,Thetajj,Svv,Lvj,Ljv,param)
%% Parameters
p             = param.p;
Ip            = param.Ip;
q             = param.q;
Iq            = param.Iq;
axi           = param.axi;
aj            = param.aj;
Ajj           = param.Ajj;
penalty       = param.penalty;
%% Eigenvalues for likelihood
disp("-->> Running prewarming, preparing constants for crossvalidation");
if(~param.use_gpu)
    [~,D]         = eig(Thetajj);
else
    gpuThetajj    = gpuArray(Thetajj);
    [~,D]         = eig(gpuThetajj);
    D             = gather(D);
    clear gpuThetajj
end
d             = diag(D);
dinv          = 1./d;
Ljj           = Ljv*Lvj;
[V,C]         = eig(Ljj,Thetajj); % Ljj = Thetajj*V*C*V^(-1)
if(~param.use_gpu)
    Vinv          = Iq/V;
else
    gpuIq         = gpuArray(Iq);
    gpuV          = gpuArray(V);
    Vinv          = gpuIq/gpuV;
    Vinv          = gather(Vinv);
    clear gpuIq gpuV
end
c             = abs(diag(C));
SigmajjLjv    = Sigmajj*Ljv;
p2Tjv         = Vinv*SigmajjLjv;
llh_norm1     = sum(abs(Ajj(:).*Thetajj(:)));
llh_norm2     = sum(sum(abs(Thetajj.*transpose(Thetajj)),2));
%% Cross-validation cycle
llh = zeros(length(grid_sigma2xi),length(grid_sigma2j));
fprintf(1,'-->> Running prewarming, crossvalidating higgs initial values: %3d%%\n',0);
count = 1;
for count_sigma2xi = 1:length(grid_sigma2xi)
    sigma2xi  = 10^(grid_sigma2xi(count_sigma2xi));
    for count_sigma2j = 1:length(grid_sigma2j)
        sigma2j  = 10^(grid_sigma2j(count_sigma2j));
        % Expectation
        c_inv    = (sigma2j/sigma2xi)./(c/sigma2xi + 1/sigma2j);
        p1Tjv    = V.*repmat(c_inv',q,1);
        Tjv      = p1Tjv*p2Tjv; % DSTF
        Tvj      = Tjv'; % Tranconjugated DSTF
        p1ThetaS = Thetajj*Tjv;
        p2ThetaS = Svv*Tvj;
        ThetaS   = p1ThetaS*p2ThetaS;
        Txiv     = (Ip - Lvj*Tjv); % DRTF
        Tvxi     = Txiv'; % Tranconjugated DRTF
        Sxixi    = Txiv*Svv*Tvxi; % Residual Empirical Covariance (REC)
        Sxixi    = (Sxixi + Sxixi')/2;
        % Likelihood
        llh_pst  = sum(log(c_inv)) + q*log(sigma2j) + sum(log(dinv));
        % llh_pst           = 0;
        llh_xixi = p*log(1/sigma2xi) - (1/sigma2xi)*sum(abs(diag(Sxixi))) - p*axi/sigma2xi;
        if penalty == 0 % naive
            llh_jj  = q*log(1/sigma2j) + sum(log(d)) - (1/sigma2j)*sum(abs(diag(ThetaS)));
        elseif penalty == 1 % lasso
            llh_jj  = q*log(1/sigma2j) + sum(log(d)) - (1/sigma2j)*sum(abs(diag(ThetaS))) - aj*(1/sigma2j)*llh_norm1;
        elseif penalty == 2 % ridge
            llh_jj  = q*log(1/sigma2j) + sum(log(d)) - (1/sigma2j)*sum(abs(diag(ThetaS))) - (aj^2/2)*(1/sigma2j)^2*llh_norm2;
        end
        llh(count_sigma2xi,count_sigma2j) = llh_pst + llh_xixi + llh_jj;
        fprintf(1,'\b\b\b\b%3.0f%%',(count/(length(grid_sigma2xi)*length(grid_sigma2j)))*100);
        count = count + 1;
    end
end
fprintf(1,'\n');
end