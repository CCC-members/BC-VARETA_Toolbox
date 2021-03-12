function llh = higgs_crossvalidation(grid_sigma2xi,grid_sigma2j,Sigmajj,Thetajj,Svv,Lvj,param)
%% Parameters
p              = param.p;
Ip             = param.Ip;
q              = param.q;
Iq             = param.Iq;
axi            = param.axi;
aj             = param.aj;
Ajj            = param.Ajj;
penalty        = param.penalty;
eigreg         = param.eigreg;
%% Eigenvalues for likelihood
disp("-->> Running prewarming, preparing constants for crossvalidation");
SigmajjLjv     = Sigmajj.X*Lvj';
SigmajjLjj     = SigmajjLjv*Lvj;
if(~param.use_gpu)
    [U,D]      = eig(SigmajjLjj); % Ljj = Thetajj*U*D*Uinv
    Uinv       = Iq/U;
    d          = diag(D);
    clear SigmajjLjj D
else
    gpuSigmajjLjj = gpuArray(SigmajjLjj);
    [U,D]      = eig(gpuSigmajjLjj); % Ljj = Thetajj*U*D*Uinv
    Uinv       = Iq/U;
    d          = diag(D);
    [U,d,Uinv] = gather(U,d,Uinv);
    clear gpuSigmajjLjj SigmajjLjj D
end
dmax          = max(d);
dmin          = min(d);
if (dmin == 0 || dmin/dmax < eigreg)
    d         = (d + eigreg*dmax)/(1 + eigreg);
end
p2Tjv         = Uinv*SigmajjLjv;
llh_norm1     = sum(abs(Ajj(:).*Thetajj.X(:)));
llh_norm2     = abs(sum(sum(Thetajj.X.*transpose(Thetajj.X)),2));
%% Cross-validation cycle
llh = zeros(length(grid_sigma2xi),length(grid_sigma2j));
fprintf(1,'-->> Running prewarming, crossvalidating higgs initial values: %3d%%\n',0);
count = 1;
for count_sigma2xi = 1:length(grid_sigma2xi)
    sigma2xi = 10^(grid_sigma2xi(count_sigma2xi));
    for count_sigma2j = 1:length(grid_sigma2j)
        sigma2j = 10^(grid_sigma2j(count_sigma2j));
        % Expectation
        d_inv        = 1./(d/sigma2xi + 1/sigma2j);
        p1Tjv        = U*spdiags(d_inv/sigma2xi,0,q,q);
        Tjv          = p1Tjv*p2Tjv; % DSTF
        p1ThetaS     = Thetajj.X*Tjv;
        p2ThetaS     = Svv*Tjv';
        Txiv         = (Ip - Lvj*Tjv); % DRTF
        Tvxi         = Txiv'; % Tranconjugated DRTF
        Sxixi        = Txiv*Svv*Tvxi; % Residual Empirical Covariance (REC)
        Sxixi        = (Sxixi + Sxixi')/2;
        % Likelihood
        llh_xixi     = p*log(1/sigma2xi) - (1/sigma2xi)*sum(abs(diag(Sxixi))) - p*axi/sigma2xi;
        llh_jj_pst   = sum(log(d_inv)) + sum(log(Sigmajj.d));
        if penalty == 0 % naive
            llh_jj   = q*log(1/sigma2j) + sum(log(Thetajj.d)) - (1/sigma2j)*sum(abs(sum(p1ThetaS.*transpose(p2ThetaS),2)));
        elseif penalty == 1 % lasso
            llh_jj   = q*log(1/sigma2j) + sum(log(Thetajj.d)) - (1/sigma2j)*sum(abs(sum(p1ThetaS.*transpose(p2ThetaS),2))) - aj*(1/sigma2j)*llh_norm1;
        elseif penalty == 2 % ridge
            llh_jj   = q*log(1/sigma2j) + sum(log(Thetajj.d)) - (1/sigma2j)*sum(abs(sum(p1ThetaS.*transpose(p2ThetaS),2))) - (aj^2/2)*(1/sigma2j)^2*llh_norm2;
        end
        llh(count_sigma2xi,count_sigma2j) = llh_xixi + llh_jj_pst + llh_jj;
        fprintf(1,'\b\b\b\b%3.0f%%',(count/(length(grid_sigma2xi)*length(grid_sigma2j)))*100);
        count = count + 1;
    end
end
fprintf(1,'\n');
end