function [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_post,Tjv] = higgs_expectation(Svv,Lvj,sigma2xi,Sigmajj,param)
Ip                = param.Ip;

%% Source Posterior Covariance (SPC)
SigmajjLjv        = Sigmajj.X*Lvj'; % SEC*SDTF'
Sigmajj_post      = higgs_eigendecomposition(Sigmajj.X - SigmajjLjv*(Ip/(Lvj*SigmajjLjv+diag(sigma2xi)))*SigmajjLjv',param); % SPC by Woodbury formula

%% Data to Source Transfer Function (DSTF)
Tjv               = Sigmajj_post.X*Lvj'*diag(1./sigma2xi); % DSTF

%% Source Empirical Covariance (ESEC)
Sjj               = Tjv*Svv*Tjv'; % Source Empirical Covariance (SEC)

%% Effective Source Empirical Covariance (ESEC)
Psijj             = higgs_eigendecomposition(Sigmajj_post.X + Sjj,param); % ESEC    

%% Data to Residuals Transfer Function (DRTF)
Txiv              = (Ip - Lvj*Tjv); % DRTF

%% Residuals Posterior Covariance (SPC)
Sigmaxixi_post    =  Lvj*Sigmajj_post.X*Lvj';

%% Residuals Empirical Covariance (REC)
Sxixi             = Txiv*Svv*Txiv'; % Residual Empirical Covariance (REC)

%% Effective Residual Empirical Covariance (EREC)
Psixixi           = Sigmaxixi_post + Sxixi; % EREC

end