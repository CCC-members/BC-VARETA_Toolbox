function [Sxixi,Psixixi,Sjj,Psijj,Sigmajj_post,Tjv] = higgs_expectation(Svv,Lvj,Ljv,sigma2xi,Sigmajj,param)
Ip                = param.Ip;
%% Source Posterior Covariance (SPC)
SigmajjLjv        = Sigmajj*Ljv; % SEC*SDTF'
LvjSigmajj        = SigmajjLjv'; % Tranconjugated SEC*SDTF'
Sigmajj_post      = Sigmajj - SigmajjLjv*(Ip/(Lvj*SigmajjLjv+diag(sigma2xi)))*LvjSigmajj; % SPC by Woodbury formula
Sigmajj_post      = (Sigmajj_post + Sigmajj_post')/2;
clearvars SigmajjLjv LvjSigmajj Sigmajj;
%%
%% Data to Source Transfer Function (DSTF)
Tjv               = Sigmajj_post*Ljv*diag(1./sigma2xi); % DSTF
Tvj               = Tjv'; % Tranconjugated DSTF
%%
%% Source Empirical Covariance (ESEC)
Sjj               = Tjv*Svv*Tvj; % Source Empirical Covariance (SEC)
Sjj               = (Sjj + Sjj')/2;
%%
%% Effective Source Empirical Covariance (ESEC)
Psijj             = Sigmajj_post + Sjj; % ESEC
%%
%% Data to Residuals Transfer Function (DRTF)
Txiv              = (Ip - Lvj*Tjv); % DRTF
Tvxi              = Txiv'; % Tranconjugated DRTF
%%
%% Residuals Posterior Covariance (SPC)
Sigmaxixi_post    =  Lvj*Sigmajj_post*Ljv;
clearvars Lvj Ljv;
Sigmaxixi_post    = (Sigmaxixi_post + Sigmaxixi_post')/2;
%%
%% Residuals Empirical Covariance (REC)
Sxixi             = Txiv*Svv*Tvxi; % Residual Empirical Covariance (REC)
clearvars Txiv Svv Tvxi;
Sxixi             = (Sxixi + Sxixi')/2;
%% Effective Residual Empirical Covariance (EREC)
Psixixi           = Sigmaxixi_post + Sxixi; % EREC
end