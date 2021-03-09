function [llh] = higgs_likelihood(Sxixi,sigma2xi,Sjj,Thetajj,Sigmajj_post,param)
axi               = param.axi;
aj                = param.aj;
Ajj               = param.Ajj;
penalty           = param.penalty;
%%
if(~param.use_gpu)
    [~,D]             = eig(Sigmajj_post);
    d_post            = diag(D);
    [~,D]             = eig(Thetajj);
    d                 = diag(D);
else
    gpuSigmajj_pst    = gpuArray(Sigmajj_post);
    [~,D]             = eig(gpuSigmajj_pst);
    clear gpuSigmajj_pst
    d_post            = diag(D);
    d_post            = gather(d_post);
    gpuThetajj        = gpuArray(Thetajj);
    [~,D]             = eig(gpuThetajj);
    clear gpuThetajj
    d                 = diag(D);
    d                 = gather(d);
end
llh_pst           = sum(log(d_post));
% llh_pst           = 0;
llh_xixi          = sum(log(1./sigma2xi)) - sum(abs(diag(diag(1./sigma2xi)*Sxixi))) - axi*sum(1./sigma2xi);
if penalty == 0 % naive
    llh_jj                = sum(log(d)) - sum(abs(sum(Thetajj.*transpose(Sjj),2)));
elseif penalty == 1 % lasso
    llh_jj                = sum(log(d)) - sum(abs(sum(Thetajj.*transpose(Sjj),2))) - aj*sum(abs(Ajj(:).*Thetajj(:)));
elseif penalty == 2 % ridge
    llh_jj                = sum(log(d)) - sum(abs(sum(Thetajj.*transpose(Sjj),2))) - (aj^2/2)*sum(abs(diag(Thetajj*Thetajj')));
end
llh               = llh_pst + llh_xixi + llh_jj;
end