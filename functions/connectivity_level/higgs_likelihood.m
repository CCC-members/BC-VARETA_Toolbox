function [llh] = higgs_likelihood(Sxixi,theta2xi,Sjj,Thetajj,Sigmajj_post,param)
axi                   = param.axi;
aj                    = param.aj;
Ajj                   = param.Ajj;
penalty               = param.penalty;
%%
llh_xixi              = sum(log(theta2xi)) - abs(sum(theta2xi.*diag(Sxixi))) - axi*sum(theta2xi);
llh_jj_pst            = sum(log(Sigmajj_post.d));
if penalty == 0 % naive
    llh_jj            = sum(log(Thetajj.d)) - sum(abs(sum(Thetajj.X.*transpose(Sjj),2)));
elseif penalty == 1 % lasso
    llh_jj            = sum(log(Thetajj.d)) - sum(abs(sum(Thetajj.X.*transpose(Sjj),2))) - aj*sum(abs(Ajj(:).*Thetajj.X(:)));
elseif penalty == 2 % ridge
    llh_jj            = sum(log(Thetajj.d)) - sum(abs(sum(Thetajj.X.*transpose(Sjj),2))) - (aj^2/2)*sum(abs(sum(Thetajj.X.*Thetajj.X,2)));
end
llh                   = llh_xixi + llh_jj_pst + llh_jj;
end