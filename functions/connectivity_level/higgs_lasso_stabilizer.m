function [rth] = higgs_lasso_stabilizer(Psijj,Thetajj,Thetajj_unb,Thetajj_var,Svv,Lvj,sigma2xi0,Sigmajj0,llh0,param)
%% try likelihood evolution for the mask given by all thresholds
aj        = param.aj;
rth1      = param.rth1;
rth2      = param.rth2;
m         = param.m;
ntry      = param.ntry;
rth_grid  = rth1:0.1:rth2;
run_bash_mode  = param.run_bash_mode;
if(~run_bash_mode)
    process_waitbar = waitbar(0,'Please wait...');
end
fprintf(1,'-->> Computing optimal Ryleigh threshold: %3d%%\n',0);
if ntry == 0
    %% check only partial likelihood
    llhjj_grid       = zeros(length(rth_grid),1);
    for rth_count = 1:length(rth_grid)
        %% Zeros mask
        rth                   = rth_grid(rth_count);
        mask                  = find(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var))));
        Thetajj_mask          = Thetajj.X;
        Thetajj_mask(mask)    = 0;
        [Thetajj_mask]        = higgs_eigendecomposition(Thetajj_mask,param);
        llhjj_grid(rth_count) = sum(log(abs(Thetajj_mask.d))) - sum(abs(sum(Thetajj_mask.X.*transpose(Psijj.X),2))) - aj*sum(abs(Thetajj_mask.X(:)));
        fprintf(1,'\b\b\b\b%3.0f%%',(rth_count/length(rth_grid))*100 - 1);
        if(~run_bash_mode)
            waitbar(rth_count/length(rth_grid),process_waitbar,strcat("Computing optimal Ryleigh threshold: ",num2str(fix((rth_count/length(rth_grid))*100)-1),"%"));
        end
    end
    [val,id_rth]              = max(llhjj_grid);
    id_rth                    = min(id_rth);
    rth                       = rth_grid(id_rth);
else
    %% check global likelihood
    llh_grid  = zeros(length(rth_grid),ntry);
    for rth_count = 1:length(rth_grid)
        %% Zeros mask
        rth                     = rth_grid(rth_count);
        mask                    = find(abs(Thetajj_unb) < (rth/sqrt(m))*(Thetajj_var - diag(diag(Thetajj_var))));
        %% Perform ntry iterations
        [llh_grid(rth_count,:)] = higgs_try_likelihood(mask,Svv,Lvj,sigma2xi0,Sigmajj0,llh0,param);
        fprintf(1,'\b\b\b\b%3.0f%%',(rth_count/length(rth_grid))*100 - 1);
        if(~run_bash_mode)
            waitbar(rth_count/length(rth_grid),process_waitbar,strcat("Computing optimal Ryleigh threshold: ",num2str(fix((rth_count/length(rth_grid))*100)-1),"%"));
        end
    end
    [rth]                       = higgs_check_likelihood_trend(llh0,llh_grid,rth_grid,ntry);
end
fprintf(1,'\b\b\b\b%3.0f%%',100);
fprintf(1,'\n');
if(~run_bash_mode)
    waitbar(1,process_waitbar,strcat("Computing optimal Ryleigh threshold: ",num2str(100),"%"));
    delete(process_waitbar)
end
end