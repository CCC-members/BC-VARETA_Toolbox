function [sigma2j,sigma2j_post,Tjv,Svv,LvjW,scaleJ,scaleLvj,gamma_grid,gcv] = lcmv(Svv,Lvj,param)

%% Standardization
W             = param.W;
Winv          = param.Winv;
LvjW          = Lvj*W;
WLjv          = LvjW';
p             = size(LvjW,1);
Ip            = eye(p);
q             = size(LvjW,2);
sigma2j_post  = zeros(q,1);
gamma1        = param.gamma1;
gamma2        = param.gamma2;
delta_gamma   = param.delta_gamma;
gamma_grid    = gamma1:delta_gamma:gamma2;
gcv           = zeros(length(gamma_grid),1);
count         = 1;
IsField       = param.field;
flag          = param.flag;
run_bash_mode = param.run_bash_mode;

%%
scaleLvj      = sqrt(sum(abs(diag(LvjW*WLjv)))/p);
LvjW          = LvjW/scaleLvj;
LvjW3D         = LvjW;
scaleJ        = (sum(abs(diag(Svv)))/p);
Svv           = Svv/scaleJ;
if(isequal(IsField,3))
    LvjW3D     = reshape(LvjW3D,p,3,q/3);
    LvjW3D     = permute(LvjW3D,[1 3 2]);
end
%%
if(~run_bash_mode)
    process_waitbar = waitbar(0,'Please wait...','windowstyle', 'modal');
    frames = java.awt.Frame.getFrames();
    frames(end).setAlwaysOnTop(1)
end
disp(flag);
fprintf(1,strcat("-->> Running LCMV (",param.str_band,') process: %3d%%\n'),0);
for gamma = gamma_grid        
    [Tvj,~,Wout]       = mkfilt_lcmv(LvjW3D,10^gamma);
    if(isequal(IsField,3))
        Tvj   = permute(Tvj,[1 3 2]);
        Tvj   = reshape(Tvj,p,q);
    end
    Tjv              = Tvj';
    Txiv             = Ip - LvjW*Tjv;
    gcv(count)       = (1/p)*sum(abs(diag(Txiv*Svv*Txiv')))/((1/p)*sum(abs(diag(Txiv))))^2;
    fprintf(1,'\b\b\b\b%3.0f%%',(count)/(length(gamma_grid))*100-1);
    if(~run_bash_mode)
        text = replace(param.str_band,'_','-');
        waitbar((count)/(length(gamma_grid)),process_waitbar,...
            strcat("Running LCMV (",text,") process: ",num2str(fix((count)/(length(gamma_grid))*100)-1),"%"));
    end
    count            = count + 1;
end
clearvars Txiv;
%%
[gcv_opt,idx_gamma] = min(gcv);
gamma               = gamma_grid(idx_gamma);
[Tvj,~,Wout]        = mkfilt_lcmv(LvjW3D,10^gamma);
clearvars LvjW3D;
if(isequal(IsField,3))
    Tvj   = permute(Tvj,[1 3 2]);
    Tvj   = reshape(Tvj,p,q);   
end
Tvj               = Tvj*Winv;
clearvars Winv;
Tjv               = Tvj';
Sjj               = Tjv*Svv*Tvj;
clearvars Tvj;
Sjj               = (Sjj + Sjj')/2;
sigma2j           = diag(Sjj);
clearvars Sjj;
Wout              = reshape(Wout,1,1,length(Wout));
if(isequal(IsField,3))
    Wout3D        = zeros(3,3,length(Wout));
    for i=1:3
       Wout3D(i,i,:) =  Wout;
    end
    Wout = Wout3D;
    clearvars Wout3D;
end
Wout              = num2cell(Wout,[1,2]);
%% Destandardization using a final pass solution
sigma2jW          = blkdiag(Wout{:})*W';
sigma2jWLjv       = blkdiag(Wout{:})*WLjv;
clearvars WLjv Wout;
sigma2j_post0     = (W*sigma2jWLjv)/(LvjW*sigma2jWLjv+10^gamma*Ip);
LvjWsigma2j       = sigma2jWLjv';
clearvars sigma2jWLjv;
LvjWsigma2jW      = LvjWsigma2j*W';
clearvars LvjWsigma2j;
% Only save the diagonals of the Posterior Covariance
for count_gen = 1:length(Lvj)
    sigma2j_post(count_gen) = W(count_gen,:)*sigma2jW(:,count_gen) - sigma2j_post0(count_gen,:)*LvjWsigma2jW(:,count_gen);
end
if(~run_bash_mode)
    text = replace(param.str_band,'_','-');
    waitbar(1,process_waitbar,strcat("Running LCMV (",text,") process: 100%"));
end
fprintf(1,'\b\b\b\b%3.0f%%',100);
fprintf(1,'\n');
if(~run_bash_mode && exist('process_waitbar','var'))
    delete(process_waitbar);
end

end