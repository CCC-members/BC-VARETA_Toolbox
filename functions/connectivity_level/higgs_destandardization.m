function [Thetajj,s2j,Tjv] = higgs_destandardization(Thetajj,Svv,Tjv,Winv,W,indms,IsField)
disp("-->> Running higgs destandardization.");

if(getGlobalGuimode())
    process_waitbar = waitbar(0,'Running higgs destandardization.','windowstyle', 'modal');
    frames = java.awt.Frame.getFrames();
    frames(end).setAlwaysOnTop(1)
end
Tjv_sp          = sparse(zeros(size(W,2),size(Tjv,2)));
Tjv_sp(indms,:) = Tjv;
Tjv             = W*Tjv_sp;
Tvj             = Tjv';
TjvSvv          = Tjv*Svv;
Tjv             = full(Tjv);
s2j             = zeros(size(W,1),1);
if(getGlobalGuimode())
        waitbar(0.20,process_waitbar,strcat("Running higgs destandardization: 20%"));
end
if IsField == 2
    indms3D      = [3*indms-2 3*indms-1 3*indms];
    indms3D      = indms3D';
    indms3D      = indms3D(:);
    for gen = indms3D'
        s2j(gen) = TjvSvv(gen,:)*Tvj(:,gen);
    end
    s2j          = sum(reshape(s2j,3,size(W,2)),1)';
else
    for gen = indms'
        s2j(gen) = TjvSvv(gen,:)*Tvj(:,gen);
    end
end
if(getGlobalGuimode())
        waitbar(0.40,process_waitbar,strcat("Running higgs destandardization: 40%"));
end
Thetajj_sp                   = sparse(zeros(size(W,2)));
Thetajj_sp(indms,indms)      = Thetajj;
Thetajj                      = Winv'*Thetajj_sp;
Thetajj                      = Thetajj*Winv;
Thetajj                      = full(Thetajj);
if(getGlobalGuimode())
        waitbar(0.75,process_waitbar,strcat("Running higgs destandardization: 75%"));
end
if IsField == 2
    theta2j                  = diag(Thetajj);
    Thetajj                  = Thetajj - diag(theta2j);    
    Thetajj                  = sum(reshape(Thetajj,3*size(W,2),3,size(W,2)),2);
    Thetajj                  = squeeze(Thetajj);
    Thetajj                  = sum(reshape(Thetajj,3,size(W,2),size(W,2)),1);
    Thetajj                  = squeeze(Thetajj);
    theta2j                  = sum(reshape(theta2j,3,size(W,2)),1)';
    Thetajj                  = Thetajj - diag(diag(Thetajj)) + diag(theta2j);
end
if(getGlobalGuimode())
        waitbar(0.95,process_waitbar,strcat("Running higgs destandardization: 95%"));
end
Thetajj                      = Thetajj(indms,indms);
Thetajj                      = (Thetajj + Thetajj')/2;
if(getGlobalGuimode() && exist('process_waitbar','var'))
    waitbar(1,process_waitbar,strcat("Running higgs destandardization: 100%"));
    delete(process_waitbar)
end
end
