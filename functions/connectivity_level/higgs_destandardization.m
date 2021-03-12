function [Thetajj,s2j,Tjv] = higgs_destandardization(Thetajj,Svv,Tjv,Winv,W,indms,IsField)
disp("-->> Running higgs destandardization.");
Tjv_sp          = sparse(zeros(size(W,2),size(Tjv,2)));
Tjv_sp(indms,:) = Tjv;
Tjv             = W*Tjv_sp;
Tvj             = Tjv';
TjvSvv          = Tjv*Svv;
Tjv             = full(Tjv);
s2j             = zeros(size(W,1),1);
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
Thetajj_sp                   = sparse(zeros(size(W,2)));
Thetajj_sp(indms,indms)      = Thetajj;
Thetajj                      = Winv'*Thetajj_sp;
Thetajj                      = Thetajj*Winv;
Thetajj                      = full(Thetajj);
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
Thetajj                      = Thetajj(indms,indms);
Thetajj                      = (Thetajj + Thetajj')/2;
end
