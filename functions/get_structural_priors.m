function subject = get_structural_priors(subject)

%% Getting Homogeneus field
disp("BC-V-->> Structural information")
disp('---------------------------------------------------------------------');
Cortex      = subject.Scortex;
Channel     = subject.Cdata.Channel;
ChanLoc     = [Channel.Loc]';
Shead       = subject.Shead;
Ng          = length(Cortex.Vertices);
Ns          = length(ChanLoc);
sGain       = zeros(Ns,3,Ng);
if(getGlobalGuimode())
    process_waitbar = waitbar(0,'Please wait...');
end
fprintf(1,strcat('-->> Getting sensor priors process: %3d%%\n'),0);
for s = 1:Ns
    for g = 1:Ng
        sGain(s,:,g) = (1/(4*pi))*(ChanLoc(s,:) - Cortex.Vertices(g,:))/...
            (sum((ChanLoc(s,:) - Cortex.Vertices(g,:)).^2)).^(3/2);
        if(getGlobalGuimode())
            waitbar(((s-1)*Ng+g)/(Ns*Ng)*100,process_waitbar,strcat("Getting sensor priors process: ",num2str(fix(((s-1)*Ng+g)/(Ns*Ng)*100)),"%"));
        end
        fprintf(1,'\b\b\b\b%3.0f%%',((s-1)*Ng+g)/(Ns*Ng)*100);
    end
end
subject.Cdata.sGain = sGain;
[U,D,V]             = svd(reshape(sGain,Ns,3*Ng),'econ');
sGain_pinv          = V*diag(diag(D).^(-1))*U';
sGain_pinv          = reshape(sGain_pinv,3,Ng,Ns);
fprintf(1,'\n');
if(getGlobalGuimode())
    delete(process_waitbar);
end

Nh = length(Shead.Vertices);
hsGain = zeros(Nh,Ns);
hGain = zeros(3,Ng);
if(getGlobalGuimode())
    process_waitbar = waitbar(0,'Please wait...');
end
fprintf(1,strcat('-->> Getting structural priors process: %3d%%\n'),0);
for h = 1:Nh
    for g = 1:Ng
        hGain(:,g) = (1/(4*pi))*(Shead.Vertices(h,:) - Cortex.Vertices(g,:))/...
            (sum((Shead.Vertices(h,:) - Cortex.Vertices(g,:)).^2)).^(3/2);
    end
    for s = 1:Ns
        hsGain(h,s) = sum(hGain.*sGain_pinv(:,:,s),"all");
    end
    if(getGlobalGuimode())
        waitbar((h)/(Nh)*100,process_waitbar,strcat("Getting structural priors process: ",num2str(fix((h)/(Nh)*100)),"%"));
    end
    fprintf(1,'\b\b\b\b%3.0f%%',(h)/(Nh)*100);
end

subject.Shead.hsGain = hsGain;
fprintf(1,'\n');
if(getGlobalGuimode())
    delete(process_waitbar);
end

end

