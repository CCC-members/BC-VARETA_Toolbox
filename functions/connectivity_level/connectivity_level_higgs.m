function [subject,properties,outputs] = connectivity_level_higgs(subject,properties)

% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019


% Updates
% - Ariosky Areces Gonzalez

% Date: Febrery 2, 2021

%%
%% Preparing params
%%
Ke                      = subject.Ke;
W                       = subject.W;
Winv                    = subject.Winv;
sub_to_FSAve            = subject.sub_to_FSAve;

%%
%% Sensor level Outputs
%%
sensor_level_out        = subject.sensor_level_out;
Svv                     = sensor_level_out.Svv;
Nseg                    = sensor_level_out.Nseg;

%%
%% Activation level Outputs
%%
activation_level_out    = subject.activation_level_out;
actv_method             = activation_level_out.method;
actv_methods            = properties.activation_params.methods;
for i=1:length(actv_methods)
    if(isequal(actv_methods{i}.method,actv_method))
        actv_th         = actv_methods{i}.threshold.value;            
    end
end

%%
%% HIGGS parameters
%%
connectivity_params     = properties.connectivity_params;
higgs_th                = connectivity_params.higgs_th.value;
IsCurv                  = connectivity_params.IsCurv.value; % 0 (no compensation) 1 (giri and sulci curvature compensation)
IsField                 = connectivity_params.IsField.value; % 1 (projected Lead Field) 3 (3D Lead Field)

param.use_gpu           = properties.general_params.use_gpu.value;
m                       = Nseg;
param.m                 = m;
param.nu                = m;
p                       = length(Svv);
param.p                 = p;
param.Ip                = eye(p);
param.Op                = ones(p,1);
param.Axixi             = eye(p);
param.Axixi_inv         = eye(p);
param.axi               = connectivity_params.axi.value;
param.maxiter_outer     = connectivity_params.maxiter_outer.value;
param.maxiter_inner     = connectivity_params.maxiter_inner.value;
param.ntry              = connectivity_params.ntry.value;
param.prew              = connectivity_params.prew.value;
param.penalty           = connectivity_params.penalty.value;
param.rth1              = connectivity_params.rth1.value;
param.rth2              = connectivity_params.rth2.value;
param.eigreg            = connectivity_params.eigreg.value;
param.str_band          = str_band;
param.W                 = W;
param.parcellation      = subject.parcellation;

%%
%% Checking activation and Connectivity Threshold
%%
disp('BC-V-->> Connectivity leakage module...');
if(isequal(actv_th,higgs_th))
    indms               = activation_level_out.indms;
else
    if IsCurv == 0
        s2j             = activation_level_out.s2j;
        sigma2j         = activation_level_out.sigma2j_post;
        if IsField == 2 || IsField == 3
            s2j         = sum(reshape(abs(s2j),3,length(Ke)/3),1)';
            sigma2j     = sum(reshape(abs(sigma2j),3,length(Ke)/3),1)';
        end
        stat            = sqrt(s2j./sigma2j);
        indms           = find(stat > higgs_th);
    elseif IsCurv == 1
        s2j             = activation_level_out.s2j;
        s2j_giri        = s2j(:,1);
        s2j_sulc        = s2j(:,2);
        sigma2j         = activation_level_out.sigma2j;
        sigma2j_giri    = sigma2j(:,1);
        sigma2j_sulc    = sigma2j(:,2);
        if IsField == 2 || IsField == 3
            s2j_giri    = sum(reshape(abs(s2j_giri),3,length(Ke)/3),1)';
            sigma2j_giri= sum(reshape(abs(sigma2j_giri),3,length(Ke)/3),1)';
            s2j_sulc    = sum(reshape(abs(s2j_sulc),3,length(Ke)/3),1)';
            sigma2j_sulc= sum(reshape(abs(sigma2j_sulc),3,length(Ke)/3),1)';
        end
        stat_giri       = sqrt(s2j_giri./sigma2j_giri);
        indms_giri      = find(stat_giri > higgs_th);
        stat_sulc       = sqrt(s2j_sulc./sigma2j_sulc);
        indms_sulc      = find(stat_sulc > higgs_th);        
        indms           = unique([indms_giri;indms_sulc]);
        clearvars sigma2j_giri sigma2j_sulc;
        clearvars stat_giri stat_sulc;
    end
end

%%
%% Connectivity Leakage Module
%%
q                       = length(indms);
param.q                 = q;
param.Iq                = eye(q);
param.Oq                = ones(q,1);
aj                      = sqrt(log(q)/m);
Ajj_diag                = 0;
Ajj_ndiag               = 1;
Ajj                     = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj                = aj;
param.Ajj               = Ajj;
if IsCurv == 0
    Ke                  = Ke*W;
    [Thetajj,Tjv,llh]   = higgs(Svv,Ke(:,indms),param);
    [Thetajj,s2j,Tjv]   = higgs_destandardization(Thetajj,Svv,Tjv,Winv,W,indms,IsField);
elseif IsCurv == 1
    Ke_giri             = subject.Ke_giri;
    Ke_sulc             = subject.Ke_sulc;
    Ke_giri             = Ke_giri*W;
    Ke_sulc             = Ke_sulc*W;
    [Thetajj_sulc,Tjv_sulc,llh_sulc,Sjj_sulc,Psijj_sulc,Sigmajj_sulc] = higgs(Svv,Ke_sulc(:,indms),param);
    [Thetajj_giri,Tjv_giri,llh_giri,Sjj_giri,Psijj_giri,Sigmajj_giri] = higgs(Svv,Ke_giri(:,indms),param);
    [Thetajj_giri,s2j_giri,Tjv_giri] = higgs_destandardization(Thetajj_giri,Svv,Tjv_giri,Winv,W,indms,IsField);
    [Thetajj_sulc,s2j_sulc,Tjv_sulc] = higgs_destandardization(Thetajj_sulc,Svv,Tjv_sulc,Winv,W,indms,IsField);
    Thetajj             = (Thetajj_giri + Thetajj_sulc)/2;
    Sjj                 = (Sjj_giri + Sjj_sulc)/2;
    Psijj               = (Psijj_giri + Psijj_sulc)/2;
    Sigmajj             = (Sigmajj_giri + Sigmajj_sulc)/2;
    clearvars Thetajj_sulc Thetajj_giri Ke_sulc Ke_giri Sjj_giri Sjj_sulc Psijj_giri Psijj_sulc Sigmajj_giri Sigmajj_sulc;
    s2j                 = (s2j_giri + s2j_sulc)/2;
    llh                 = [llh_giri llh_sulc];
    Tjv                 = cat(3,Tjv_giri,Tjv_sulc);
    clearvars Tjv_giri Tjv_sulc s2j_giri s2j_sulc;
end
clearvars param W Winv;

%% Ordering results by FSAve indices
disp("-->> Ordering connectivity results by FSAve indices");
Msub_to_FSAve   = zeros(length(sub_to_FSAve),size(Ke,2));
for h = 1:length(sub_to_FSAve)
    indices = sub_to_FSAve(h,:);
    Msub_to_FSAve(h,[indices(1) indices(2) indices(3)]) = 1/3;
end

Msub_to_FSAve           = sparse(Msub_to_FSAve);
Thetajj_FSAve           = zeros(size(Ke,2));
Thetajj_FSAve(indms,indms) = Thetajj;
Thetajj_FSAve           = sparse(Thetajj_FSAve);
Thetajj_FSAve           = Msub_to_FSAve*Thetajj_FSAve;
Thetajj_FSAve           = Thetajj_FSAve*Thetajj_FSAve';
Thetajj_FSAve           = full(Thetajj_FSAve);
indms_FSAve             = find(diag(Thetajj_FSAve) > 0);
Thetajj_FSAve           = Thetajj_FSAve(indms_FSAve,indms_FSAve);

%% Outputs variables
outputs.Thetajj         = Thetajj;
outputs.s2j             = s2j;
outputs.Tjv             = Tjv;
outputs.llh             = llh;
outputs.Svv             = Svv;
outputs.Thetajj_FSAve   = Thetajj_FSAve;
outputs.indms_FSAve     = indms_FSAve;
outputs.Sjj             = Sjj;
outputs.Psijj           = Psijj;
outputs.Sigmajj         = Sigmajj;
end

