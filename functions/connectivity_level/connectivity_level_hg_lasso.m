function [subject,properties,outputs] = connectivity_level_hg_lasso(subject,properties)

% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019

% Updates
% - Ariosky Areces Gonzalez

% Date: March 20, 2019

GridOrient          = subject.GridOrient;
GridAtlas           = subject.GridAtlas;
Nseg                = properties.sensor_level_out.Nseg;
peak_pos            = properties.sensor_level_out.peak_pos;
band                = properties.sensor_level_out.band;
Svv                 = properties.sensor_level_out.Svv;
indms               = properties.activation_level_out.indms;
Tjv                 = properties.activation_level_out.T;
scaleSvv            = properties.activation_level_out.scaleSvv;
scaleKe             = properties.activation_level_out.scaleKe;
activation_params   = properties.activation_params;
connectivity_params = properties.connectivity_params;
IsCurv              = activation_params.IsCurv.value;  % 0 (no compensation) 1 (giri and sulci curvature compensation)
IsField             = activation_params.IsField.value; % 1 (projected Lead Field) 3 (3D Lead Field)

%%
%% HG-LASSO parameters
%%
param.use_gpu         = properties.general_params.use_gpu.value;
m                     = length(peak_pos)*Nseg;
param.m               = m;
param.nu              = m;
p                     = length(Svv);
param.p               = p;
param.Ip              = eye(p);
param.Op              = ones(p,1);
param.Axixi           = eye(p);
param.Axixi_inv       = eye(p);
param.axi             = connectivity_params.axi.value;
param.maxiter_outer   = connectivity_params.maxiter_outer.value;
param.maxiter_inner   = connectivity_params.maxiter_inner.value;
param.ntry            = connectivity_params.ntry.value;
param.prew            = connectivity_params.prew.value;
param.penalty         = connectivity_params.penalty.value;
param.rth1            = connectivity_params.rth1.value;
param.rth2            = connectivity_params.rth2.value;
param.eigreg          = connectivity_params.eigreg.value;
param.str_band        = str_band;

%% Connectivity Leakage Module
disp('BC-V-->> Connectivity leakage module...');
q                     = length(indms);
param.q               = q;
param.Iq              = eye(q);
param.Oq              = ones(q,1);
aj                    = sqrt(log(q)/m);
Ajj_diag              = 0;
Ajj_ndiag             = 1;
Ajj                   = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj              = aj;
param.Ajj             = Ajj;

%% Calling two step lasso
if (IsCurv == 0)
    if (IsField == 2) || (IsField == 3)
        Tvj = transpose(Tjv);
        Tvj = bst_gain_orient(Tvj, GridOrient, GridAtlas);
        Tjv = transpose(Tvj);
    end
    Sjj      = Tjv(indms,:)*Svv*Tjv(indms,:)';
    Sjj      = (Sjj + Sjj')/2;
    Sjj      = Sjj*sqrt(scaleSvv)/(scaleKe);
elseif (IsCurv == 1)
    Tjv_giri     = squeeze(Tjv(:,:,1));
    Tjv_sulc     = squeeze(Tjv(:,:,2));
    if (IsField == 2) || (IsField == 3)
        Tjv_giri = transpose(Tjv_giri);
        Tjv_giri = bst_gain_orient(Tjv_giri, GridOrient, GridAtlas);
        Tjv_giri = transpose(Tjv_giri);
        Tjv_sulc = transpose(Tjv_sulc);
        Tjv_sulc = bst_gain_orient(Tjv_sulc, GridOrient, GridAtlas);
        Tjv_sulc = transpose(Tjv_sulc);
    end
    Sjj_giri = Tjv_giri(indms,:)*Svv*Tjv_giri(indms,:)';
    Sjj_sulc = Tjv_sulc(indms,:)*Svv*Tjv_sulc(indms,:)';
    Sjj      = (Sjj_giri + Sjj_sulc)/2;
    Sjj      = (Sjj + Sjj')/2;
    Sjj      = Sjj*sqrt(scaleSvv(1)*scaleSvv(2))/(scaleKe(1)*scaleKe(2));
end
[Thetajj,Sigmajj]     = twostep_lasso_caller(Sjj,param);

outputs.Thetajj = Thetajj;
outputs.Sjj = Sjj;
outputs.Sigmajj = Sigmajj;
end

