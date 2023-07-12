function [subject,properties,outputs] = activation_level_eloreta(subject,properties)

% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019


% Updates
% - Ariosky Areces Gonzalez

% Date: Febrery 2, 2021

%%
%% BC-VARETA Toolbox...
%% Preparing params
Lvj           = subject.Ke;
W             = subject.W;
sub_to_FSAve  = subject.sub_to_FSAve;

%%
%% Sensor level Outputs
%%
sensor_level_out    = subject.sensor_level_out;
Svv                 = sensor_level_out.Svv;
band                = sensor_level_out.band;
str_band            = band.str_band;

%%
%% eLORETA activation parameters
%%
activation_params   = properties.activation_params.methods{2};
gamma1              = activation_params.gamma1.value;
gamma2              = activation_params.gamma2.value;
delta_gamma         = activation_params.delta_gamma.value;
threshold           = activation_params.threshold.value;
IsCurv              = activation_params.IsCurv.value; % 0 (no compensation) 1 (giri and sulci curvature compensation)
IsField             = activation_params.IsField.value; % 1 (projected Lead Field) 3 (3D Lead Field)

%%
%% Activation Leakage Module spectral eLORETA
%%
flag = "-->> Running source activation level.";
disp('BC-V-->> eLORETA activation leakage module.');

%% eLORETA params
param.field         = IsField;
param.gamma1        = gamma1;
param.gamma2        = gamma2;
param.delta_gamma   = delta_gamma;
param.flag          = flag;
param.str_band      = str_band;
param.W             = W;
param.Winv          = subject.Winv;

%%
if IsCurv == 0    
    [s2j,sigma2j_post,T,~,~,scaleSvv,scaleKe] = eloreta(Svv,Lvj,param);
    clearvars param Svv;
    if IsField == 2 || IsField == 3
        s2j               = sqrt(sum(reshape(abs(s2j),3,length(Lvj)/3),1))';
        sigma2j_post      = sqrt(sum(reshape(abs(sigma2j_post),3,length(Lvj)/3),1))';
    end
    clearvars Ke;
    stat                  = s2j./sigma2j_post;
    indms                 = find(stat > threshold);
    J                     = s2j;
    J                     = J*scaleSvv/scaleKe^2;
    Jsp                   = zeros(length(stat),1);
    Jsp(indms)            = J(indms);
elseif IsCurv == 1
    param.flag = strcat(flag," Giri compensation");
    [s2j_giri,sigma2j_post_giri,Tgiri,~,~,scaleSvv_giri,scaleKe_giri] = eloreta(Svv,subject.Ke_giri,param);
    param.flag = strcat(flag," Sulci compensation");
    [s2j_sulc,sigma2j_post_sulc,Tsulc,~,~,scaleSvv_sulc,scaleKe_sulc] = eloreta(Svv,subject.Ke_sulc,param);
    clearvars param Svv;
    disp("-->> Applying giri and sulci compensation.");
    if IsField == 2 || IsField == 3
        s2j_giri               = sqrt(sum(reshape(abs(s2j_giri),3,length(Lvj)/3),1))';
        sigma2j_post_giri      = sqrt(sum(reshape(abs(sigma2j_post_giri),3,length(Lvj)/3),1))';
        s2j_sulc               = sqrt(sum(reshape(abs(s2j_sulc),3,length(Lvj)/3),1))';
        sigma2j_post_sulc      = sqrt(sum(reshape(abs(sigma2j_post_sulc),3,length(Lvj)/3),1))';
    end
    clearvars Ke;
    stat_giri             = s2j_giri./sigma2j_post_giri;
    indms_giri            = find(stat_giri > threshold);
    stat_sulc             = s2j_sulc./sigma2j_post_sulc;
    indms_sulc            = find(stat_sulc > threshold);    
    s2j                   = [s2j_giri s2j_sulc];
    sigma2j_post          = [sigma2j_post_giri sigma2j_post_sulc];
    clearvars sigma2j_post_giri sigma2j_post_sulc;
    scaleSvv              = [scaleSvv_giri scaleSvv_sulc];
    scaleKe               = [scaleKe_giri scaleKe_sulc];
    stat                  = [stat_giri stat_sulc];
    clearvars stat_giri stat_sulc;
    indms                 = unique([indms_giri;indms_sulc]);
    clearvars indms_giri indms_sulc;
    J                     = (s2j_giri + s2j_sulc)/2;
    clearvars s2j_giri s2j_sulc;
    J                     = J*sqrt(scaleSvv_giri*scaleSvv_sulc)/(scaleKe_giri*scaleKe_sulc);
    clearvars scaleSvv_giri scaleSvv_sulc;
    clearvars scaleKe_giri scaleKe_sulc;
    Jsp                   = zeros(length(stat),1);
    Jsp(indms)            = J(indms);
    T                     = cat(3,Tgiri,Tsulc);
    clearvars Tgiri Tsulc;
end
% Ordering results by FSAve indices
J_FSAve   = zeros(length(J),1);
Jsp_FSAve = zeros(length(J),1);
for h=1:length(sub_to_FSAve)
    indices           = sub_to_FSAve(h,:);    
    J_FSAve(h)        = (J(indices(1))+J(indices(2))+J(indices(3)))/3;
    Jsp_FSAve(h)      = (Jsp(indices(1))+Jsp(indices(2))+Jsp(indices(3)))/3; 
end
outputs.s2j = s2j;
outputs.sigma2j_post = sigma2j_post;
outputs.T = T;
outputs.scaleSvv = scaleSvv;
outputs.scaleKe = scaleKe;
outputs.stat = stat;
outputs.J = J;
outputs.Jsp = Jsp;
outputs.indms = indms;
outputs.J_FSAve = J_FSAve;
outputs.Jsp_FSAve = Jsp_FSAve;
end

