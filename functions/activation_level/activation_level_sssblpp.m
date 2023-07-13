function [subject,properties,outputs] = activation_level_sssblpp(subject,properties)
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
Nseg                = sensor_level_out.Nseg;
Svv                 = sensor_level_out.Svv;
band                = sensor_level_out.band;
str_band            = band.str_band;

%%
%% sSSBL++ activation parameters
%%
activation_params   = properties.activation_params;
method              = activation_params.methods{1};
threshold           = method.threshold.value;
IsField             = activation_params.IsField.value; % 1 (projected Lead Field) 3 (3D Lead Field)
IsCurv              = activation_params.IsCurv.value; % 0 (no compensation) 1 (giri and sulci curvature compensation)

%%
%% Activation Leakage Module spectral enet-ssbl
%%
disp('BC-V-->> sSSBL++ activation leakage module.');
flag = "-->> Running source activation level.";
param.Nsamp         = Nseg;
param.str_band      = str_band;
param.W             = W;
param.parcellation  = subject.parcellation;

if IsCurv == 0
    param.flag          = flag;
    [s2j,sigma2j,T,~,scaleSvv,scaleKe] = sSSBLpp(Svv,Lvj,param);
    clearvars param Svv;
    if IsField == 2 || IsField == 3
        s2j                 = sum(reshape(abs(s2j),3,length(Lvj)/3),1)';
    end
    clearvars Ke;    
    stat                    = sqrt(2)*s2j./sqrt(var(s2j));
    indms                   = find(stat > threshold);  
    J                       = s2j;
    J                       = J*scaleSvv/scaleKe^2;
    Jsp                     = zeros(length(stat),1);
    Jsp(indms)              = J(indms);
elseif IsCurv == 1
    param.flag    = strcat(flag," Giri compensation");
    [s2j_giri,sigma2j_giri,Tgiri,~,scaleSvv_giri,scaleKe_giri] = sSSBLpp(Svv,subject.Ke_giri,param);
    param.flag              = strcat(flag," Sulci compensation");
    [s2j_sulc,sigma2j_sulc,Tsulc,~,scaleSvv_sulc,scaleKe_sulc] = sSSBLpp(Svv,subject.Ke_sulc,param);
    clearvars param Svv;
    if IsField == 2 || IsField == 3
        s2j_giri            = sum(reshape(abs(s2j_giri),3,length(Lvj)/3),1)';
        s2j_sulc            = sum(reshape(abs(s2j_sulc),3,length(Lvj)/3),1)';
    end
    clearvars Ke;
    stat_giri               = sqrt(2)*s2j_giri/sqrt(var(s2j_giri));
    indms_giri              = find(stat_giri > threshold);
    stat_sulc               = sqrt(2)*s2j_sulc/sqrt(var(s2j_sulc));
    indms_sulc              = find(stat_sulc > threshold);   
    s2j                     = [s2j_giri s2j_sulc];
    sigma2j                 = [sigma2j_giri sigma2j_sulc];
    clearvars sigma2j_post_giri sigma2j_post_sulc;
    scaleSvv                = [scaleSvv_giri scaleSvv_sulc];
    scaleKe                 = [scaleKe_giri scaleKe_sulc];
    stat                    = [stat_giri stat_sulc];
    clearvars stat_giri stat_sulc;
    indms                   = unique([indms_giri;indms_sulc]);
    clearvars indms_sulc indms_giri;
    J                       = (s2j_giri + s2j_sulc)/2;
    clearvars s2j_giri s2j_sulc;
    J                       = J*sqrt(scaleSvv_giri*scaleSvv_sulc)/(scaleKe_giri*scaleKe_sulc);
    clearvars scaleSvv_giri scaleSvv_sulc;
    clearvars scaleKe_giri scaleKe_sulc;
    Jsp                     = zeros(length(stat),1);
    Jsp(indms)              = J(indms);
    T                       = cat(3,Tgiri,Tsulc);  
    clearvars Tgiri Tsulc;
end
% Ordering results by FSAve indices
J_FSAve                     = zeros(length(J),1);
Jsp_FSAve                   = zeros(length(J),1);
for h=1:length(sub_to_FSAve)
    indices                 = sub_to_FSAve(h,:);    
    J_FSAve(h)              = (J(indices(1))+J(indices(2))+J(indices(3)))/3;
    Jsp_FSAve(h)            = (Jsp(indices(1))+Jsp(indices(2))+Jsp(indices(3)))/3; 
end
outputs.s2j                 = s2j;
outputs.sigma2j             = sigma2j;
outputs.T                   = T;
outputs.scaleSvv            = scaleSvv;
outputs.scaleKe             = scaleKe;
outputs.stat                = stat;
outputs.J                   = J;
outputs.Jsp                 = Jsp;
outputs.indms               = indms;
outputs.J_FSAve             = J_FSAve;
outputs.Jsp_FSAve           = Jsp_FSAve;
end
