function [subject,properties] = activation_level_sssblpp(subject,properties)
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
cmap_a        = properties.cmap_a;
cmap          = properties.cmap;
Sc            = subject.Scortex;
sub_to_FSAve  = subject.sub_to_FSAve;
pathname      = subject.pathname;

%%
%% Sensor level Outputs
%%
sensor_level_out    = subject.sensor_level_out;
Nseg                = sensor_level_out.Nseg;
peak_pos            = sensor_level_out.peak_pos;
Svv                 = sensor_level_out.Svv;
band                = sensor_level_out.band;
str_band            = band.str_band;

%%
%% sSSBL++ activation parameters
%%
activation_params   = properties.activation_params;
sssblpp_th          = activation_params.sssblpp_th.value;
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
        s2j               = sum(reshape(abs(s2j),3,length(Lvj)/3),1)';
    end
    clearvars Ke;
    figure_name = strcat('BC-VARETA-stat - ',str_band);    
    figure_BC_VARETA_stat = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on;    
    stat             = sqrt(2)*s2j./sqrt(var(s2j));
    indms            = find(stat > sssblpp_th);
    h                = histogram(stat,1000,'Normalization','pdf');
    bins             = h.BinEdges;
    pdf              = (1/(sqrt(2*pi)))*(bins.^(-1/2)).*exp(-bins/2);
    hold on
    plot(bins,pdf,'LineWidth',2,'Color','r')
    legend('Empirical','Chi2');
    title('BC-VARETA-stat','Color','k','FontSize',16);
    disp('-->> Saving figure');
    file_name = strcat('BC_VARETA_stat','_',str_band,'.fig');
    saveas(figure_BC_VARETA_stat,fullfile(pathname,file_name));
    pause(1e-12);
    close(figure_BC_VARETA_stat);
    
    J                     = s2j;
    J                     = J*scaleSvv/scaleKe^2;
    Jsp                   = zeros(length(stat),1);
    Jsp(indms)            = J(indms);
elseif IsCurv == 1
    param.flag    = strcat(flag," Giri compensation");
    [s2j_giri,sigma2j_giri,Tgiri,~,scaleSvv_giri,scaleKe_giri] = sSSBLpp(Svv,subject.Ke_giri,param);
    param.flag    = strcat(flag," Sulci compensation");
    [s2j_sulc,sigma2j_sulc,Tsulc,~,scaleSvv_sulc,scaleKe_sulc] = sSSBLpp(Svv,subject.Ke_sulc,param);
    clearvars param Svv;
    if IsField == 2 || IsField == 3
        s2j_giri               = sum(reshape(abs(s2j_giri),3,length(Lvj)/3),1)';
        s2j_sulc               = sum(reshape(abs(s2j_sulc),3,length(Lvj)/3),1)';
    end
    clearvars Ke;
    figure_name = strcat('BC-VARETA-stat - ',str_band);
    
    figure_BC_VARETA_stat = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on;
    
    subplot(1,2,1);
    stat_giri             = sqrt(2)*s2j_giri/sqrt(var(s2j_giri));
    indms_giri            = find(stat_giri > sssblpp_th);
    h                     = histogram(stat_giri,1000,'Normalization','pdf');
    bins                  = h.BinEdges;
    pdf                   = (1/(sqrt(2*pi)))*(bins.^(-1/2)).*exp(-bins/2);
    hold on
    plot(bins,pdf,'LineWidth',2,'Color','r')
    legend('Empirical','Chi2');
    xlabel('amplitude');
    ylabel('probabilities');
    title('BC-VARETA-stat-giri','Color','k','FontSize',16);
    subplot(1,2,2);
    stat_sulc             = sqrt(2)*s2j_sulc/sqrt(var(s2j_sulc));
    indms_sulc            = find(stat_sulc > sssblpp_th);
    h = histogram(stat_sulc,1000,'Normalization','pdf');
    bins = h.BinEdges;
    pdf = (1/(sqrt(2*pi)))*(bins.^(-1/2)).*exp(-bins/2);
    hold on
    plot(bins,pdf,'LineWidth',2,'Color','r')
    legend('Empirical','Chi2');
    xlabel('amplitude');
    ylabel('probabilities');
    title('BC-VARETA-stat-sulc','Color','k','FontSize',16);
    disp('-->> Saving figure');
    file_name = strcat('BC_VARETA_stat','_',str_band,'.fig');
    saveas(figure_BC_VARETA_stat,fullfile(pathname,file_name));
    pause(1e-12);   
    close(figure_BC_VARETA_stat);
    
    s2j                   = [s2j_giri s2j_sulc];
    sigma2j               = [sigma2j_giri sigma2j_sulc];
    clearvars sigma2j_post_giri sigma2j_post_sulc;
    scaleSvv              = [scaleSvv_giri scaleSvv_sulc];
    scaleKe               = [scaleKe_giri scaleKe_sulc];
    stat                  = [stat_giri stat_sulc];
    clearvars stat_giri stat_sulc;
    indms                 = unique([indms_giri;indms_sulc]);
    clearvars indms_sulc indms_giri;
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

%%
%% Plotting results
%%
sources_iv          = sqrt(abs(J));
sources_iv          = sources_iv/max(sources_iv(:));

figure_name = strcat('BC-VARETA-activation - ',str_band);
figure_BC_VARETA1 = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on;
smoothValue          = 0.66;
SurfSmoothIterations = 10;
Vertices             = tess_smooth(Sc.Vertices, smoothValue, SurfSmoothIterations, Sc.VertConn, 1);

define_ico(figure_BC_VARETA1);
patch('Faces',Sc.Faces,'Vertices',Vertices,'FaceVertexCData',Sc.SulciMap*0.06+...
    log(1+sources_iv),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
set(gca,'xcolor','w','ycolor','w','zcolor','w');
az = 0; el = 0;
view(az,el);
rotate3d on;
colormap(gca,cmap);
title('BC-VARETA-activation','Color','k','FontSize',16);
axis square;
disp('-->> Saving figure');
file_name = strcat('BC_VARETA_activation','_',str_band,'.fig');
saveas(figure_BC_VARETA1,fullfile(pathname,file_name));

close(figure_BC_VARETA1);

%% Saving files
disp('-->> Saving file')
subject.file_name = strcat('MEEG_source_',str_band,'.mat');
disp(strcat("File: ", subject.file_name));
parsave(fullfile(subject.pathname ,subject.file_name ),s2j,sigma2j,T,scaleSvv,scaleKe,stat,J,Jsp,indms,J_FSAve,Jsp_FSAve);


end
