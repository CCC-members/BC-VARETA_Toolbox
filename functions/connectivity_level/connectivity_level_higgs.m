function [Thetajj,s2j,Tjv,llh] = connectivity_level_higgs(subject,properties)

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
Ke                  = subject.Ke;
W                   = subject.W;
Winv                = subject.Winv;
Cdata               = subject.Cdata;
cmap                = load(properties.general_params.colormap_path);
cmap_a              = cmap.cmap_a;
cmap_c              = cmap.cmap_c;
cmap                = cmap.cmap;
Sc                  = subject.Sc;
Atlas               = Sc.Atlas(Sc.iAtlas).Scouts;
run_bash_mode       = properties.run_bash_mode.value;

%%
%% Sensor level Outputs
%%
sensor_level_out    = properties.sensor_level_out;
Svv                 = sensor_level_out.Svv;
Nseg                = sensor_level_out.Nseg;
peak_pos            = sensor_level_out.peak_pos;
band                = sensor_level_out.band;

%%
%% Activation level Outputs
%%
activation_level_out    = properties.activation_level_out;
% s2j                     = activation_level_out.s2j;
% sigma2j_post            = activation_level_out.sigma2j_post;
   
%%
%% HIGGS parameters
%%
connectivity_params = properties.connectivity_params;
higgs_th            = connectivity_params.higgs_th.value;
IsCurv              = connectivity_params.IsCurv.value; % 0 (no compensation) 1 (giri and sulci curvature compensation)
IsField             = connectivity_params.IsField.value; % 1 (projected Lead Field) 3 (3D Lead Field)

%% Defining path
disp('=================================================================');
if(isfield(band,'f_bin'))
    disp(strcat( 'BC-V-->> Connectivity level for frequency band: (' , band.name , ') bin ->>>' , string(band.f_bin), 'Hz') );
    str_band =  strcat( band.name,'_',string(band.f_bin),'Hz');
else
    disp(strcat( 'BC-V-->> Connectivity level for frequency band: (' , band.name , ') ' , string(band.f_start), 'Hz-->' , string(band.f_end) , 'Hz') );
    str_band =  strcat( band.name,'_',string(band.f_start),'Hz_',string(band.f_end),'Hz');
end
text_level = 'Connectivity_level';
if(properties.general_params.run_by_trial.value) 
    trial_name = properties.trial_name;
    pathname = fullfile(subject.subject_path,trial_name,text_level,'HIGGS',band.name);
else
    pathname = fullfile(subject.subject_path,text_level,'HIGGS',band.name);
end
if(~isfolder(pathname))
    mkdir(pathname);
end

%%
%% HIGGS parameters
%%
param.use_gpu         = properties.run_bash_mode.use_gpu;
m                     = Nseg;
param.run_bash_mode   = run_bash_mode;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if IsCurv == 0    
%     if IsField == 2 || IsField == 3
%         s2j               = sum(reshape(abs(s2j),3,length(Ke)/3),1)';
%         sigma2j_post      = sum(reshape(abs(sigma2j_post),3,length(Ke)/3),1)';
%     end
%     stat                  = sqrt(s2j./sigma2j_post);
%     indms                 = find(stat > higgs_th);    
% elseif IsCurv == 1
%     Ke_giri                          = subject.Ke_giri;
%     Ke_sulc                          = subject.Ke_sulc;
%     Ke_giri                          = Ke_giri*W;
%     Ke_sulc                          = Ke_sulc*W;
%     if IsField == 2 || IsField == 3
%         s2j_giri               = sum(reshape(abs(s2j_giri),3,length(Ke)/3),1)';
%         sigma2j_post_giri      = sum(reshape(abs(sigma2j_post_giri),3,length(Ke)/3),1)';
%         s2j_sulc               = sum(reshape(abs(s2j_sulc),3,length(Ke)/3),1)';
%         sigma2j_post_sulc      = sum(reshape(abs(sigma2j_post_sulc),3,length(Ke)/3),1)';
%     end
%     stat_giri             = sqrt(s2j_giri./sigma2j_post_giri);
%     indms_giri            = find(stat_giri > sssblpp_th);
%     stat_sulc             = sqrt(s2j_sulc./sigma2j_post_sulc);
%     indms_sulc            = find(stat_sulc > sssblpp_th);
%     s2j                   = [s2j_giri s2j_sulc];
%     sigma2j_post          = [sigma2j_post_giri sigma2j_post_sulc];
%     clearvars sigma2j_post_giri sigma2j_post_sulc;
%     scaleSvv              = [scaleSvv_giri scaleSvv_sulc];
%     scaleKe               = [scaleKe_giri scaleKe_sulc];
%     stat                  = [stat_giri stat_sulc];
%     clearvars stat_giri stat_sulc;
%     indms                 = unique([indms_giri;indms_sulc]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indms                 = activation_level_out.indms;
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
if IsCurv == 0
    Ke                               = Ke*W;   
    [Thetajj,Tjv,llh]                = higgs(Svv,Ke(:,indms),param);
    [Thetajj,s2j,Tjv]                = higgs_destandardization(Thetajj,Svv,Tjv,Winv,W,indms,IsField);
elseif IsCurv == 1
    Ke_giri                          = subject.Ke_giri;
    Ke_sulc                          = subject.Ke_sulc;
    Ke_giri                          = Ke_giri*W;
    Ke_sulc                          = Ke_sulc*W;
    [Thetajj_sulc,Tjv_sulc,llh_sulc] = higgs(Svv,Ke_sulc(:,indms),param);
    [Thetajj_giri,Tjv_giri,llh_giri] = higgs(Svv,Ke_giri(:,indms),param);
    [Thetajj_giri,s2j_giri,Tjv_giri] = higgs_destandardization(Thetajj_giri,Svv,Tjv_giri,Winv,W,indms,IsField);
    [Thetajj_sulc,s2j_sulc,Tjv_sulc] = higgs_destandardization(Thetajj_sulc,Svv,Tjv_sulc,Winv,W,indms,IsField);
    Thetajj                          = (Thetajj_giri + Thetajj_sulc)/2;
    s2j                              = (s2j_giri + s2j_sulc)/2;
    llh                              = [llh_giri llh_sulc];
    Tjv                              = cat(3,Tjv_giri,Tjv_sulc);
end

%%
%% Plotting results
%%
J                   = s2j;
sources_iv          = zeros(length(J),1);
sources_iv(indms)   = sqrt(abs(J(indms)));
sources_iv          = sources_iv/max(sources_iv(:));

figure_name = strcat('BC-VARETA-activation - ',str_band);
if(properties.run_bash_mode.disabled_graphics)
    figure_BC_VARETA1 = figure('Color','w','Name',figure_name,'NumberTitle','off','visible','off'); hold on;
else
    figure_BC_VARETA1 = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on;
end
define_ico(figure_BC_VARETA1);
patch('Faces',Sc.Faces,'Vertices',Sc.Vertices,'FaceVertexCData',sources_iv,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
set(gca,'Color','w');
az = 0; el = 0;
view(az,el);
colormap(gca,cmap_a);
title('BC-VARETA-activation','Color','k','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_activation','_',str_band,'.fig');
saveas(figure_BC_VARETA1,fullfile(pathname,file_name));

close(figure_BC_VARETA1);

pause(1e-12);

%%
%% Plotting results

temp_iv    = abs(s2j(indms));
connect_iv = abs(Thetajj);
temp       = abs(connect_iv);
temp_diag  = diag(diag(temp));
temp_ndiag = temp - temp_diag;
temp_ndiag = temp_ndiag/max(temp_ndiag(:));
temp_diag  = diag(temp_iv);
temp_diag  = temp_diag/max(temp_diag(:));
temp_diag  = diag(diag(temp_diag) + 1);
temp_comp  = temp_diag + temp_ndiag;
label_gen = [];
for ii = 1:length(indms)
    label_gen{ii} = num2str(ii);
end

figure_name = strcat('BC-VARETA-node-wise-conn - ',str_band);
if(properties.run_bash_mode.disabled_graphics)
    figure_BC_VARETA2 = figure('Color','w','Name',figure_name,'NumberTitle','off','visible','off');
else
    figure_BC_VARETA2 = figure('Color','w','Name',figure_name,'NumberTitle','off');
end
define_ico(figure_BC_VARETA2);
imagesc(temp_comp);
set(gca,'Color','w','XColor','k','YColor','k','ZColor','k',...
    'XTick',1:length(indms),'YTick',1:length(indms),...
    'XTickLabel',label_gen,'XTickLabelRotation',90,...
    'YTickLabel',label_gen,'YTickLabelRotation',0);
xlabel('sources','Color','k');
ylabel('sources','Color','k');
colorbar;
colormap(gca,cmap_c);
axis square;
title('BC-VARETA-node-wise-conn','Color','k','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_node_wise_conn','_',str_band,'.fig');
saveas(figure_BC_VARETA2,fullfile(pathname,file_name));

close(figure_BC_VARETA2);

%% Roi analysis
Thetajj_full              = zeros(length(Ke));
Sjj_full                  = zeros(length(Ke));
Thetajj_full(indms,indms) = Thetajj;
Sjj_full(indms,indms)     = diag(temp_iv);
atlas_label               = cell(1,length(Atlas));
conn_roi                  = zeros(length(Atlas));
act_roi                   = zeros(length(Atlas),1);
for roi1 = 1:length(Atlas)
    for roi2 = 1:length(Atlas)
        conn_tmp             = Thetajj_full(Atlas(roi1).Vertices,Atlas(roi2).Vertices);
        conn_tmp             = mean(abs(conn_tmp(:)));
        conn_roi(roi1,roi2)  = conn_tmp;
    end
    atlas_label{roi1} = Atlas(roi1).Label;
end

for roi1 = 1:length(Atlas)
    act_tmp              = diag(Sjj_full(Atlas(roi1).Vertices,Atlas(roi1).Vertices));
    act_tmp              = mean(abs(act_tmp));
    act_roi(roi1)        = act_tmp;
end
act_roi    = diag(act_roi);
temp_iv    = abs(act_roi);
connect_iv = abs(conn_roi);
temp       = abs(connect_iv);
temp_diag  = diag(diag(temp));
temp_ndiag = temp-temp_diag;
temp_ndiag = temp_ndiag/max(temp_ndiag(:));
temp_diag  = diag(abs(diag(temp_iv)));
temp_diag  = temp_diag/max(temp_diag(:));
temp_diag  = diag(diag(temp_diag)+1);
temp_comp  = temp_diag+temp_ndiag;

figure_name = strcat('BC-VARETA-roi-conn - ',str_band);
if(properties.run_bash_mode.disabled_graphics)
    figure_BC_VARETA3 = figure('Color','w','Name',figure_name,'NumberTitle','off','visible','off');
else
    figure_BC_VARETA3 = figure('Color','w','Name',figure_name,'NumberTitle','off');
end
define_ico(figure_BC_VARETA3);
imagesc(temp_comp);
set(gca,'Color','w','XColor','k','YColor','k','ZColor','k',...
    'XTick',1:length(Atlas),'YTick',1:length(Atlas),...
    'XTickLabel',atlas_label,'XTickLabelRotation',90,...
    'YTickLabel',atlas_label,'YTickLabelRotation',0);
xlabel('sources','Color','k');
ylabel('sources','Color','k');
colorbar;
colormap(gca,cmap_c);
axis square;
title('BC-VARETA-roi-conn','Color','k','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_roi_conn','_',str_band,'.fig');
saveas(figure_BC_VARETA3,fullfile(pathname,file_name));

close(figure_BC_VARETA3);


%% Saving files
disp('-->> Saving file.')
disp(strcat("Path: ",pathname));
file_name = strcat('MEEG_source_',str_band,'.mat');
disp(strcat("File: ", file_name));
parsave(fullfile(pathname ,file_name ),Thetajj,s2j,Tjv,llh,Svv,W,indms);


reference_path = strsplit(pathname,subject.name);
if(properties.general_params.run_by_trial.value) 
    if(properties.general_params.run_frequency_bin.value)
        properties.BC_V_info.connectivity_level.(trial_name).(band.name).(band.f_bin).name = file_name;
        properties.BC_V_info.connectivity_level.(trial_name).(band.name).(band.f_bin).ref_path = reference_path{2};
    else
        properties.BC_V_info.connectivity_level.(trial_name).(band.name).name = file_name;
        properties.BC_V_info.connectivity_level.(trial_name).(band.name).ref_path = reference_path{2};
    end
else
    if(properties.general_params.run_frequency_bin.value)
        properties.BC_V_info.connectivity_level.(band.name).(band.f_bin).name = file_name;
        properties.BC_V_info.connectivity_level.(band.name).(band.f_bin).ref_path = reference_path{2};
    else
        properties.BC_V_info.connectivity_level.(band.name).name = file_name;
        properties.BC_V_info.connectivity_level.(band.name).ref_path = reference_path{2};
    end
end

pause(1e-12);

end

