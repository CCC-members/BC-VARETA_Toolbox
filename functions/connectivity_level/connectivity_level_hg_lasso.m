function [Thetajj,Sjj,Sigmajj] = connectivity_level_hg_lasso(subject,properties)

% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019

% Updates
% - Ariosky Areces Gonzalez

% Date: March 20, 2019

Ke                  = subject.Ke;
Cdata               = subject.Cdata;
Sh                  = subject.Sh;
cmap                = load(properties.general_params.colormap_path);
cmap_a              = cmap.cmap_a;
cmap_c              = cmap.cmap_c;
Sc                  = subject.Sc;
GridOrient          = subject.GridOrient;
GridAtlas           = subject.GridAtlas;
Atlas               = Sc.Atlas(Sc.iAtlas).Scouts;
Nseg                = properties.sensor_level_out.Nseg;
peak_pos            = properties.sensor_level_out.peak_pos;
band                = properties.sensor_level_out.band;
run_bash_mode       = properties.run_bash_mode.value;
Svv                 = properties.sensor_level_out.Svv;
indms               = properties.activation_level_out.indms;
Tjv                 = properties.activation_level_out.T;
scaleSvv            = properties.activation_level_out.scaleSvv;
scaleKe             = properties.activation_level_out.scaleKe;
activation_params   = properties.activation_params;
connectivity_params = properties.connectivity_params;
IsCurv              = activation_params.IsCurv.value;  % 0 (no compensation) 1 (giri and sulci curvature compensation)
IsField             = activation_params.IsField.value; % 1 (projected Lead Field) 3 (3D Lead Field)

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
    pathname = fullfile(subject.subject_path,trial_name,text_level,'sSSBLpp',band.name);
else
    pathname = fullfile(subject.subject_path,text_level,'sSSBLpp',band.name);
end
if(~isfolder(pathname))
    mkdir(pathname);
end

%%
%% bc-vareta toolbox...
%% HG-LASSO parameters
%%

param.use_gpu         = properties.run_bash_mode.use_gpu;
m                     = length(peak_pos)*Nseg;
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
%%
%% Plotting results
temp_iv    = abs(Sjj);
connect_iv = abs(Thetajj);
temp       = abs(connect_iv);
temp_diag  = diag(diag(temp));
temp_ndiag = temp-temp_diag;
temp_ndiag = temp_ndiag/max(temp_ndiag(:));
temp_diag  = diag(abs(diag(temp_iv)));
temp_diag  = temp_diag/max(temp_diag(:));
temp_diag  = diag(diag(temp_diag)+1);
temp_comp  = temp_diag+temp_ndiag;
label_gen = [];
for ii = 1:length(indms)
    label_gen{ii} = num2str(ii);
end

figure_name = strcat('BC-VARETA-node-wise-conn - ',str_band);
if(properties.run_bash_mode.disabled_graphics)
    figure_BC_VARETA2 = figure('Color','k','Name',figure_name,'NumberTitle','off','visible','off');
else
    figure_BC_VARETA2 = figure('Color','k','Name',figure_name,'NumberTitle','off');
end
define_ico(figure_BC_VARETA2);
imagesc(temp_comp);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length(indms),'YTick',1:length(indms),...
    'XTickLabel',label_gen,'XTickLabelRotation',90,...
    'YTickLabel',label_gen,'YTickLabelRotation',0);
xlabel('sources','Color','w');
ylabel('sources','Color','w');
colorbar;
colormap(gca,cmap_c);
axis square;
title('BC-VARETA-node-wise-conn','Color','w','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_node_wise_conn','_',str_band,'.fig');
saveas(figure_BC_VARETA2,fullfile(pathname,file_name));

close(figure_BC_VARETA2);

%% Roi analysis
Thetajj_full              = zeros(length(Ke)/3);
Sjj_full                  = zeros(length(Ke)/3);
Thetajj_full(indms,indms) = Thetajj;
Sjj_full(indms,indms)     = Sjj;
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
    figure_BC_VARETA3 = figure('Color','k','Name',figure_name,'NumberTitle','off','visible','off');
else
    figure_BC_VARETA3 = figure('Color','k','Name',figure_name,'NumberTitle','off');
end
define_ico(figure_BC_VARETA3);
imagesc(temp_comp);
set(gca,'Color','k','XColor','w','YColor','w','ZColor','w',...
    'XTick',1:length(Atlas),'YTick',1:length(Atlas),...
    'XTickLabel',atlas_label,'XTickLabelRotation',90,...
    'YTickLabel',atlas_label,'YTickLabelRotation',0);
xlabel('sources','Color','w');
ylabel('sources','Color','w');
colorbar;
colormap(gca,cmap_c);
axis square;
title('BC-VARETA-roi-conn','Color','w','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('BC_VARETA_roi_conn','_',str_band,'.fig');
saveas(figure_BC_VARETA3,fullfile(pathname,file_name));

close(figure_BC_VARETA3);


%% Saving files
disp('-->> Saving file.')
disp(strcat("Path: ",pathname));
file_name = strcat('MEEG_source_',str_band,'.mat');
disp(strcat("File: ", file_name));
parsave(fullfile(pathname ,file_name ),Thetajj,Sjj,Sigmajj);


reference_path = strsplit(properties.pathname,subject.name);
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

