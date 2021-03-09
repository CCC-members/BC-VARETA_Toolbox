function [subject,properties] = sensor_level_analysis(Svv_channel,PSD,Nf,F,subject,properties)

disp('=================================================================');
band = properties.band;
if(properties.general_params.run_frequency_bin.value)
    disp(strcat( 'BC-V-->> Sensor level for frequency band: (' , band.name , ') bin ->>>' , string(band.f_bin), 'Hz') );
    properties.str_band =  strcat( band.name,'_',string(band.f_bin),'Hz');
else
    disp(strcat( 'BC-V-->> Sensor level for frequency band: (' , band.name , ') ' , string(band.f_start), 'Hz-->' , string(band.f_end) , 'Hz') );
    properties.str_band =  strcat( band.name,'_',string(band.f_start),'Hz_',string(band.f_end),'Hz');
end
text_level      = 'Sensor_level';
if(properties.general_params.run_by_trial.value)
    trial_name  = properties.trial_name;
    pathname    = fullfile(subject.subject_path,trial_name,text_level,band.name);
else
    pathname    = fullfile(subject.subject_path,text_level,band.name);
end
if(~isfolder(pathname))
    mkdir(pathname);
end
properties.pathname = pathname;

%%
%% Get band
%%
[Svv,subject,properties] = get_band(Svv_channel,PSD,Nf,F,subject,properties);

%%
%% Preparing params
%%

Ke           = subject.Ke;
Cdata        = subject.Cdata;
Sh           = subject.Sh;
cmap_a       = properties.cmap_a;
cmap_c       = properties.cmap_c;
str_band     = properties.str_band;
Nseg         = properties.Nseg;

%%
%% Test
%%

figure_name = strcat('Scalp 2D - ',str_band);
if(properties.run_bash_mode.disabled_graphics)
    figure_scalp_2D = figure('Color','w','Name',figure_name,'NumberTitle','off','visible','off'); hold on; set(gca,'Color','w');
else
    figure_scalp_2D = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on; set(gca,'Color','w');
end
define_ico(figure_scalp_2D);

Channel          = Cdata.Channel;
%%
elec_data               = [];
elec_data.pos           = zeros(length(Channel),3);
for ii = 1:length(Channel)
    elec_data.lbl{ii}   = Channel(ii).Name;
    temp                = Channel(ii).Loc;
    elec_data.pos(ii,:) = mean(temp,2);
end
elec_data.label         = elec_data.lbl;
elec_data.elecpos       = elec_data.pos;
elec_data.unit          = 'mm';


temp    = diag(Svv);
temp    = abs(temp)/max(abs(temp(:)));
cfg     = [];
topo    = [];
if(isequal(subject.modality,'MEG'))
    %% MEG topography
    cfg.layout          = '4D248_helmet.mat';
    cfg.channel         = 'meg';
    cfg.markers         = '.';
    cfg.markersymbol    = '.';
    cfg.colormap        = cmap_a;
    cfg.markersize      = 3;
    cfg.markercolor     = [1 1 1];
    topo.sens           = elec_data;
    topo.tra            = elec_data.pos;
    topo.coilpos        = elec_data.pos;
    topo.label          = elec_data.lbl';
    topo.dimord         = 'chan_freq';
    topo.freq           = 1;
    topo.powspctrm      = temp;    
else
    %% EEG topography
    cfg.marker          = '';
    cfg.layout          = properties.sensor_params.fieldtrip.layout.value;   
    cfg.channel         = 'eeg';
    cfg.markersymbol    = '.';
    cfg.colormap        = cmap_a;
    cfg.markersize      = 3;
    cfg.markercolor     = [1 1 1];
    topo.elec           = elec_data;
    topo.label          = elec_data.lbl;
    topo.dimord         = 'chan_freq';
    topo.freq           = 1;
    topo.powspctrm      = temp;
end

ft_topoplotTFR(cfg,topo);
title(['MEG' ' ' band.name ' ' 'topography'])

disp('-->> Saving figure');
file_name = strcat('Scalp_2D','_',str_band,'.fig');
saveas(figure_scalp_2D,fullfile(properties.pathname,file_name));

close(figure_scalp_2D);

%%
%%
%%
%%


%%
%% topography...
%%
Nelec = size(Ke,1);
Svv_inv = sqrtm(Svv*Svv+4*eye(Nelec))-Svv;
Loc = [Cdata.Channel.Loc];
if(isequal(subject.modality,'MEG'))
    Loc = squeeze(mean(reshape(Loc,3,4,Nelec),2));
end
for ii = 1:length(Loc)
    X(ii) = Loc(1,ii);
    Y(ii) = Loc(2,ii);
    Z(ii) = Loc(3,ii);
end
C = abs(diag(Svv));
C = C/max(C);
C(C<0.01) = 0;

figure_name = strcat('Scalp 3D - ',str_band);
if(properties.run_bash_mode.disabled_graphics)
    figure_scalp_3D = figure('Color','w','Name',figure_name,'NumberTitle','off','visible','off'); hold on; set(gca,'Color','w');
else
    figure_scalp_3D = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on; set(gca,'Color','w');
end
define_ico(figure_scalp_3D);
scatter3(X,Y,Z,100,C.^1,'filled');
patch('Faces',Sh.Faces,'Vertices',Sh.Vertices,'FaceVertexCData',0.01*(ones(length(Sh.Vertices),1)),'FaceColor','interp','EdgeColor','none','FaceAlpha',.35);
colormap(gca,cmap_a);
az = 0; el = 0;
view(az, el);
title('Scalp','Color','k','FontSize',16);
axis equal;
axis off;
disp('-->> Saving figure');
file_name = strcat('Scalp_3D','_',str_band,'.fig');
saveas(figure_scalp_3D,fullfile(properties.pathname,file_name));

close(figure_scalp_3D);

%%
%% inverse covariance matrix...
%%
temp_diag  = diag(diag(abs(Svv_inv)));
temp_ndiag = abs(Svv_inv)-temp_diag;
temp_ndiag = temp_ndiag/max(temp_ndiag(:));
temp_diag  = diag(abs(diag(Svv)));
temp_diag  = temp_diag/max(temp_diag(:));
temp_diag  = diag(diag(temp_diag)+1);
temp_comp  = temp_diag+temp_ndiag;

figure_name = strcat('Scalp - ',str_band);
if(properties.run_bash_mode.disabled_graphics)
    figure_scalp_electrodes = figure('Color','w','Name',figure_name,'NumberTitle','off','visible','off');
else
    figure_scalp_electrodes = figure('Color','w','Name',figure_name,'NumberTitle','off');
end
define_ico(figure_scalp_electrodes);
imagesc(temp_comp);
set(gca,'Color','w','XColor','k','YColor','k','ZColor','k',...
    'XTick',1:length(Loc),'YTick',1:length(Loc),...
    'XTickLabel',{Cdata.Channel.Name},'XTickLabelRotation',90,...
    'YTickLabel',{Cdata.Channel.Name},'YTickLabelRotation',0);

xlabel('electrodes','Color','k');
ylabel('electrodes','Color','k');
colormap(gca,cmap_c);
colorbar;
axis square;
title('Scalp','Color','k','FontSize',16);

disp('-->> Saving figure');
file_name = strcat('Covariance_Matrix','_',str_band,'.fig');
saveas(figure_scalp_electrodes,fullfile(properties.pathname,file_name));

close(figure_scalp_electrodes);

% subject.Svv = Svv;

%% Saving files
disp('-->> Saving file')
file_name = strcat('Sensor_level_',str_band,'.mat');
disp(strcat("File: ", file_name));
peak_pos = properties.peak_pos;
parsave(fullfile(properties.pathname ,file_name ),Svv,peak_pos,Nseg,band);

reference_path = strsplit(properties.pathname,subject.name);
if(properties.general_params.run_by_trial.value)
    if(properties.general_params.run_frequency_bin.value)
        f_bin = replace(num2str(band.f_bin),'.','_');
        f_bin = strcat(band.name,'_',f_bin);
        properties.BC_V_info.(trial_name).sensor_level.(band.name).(f_bin).name         = file_name;
        properties.BC_V_info.(trial_name).sensor_level.(band.name).(f_bin).ref_path     = reference_path{2};
    else
        properties.BC_V_info.(properties.trial_name).sensor_level.(band.name).name      = file_name;
        properties.BC_V_info.(properties.trial_name).sensor_level.(band.name).ref_path  = reference_path{2};
    end
else
    if(properties.general_params.run_frequency_bin.value)
        f_bin = replace(num2str(band.f_bin),'.','_');
        f_bin = strcat(band.name,'_',f_bin);
        properties.BC_V_info.sensor_level.(band.name).(f_bin).name      = file_name;
        properties.BC_V_info.sensor_level.(band.name).(f_bin).ref_path  = reference_path{2};
    else
        properties.BC_V_info.sensor_level.(band.name).name              = file_name;
        properties.BC_V_info.sensor_level.(band.name).ref_path          = reference_path{2};
    end
end
pause(1e-12);

end

