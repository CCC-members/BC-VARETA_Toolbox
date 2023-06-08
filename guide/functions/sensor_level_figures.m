

%% Adding fieltrip external functions
f_path          = mfilename('fullpath');
[ref_path,~,~]  = fileparts(fileparts(f_path));
addpath(genpath(fullfile(ref_path,'external/fieldtrip')));
ft_defaults


%%
%% Plot Power Spectral Density
%%

PSD_log = 10*log10(abs(PSD));
min_psd = min(PSD_log(:));
max_psd = max(PSD_log(:));
plot_peak = min_psd*ones(Nf,1);
peak_pos = nf1:nf2;
properties.peak_pos = peak_pos;
Svv = mean(Svv(:,:,peak_pos),3);
fig_title = strcat("Power Spectral Density - ", band.name,'_',string(band.f_start),'Hz-',string(band.f_end),'Hz');
if(~isfile(fullfile(properties.pathname,strcat(fig_title,'.fig'))))
    % just for the plot
    [f1,nf1] = min(abs(F - band.f_start));
    [f2,nf2] = min(abs(F - band.f_end));
    plot_peak(nf1:nf2) = max_psd;
    
    
    figure_band = figure('Color','w','Name',fig_title,'NumberTitle','off');
    
    %     define_ico(figure_band);
    hold on;
    plot(F,PSD_log);
    plot(F,plot_peak,'--k');
    set(gca,'Color','w','XColor','k','YColor','k');
    ylabel('PSD (dB)','Color','k');
    xlabel('Freq. (Hz)','Color','k');
    title(strcat("Power Spectral Density - ",band.name, " band"),'Color','k');
    
    text_cross = strcat(string(band.f_start),"Hz - ", string(band.f_end), "Hz");
    text(band.f_end,max_psd*0.9,text_cross,'Color','k','FontSize',12,'HorizontalAlignment','center');
    
    pause(1e-10);
    
    % Saving figure band
    saveas(figure_band,fullfile(properties.pathname,strcat(fig_title,'.fig')));
    close(figure_band);
end

%% Sensor level figures

figure_name = strcat('Scalp 2D - ',band.str_band);
figure_scalp_2D = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on; set(gca,'Color','w');
define_ico(figure_scalp_2D);

Channel                 = Cdata.Channel;
elec_data               = [];
elec_data.pos           = zeros(length(Channel),3);
for ii = 1:length(Channel)
    elec_data.lbl{ii}   = Channel(ii).Name;
    loc                 = Channel(ii).Loc;
    elec_data.pos(ii,:) = mean(loc,2);
end
elec_data.label         = elec_data.lbl;
elec_data.elecpos       = elec_data.pos;
elec_data.unit          = 'mm';


SvvDiag                 = diag(Svv);
SvvDiag                 = abs(SvvDiag)/max(abs(SvvDiag(:)));
cfg                     = [];
topo                    = [];
if(isequal(subject.modality,'MEG'))
    %% MEG topography
    if(isequal(properties.sensor_params.fieldtrip.layout.value,'4D248_helmet.mat'))
        cfg.layout          = properties.sensor_params.fieldtrip.layout.value;
        cfg.channel         = 'meg';
        cfg.markers         = '.';
        cfg.markersymbol    = '.';
        cfg.colormap        = cmap_a;
        cfg.markersize      = 3;
        cfg.markercolor     = [1 1 1];
    elseif(isequal(properties.sensor_params.fieldtrip.layout.value,'another'))

    end
    topo.sens           = elec_data;
    topo.tra            = elec_data.pos;
    topo.coilpos        = elec_data.pos;
    topo.label          = elec_data.lbl';
    topo.dimord         = 'chan_freq';
    topo.freq           = 1;
    topo.powspctrm      = SvvDiag;
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
    topo.powspctrm      = SvvDiag;
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
Nelec = size(Lvj,1);
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
figure_scalp_3D = figure('Color','w','Name',figure_name,'NumberTitle','off'); hold on; set(gca,'Color','w');
define_ico(figure_scalp_3D);
scatter3(X,Y,Z,100,C.^1,'filled');
patch('Faces',Sh.Faces,'Vertices',Sh.Vertices,'FaceVertexCData',0.01*(ones(length(Sh.Vertices),1)),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
colormap(gca,cmap_a);
az = 0; el = 0;
view(az, el);
rotate3d on;
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
figure_scalp_electrodes = figure('Color','w','Name',figure_name,'NumberTitle','off');
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