function PlotSensorTopography(reportPath,sensorData,Cdata)

if(~isfolder(reportPath))
    mkdir(reportPath);
end

template = load('axes.mat');
currentAxes = template.axes;
colorMap = load('tools/mycolormap_brain_basic_conn.mat');
fig = figure("Name","Sensor Topography","Color","w", 'Position', get(0, 'Screensize'),'Visible','off');
set(currentAxes,"Parent",fig);
Svv = sensorData.Svv;
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

% if(isequal(subject.modality,'MEG'))
%     %% MEG topography
%     if(isequal(properties.sensor_params.fieldtrip.layout.value,'4D248_helmet.mat'))
%         cfg.layout          = properties.sensor_params.fieldtrip.layout.value;
%         cfg.channel         = 'meg';
%         cfg.markers         = '.';
%         cfg.markersymbol    = '.';
%         cfg.colormap        = cmap_a;
%         cfg.markersize      = 3;
%         cfg.markercolor     = [1 1 1];
%     elseif(isequal(properties.sensor_params.fieldtrip.layout.value,'another'))
%
%     end
%     topo.sens           = elec_data;
%     topo.tra            = elec_data.pos;
%     topo.coilpos        = elec_data.pos;
%     topo.label          = elec_data.lbl';
%     topo.dimord         = 'chan_freq';
%     topo.freq           = 1;
%     topo.powspctrm      = SvvDiag;
% else
%% EEG topography
cfg.marker          = 'on';
cfg.layout          = 'EEG1005.lay';
cfg.channel         = elec_data.label';
cfg.markersymbol    = '.';
cfg.colormap        = colorMap.cmap_a;
cfg.markersize      = 20;
cfg.markercolor     = 'g';
cfg.comment         = 'no';
cfg.interactive     = 'no';
cfg.colorbar        = 'no';
% Topo
topo.elec           = elec_data;
topo.label          = elec_data.lbl;
topo.dimord         = 'chan_freq';
topo.freq           = 1;
topo.powspctrm      = SvvDiag;
% end
ft_topoplotTFR(cfg,topo);
view(currentAxes,0,0);
view(currentAxes,0,90);

export_fig(fullfile(reportPath,strcat('sensorTopography_',sensorData.band.str_band)),'-transparent','-png');
close(fig);
end

