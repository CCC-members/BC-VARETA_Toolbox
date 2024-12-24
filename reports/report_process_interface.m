function result = report_process_interface(dataset,varargin)

for k = 2 : nargin
    eval([inputname(k) '=  varargin{k-1};']);
end


result = "success";
reportPath = fullfile(pwd,'tmp','report','subject');
if(~isfolder(reportPath))
    mkdir(reportPath);
end

switch clasif
    case 'subject'
        subjec_info = nodeData;
        disp(strcat("-->> Exportin report for subject: ", subjec_info.SubID));
        for i=1:length(subjec_info.common)
            common_file = subjec_info.common(i);
            comment = strrep(common_file.Name,'.mat','');
            common.(comment) = load(fullfile(dataset.Path, subjec_info.SubID, common_file.Ref_path,common_file.Name));  
        end
        prop.general_params  = dataset.general_params;

        % Structural processing
        PlotHeadModel(fullfile(reportPath,'structural'),common.Sscalp,common.Souterskull,common.Sinnerskull,common.Cortex);
        PlotSChannel(fullfile(reportPath,'structural'),common.Sscalp,common.Channels);
        PlotLeadfield(fullfile(reportPath,'structural'),common.Leadfield,common.Channels,common.Sscalp,common.Cortex);
        PlotSCortex(fullfile(reportPath,'structural'),common.Cortex);        

        % Sensor PSD
        PlotSensorSpectrum(fullfile(reportPath,'sensor'),dataset, nodeData);

        % Sensor level
        if(isequal(dataset.general_params.analysis_level.value,'all') ...
                || isequal(dataset.general_params.analysis_level.value,'1') ...
                || isequal(dataset.general_params.analysis_level.value,'12'))
            for i=1:length(subjec_info.sensor_level)
                sensor_file = subjec_info.sensor_level(i);
                sensor = load(fullfile(dataset.Path,subjec_info.SubID, sensor_file.Ref_path,sensor_file.Name));
                PlotSensorTopography(fullfile(reportPath,'sensor'),sensor,common.Channels);
                freqs(i) = sensor.band;
            end
             prop.sensor_params  = dataset.sensor_params;
        end
        
        % Activation level
        if(isequal(dataset.general_params.analysis_level.value,'all') ...
                || isequal(dataset.general_params.analysis_level.value,'2') ...
                || isequal(dataset.general_params.analysis_level.value,'12') ...
                || isequal(dataset.general_params.analysis_level.value,'23'))
            for i=1:length(subjec_info.activation_level)
                activation_file = subjec_info.activation_level(i);
                activation = load(fullfile(dataset.Path,subjec_info.SubID, activation_file.Ref_path,activation_file.Name));
                % J3D(:,i) = activation.J;
                PlotSourceTopography(fullfile(reportPath,'activation'),freqs(i),activation,common.Cortex);
            end
             prop.activation_params  = dataset.activation_params;
        end
        % Source PSD
        % PlotSourceSpectrum(fullfile(reportPath,'activation'), dataset, J3D);

       
        createAdvancedReport(userData, reportPath, 'Subject', prop);


    case 'dataset'

    case 'Sensor'

    case 'Activation'

    case 'Connectivity'

end
end

