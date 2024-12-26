function result = report_process_interface(dataset,varargin)

for k = 2 : nargin
    eval([inputname(k) '=  varargin{k-1};']);
end


result = "success";
reportPath = fullfile(pwd,'tmp','report','subject');
if(~isfolder(reportPath))
    mkdir(reportPath);
end


subjec_info = nodeData;
disp(strcat("-->> Exporting report for subject: ", subjec_info.SubID));
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
        || isequal(dataset.general_params.analysis_level.value,'12') ...
        && (isequal(clasif,'subject') || isequal(clasif,'sensor_level')))
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
        || isequal(dataset.general_params.analysis_level.value,'23') ...
        && (isequal(clasif,'subject') || isequal(clasif,'activation_level')))
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


%% 
%% Defining Participant information
%%
if(isfield(userData,'Name'))
    Partic_info.Name = userData.Name;
else
        Partic_info.Name = 'Anonymized';
end
Partic_info.SubID = userData.SubID;
if(isfield(userData, 'Sex'))
    Partic_info.Gender = userData.Sex;
elseif(isfield(userData, 'Gender'))
    Partic_info.Gender = userData.Gender;
else
    Partic_info.Gender = 'Undefinited';
end
if(isfield(userData,'Age'))
    Partic_info.Age = userData.Age;
else
    Partic_info.Age = 'Undefinited';
end
if(isfield(userData,'Condition'))
        Partic_info.Condition =  userData.Condition;
elseif(isfield(userData,'Task'))
    Partic_info.Condition =  userData.Task;
else
     Partic_info.Condition = 'Undefinited';
end

outputFile = createAdvancedReport(Partic_info, reportPath, 'Subject', prop);
movefile(outputFile,fullfile(dataset.Path, subjec_info.SubID),"f");
rmdir(reportPath, 's');

end

