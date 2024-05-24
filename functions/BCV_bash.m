function BCV_bash(varargin)
%% BC_VARETA_bash Summary of this function goes here
%   Detailed explanation goes here
%
%
% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa
%
% Date: March 18, 2019
%
% Updates:
% - Ariosky Areces Gonzalez
%
% Date: March 26, 2019
%%

if(isequal(nargin,2))
    idnode      = varargin{1};
    total_node  = varargin{2};
    if(~isnumeric(idnode) || ~isnumeric(total_node))
        fprintf(2,"\n ->> Error: The selected node and count of nodes have to be numbers \n");
        return;
    end
else
    idnode      = 1;
    total_node  = 1;
end
disp(strcat("-->> Working in instance: ",num2str(idnode)));
disp('---------------------------------------------------------------------');

%% Adding paths
addpath(genpath('functions'));
addpath(('external'));
addpath(('external/osl_core'));
addpath(genpath('external/MEG-ROI-nets'));

%%
%% Defining parameters
%%
properties                  = get_properties();
if(isequal(properties,'canceled'))
    return;
end
[properties,status,reject]  = check_properties(properties);
if(~status)
    fprintf(2,strcat('\nBC-V-->> Error: The current configuration files are wrong \n'));
    disp('Please check the configuration files.');
    return;
end

root_path                   = properties.general_params.bcv_workspace.BCV_input_dir;
CiftiStorm                  = jsondecode(fileread(fullfile(root_path,'CiftiStorm.json')));
if(isfield(CiftiStorm,'Template'))
    sufix = CiftiStorm.Template.SubID;
end
subjects                    = CiftiStorm.Participants;

%%
%% Multi-node process ( splitting data subjects by each node)
%%
subjects                    = multinode_subjects(subjects,idnode,total_node);

%%
%% Starting subjects analysis
%%
if(isempty(subjects))
    fprintf(2,strcat('\nBC-V-->> Error: The folder structure: \n'));
    disp(root_path);
    fprintf(2,strcat('BC-V-->> Error: Do not contain any subject information file.\n'));
    disp("Please verify the configuration of the input data and start the process again.");
    return;
end
[properties]                                = define_frequency_bands(properties);
analysis_level                              = properties.general_params.analysis_level.value;

%%
%% Creating Dataset
%%
if(isequal(total_node,1))
    BCV_filename = 'BC_VARETA.json';
else
    BCV_filename = strcat('BC_VARETA_node-',num2str(idnode),'.json');
end
if(exist('sufix','var'))
    folder = strcat('bcvareta-',sufix);
else
    folder = 'bcvareta';
end
BCV_file = fullfile(properties.general_params.bcv_workspace.BCV_work_dir,folder,properties.general_params.dataset.task.value,BCV_filename);
if(isfile(BCV_file))
    BC_VARETA = jsondecode(fileread(BCV_file));
else
    BC_VARETA.Name                              = properties.general_params.dataset.Name;
    BC_VARETA.Description                       = properties.general_params.dataset.Description;
    BC_VARETA.Task                              = properties.general_params.dataset.task.value;
    BC_VARETA.Status                            = "Processing";
    BC_VARETA.general_params                    = properties.general_params;
    BC_VARETA.sensor_params                     = properties.sensor_params;
    BC_VARETA.activation_params                 = properties.activation_params;
    BC_VARETA.connectivity_params               = properties.connectivity_params;
    BC_VARETA.Participants                      = [];
end

%%
%% Starting analysis
%%
for i=1:length(subjects)
    if(isfile(BCV_file))
        BC_VARETA = jsondecode(fileread(BCV_file));
    end
    participant                            = subjects(i);
    if(~isempty(BC_VARETA.Participants))
        iPart = find(ismember({BC_VARETA.Participants.SubID},participant.SubID),1);
        if(~isempty(iPart) && isequal(BC_VARETA.Participants(iPart).Status,"Completed"))
            continue;
        end
    end
    [subject,checked,error_msg_array]       = checked_subject_data(participant,properties);
    if(~checked)
        BC_VARETA.Participants(i).SubID     = subject.name;
        BC_VARETA.Participants(i).Status    = "Rejected";
        BC_VARETA.Participants(i).FileInfo  = "";
        fprintf(2,strcat('\nBC-V-->> Error: The folder structure for subject: ',subject.name,' \n'));
        fprintf(2,strcat('BC-V-->> Have the folows errors.\n'));
        for j=1:length(error_msg_array)
            fprintf(2,strcat('BC-V-->>' ,error_msg_array(j), '.\n'));
            BC_VARETA.Participants(i).Error(j)  = error_msg_array(j);
        end
        fprintf(2,strcat('BC-V-->> Jump to an other subject.\n'));
        continue;
    end
   
    subject.subject_path                    = fullfile(properties.general_params.bcv_workspace.BCV_work_dir,folder,BC_VARETA.Task,subject.name);
    if(~isfolder(subject.subject_path))
        mkdir(subject.subject_path);
    end

    %%
    %% Data analysis for sensor level
    %%
    if(ismember(analysis_level,{'1','12','all'}))
        [subject,status]                    = check_BC_V_info(properties,subject,1);
        if(status)
            subject                         = get_structural_priors(subject);
            subject                         = BC_V_save(properties,subject,'common');
            if(properties.general_params.run_by_trial.value)
                data                        = subject.MEEG.data;
                for m=1:length(data)
                    properties.trial_name   = ['trial_',num2str(m)];
                    disp("=====================================================================");
                    disp(strcat("BC-V-->> Processing trial: ",properties.trial_name));
                    disp("=====================================================================");
                    subject.MEEG.data       = data{1,m};
                    [subject,properties]    = sensor_level_analysis(subject,properties);
                    disp("=====================================================================");
                    subject                         = BC_V_save(properties,subject,'level1');
                end
                subject.MEEG.data           = data;
            else
                [subject,properties]        = sensor_level_analysis(subject,properties);
                disp("=====================================================================");
                subject                     = BC_V_save(properties,subject,'level1');
            end
            
        end
    end

    %%
    %% Data analysis for activation level
    %%
    if(ismember(analysis_level,{'2','12','23','all'}))
        [subject,status]                    = check_BC_V_info(properties,subject,2);
        if(status)
            if(properties.general_params.run_by_trial.value)
                data                        = subject.MEEG.data;
                for m=1:length(data)
                    properties.trial_name   = ['trial_',num2str(m)];
                    disp("=====================================================================");
                    disp(strcat("BC-V-->> Processing trial: ",properties.trial_name));
                    disp("=====================================================================");
                    subject.data            = data{1,m};
                    [subject,properties]    = activation_level_interface(subject,properties);
                end
                subject.MEEG.data                = data;
            else
                [subject,properties]        = activation_level_interface(subject,properties);
            end
            disp("=====================================================================");
            subject                         = BC_V_save(properties,subject,'level2');
        end
    end

    %%
    %% Data analysis for connectivity level
    %%
    if(ismember(analysis_level,{'3','23','all'}))
        [subject,status]                    = check_BC_V_info(properties,subject,3);
        if(status)
            if(properties.general_params.run_by_trial.value)
                data                        = subject.MEEG.data;
                for m=1:length(data)
                    properties.trial_name   = ['trial_',num2str(m)];
                    disp("=====================================================================");
                    disp(strcat("BC-V-->> Processing trial: ",properties.trial_name));
                    disp("=====================================================================");
                    subject.data            = data{1,m};
                    [subject,properties]    = connectivity_level_interface(subject,properties);
                end
                subject.MEEG.data                = data;
            else
                [subject,properties]        = connectivity_level_interface(subject,properties);
            end
            disp("=====================================================================");
            subject                         = BC_V_save(properties,subject,'level3');
        end
    end

    %%
    %% Subject Process status
    %%
    BC_VARETA.Participants(i).SubID         = subject.name;
    BC_VARETA.Participants(i).Status        = "Completed";
    BC_VARETA.Participants(i).FileInfo      = strcat(subject.name,".json");
    BC_VARETA.Participants(i).Error         = [];
    saveJSON(BC_VARETA,BCV_file);
end

%%
%% Saving BC-VARETA file
%%
BC_VARETA = jsondecode(fileread(BCV_file));
BC_VARETA.Status = "Completed";
saveJSON(BC_VARETA,BCV_file);

%%
%% Checking Multi-instance processing
%%
BC_VARETA_files = dir(fullfile(properties.general_params.bcv_workspace.BCV_work_dir,folder,BC_VARETA.Task));
BC_VARETA_files([BC_VARETA_files.isdir]==1) = [];
if(isequal(length(BC_VARETA_files),total_node) && ~isequal(total_node,1))
    complete = true;
    BC_VARETA = jsondecode(fileread(fullfile(BC_VARETA_files(1).folder,BC_VARETA_files(1).name)));
    if(~isequal(BC_VARETA.Status,"Completed"))
        return;
    end
    for i=2:length(BC_VARETA_files)
        BC_VARETA_part = jsondecode(fileread(fullfile(BC_VARETA_files(i).folder,BC_VARETA_files(i).name)));
        if(isequal(BC_VARETA_part.Status,"Completed"))
            Participants = BC_VARETA_part.Participants;
            BC_VARETA.Participants = [BC_VARETA.Participants; Participants];
        else
            complete = false;
            break;
        end
    end
    if(complete)
        BCV_file = fullfile(properties.general_params.bcv_workspace.BCV_work_dir,folder,BC_VARETA.Task,'BC_VARETA.json');
        saveJSON(BC_VARETA,BCV_file);
    end
end

%%
%% Create BC-VARETA Dataset and import
%%
BCV_file = fullfile(properties.general_params.bcv_workspace.BCV_work_dir,folder,BC_VARETA.Task,'BC_VARETA.json');
if(isfile(BCV_file))
    % Loading existed Datasets
    homedir = char(java.lang.System.getProperty('user.home'));
    BCVdir  = fullfile(homedir,"BC_VARETA");
    if(~isfolder(BCVdir))
        mkdir(BCVdir);
    end
    Datasets_file = fullfile(BCVdir,"Datasets.json");
    if(isfile(Datasets_file))
        TempDatasets = jsondecode(fileread(Datasets_file));
        if(~isempty(TempDatasets))
            Datasets = TempDatasets;
        else
            Datasets = struct([]);
        end
    else
        Datasets = struct([]);
    end

    % Including new dataset
    dataset_file            = fullfile(properties.general_params.bcv_workspace.BCV_work_dir,folder,BC_VARETA.Task,BCV_filename);
    BC_VARETA_info          = jsondecode(fileread(dataset_file));
    BC_VARETA_info.Path     = fullfile(properties.general_params.bcv_workspace.BCV_work_dir,folder,BC_VARETA_info.Task);
    TempUUID                = java.util.UUID.randomUUID;
    BC_VARETA_info.UUID     = char(TempUUID.toString);
    if(isempty(Datasets))
        Datasets            = jsondecode(fileread(BCV_file));
    else
        Datasets(end + 1)   = jsondecode(fileread(BCV_file));
    end
    disp("-->> Dataset saved");
    saveJSON(Datasets,Datasets_file);
end
end
