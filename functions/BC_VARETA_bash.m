function BC_VARETA_bash(varargin)
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
    idnode = varargin{1};
    total_node = varargin{2};
    if(~isnumeric(idnode) || ~isnumeric(total_node))
        fprintf(2,"\n ->> Error: The selected node and count of nodes have to be numbers \n");
        return;
    end
else
    idnode = 1;
    total_node = 1;
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
properties                      = get_properties();
if(isequal(properties,'canceled'))
    return;
end
[status,reject_subjects]        = check_properties(properties);
if(~status)
    fprintf(2,strcat('\nBC-V-->> Error: The current configuration files are wrong \n'));
    disp('Please check the configuration files.');
    return;
end
properties.general_params       = properties.general_params.params;
properties.sensor_params        = properties.sensor_params.params;
properties.activation_params    = properties.activation_params.params;
properties.connectivity_params  = properties.connectivity_params.params;
root_path                       = properties.general_params.bcv_workspace.BCV_input_dir;
subjects                        = dir(fullfile(root_path,'**','subject.mat'));

%%
%% Multi-node process ( splitting data subjects by each node)
%%
subjects                        = multinode_subjects(subjects,idnode,total_node);

%%
%% Starting subjects analysis
%%
if(~isempty(subjects))
    [properties]                                                = define_frequency_bands(properties);
    color_map                                                   = load(properties.general_params.colormap_path);
    properties.cmap                                             = color_map.cmap;
    properties.cmap_a                                           = color_map.cmap_a;
    properties.cmap_c                                           = color_map.cmap_c;
    
    %% Creating Dataset
    BC_VARETA.Name                                              = properties.general_params.dataset.Name;
    BC_VARETA.Description                                       = properties.general_params.dataset.Description;
    BC_VARETA.general_params                                    = properties.general_params;
    BC_VARETA.sensor_params                                     = properties.sensor_params;
    BC_VARETA.activation_params                                 = properties.activation_params;
    BC_VARETA.connectivity_params                               = properties.connectivity_params;
    BC_VARETA.Participants                                      = [];
    
    %%
    %% Starting analysis
    %%
    for i=1:length(subjects)
        subject_file                                            = subjects(i);
        [subject,checked,error_msg_array]                       = checked_subject_data(subject_file,properties);
        if(checked)
            subject.subject_path                                = fullfile(properties.general_params.bcv_workspace.BCV_work_dir,subject.name);
            if(~isfolder(subject.subject_path))
                mkdir(subject.subject_path);
            end         

            %%
            %% Data analysis for sensor level
            %%
            if(isequal(properties.general_params.analysis_level.value,'1')...
                    || isequal(properties.general_params.analysis_level.value,'12')...
                    || isequal(properties.general_params.analysis_level.value,'all'))
                [subject,status]                                = check_BC_V_info(properties,subject,1);
                if(status)
                    % Saving data
                    subject                                     = BC_V_save(properties,subject,'common');                    
                   
                    if(properties.general_params.run_by_trial.value)
                        data                                    = subject.MEEG.data;
                        for m=1:length(data)
                            properties.trial_name               = ['trial_',num2str(m)];
                            subject.MEEG.data                   = data{1,m};
                            [subject,properties]                = sensor_level_analysis(subject,properties);
                        end
                        subject.MEEG.data                       = data;
                    else
                        [subject,properties]                    = sensor_level_analysis(subject,properties);
                    end
                    disp('=================================================================');
                    subject                                     = BC_V_save(properties,subject,'level1');
                    
                end
            end
            %%
            %% Data analysis for activation level
            %%
            if(isequal(properties.general_params.analysis_level.value,'2')...
                    || isequal(properties.general_params.analysis_level.value,'12')...
                    || isequal(properties.general_params.analysis_level.value,'23')...
                    || isequal(properties.general_params.analysis_level.value,'all'))
                [subject,status]                                = check_BC_V_info(properties,subject,2);
                if(status)
                    if((isfield(subject.BC_V_info,'sensor_level')))
                        if(properties.general_params.run_by_trial.value)
                            data                                = subject.data;
                            for m=1:length(data)
                                properties.trial_name           = ['trial_',num2str(m)];
                                subject.data                    = data{1,m};
                                [subject,properties]            = activation_level_interface(subject,properties);
                            end
                            subject.data                        = data;
                        else
                            [subject,properties]                = activation_level_interface(subject,properties);
                        end
                        disp('=================================================================');
                        subject                                 = BC_V_save(properties,subject,'level2');
                        
                    else
                        fprintf(2,strcat('\nBC-V-->> Error: Do not process activation level for subject: \n'));
                        disp(subject.name);
                        fprintf(2,strcat('BC-V-->> Error: This subject do not countain the sensor process output.\n'));
                        disp("Please, run first the sensor process.");
                    end                    
                end
            end
            %%
            %% Data analysis for connectivity level
            %%
            if(isequal(properties.general_params.analysis_level.value,'3')...
                    || isequal(properties.general_params.analysis_level.value,'23')...
                    || isequal(properties.general_params.analysis_level.value,'all'))
                [subject,status]                                = check_BC_V_info(properties,subject,3);
                if(status)
                    if(properties.general_params.run_by_trial.value)
                        data                                    = subject.data;
                        for m=1:length(data)
                            properties.trial_name               = ['trial_',num2str(m)];
                            subject.data                        = data{1,m};
                            [subject,properties]                = connectivity_level_interface(subject,properties);
                        end
                        subject.data                            = data;
                    else
                        [subject,properties]                    = connectivity_level_interface(subject,properties);
                    end 
                    disp('=================================================================');
                    subject                                     = BC_V_save(properties,subject,'level3');                    
                else
                    fprintf(2,strcat('\nBC-V-->> Error: Do not process connectivity level for subject: \n'));
                    disp(subject.name);
                    fprintf(2,strcat('BC-V-->> Error: This subject do not countain the sensor and activation process output.\n'));
                    disp("Please, run first the sensor and activation process.");
                end
            end
            BC_VARETA.Participants(i).SubID                     = subject.name;
            BC_VARETA.Participants(i).Status                    = "Completed";
            BC_VARETA.Participants(i).FileInfo                  = "BC_V_info.mat";
            BC_VARETA.Participants(i).Error                     = [];
        else
            BC_VARETA.Participants(i).SubID                     = subject.name;
            BC_VARETA.Participants(i).Status                    = "Rejected";
            BC_VARETA.Participants(i).FileInfo                  = "";
            fprintf(2,strcat('\nBC-V-->> Error: The folder structure for subject: ',subject.name,' \n'));
            fprintf(2,strcat('BC-V-->> Have the folows errors.\n'));
            for j=1:length(error_msg_array)
                fprintf(2,strcat('BC-V-->>' ,error_msg_array(j), '.\n'));
                BC_VARETA.Participants(i).Error(j)              = error_msg_array(j);
            end
            fprintf(2,strcat('BC-V-->> Jump to an other subject.\n'));
            continue;
        end        
    end
    save(fullfile(properties.general_params.bcv_workspace.BCV_work_dir ,'BC_VARETA.mat'),'-struct','BC_VARETA');
else
    fprintf(2,strcat('\nBC-V-->> Error: The folder structure: \n'));
    disp(root_path);
    fprintf(2,strcat('BC-V-->> Error: Do not contain any subject information file.\n'));
    disp("Please verify the configuration of the input data and start the process again.");
    return;
end
end
