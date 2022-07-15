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
    count_node = varargin{2};
    if(~isnumeric(idnode) || ~isnumeric(count_node))
        fprintf(2,"\n ->> Error: The selected node and count of nodes have to be numbers \n");
        return;
    end
else
    idnode = 1;
    count_node = 1;
end
disp(strcat("-->> Working in instance: ",num2str(idnode)));
disp('---------------------------------------------------------------------');

%%

addpath(genpath('functions'));
addpath(('external'));
addpath(('external/osl_core'));
addpath(genpath('external/MEG-ROI-nets'));
% process_waitbar = waitbar(0,'Please wait...');
%%
%% Defining parameters
%%
properties                      = get_properties();
if(isequal(properties,'canceled'))
    return;
end
status                          = check_properties(properties);
if(~status)
    fprintf(2,strcat('\nBC-V-->> Error: The current configuration files are wrong \n'));
    disp('Please check the configuration files.');
    return;
end
properties.general_params       = properties.general_params.params;
properties.sensor_params        = properties.sensor_params.params;
properties.activation_params    = properties.activation_params.params;
properties.connectivity_params  = properties.connectivity_params.params;
properties.spectral_params      = properties.spectral_params.params;
root_path                       = properties.general_params.bcv_workspace.BCV_input_dir;
subjects                        = dir(fullfile(root_path,'**','subject.mat'));
predef_folder                   = properties.run_bash_mode.predefinition_params;
%%
%% Multi-node process ( dividing data subjects by each node)
%%
if(length(subjects)< idnode)
    warning("BC-V-->> There is not enough data to run in this node");
    disp(strcat("BC-V-->> Closing process in nodo: ",num2str(idnode),"."));
    return;
end
sub_count                       = fix(length(subjects)/count_node);
rest_sub                        = mod(length(subjects),count_node);
start_ind                       = idnode * sub_count - sub_count + 1;
% end_ind = idnode * sub_count;
if(start_ind>length(subjects))
    fprintf(2,strcat('\nBC-V-->> Error: The follow folder: \n'));
    disp(root_path);
    fprintf(2,strcat('BC-V-->> Error: Does not contain enough data to process on this node.\n'));
    fprintf(2,strcat('BC-V-->> Error: Or Does not contain any subject information file.\n'));
    disp("Please verify the configuration of the input data and start the process again.");
    return;
end
if(~isequal(rest_sub,0))
    if(idnode<=rest_sub)
        start_ind               = start_ind + idnode - 1;
        end_ind                 = start_ind + sub_count;
    else
        start_ind               = start_ind + rest_sub;
        end_ind                 = start_ind + sub_count - 1;
    end
else
    end_ind                     = start_ind + sub_count - 1;
end
if(sub_count == 0)
    start_ind                   = idnode;
    end_ind                     = idnode;
end
subjects                        = subjects([start_ind:end_ind]);

%%
%% Starting subjects analysis
%%
if(~isempty(subjects))
    [properties]                                            = define_frequency_bands(properties);
    color_map                                               = load(properties.general_params.colormap_path);
    properties.cmap                                         = color_map.cmap;
    properties.cmap_a                                       = color_map.cmap_a;
    properties.cmap_c                                       = color_map.cmap_c;
    %% Starting analysis
    for i=1:length(subjects)
        subject_file                                        = subjects(i);
        [subject,checked,error_msg_array]                   = checked_subject_data(subject_file,properties);
        if(checked)
            if(isequal(properties.general_params.bcv_workspace.BCV_work_dir,'local'))
                subject.subject_path                        = fullfile(subject_file.folder,'BC-V_Result');
            else
                subject.subject_path                        = fullfile(properties.general_params.bcv_workspace.BCV_work_dir,subject.name);
            end
            if(~isfolder(subject.subject_path))
                mkdir(subject.subject_path);
            end
            
            %%
            %% Data analysis for sensor level
            %%
            if(isequal(properties.general_params.analysis_level.value,'1')...
                    || isequal(properties.general_params.analysis_level.value,'12')...
                    || isequal(properties.general_params.analysis_level.value,'all'))
                [properties,canceled]                       = check_BC_V_info(properties,subject,1);
                if(canceled)
                    continue;
                end
                %% Saving general variables for sensor level                
                pathname_common                             = fullfile(subject.subject_path,'Common');
                if(~isfolder(pathname_common))
                    mkdir(pathname_common);
                end
                atlas = subject.Scortex.Atlas(subject.Scortex.iAtlas);
                cortex = subject.Scortex;
                
                file_name                                   = strcat('Atlas.mat');
                disp(strcat("File: ", file_name));                                
                parsave(fullfile(pathname_common ,file_name ),atlas);
                reference_path                              = strsplit(pathname_common,subject.name);
                properties.BC_V_info.common.Comment         = 'Common Atlas';
                properties.BC_V_info.common.Ref_path        = strrep(reference_path{2},'\','/');
                properties.BC_V_info.common.Name            = file_name;
                
                file_name                                   = strcat('Cortex.mat');
                disp(strcat("File: ", file_name));                                
                parsave(fullfile(pathname_common ,file_name ),cortex);
                reference_path                              = strsplit(pathname_common,subject.name);
                properties.BC_V_info.common.Comment         = 'Common Cortex';
                properties.BC_V_info.common.Ref_path        = strrep(reference_path{2},'\','/');
                properties.BC_V_info.common.Name            = file_name;               
                
                properties.analysis_level                   = 1;
                if(properties.general_params.run_by_trial.value)
                    %% Data analysis for sensor and activation level by trials
                    if(iscell(subject.data))
                        data = subject.data;
                        for m=1:length(data)
                            properties.trial_name           = ['trial_',num2str(m)];
                            subject.data                    = data{1,m};
                            [subject,properties]            = data_analysis(subject,properties);
                        end
                        subject.data = data;
                    else
                        %% Data analysis for sensor and activation level by complete data
                        [subject,properties]                = data_analysis(subject,properties);
                    end
                else
                    %% Data analysis for sensor and activation level by complete data
                    [subject,properties]                    = data_analysis(subject,properties);
                end
                disp('=================================================================');
                disp('-->> Saving BC-VARETA Information file.')
                BC_V_info                                   = properties.BC_V_info;
                BC_V_info.subjectID                         = subject.name;
                properties.level_process.name               = 'Sensor_level';
                properties.level_process.completed          = true;
                BC_V_info.properties.general_params         = properties.general_params;
                BC_V_info.properties.spectral_params        = properties.spectral_params;
                disp(strcat("File: ", "BC_V_info.mat"));
                save(fullfile(subject.subject_path ,'BC_V_info.mat'),'-struct','BC_V_info');
            end
            %%
            %% Data analysis for activation level
            %%
            if(isequal(properties.general_params.analysis_level.value,'2')...
                    || isequal(properties.general_params.analysis_level.value,'12')...
                    || isequal(properties.general_params.analysis_level.value,'23')...
                    || isequal(properties.general_params.analysis_level.value,'all'))
                [properties,canceled]                   = check_BC_V_info(properties,subject,2);
                if(canceled)
                    continue;
                end                
                properties.analysis_level               = 2;
                if(properties.general_params.run_by_trial.value)
                    if((isfield(properties.BC_V_info,'sensor_level')))
                        %% Data analysis for sensor and activation level by trials
                        if(iscell(subject.data))
                            data = subject.data;
                            for m=1:length(data)
                                properties.trial_name   = ['trial_',num2str(m)];
                                subject.data            = data{1,m};
                                [subject,properties]    = get_activation_priors(subject,properties);
                                [subject,properties]    = data_analysis(subject,properties);
                            end
                            subject.data = data;
                        else
                            %% Data analysis for sensor and activation level by complete data
                            [subject,properties]        = get_activation_priors(subject,properties);
                            [subject,properties]        = data_analysis(subject,properties);
                        end
                    else
                        fprintf(2,strcat('\nBC-V-->> Error: Do not process activation level for subject: \n'));
                        disp(subject.name);
                        fprintf(2,strcat('BC-V-->> Error: This subject do not countain the sensor process output.\n'));
                        disp("Please, run first the sensor process.");
                        continue;
                    end
                else
                    if((isfield(properties.BC_V_info,'sensor_level')))
                        %% Data analysis for sensor and activation level by complete data
                        [subject,properties]            = get_activation_priors(subject,properties);
                        [subject,properties]            = data_analysis(subject,properties);
                    else
                        fprintf(2,strcat('\nBC-V-->> Error: Do not process activation level for subject: \n'));
                        disp(subject.name);
                        fprintf(2,strcat('BC-V-->> Error: This subject do not countain the sensor process output.\n'));
                        disp("Please, run first the sensor process.");
                        continue;
                    end
                end
                disp('-->> Saving file.')
                BC_V_info                               = properties.BC_V_info;
                properties.level_process.name           = 'Activation_level';
                properties.level_process.completed      = true;
                BC_V_info.properties.activation_params  = properties.activation_params; 
                disp(strcat("File: ", "BC_V_info.mat"));
                save(fullfile(subject.subject_path ,'BC_V_info.mat'),'-struct','BC_V_info');
            end
            %%
            %% Data analysis for connectivity level
            %%
            if(isequal(properties.general_params.analysis_level.value,'3')...
                    || isequal(properties.general_params.analysis_level.value,'23')...
                    || isequal(properties.general_params.analysis_level.value,'all'))
                [properties,canceled]                   = check_BC_V_info(properties,subject,3);
                if(canceled)
                    continue;
                end                
                properties.analysis_level               = 3;
                if(properties.general_params.run_by_trial.value)
                    if(isfield(properties.BC_V_info,'sensor_level') && isfield(properties.BC_V_info,'activation_level'))
                        %% Data analysis for connectivity level by trials
                        if(iscell(subject.data))
                            data = subject.data;
                            for m=1:length(data)
                                properties.trial_name   = ['trial_',num2str(m)];
                                subject.data            = data{1,m};
                                [subject,properties]    = get_connectivity_priors(subject,properties);
                                [subject,properties]    = data_analysis(subject,properties);
                            end
                            subject.data = data;
                        else
                            %% Data analysis for connectivity level by complete data
                            [subject,properties]        = get_connectivity_priors(subject,properties);
                            [subject,properties]        = data_analysis(subject,properties);
                        end
                    else
                        fprintf(2,strcat('\nBC-V-->> Error: Do not process connectivity level for subject: \n'));
                        disp(subject.name);
                        fprintf(2,strcat('BC-V-->> Error: This subject do not countain the sensor and activation process output.\n'));
                        disp("Please, run first the sensor and activation process.");
                        continue;
                    end
                else
                    if((isfield(properties.BC_V_info,'sensor_level') && isfield(properties.BC_V_info,'activation_level')))
                        %% Data analysis for connectivity level by complete data
                        [subject,properties]                = get_connectivity_priors(subject,properties);
                        [subject,properties]                = data_analysis(subject,properties);
                    else
                        fprintf(2,strcat('\nBC-V-->> Error: Do not process connectivity level for subject: \n'));
                        disp(subject.name);
                        fprintf(2,strcat('BC-V-->> Error: This subject do not countain the sensor and activation process output.\n'));
                        disp("Please, run first the sensor and activation process.");
                        continue;
                    end
                end
                disp('-->> Saving file.')
                BC_V_info                                   = properties.BC_V_info;
                properties.level_process.name               = 'Connectivity_level';
                properties.level_process.completed          = true;
                BC_V_info.properties.connectivity_params    = properties.connectivity_params;               
                disp(strcat("File: ", "BC_V_info.mat"));
                save(fullfile(subject.subject_path ,'BC_V_info.mat'),'-struct','BC_V_info');
            end
        else
            fprintf(2,strcat('\nBC-V-->> Error: The folder structure for subject: ',subject.name,' \n'));
            fprintf(2,strcat('BC-V-->> Have the folows errors.\n'));
            for j=1:length(error_msg_array)
                fprintf(2,strcat('BC-V-->>' ,error_msg_array(j), '.\n'));
            end
            fprintf(2,strcat('BC-V-->> Jump to an other subject.\n'));
            continue;
        end
    end
else
    fprintf(2,strcat('\nBC-V-->> Error: The folder structure: \n'));
    disp(root_path);
    fprintf(2,strcat('BC-V-->> Error: Do not contain any subject information file.\n'));
    disp("Please verify the configuration of the input data and start the process again.");
    return;
end

end