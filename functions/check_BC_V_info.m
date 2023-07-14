function [subject,status] = check_BC_V_info(properties,subject,level)
status = true;
subject_path = subject.subject_path;
%%
%% Get information file for previous analysis
%%
if(isequal(level,1))
    BC_V_info = struct;
    BC_V_info.SubID = subject.name;
    files = dir(subject_path);
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    for k =3:length(subFolders)
        rmdir(fullfile(subFolders(k).folder,subFolders(k).name),'s');
    end
end
if(isequal(level,2))
    BC_V_file = dir(fullfile(subject_path ,'BC_V_info.mat'));
    if(isempty(BC_V_file))
        fprintf(2,strcat('\nBC-V-->> Error: The BC-V info file is not correct for subject: ',subject.name,' \n'));
        fprintf(2,strcat('BC-V-->> You need to exjecute the sensor lever first.\n'));
        fprintf(2,strcat('BC-V-->> Jump to an other subject.\n'));
        status = false;
        return;
    end
    BC_V_info = load(fullfile(BC_V_file.folder,BC_V_file.name));
    if((~isfield(BC_V_info,'sensor_level')))
        fprintf(2,strcat('\nBC-V-->> Error: Do not process activation level for subject: \n'));
        disp(subject.name);
        fprintf(2,strcat('BC-V-->> Error: This subject do not countain the sensor process output.\n'));
        disp("Please, run first the sensor process.");
        status = false;
        return;
    end    
    if(properties.general_params.run_by_trial.value)
        trials = fieldnames(BC_V_info);
        for m=1:length(trials)
            trial_name = trials{m};
            trial = BC_V_info.(trial_name);
            if( isfield(trial,'activation_level'))
                BC_V_info.(trial_name) = rmfield(trial,'activation_level');
            end
            if( isfield(trial,'connectivity_level'))
                BC_V_info.(trial_name) = rmfield(trial,'connectivity_level');
            end
            if( exist(fullfile(subject_path,trial_name,'Activation_level'),'dir'))
                rmdir(fullfile(subject_path,trial_name,'Activation_level'),'s');
            end
            if( exist(fullfile(subject_path,trial_name,'Connectivity_level'),'dir'))
                rmdir(fullfile(subject_path,trial_name,'Connectivity_level'),'s');
            end
        end
    else
        if( isfield(BC_V_info,'activation_level'))
            BC_V_info = rmfield(BC_V_info,'activation_level');
        end
        if( isfield(BC_V_info,'connectivity_level'))
            BC_V_info = rmfield(BC_V_info,'connectivity_level');
        end
        if( exist(fullfile(subject_path,'Activation_level'),'dir'))
            rmdir(fullfile(subject_path,'Activation_level'),'s');
        end
        if( exist(fullfile(subject_path,'Connectivity_level'),'dir'))
            rmdir(fullfile(subject_path,'Connectivity_level'),'s');
        end
    end
end
if(isequal(level,3))
    BC_V_file = dir(fullfile(subject_path ,'BC_V_info.mat'));
    if(isempty(BC_V_file))
        fprintf(2,strcat('\nBC-V-->> Error: The BC-V info file is not correct for subject: ',subject.name,' \n'));
        fprintf(2,strcat('BC-V-->> You need to exjecute the activation lever first.\n'));
        fprintf(2,strcat('BC-V-->> Jump to an other subject.\n'));
        status = false;
        return;
    end
    if((~isfield(BC_V_info,'sensor_level')))
        fprintf(2,strcat('\nBC-V-->> Error: Do not process activation level for subject: \n'));
        disp(subject.name);
        fprintf(2,strcat('BC-V-->> Error: This subject do not countain the sensor process output.\n'));
        disp("Please, run first the sensor process.");
        status = false;
        return;
    end
    if((~isfield(BC_V_info,'activation_level')))
        fprintf(2,strcat('\nBC-V-->> Error: Do not process activation level for subject: \n'));
        disp(subject.name);
        fprintf(2,strcat('BC-V-->> Error: This subject do not countain the activation process output.\n'));
        disp("Please, run first the activation process.");
        status = false;
        return;
    end
    BC_V_info = load(fullfile(BC_V_file.folder,BC_V_file.name));
    if(properties.general_params.run_by_trial.value)
        trials = fieldnames(BC_V_info);
        for m=1:length(trials)
            trial_name = trials{m};
            trial = BC_V_info.(trial_name);
            if( isfield(trial,'connectivity_level'))
                BC_V_info.(trial_name) = rmfield(trial,'connectivity_level');
            end
            if( exist(fullfile(subject_path,trial_name,'Connectivity_level'),'dir'))
                rmdir(fullfile(subject_path,trial_name,'Connectivity_level'),'s');
            end
        end
    else
        if( isfield(BC_V_info,'connectivity_level'))
            BC_V_info = rmfield(BC_V_info,'connectivity_level');
        end
        if( exist(fullfile(subject_path,'Connectivity_level'),'dir'))
            rmdir(fullfile(subject_path,'Connectivity_level'),'s');
        end
    end
end
subject.BC_V_info = BC_V_info;
end


