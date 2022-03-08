function [properties,canceled] = check_BC_V_info(properties,subject,level)
canceled = false;
subject_path = subject.subject_path;
%%
%% Get information file for previous analysis
%%
BC_V_file = dir(fullfile(subject_path ,'BC_V_info.mat'));
if(~isempty(BC_V_file))
    BC_V_info = load(fullfile(BC_V_file.folder,BC_V_file.name));
    reseted = false;
    if(~isequal(properties.general_params.run_by_trial.value,BC_V_info.properties.general_params.run_by_trial.value))
        BC_V_info = struct;
        files = dir(subject_path);
        dirFlags = [files.isdir];
        subFolders = files(dirFlags);
        for k =3:length(subFolders)            
                rmdir(fullfile(subFolders(k).folder,subFolders(k).name),'s');            
        end
        reseted = true;
    end
    if(~reseted && ~isequal(properties.general_params.run_frequency_bin.value,BC_V_info.properties.general_params.run_frequency_bin.value))
        BC_V_info = struct;
        files = dir(subject_path);
        dirFlags = [files.isdir];
        subFolders = files(dirFlags);
        for k =3:length(subFolders)
            rmdir(fullfile(subFolders(k).folder,subFolders(k).name),'s');
        end
        reseted = true;
    end
    properties.BC_V_info = BC_V_info;
end
if(isequal(level,1))
    if(isempty(BC_V_file) || properties.general_params.analysis_level.reset_all)
        BC_V_info = struct;
        files = dir(subject_path);
        dirFlags = [files.isdir];
        subFolders = files(dirFlags);
        for k =3:length(subFolders)
            rmdir(fullfile(subFolders(k).folder,subFolders(k).name),'s');
        end
    else
        
    end
    properties.BC_V_info = BC_V_info;
end
if(isequal(level,2))
    if(isempty(BC_V_file))
        fprintf(2,strcat('\nBC-V-->> Error: The BC-V info file is not correct for subject: ',subject.name,' \n'));
        fprintf(2,strcat('BC-V-->> You need to exjecute the sensor lever first.\n'));
        fprintf(2,strcat('BC-V-->> Jump to an other subject.\n'));
        canceled = true;
    else
        if(properties.general_params.analysis_level.reset_all)
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
        else
            for m=1:length(BC_V_info.properties.spectral_params.frequencies)
                BC_V_info.properties.spectral_params.frequencies(m).run = (BC_V_info.properties.spectral_params.frequencies(m).run || properties.spectral_params.frequencies(m).run);
            end
            if(isfield(BC_V_info.properties,'activation_params'))
                for m=1:length(BC_V_info.properties.activation_params.methods)
                    method_name = fieldnames(BC_V_info.properties.activation_params.methods{m});
                    method_name = method_name{1};
                    BC_V_info.properties.activation_params.methods{m}.(method_name).run = (BC_V_info.properties.activation_params.methods{m}.(method_name).run || properties.activation_params.methods{m}.(method_name).run);
                end
            end
        end
    end
    properties.BC_V_info = BC_V_info;
end
if(isequal(level,3))
    if(isempty(BC_V_file))
        fprintf(2,strcat('\nBC-V-->> Error: The BC-V info file is not correct for subject: ',subject.name,' \n'));
        fprintf(2,strcat('BC-V-->> You need to exjecute the activation lever first.\n'));
        fprintf(2,strcat('BC-V-->> Jump to an other subject.\n'));
        canceled = true;
    else
        if(properties.general_params.analysis_level.reset_all)
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
        else
            for m=1:length(BC_V_info.properties.spectral_params.frequencies)
                BC_V_info.properties.spectral_params.frequencies(m).run = (BC_V_info.properties.spectral_params.frequencies(m).run || properties.spectral_params.frequencies(m).run);
            end
            if(isfield(BC_V_info.properties,'connectivity_params'))
                for m=1:length(BC_V_info.properties.connectivity_params.methods)
                    method_name = fieldnames(BC_V_info.properties.connectivity_params.methods{m});
                    method_name = method_name{1};
                    BC_V_info.properties.connectivity_params.methods{m}.(method_name).run = (BC_V_info.properties.connectivity_params.methods{m}.(method_name).run || properties.connectivity_params.methods{m}.(method_name).run);
                end
            end
        end
    end
    properties.BC_V_info = BC_V_info;
end

end


