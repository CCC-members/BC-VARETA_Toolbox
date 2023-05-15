function [subject,properties] = activation_level_interface(subject,properties)

%% Get Activation priors
[subject,properties]    = get_activation_priors(subject,properties);

%% Starting Activation Analysis
file_name                                       = strcat('W.mat');
W                                               = subject.W;

if(properties.general_params.run_by_trial.value)
    trial_name = properties.trial_name;
    pathname                                    = fullfile(subject.subject_path,trial_name,'Generals','Structural','sSSBL');
    reference_path                              = strsplit(pathname,subject.name);
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    sensor_level                                = subject.BC_V_info.sensor_level(contains({subject.BC_V_info.sensor_level.Ref_path},properties.trial_name));
else
    pathname                                    = fullfile(subject.subject_path,'Generals','Structural','sSSBL');
    reference_path                              = strsplit(pathname,subject.name);
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    sensor_level                                = subject.BC_V_info.sensor_level;
end
subject.BC_V_info.generals(2).Comment        = 'Generals';
subject.BC_V_info.generals(2).Ref_path       = strrep(reference_path{2},'\','/');
subject.BC_V_info.generals(2).Name           = file_name;
disp(strcat("File: ", file_name));
if(getGlobalGuimode)
    dlg = msgbox('Save operation in progress...');
end
parsave(fullfile(pathname ,file_name ),W);
if ishghandle(dlg)
    delete(dlg);
end

for f=1:length(sensor_level)
    ref_path                                    = sensor_level(f).Ref_path;
    file_name                                   = sensor_level(f).Name;
    subject.sensor_level_out                    = load(fullfile(subject.subject_path,ref_path,file_name));
    band                                        = subject.sensor_level_out.band;
    %%
    %% Defining path
    %%
    disp('=================================================================');
    disp(strcat("BC-V-->> Activation level for frequency band:", band.str_band));
    text_level = 'Activation_level';

    %%
    %% Band Analysis, activation level
    %%
    for m=1:length(properties.activation_params.methods)
        analysis_method = properties.activation_params.methods{m};
        fields          = fieldnames(analysis_method);
        method_name     = fields{1};
        if(analysis_method.(method_name).run)
            disp('-----------------------------------------------------------------');
            disp(strcat("-->> Start time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
            switch method_name
                case 'sssblpp'
                    if(properties.general_params.run_by_trial.value)
                        trial_name                                                  = properties.trial_name;
                        subject.pathname                                         = fullfile(subject.subject_path,trial_name,text_level,'sSSBLpp',band.name);
                    else
                        subject.pathname                                         = fullfile(subject.subject_path,text_level,'sSSBLpp',band.name);
                    end
                    if(~isfolder(subject.pathname))
                        mkdir(subject.pathname);
                    end
                    properties.activation_params.sssblpp_th                         = analysis_method.(method_name).sssblpp_th;
                    [subject,properties]                                     = activation_level_sssblpp(subject,properties);
                case 'eloreta'
                    if(properties.general_params.run_by_trial.value)
                        trial_name                                                  = properties.trial_name;
                        subject.pathname                                         = fullfile(subject.subject_path,trial_name,text_level,'eLORETA',band.name);
                    else
                        subject.pathname                                         = fullfile(subject.subject_path,text_level,'eLORETA',band.name);
                    end
                    if(~isfolder(subject.pathname))
                        mkdir(subject.pathname);
                    end
                    properties.activation_params.gamma1                             = analysis_method.(method_name).gamma1;
                    properties.activation_params.gamma2                             = analysis_method.(method_name).gamma2;
                    properties.activation_params.delta_gamma                        = analysis_method.(method_name).delta_gamma;
                    properties.activation_params.eloreta_th                         = analysis_method.(method_name).eloreta_th;
                    [subject,properties]                                     = activation_level_eloreta(subject,properties);
                case 'lcmv'
                    if(properties.general_params.run_by_trial.value)
                        trial_name                                                  = properties.trial_name;
                        subject.pathname                                         = fullfile(subject.subject_path,trial_name,text_level,'LCMV',band.name);
                    else
                        subject.pathname                                         = fullfile(subject.subject_path,text_level,'LCMV',band.name);
                    end
                    if(~isfolder(subject.pathname))
                        mkdir(subject.pathname);
                    end
                    properties.activation_params.gamma1                             = analysis_method.(method_name).gamma1;
                    properties.activation_params.gamma2                             = analysis_method.(method_name).gamma2;
                    properties.activation_params.delta_gamma                        = analysis_method.(method_name).delta_gamma;
                    properties.activation_params.lcmv_th                            = analysis_method.(method_name).lcmv_th;
                    [subject,properties]                                     = activation_level_lcmv(subject,properties);
            end

            disp(strcat("-->> End time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));

            reference_path = strsplit(subject.pathname,subject.name);
            if(~isfield(subject.BC_V_info,'activation_level'))
                iter = 1;
            else
                iter = length(subject.BC_V_info.activation_level) + 1;
            end
            subject.BC_V_info.activation_level(iter).Comment     = 'Activation_level';
            subject.BC_V_info.activation_level(iter).Band        = band.name;
            subject.BC_V_info.activation_level(iter).Method      = lower(method_name);
            subject.BC_V_info.activation_level(iter).Freq        = char(band.str_band);
            subject.BC_V_info.activation_level(iter).Ref_path    = strrep(reference_path{2},'\','/');
            subject.BC_V_info.activation_level(iter).Name        = subject.file_name;
        end
    end
end
end

