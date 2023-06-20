function [subject,properties] = connectivity_level_interface(subject,properties)

%% Getting connectivity priors
[subject,properties]                            = get_connectivity_priors(subject,properties);

%%
%% Getting system response
%%
if(properties.general_params.system_response.value)
    [syst_resp_out]                             = get_system_response(subject,properties);
end
%%
%% Band Analysis, connectivity level
%%
% Saving general variables for connectivity level
file_name                                       = strcat('W.mat');
W                                               = subject.W;
if(properties.general_params.run_by_trial.value)
    trial_name                                  = properties.trial_name;
    pathname                                    = fullfile(subject.subject_path,trial_name,'Generals','Structural','HIGGS');
    reference_path                              = strsplit(pathname,subject.name);
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    sensor_level                                = subject.BC_V_info.sensor_level(contains({subject.BC_V_info.sensor_level.Ref_path},properties.trial_name));
    activation_level                            = subject.BC_V_info.activation_level(contains({subject.BC_V_info.activation_level.Ref_path},properties.trial_name));
else
    pathname                                    = fullfile(subject.subject_path,'Generals','Structural','HIGGS');
    reference_path                              = strsplit(pathname,subject.name);
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    sensor_level                                = subject.BC_V_info.sensor_level;
    activation_level                            = subject.BC_V_info.activation_level;
end
if(~isfield(subject.BC_V_info,'generals'))
    iter = 1;
else
    iter = length(subject.BC_V_info.generals) + 1;
end
subject.BC_V_info.generals(iter).Comment     = 'Generals';
subject.BC_V_info.generals(iter).Ref_path    = strrep(reference_path{2},'\','/');
subject.BC_V_info.generals(iter).Name        = file_name;
disp(strcat("File: ", file_name));
parsave(fullfile(pathname ,file_name ),W);
pos = 1;
for i=1:length(activation_level)
    activ_file                                   = activation_level(i);
    sensor_file                                  = sensor_level(contains({sensor_level.Freq},activ_file.Freq));
    subject.sensor_level_out                     = load(fullfile(subject.subject_path,sensor_file.Ref_path,sensor_file.Name));
    subject.activation_level_out                 = load(fullfile(subject.subject_path,activ_file.Ref_path,activ_file.Name));
    if(properties.general_params.system_response.value)
        if(properties.general_params.run_by_trial.value)
            subject.activation_level_out.indms   = syst_resp_out.(trial_name).(activ_file.Method).indms;
            subject.activation_level_out.stat    = syst_resp_out.(trial_name).(activ_file.Method).ave_stat;
            subject.activation_level_out.method  = activ_file.Method;
        else
            subject.activation_level_out.indms   = syst_resp_out.(activ_file.Method).indms;
            subject.activation_level_out.stat    = syst_resp_out.(activ_file.Method).ave_stat;
            subject.activation_level_out.method  = activ_file.Method;
        end
    end
    band  = subject.sensor_level_out.band;
    
    %%
    %% Defining path
    %%
    disp('=================================================================');
    disp(strcat( "BC-V-->> Connectivity level for frequency band: ",band.str_band ));
    text_level = 'Connectivity_level';

    %%
    %% Band Analysis, connectivity level
    %%
    for m=1:length(properties.connectivity_params.methods)
        analysis_method                                             = properties.connectivity_params.methods{m};        
        method                                                      = analysis_method.method;  
        if(analysis_method.run)
            disp('-----------------------------------------------------------------');
            disp(strcat("-->> Start time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
            switch m
                case 1
                    if(properties.general_params.run_by_trial.value)
                        trial_name                                  = properties.trial_name;
                        subject.pathname                         = fullfile(subject.subject_path,trial_name,text_level,'HIGGS',band.name);
                    else
                        subject.pathname                         = fullfile(subject.subject_path,text_level,'HIGGS',band.name);
                    end
                    if(~isfolder(subject.pathname))
                        mkdir(subject.pathname);
                    end
                   
                    [Thetajj,s2j,Tjv,llh,properties]                = connectivity_level_higgs(subject,properties);
                case 2
                    if(properties.general_params.run_by_trial.value)
                        trial_name                                  = properties.trial_name;
                        subject.pathname                         = fullfile(subject.subject_path,trial_name,text_level,'HG_LASSO',band.name);
                    else
                        subject.pathname                         = fullfile(subject.subject_path,text_level,'HG_LASSO',band.name);
                    end
                    if(~isfolder(subject.pathname))
                        mkdir(subject.pathname);
                    end
                    properties.connectivity_params.hg_lasso_th      = analysis_method.threshold;
                    [Thetajj,Sjj,Sigmajj]                           = connectivity_level_hg_lasso(subject,properties);
            end
            disp(strcat("-->> End time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
            subject                                             = BC_V_save(properties,subject,'activation',method,outputs,pos,band);
            pos                                                 = pos + 1;

            reference_path = strsplit(subject.pathname,subject.name);
            if(~isfield(subject.BC_V_info,'connectivity_level'))
                iter = 1;
            else
                iter = length(subject.BC_V_info.connectivity_level) + 1;
            end
            subject.BC_V_info.connectivity_level(iter).Comment   = 'Connectivity_level';
            subject.BC_V_info.connectivity_level(iter).Band      = band.name;
            subject.BC_V_info.connectivity_level(iter).Method    = lower(method_name);
            subject.BC_V_info.connectivity_level(iter).Freq      = band.str_band;
            subject.BC_V_info.connectivity_level(iter).Ref_path  = strrep(reference_path{2},'\','/');
            subject.BC_V_info.connectivity_level(iter).Name      = subject.file_name;
        end
    end 
end
end

