function [subject,properties] = connectivity_level_interface(subject,properties)

%% Getting connectivity priors
[subject,properties]                                = get_connectivity_priors(subject,properties);
subject                                             = BC_V_save(properties,subject,'c_priors');
%%
%% Getting system response
%%
if(properties.general_params.system_response.value)
    [syst_resp_out]                                 = get_system_response(subject,properties);
end
%%
%% Band Analysis, connectivity level
%% 
if(properties.general_params.run_by_trial.value)    
    sensor_level                                    = subject.BC_V_info.sensor_level(contains({subject.BC_V_info.sensor_level.Ref_path},properties.trial_name));
    activation_level                                = subject.BC_V_info.activation_level(contains({subject.BC_V_info.activation_level.Ref_path},properties.trial_name));
else   
    sensor_level                                    = subject.BC_V_info.sensor_level;
    activation_level                                = subject.BC_V_info.activation_level;
end
pos = 1;
for i=1:length(activation_level)
    activ_file                                      = activation_level(i);
    sensor_file                                     = sensor_level(contains({sensor_level.Freq},activ_file.Freq));
    subject.sensor_level_out                        = load(fullfile(subject.subject_path,sensor_file.Ref_path,sensor_file.Name));
    subject.activation_level_out                    = load(fullfile(subject.subject_path,activ_file.Ref_path,activ_file.Name));
    if(properties.general_params.system_response.value)
        if(properties.general_params.run_by_trial.value)
            subject.activation_level_out.indms      = syst_resp_out.(trial_name).(activ_file.Method).indms;
            subject.activation_level_out.stat       = syst_resp_out.(trial_name).(activ_file.Method).ave_stat;
            subject.activation_level_out.method     = activ_file.Method;
        else
            subject.activation_level_out.indms      = syst_resp_out.(activ_file.Method).indms;
            subject.activation_level_out.stat       = syst_resp_out.(activ_file.Method).ave_stat;
            subject.activation_level_out.method     = activ_file.Method;
        end
    end
    band  = subject.sensor_level_out.band;  
    disp("=====================================================================");
    disp(strcat( "BC-V-->> Connectivity level for frequency band: ",band.str_band ));

    %%
    %% Band Analysis, connectivity level
    %%
    for m=1:length(properties.connectivity_params.methods)
        analysis_method                             = properties.connectivity_params.methods{m};        
        method                                      = analysis_method.method;  
        if(analysis_method.run)
            disp('-----------------------------------------------------------------');
            disp(strcat("-->> Start time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
            switch m
                case 1                    
                    [subject,properties,outputs]    = connectivity_level_higgs(subject,properties);
                case 2                    
                    [subject,properties,outputs]    = connectivity_level_hg_lasso(subject,properties);
            end
            subject                                 = BC_V_save(properties,subject,'connectivity',method,outputs,pos,band);
            pos                                     = pos + 1;
            disp(strcat("-->> End time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
        end
    end 
end
end

