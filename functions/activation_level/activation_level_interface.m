function [subject,properties] = activation_level_interface(subject,properties)

band  = properties.sensor_level_out.band;
%%
%% Defining path
%%
disp('=================================================================');
if(isfield(band,'f_bin'))    
    disp(strcat( 'BC-V-->> Activation level for frequency band: (' , band.name , ') bin ->>>' , string(band.f_bin), 'Hz') );    
    properties.str_band =  strcat( band.name,'_',string(band.f_bin),'Hz');
else
    disp(strcat( 'BC-V-->> Activation level for frequency band: (' , band.name , ') ' , string(band.f_start), 'Hz-->' , string(band.f_end) , 'Hz') );
    properties.str_band =  strcat( band.name,'_',string(band.f_start),'Hz_',string(band.f_end),'Hz');
end
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
                if(properties.BC_V_info.properties.general_params.run_by_trial.value)
                    trial_name                                                  = properties.trial_name;
                    properties.pathname                                         = fullfile(subject.subject_path,trial_name,text_level,'sSSBLpp',band.name);
                else
                    properties.pathname                                         = fullfile(subject.subject_path,text_level,'sSSBLpp',band.name);
                end
                if(~isfolder(properties.pathname))
                    mkdir(properties.pathname);
                end
                properties.activation_params.sssblpp_th                         = analysis_method.(method_name).sssblpp_th;
                [stat,J,T,indms,properties]                                     = activation_level_sssblpp(subject,properties);
            case 'eloreta'
                if(properties.general_params.run_by_trial.value)
                    trial_name                                                  = properties.trial_name;
                    properties.pathname                                         = fullfile(subject.subject_path,trial_name,text_level,'eLORETA',band.name);
                else
                    properties.pathname                                         = fullfile(subject.subject_path,text_level,'eLORETA',band.name);
                end
                if(~isfolder(properties.pathname))
                    mkdir(properties.pathname);
                end
                properties.activation_params.gamma1                             = analysis_method.(method_name).gamma1;
                properties.activation_params.gamma2                             = analysis_method.(method_name).gamma2;
                properties.activation_params.delta_gamma                        = analysis_method.(method_name).delta_gamma;
                properties.activation_params.eloreta_th                         = analysis_method.(method_name).eloreta_th;
                [stat,J,T,indms,properties]                                     = activation_level_eloreta(subject,properties);
            case 'lcmv'                
                if(properties.general_params.run_by_trial.value)
                    trial_name                                                  = properties.trial_name;
                    properties.pathname                                         = fullfile(subject.subject_path,trial_name,text_level,'LCMV',band.name);
                else
                    properties.pathname                                         = fullfile(subject.subject_path,text_level,'LCMV',band.name);
                end
                if(~isfolder(properties.pathname))
                    mkdir(properties.pathname);
                end
                properties.activation_params.gamma1                             = analysis_method.(method_name).gamma1;
                properties.activation_params.gamma2                             = analysis_method.(method_name).gamma2;
                properties.activation_params.delta_gamma                        = analysis_method.(method_name).delta_gamma;
                properties.activation_params.lcmv_th                            = analysis_method.(method_name).lcmv_th;
                [stat,J,T,indms,properties]                                     = activation_level_lcmv(subject,properties);
        end
        
        disp(strcat("-->> End time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
        
        reference_path = strsplit(properties.pathname,subject.name);
        if(~isfield(properties.BC_V_info,'activation_level'))
            iter = 1;
        else
            iter = length(properties.BC_V_info.activation_level) + 1;
        end
        properties.BC_V_info.activation_level(iter).Comment     = 'Activation_level';
        [base_path,band_name,~]                                 = fileparts(reference_path{2});
        properties.BC_V_info.activation_level(iter).Band        = band_name;
        [~,method,~]                                            = fileparts(base_path);
        properties.BC_V_info.activation_level(iter).Method      = lower(method);
        properties.BC_V_info.activation_level(iter).Freq        = char(properties.str_band);
        properties.BC_V_info.activation_level(iter).Ref_path    = reference_path{2};
        properties.BC_V_info.activation_level(iter).Name        = properties.file_name;
    end
end

end

