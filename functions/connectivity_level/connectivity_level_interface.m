function [subject,properties] = connectivity_level_interface(subject,properties)

band  = properties.sensor_level_out.band;
%%
%% Defining path
%%
disp('=================================================================');
if(isfield(band,'f_bin'))    
    disp(strcat( 'BC-V-->> Connectivity level for frequency band: (' , band.name , ') bin ->>>' , string(band.f_bin), 'Hz') );
    properties.str_band =  strcat( band.name,'_',string(band.f_bin),'Hz');
else
    disp(strcat( 'BC-V-->> Connectivity level for frequency band: (' , band.name , ') ' , string(band.f_start), 'Hz-->' , string(band.f_end) , 'Hz') );
    properties.str_band =  strcat( band.name,'_',string(band.f_start),'Hz_',string(band.f_end),'Hz');
end
text_level = 'Connectivity_level'; 

%%
%% Band Analysis, connectivity level
%%
for m=1:length(properties.connectivity_params.methods)
    analysis_method = properties.connectivity_params.methods{m};
    fields          = fieldnames(analysis_method);
    method_name     = fields{1};
    if(analysis_method.(method_name).run)
        disp('-----------------------------------------------------------------');
        disp(strcat("-->> Start time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
        switch method_name
            case 'higgs'
                if(properties.BC_V_info.properties.general_params.run_by_trial.value)
                    trial_name          = properties.trial_name;
                    properties.pathname = fullfile(subject.subject_path,trial_name,text_level,'HIGGS',band.name);
                else
                    properties.pathname = fullfile(subject.subject_path,text_level,'HIGGS',band.name);
                end
                if(~isfolder(properties.pathname))
                    mkdir(properties.pathname);
                end
                properties.connectivity_params.higgs_th = analysis_method.(method_name).higgs_th;
                [Thetajj,s2j,Tjv,llh]                   = connectivity_level_higgs(subject,properties);
            case 'hg_lasso'
                if(properties.BC_V_info.properties.general_params.run_by_trial.value)
                    trial_name          = properties.trial_name;
                    properties.pathname = fullfile(subject.subject_path,trial_name,text_level,'HG_LASSO',band.name);
                else
                    properties.pathname = fullfile(subject.subject_path,text_level,'HG_LASSO',band.name);
                end
                if(~isfolder(properties.pathname))
                    mkdir(properties.pathname);
                end
                properties.connectivity_params.hg_lasso_th  = analysis_method.(method_name).hg_lasso_th;
                [Thetajj,Sjj,Sigmajj]                       = connectivity_level_hg_lasso(subject,properties);
        end
        disp(strcat("-->> End time: ",datestr(now,'mmmm dd, yyyy HH:MM:SS AM')));
        
          reference_path = strsplit(properties.pathname,subject.name);
        if(properties.BC_V_info.properties.general_params.run_by_trial.value)
            if(properties.BC_V_info.properties.general_params.run_frequency_bin.value)
                f_bin                                                                                               = replace(num2str(band.f_bin),'.','_');
                f_bin                                                                                               = strcat(band.name,'_',f_bin);
                properties.BC_V_info.(trial_name).connectivity_level.(method_name).(band.name).(f_bin).name         = properties.file_name;
                properties.BC_V_info.(trial_name).connectivity_level.(method_name).(band.name).(f_bin).ref_path     = reference_path{2};
            else
                properties.BC_V_info.(trial_name).connectivity_level.(method_name).(band.name).name                 = properties.file_name;
                properties.BC_V_info.(trial_name).connectivity_level.(method_name).(band.name).ref_path             = reference_path{2};
            end
        else
            if(properties.BC_V_info.properties.general_params.run_frequency_bin.value)
                f_bin                                                                               = replace(num2str(band.f_bin),'.','_');
                f_bin                                                                               = strcat(band.name,'_',f_bin);
                properties.BC_V_info.connectivity_level.(method_name).(band.name).(f_bin).name      = properties.file_name;
                properties.BC_V_info.connectivity_level.(method_name).(band.name).(f_bin).ref_path  = reference_path{2};
            else
                properties.BC_V_info.connectivity_level.(method_name).(band.name).name              = properties.file_name;
                properties.BC_V_info.connectivity_level.(method_name).(band.name).ref_path          = reference_path{2};
            end
        end
    end
end

end

