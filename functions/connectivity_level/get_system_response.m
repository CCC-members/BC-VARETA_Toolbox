function [syst_resp_out] = get_system_response(subject,properties)

IsCurv        = properties.activation_params.IsCurv.value;
BC_V_info   = subject.BC_V_info;
act_methods = properties.activation_params.methods;
if(properties.general_params.run_by_trial.value) 
    trial_name  = properties.trial_name;
    variant     = 'trial';   
else
    variant     = 'not_trial';    
end
for i=1:length(act_methods)
    act_method  = act_methods{i};
    method_name = fieldnames(act_method);
    method_name = method_name{1};
    if(act_methods{i}.(method_name).run)
        switch variant
            case 'not_trial'
                activ_files         = BC_V_info.activation_level(contains(lower({BC_V_info.activation_level.Ref_path}),method_name));
                count_bin           = 0;
                for j=1:length(activ_files)
                    activ_file      = activ_files(j);
                    ref_path        = activ_file.Ref_path;
                    file_name       = activ_file.Name;
                    result_file     = fullfile(subject.subject_path,ref_path,file_name);
                    act_result      = load(result_file);
                    if(isequal(j,1))
                        ave_stat    = act_result.stat;
                    else
                        ave_stat    = ave_stat + act_result.stat;
                    end
                    count_bin       = count_bin + 1;
                end
            case 'trial'
                trial_files         = BC_V_info.activation_level(contains({BC_V_info.activation_level.Ref_path},trial_name));
                activ_files         = BC_V_info.activation_level(contains(lower({trial_files.Ref_path}),method_name));
                count_bin           = 0;
                for j=1:length(activ_files)
                    activ_file      = activ_files(j);
                    ref_path        = activ_file.Ref_path;
                    file_name       = activ_file.Name;
                    result_file     = fullfile(subject.subject_path,ref_path,file_name);
                    act_result      = load(result_file);
                    if(isequal(j,1))
                        ave_stat    = act_result.stat;
                    else
                        ave_stat    = ave_stat + act_result.stat;
                    end
                    count_bin       = count_bin + 1;
                end            
        end        
        ave_stat                                = ave_stat/count_bin;
        if(properties.general_params.run_by_trial.value)
            syst_resp_out.(trial_name).(method_name).ave_stat    = ave_stat;
        else
            syst_resp_out.(method_name).ave_stat    = ave_stat;
        end
        param_name                              = strcat(method_name,'_th');
        if IsCurv == 0
            indms                                                       = find(ave_stat > properties.activation_params.methods{i}.(method_name).(param_name).value);
            if(properties.general_params.run_by_trial.value)
                syst_resp_out.(trial_name).(method_name).indms          = indms;
            else
                syst_resp_out.(method_name).indms                       = indms;
            end
        else
            indms_giri                                                  = find(ave_stat(:,1) > properties.activation_params.methods{i}.(method_name).(param_name).value);
            indms_sulc                                                  = find(ave_stat(:,2) > properties.activation_params.methods{i}.(method_name).(param_name).value);
            indms                                                       = unique([indms_giri;indms_sulc]);
            if(properties.general_params.run_by_trial.value)
                syst_resp_out.(trial_name).(method_name).indms_giri     = indms_giri;
                syst_resp_out.(trial_name).(method_name).indms_sulc     = indms_sulc;
                syst_resp_out.(trial_name).(method_name).indms          = indms;
            else
                syst_resp_out.(method_name).indms_giri                  = indms_giri;
                syst_resp_out.(method_name).indms_sulc                  = indms_sulc;
                syst_resp_out.(method_name).indms                       = indms;
            end            
        end
        disp(strcat("BC-V-->> The average number of responsive sources across frequencies detected by ",method_name," screening was: " ,num2str(length(indms))));
    end
end
end

