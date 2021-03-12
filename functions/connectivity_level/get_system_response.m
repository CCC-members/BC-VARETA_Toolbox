function [syst_resp_out] = get_system_response(subject,properties)

IsCurv        = properties.activation_params.IsCurv.value;
BC_V_info   = properties.BC_V_info;
act_methods = properties.activation_params.methods;
if(properties.general_params.run_by_trial.value) 
    trial_name  = properties.trial_name;
    variant     = 'trial';
    if(properties.general_params.run_frequency_bin.value)
        variant = 'trial_bin';
    end
else
    variant     = 'not_trial';
    if(properties.general_params.run_frequency_bin.value)
        variant = 'not_trial_bin';
    end
end

for i=1:length(act_methods)
    act_method  = act_methods{i};
    method_name = fieldnames(act_method);
    method_name = method_name{1};
    if(act_methods{i}.(method_name).run)
        switch variant
            case 'not_trial'
                bands_runned = fieldnames(BC_V_info.activation_level.(method_name));
                count_bin = 0;
                for j=1:length(bands_runned)
                    band_name = bands_runned{j};
                    file_name = BC_V_info.activation_level.(method_name).(band_name).name;
                    ref_path = BC_V_info.activation_level.(method_name).(band_name).ref_path;
                    result_file = fullfile(subject.subject_path,ref_path,file_name);
                    act_result = load(result_file);
                    if(isequal(j,1))
                        ave_stat = act_result.stat;
                    else
                        ave_stat = ave_stat + act_result.stat;
                    end
                    count_bin = count_bin + 1;
                end
            case 'trial'
                bands_runned = fieldnames(BC_V_info.(trial_name).activation_level.(method_name));
                count_bin = 0;
                for j=1:length(bands_runned)
                    band_name = bands_runned{j};
                    file_name = BC_V_info.(trial_name).activation_level.(method_name).(band_name).name;
                    ref_path = BC_V_info.(trial_name).activation_level.(method_name).(band_name).ref_path;
                    result_file = fullfile(subject.subject_path,ref_path,file_name);
                    act_result = load(result_file);
                    if(isequal(j,1))
                        ave_stat = act_result.stat;
                    else
                        ave_stat = ave_stat + act_result.stat;
                    end
                    count_bin = count_bin + 1;
                end
            case 'not_trial_bin'
                bands_runned = fieldnames(BC_V_info.activation_level.(method_name));
                count_bin = 0;
                for j=1:length(bands_runned)
                    band_name = bands_runned{j};
                    bin_files = fieldnames(BC_V_info.activation_level.(method_name).(band_name));
                    for k=1:length(bin_files)
                        bin = bin_files{k};
                        file_name = BC_V_info.activation_level.(method_name).(band_name).(bin).name;
                        ref_path = BC_V_info.activation_level.(method_name).(band_name).(bin).ref_path;
                        result_file = fullfile(subject.subject_path,ref_path,file_name);
                        act_result = load(result_file);
                        if(isequal(j,1) && isequal(k,1))
                            ave_stat = act_result.stat;
                        else
                            ave_stat = ave_stat + act_result.stat;
                        end
                        count_bin = count_bin + 1;
                    end
                end
            case 'trial_bin'
                bands_runned = fieldnames(BC_V_info.(trial_name).activation_level.(method_name));
                count_bin = 0;
                for j=1:length(bands_runned)
                    band_name = bands_runned{j};
                    bin_files = fieldnames(BC_V_info.(trial_name).activation_level.(method_name).(band_name));
                    for k=1:length(bin_files)
                        bin = bin_files{k};
                        file_name = BC_V_info.(trial_name).activation_level.(method_name).(band_name).(bin).name;
                        ref_path = BC_V_info.(trial_name).activation_level.(method_name).(band_name).(bin).ref_path;
                        result_file = fullfile(subject.subject_path,ref_path,file_name);
                        act_result = load(result_file);
                        if(isequal(j,1) && isequal(k,1))
                            ave_stat = act_result.stat;
                        else
                            ave_stat = ave_stat + act_result.stat;
                        end
                        count_bin = count_bin + 1;
                    end
                end
        end
        ave_stat = ave_stat/count_bin;
        syst_resp_out.(method_name).ave_stat = ave_stat;
        param_name = strcat(method_name,'_th');
        if IsCurv == 0
            indms                                    = find(ave_stat > properties.activation_params.methods{i}.(method_name).(param_name).value);
            syst_resp_out.(method_name).indms        = indms;
        else
            indms_giri                               = find(ave_stat(:,1) > properties.activation_params.methods{i}.(method_name).(param_name).value);
            indms_sulc                               = find(ave_stat(:,2) > properties.activation_params.methods{i}.(method_name).(param_name).value);
            indms                                    = unique([indms_giri;indms_sulc]);
            syst_resp_out.(method_name).indms_giri   = indms_giri;
            syst_resp_out.(method_name).indms_sulc   = indms_sulc;
            syst_resp_out.(method_name).indms        = indms;
        end
        disp(strcat("BC-V-->> The average number of responsive sources across frequencies detected by ",method_name," screening was: " ,num2str(length(indms))));
    end
end
end

