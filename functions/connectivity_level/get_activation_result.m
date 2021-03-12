function [subject,properties] = get_activation_result(subject,properties)


BC_V_info  = properties.BC_V_info;
band       = properties.band;
act_method = properties.act_method;

if(properties.general_params.run_by_trial.value) 
    trial_name = properties.trial_name;
    variant = 'trial';
    if(properties.general_params.run_frequency_bin.value)
        variant = 'trial_bin';
        bin     = band.f_bin;
    end
else
    variant = 'not_trial';
    if(properties.general_params.run_frequency_bin.value)
        variant = 'not_trial_bin';
        bin     = band.f_bin;
    end
end
switch variant
    case 'not_trial'
        ref_path    = BC_V_info.sensor_level.(band.name).ref_path;
        file_name   = BC_V_info.sensor_level.(band.name).name;
        s_ref_file  = fullfile(ref_path,file_name);
        ref_path    = BC_V_info.activation_level.(act_method).(band.name).ref_path;
        file_name   = BC_V_info.activation_level.(act_method).(band.name).name;
        a_ref_file  = fullfile(ref_path,file_name);
    case 'trial'
        ref_path    = BC_V_info.(trial_name).sensor_level.(band.name).ref_path;
        file_name   = BC_V_info.(trial_name).sensor_level.(band.name).name;
        s_ref_file  = fullfile(ref_path,file_name);
        ref_path    = BC_V_info.(trial_name).activation_level.(act_method).(band.name).ref_path;
        file_name   = BC_V_info.(trial_name).activation_level.(act_method).(band.name).name;
        a_ref_file  = fullfile(ref_path,file_name);
    case 'not_trial_bin'
        ref_path    = BC_V_info.sensor_level.(band.name).(bin).ref_path;
        file_name   = BC_V_info.sensor_level.(band.name).(bin).name;
        s_ref_file  = fullfile(ref_path,file_name);
        ref_path    = BC_V_info.activation_level.(act_method).(band.name).(bin).ref_path;
        file_name   = BC_V_info.activation_level.(act_method).(band.name).(bin).name;
        a_ref_file  = fullfile(ref_path,file_name);
    case 'trial_bin'
        ref_path    = BC_V_info.(trial_name).sensor_level.(band.name).(bin).ref_path;
        file_name   = BC_V_info.(trial_name).sensor_level.(band.name).(bin).name;
        s_ref_file  = fullfile(ref_path,file_name);
        ref_path    = BC_V_info.(trial_name).activation_level.(act_method).(band.name).(bin).ref_path;
        file_name   = BC_V_info.(trial_name).activation_level.(act_method).(band.name).(bin).name;
        a_ref_file  = fullfile(ref_path,file_name);
end

if(isequal(properties.general_params.BCV_work_dir,'local'))
    [subject_path,~] = fileparts(subject.subject_path);
else
    subject_path = subject.subject_path;
end
sensor_out = load(fullfile(subject_path,s_ref_file));
activation_out = load(fullfile(subject_path,a_ref_file));

properties.peak_pos     = sensor_out.peak_pos;
properties.str_band     = sensor_out.str_band;
properties.Nseg         = sensor_out.Nseg;
properties.Svv          = sensor_out.Svv;
properties.indms        = activation_out.indms;


end

