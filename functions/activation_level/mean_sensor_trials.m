function [subject,properties] = mean_sensor_trials(subject,properties)


sensor_trials                = subject.BC_V_info.trials;
sensor_tensor                = {};
for t=1:length(sensor_trials)
    bands = sensor_trials(t).sensor_level;
    for f=1:length(bands)
        sensor = bands(f);
        ref_path                    = sensor.Ref_path;
        file_name                   = sensor.Name;
        sensor_level    = load(fullfile(subject.subject_path,ref_path,file_name));
        sensor_tensor{t,f} = sensor_level.Svv;
    end
end
no_trials = size(sensor_tensor,1);
no_freqs = size(sensor_tensor,2);
no_channels = size(sensor_tensor{1,1},1);

sensor_tensor_out = zeros(no_trials,no_freqs,no_channels,no_channels);
for t=1:no_trials
    for f=1:no_freqs
        sensor_tensor_out(t,f,:,:) = sensor_tensor{t,f};
    end
end
sensor_out = mean(sensor_tensor_out,1);
subject    = BC_V_save(properties,subject,'sensor_mean',Svv,peak_pos,Nseg,band,pos, trial_info);

end
