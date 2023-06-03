function [subject,properties] = sensor_level_analysis(subject,properties)

text_level                  = 'Sensor_level';

%%
%% Estimating cross-spectra
%%
disp('BC-V-->> Estimating cross-spectra for M/EEG data.');
[Svv_channel,~,PSD,Nseg]    = cross_spectra(subject, properties);

%%
%% Saving general variables for sensor level
%%
file_name      = strcat('Data_spectrum.mat');
disp(strcat("File: ", file_name));
if(properties.general_params.run_by_trial.value)
    trial_name  = properties.trial_name;
    pathname              = fullfile(subject.subject_path,trial_name,text_level);
    pathname_generals     = fullfile(subject.subject_path,trial_name,'Generals');
    pathname_funtional    = fullfile(subject.subject_path,trial_name,'Generals','Functional');
else
    pathname              = fullfile(subject.subject_path,text_level);
    pathname_generals     = fullfile(subject.subject_path,'Generals');
    pathname_funtional    = fullfile(subject.subject_path,'Generals','Functional');
end
if(~isfolder(pathname))
    mkdir(pathname);
    mkdir(pathname_generals);
    mkdir(pathname_funtional);
end
if(~properties.general_params.run_frequency_bin.value &&  properties.general_params.run_frequency_bin.band_mean)
    parsave(fullfile(pathname_funtional ,file_name ),Svv_channel);
else
    parsave(fullfile(pathname_funtional ,file_name ),Svv_channel,PSD);
end
reference_path                          = strsplit(pathname_funtional,subject.name);
subject.BC_V_info.generals(1).Comment   = 'Generals';
subject.BC_V_info.generals(1).Ref_path  = strrep(reference_path{2},'\','/');
subject.BC_V_info.generals(1).Name      = file_name;

%% Sensor analysis
for h=1:length(properties.sensor_params.frequencies)
    band                                            = properties.sensor_params.frequencies(h);
    if(band.run)
        % Get band
        if(~isempty(PSD))
            [Svv,peak_pos,subject,properties]       = get_band(Svv_channel,PSD,band,subject,properties);
        else
            Svv = Svv_channel(:,:,h);
            peak_pos = h;
        end
        disp('=================================================================');
        disp(strcat( "BC-V-->> Sensor level for frequency band: " , band.str_band));
        
        if(properties.general_params.run_by_trial.value)
            trial_name  = properties.trial_name;
            pathname    = fullfile(subject.subject_path,trial_name,text_level,band.name);
        else
            pathname    = fullfile(subject.subject_path,text_level,band.name);
        end
        if(~isfolder(pathname))
            mkdir(pathname);
        end

        %% Saving files
        disp('-->> Saving file')
        file_name                               = strcat('Sensor_level_',band.str_band,'.mat');
        disp(strcat("File: ", file_name));
        parsave(fullfile(pathname ,file_name ),Svv,peak_pos,Nseg,band);
        reference_path                          = strsplit(pathname,subject.name);
        subject.BC_V_info.sensor_level(h).Comment       = 'Sensor_level';
        [~,band_name,~]                         = fileparts(reference_path{2});
        subject.BC_V_info.sensor_level(h).Band          = band_name;
        subject.BC_V_info.sensor_level(h).Freq          = char(band.str_band);
        subject.BC_V_info.sensor_level(h).Ref_path      = strrep(reference_path{2},'\','/');
        subject.BC_V_info.sensor_level(h).Name          = file_name;
    end
end
end

