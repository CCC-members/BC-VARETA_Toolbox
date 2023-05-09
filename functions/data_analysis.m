function [subject,properties] = data_analysis(subject,properties)
% try
%%
%% Running sensor level
%%
if(isequal(properties.analysis_level,1))
   
    %%
    %% Estimating cross-spectra
    %%
    disp('BC-V-->> Estimating cross-spectra for M/EEG data.');
    [Svv_channel,~,PSD,Nseg] = cross_spectra(subject, properties);
    
    %% Adding fieltrip external functions
    f_path          = mfilename('fullpath');
    [ref_path,~,~]  = fileparts(fileparts(f_path));
    addpath(genpath(fullfile(ref_path,'external/fieldtrip')));
    ft_defaults
    
    %% Saving general variables for sensor level
    file_name      = strcat('Data_spectrum.mat');
    disp(strcat("File: ", file_name));
    if(properties.general_params.run_by_trial.value)
        trial_name  = properties.trial_name;
        pathname              = fullfile(subject.subject_path,trial_name,'Sensor_level');
        pathname_generals     = fullfile(subject.subject_path,trial_name,'Generals');
        pathname_funtional    = fullfile(subject.subject_path,trial_name,'Generals','Funtional');
    else
        pathname              = fullfile(subject.subject_path,'Sensor_level');
        pathname_generals     = fullfile(subject.subject_path,'Generals');
        pathname_funtional    = fullfile(subject.subject_path,'Generals','Funtional');
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
    if(~isfield(properties.BC_V_info,'generals'))
        iter = 1;
    else
        iter = length(properties.BC_V_info.generals) + 1;
    end
    properties.BC_V_info.generals(iter).Comment   = 'Generals';
    properties.BC_V_info.generals(iter).Ref_path  = strrep(reference_path{2},'\','/');
    properties.BC_V_info.generals(iter).Name      = file_name;
    
    %% Sensor analysis
    properties.Nseg = Nseg;
    for h=1:length(properties.sensor_params.frequencies)
        band                                = properties.sensor_params.frequencies(h);
        if(band.run)
            % Get band
            if(~isempty(PSD))
                [Svv,subject,properties]    = get_band(Svv_channel,PSD,band,subject,properties);
            else
                Svv = Svv_channel(:,:,h);
                properties.peak_pos = h;
            end            
            [subject,properties]            = sensor_level_analysis(Svv,band,subject,properties);
        end
    end
    return;
end
%%
%% Running activation level
%%
if(isequal(properties.analysis_level,2))
    % Saving general variables for activation level
    file_name                                       = strcat('W.mat');
    W                                               = subject.W;
    
    if(properties.general_params.run_by_trial.value)
        trial_name = properties.trial_name;
        pathname                                    = fullfile(subject.subject_path,trial_name,'Generals','Structural','sSSBL');
        reference_path                              = strsplit(pathname,subject.name);
        if(~isfolder(pathname))
            mkdir(pathname);
        end
        sensor_level                                = properties.BC_V_info.sensor_level(contains({properties.BC_V_info.sensor_level.Ref_path},properties.trial_name));
    else
        pathname                                    = fullfile(subject.subject_path,'Generals','Structural','sSSBL');
        reference_path                              = strsplit(pathname,subject.name);
        if(~isfolder(pathname))
            mkdir(pathname);
        end
        sensor_level                                = properties.BC_V_info.sensor_level;
    end
    if(~isfield(properties.BC_V_info,'generals'))
        iter                                        = 1;
    else
        iter                                        = length(properties.BC_V_info.generals) + 1;
    end
    properties.BC_V_info.generals(iter).Comment     = 'Generals';
    properties.BC_V_info.generals(iter).Ref_path    = strrep(reference_path{2},'\','/');
    properties.BC_V_info.generals(iter).Name        = file_name;    
    disp(strcat("File: ", file_name));
    parsave(fullfile(pathname ,file_name ),W);
    
    for h=1:length(sensor_level)
        ref_path                                    = sensor_level(h).Ref_path;
        file_name                                   = sensor_level(h).Name;
        properties.sensor_level_out                 = load(fullfile(subject.subject_path,ref_path,file_name));
        [subject,properties]                        = activation_level_interface(subject,properties);
    end
    return;
end
%%
%% Running connectivity level
%%
if(isequal(properties.analysis_level,3))
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
        sensor_level                                = properties.BC_V_info.sensor_level(contains({properties.BC_V_info.sensor_level.Ref_path},properties.trial_name));
        activation_level                            = properties.BC_V_info.activation_level(contains({properties.BC_V_info.activation_level.Ref_path},properties.trial_name));
    else
        pathname                                    = fullfile(subject.subject_path,'Generals','Structural','HIGGS');
        reference_path                              = strsplit(pathname,subject.name);
        if(~isfolder(pathname))
            mkdir(pathname);
        end
        sensor_level                                = properties.BC_V_info.sensor_level;
        activation_level                            = properties.BC_V_info.activation_level;
    end
    if(~isfield(properties.BC_V_info,'generals'))
        iter = 1;
    else
        iter = length(properties.BC_V_info.generals) + 1;
    end
    properties.BC_V_info.generals(iter).Comment     = 'Generals';
    properties.BC_V_info.generals(iter).Ref_path    = strrep(reference_path{2},'\','/');
    properties.BC_V_info.generals(iter).Name        = file_name;
    disp(strcat("File: ", file_name));
    parsave(fullfile(pathname ,file_name ),W);
    
    for i=1:length(activation_level)
        activ_file                                      = activation_level(i);
        sensor_file                                     = sensor_level(contains({sensor_level.Freq},activ_file.Freq));        
        properties.sensor_level_out                     = load(fullfile(subject.subject_path,sensor_file.Ref_path,sensor_file.Name));
        properties.activation_level_out                 = load(fullfile(subject.subject_path,activ_file.Ref_path,activ_file.Name));
        if(properties.general_params.system_response.value)
            if(properties.general_params.run_by_trial.value)
                properties.activation_level_out.indms   = syst_resp_out.(trial_name).(activ_file.Method).indms;
                properties.activation_level_out.stat    = syst_resp_out.(trial_name).(activ_file.Method).ave_stat;
                properties.activation_level_out.method  = activ_file.Method;
            else
                properties.activation_level_out.indms   = syst_resp_out.(activ_file.Method).indms;
                properties.activation_level_out.stat    = syst_resp_out.(activ_file.Method).ave_stat;
                properties.activation_level_out.method  = activ_file.Method;
            end
        end
        [subject,properties]                            = connectivity_level_interface(subject,properties);
        
    end
    return;
end
% catch ME
%     if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
%       msg = ['Dimension mismatch occurred: First argument has ', ...
%             num2str(size(A,2)),' columns while second has ', ...
%             num2str(size(B,2)),' columns.'];
%         causeException = MException('MATLAB:myCode:dimensions',msg);
%         ME = addCause(ME,causeException);
%    end
%    rethrow(ME)
% end

end

