function [subject,properties] = data_analysis(subject,properties)
% try
%%
%% Running sensor level
%%
if(isequal(properties.analysis_level,1))    
    %%
    %% Preparing params
    %%
    data    = subject.data; 
    Fs      = properties.spectral_params.samp_freq.value;       % sampling frequency
    Fm      = properties.spectral_params.max_freq.value;        % maximum frequency
    deltaf  = properties.spectral_params.freq_resol.value;      % frequency resolution
    varf    = properties.spectral_params.freq_gfiltvar.value;   % gaussian filter variance
    Nw      = properties.spectral_params.win_order.value;       % Slepian windows
    
    %%
    %% Estimating cross-spectra
    %%
    disp('BC-V-->> Estimating cross-spectra for M/EEG data.');
    [Svv_channel,~,PSD,Nf,F,Nseg]           = cross_spectra(data,Fs,Fm,deltaf,subject.Ke,varf,Nw,'app_properties',properties);
    
    %% Adding fieltrip external functions
    addpath(genpath('external/fieldtrip'));
    ft_defaults
       
    %% Saving general variables for sensor level
    disp('-->> Saving file')
    file_name                               = strcat('Data_spectrum.mat');
    disp(strcat("File: ", file_name));
    pathname                                = fullfile(subject.subject_path,'Sensor_level');
    pathname_generals                       = fullfile(subject.subject_path,'Generals');
    pathname_funtional                      = fullfile(subject.subject_path,'Generals','Funtional');
    if(~isfolder(pathname))
        mkdir(pathname);
        mkdir(pathname_generals);
        mkdir(pathname_funtional);
    end
    parsave(fullfile(pathname_funtional ,file_name ),Svv_channel);
    reference_path                          = strsplit(pathname_funtional,subject.name);
    if(properties.general_params.run_by_trial.value) 
        trial_name                          = properties.trial_name;
        properties.BC_V_info.(trial_name).generals.funtional.spectrum_data.name         = file_name;
        properties.BC_V_info.(trial_name).generals.funtional.spectrum_data.ref_path     = reference_path{2};        
    else        
        properties.BC_V_info.generals.funtional.spectrum_data.name                      = file_name;
        properties.BC_V_info.generals.funtional.spectrum_data.ref_path                  = reference_path{2};        
    end
    properties.Nseg = Nseg;
    for h=1:length(properties.spectral_params.frequencies)
        band                                = properties.spectral_params.frequencies(h);
        if(band.run)            
            properties.band                 = band;
            [subject,properties]            = sensor_level_analysis(Svv_channel,PSD,Nf,F,subject,properties);
        end
    end
end
%%
%% Running activation level
%%
if(isequal(properties.analysis_level,2))
    % Saving general variables for activation level
    disp('-->> Saving file')
    file_name                                                               = strcat('W.mat');
    disp(strcat("File: ", file_name));
    pathname                                                                = fullfile(subject.subject_path,'Generals','Structural');
    if(~isfolder(pathname))
        mkdir(pathname);
    end
    W                                                                       = subject.W;
    parsave(fullfile(pathname ,file_name ),W);
    reference_path                                                          = strsplit(pathname,subject.name);
    if(properties.general_params.run_by_trial.value) 
        trial_name = properties.trial_name;
        properties.BC_V_info.(trial_name).generals.structural.W.name        = file_name;
        properties.BC_V_info.(trial_name).generals.structural.W.ref_path    = reference_path{2};        
    else        
        properties.BC_V_info.generals.structural.W.name                     = file_name;
        properties.BC_V_info.generals.structural.W.ref_path                 = reference_path{2};        
    end
    if(properties.general_params.run_by_trial.value)
        sensor_level                                                        = properties.BC_V_info.(properties.trial_name).sensor_level;
    else
        sensor_level                                                        = properties.BC_V_info.sensor_level;
    end    
    band_fields = fieldnames(sensor_level);
    for h=1:length(band_fields)
        band_name = band_fields{h};
        if(properties.general_params.run_frequency_bin.value)
            bin_fields                                                      = fieldnames(sensor_level.(band_name));
            for k=1:length(bin_fields)
                bin_name                                                    = bin_fields{k};
                properties.sensor_level_out                                 = load(fullfile(subject.subject_path,sensor_level.(band_name).(bin_name).ref_path,sensor_level.(band_name).(bin_name).name));
                [subject,properties]                                        = activation_level_interface(subject,properties);
            end
        else
            properties.sensor_level_out                                     = load(fullfile(subject.subject_path,sensor_level.(band_name).ref_path,sensor_level.(band_name).name));
            [subject,properties]                                            = activation_level_interface(subject,properties);
        end
    end
end
%%
%% Running connectivity level
%%
if(isequal(properties.analysis_level,3))    
    %%
    %% Getting system response
    %%
    if(properties.general_params.system_response.value)
        [syst_resp_out]                                         = get_system_response(subject,properties);
    end
    %%
    %% Band Analysis, connectivity level
    %%    
    if(properties.general_params.run_by_trial.value)
        sensor_level                                            = properties.BC_V_info.(properties.trial_name).sensor_level;
        activation_level                                        = properties.BC_V_info.(properties.trial_name).activation_level;
    else
        sensor_level                                            = properties.BC_V_info.sensor_level;
        activation_level                                        = properties.BC_V_info.activation_level;
    end
    activ_fields = fieldnames(activation_level);
    for i=1:length(activ_fields)
        act_method                                              = activ_fields{i};        
        band_fields                                             = fieldnames(sensor_level);
        for h=1:length(band_fields)
            band_name                                           = band_fields{h};
            if(properties.general_params.run_frequency_bin.value)
                bin_fields                                      = fieldnames(sensor_level.(band_name));
                for k=1:length(bin_fields)
                    bin_name                                    = bin_fields{k};
                    properties.sensor_level_out                 = load(fullfile(subject.subject_path,sensor_level.(band_name).(bin_name).ref_path,sensor_level.(band_name).(bin_name).name));
                    properties.activation_level_out             = load(fullfile(subject.subject_path,activation_level.(act_method).(band_name).(bin_name).ref_path,activation_level.(act_method).(band_name).(bin_name).name));
                    if(properties.general_params.system_response.value)
                        properties.activation_level_out.indms   = syst_resp_out.(act_method).indms;
                        properties.activation_level_out.stat    = syst_resp_out.(act_method).ave_stat;
                    end
                    [subject,properties]                        = connectivity_level_interface(subject,properties);
                end
            else
                properties.sensor_level_out                     = load(fullfile(subject.subject_path,sensor_level.(band_name).ref_path,sensor_level.(band_name).name));
                properties.activation_level_out                 = load(fullfile(subject.subject_path,activation_level.(act_method).(band_name).ref_path,activation_level.(act_method).(band_name).name));
                if(properties.general_params.system_response.value)
                    properties.activation_level_out.indms       = syst_resp_out.(act_method).indms;
                    properties.activation_level_out.stat        = syst_resp_out.(act_method).ave_stat;
                end
                [subject,properties]                            = connectivity_level_interface(subject,properties);
            end
        end
    end
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

