function [subject,properties] = sensor_level_analysis(subject,properties)

%%
%% Preparing params
%%
Lvj             = subject.Headmodel.Gain;
Fs              = properties.sensor_params.samp_freq.value;         % sampling frequency
Fmax            = properties.sensor_params.max_freq.value;          % maximum frequency
deltaf          = properties.sensor_params.freq_resol.value;        % frequency resolution
F               = 0:deltaf:Fmax;                                    % frequency vector
if(isfield(subject.MEEG,'data'))
    data        = subject.MEEG.data;    
    varf        = properties.sensor_params.freq_gfiltvar.value;     % gaussian filter variance
    Nw          = properties.sensor_params.win_order.value;         % Slepian windows
    method      = properties.sensor_params.method.value;            % Spectrum method
    %%
    %% Estimates the Cross Spectrum of the input M/EEG data
    %%
    disp("=====================================================================");
    disp('BC-V-->> Estimating cross-spectra for M/EEG data.');
    switch lower(method)
        case "hilbert"
            if(~properties.general_params.run_frequency_bin.value &&  properties.general_params.run_frequency_bin.band_mean)
                [Svv_channel,Nf,Nseg,PSD]   = xspectrum_band(data, Fs, deltaf, varf, Nw, properties);
            else
                [Svv_channel,Nf,Nseg,PSD]   = xspectrum(data, Fs, Fmax, deltaf, varf, Nw, properties);
            end
        case "thomson"
            [Svv_channel,Nf,Nseg,PSD]       = xspectrum_thomson(data, Fs, Fmax, deltaf, properties);
    end

    %% Applying average reference
    disp('-->> Applying average reference.');
    for jj = 1:Nf
        [Svv_channel(:,:,jj),Lvj]           = applying_reference(Svv_channel(:,:,jj),Lvj);    % applying average reference...
    end
else
    Svv_channel                             = subject.MEEG.CrossM;
    PSD                                     = subject.MEEG.Spec;
    Nseg                                    = [];
end
subject                                     = BC_V_save(properties,subject,'fuctional',Svv_channel,PSD);

%% Sensor analysis
for pos=1:length(properties.sensor_params.frequencies)
    band                                = properties.sensor_params.frequencies(pos);
    disp("=====================================================================");
    disp(strcat( "BC-V-->> Sensor level for frequency band: " , band.str_band));
    if(~isempty(PSD) && isfield(subject.MEEG,'data'))
        if(isfield(band,'f_bin'))
            [f1,nf1]                = min(abs(F - band.f_bin));
            [f2,nf2]                = min(abs(F - band.f_bin));
        else
            [f1,nf1]                = min(abs(F - band.f_start));
            [f2,nf2]                = min(abs(F - band.f_end));
        end
        peak_pos                    = nf1:nf2;
        Svv                         = mean(Svv_channel(:,:,peak_pos),3);
    else
        Svv                         = Svv_channel(:,:,pos);
        peak_pos                    = pos;
    end
    subject                         = BC_V_save(properties,subject,'sensor',Svv,peak_pos,Nseg,band,pos);

end
end

