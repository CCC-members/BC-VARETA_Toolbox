function [Svv_channel,Lvj,PSD,Nseg] = cross_spectra(subject, properties)


% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019

% Updates
% - Ariosky Areces Gonzalez

% Date: January 30, 2021

%%
%% Preparing params
%%
data    = subject.MEEG.data;
Lvj     = subject.Headmodel.Gain;
Fs      = properties.sensor_params.samp_freq.value;       % sampling frequency
Fmax    = properties.sensor_params.max_freq.value;        % maximum frequency
deltaf  = properties.sensor_params.freq_resol.value;      % frequency resolution
varf    = properties.sensor_params.freq_gfiltvar.value;   % gaussian filter variance
Nw      = properties.sensor_params.win_order.value;       % Slepian windows

%% estimating cross-spectra...

% Remove fieldtrip path for override functions 
warning off;
rmpath(genpath(fullfile('external/fieldtrip')));
warning on;
% estimates the Cross Spectrum of the input M/EEG data
if(isequal(lower(properties.sensor_params.method.value),"hilbert"))
    if(~properties.general_params.run_frequency_bin.value &&  properties.general_params.run_frequency_bin.band_mean)
        [Svv_channel,Nf,Nseg] = xspectrum_band(data, Fs, deltaf, varf, Nw, properties);
        PSD = [];
    else
        [Svv_channel,Nf,Nseg,PSD] = xspectrum(data, Fs, Fmax, deltaf, varf, Nw, properties);
    end
elseif(isequal(lower(properties.sensor_params.method.value),"thomson"))
    [Svv_channel,Nf,Nseg,PSD] = xspectrum_thomson(data, Fs, Fmax, deltaf, properties);
else

end
disp('-->> Applying average reference.');
for jj = 1:Nf
    [Svv_channel(:,:,jj),Lvj] = applying_reference(Svv_channel(:,:,jj),Lvj);    % applying average reference...
end
end