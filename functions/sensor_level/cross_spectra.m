function [Svv_channel,K_vK,PSD,Nseg] = cross_spectra(subject, properties)


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
data    = subject.data;
K_vK    = subject.Ke;
Fs      = properties.spectral_params.samp_freq.value;       % sampling frequency
Fmax    = properties.spectral_params.max_freq.value;        % maximum frequency
deltaf  = properties.spectral_params.freq_resol.value;      % frequency resolution
varf    = properties.spectral_params.freq_gfiltvar.value;   % gaussian filter variance
Nw      = properties.spectral_params.win_order.value;       % Slepian windows

%% estimating cross-spectra...

% Remove fieldtrip path for override functions 
warning off;
rmpath(genpath(fullfile('external/fieldtrip')));
warning on;
% estimates the Cross Spectrum of the input M/EEG data
if(isequal(properties.sensor_params.method.value,"hilbert"))
if(~properties.general_params.run_frequency_bin.value &&  properties.general_params.run_frequency_bin.band_mean)    
    [Svv_channel,Nf,Nseg] = xspectrum_band(data, Fs, deltaf, varf, Nw, properties);
    PSD = [];
else
    [Svv_channel,Nf,Nseg,PSD] = xspectrum(data, Fs, Fmax, deltaf, varf, Nw, properties); 
end
elseif(isequal(properties.sensor_params.method.value,"thomson"))
    [Svv_channel,Nf,Nseg,PSD] = xspectrum_thomson(data, Fs, Fmax, deltaf, properties); 
else
    
end
disp('-->> Applying average reference.');
for jj = 1:Nf
    [Svv_channel(:,:,jj),K_vK] = applying_reference(Svv_channel(:,:,jj),K_vK);    % applying average reference...
end


end