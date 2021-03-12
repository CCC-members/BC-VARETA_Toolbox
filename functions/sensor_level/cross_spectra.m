function [Svv_channel,K_vK,PSD,Nf,F,Nseg] = cross_spectra(data, Fs, Fm, deltaf, K_vK, varf, Nw, varargin)


% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019

% Updates
% - Ariosky Areces Gonzalez

% Date: January 30, 2021

%% estimating cross-spectra...
for i=1:2:length(varargin)
    eval([varargin{i} '=  varargin{(i+1)};'])
end
% Remove fieldtrip path for override functions 
warning off;
rmpath(genpath(fullfile('external/fieldtrip')));
warning on;
% estimates the Cross Spectrum of the input M/EEG data
if(exist('app_properties','var'))    
    [Svv_channel,F,Nseg,PSD] = xspectrum(data, Fs, Fm , deltaf, varf, Nw, 'app_properties', app_properties); 
else
    [Svv_channel,F,Nseg,PSD] = xspectrum(data, Fs, Fm , deltaf, varf, Nw);
end
disp('-->> Applying average reference.');
Nf = length(F);
for jj = 1:Nf
    [Svv_channel(:,:,jj),K_vK] = applying_reference(Svv_channel(:,:,jj),K_vK);    % applying average reference...
end

end

