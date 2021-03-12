function [Svv,F,Ns,PSD] = xspectrum(data, Fs, Fm, deltaf, varf, Nw, varargin)
% xspectrum estimates the Cross Spectrum of the input M/EEG data by Hilbert Transform
%%
% =============================================================================
% This function is part of the BC-VARETA toolbox:
% https://github.com/CCC-members/BC-VARETA_Toolbox
% =============================================================================@
%
%% Authors:
%   Pedro A. Valdes-Sosa, 2010-2021
%   Alberto Taboada-Crispi, 2016
%   Deirel Paz-Linares, 2017-2021
%   Eduardo Gonzalez-Moreira, 2017-2018
%   Ariosky Areces-Gonzalez, 2018-2021
%
%
%% Inputs:
%    data     = M/EEG data matrix, in which every row is a channel
%    Fs       = sampling frequency (in Hz)
%    Fm       = maximun frequency (in Hz) in the estimated spectrum
%    deltaf   = frequency resolution
%    use_gpu  = key:'use_gpu' value: true or false
%% Outputs:
%    PSD      = estimated power spectral density of input EEG data
%    Svv      = estimated cross spectrum of input EEG data
%    Ns       = number of segments in which the EEG signal is wrapped
%    Svv      =
%
%
%**************************************************************************
%% Initialization oF variables...
for i=1:2:length(varargin)
    eval([varargin{i} '=  varargin{(i+1)};'])
end
use_seg     = true;
[Nc,Nt]     = size(data);
deltat      = 1/(2*varf); % sliding window for hilbert envelope (adjusted by the response of the gaussian filter)
Nt_seg      = floor(Fs/deltaf); % number of time points in a segments
Nseg        = fix(Nt/Nt_seg); % number of segments 
data        = data(:,1:(Nseg*Nt_seg));
data        = reshape(data,Nc,Nt_seg,Nseg);
data        = repmat(data,1,1,1,2*Nw);
e           = dpss(Nt_seg,Nw);                                % discrete prolate spheroidal (Slepian) sequences
e           = reshape(e,[1,Nt_seg,1,2*Nw]);
e           = repmat(e,Nc,1,Nseg,1);
data        = data.*e;
Nseg        = Nseg*2*Nw; % augmented number of segments 
data        = reshape(data,Nc,Nt_seg,Nseg);
Ns          = Nseg*Nt_seg; % sample number
F           = 0:deltaf:Fm; % frequency vector
Nf          = length(F);
Svv         = zeros(Nc,Nc,Nf);
if(~exist('app_properties','var'))
    app_properties.run_bash_mode.use_gpu = false;
end
use_gpu     = app_properties.run_bash_mode.use_gpu;
use_gpu = false;
if(use_gpu)
    Svv_gpu = gpuArray(Svv);
end
%% Estimation of the Cross Spectrum...
fprintf(1,strcat('-->> Computing cross-spectra: %3d%%\n'),0);
if(use_gpu)
    W = gather(fft(gpuArray(data),[],2));
else
    W = fft(data,[],2);
end
W   = cat(2,W,conj(flip(W,2)));
for freq = 1:Nf
    deltaFilt               = filt_data(W,Fs,F(freq),varf,'use_gpu',use_gpu);
    dataEnv                 = envelope_data(deltaFilt,Fs,deltat,'use_seg',use_seg,'use_gpu',use_gpu);
    clearvars deltaFilt;
    if(use_gpu)
        dataEnv             = reshape(dataEnv,Nc,Nt_seg*Nseg).';
        Svv_gpu(:,:,freq)   = cov(dataEnv);
        Svv_gpu(:,:,freq)   = (Svv_gpu(:,:,freq) + Svv_gpu(:,:,freq)')/2;
    else
        dataEnv             = reshape(dataEnv,Nc,Nt_seg*Nseg).';
        Svv(:,:,freq)       = cov(dataEnv);
        Svv(:,:,freq)       = (Svv(:,:,freq) + Svv(:,:,freq)')/2;
    end
    clearvars dataEnv;
    fprintf(1,'\b\b\b\b%3.0f%%',(freq/Nf)*100);
end
fprintf(1,'\n');
if(use_gpu)
    Svv = gather(Svv_gpu);
end

%% Estimation of Power Spectral Density (PSD)...
PSD = zeros(Nc,Nf);
for freq = 1:Nf
    PSD(:,freq) = diag(squeeze(abs(Svv(:,:,freq))));
end

end