function [Svv,Nf,Ns,PSD] = xspectrum_thomson(data,Fs,Fm,deltaf,properties)
% xspectrum estimates the Cross Spectrum of the input M/EEG data
% Inputs:
%    data     = M/EEG data matrix, in which every row is a channel
%    Fs       = sampling frequency (in Hz)
%    Fm       = maximun frequency (in Hz) in the estimated spectrum
%    deltaf   = frequency resolution
% Outputs:
%    PSD      = estimated power spectral density of input EEG data
%    Svv      = estimated cross spectrum of input EEG data
%    Ns       = number of segments in which the EEG signal is wrapped
%
%
%% 
% =============================================================================
% This function is part of the BC-VARETA toolbox:
% https://github.com/egmoreira/BC-VARETA-toolbox
% =============================================================================@
%
% Authors:
% Pedro A. Valdes-Sosa, 2010-2018
% Alberto Taboada-Crispi, 2016
% Deirel Paz-Linares, 2017-2018
% Eduardo Gonzalez-Moreira, 2017-2018
%
%**************************************************************************
%% Initialization oF variables...
NFFT     = round(Fs/deltaf);                            % number of time points per window
Nw       = 1;                                           % number of windows for Thomson spectral estimate
F        = 0:deltaf:Fm;                                 % frequency vector
%% Estimation of the Cross Spectrum...
e       = dpss(NFFT,Nw);                                % discrete prolate spheroidal (Slepian) sequences
e       = reshape(e,[1,NFFT,2*Nw]);
[Nc,Ns] = size(data);                                   % number of channels (rows) and samples (columns)
Ns      = fix(Ns/NFFT);                                 % number of segments in which the EEG signal is wrapped
Ns      = max(1,Ns);

if NFFT > Ns
    data = [data zeros(Nc,NFFT-Ns)];                    % zero padding
end
fprintf(1,strcat('-->> Computing cross-spectra: %3d%%\n'),0);
data(:,Ns*NFFT+1:end)   = [];                             % discards samples after Ns*NFFT
data                    = reshape(data,Nc,NFFT,Ns);       % 'resized' EEG data
Nf                      = length(F);                      % length of F vector
Svv                     = zeros(Nc,Nc,Nf);                % allocated matrix for the cross spectrum
if(getGlobalGuimode())
    process_waitbar = waitbar(0,'Please wait...','windowstyle', 'modal');
    frames = java.awt.Frame.getFrames();
    frames(end).setAlwaysOnTop(1);
end
for k = 1:Ns
    w = data(:,:,k);                                    % k-th window
    W = repmat(w,[1,1,2*Nw]).*repmat(e,[Nc,1,1]);       % multiplied by Slepian seq
    W = fft(W,[],2);                                    % FFT
    W = W(:,1:Nf,:);                                    % pruning values of the FFT
    for i=1:Nf
        Svv(:,:,i) = Svv(:,:,i)+cov(squeeze(W(:,i,:)).',1);
    end
    fprintf(1,'\b\b\b\b%3.0f%%',(k/Ns)*100);
    if(getGlobalGuimode())
        waitbar((k)/(Ns),process_waitbar,strcat("Computing cross-spectra: ",num2str(fix((k/Ns)*100)),"%"));
    end
end
fprintf(1,'\n');
Svv = Svv/Ns;                                           % normalizing
%% Estimation of Power Spectral Density (PSD)...
PSD = zeros(Nc,Nf);                                     
for freq = 1:Nf
    PSD(:,freq) = diag(squeeze(abs(Svv(:,:,freq))));
end
if(getGlobalGuimode())
    delete(process_waitbar);
end
end