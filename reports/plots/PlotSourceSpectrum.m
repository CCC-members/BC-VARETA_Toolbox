function PlotSourceSpectrum(reportPath,dataset, J3D)

if(~isfolder(reportPath))
    mkdir(reportPath);
end

frequencies  = dataset.sensor_params.frequencies;
if(isfield(frequencies,'f_bin'))
    freqs = [frequencies.f_bin];
else
    maxFreq = dataset.sensor_params.max_freq.value;
    freqRes = dataset.sensor_params.freq_resol.value;
    freqs = 0:freqRes:maxFreq;
end

fig = figure("Name","Sensor PSD","Color","w", 'Position', get(0, 'Screensize'), 'Visible','on');

PSD_log = 10*log10(abs(J3D));
plot(freqs,PSD_log);
xlabel('Freq. (Hz)');
ylabel('PSD (dB)');
title('Power Spectral Density');

export_fig(fullfile(reportPath,'sensorPSD'),'-transparent','-png');
close(fig);
end
