function PlotSensorSpectrum(reportPath,dataset, nodeData)

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
spectrum = load(fullfile(dataset.Path,nodeData.SubID,'Generals/Functional/Data_spectrum.mat'));
fig = figure("Name","Sensor PSD","Color","w", 'Position', get(0, 'Screensize'), 'Visible','off');

PSD_log = 10*log10(abs(spectrum.PSD));
plot(freqs,PSD_log);
xlabel('Freq. (Hz)');
ylabel('PSD (dB)');
title('Power Spectral Density');

export_fig(fullfile(reportPath,'sensorPSD'),'-transparent','-png');
close(fig);
end

