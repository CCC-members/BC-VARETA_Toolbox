function [Svv,peak_pos,subject,properties] = get_band(Svv,PSD,band,subject,properties)
%GET_BAND Summary of this function goes here
%   Detailed explanation goes here


% Authors:
% - Deirel Paz Linares
% - Eduardo Gonzalez Moreira
% - Pedro A. Valdes Sosa

% Date: March 16, 2019

% Updates
% - Ariosky Areces Gonzalez

% Date: March 20, 2019

%%
%% Preparing params
%%
Fm      = properties.sensor_params.max_freq.value;
deltaf  = properties.sensor_params.freq_resol.value;      % frequency resolution
F       = 0:deltaf:Fm; % frequency vector
Nf      = length(F);

disp(strcat("BC-V-->> Sensor level for frequency band: " , band.str_band));
text_level      = 'Sensor_level';
if(properties.general_params.run_by_trial.value)
    trial_name  = properties.trial_name;
    pathname    = fullfile(subject.subject_path,trial_name,text_level,band.name);
else
    pathname    = fullfile(subject.subject_path,text_level,band.name);
end
if(~isfolder(pathname))
    mkdir(pathname);
end
properties.pathname = pathname;

if(isfield(band,'f_bin'))
    [f1,nf1] = min(abs(F - band.f_bin));
    [f2,nf2] = min(abs(F - band.f_bin));
else
    [f1,nf1] = min(abs(F - band.f_start));
    [f2,nf2] = min(abs(F - band.f_end));
end

%%
%% Plot Power Spectral Density
%%

PSD_log = 10*log10(abs(PSD));
min_psd = min(PSD_log(:));
max_psd = max(PSD_log(:));
plot_peak = min_psd*ones(Nf,1);
peak_pos = nf1:nf2;
properties.peak_pos = peak_pos;
Svv = mean(Svv(:,:,peak_pos),3);
fig_title = strcat("Power Spectral Density - ", band.name,'_',string(band.f_start),'Hz-',string(band.f_end),'Hz');
if(~isfile(fullfile(properties.pathname,strcat(fig_title,'.fig'))))
    % just for the plot
    [f1,nf1] = min(abs(F - band.f_start));
    [f2,nf2] = min(abs(F - band.f_end));
    plot_peak(nf1:nf2) = max_psd;
    
    
    figure_band = figure('Color','w','Name',fig_title,'NumberTitle','off');
    
    %     define_ico(figure_band);
    hold on;
    plot(F,PSD_log);
    plot(F,plot_peak,'--k');
    set(gca,'Color','w','XColor','k','YColor','k');
    ylabel('PSD (dB)','Color','k');
    xlabel('Freq. (Hz)','Color','k');
    title(strcat("Power Spectral Density - ",band.name, " band"),'Color','k');
    
    text_cross = strcat(string(band.f_start),"Hz - ", string(band.f_end), "Hz");
    text(band.f_end,max_psd*0.9,text_cross,'Color','k','FontSize',12,'HorizontalAlignment','center');
    
    pause(1e-10);
    
    % Saving figure band
    saveas(figure_band,fullfile(properties.pathname,strcat(fig_title,'.fig')));
    close(figure_band);
end

end

