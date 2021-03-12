function [Svv,subject,properties] = get_band(Svv,PSD,Nf,F,subject,properties)
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


%% Preparing params
band = properties.band;

%%
if(isequal(subject.modality,'EEG'))
    PSD_log = 10*log10(abs(PSD));
else
    PSD_log = abs(PSD);
end
min_psd = min(PSD_log(:));
max_psd = max(PSD_log(:));
plot_peak = min_psd*ones(Nf,1);
if(isfield(band,'f_bin'))
    [f1,nf1] = min(abs(F - band.f_bin));
    [f2,nf2] = min(abs(F - band.f_bin));
else
    [f1,nf1] = min(abs(F - band.f_start));
    [f2,nf2] = min(abs(F - band.f_end));
end

peak_pos = nf1:nf2;
properties.peak_pos = peak_pos;
Svv = mean(Svv(:,:,peak_pos),3);

fig_title = strcat("Power Spectral Density - ", band.name,'_',string(band.f_start),'Hz-',string(band.f_end),'Hz');
if(~isfile(fullfile(properties.pathname,strcat(fig_title,'.fig'))))
    % just for the plot
    [f1,nf1] = min(abs(F - band.f_start));
    [f2,nf2] = min(abs(F - band.f_end));
    plot_peak(nf1:nf2) = max_psd;
    
    if(properties.run_bash_mode.disabled_graphics)
        figure_band = figure('Color','w','Name',fig_title,'NumberTitle','off','visible','off');
    else
        figure_band = figure('Color','w','Name',fig_title,'NumberTitle','off');
    end
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

