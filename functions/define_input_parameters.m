function [result,properties] = define_input_parameters(properties)
result = 1;
%%
%% Defining subject folder
%%
tittle = "Please select the root subjects folder to analyze.";
root_path = uigetdir('tittle',tittle);
if(root_path==0)
    result = 'canceled';
    return;
end
properties.general_params.BCV_input_dir = root_path;
saveJSON(properties,fullfile('bcv_properties','bcv_properties.json'));

subjects = dir(fullfile(root_path,'**','subject.mat'));
if(~isempty(subjects))
    if(length(subjects)>1)
        answer = questdlg({'There is many subject on slected folder',' Do you want to process all subjects?'}, ...
            'Information', ...
            'Yes','No','default');
        % Handle response
        switch answer
            case 'No'
                result = 'canceled';
                return;
            case ''
                return;
        end
    end
    for i=1:length(subjects)
        subject_file = subjects(i);
        [subject,checked,error_msg_array] = ischecked_subject_data(subject_file);
        if(checked)
            break;
        else
            fprintf(2,strcat('\nBC-V-->> Error: The selected subject:',subject_file.name,' do not have a correcxt structure.\n'));
            fprintf(2,strcat('\nBC-V-->> Error: Finding other subject in folder.\n'));
        end
    end
else
    fprintf(2,strcat('\nBC-V-->> Error: The selected folder do not have any subject to process.\n'));
    disp(root_path);
    fprintf(2,strcat('BC-V-->> Error: Please other folder and start the process again.\n'));
    result = 'canceled';
    return;
end

%%
%% Defining data params
%%
guiHandle = freqresol_maxfreq_samplfreq_guide;
disp('-----Waiting for Windows number and frequency''s resolution------');
uiwait(guiHandle.SpectralpropertiesUIFigure);
if( isvalid(guiHandle) && ~guiHandle.canceled)
    delete(guiHandle);
else
    fprintf(2,'-----------Canceled by User------------');
    delete(guiHandle);
    result = 'canceled';
    return;
end

%%
%%  Ploting cross-spectra first subject
%%
properties = jsondecode(fileread(fullfile('bcv_properties','bcv_properties.json')));
disp('>> Estimating cross-spectra for M/EEG data...');
Fs = properties.samp_freq.value; % sampling frequency
Fm = properties.max_freq.value; % maximum frequency
deltaf = properties.freq_resol.value; % frequency resolution
Nw = properties.win_order.value;
[Svv_channel,Ke,PSD,Nf,F,Nseg] = cross_spectra(subject.data,Fs,Fm,deltaf,subject.Ke,Nw);


if(isequal(subject.modality,'EEG'))
    PSD_log = 10*log10(abs(PSD));
else
    PSD_log = abs(PSD);
end
min_psd = min(PSD_log(:));

plot_peak = min_psd*ones(Nf,1);

figure_cross = figure('Color','w','Name','Power Spectral Density','NumberTitle','off');
define_ico(figure_cross);
hold on;
plot(F,PSD_log);
plot(F,plot_peak,'--k');
set(gca,'Color','w','XColor','k','YColor','k');
ylabel('PSD (dB)','Color','k');
xlabel('Freq. (Hz)','Color','k');
title('Power Spectral Density','Color','k');
pause(1e-10);

%%
%% Defining frequency bands
%%
guiHandle = frequency_bands_guide;
disp('------Waitintg for frequency_bands------');
uiwait(guiHandle.UIFigure);
%waitfor(guiHandle);
if(isvalid(guiHandle) && ~guiHandle.canceled)
    disp('finishing frequencies_band...');
    delete(guiHandle);
else
    fprintf(2,'-----------Canceled by User------------\n');
    delete(guiHandle);
    result = 'canceled';
    if(exist('figure_cross','var'))
        delete(figure_cross);
    end
    return;
end
if(exist('figure_cross','var'))
    delete(figure_cross);
end

%%
%%  Defining HIGGS params
%%
guiHandle = hhgm_params_guide;
disp('-----Waiting for H-HHGM Parameters------');
uiwait(guiHandle.HHHGMParametersUIFigure);
if( isvalid(guiHandle) && ~guiHandle.canceled)
    delete(guiHandle);
else
    fprintf(2,'-----------Canceled by User------------');
    delete(guiHandle);
    result = 'canceled';
    return;
end

properties.general_params       = jsondecode(fileread(fullfile('bcv_properties','bcv_general_params.json')));
properties.sensor_params        = jsondecode(fileread(fullfile('bcv_properties','bcv_sensor_params.json')));
properties.activation_params    = jsondecode(fileread(fullfile('bcv_properties','bcv_activation_params.json')));
properties.connectivity_params  = jsondecode(fileread(fullfile('bcv_properties','bcv_connectivity_params.json')));


end

