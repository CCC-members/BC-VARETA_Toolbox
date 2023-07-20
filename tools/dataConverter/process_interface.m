function [process_error] = process_interface(properties, reject_subjects)
% Description here
%
%
%
% Author:
% - Ariosky Areces Gonzalez
% - Deirel Paz Linares
%%

%%
%% Preparing selected protocol
%%
process_error           = [];
modality                = properties.general_params.modality;
subjects_process_error  = [];
subjects_processed      = [];
report_path             = properties.general_params.reports.output_path;
general_params          = properties.general_params;
preprocessed_params     = properties.prep_data_params;
channel_label_file      = preprocessed_params.channel_label_file;
clean_config            = preprocessed_params.clean_data;
data_params             = preprocessed_params.data_config;

disp(strcat('-->> Data Source:  ', data_params.base_path ));
bcv_path = data_params.base_path;
subjects = dir(bcv_path);
subjects(ismember( {subjects.name}, {'.', '..','derivatives'})) = [];  %remove . and ..
subjects([subjects.isdir] == 0) = [];  %remove . and ..
subjects(ismember( {subjects.name}, reject_subjects)) = [];

for i=1:length(subjects)
    subject_name = subjects(i).name;
    subID = subject_name;
    disp(strcat("-->> Processing subject: ", subID));
    disp('==========================================================================');
    
    %%
    %%  Starting report
    %%
    % Creating a report
    report_name = subID;
    f_report('New');
    % Add a title to the report
    f_report('Title',strcat("EEG preprocessing: ",subID));
    % Add a title to the report
    f_report('Header',strcat("EEG preprocessing: ",subID));
         
    sub_report = fullfile(report_path,subID);
    if(~isfolder(sub_report))
        mkdir(sub_report);
    end
    
    %%
    %% Processing MEG/EEG file
    %%    
    data_path = dir(fullfile(data_params.base_path,subID,'**',strrep(data_params.file_location,'SubID',subID)));
    data_path = fullfile(data_path.folder,data_path.name);
    disp ("-->> Genering MEG/EEG file");
    if(isequal(modality,'EEG'))
        data_type    = data_params.format;
        if(~isequal(channel_label_file,"none") && ~isempty(channel_label_file))
            user_labels = jsondecode(fileread(channel_label_file));
        else
            user_labels = [];
        end
        if(~isequal(data_params.electrodes_file,"none") && ~isempty(data_params.electrodes_file))
            filepath = strrep(data_params.electrodes_file,'SubID',subID);
            base_path =  data_params.base_path;
            electrodes_file = dir(fullfile(base_path,subID,'**',filepath));
            electrodes_file = fullfile(electrodes_file.folder,electrodes_file.name);
            if(isfile(electrodes_file))
                electrodes = tsvread(electrodes_file);
                user_labels = electrodes.name;
            end
        end
        if(~isequal(data_params.derivatives_file,"none") && ~isempty(data_params.derivatives_file))
            file_name = strrep(data_params.derivatives_file,'SubID',subID);
            derivatives_file = dir(fullfile(data_params.base_path,'derivatives','**',subID,'eeg',file_name));
            derivatives_file = fullfile(derivatives_file.folder,derivatives_file.name);
            if(isfile(derivatives_file))
                derivatives = tsvread(derivatives_file);
            else
                derivatives = [];
            end
        else
            derivatives = [];
        end
        
        clean               = clean_config.run;
        toolbox_path        = clean_config.toolbox_path;
        min_freq            = clean_config.min_freq;
        max_freq            = clean_config.max_freq;
        chan_action         = clean_config.rej_or_interp_chan.action;
        select_events       = clean_config.select_events;
        clean_art_params    = clean_config.clean_artifacts;
        decompose_ica       = clean_config.decompose_ica;  
        freq_list           = clean_config.freq_list';
        electrode_file      = clean_config.electrode_file;
        
        MEEGs               = eeglab_preproc(subID, data_path, data_type,toolbox_path,...
        'clean', clean, 'verbosity', true, 'max_freq', max_freq,'min_freq', min_freq, ...
        'electrode_file', electrode_file, 'labels', user_labels, 'select_events', ...
        select_events, 'derivatives', derivatives, 'freq_list', freq_list,'save_path',...
        sub_report, 'chan_action', chan_action,'clean_art_params', clean_art_params,...
        'decompose_ica', decompose_ica);
    elseif(isequal(modality,'MEG'))
        MEEGs = import_meg_format(subID, preprocessed_params, data_path);
    else
        MEEGs = load(data_path);
        MEEGs = MEEGs.data_struct;   
    end
    
    for j=1:length(MEEGs)
        MEEG                    = MEEGs(j);
        %%
        %% Filter Channels and LeadField by Preprocessed MEEG
        %%
        base_path               = fullfile(general_params.bcv_config.base_path);
        if(general_params.bcv_config.anat_template.use_template)
            template_name       = general_params.bcv_config.anat_template.template_name;
            bcv_path            = fullfile(base_path,template_name);
        else
            bcv_path            = fullfile(base_path,subID);
        end
        subject_info            = load(fullfile(bcv_path,'subject.mat'));
        Cdata                   = load(fullfile(bcv_path,subject_info.channel_dir));
        HeadModels              = load(fullfile(bcv_path,subject_info.leadfield_dir));
        if(isequal(modality,'EEG'))
            labels              = {MEEG.chanlocs(:).labels};
        elseif(isequal(modality,'MEG'))
            labels              = MEEG.labels;
        else
            labels              = MEEG.dnames;
        end
        [Cdata, HeadModel]      = filter_structural_result_by_preproc_data(labels, Cdata, HeadModels.HeadModel);
        HeadModels.HeadModel    = HeadModel;
        
        %%
        %% Exporting files
        %%
        if(general_params.bcv_config.anat_template.use_template)
            Shead               = load(fullfile(bcv_path,subject_info.scalp_dir));
            Sout                = load(fullfile(bcv_path,subject_info.outerskull_dir));
            Sinn                = load(fullfile(bcv_path,subject_info.innerskull_dir));
            Scortex             = load(fullfile(bcv_path,subject_info.surf_dir));
            action              = 'new';
            save_output_files(action, base_path, subject_info, HeadModels, Cdata, MEEG, Shead, Sout, Sinn, Scortex);
        else
            if(~isfield(MEEG,'event_name'))
                action              = 'update';
                save_output_files(action, base_path, subject_info, HeadModels, Cdata, MEEG);
            else
                Shead               = load(fullfile(bcv_path,subject_info.scalp_dir));
                Sout                = load(fullfile(bcv_path,subject_info.outerskull_dir));
                Sinn                = load(fullfile(bcv_path,subject_info.innerskull_dir));
                Scortex             = load(fullfile(bcv_path,subject_info.surf_dir));
                action              = 'event';
                save_output_files(action, base_path, subject_info, HeadModels, Cdata, MEEG, Shead, Sout, Sinn, Scortex);
            end
        end
        disp("---------------------------------------------------------------------");
    end
    
    % Add the footer info to the report
    footer_title = 'Organization';
    text = 'All rights reserved';
    copyright = '@copy CC-Lab';
    contact = 'cc-lab@neuroinformatics-collaboratory.org';
    references = {'https://github.com/CCC-members', 'https://www.neuroinformatics-collaboratory.org/'};
    
    f_report('Footer', footer_title, text, 'ref', references, 'copyright', copyright, 'contactus', contact);
    disp('-->> Saving report.');
    % Export the report
    disp('-->> Exporting report.');
    FileFormat = 'html';
    f_report('Export',sub_report, report_name, FileFormat);
    disp('==========================================================================');
end
end

