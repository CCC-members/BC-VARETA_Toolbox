function [subject,checked,error_msg_array] = checked_subject_data(subject_file_info,properties)

base_path = properties.general_params.bcv_workspace.BCV_input_dir;
checked         = true;
subject         = struct;
subject_info    = jsondecode(fileread(fullfile(base_path,subject_file_info.SubID,strcat(subject_file_info.SubID,'.json'))));
subject.name    = subject_info.name;
subject_dir        = fullfile(base_path,subject_file_info.SubID);
error_msg_array = [];

if(~isequal(subject_file_info.Status,'Completed'))
    error_msg       = strcat("The subject information is not completed");
    error_msg_array = [error_msg_array; error_msg];
    checked         = false;
    return;
end

if(~isfield(subject_info,'name')...
        || ~isfield(subject_info,'meeg_dir')...
        || ~isfield(subject_info,'leadfield_dir')...
        || ~isfield(subject_info,'sourcemodel_dir')...
        || ~isfield(subject_info,'channel_dir')...
        || ~isfield(subject_info,'headmodel_dir'))
    checked = false;
    return;
else
    for i=1:length(subject_info.meeg_dir)
        if(contains(subject_info.meeg_dir{i},properties.general_params.dataset.task.value))
            subject_info.meeg_dir = subject_info.meeg_dir{i};
            subject_info.meeg_dir       = replace(subject_info.meeg_dir,'\','/');
            break;
        end
    end
    subject_info.channel_dir    = replace(subject_info.channel_dir,'\','/');
    subject_info.leadfield_dir.AQCI  = replace(subject_info.leadfield_dir.AQCI,'\','/');
    subject_info.leadfield_dir.leadfield  = replace(subject_info.leadfield_dir.leadfield,'\','/');
    subject_info.headmodel_dir.scalp       = replace(subject_info.headmodel_dir.scalp,'\','/');
    subject_info.headmodel_dir.outerskull       = replace(subject_info.headmodel_dir.outerskull,'\','/');
    subject_info.headmodel_dir.innerskull       = replace(subject_info.headmodel_dir.innerskull,'\','/');
    subject_info.sourcemodel_dir      = replace(subject_info.sourcemodel_dir,'\','/');
end
disp("=====================================================================");
disp("=====================================================================");
disp(strcat('BC-V-->>Processing subject:',subject_info.name));
disp("=====================================================================");
disp("=====================================================================");


if(~isfile(fullfile(subject_dir,subject_info.meeg_dir)))
    error_msg       = strcat("The meeg file do not exist");
    error_msg_array = [error_msg_array; error_msg];
    checked = false;
end
if(~isfile(fullfile(subject_dir,subject_info.leadfield_dir.leadfield)))
    error_msg       = strcat("The leadfield file do not exist");
    error_msg_array = [error_msg_array; error_msg];
    checked = false;
end
if(~isfile(fullfile(subject_dir,subject_info.channel_dir)))
    error_msg       = strcat("The channel file do not exist");
    error_msg_array = [error_msg_array; error_msg];
    checked = false;
end
if(~isfile(fullfile(subject_dir,subject_info.sourcemodel_dir)))
    error_msg       = strcat("The surface file do not exist");
    error_msg_array = [error_msg_array; error_msg];
    checked = false;
end
if(~isfile(fullfile(subject_dir,subject_info.headmodel_dir.scalp)))
    error_msg       = strcat("The scalp file do not exist");
    error_msg_array = [error_msg_array; error_msg];
    checked = false;
end
if(~isfile(fullfile(subject_dir,subject_info.headmodel_dir.innerskull)))
    error_msg       = strcat("The innerskull file do not exist");
    error_msg_array = [error_msg_array; error_msg];
    checked = false;
end
if(~isfile(fullfile(subject_dir,subject_info.headmodel_dir.outerskull)))
    error_msg       = strcat("The outerskull file do not exist");
    error_msg_array = [error_msg_array; error_msg];
    checked = false;
end

if(isempty(error_msg_array))
    Scortex     = load(fullfile(subject_dir,subject_info.sourcemodel_dir));
    Cdata       = load(fullfile(subject_dir,subject_info.channel_dir));
    Shead       = load(fullfile(subject_dir,subject_info.headmodel_dir.scalp));
    Sout        = load(fullfile(subject_dir,subject_info.headmodel_dir.outerskull));
    Sinn        = load(fullfile(subject_dir,subject_info.headmodel_dir.innerskull));
    Headmodels  = load(fullfile(subject_dir,subject_info.leadfield_dir.leadfield));
    try
        MEEG        = load(fullfile(subject_dir,subject_info.meeg_dir));
    catch
        error_msg       = strcat("Errer loading EEG file");
        error_msg_array = [error_msg_array; error_msg];
        checked = false;
        return;
    end
    if(properties.general_params.run_by_trial.value && ~iscell(MEEG.data))
        error_msg       = strcat("Process run by trial and the data format is not a cellarray");
        error_msg_array = [error_msg_array; error_msg];
        checked = false;
        return;
    end

    subject.name            = subject_info.name;
    subject.modality        = 'EEG';
    subject.Headmodel       = Headmodels.HeadModel(Headmodels.iHeadModel);
    subject.Scortex         = Scortex.Sc(Scortex.iCortex);
    if(isfield(Scortex,'sub_to_FSAve'))
        subject.sub_to_FSAve    = Scortex.sub_to_FSAve;
    else
        subject.sub_to_FSAve    = [];
    end
    subject.Shead           = Shead;
    subject.Cdata           = Cdata;
    subject.Sinn            = Sinn;
    subject.Sout            = Sout;
    subject.MEEG            = MEEG;
end
end